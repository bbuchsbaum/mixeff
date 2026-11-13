boundary_policy_data <- function() {
  dat <- read.csv(test_path("..", "fixtures", "boundary_collapse.csv"))
  dat$subject <- factor(dat$subject)
  dat
}

test_that("all LMM inference surfaces preserve conditional results and interpretation", {
  dat <- boundary_policy_data()
  reduced <- lmm(y ~ x + (0 + x | subject), dat, control = mm_control(verbose = -1))
  for (formula in list(y ~ x + (1 + x || subject),
                       y ~ x + (1 | subject) + (0 + x | subject))) {
    fit <- lmm(formula, dat, control = mm_control(verbose = -1))
    baseline <- list(beta = fixef(fit), theta = getME(fit, "theta"),
                     ll = logLik(fit), artifact = fit$artifact,
                     auto = summary(fit)$inference$table)
    for (method in c("satterthwaite", "kenward_roger")) {
      expected <- contrast(reduced, c(0, 1), method = method)$table
      ct <- contrast(fit, c(0, 1), method = method, boundary_policy = "active_set")
      te <- test_effect(fit, "x", method = method, boundary_policy = "active_set")
      it <- inference_table(fit, method = method, boundary_policy = "active_set")
      sm <- summary(fit, method = method, boundary_policy = "active_set")
      for (obj in list(ct, te, it, sm$inference)) {
        expect_identical(obj$boundary_policy, "active_set")
        rows <- obj$table
        expect_true(all(rows$boundary_policy == "active_set"))
        expect_true(all(rows$conditioning_applied))
        expect_true(all(rows$reliability == "low"))
        expect_true(all(rows$status == "available"))
        expect_true(all(rows$method == method))
        expect_true(all(lengths(rows$fixed_components) == 1L))
        expect_true(all(vapply(rows$fixed_theta_indices, identical, logical(1), 0L)))
        expect_true(all(vapply(rows$details, function(d) {
          isTRUE(d$boundary$conditioning_applied) &&
            identical(d$boundary$policy, "active_set") &&
            identical(d$boundary$requested_method, method)
        }, logical(1))))
      }
      expect_lt(abs(ct$table$df - expected$df), 0.02)
      expect_lt(abs(ct$table$std_error - expected$std_error), 1e-6)
      expect_lt(abs(ct$table$p_value - expected$p_value), 1e-5)
      # Independent lmerTest 3.2.1 / pbkrtest 0.5.5 reference supplied upstream.
      oracle <- if (method == "satterthwaite") {
        c(df = 37.96719309, se = 0.03930555189, p = 0.01787728420)
      } else c(df = 31.34118800, se = 0.03941000993, p = 0.01919515610)
      expect_lt(abs(ct$table$df - oracle[["df"]]), 0.02)
      expect_lt(abs(ct$table$std_error - oracle[["se"]]), 1e-6)
      expect_lt(abs(ct$table$p_value - oracle[["p"]]), 1e-5)
      expect_equal(te$table$p_value, ct$table$p_value, tolerance = 1e-10)
      expect_equal(it$table$p_value[it$table$term == "x"], ct$table$p_value)
      expect_equal(sm$inference$table, it$table)
      printed <- paste(capture.output(print(sm)), collapse = "\n")
      expect_match(printed, "conditioning applied: yes", fixed = TRUE)
      expect_match(printed, "selection uncertainty is not included", fixed = TRUE)
      expect_match(printed, "Reliability: low", fixed = TRUE)
      default <- contrast(fit, c(0, 1), method = method)
      expect_true(is.na(default$table$p_value))
      expect_false(default$table$conditioning_applied)
      expect_identical(default$boundary_policy, "default")
    }
    expect_identical(fixef(fit), baseline$beta)
    expect_identical(getME(fit, "theta"), baseline$theta)
    expect_identical(logLik(fit), baseline$ll)
    expect_identical(fit$artifact, baseline$artifact)
    expect_identical(summary(fit)$inference$table, baseline$auto)
    expect_identical(inference_table(fit)$table,
                     inference_table(fit, boundary_policy = "default")$table)
  }
})

test_that("interior active-set requests retain ordinary finite-sample inference", {
  fit <- lmm(y ~ x + (0 + x | subject), boundary_policy_data(),
             control = mm_control(verbose = -1))
  expect_identical(fit_status(fit), "converged_interior")
  for (method in c("satterthwaite", "kenward_roger")) {
    active <- contrast(fit, c(0, 1), method = method, boundary_policy = "active_set")$table
    default <- contrast(fit, c(0, 1), method = method)$table
    fields <- c("estimate", "std_error", "df", "statistic", "p_value", "method",
                "status", "reliability", "reliability_reason")
    expect_equal(active[, fields], default[, fields])
    # JSON object member order is not part of the estimability contract.
    a <- active$estimability[[1L]]
    d <- default$estimability[[1L]]
    expect_identical(a$estimand, d$estimand)
    expect_setequal(names(a$assessment), names(d$assessment))
    expect_identical(a$assessment[names(d$assessment)], d$assessment)
    expect_false(active$conditioning_applied)
    expect_length(active$fixed_components[[1L]], 0)
  }
})

test_that("joint conditional term tests match an independent linear-model F test", {
  dat <- data.frame(y = rep(c(-2, -1, 0, 1, 2, 3), 10),
                    g = factor(rep(1:10, each = 6)),
                    trt = factor(rep(c("a", "b", "c"), 20)))
  fit <- lmm(y ~ trt + (1 | g), dat, control = mm_control(verbose = -1))
  ref <- anova(lm(y ~ trt, dat))
  for (method in c("satterthwaite", "kenward_roger")) {
    row <- test_effect(fit, "trt", method = method, boundary_policy = "active_set")$table
    expect_true(row$conditioning_applied)
    expect_identical(row$statistic_name, "f")
    expect_equal(row$num_df, 2)
    expect_lt(abs(row$den_df - 57), 1e-4)
    expect_lt(abs(row$statistic - ref$`F value`[1]), 1e-5)
    expect_lt(abs(row$p_value - ref$`Pr(>F)`[1]), 1e-6)
    expect_identical(row$details[[1L]]$contrast_family$restriction_rows, 2L)
    if (method == "kenward_roger") {
      expect_lt(abs(row$details[[1L]]$kenward_roger$f_scaling - 1), 1e-6)
    }
  }
})

test_that("unsupported boundary requests preserve refusal reasons without fallback", {
  dat <- boundary_policy_data()
  fits <- list(
    correlated = lmm(y ~ x + (1 + x | subject), dat, control = mm_control(verbose = -1)),
    ml = lmm(y ~ x + (1 + x || subject), dat, REML = FALSE,
             control = mm_control(verbose = -1)),
    weighted = lmm(y ~ x + (1 + x || subject), dat, weights = rep(1, nrow(dat)),
                   control = mm_control(verbose = -1))
  )
  reasons <- c(correlated = "independent scalar or diagonal", ml = "REML",
               weighted = "unweighted")
  for (name in names(fits)) {
    fit <- fits[[name]]
    for (method in c("satterthwaite", "kenward_roger")) {
      rows <- list(
        contrast(fit, c(0, 1), method = method, boundary_policy = "active_set")$table,
        test_effect(fit, "x", method = method, boundary_policy = "active_set")$table,
        inference_table(fit, method = method, boundary_policy = "active_set")$table,
        summary(fit, method = method, boundary_policy = "active_set")$inference$table
      )
      for (row in rows) {
        expect_true(all(is.na(row$p_value)))
        expect_false(any(row$conditioning_applied))
        expect_true(all(row$boundary_policy == "active_set"))
        expect_true(all(row$reason_code == "boundary_active_set_unavailable"))
        expect_true(all(grepl(reasons[[name]], row$reason, fixed = TRUE)))
        expect_false(any(row$method == "asymptotic_wald_z"))
      }
    }
  }
})

test_that("active-set policy survives saving results and fitted-state restoration", {
  fit <- lmm(y ~ x + (1 + x || subject), boundary_policy_data(),
             control = mm_control(verbose = -1))
  expected <- summary(fit, method = "satterthwaite", boundary_policy = "active_set")
  restored <- unserialize(serialize(fit, NULL))
  actual <- summary(restored, method = "satterthwaite", boundary_policy = "active_set")
  expect_equal(actual$inference$table, expected$inference$table)
  saved <- unserialize(serialize(expected, NULL))
  # Formula environments may contain process-local pointers; inference results
  # and the report must survive serialization independently of that identity.
  expect_identical(saved$inference, expected$inference)
  expect_identical(saved$boundary_policy, expected$boundary_policy)
  expect_identical(capture.output(print(saved)), capture.output(print(expected)))
})

test_that("active-set policy requires explicit finite-sample methods and LMMs", {
  dat <- data.frame(y = rep(c(-2, -1, 1, 2), 10), g = factor(rep(1:10, each = 4)))
  fit <- lmm(y ~ 1 + (1 | g), dat, control = mm_control(verbose = -1))
  for (method in c("auto", "asymptotic", "none", "bootstrap")) {
    expect_error(contrast(fit, 1, method = method, boundary_policy = "active_set"),
                 class = "mm_arg_error")
    expect_error(test_effect(fit, "(Intercept)", method = method,
                             boundary_policy = "active_set"), class = "mm_arg_error")
    expect_error(summary(fit, method = method, boundary_policy = "active_set"),
                 class = "mm_arg_error")
    if (method != "bootstrap") {
      expect_error(inference_table(fit, method = method, boundary_policy = "active_set"),
                   class = "mm_arg_error")
    }
  }
  expect_error(summary(fit, tests = "none", method = "satterthwaite",
                       boundary_policy = "active_set"), class = "mm_arg_error")
  for (policy in list("unknown", NA_character_, NULL, c("default", "active_set"))) {
    expect_error(contrast(fit, 1, boundary_policy = policy), class = "mm_arg_error")
  }
  gf <- glmm(I(y > 0) ~ 1 + (1 | g), dat, family = binomial(),
              control = mm_control(verbose = -1))
  for (call in list(quote(contrast(gf, 1, boundary_policy = "active_set")),
                    quote(test_effect(gf, "(Intercept)", boundary_policy = "active_set")),
                    quote(summary(gf, boundary_policy = "active_set")),
                    quote(inference_table(gf, boundary_policy = "active_set")))) {
    expect_error(eval(call), class = "mm_inference_unavailable")
  }
})

test_that("inference menu exposes policy and conditional interpretation without certifying it", {
  fit <- lmm(y ~ x + (1 + x || subject), boundary_policy_data(),
             control = mm_control(verbose = -1))
  options <- inference_options(fit)
  active <- subset(options$table, boundary_policy == "active_set")
  expect_setequal(active$method, c("satterthwaite", "kenward_roger"))
  expect_true(all(active$expected_status == "requires_check"))
  expect_false(any(active$current))
  expect_true(all(grepl("low reliability", active$notes, fixed = TRUE)))
  expect_true(all(grepl("boundary_policy = 'active_set'", active$r_verb, fixed = TRUE)))
})
