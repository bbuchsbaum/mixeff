mk_phase4_fit <- function(seed = 40L, slope = TRUE, reml = TRUE) {
  set.seed(seed)
  n_subjects <- 8L
  n_per <- 5L
  subject <- factor(rep(seq_len(n_subjects), each = n_per))
  x <- rep(seq_len(n_per) - 1, n_subjects)
  z <- rep(c(0, 1), length.out = length(x))
  b0 <- rnorm(n_subjects, sd = 0.5)
  y <- 1 + 0.35 * x + 0.2 * z + b0[as.integer(subject)] +
    rnorm(length(x), sd = 0.25)
  df <- data.frame(y = y, x = x, z = z, subject = subject)
  formula <- if (slope) y ~ x + z + (1 | subject) else y ~ x + (1 | subject)
  lmm(formula, df, REML = reml, control = mm_control(verbose = -1))
}

mk_cbpp_glmm_fit <- function() {
  skip_if_not_installed("lme4")
  env <- new.env(parent = emptyenv())
  utils::data("cbpp", package = "lme4", envir = env)
  cbpp <- get("cbpp", envir = env, inherits = FALSE)
  rows <- lapply(seq_len(nrow(cbpp)), function(i) {
    successes <- cbpp$incidence[[i]]
    failures <- cbpp$size[[i]] - cbpp$incidence[[i]]
    data.frame(
      y = c(rep(1, successes), rep(0, failures)),
      period = cbpp$period[[i]],
      herd = cbpp$herd[[i]]
    )
  })
  df <- droplevels(do.call(rbind, rows))

  list(
    fit = glmm(y ~ period + (1 | herd), df, family = binomial(),
               method = "pirls_profiled", control = mm_control(verbose = -1)),
    data = df
  )
}

mk_glmm_contract_data <- function(family) {
  set.seed(416)
  n_groups <- 8L
  n_per_group <- 8L
  subject <- factor(rep(seq_len(n_groups), each = n_per_group))
  x <- rep(seq(-1, 1, length.out = n_per_group), n_groups)
  b0 <- rnorm(n_groups, sd = 0.25)[as.integer(subject)]
  eta <- 0.4 + 0.3 * x + b0

  y <- switch(
    family,
    binomial = stats::rbinom(length(x), size = 1L, prob = stats::plogis(eta)),
    poisson = stats::rpois(length(x), lambda = exp(eta)),
    gamma = stats::rgamma(length(x), shape = 8, rate = 8 / exp(eta)),
    stop(sprintf("unsupported test family `%s`", family), call. = FALSE)
  )
  data.frame(y = y, x = x, subject = subject)
}

test_that("glmm() fits expanded cbpp binomial smoke through profiled PIRLS", {
  pair <- mk_cbpp_glmm_fit()
  fit <- pair$fit
  df <- pair$data

  expect_s3_class(fit, "mm_glmm")
  expect_s3_class(fit, "mm_fit")
  expect_identical(fit$family$family, "binomial")
  expect_identical(fit$family$link, "logit")
  expect_identical(fit$method, "pirls_profiled")
  expect_equal(nobs(fit), nrow(df))
  expect_true(all(is.finite(fixef(fit))))
  expect_true(all(is.finite(fit$theta)))
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_true(is.finite(AIC(fit)))
  expect_s3_class(summary(fit), "summary.mm_glmm")
  expect_true(all(is.finite(diag(vcov(fit)))))
  expect_equal(nrow(model.matrix(fit)), nrow(df))
  ## The default fast-PIRLS path is NOT certified for fixed-effect Wald
  ## inference: summary() flags the covariance status unsupported and withholds
  ## SE/z/p (no fake certainty) rather than reporting certified glmer-style
  ## columns. Certified Wald requires method = "joint_laplace" (see
  ## test-glmm-summary-tests.R / test-glmm-confint.R).
  s_coef <- summary(fit, tests = "coefficients")
  expect_s3_class(s_coef, "summary.mm_glmm")
  expect_identical(s_coef$vcov_status$status, "unsupported")
  expect_true("Estimate" %in% colnames(s_coef$coefficients))
  expect_true(all(is.na(s_coef$coefficients[, "p.value"])))
})

test_that("glmm() family/link surface matches the upstream support contract", {
  expected <- data.frame(
    family = c(rep("binomial", 3), rep("poisson", 2), "Gamma"),
    link = c("logit", "probit", "cloglog", "log", "sqrt", "log"),
    stringsAsFactors = FALSE
  )
  supported <- mm_glmm_supported_family_link_table()
  expect_equal(supported[order(supported$family, supported$link), ],
               expected[order(expected$family, expected$link), ],
               ignore_attr = TRUE)

  cases <- list(
    list(family = binomial(link = "logit"), engine_family = "binomial"),
    list(family = binomial(link = "probit"), engine_family = "binomial"),
    list(family = binomial(link = "cloglog"), engine_family = "binomial"),
    list(family = poisson(link = "log"), engine_family = "poisson"),
    list(family = poisson(link = "sqrt"), engine_family = "poisson"),
    list(family = Gamma(link = "log"), engine_family = "gamma")
  )

  for (case in cases) {
    data_family <- case$engine_family
    data_family <- if (identical(data_family, "gamma")) "gamma" else data_family
    df <- mk_glmm_contract_data(data_family)
    fit <- glmm(y ~ x + (1 | subject), df, family = case$family,
                method = "pirls_profiled", control = mm_control(verbose = -1))
    expect_s3_class(fit, "mm_glmm")
    expect_identical(fit$family$family, case$engine_family)
    expect_identical(fit$family$link, case$family$link)
    expect_true(is.finite(as.numeric(logLik(fit))))
    expect_true(all(is.finite(fixef(fit))))
  }
})

test_that("glmm() refuses off-contract family/link pairs with a stable reason code", {
  df <- data.frame(
    y = c(0, 1, 0, 1, 1, 0),
    x = c(0, 0, 1, 1, 0, 1),
    subject = factor(rep(1:3, each = 2))
  )
  bad_cases <- list(
    gaussian_identity = gaussian(),
    binomial_log = binomial(link = "log"),
    poisson_identity = poisson(link = "identity"),
    gamma_inverse = Gamma(link = "inverse"),
    inverse_gaussian_log = inverse.gaussian(link = "log")
  )

  for (family in bad_cases) {
    err <- tryCatch(
      glmm(y ~ x + (1 | subject), df, family = family,
           control = mm_control(verbose = -1)),
      mm_inference_unavailable = function(cnd) cnd
    )
    expect_s3_class(err, "mm_inference_unavailable")
    expect_identical(err$reason_code, "unsupported_glmm_family_link")
    expect_identical(err$family, family$family)
    expect_identical(err$link, family$link)
    expect_true(all(c("family", "link") %in% names(err$supported)))
  }
})

test_that("mm_glmm revive preserves durable extractor behavior", {
  fit <- mk_cbpp_glmm_fit()$fit
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- revive(readRDS(path))

  expect_s3_class(restored, "mm_glmm")
  expect_false(fit_handle_alive(restored))
  expect_true(is.environment(restored$lazy_cache))
  expect_equal(fixef(restored), fixef(fit))
  expect_equal(ranef(restored), ranef(fit))
  expect_equal(VarCorr(restored)$table, VarCorr(fit)$table)
  # GLMM predict() is now available (R-side plug-in); revive must preserve it,
  # so the revived fit reproduces the original predictions.
  expect_equal(stats::predict(restored), stats::predict(fit))
})

test_that("glmm() validates family metadata and labelled joint backend", {
  df <- data.frame(
    y = c(0, 1, 0, 1, 1, 0),
    x = c(0, 0, 1, 1, 0, 1),
    subject = factor(rep(1:3, each = 2))
  )
  joint <- glmm(y ~ x + (1 | subject), df, family = binomial(),
                method = "joint_laplace", nAGQ = 1L,
                control = mm_control(verbose = -1))
  expect_s3_class(joint, "mm_glmm")
  expect_equal(joint$method, "joint_laplace")
  expect_equal(joint$nAGQ, 1L)
  expect_true(is.finite(as.numeric(logLik(joint))))
  expect_error(
    glmm(y ~ x + (1 | subject), df, family = binomial(),
         method = "joint_laplace", nAGQ = 3L,
         control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  unsupported <- tryCatch(
    glmm(y ~ x + (1 | subject), df, family = gaussian(),
         control = mm_control(verbose = -1)),
    mm_inference_unavailable = function(cnd) cnd
  )
  expect_s3_class(unsupported, "mm_inference_unavailable")
  expect_identical(unsupported$reason_code, "unsupported_glmm_family_link")
})

test_that("simulate.mm_lmm is reproducible and refit() uses the stored model", {
  fit <- mk_phase4_fit()
  sims1 <- simulate(fit, nsim = 2, seed = 100)
  sims2 <- simulate(fit, nsim = 2, seed = 100)

  expect_s3_class(sims1, "data.frame")
  expect_equal(dim(sims1), c(nobs(fit), 2L))
  expect_equal(sims1, sims2)
  expect_identical(attr(sims1, "mm_method"), "r_side_gaussian_parametric")

  ref <- refit(fit, sims1[[1L]])
  expect_s3_class(ref, "mm_lmm")
  expect_equal(nobs(ref), nobs(fit))
  expect_identical(formula(ref), formula(fit))
})

test_that("compare() and multi-model anova() refit REML fits by ML for comparison", {
  full <- mk_phase4_fit(slope = TRUE)
  reduced <- mk_phase4_fit(slope = FALSE)

  cmp <- compare(reduced, full)
  av <- stats::anova(reduced, full)

  expect_s3_class(cmp, "mm_model_comparison")
  expect_equal(nrow(cmp$table), 2L)
  expect_equal(nrow(cmp$ledger), 2L)
  expect_true(all(c("comparison_id", "formula", "fit_method", "refit",
                    "comparison_method", "fit_status", "validity_status",
                    "source") %in% names(cmp$ledger)))
  expect_identical(cmp$ledger$formula, cmp$table$formula)
  expect_true(all(cmp$ledger$fit_method == "ML"))
  expect_true(any(cmp$ledger$refit))
  expect_true(all(!cmp$table$REML))
  expect_true(any(cmp$table$refit))
  expect_true(is.finite(cmp$table$LRT[[2L]]))
  expect_true(cmp$table$lrt_available[[2L]])
  expect_identical(cmp$ledger$comparison_method[[2L]], cmp$table$method[[2L]])
  expect_identical(cmp$ledger$validity_status[[2L]], cmp$table$status[[2L]])
  expect_identical(cmp$table$comparison_class[[2L]], "nested_fixed_effects")
  expect_s3_class(av, "mm_model_comparison")
  expect_error(
    compare(reduced, full, refit_for_comparison = "error"),
    class = "mm_inference_unavailable"
  )
})

test_that("compare() consumes Rust validity rows for non-nested models", {
  set.seed(409)
  subject <- factor(rep(seq_len(8L), each = 4L))
  d <- data.frame(
    y = rnorm(32),
    x = rep(seq_len(4L), 8L),
    z = rep(c(0, 1), length.out = 32L),
    subject = subject
  )
  fit_x <- lmm(y ~ x + (1 | subject), d, REML = FALSE,
               control = mm_control(verbose = -1))
  fit_z <- lmm(y ~ z + (1 | subject), d, REML = FALSE,
               control = mm_control(verbose = -1))

  cmp <- compare(fit_x, fit_z)

  expect_s3_class(cmp, "mm_model_comparison")
  expect_true(all(is.finite(cmp$table$AIC)))
  expect_true(all(is.finite(cmp$table$BIC)))
  expect_false(cmp$table$lrt_available[[2L]])
  expect_true(is.na(cmp$table$LRT[[2L]]))
  expect_true(is.na(cmp$table$p_value[[2L]]))
  expect_identical(cmp$table$reason_code[[2L]], "non_nested_models_lrt_invalid")
  expect_identical(cmp$table$comparison_class[[2L]], "non_nested_fixed_effects")
  expect_identical(cmp$ledger$reason_code[[2L]], "non_nested_models_lrt_invalid")
  expect_identical(cmp$ledger$validity_status[[2L]], cmp$table$status[[2L]])
  expect_true(nzchar(cmp$ledger$reason[[2L]]))
  expect_error(
    compare(fit_x, fit_z, method = "lrt"),
    class = "mm_inference_unavailable"
  )
})

test_that("drop1() preserves random effects and reports deletion rows", {
  fit <- mk_phase4_fit(slope = TRUE)
  d <- stats::drop1(fit, test = "Chisq")

  expect_s3_class(d, "mm_drop1")
  expect_true(all(c("dropped", "formula", "LRT", "p_value") %in% names(d$table)))
  expect_equal(nrow(d$ledger), nrow(d$table))
  expect_true(all(c("dropped", "formula", "reference_formula",
                    "comparison_method", "statistic", "status",
                    "fit_status") %in% names(d$ledger)))
  expect_true("x" %in% d$table$dropped)
  expect_true("z" %in% d$table$dropped)
  expect_true(all(grepl("(1 | subject)", d$table$formula, fixed = TRUE)))
  expect_true(all(d$ledger$reference_formula == deparse1(formula(fit))))
})

test_that("parametric bootstrap comparison runs on a tiny nsim", {
  # parametric_bootstrap() routes through the certified Rust LRT, which
  # requires ML fits; compare() refits REML -> ML on its own.
  full <- mk_phase4_fit(slope = TRUE, reml = FALSE)
  reduced <- mk_phase4_fit(slope = FALSE, reml = FALSE)

  boot <- parametric_bootstrap(reduced, full, nsim = 2, seed = 101)
  cmp <- compare(reduced, full, method = "bootstrap", nsim = 2, seed = 101)

  expect_s3_class(boot, "mm_parametric_bootstrap")
  # Engine either certifies a p-value or refuses with a structured reason;
  # both are contract-honest. No bare mean() p-value is fabricated.
  expect_true(boot$status %in% c("available", "not_assessed"))
  if (identical(boot$status, "available")) {
    expect_true(is.finite(boot$p_value))
  } else {
    expect_true(is.na(boot$p_value))
    expect_false(is.na(boot$reason))
  }
  expect_true(is.numeric(boot$simulated))
  expect_true("successful_replicates" %in% names(boot))
  expect_s3_class(cmp$bootstrap, "mm_parametric_bootstrap")
  expect_true(cmp$table$status[[2L]] %in% c("available", "not_assessed"))
  expect_identical(cmp$table$method[[2L]], "parametric_bootstrap_lrt")
})

test_that("parametric bootstrap LRT refuses non-nested or mismatched fit data", {
  set.seed(401)
  d1 <- data.frame(
    y = rnorm(24),
    x = rnorm(24),
    z = rnorm(24),
    subject = factor(rep(1:6, each = 4))
  )
  d2 <- transform(d1, y = rnorm(24, mean = 10))

  reduced <- lmm(y ~ x + (1 | subject), d1, REML = FALSE,
                 control = mm_control(verbose = -1))
  full_other_response <- lmm(y ~ x + z + (1 | subject), d2, REML = FALSE,
                             control = mm_control(verbose = -1))
  err <- tryCatch(
    parametric_bootstrap(reduced, full_other_response, nsim = 2, seed = 1),
    mm_arg_error = function(cnd) cnd
  )
  expect_s3_class(err, "mm_arg_error")
  expect_identical(err$reason_code, "bootstrap_lrt_requires_same_observations")

  same_df_non_nested <- lmm(y ~ z + (1 | subject), d1, REML = FALSE,
                            control = mm_control(verbose = -1))
  err <- tryCatch(
    compare(reduced, same_df_non_nested, method = "bootstrap", nsim = 2),
    mm_arg_error = function(cnd) cnd
  )
  expect_s3_class(err, "mm_arg_error")
  expect_identical(err$reason_code, "bootstrap_lrt_requires_nested_models")
})

test_that("compare(method = 'bootstrap') validates nsim before dispatch", {
  full <- mk_phase4_fit(slope = TRUE, reml = FALSE)
  reduced <- mk_phase4_fit(slope = FALSE, reml = FALSE)

  expect_error(
    compare(reduced, full, method = "bootstrap", nsim = NA),
    class = "mm_arg_error"
  )
  expect_error(
    compare(reduced, full, method = "bootstrap", nsim = -1),
    class = "mm_arg_error"
  )
})

test_that("parametric_bootstrap() refuses REML fits with a typed contract reason", {
  full <- mk_phase4_fit(slope = TRUE, reml = TRUE)
  reduced <- mk_phase4_fit(slope = FALSE, reml = TRUE)
  err <- tryCatch(
    parametric_bootstrap(reduced, full, nsim = 2, seed = 1),
    mm_inference_unavailable = function(cnd) cnd
  )
  expect_s3_class(err, "mm_inference_unavailable")
  expect_identical(err$reason_code, "bootstrap_lrt_requires_ml")
})

test_that("manifest advertises shipped simulation and inference surfaces", {
  cap <- mm_formula_manifest()$capabilities
  expect_true(cap$simulate)
  expect_true(cap$inference)
  expect_true(cap$fit_glmm)
})
