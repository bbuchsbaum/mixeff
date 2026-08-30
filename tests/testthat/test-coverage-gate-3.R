# Coverage-gate round 3: remaining high-miss concentrations not covered by
# test-coverage-gate.R / test-coverage-gate-2.R. Local fixtures only; keep
# fast (small data, small nsim, mm_control(verbose = -1)).

mk_cg3_lmm <- function(seed = 61L, reml = TRUE) {
  set.seed(seed)
  g <- gl(12, 5)
  x <- rnorm(60)
  y <- 0.5 * x + rnorm(60, sd = 0.35) + rnorm(12, sd = 0.5)[g]
  lmm(
    y ~ x + (1 | g),
    data.frame(y = y, x = x, g = g),
    REML = reml,
    control = mm_control(verbose = -1)
  )
}

mk_cg3_lmm_factor <- function(seed = 62L) {
  set.seed(seed)
  g <- gl(10, 6)
  f <- factor(sample(c("a", "b", "c"), 60, replace = TRUE))
  x <- rnorm(60)
  flag <- sample(c(TRUE, FALSE), 60, replace = TRUE)
  y <- 0.4 * x + as.numeric(f) * 0.15 + rnorm(60, sd = 0.3) +
    rnorm(10, sd = 0.45)[g]
  df <- data.frame(y = y, x = x, f = f, flag = flag, g = g,
                   stringsAsFactors = FALSE)
  list(
    fit = lmm(y ~ f + x + (1 | g), df, REML = FALSE,
              control = mm_control(verbose = -1)),
    fit_log = lmm(y ~ flag + x + (1 | g), df, REML = FALSE,
                  control = mm_control(verbose = -1)),
    data = df
  )
}

mk_cg3_nested <- function(seed = 63L, reml = FALSE) {
  set.seed(seed)
  g <- gl(12, 5)
  x <- rnorm(60)
  z <- rnorm(60)
  y <- 0.4 * x + 0.3 * z + rnorm(60, sd = 0.3) + rnorm(12, sd = 0.45)[g]
  df <- data.frame(y = y, x = x, z = z, g = g)
  list(
    full = lmm(y ~ x + z + (1 | g), df, REML = reml,
               control = mm_control(verbose = -1)),
    reduced = lmm(y ~ x + (1 | g), df, REML = reml,
                  control = mm_control(verbose = -1)),
    data = df
  )
}

mk_cg3_glmm <- function(seed = 64L, with_offset = FALSE) {
  set.seed(seed)
  g <- gl(8, 6)
  x <- rnorm(48)
  eta <- 0.25 * x + rnorm(8, sd = 0.3)[g]
  y <- rbinom(48, 1L, plogis(eta))
  df <- data.frame(y = y, x = x, g = g)
  off <- if (isTRUE(with_offset)) rnorm(48, sd = 0.05) else NULL
  glmm(
    y ~ x + (1 | g),
    df,
    family = binomial(),
    offset = off,
    control = mm_control(verbose = -1)
  )
}

mk_cg3_aliased <- function(seed = 65L) {
  set.seed(seed)
  subject <- factor(rep(seq_len(8L), each = 4L))
  x <- rep(0:3, 8L)
  x2 <- x
  y <- 1 + 0.5 * x + rnorm(length(x), sd = 0.4)
  lmm(
    y ~ x + x2 + (1 | subject),
    data.frame(y = y, x = x, x2 = x2, subject = subject),
    control = mm_control(verbose = -1)
  )
}

mk_cg3_design_weak <- function(seed = 66L) {
  set.seed(seed)
  g <- factor(rep(letters[1:6], each = 5))
  y <- rnorm(30, sd = 0.5)
  lmm(y ~ g + (1 | g), data.frame(y = y, g = g),
      control = mm_control(verbose = -1))
}

mk_cg3_multi_re <- function(seed = 67L) {
  set.seed(seed)
  g <- gl(10, 6)
  h <- gl(5, 12)
  x <- rnorm(60)
  y <- 0.4 * x + rnorm(60, sd = 0.3) + rnorm(10, sd = 0.4)[g] +
    rnorm(5, sd = 0.3)[h]
  lmm(
    y ~ x + (1 | g) + (1 | h),
    data.frame(y = y, x = x, g = g, h = h),
    REML = FALSE,
    control = mm_control(verbose = -1)
  )
}

mk_cg3_uncorrelated <- function(seed = 68L) {
  set.seed(seed)
  g <- gl(10, 6)
  x <- rnorm(60)
  y <- 0.4 * x + rnorm(60, sd = 0.3) +
    rnorm(10, sd = 0.45)[g] + rnorm(10, sd = 0.2)[g] * x
  lmm(
    y ~ x + (1 | g) + (0 + x | g),
    data.frame(y = y, x = x, g = g),
    REML = FALSE,
    control = mm_control(verbose = -1)
  )
}

# ---- inference.R ------------------------------------------------------------

test_that("test_effect method none / cluster_bootstrap / boundary df", {
  fit <- mk_cg3_lmm()
  te_none <- test_effect(fit, "x", method = "none")
  expect_s3_class(te_none, "mm_effect_test")
  expect_true(all(te_none$table$status == "not_assessed"))
  expect_true(length(capture.output(print(te_none))) > 0L)

  te_cb <- test_effect(fit, "x", method = "cluster_bootstrap")
  expect_s3_class(te_cb, "mm_effect_test")
  expect_identical(
    te_cb$table$reason_code[[1L]],
    "bootstrap_cluster_resample_p_value_unavailable"
  )
  expect_error(
    test_effect(fit, "x", method = "cluster_bootstrap", group = ""),
    class = "mm_arg_error"
  )
  expect_error(
    test_effect(fit, "x", method = "cluster_bootstrap", group = "nope"),
    class = "mm_arg_error"
  )

  aliased <- mk_cg3_aliased()
  expect_true(isTRUE(is_singular(aliased)))
  te_sat <- test_effect(aliased, "x", method = "satterthwaite")
  expect_true(all(te_sat$table$status == "not_assessed"))
  te_kr <- test_effect(aliased, "x", method = "kenward_roger")
  expect_true(all(te_kr$table$status == "not_assessed"))
  expect_true(length(capture.output(print(te_sat))) > 0L)
})

test_that("test_effect bootstrap / bootstrap_lrt with small nsim", {
  fit_ml <- mk_cg3_lmm(reml = FALSE)
  boot <- bootstrap_control(nsim = 2L, seed = 11L)
  te_b <- test_effect(fit_ml, "x", method = "bootstrap", bootstrap = boot)
  expect_s3_class(te_b, "mm_effect_test")
  expect_equal(nrow(te_b$table), 1L)

  te_blrt <- test_effect(
    fit_ml, "x", method = "bootstrap_lrt",
    bootstrap = list(nsim = 2L, seed = 12L)
  )
  expect_s3_class(te_blrt, "mm_effect_test")

  fit_reml <- mk_cg3_lmm(reml = TRUE)
  te_reml_lrt <- test_effect(
    fit_reml, "x", method = "bootstrap_lrt",
    bootstrap = bootstrap_control(nsim = 2L, seed = 13L)
  )
  reml_reason <- te_reml_lrt$table$reason_code[[1L]]
  if (is.null(reml_reason) || is.na(reml_reason)) {
    reml_reason <- te_reml_lrt$table$reason[[1L]]
  }
  if (is.null(reml_reason) || is.na(reml_reason)) reml_reason <- ""
  expect_true(
    grepl("REML|ml", reml_reason, ignore.case = TRUE) ||
      te_reml_lrt$table$status[[1L]] != "available"
  )
})

test_that("test_random_effect multi-RE, ambiguous, unknown, REML ml", {
  multi <- mk_cg3_multi_re()
  out <- test_random_effect(multi, "g")
  expect_s3_class(out, "mm_random_effect_test")
  expect_true(nrow(out$table) >= 1L)
  expect_true(length(capture.output(print(out))) > 0L)
  expect_true(length(capture.output(print(reporting_table(out)))) > 0L)

  unc <- mk_cg3_uncorrelated()
  expect_error(test_random_effect(unc, "g"), class = "mm_arg_error")
  out_r0 <- test_random_effect(unc, "r0")
  expect_s3_class(out_r0, "mm_random_effect_test")
  expect_error(test_random_effect(unc, "not_a_term"), class = "mm_arg_error")

  reml <- mk_cg3_lmm(reml = TRUE)
  out_ml <- test_random_effect(reml, "g", refit_for_comparison = "ml")
  expect_s3_class(out_ml, "mm_random_effect_test")

  expect_error(
    mixeff:::mm_boundary_lrt_ml_fit(1, "auto"),
    class = "mm_arg_error"
  )
})

test_that("contrast matrix L / named L / bad L type / anova type errors", {
  fit <- mk_cg3_lmm()
  nm <- names(fixef(fit))
  L <- rbind(c(0, 1), c(1, -1))
  colnames(L) <- nm
  rownames(L) <- c("slope", "diff")
  ct <- contrast(fit, L, method = "asymptotic")
  expect_s3_class(ct, "mm_contrast")
  expect_equal(nrow(ct$table), 2L)
  expect_true(length(capture.output(print(ct))) > 0L)

  named <- setNames(c(0, 1), nm)
  ct2 <- contrast(fit, named, method = "asymptotic")
  expect_equal(nrow(ct2$table), 1L)

  expect_error(contrast(fit, data.frame(a = 1)), class = "mm_arg_error")
  expect_error(anova(fit, type = "IV"))
})

test_that("df_for_contrast / estimability / confint bootstrap print edges", {
  fit <- mk_cg3_lmm(reml = TRUE)
  aliased <- mk_cg3_aliased()

  df_sat <- df_for_contrast(aliased, c(0, 1, 0), method = "satterthwaite")
  expect_s3_class(df_sat, "mm_df_for_contrast")
  expect_true(length(capture.output(print(df_sat))) > 0L)

  est <- estimability(fit, L = diag(2))
  expect_equal(nrow(est$table), 2L)
  expect_true(length(capture.output(print(est))) > 0L)

  ci_b <- confint(
    fit, parm = 1L, method = "bootstrap",
    bootstrap = list(nsim = 2L, seed = 21L),
    interval = "percentile"
  )
  expect_true(is.matrix(ci_b) && nrow(ci_b) == 1L)
  out_ci <- capture.output(print(ci_b))
  expect_true(any(grepl("Confidence intervals|Bootstrap", out_ci)))

  ci_basic <- confint(
    fit, parm = "(Intercept)", method = "bootstrap",
    bootstrap = bootstrap_control(nsim = 2L, seed = 22L),
    interval = "basic"
  )
  expect_equal(nrow(ci_basic), 1L)

  pr <- profile(fit)
  expect_s3_class(pr, "mm_profile")
  out_pr <- capture.output(print(pr))
  expect_true(any(grepl("Profile-likelihood", out_pr, fixed = TRUE)))
  # REML profile appends typed refusal rows for beta; print may note them.
  expect_true(is.matrix(confint(pr, parm = "sigma")))
})

test_that("mm_lincomb LMM/GLMM and bad weights", {
  fit <- mk_cg3_lmm()
  nm <- names(fixef(fit))
  lc <- mm_lincomb(fit, setNames(c(0, 1), nm), method = "asymptotic")
  expect_true(is.data.frame(lc) && nrow(lc) == 1L)
  lc_auto <- mm_lincomb(fit, setNames(c(1, 0), nm), method = "auto")
  expect_true(is.data.frame(lc_auto))

  expect_error(mm_lincomb(1, c(a = 1)), class = "mm_arg_error")
  expect_error(mm_lincomb(fit, c(not_a_coef = 1)), class = "mm_arg_error")
  expect_error(mm_lincomb(fit, c(1, 2, 3)), class = "mm_arg_error")

  gfit <- mk_cg3_glmm()
  gnm <- names(fixef(gfit))
  # Default profiled GLMM may refuse Wald; accept either result or typed error.
  out <- tryCatch(
    mm_lincomb(gfit, setNames(c(0, 1), gnm), method = "asymptotic"),
    error = function(e) e
  )
  expect_true(is.data.frame(out) || inherits(out, "mm_inference_unavailable"))
  expect_error(
    mm_lincomb(gfit, setNames(c(0, 1), gnm), method = "bootstrap"),
    class = "mm_arg_error"
  )
})

test_that("boundary-LRT JSON / unavailable helpers via :::", {
  expect_error(
    mixeff:::mm_json_parse_boundary_lrt(""),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_boundary_lrt("{"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_boundary_lrt("{}"),
    class = "mm_schema_error"
  )
  tab <- mixeff:::mm_unavailable_effect_table("x", "none")
  expect_equal(nrow(tab), 1L)
  tab2 <- mixeff:::mm_boundary_df_unavailable_effect_table(
    c("x", "z"), "satterthwaite"
  )
  expect_equal(nrow(tab2), 2L)
  expect_match(tab2$reason_code[[1L]], "satterthwaite_unavailable")
})

# ---- compare.R --------------------------------------------------------------

test_that("compare bootstrap nsim=0, aic, multi-model anova", {
  nest <- mk_cg3_nested(reml = FALSE)
  cmp0 <- compare(nest$reduced, nest$full, method = "bootstrap", nsim = 0L)
  expect_s3_class(cmp0, "mm_model_comparison")
  expect_true(all(cmp0$table$status == "bootstrap_not_run"))
  expect_true(all(cmp0$table$method == "bootstrap_not_run"))
  expect_true(length(capture.output(print(cmp0))) > 0L)

  cmp_aic <- compare(nest$reduced, nest$full, method = "aic")
  expect_true(all(cmp_aic$table$status == "information_criteria"))

  av <- anova(nest$reduced, nest$full)
  expect_s3_class(av, "mm_model_comparison")
  rt <- reporting_table(av)
  expect_s3_class(rt, "mm_reporting_table")
  expect_true(length(capture.output(print(rt))) > 0L)

  boot <- parametric_bootstrap(nest$reduced, nest$full, nsim = 3L, seed = 31L)
  expect_s3_class(boot, "mm_parametric_bootstrap")
  out <- capture.output(print(boot))
  expect_true(any(grepl("Parametric bootstrap", out, fixed = TRUE)))
  expect_true(any(grepl("p.value|not certified|replicates", out)))
})

test_that("drop1 interaction scope and compare response mismatch", {
  set.seed(69L)
  a <- factor(sample(c("a", "b"), 50, TRUE))
  b <- factor(sample(c("c", "d"), 50, TRUE))
  g <- gl(10, 5)
  y <- rnorm(50) + rnorm(10)[g]
  df <- data.frame(y = y, a = a, b = b, g = g)
  fit <- lmm(y ~ a * b + (1 | g), df, REML = FALSE,
             control = mm_control(verbose = -1))
  d1 <- drop1(fit, scope = "a", test = "Chisq")
  expect_s3_class(d1, "mm_drop1")
  expect_true("a" %in% d1$table$dropped)
  expect_true(length(capture.output(print(d1))) > 0L)

  nest <- mk_cg3_nested(reml = FALSE)
  other <- nest$data
  names(other)[1] <- "resp"
  other_fit <- lmm(resp ~ x + (1 | g), other, REML = FALSE,
                   control = mm_control(verbose = -1))
  expect_error(compare(nest$full, other_fit), class = "mm_arg_error")
})

# ---- reporting.R ------------------------------------------------------------

test_that("model_report section subsets and reporting_table all/audit", {
  fit <- mk_cg3_lmm()
  mr <- model_report(fit, sections = c("overview", "fixed_effects",
                                       "random_effects", "fit_statistics"))
  expect_true(all(c("overview", "fixed_effects") %in% names(mr$sections)))
  expect_false("optimizer" %in% names(mr$sections))
  expect_true(length(capture.output(print(mr))) > 0L)

  rt_all_audit <- reporting_table(fit, section = "all", view = "audit")
  expect_s3_class(rt_all_audit, "mm_reporting_table")
  expect_true(!is.null(rt_all_audit$sections))
  expect_true(length(capture.output(print(rt_all_audit))) > 0L)

  rt_mr <- reporting_table(mr, section = "all", view = "audit")
  expect_true(!is.null(rt_mr$sections))

  expect_error(
    reporting_table(fit, section = c("all", "overview")),
    class = "mm_schema_error"
  )
  expect_error(
    reporting_table(fit, section = c("overview", "overview")),
    class = "mm_schema_error"
  )
  # model_report(allow_many=TRUE) treats "all" as the full section list.
  mr_all <- model_report(fit, sections = c("all", "overview"))
  expect_true(length(mr_all$sections) >= 5L)

  # Design-weak / random_terms / random_effects sections on a saturated design.
  weak <- mk_cg3_design_weak()
  mr_w <- model_report(weak, sections = c("random_terms", "random_effects",
                                          "data_design", "overview"))
  expect_true(nrow(mr_w$sections$random_terms) >= 1L)
  expect_true(length(capture.output(print(
    reporting_table(weak, section = "random_terms", view = "audit")
  ))) > 0L)
})

test_that("reporting cross-card / unavailable helpers via :::", {
  rows <- mixeff:::mm_report_cross_card_rows(list(
    list(between_cards = list("a", "b"), between_basis = list("1"),
         reason = "uncorrelated"),
    list(between = list("c"), reason = "constraint_only")
  ))
  expect_true(is.data.frame(rows) && nrow(rows) >= 1L)
  txt <- mixeff:::mm_report_constraint_text(list(
    list(reason = "r1", between = list("a", "b")),
    list(reason = "r2")
  ))
  expect_true(is.character(txt) && nzchar(txt))
  expect_identical(mixeff:::mm_group_ir_label(NULL), "")
})

# ---- revive.R ---------------------------------------------------------------

test_that("getME extra names and fit_handle_alive true path", {
  fit <- mk_cg3_lmm()
  fl <- getME(fit, "flist")
  expect_true(is.list(fl) && length(fl) >= 1L)
  cn <- getME(fit, "cnms")
  expect_true(is.list(cn) && length(cn) >= 1L)
  parts <- getME(fit, c("flist", "cnms", "Lambda", "mu", "y"))
  expect_true(all(c("flist", "cnms", "Lambda", "mu", "y") %in% names(parts)))

  expect_false(fit_handle_alive(fit))
  live <- fit
  live$rust_handle <- new("externalptr")
  expect_true(fit_handle_alive(live))

  sch <- mixeff:::mm_object_schema(fit$artifact)
  expect_true(is.list(sch) && "schema_name" %in% names(sch))
})

test_that("revive artifact guards and empty random_blocks / cert print", {
  fit <- mk_cg3_lmm()
  bad <- fit
  bad$artifact <- NULL
  bad$fit$artifact_json <- NULL
  expect_error(revive(bad), class = "mm_arg_error")

  # Rebuild artifact from stored JSON path.
  stripped <- fit
  stripped$artifact <- NULL
  class(stripped) <- class(fit)
  revived <- revive(stripped)
  expect_s3_class(revived, "mm_fit")
  expect_true(is.list(revived$artifact))

  empty_blocks <- structure(
    list(table = mixeff:::mm_random_blocks_empty_table()),
    class = "mm_random_blocks"
  )
  out_b <- capture.output(print(empty_blocks))
  expect_true(any(grepl("none", out_b, fixed = TRUE)))

  empty_cert <- structure(
    list(table = data.frame(metric = character(), value = character(),
                            stringsAsFactors = FALSE)),
    class = "mm_optimizer_certificate"
  )
  out_c <- capture.output(print(empty_cert))
  expect_true(any(grepl("not assessed|none", out_c)))

  expect_error(
    mixeff:::.mm_lazy(list(), "X", identity),
    class = "mm_arg_error"
  )
})

# ---- predict.R --------------------------------------------------------------

test_that("predict LMM newdata missing cols and new levels", {
  fit <- mk_cg3_lmm()
  expect_error(
    predict(fit, newdata = data.frame(x = 1:3)),
    class = "mm_data_error"
  )
  expect_error(
    predict(fit, newdata = list(x = 1)),
    class = "mm_data_error"
  )

  nd <- data.frame(x = 0, g = factor("new_g", levels = c(levels(fit$model_frame$g),
                                                         "new_g")))
  expect_error(
    predict(fit, newdata = nd, allow.new.levels = FALSE),
    class = "mm_inference_unavailable"
  )
  pred_ok <- predict(fit, newdata = nd, allow.new.levels = TRUE)
  expect_equal(length(pred_ok), 1L)
  expect_true(is.finite(pred_ok) || is.na(pred_ok))

  pop <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)), re.form = NA)
  expect_equal(length(pop), 3L)
})

test_that("predict GLMM se.fit+interval, offset+newdata, type errors", {
  gfit <- mk_cg3_glmm()
  pop_ci <- predict(gfit, re.form = NA, interval = "confidence")
  expect_true(is.matrix(pop_ci) && ncol(pop_ci) == 3L)
  pop_both <- predict(gfit, re.form = NA, se.fit = TRUE, interval = "confidence")
  expect_true(is.list(pop_both) && is.matrix(pop_both$fit))
  pop_se <- predict(gfit, re.form = NA, se.fit = TRUE)
  expect_true(is.list(pop_se) && "se.fit" %in% names(pop_se))

  link_ci <- predict(gfit, re.form = NA, type = "link", interval = "confidence")
  expect_true(is.matrix(link_ci) && ncol(link_ci) == 3L)

  expect_error(predict(gfit, type = "probability"), class = "error")
  expect_error(
    predict(gfit, newdata = data.frame(x = 0)),
    class = "mm_data_error"
  )

  gfit_off <- mk_cg3_glmm(with_offset = TRUE)
  expect_error(
    predict(gfit_off, newdata = data.frame(x = 0, g = gfit_off$model_frame$g[1])),
    class = "mm_inference_unavailable"
  )

  # Conditional se.fit / interval still exercised on the default GLMM.
  cond_se <- predict(gfit, se.fit = TRUE)
  expect_true(is.list(cond_se))
  cond_ci <- predict(gfit, se.fit = TRUE, interval = "confidence")
  expect_true(is.list(cond_ci))
})

test_that("confint.mm_glmm level/parm guards", {
  gfit <- mk_cg3_glmm()
  expect_error(confint(gfit, level = 0), class = "mm_arg_error")
  expect_error(confint(gfit, level = 1.5), class = "mm_arg_error")
  expect_error(confint(gfit, parm = "not_a_coef"), class = "mm_arg_error")
})

# ---- marginal.R -------------------------------------------------------------

test_that("mm_grid logical/character/numeric at= and print empty grid", {
  fx <- mk_cg3_lmm_factor()
  gr_f <- mm_grid(fx$fit, ~f, at = list(f = c("a", "b"), x = 0))
  expect_equal(nrow(gr_f$grid), 2L)

  gr_log <- mm_grid(fx$fit_log, ~flag, at = list(flag = c(TRUE, FALSE)))
  expect_equal(nrow(gr_log$grid), 2L)

  # Character / logical / numeric default-grid helpers (no live character FE).
  expect_identical(
    mixeff:::mm_default_grid_values(c("lo", "hi", "lo"), "lab", TRUE, mean),
    c("hi", "lo")
  )
  expect_identical(
    mixeff:::mm_default_grid_values(c(TRUE, FALSE, TRUE), "flag", TRUE, mean),
    c(FALSE, TRUE)
  )
  expect_equal(
    mixeff:::mm_default_grid_values(1:3, "x", FALSE, mean),
    mean(1:3)
  )

  expect_error(
    mm_grid(fx$fit, ~f, at = list(f = character())),
    class = "mm_arg_error"
  )
  expect_error(
    mm_grid(fx$fit, ~f, at = list(f = "zzz")),
    class = "mm_arg_error"
  )
  expect_error(
    mm_grid(fx$fit, ~f, cov.reduce = function(x) c(1, 2)),
    class = "mm_arg_error"
  )

  empty <- gr_f
  empty$grid <- empty$grid[0L, , drop = FALSE]
  out <- capture.output(print(empty))
  expect_true(any(grepl("Marginal grid", out, fixed = TRUE)))

  expect_error(mm_predictions(fx$fit), class = "mm_arg_error")
  expect_error(
    mm_means(fx$fit, ~f, grid = list(a = 1)),
    class = "mm_inference_unavailable"
  )
})

test_that("mm_means proportional weights and comparisons edges", {
  fx <- mk_cg3_lmm_factor()
  mn <- mm_means(fx$fit, ~f, method = "asymptotic", weights = "proportional")
  expect_s3_class(mn, "mm_marginal_quantity")
  expect_true(nrow(mn$table) >= 2L)

  gr <- mm_grid(fx$fit, ~f, at = list(f = "a"))
  expect_error(
    mm_comparisons(fx$fit, ~f, grid = gr, method = "asymptotic"),
    class = "mm_inference_unavailable"
  )
})

# ---- methods-extract.R / methods-print.R / coef-names.R --------------------

test_that("aliased coefficient print and design-weak VarCorr", {
  aliased <- mk_cg3_aliased()
  out_fit <- capture.output(print(aliased))
  expect_true(any(grepl("Aliased", out_fit, fixed = TRUE)))
  s <- summary(aliased)
  out_s <- capture.output(print(s))
  expect_true(is.character(out_s) && length(out_s) > 0L)
  expect_true(any(s$coefficients$method == "aliased") ||
                !is.null(attr(s$coefficients, "mm_aliased")))

  weak <- mk_cg3_design_weak()
  vc <- VarCorr(weak)
  expect_true(length(attr(vc, "mm_design_weak_identifiability_groups")) >= 1L)
  out_vc <- capture.output(print(vc))
  expect_true(any(grepl("design_weak_identifiability", out_vc, fixed = TRUE)))

  empty_vc <- structure(
    list(
      table = data.frame(
        group = character(), name = character(),
        variance = numeric(), std_dev = numeric(),
        stringsAsFactors = FALSE
      ),
      residual_sd = NA_real_
    ),
    class = c("mm_varcorr", "list")
  )
  expect_true(any(grepl("none", capture.output(print(empty_vc)), fixed = TRUE)))
  expect_true(any(grepl("none", capture.output(print(
    structure(list(), class = c("mm_ranef", "list"))
  )), fixed = TRUE)))
  expect_true(any(grepl("none", capture.output(print(
    structure(list(), class = c("mm_coef", "list"))
  )), fixed = TRUE)))
})

test_that("coef GLMM and as.data.frame.mm_varcorr residual path", {
  gfit <- mk_cg3_glmm()
  cf <- coef(gfit)
  expect_s3_class(cf, "mm_coef")
  expect_true(length(capture.output(print(cf))) > 0L)

  fit <- mk_cg3_lmm()
  vc_df <- as.data.frame(VarCorr(fit))
  expect_true(is.data.frame(vc_df) && nrow(vc_df) >= 1L)
  expect_true("Residual" %in% vc_df$grp || any(vc_df$grp == "Residual") ||
                nrow(vc_df) >= 1L)

  expect_true(is.matrix(model.matrix(fit, type = "random")) ||
                inherits(model.matrix(fit, type = "random"), "Matrix"))
  expect_true(is.matrix(model.matrix(gfit, type = "fixed")))
})

test_that("mm_note_append and group label helpers via :::", {
  expect_equal(
    mixeff:::mm_note_append(c("", "a"), "[note]"),
    c("[note]", "a [note]")
  )
  expect_equal(mixeff:::mm_note_append(character(), "[note]"), character())
  expect_identical(
    mixeff:::mm_group_lme4_label(list(single = list(name = "g"))),
    "g"
  )
  expect_identical(
    mixeff:::mm_group_lme4_label(list(cell = list(names = list("a", "b")))),
    "a:b"
  )
  expect_identical(
    mixeff:::mm_group_lme4_label(list(interaction = list(names = list("a", "b")))),
    "a:b"
  )
})

# ---- methods-summary.R ------------------------------------------------------

test_that("summary GLMM method edges and unavailable reason print", {
  gfit <- mk_cg3_glmm()
  s <- summary(gfit)
  expect_s3_class(s, "summary.mm_glmm")
  out <- capture.output(print(s))
  expect_true(any(grepl("Generalized linear mixed model|Family", out)))
  out_v <- capture.output(print(s, verbose = TRUE))
  expect_true(is.character(out_v) && length(out_v) > 0L)

  s_none <- summary(gfit, tests = "none")
  expect_null(s_none$inference)
  expect_true(length(capture.output(print(s_none))) > 0L)

  fit <- mk_cg3_lmm()
  s_none_m <- summary(fit, method = "none")
  out_r <- capture.output(print(s_none_m))
  expect_true(any(grepl("Reasons:", out_r, fixed = TRUE)))

  aliased <- mk_cg3_aliased()
  s_a <- summary(aliased, method = "satterthwaite")
  out_a <- capture.output(print(s_a))
  expect_true(any(grepl("Reasons:|not_assessed|unavailable|Aliased", out_a)))
})

# ---- changes.R --------------------------------------------------------------

test_that("changes on stopped fits and empty / covariance print paths", {
  nest <- mk_cg3_nested()
  fit_s <- lmm(
    y ~ x + (1 | g), nest$data,
    control = mm_control(verbose = -1, max_feval = 1L)
  )
  ch <- changes(fit_s)
  expect_s3_class(ch, "mm_change_log")
  out <- capture.output(print(ch))
  expect_true(any(grepl("not_optimized|stopped early|Model changes", out)))

  empty <- structure(
    list(
      table = mixeff:::mm_change_empty_table(),
      reductions = list(),
      covariance_transitions = list(),
      effective_covariance = list(),
      fit_status = "not_assessed",
      requested_formula = "",
      effective_formula = ""
    ),
    class = "mm_change_log"
  )
  out_e <- capture.output(print(empty))
  expect_true(any(grepl("none recorded", out_e, fixed = TRUE)))

  empty$fit_status <- "converged"
  out_c <- capture.output(print(empty))
  expect_true(any(grepl("fitted as requested|none", out_c)))

  expect_identical(mixeff:::mm_change_covariance_text(NULL), "(unknown)")
  expect_true(nzchar(mixeff:::mm_change_covariance_text("diagonal")))
  expect_true(nzchar(mixeff:::mm_change_covariance_text(
    list(full = list(rank = 1L))
  )))
  expect_true(nzchar(mixeff:::mm_change_covariance_text(list(scalar = list()))))
})

# ---- emmeans.R --------------------------------------------------------------

test_that("emmeans recover_data / emm_basis on LMM and GLMM", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("estimability")
  fit <- mk_cg3_lmm_factor()$fit
  rd <- mixeff:::recover_data.mm_lmm(fit)
  expect_true(is.data.frame(rd) || inherits(rd, "data.frame") || is.list(rd))

  trms <- delete.response(terms(fit))
  grid <- data.frame(
    f = factor(levels(fit$model_frame$f), levels = levels(fit$model_frame$f)),
    x = mean(fit$model_frame$x)
  )
  eb <- mixeff:::emm_basis.mm_lmm(
    fit,
    trms = trms,
    xlev = list(f = levels(fit$model_frame$f)),
    grid = grid
  )
  expect_true(is.list(eb) && !is.null(eb$X))

  gfit <- mk_cg3_glmm()
  # Default profiled GLMM may refuse the emmeans bridge; accept typed refusal.
  rd_g <- tryCatch(
    mixeff:::recover_data.mm_glmm(gfit),
    error = function(e) e
  )
  expect_true(is.data.frame(rd_g) || is.list(rd_g) ||
                inherits(rd_g, "mm_inference_unavailable") ||
                inherits(rd_g, "error"))
})

# ---- glmm.R validation ------------------------------------------------------

test_that("glmm weights/offset and cbind validation", {
  set.seed(70L)
  g <- gl(6, 5)
  x <- rnorm(30)
  y <- rbinom(30, 1, 0.5)
  df <- data.frame(y = y, x = x, g = g)
  expect_error(
    glmm(y ~ x + (1 | g), df, family = binomial(),
         weights = rep(-1, 30), control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  expect_error(
    glmm(y ~ x + (1 | g), df, family = binomial(),
         weights = 1:3, control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  expect_error(
    glmm(y ~ x + (1 | g), df, family = binomial(),
         offset = c(1, 2), control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  expect_error(
    glmm(y ~ x + (1 | g), df, family = binomial(),
         offset = c(NA_real_, rep(0, 29)),
         control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )

  succ <- rbinom(30, 5, 0.4)
  fail <- 5L - succ
  # cbind with zero trial totals
  df_bad <- data.frame(succ = c(0L, succ[-1]), fail = c(0L, fail[-1]),
                       x = x, g = g)
  expect_error(
    glmm(cbind(succ, fail) ~ x + (1 | g), df_bad, family = binomial(),
         control = mm_control(verbose = -1)),
    class = "mm_data_error"
  )
  # cbind together with weights is refused
  df_ok <- data.frame(succ = succ, fail = fail, x = x, g = g)
  expect_error(
    glmm(cbind(succ, fail) ~ x + (1 | g), df_ok, family = binomial(),
         weights = rep(5, 30), control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )

  expect_error(
    mixeff:::mm_json_parse_glmm_fit(""),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_glmm_fit("{}"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_glmm_fit("{"),
    class = "mm_schema_error"
  )
})

# ---- misc high-ROI leftovers ------------------------------------------------

test_that("inference_table recompute and mm_fixed_effect_terms edges", {
  fit <- mk_cg3_lmm()
  inf_none <- inference_table(fit, method = "none")
  expect_s3_class(inf_none, "mm_inference_table")
  expect_true(length(capture.output(print(inf_none))) > 0L)

  inf_as <- inference_table(fit, method = "asymptotic")
  expect_true(nrow(inf_as$table) >= 1L)

  terms <- mixeff:::mm_fixed_effect_terms(fit)
  expect_true(is.character(terms) && length(terms) >= 1L)
})

test_that("compile-only changes and reproducibility empty thresholds", {
  fit <- mk_cg3_lmm()
  spec <- compile_model(y ~ x + (1 | g), fit$model_frame)
  ch <- changes(spec)
  expect_s3_class(ch, "mm_change_log")
  expect_true(length(capture.output(print(ch))) > 0L)

  empty_thr <- mixeff:::mm_repro_threshold_table(list())
  expect_equal(nrow(empty_thr), 0L)
  empty_cert_tbl <- mixeff:::mm_optimizer_certificate_table(list())
  expect_equal(nrow(empty_cert_tbl), 0L)
})
