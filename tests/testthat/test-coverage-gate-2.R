# Coverage-gate round 2: high-ROI paths not covered by test-coverage-gate.R.
# Local fixtures only (helpers from other test-*.R are not shared).
# Keep fast: small data, small nsim, mm_control(verbose = -1).

mk_cg2_lmm <- function(seed = 51L, reml = TRUE, factor_x = FALSE) {
  set.seed(seed)
  g <- gl(12, 5)
  x <- rnorm(60)
  f <- factor(sample(c("a", "b", "c"), 60, replace = TRUE))
  y <- 0.5 * x + as.numeric(f) * 0.2 + rnorm(60, sd = 0.35) +
    rnorm(12, sd = 0.5)[g]
  df <- data.frame(y = y, x = x, f = f, g = g)
  form <- if (isTRUE(factor_x)) y ~ f + (1 | g) else y ~ x + f + (1 | g)
  lmm(form, df, REML = reml, control = mm_control(verbose = -1))
}

mk_cg2_lmm_simple <- function(seed = 52L, reml = TRUE) {
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

mk_cg2_glmm <- function(seed = 53L) {
  set.seed(seed)
  g <- gl(8, 6)
  x <- rnorm(48)
  eta <- 0.2 * x + rnorm(8, sd = 0.3)[g]
  y <- rbinom(48, 1L, plogis(eta))
  glmm(
    y ~ x + (1 | g),
    data.frame(y = y, x = x, g = g),
    family = binomial(),
    control = mm_control(verbose = -1)
  )
}

mk_cg2_intercept_only <- function(seed = 54L) {
  set.seed(seed)
  g <- gl(10, 4)
  y <- rnorm(40, sd = 0.4) + rnorm(10, sd = 0.5)[g]
  lmm(
    y ~ 1 + (1 | g),
    data.frame(y = y, g = g),
    REML = FALSE,
    control = mm_control(verbose = -1)
  )
}

# ---- summary / print.summary ------------------------------------------------

test_that("print.summary.mm_lmm covers coefficients, none, reason, verbose", {
  fit <- mk_cg2_lmm_simple()
  s <- summary(fit)
  out <- capture.output(print(s))
  expect_true(any(grepl("Linear mixed model", out, fixed = TRUE)))
  expect_true(any(grepl("Fixed effects:", out, fixed = TRUE)))
  expect_true(any(grepl("Inference status:", out, fixed = TRUE)))

  out_v <- capture.output(print(s, verbose = TRUE))
  expect_true(any(grepl("Inference rows are supplied by Rust", out_v, fixed = TRUE)))

  s_none_tests <- summary(fit, tests = "none")
  expect_null(s_none_tests$inference)
  expect_true(length(capture.output(print(s_none_tests))) > 0L)

  s_method_none <- summary(fit, method = "none")
  out_r <- capture.output(print(s_method_none))
  expect_true(any(grepl("Reasons:", out_r, fixed = TRUE)))
  expect_true(any(grepl("inference_not_requested", out_r, fixed = TRUE)))
})

test_that("summary.mm_glmm and print.summary.mm_glmm cover default and none", {
  fit <- mk_cg2_glmm()
  s <- summary(fit)
  expect_s3_class(s, "summary.mm_glmm")
  out <- capture.output(print(s))
  expect_true(any(grepl("Generalized linear mixed model", out, fixed = TRUE)))
  expect_true(any(grepl("Family/link:", out, fixed = TRUE)))
  expect_true(any(grepl("Fixed effects:", out, fixed = TRUE)))

  s_none <- summary(fit, tests = "none")
  expect_null(s_none$inference)
  expect_true(length(capture.output(print(s_none))) > 0L)

  out_v <- capture.output(print(s, verbose = TRUE))
  expect_true(is.character(out_v) && length(out_v) > 0L)
})

# ---- predict guards and se.fit / interval shapes ----------------------------

test_that("predict.mm_lmm rejects bad allow.new.levels and unsupported re.form", {
  fit <- mk_cg2_lmm_simple()
  expect_error(predict(fit, allow.new.levels = NA), class = "mm_arg_error")
  expect_error(predict(fit, allow.new.levels = c(TRUE, FALSE)),
               class = "mm_arg_error")
  expect_error(predict(fit, re.form = ~(1 | g)),
               class = "mm_inference_unavailable")
  expect_error(predict(fit, random.only = TRUE), class = "mm_arg_error")
})

test_that("predict.mm_lmm se.fit and interval shapes (conditional + population)", {
  fit <- mk_cg2_lmm_simple()

  se_only <- predict(fit, se.fit = TRUE)
  expect_true(is.list(se_only) && all(c("fit", "se.fit") %in% names(se_only)))
  expect_equal(length(se_only$fit), nobs(fit))

  ci <- predict(fit, interval = "confidence")
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 3L)
  expect_identical(attr(ci, "interval"), "confidence")

  both <- predict(fit, se.fit = TRUE, interval = "confidence")
  expect_true(is.list(both) && all(c("fit", "se.fit") %in% names(both)))
  expect_true(is.matrix(both$fit) && ncol(both$fit) == 3L)

  pred_int <- predict(fit, interval = "prediction")
  expect_true(is.matrix(pred_int) && ncol(pred_int) == 3L)
  expect_identical(attr(pred_int, "interval"), "prediction")

  pop_se <- predict(fit, re.form = NA, se.fit = TRUE)
  expect_true(is.list(pop_se) && all(c("fit", "se.fit") %in% names(pop_se)))

  pop_ci <- predict(fit, re.form = NA, interval = "confidence")
  expect_true(is.matrix(pop_ci) && ncol(pop_ci) == 3L)

  pop_both <- predict(fit, re.form = NA, se.fit = TRUE, interval = "confidence")
  expect_true(is.list(pop_both) && is.matrix(pop_both$fit))

  pop_pred <- predict(fit, re.form = NA, interval = "prediction")
  expect_true(is.matrix(pop_pred) && ncol(pop_pred) == 3L)

  nd <- fit$model_frame[1:3, , drop = FALSE]
  nd_out <- predict(fit, newdata = nd, re.form = NA, se.fit = TRUE)
  expect_equal(length(nd_out$fit), 3L)
})

test_that("predict.mm_glmm guards and basic response/link paths", {
  fit <- mk_cg2_glmm()
  expect_error(predict(fit, allow.new.levels = NA), class = "mm_arg_error")
  expect_error(predict(fit, re.form = ~(1 | g)),
               class = "mm_inference_unavailable")
  expect_error(
    predict(fit, interval = "prediction", type = "link"),
    class = "mm_inference_unavailable"
  )
  expect_error(
    predict(fit, interval = "prediction", re.form = NA),
    class = "mm_inference_unavailable"
  )
  # Bernoulli (0/1 binomial) prediction intervals are available on response scale.
  pi <- predict(fit, interval = "prediction")
  expect_true(is.matrix(pi) && ncol(pi) == 3L)
  expect_identical(attr(pi, "interval"), "prediction")

  pr <- predict(fit)
  expect_equal(length(pr), nobs(fit))
  pl <- predict(fit, type = "link")
  expect_equal(length(pl), nobs(fit))
  pop <- predict(fit, re.form = NA)
  expect_equal(length(pop), nobs(fit))
  expect_equal(length(fitted(fit)), nobs(fit))
  expect_equal(length(residuals(fit)), nobs(fit))
})

test_that("residuals.mm_lmm pearson/scaled and fitted names", {
  fit <- mk_cg2_lmm_simple()
  r0 <- residuals(fit)
  rp <- residuals(fit, type = "pearson")
  rs <- residuals(fit, scaled = TRUE)
  expect_equal(length(r0), nobs(fit))
  expect_equal(length(rp), nobs(fit))
  expect_false(isTRUE(all.equal(unname(r0), unname(rp))))
  expect_equal(length(rs), nobs(fit))
  expect_equal(length(fitted(fit)), nobs(fit))
})

# ---- mm_grid / mm_means / mm_predictions / mm_comparisons -------------------

test_that("mm_grid happy path and arg errors", {
  fit <- mk_cg2_lmm(factor_x = TRUE)
  gr <- mm_grid(fit, ~f)
  expect_s3_class(gr, "mm_grid")
  expect_true(nrow(gr$grid) >= 2L)
  out <- capture.output(print(gr))
  expect_true(is.character(out) && length(out) > 0L)

  gr_at <- mm_grid(fit, ~f, at = list(f = c("a", "b")))
  expect_equal(nrow(gr_at$grid), 2L)

  expect_error(mm_grid(fit, ~not_a_var), class = "mm_arg_error")
  expect_error(mm_grid(fit, ~f, at = 1:3), class = "mm_arg_error")
  expect_error(mm_grid(fit, ~f, at = list(1, 2)), class = "mm_arg_error")
})

test_that("mm_means and mm_predictions happy + level errors", {
  fit <- mk_cg2_lmm(factor_x = TRUE)
  mn <- mm_means(fit, ~f, method = "asymptotic")
  expect_s3_class(mn, "mm_marginal_quantity")
  expect_true(nrow(mn$table) >= 2L)
  expect_true(length(capture.output(print(mn))) > 0L)

  pr <- mm_predictions(fit, specs = ~f, method = "asymptotic")
  expect_s3_class(pr, "mm_marginal_quantity")
  expect_true(nrow(pr$table) >= 2L)
  expect_true(length(capture.output(print(pr))) > 0L)

  expect_error(mm_means(fit, ~f, level = 0), class = "mm_arg_error")
  expect_error(mm_means(fit, ~f, level = 1), class = "mm_arg_error")
  expect_error(mm_predictions(fit, specs = ~f, level = c(0.9, 0.95)),
               class = "mm_arg_error")
})

test_that("mm_comparisons difference path and unsupported comparison", {
  fit <- mk_cg2_lmm(factor_x = TRUE)
  cmp <- mm_comparisons(fit, ~f, method = "asymptotic")
  expect_s3_class(cmp, "mm_marginal_quantity")
  expect_true(nrow(cmp$table) >= 1L)
  expect_true(length(capture.output(print(cmp))) > 0L)

  expect_error(
    mm_comparisons(fit, ~f, comparison = "ratio"),
    class = "mm_inference_unavailable"
  )
})

# ---- getME / extractors -----------------------------------------------------

test_that("getME common names, multi-name, and unknown name", {
  fit <- mk_cg2_lmm_simple()
  expect_true(is.matrix(getME(fit, "X")))
  expect_true(inherits(getME(fit, "Z"), "Matrix") || is.matrix(getME(fit, "Z")))
  expect_true(inherits(getME(fit, "Zt"), "Matrix") || is.matrix(getME(fit, "Zt")))
  Lam <- getME(fit, "Lambda")
  expect_true(!is.null(Lam))
  expect_true(!is.null(getME(fit, "Lambdat")))
  expect_equal(length(getME(fit, "theta")), length(fit$theta))
  expect_equal(unname(getME(fit, "beta")), unname(fixef(fit)))
  expect_equal(unname(getME(fit, "fixef")), unname(fixef(fit)))
  expect_equal(length(getME(fit, "y")), nobs(fit))
  expect_equal(length(getME(fit, "mu")), nobs(fit))
  expect_true(is.list(getME(fit, "flist")))
  expect_true(is.list(getME(fit, "cnms")))

  parts <- getME(fit, c("theta", "beta", "y"))
  expect_true(is.list(parts) && all(c("theta", "beta", "y") %in% names(parts)))

  expect_error(getME(fit, character()), class = "mm_arg_error")
  expect_error(getME(fit, 1L), class = "mm_arg_error")
  expect_error(getME(fit, "bogus_component"), class = "mm_arg_error")
})

test_that("model.matrix, vcov, ngrps, sigma, logLik, formula, nobs, terms", {
  fit <- mk_cg2_lmm_simple()
  expect_equal(ncol(model.matrix(fit, type = "fixed")), length(fixef(fit)))
  Z <- model.matrix(fit, type = "random")
  expect_true(nrow(Z) == nobs(fit) || ncol(Z) == nobs(fit))

  V <- vcov(fit)
  expect_true(is.matrix(unclass(V)))
  V_corr <- vcov(fit, correlation = TRUE)
  expect_true(!is.null(attr(V_corr, "correlation")))
  Vt <- vcov(fit, type = "theta")
  expect_identical(attr(Vt, "mm_unavailable_reason"), "theta_covariance_unavailable")

  expect_equal(unname(ngrps(fit)), 12L)
  expect_error(ngrps(1), class = "mm_arg_error")
  expect_true(is.finite(sigma(fit)))
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_equal(nobs(fit), 60L)
  expect_true(inherits(formula(fit), "formula"))
  expect_true(inherits(terms(fit), "terms"))
  expect_equal(nrow(model.frame(fit)), nobs(fit))
  expect_true(is.null(weights(fit)) || is.numeric(weights(fit)))
  aic <- extractAIC(fit)
  expect_equal(length(aic), 2L)
})

test_that("AIC/BIC refuse multi-object comparison", {
  fit <- mk_cg2_lmm_simple()
  expect_error(AIC(fit, fit), class = "mm_inference_unavailable")
  expect_error(BIC(fit, fit), class = "mm_inference_unavailable")
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

test_that("as.data.frame.mm_varcorr includes residual and correlation rows", {
  set.seed(55L)
  g <- gl(10, 6)
  x <- rnorm(60)
  y <- 0.4 * x + rnorm(60, sd = 0.35) +
    rnorm(10, sd = 0.5)[g] + rnorm(10, sd = 0.25)[g] * x
  fit <- lmm(
    y ~ x + (1 + x | g),
    data.frame(y = y, x = x, g = g),
    control = mm_control(verbose = -1)
  )
  vc <- as.data.frame(VarCorr(fit))
  expect_true(nrow(vc) >= 3L)
  expect_true("Residual" %in% vc$grp)
  expect_true(any(!is.na(vc$var2)))
})

# ---- changes / model_report / parameterization ------------------------------

test_that("changes() and print.mm_change_log work on fits", {
  fit <- mk_cg2_lmm_simple()
  ch <- changes(fit)
  expect_s3_class(ch, "mm_change_log")
  expect_true(is.data.frame(ch$table))
  out <- capture.output(print(ch))
  expect_true(any(grepl("Model changes:", out, fixed = TRUE)))
  expect_true(any(grepl("Stage-by-stage", out, fixed = TRUE)))

  spec <- compile_model(y ~ x + (1 | g), fit$model_frame)
  ch_spec <- changes(spec)
  expect_s3_class(ch_spec, "mm_change_log")
  expect_true(length(capture.output(print(ch_spec))) > 0L)
})

test_that("model_report and reporting_table section/view paths", {
  fit <- mk_cg2_lmm_simple()
  mr <- model_report(fit)
  expect_s3_class(mr, "mm_model_report")
  expect_true(length(mr$sections) >= 5L)
  out <- capture.output(print(mr))
  expect_true(is.character(out) && length(out) > 0L)

  mr_sub <- model_report(fit, sections = c("overview", "fixed_effects"))
  expect_true(all(c("overview", "fixed_effects") %in% names(mr_sub$sections)))

  rt <- reporting_table(fit, section = "overview", view = "audit")
  expect_s3_class(rt, "mm_reporting_table")
  expect_true(length(capture.output(print(rt))) > 0L)

  rt2 <- reporting_table(mr, section = "fit_statistics")
  expect_s3_class(rt2, "mm_reporting_table")

  expect_error(model_report(fit, sections = 1L), class = "mm_schema_error")
  expect_error(
    reporting_table(fit, section = c("overview", "fixed_effects")),
    class = "mm_schema_error"
  )
  expect_error(
    model_report(fit, sections = "not_a_section"),
    class = "mm_schema_error"
  )
})

test_that("parameterization prints theta map", {
  fit <- mk_cg2_lmm_simple()
  pm <- parameterization(fit)
  expect_s3_class(pm, "mm_theta_map")
  expect_true(nrow(pm$table) >= 1L)
  out <- capture.output(print(pm))
  expect_true(any(grepl("Covariance parameterization:", out, fixed = TRUE)))
})

# ---- broom ------------------------------------------------------------------

test_that("broom tidy/glance/augment work via package methods", {
  skip_if_not_installed("generics")
  fit <- mk_cg2_lmm_simple()
  gfit <- mk_cg2_glmm()

  tidy_fn <- getExportedValue("generics", "tidy")
  glance_fn <- getExportedValue("generics", "glance")
  augment_fn <- getExportedValue("generics", "augment")

  t1 <- tidy_fn(fit, effects = c("fixed", "ran_pars", "ran_vals"),
                conf.int = TRUE)
  expect_true(nrow(t1) >= 3L)
  expect_true(all(c("fixed", "ran_pars", "ran_vals") %in% t1$effect))
  expect_true(all(c("conf.low", "conf.high") %in% names(t1)))

  g1 <- glance_fn(fit)
  expect_equal(nrow(g1), 1L)
  expect_true(all(c("nobs", "AIC", "BIC", "logLik") %in% names(g1)))

  a1 <- augment_fn(fit)
  expect_equal(nrow(a1), nobs(fit))
  expect_true(all(c(".fitted", ".resid") %in% names(a1)))

  tg <- tidy_fn(gfit, effects = "fixed", conf.int = TRUE)
  expect_equal(nrow(tg), length(fixef(gfit)))
  expect_equal(nrow(glance_fn(gfit)), 1L)
  expect_equal(nrow(augment_fn(gfit)), nobs(gfit))
})

test_that("broom methods are reachable without generics installed", {
  fit <- mk_cg2_lmm_simple()
  t1 <- mixeff:::tidy.mm_lmm(fit, effects = c("fixed", "ran_pars"),
                             conf.int = TRUE)
  expect_true(nrow(t1) >= 2L)
  expect_true("conf.low" %in% names(t1))
  t2 <- mixeff:::tidy.mm_lmm(fit, effects = "ran_vals")
  expect_true(nrow(t2) >= 1L)
  expect_equal(nrow(mixeff:::glance.mm_lmm(fit)), 1L)
  expect_equal(nrow(mixeff:::augment.mm_lmm(fit)), nobs(fit))

  gfit <- mk_cg2_glmm()
  tg <- mixeff:::tidy.mm_glmm(gfit, effects = "fixed", conf.int = TRUE)
  expect_equal(nrow(tg), length(fixef(gfit)))
  expect_equal(nrow(mixeff:::glance.mm_glmm(gfit)), 1L)
  expect_equal(nrow(mixeff:::augment.mm_glmm(gfit)), nobs(gfit))
})

# ---- random_options / compare_covariance ------------------------------------

test_that("random_options and compare_covariance print nearby map", {
  fit <- mk_cg2_lmm_simple()
  spec <- compile_model(y ~ x + (1 | g), fit$model_frame)
  opts <- random_options(spec, "g", "x")
  expect_s3_class(opts, "mm_random_options")
  expect_true(nrow(opts$options) >= 1L)
  out <- capture.output(print(opts))
  expect_true(any(grepl("Random-effect options", out, fixed = TRUE)))

  expect_error(random_options(spec, ""), class = "mm_arg_error")
  expect_error(random_options(spec, "not_a_group"), class = "mm_schema_error")

  cmp <- compare_covariance(spec)
  expect_s3_class(cmp, "mm_compare_covariance")
  expect_true(nrow(cmp$table) >= 1L)
  out_c <- capture.output(print(cmp))
  expect_true(any(grepl("Covariance comparison:", out_c, fixed = TRUE)))
})

# ---- inference: estimability, df, confint, profile, random effect -----------

test_that("estimability() and print cover default and L paths", {
  fit <- mk_cg2_lmm_simple()
  est <- estimability(fit)
  expect_s3_class(est, "mm_estimability")
  expect_true(nrow(est$table) >= 2L)
  out <- capture.output(print(est))
  expect_true(any(grepl("Estimability:", out, fixed = TRUE)))

  est_L <- estimability(fit, L = c(0, 1))
  expect_equal(nrow(est_L$table), 1L)
  expect_true(length(capture.output(print(est_L))) > 0L)
})

test_that("df_for_contrast() prints none and asymptotic paths", {
  fit <- mk_cg2_lmm_simple()
  df_none <- df_for_contrast(fit, c(0, 1), method = "none")
  expect_s3_class(df_none, "mm_df_for_contrast")
  out <- capture.output(print(df_none))
  expect_true(any(grepl("Degrees of freedom", out, fixed = TRUE)))

  df_as <- df_for_contrast(fit, c(0, 1), method = "asymptotic")
  expect_s3_class(df_as, "mm_df_for_contrast")
  expect_true(length(capture.output(print(df_as))) > 0L)

  df_auto <- df_for_contrast(fit, diag(2), method = "auto")
  expect_equal(nrow(df_auto$table), 2L)
})

test_that("confint parm numeric index and bootstrap/profile edges", {
  fit <- mk_cg2_lmm_simple(reml = TRUE)
  ci1 <- confint(fit, parm = 1L, method = "asymptotic")
  expect_s3_class(ci1, "mm_confint")
  expect_equal(nrow(ci1), 1L)
  expect_true(length(capture.output(print(ci1))) > 0L)

  ci2 <- confint(fit, parm = 2L, method = "wald")
  expect_equal(nrow(ci2), 1L)

  ci_b <- confint(
    fit, parm = 1L, method = "bootstrap",
    bootstrap = bootstrap_control(nsim = 2L, seed = 1L)
  )
  expect_true(is.matrix(ci_b) && nrow(ci_b) == 1L)

  ci_p <- confint(fit, method = "profile")
  expect_true(is.matrix(ci_p) || inherits(ci_p, "mm_confint"))
  expect_true(length(capture.output(print(ci_p))) > 0L)
})

test_that("profile() print and confint.mm_profile edges", {
  fit <- mk_cg2_lmm_simple(reml = TRUE)
  pr <- profile(fit)
  expect_s3_class(pr, "mm_profile")
  out <- capture.output(print(pr))
  expect_true(any(grepl("Profile-likelihood", out, fixed = TRUE)))

  mat <- confint(pr)
  expect_true(is.matrix(mat) && ncol(mat) == 2L)
  mat2 <- confint(pr, parm = "sigma")
  expect_equal(nrow(mat2), 1L)
  expect_error(confint(pr, level = 0.5), class = "mm_arg_error")
})

test_that("test_random_effect() guards, ML path, and print", {
  fit_ml <- mk_cg2_lmm_simple(reml = FALSE)
  out <- test_random_effect(fit_ml, "g")
  expect_s3_class(out, "mm_random_effect_test")
  expect_true(nrow(out$table) >= 1L)
  printed <- capture.output(print(out))
  expect_true(any(grepl("Random-effect variance-component test", printed,
                        fixed = TRUE)))

  expect_error(test_random_effect(fit_ml, ""), class = "mm_arg_error")
  expect_error(test_random_effect(fit_ml, character()), class = "mm_arg_error")
  expect_error(test_random_effect(fit_ml, 1L), class = "mm_arg_error")

  fit_reml <- mk_cg2_lmm_simple(reml = TRUE)
  expect_error(
    test_random_effect(fit_reml, "g", refit_for_comparison = "error"),
    class = "mm_inference_unavailable"
  )
  # auto refit from REML -> ML
  out_auto <- test_random_effect(fit_reml, "g", refit_for_comparison = "auto")
  expect_s3_class(out_auto, "mm_random_effect_test")
})

test_that("mm_lincomb and contrast asymptotic print", {
  fit <- mk_cg2_lmm_simple()
  w <- setNames(c(0, 1), names(fixef(fit)))
  lc <- mm_lincomb(fit, w)
  expect_true(!is.null(lc$estimate) || is.numeric(lc) || is.list(lc))
  ct <- contrast(fit, c(0, 1), method = "asymptotic")
  expect_true(length(capture.output(print(ct))) > 0L)
})

# ---- compare / drop1 / anova ------------------------------------------------

test_that("drop1 on intercept-only yields empty table", {
  fit0 <- mk_cg2_intercept_only()
  d1 <- drop1(fit0, test = "Chisq")
  expect_s3_class(d1, "mm_drop1")
  expect_equal(nrow(d1$table), 0L)
  out <- capture.output(print(d1))
  expect_true(any(grepl("Single-term deletion", out, fixed = TRUE)))
})

test_that("anova types I/II/III/block and method none", {
  fit <- mk_cg2_lmm()
  a1 <- anova(fit, type = "I", method = "asymptotic")
  expect_s3_class(a1, "mm_anova")
  expect_true(nrow(a1$table) >= 1L)

  a2 <- anova(fit, type = "II", method = "asymptotic")
  expect_true(nrow(a2$table) >= 1L)

  a3 <- anova(fit, type = "III", method = "asymptotic")
  expect_true(nrow(a3$table) >= 1L)

  ab <- anova(fit, type = "block", method = "asymptotic")
  expect_true(nrow(ab$table) >= 1L)

  an <- anova(fit, type = "III", method = "none")
  expect_true(nrow(an$table) >= 1L)
  expect_true(length(capture.output(print(a3))) > 0L)
})

test_that("drop1 with explicit scope and test none", {
  fit <- mk_cg2_lmm(reml = FALSE)
  d1 <- drop1(fit, scope = "x", test = "none")
  expect_s3_class(d1, "mm_drop1")
  expect_true(nrow(d1$table) >= 1L)
  expect_true("x" %in% d1$table$dropped)
})

# ---- revive / lazy cache / inference_table print ----------------------------

test_that("revive roundtrip clears handle and rebuilds lazy cache", {
  fit <- mk_cg2_lmm_simple()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  restored <- revive(readRDS(tmp))
  expect_s3_class(restored, "mm_fit")
  expect_false(fit_handle_alive(restored))
  expect_equal(dim(getME(restored, "X")), dim(getME(fit, "X")))
  expect_equal(unname(fixef(restored)), unname(fixef(fit)))
  expect_true(length(capture.output(print(inference_table(restored)))) > 0L)
})

test_that("fit_handle_alive.default is FALSE", {
  expect_false(fit_handle_alive(1))
  expect_false(fit_handle_alive(list()))
})

# ---- json parse remaining paths ---------------------------------------------

test_that("mm_json_parse_* cover missing schema and empty inference table", {
  expect_error(
    mixeff:::mm_json_parse_artifact("{}"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_audit_report("{"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_lmm_fit(""),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_lmm_fit("{"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_lmm_fit("{}"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_fit_summary(NULL),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_fit_summary(list()),
    class = "mm_schema_error"
  )

  empty <- mixeff:::mm_json_parse_fixed_effect_inference_table(list(
    schema_name = "mixedmodels.fixed_effect_inference_table",
    schema_version = "1.1.0",
    rows = list()
  ))
  expect_equal(nrow(empty$table), 0L)
  expect_true(inherits(empty, "mm_fixed_effect_inference_table"))
})

# ---- fit-lmm / mm_control validation ----------------------------------------

test_that("lmm rejects REML=NA and bad control", {
  df <- data.frame(y = rnorm(20), x = rnorm(20), g = gl(5, 4))
  expect_error(
    lmm(y ~ x + (1 | g), df, REML = NA, control = mm_control(verbose = -1)),
    class = "mm_fit_error"
  )
  expect_error(
    lmm(y ~ x + (1 | g), df, control = "bad"),
    class = "mm_arg_error"
  )
})

test_that("mm_control validates verbose and max_feval", {
  expect_error(mm_control(verbose = NA), class = "mm_arg_error")
  expect_error(mm_control(verbose = c(0, 1)), class = "mm_arg_error")
  expect_error(mm_control(max_feval = 0), class = "mm_arg_error")
  expect_error(mm_control(max_feval = -1), class = "mm_arg_error")
  expect_error(mm_control(optimizer = 1L), class = "mm_arg_error")
  expect_error(mm_control(start = NA_real_), class = "mm_arg_error")
  expect_error(mm_control(ftol_rel = 0), class = "mm_arg_error")
  expect_error(mm_control(ftol_abs = -1), class = "mm_arg_error")
  expect_error(mm_control(xtol_rel = NA_real_), class = "mm_arg_error")
  ok <- mm_control(verbose = -1, max_feval = 10L, optimizer = "bobyqa",
                   ftol_rel = 1e-6)
  expect_s3_class(ok, "mm_control")
})

# ---- simulate / diagnostics / confint glmm ----------------------------------

test_that("simulate population path and re.form guard", {
  fit <- mk_cg2_lmm_simple()
  sims <- simulate(fit, nsim = 2L, seed = 1L, re.form = NA)
  expect_equal(ncol(sims), 2L)
  expect_equal(nrow(sims), nobs(fit))
  expect_error(simulate(fit, re.form = ~(1 | g)),
               class = "mm_inference_unavailable")
})

test_that("diagnostics filters and fit_status on fit and spec", {
  fit <- mk_cg2_lmm_simple()
  d_all <- diagnostics(fit)
  expect_s3_class(d_all, "mm_diagnostics")
  d_empty <- diagnostics(fit, stage = "does_not_exist_stage")
  expect_equal(nrow(d_empty$table), 0L)
  out <- capture.output(print(d_empty))
  expect_true(any(grepl("none", out, fixed = TRUE)))

  expect_true(is.character(fit_status(fit)) && length(fit_status(fit)) == 1L)
  spec <- compile_model(y ~ x + (1 | g), fit$model_frame)
  expect_true(is.character(fit_status(spec)))
})

test_that("confint.mm_glmm refuses profile", {
  fit <- mk_cg2_glmm()
  expect_error(confint(fit, method = "profile"),
               class = "mm_inference_unavailable")
})

# ---- GLMM extractors and VarCorr print notes --------------------------------

test_that("GLMM extractors and VarCorr/coef print", {
  fit <- mk_cg2_glmm()
  expect_equal(length(fixef(fit)), 2L)
  expect_true(is.list(ranef(fit)))
  expect_s3_class(VarCorr(fit), "mm_varcorr")
  expect_equal(nobs(fit), 48L)
  expect_true(inherits(formula(fit), "formula"))
  expect_true(inherits(terms(fit), "terms"))
  expect_equal(unname(ngrps(fit)), 8L)
  expect_true(is.matrix(unclass(vcov(fit))))
  expect_true(length(capture.output(print(VarCorr(fit)))) > 0L)
  expect_true(length(capture.output(print(coef(fit)))) > 0L)
})

# ---- reporting_table comparison objects -------------------------------------

test_that("reporting_table on drop1 and compare objects", {
  nest_df <- {
    set.seed(56L)
    g <- gl(12, 5)
    x <- rnorm(60)
    z <- rnorm(60)
    y <- 0.4 * x + 0.3 * z + rnorm(60, sd = 0.35) + rnorm(12, sd = 0.5)[g]
    data.frame(y = y, x = x, z = z, g = g)
  }
  full <- lmm(y ~ x + z + (1 | g), nest_df, REML = FALSE,
              control = mm_control(verbose = -1))
  red <- lmm(y ~ x + (1 | g), nest_df, REML = FALSE,
             control = mm_control(verbose = -1))
  cmp <- compare(red, full)
  rt <- reporting_table(cmp)
  expect_s3_class(rt, "mm_reporting_table")
  expect_true(length(capture.output(print(rt))) > 0L)

  d1 <- drop1(full, test = "Chisq")
  rt_d1 <- reporting_table(d1)
  expect_s3_class(rt_d1, "mm_reporting_table")
})

# ---- is_singular / optimizer_certificate / random_blocks --------------------

test_that("is_singular and certificate helpers on simple LMM", {
  fit <- mk_cg2_lmm_simple()
  expect_type(is_singular(fit), "logical")
  expect_error(is_singular(1), class = "mm_arg_error")
  cert <- optimizer_certificate(fit)
  expect_true(length(capture.output(print(cert))) > 0L)
  blocks <- random_blocks(fit)
  expect_true(nrow(blocks$table) >= 1L)
  expect_true(length(capture.output(print(blocks))) > 0L)
})
