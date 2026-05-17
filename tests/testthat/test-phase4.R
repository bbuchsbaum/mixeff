mk_phase4_fit <- function(seed = 40L, slope = TRUE) {
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
  lmm(formula, df, control = mm_control(verbose = -1))
}

test_that("glmm() validates family metadata and reports explicit unavailable bridge", {
  df <- data.frame(
    y = c(0, 1, 0, 1, 1, 0),
    x = c(0, 0, 1, 1, 0, 1),
    subject = factor(rep(1:3, each = 2))
  )
  err <- tryCatch(
    glmm(y ~ x + (1 | subject), df, family = binomial(),
         control = mm_control(verbose = -1)),
    mm_fit_error = function(cnd) cnd
  )
  expect_s3_class(err, "mm_fit_error")
  expect_match(conditionMessage(err), "GLMM fitting is not available", fixed = TRUE)
  expect_equal(err$metadata$family$family, "binomial")
  expect_equal(err$metadata$family$link, "logit")
  expect_error(
    glmm(y ~ x + (1 | subject), df, family = gaussian(),
         control = mm_control(verbose = -1)),
    class = "mm_fit_error"
  )
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
  expect_true(all(!cmp$table$REML))
  expect_true(any(cmp$table$refit))
  expect_true(is.finite(cmp$table$LRT[[2L]]))
  expect_s3_class(av, "mm_model_comparison")
  expect_error(
    compare(reduced, full, refit_for_comparison = "error"),
    class = "mm_inference_unavailable"
  )
})

test_that("drop1() preserves random effects and reports deletion rows", {
  fit <- mk_phase4_fit(slope = TRUE)
  d <- stats::drop1(fit, test = "Chisq")

  expect_s3_class(d, "mm_drop1")
  expect_true(all(c("dropped", "formula", "LRT", "p_value") %in% names(d$table)))
  expect_true("x" %in% d$table$dropped)
  expect_true("z" %in% d$table$dropped)
  expect_true(all(grepl("(1 | subject)", d$table$formula, fixed = TRUE)))
})

test_that("parametric bootstrap comparison runs on a tiny nsim", {
  full <- mk_phase4_fit(slope = TRUE)
  reduced <- mk_phase4_fit(slope = FALSE)

  boot <- parametric_bootstrap(reduced, full, nsim = 2, seed = 101)
  cmp <- compare(reduced, full, method = "bootstrap", nsim = 2, seed = 101)

  expect_s3_class(boot, "mm_parametric_bootstrap")
  expect_equal(length(boot$simulated), 2L)
  expect_true(is.finite(boot$p_value))
  expect_s3_class(cmp$bootstrap, "mm_parametric_bootstrap")
  expect_equal(cmp$table$status[[2L]], "parametric_bootstrap")
})

test_that("manifest advertises shipped simulation and inference surfaces", {
  cap <- mm_formula_manifest()$capabilities
  expect_true(cap$simulate)
  expect_true(cap$inference)
  expect_false(cap$fit_glmm)
})
