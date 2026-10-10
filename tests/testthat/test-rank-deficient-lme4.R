# Rank-deficient fixed effects (mixeff-rs 2873312): the engine keeps the
# earlier of two collinear columns (R's qr() rule) and lists the dropped
# ones; fixef()/vcov()/summary() follow lme4's shapes.

mm_rank_def_data <- function() {
  set.seed(3)
  d <- data.frame(g = factor(rep(1:15, each = 10)), x = rnorm(150))
  d$z <- 2 * d$x
  d$w <- rnorm(150)
  d$y <- 1 + d$x + rnorm(15)[d$g] + rnorm(150)
  d
}

test_that("collinear columns: fixef/vcov/summary match lme4", {
  skip_if_not_installed("lme4")
  d <- mm_rank_def_data()
  fit <- lmm(y ~ x + z + w + (1 | g), d, control = mm_control(verbose = -1))
  ref <- suppressMessages(lme4::lmer(y ~ x + z + w + (1 | g), d))

  expect_identical(mm_aliased_coefficients(fit), "z")
  expect_identical(names(fixef(fit)), names(lme4::fixef(ref)))
  expect_equal(fixef(fit), lme4::fixef(ref), tolerance = 1e-5)
  full <- fixef(fit, add.dropped = TRUE)
  expect_equal(full, lme4::fixef(ref, add.dropped = TRUE), tolerance = 1e-5)
  expect_true(is.na(full[["z"]]))
  # the stored coefficient vector keeps the dropped column as 0
  expect_identical(unname(fit$beta[["z"]]), 0)

  V <- vcov(fit)
  expect_identical(dimnames(V), dimnames(as.matrix(stats::vcov(ref))))
  expect_equal(unclass(V), as.matrix(stats::vcov(ref)), tolerance = 1e-3,
               ignore_attr = TRUE)

  cf <- summary(fit)$coefficients
  ref_cf <- stats::coef(summary(ref))
  expect_identical(rownames(cf), rownames(ref_cf))
  expect_equal(cf[, c("Estimate", "Std. Error", "t value")],
               ref_cf[, c("Estimate", "Std. Error", "t value")],
               tolerance = 1e-3)
  expect_identical(attr(cf, "mm_aliased"), "z")
  # standard errors stay aligned with their coefficients
  expect_true(is.na(fit$std_errors[["z"]]))
  expect_equal(fit$std_errors[c("(Intercept)", "x", "w")],
               sqrt(diag(as.matrix(stats::vcov(ref)))), tolerance = 1e-3,
               ignore_attr = TRUE)
})

test_that("predictions are unaffected by the dropped column", {
  d <- mm_rank_def_data()
  fit <- lmm(y ~ x + z + w + (1 | g), d, control = mm_control(verbose = -1))
  fit2 <- lmm(y ~ x + w + (1 | g), d, control = mm_control(verbose = -1))
  expect_equal(fitted(fit), fitted(fit2), tolerance = 1e-6)
  expect_equal(fixef(fit), fixef(fit2), tolerance = 1e-6)
})
