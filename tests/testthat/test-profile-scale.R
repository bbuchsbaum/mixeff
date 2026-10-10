# Profile likelihood must scale with the sparse factor, not with n^2.
# Before mixeff-rs#12 the fixed-effect profiles built and factored a dense
# n x n marginal covariance on every evaluation: at n = 4000 that is a
# 128 MB matrix and an O(n^3) Cholesky per step (minutes), and at
# n = 60000 it aborted R trying to allocate 28.8 GB.

test_that("confint(method = 'profile') on a crossed LMM scales to moderate n", {
  skip_on_cran()
  set.seed(42)
  ns <- 80; ni <- 50
  d <- expand.grid(s = factor(seq_len(ns)), i = factor(seq_len(ni)))
  d$x <- rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + rnorm(ns, sd = 0.8)[d$s] +
    rnorm(ni, sd = 0.5)[d$i] + rnorm(nrow(d))
  expect_identical(nrow(d), 4000L)
  fit <- lmm(y ~ x + (1 | s) + (1 | i), d, REML = FALSE,
             control = mm_control(verbose = -1))
  elapsed <- system.time(
    ci <- confint(fit, method = "profile")
  )[["elapsed"]]
  expect_lt(elapsed, 60)
  expect_identical(rownames(ci),
                   c(".sig01", ".sig02", ".sigma", "(Intercept)", "x"))
  expect_true(all(is.finite(unclass(ci)[, 1:2])))
  # every interval brackets its estimate
  est <- c(NA, NA, sigma(fit), fixef(fit))
  ok <- !is.na(est)
  expect_true(all(ci[ok, 1] < est[ok] & est[ok] < ci[ok, 2]))
  skip_if_not_installed("lme4")
  ref <- suppressMessages(confint(
    lme4::lmer(y ~ x + (1 | s) + (1 | i), d, REML = FALSE),
    method = "profile"
  ))
  expect_equal(unclass(ci)[, 1:2], ref, tolerance = 1e-3,
               ignore_attr = TRUE)
})
