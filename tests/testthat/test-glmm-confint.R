# Tests for confint.mm_glmm (asymptotic Wald intervals).

mm_confint_binom_data <- function() {
  set.seed(505)
  ng <- 15L; per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.7)[as.integer(g)]
  y <- rbinom(n, 1, plogis(-0.2 + 0.7 * x + re))
  data.frame(y = y, x = x, g = g)
}

test_that("confint.mm_glmm returns Wald intervals bracketing the estimates", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  ci <- confint(fit)
  expect_s3_class(ci, "mm_confint")
  expect_equal(rownames(ci), names(fixef(fit)))
  est <- fixef(fit)
  expect_true(all(ci[, 1] < est & est < ci[, 2]))
  expect_identical(attr(ci, "method"),
                   "wald_asymptotic_from_stored_standard_errors")
})

test_that("confint.mm_glmm matches estimate +/- z*SE and respects level/parm", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  ci <- confint(fit, parm = "x", level = 0.90)
  z <- stats::qnorm(0.95)
  est <- unname(fixef(fit)["x"])
  se <- unname(fit$std_errors["x"])
  expect_equal(nrow(ci), 1L)
  expect_equal(unname(ci[1, 1]), est - z * se, tolerance = 1e-10)
  expect_equal(unname(ci[1, 2]), est + z * se, tolerance = 1e-10)
})

test_that("confint.mm_glmm refuses profile/bootstrap with a typed reason", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  expect_error(confint(fit, method = "profile"),
               class = "mm_inference_unavailable")
  expect_error(confint(fit, method = "bootstrap"),
               class = "mm_inference_unavailable")
})

test_that("confint.mm_glmm accepts the asymptotic synonym", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  a <- confint(fit, method = "wald")
  b <- confint(fit, method = "asymptotic")
  expect_equal(unclass(a), unclass(b))
})
