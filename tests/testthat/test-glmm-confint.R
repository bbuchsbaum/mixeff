# Tests for confint.mm_glmm (asymptotic Wald intervals).

mm_confint_gamma_data <- function() {
  y <- x <- numeric()
  g <- character()
  for (group in 0:3) {
    for (obs in 0:4) {
      xv <- obs - 2
      eta <- 0.5 + 0.2 * xv + (group - 1.5) * 0.08
      y <- c(y, exp(eta) * (0.95 + 0.02 * ((group + obs) %% 3)))
      x <- c(x, xv)
      g <- c(g, paste0("g", group + 1L))
    }
  }
  data.frame(y = y, x = x, g = g)
}

test_that("confint.mm_glmm returns Wald intervals bracketing the estimates", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_gamma_data(),
              family = Gamma(link = "log"),
              method = "joint_laplace", control = mm_control(verbose = -1))
  ci <- confint(fit)
  expect_s3_class(ci, "mm_confint")
  expect_equal(rownames(ci), names(fixef(fit)))
  est <- fixef(fit)
  expect_true(all(ci[, 1] < est & est < ci[, 2]))
  expect_identical(attr(ci, "method"),
                   "wald_asymptotic_from_rust_inference_table")
  expect_identical(attr(ci, "status"), "available")
})

test_that("confint.mm_glmm matches estimate +/- z*SE and respects level/parm", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_gamma_data(),
              family = Gamma(link = "log"),
              method = "joint_laplace", control = mm_control(verbose = -1))
  ci <- confint(fit, parm = "x", level = 0.90)
  z <- stats::qnorm(0.95)
  est <- unname(fixef(fit)["x"])
  se <- summary(fit, tests = "coefficients")$coefficients["x", "Std. Error"]
  expect_equal(nrow(ci), 1L)
  expect_equal(unname(ci[1, 1]), est - z * se, tolerance = 1e-10)
  expect_equal(unname(ci[1, 2]), est + z * se, tolerance = 1e-10)
})

test_that("confint.mm_glmm refuses profile/bootstrap with a typed reason", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_gamma_data(),
              family = Gamma(link = "log"),
              method = "joint_laplace", control = mm_control(verbose = -1))
  expect_error(confint(fit, method = "profile"),
               class = "mm_inference_unavailable")
  expect_error(confint(fit, method = "bootstrap"),
               class = "mm_inference_unavailable")
})

test_that("confint.mm_glmm accepts the asymptotic synonym", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_gamma_data(),
              family = Gamma(link = "log"),
              method = "joint_laplace", control = mm_control(verbose = -1))
  a <- confint(fit, method = "wald")
  b <- confint(fit, method = "asymptotic")
  expect_equal(unclass(a), unclass(b))
})

test_that("confint.mm_glmm refuses unsupported fast-PIRLS Wald rows", {
  fit <- glmm(y ~ x + (1 | g), mm_confint_gamma_data(),
              family = Gamma(link = "log"),
              method = "pirls_profiled", control = mm_control(verbose = -1))
  expect_error(confint(fit), class = "mm_inference_unavailable")
})
