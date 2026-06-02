# Tests for population (fixed-effect) prediction standard errors and intervals.

mm_se_data <- function() {
  set.seed(321)
  ng <- 14L; per <- 9L
  g <- factor(rep(seq_len(ng), each = per)); n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 1.0)[as.integer(g)]
  y <- 0.7 + 0.4 * x + re + rnorm(n, sd = 0.8)
  data.frame(y = y, x = x, g = g)
}

test_that("population se.fit equals sqrt(diag(X V X'))", {
  df <- mm_se_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  out <- predict(fit, re.form = NA, se.fit = TRUE)
  expect_true(is.list(out) && all(c("fit", "se.fit") %in% names(out)))
  # manual reference (numeric-only fixed part: engine basis == model.matrix)
  X <- model.matrix(~ x, df)
  V <- as.matrix(unclass(vcov(fit)))
  se_manual <- sqrt(diag(X %*% V %*% t(X)))
  expect_equal(unname(out$se.fit), unname(se_manual), tolerance = 1e-8)
  expect_true(all(out$se.fit > 0))
})

test_that("population confidence interval is fit +/- z*se", {
  df <- mm_se_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  ci <- predict(fit, re.form = NA, interval = "confidence", level = 0.95)
  expect_true(is.matrix(ci))
  expect_equal(colnames(ci), c("fit", "lwr", "upr"))
  se <- predict(fit, re.form = NA, se.fit = TRUE)$se.fit
  z <- qnorm(0.975)
  expect_equal(unname(ci[, "upr"] - ci[, "fit"]), unname(z * se), tolerance = 1e-8)
  expect_true(all(ci[, "lwr"] < ci[, "fit"] & ci[, "fit"] < ci[, "upr"]))
})

test_that("prediction interval is wider than the confidence interval", {
  df <- mm_se_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  ci <- predict(fit, re.form = NA, interval = "confidence")
  pi <- predict(fit, re.form = NA, interval = "prediction")
  expect_true(all((pi[, "upr"] - pi[, "lwr"]) > (ci[, "upr"] - ci[, "lwr"])))
})

test_that("population se.fit works on newdata", {
  df <- mm_se_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  nd <- data.frame(x = c(-1, 0, 1), g = factor(c(1, 2, 3)))
  out <- predict(fit, newdata = nd, re.form = NA, se.fit = TRUE)
  expect_length(out$se.fit, 3L)
  expect_true(all(is.finite(out$se.fit) & out$se.fit > 0))
})

test_that("conditional SE/intervals are refused or NA", {
  df <- mm_se_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  # conditional (default re.form = NULL) interval is refused
  expect_error(predict(fit, interval = "confidence"),
               class = "mm_inference_unavailable")
  # conditional se.fit is NA with a reason
  out <- predict(fit, se.fit = TRUE)
  expect_true(all(is.na(out$se.fit)))
  expect_identical(attr(out$se.fit, "mm_unavailable_reason"),
                   "conditional_prediction_se_unavailable")
})
