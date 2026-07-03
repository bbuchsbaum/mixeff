# Tests for prediction standard errors and intervals: population
# (fixed-effect Wald, R-side) and conditional (engine prediction-variance
# payload including the random-effect contribution).

mm_se_data <- function() {
  set.seed(321)
  ng <- 14L
  per <- 9L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
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

test_that("conditional se.fit comes from the engine prediction variance", {
  df <- mm_se_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  out <- predict(fit, se.fit = TRUE)  # default re.form = NULL (conditional)
  expect_true(is.list(out) && all(c("fit", "se.fit") %in% names(out)))
  expect_true(all(is.finite(out$se.fit) & out$se.fit > 0))
  # the engine payload's fixed component must reproduce the population Wald SE
  pop <- predict(fit, re.form = NA, se.fit = TRUE)
  pv <- mm_lmm_prediction_variance(fit, fit$model_frame, FALSE, 0.95)
  expect_equal(sqrt(pv$fixed_variance), unname(pop$se.fit), tolerance = 1e-6)
  expect_true(all(pv$status == "available"))
})

test_that("conditional intervals bracket fit and prediction is wider", {
  df <- mm_se_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  ci <- predict(fit, interval = "confidence", level = 0.95)
  pi <- predict(fit, interval = "prediction", level = 0.95)
  expect_true(is.matrix(ci) && all(is.finite(ci)))
  expect_equal(colnames(ci), c("fit", "lwr", "upr"))
  expect_true(all(ci[, "lwr"] < ci[, "fit"] & ci[, "fit"] < ci[, "upr"]))
  se <- predict(fit, se.fit = TRUE)$se.fit
  expect_equal(unname(ci[, "upr"] - ci[, "fit"]),
               unname(qnorm(0.975) * as.numeric(se)), tolerance = 1e-6)
  expect_true(all((pi[, "upr"] - pi[, "lwr"]) > (ci[, "upr"] - ci[, "lwr"])))
})

test_that("conditional se.fit on unseen levels is NA with an engine reason", {
  df <- mm_se_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  nd <- data.frame(x = c(0, 1), g = factor(c("999", "998")))
  out <- predict(fit, newdata = nd, re.form = NULL, se.fit = TRUE,
                 allow.new.levels = TRUE)
  expect_true(all(is.finite(out$fit)))
  expect_true(all(is.na(out$se.fit)))
  reasons <- attr(out$se.fit, "mm_reason")
  expect_true(is.character(reasons) && all(!is.na(reasons)))
})
