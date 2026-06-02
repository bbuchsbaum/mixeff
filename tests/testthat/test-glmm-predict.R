# Tests for predict.mm_glmm (R-side plug-in predictions).

mm_pred_binom_data <- function() {
  set.seed(606)
  ng <- 18L; per <- 14L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.8)[as.integer(g)]
  y <- rbinom(n, 1, plogis(-0.3 + 0.8 * x + re))
  data.frame(y = y, x = x, g = g)
}

test_that("predict.mm_glmm in-sample response matches engine fitted values", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  p <- predict(fit, type = "response")
  expect_equal(unname(p), unname(fitted(fit)), tolerance = 1e-10)
  expect_true(all(p > 0 & p < 1))
})

test_that("predict.mm_glmm link/response are consistent through the link", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  eta <- predict(fit, type = "link")
  mu <- predict(fit, type = "response")
  expect_equal(unname(mu), plogis(unname(eta)), tolerance = 1e-10)
})

test_that("predict.mm_glmm population path is linkinv(X beta)", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  eta_pop <- predict(fit, re.form = NA, type = "link")
  # population linear predictor is intercept + slope * x (no RE)
  beta <- fixef(fit)
  manual <- beta[["(Intercept)"]] + beta[["x"]] * fit$model_frame$x
  expect_equal(unname(eta_pop), unname(manual), tolerance = 1e-8)
  mu_pop <- predict(fit, re.form = NA, type = "response")
  expect_equal(unname(mu_pop), plogis(unname(manual)), tolerance = 1e-10)
})

test_that("predict.mm_glmm conditional newdata reproduces in-sample fitted", {
  df <- mm_pred_binom_data()
  fit <- glmm(y ~ x + (1 + x | g), df, family = binomial(),
              control = mm_control(verbose = -1))
  # Predicting the training rows (shuffled) with re.form = NULL must reproduce
  # the engine's conditional fitted values -> validates the X beta + Z b R-side
  # reconstruction against the authoritative engine output.
  ord <- sample(nrow(df))
  nd <- df[ord, , drop = FALSE]
  p <- predict(fit, newdata = nd, re.form = NULL, type = "response")
  expect_equal(unname(p), unname(fitted(fit))[ord], tolerance = 1e-6)
})

test_that("predict.mm_glmm respects allow.new.levels for unseen groups", {
  df <- mm_pred_binom_data()
  fit <- glmm(y ~ x + (1 | g), df, family = binomial(),
              control = mm_control(verbose = -1))
  nd <- data.frame(x = c(0, 1), g = factor(c("999", "998")))
  expect_error(predict(fit, newdata = nd, re.form = NULL),
               class = "mm_inference_unavailable")
  # With allow.new.levels, unseen groups fall back to the population mean.
  p <- predict(fit, newdata = nd, re.form = NULL, allow.new.levels = TRUE,
               type = "link")
  beta <- fixef(fit)
  expect_equal(unname(p),
               beta[["(Intercept)"]] + beta[["x"]] * c(0, 1),
               tolerance = 1e-8)
})

test_that("predict.mm_glmm refuses intervals and re.form subsets", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  expect_error(predict(fit, interval = "confidence"),
               class = "mm_inference_unavailable")
  expect_error(predict(fit, re.form = ~ (1 | g)),
               class = "mm_inference_unavailable")
})

test_that("predict.mm_glmm se.fit returns NA with a reason", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  out <- predict(fit, se.fit = TRUE)
  expect_true(is.list(out) && all(c("fit", "se.fit") %in% names(out)))
  expect_true(all(is.na(out$se.fit)))
  expect_identical(attr(out, "mm_unavailable_reason"),
                   "glmm_prediction_se_unavailable")
})
