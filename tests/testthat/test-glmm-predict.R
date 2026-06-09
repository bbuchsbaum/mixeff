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

test_that("predict.mm_glmm refuses prediction intervals and re.form subsets", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  # future-observation intervals require the family variance, which the
  # engine payload does not include for GLMMs
  expect_error(predict(fit, interval = "prediction"),
               class = "mm_inference_unavailable")
  expect_error(predict(fit, re.form = ~ (1 | g)),
               class = "mm_inference_unavailable")
})

test_that("pirls conditional se.fit is withheld with the engine reason", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  out <- predict(fit, se.fit = TRUE)  # default re.form = NULL (conditional)
  expect_true(is.list(out) && all(c("fit", "se.fit") %in% names(out)))
  # the default pirls_profiled fit only carries the uncertified working-delta
  # variance (engine status "degraded"), so the SE is withheld, not reported
  expect_true(all(is.na(out$se.fit)))
  reasons <- attr(out$se.fit, "mm_reason")
  expect_true(is.character(reasons) && all(!is.na(reasons)))
  # same for the confidence interval: fit column real, bounds withheld
  ci <- predict(fit, interval = "confidence")
  expect_true(is.matrix(ci))
  expect_true(all(is.finite(ci[, "fit"])))
  expect_true(all(is.na(ci[, "lwr"])) && all(is.na(ci[, "upr"])))
})

test_that("joint_laplace conditional se.fit and CI come from the engine", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              method = "joint_laplace",
              control = mm_control(verbose = -1))
  out <- predict(fit, type = "response", se.fit = TRUE)
  expect_true(all(is.finite(out$se.fit) & out$se.fit > 0))
  # response-scale bounds are symmetric delta-method bounds, so they bracket
  # the fit but are not guaranteed to stay inside (0, 1) near the boundary
  ci <- predict(fit, type = "response", interval = "confidence")
  expect_true(is.matrix(ci) && all(is.finite(ci)))
  expect_true(all(ci[, "lwr"] < ci[, "fit"] & ci[, "fit"] < ci[, "upr"]))
  # link-scale interval half-width must be z * link-scale se
  link_se <- predict(fit, type = "link", se.fit = TRUE)$se.fit
  link_ci <- predict(fit, type = "link", interval = "confidence")
  expect_equal(unname(link_ci[, "upr"] - link_ci[, "fit"]),
               unname(qnorm(0.975) * as.numeric(link_se)), tolerance = 1e-6)
})

test_that("predict.mm_glmm population se.fit gives finite Wald SEs", {
  fit <- glmm(y ~ x + (1 | g), mm_pred_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  link <- predict(fit, re.form = NA, type = "link", se.fit = TRUE)
  expect_true(all(is.finite(link$se.fit)))
  expect_true(all(link$se.fit > 0))
  resp <- predict(fit, re.form = NA, type = "response", se.fit = TRUE)
  expect_true(all(is.finite(resp$se.fit)))
  # response-scale SE = link SE * |d mu / d eta|, which is < link SE for logit
  expect_true(all(resp$se.fit <= link$se.fit + 1e-9))
})
