# PR-4 coverage: predict() argument refusals, empty newdata, offsets on
# newdata, re.form validation, and scaled GLMM residuals.

cp4pr_ctrl <- function() mm_control(verbose = -1)

cp4pr_data <- function(seed = 451L) {
  set.seed(seed)
  g <- gl(10, 6)
  x <- rnorm(60)
  o <- runif(60, 0, 1)
  y <- 1 + 0.5 * x + o + rnorm(10, sd = 0.7)[g] + rnorm(60, sd = 0.5)
  cnt <- rpois(60, exp(0.2 + 0.3 * x + rnorm(10, sd = 0.3)[g]))
  data.frame(y = y, x = x, o = o, g = g, cnt = cnt)
}

test_that("in-sample predictions refuse a new offset", {
  d <- cp4pr_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4pr_ctrl())
  expect_error(predict(fit, offset = rep(1, 60)), class = "mm_arg_error")
  gf <- glmm(cnt ~ x + (1 | g), d, family = poisson(), control = cp4pr_ctrl())
  expect_error(predict(gf, offset = rep(1, 60)), class = "mm_arg_error")
})

test_that("newdata without complete rows predicts NA and refuses SEs", {
  d <- cp4pr_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4pr_ctrl())
  nd <- data.frame(x = c(NA, NA), g = factor(c("1", "2"), levels = levels(d$g)))
  p <- predict(fit, newdata = nd)
  expect_identical(length(p), 2L)
  expect_true(all(is.na(p)))
  expect_error(predict(fit, newdata = nd, re.form = NA, se.fit = TRUE),
               class = "mm_data_error")
  gf <- glmm(cnt ~ x + (1 | g), d, family = poisson(), control = cp4pr_ctrl())
  pg <- predict(gf, newdata = nd, type = "response")
  expect_true(all(is.na(pg)))
  expect_error(predict(gf, newdata = nd, re.form = NA, se.fit = TRUE),
               class = "mm_data_error")
})

test_that("re.form and na.action arguments are validated", {
  d <- cp4pr_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4pr_ctrl())
  expect_error(predict(fit, re.form = "g"), class = "mm_arg_error")
  expect_error(predict(fit, re.form = 3), class = "mm_arg_error")
  nd <- data.frame(x = c(0.5, NA, -1),
                   g = factor(c("1", "2", "3"), levels = levels(d$g)))
  p_omit <- predict(fit, newdata = nd, na.action = "na.omit")
  expect_identical(names(p_omit), c("1", "3"))
  p_pass <- predict(fit, newdata = nd)
  expect_equal(unname(p_pass[c(1, 3)]), unname(p_omit))
  expect_true(is.na(p_pass[[2]]))
})

test_that("offset = argument fits need a matching newdata offset", {
  skip_if_not_installed("lme4")
  d <- cp4pr_data()
  fit <- lmm(y ~ x + (1 | g), d, offset = o, control = cp4pr_ctrl())
  nd <- d[1:4, ]
  expect_error(predict(fit, newdata = nd, offset = 1:3), class = "mm_arg_error")
  expect_error(predict(fit, newdata = nd, offset = letters[1:4]),
               class = "mm_arg_error")
  p <- predict(fit, newdata = nd, offset = o)
  ref <- lme4::lmer(y ~ x + (1 | g), d, offset = o)
  # lme4 drops an `offset =` argument on newdata; adding it back must agree.
  expect_equal(unname(p), unname(predict(ref, newdata = nd)) + nd$o,
               tolerance = 1e-4)
  expect_equal(unname(p), unname(fitted(fit)[1:4]), tolerance = 1e-6)
  # A legacy fit object (no `offset_arg` slot) reads the offset from `$offset`.
  legacy <- fit
  legacy$offset_arg <- NULL
  expect_false("offset_arg" %in% names(legacy))
  expect_equal(mixeff:::mm_fit_offset_arg(legacy), legacy$offset)
})

test_that("prediction maps skip absent components", {
  res <- list(fit = c(a = 1, b = 2), se.fit = NULL)
  out <- mixeff:::mm_prediction_map(res, function(x) x * 10)
  expect_identical(out$fit, c(a = 10, b = 20))
  expect_null(out$se.fit)
})

test_that("scaled GLMM residuals divide by sigma", {
  set.seed(452)
  g <- gl(8, 8)
  x <- rnorm(64)
  mu <- exp(0.5 + 0.3 * x + rnorm(8, sd = 0.3)[g])
  d <- data.frame(y = rgamma(64, shape = 4, rate = 4 / mu), x = x, g = g)
  fit <- glmm(y ~ x + (1 | g), d, family = Gamma(link = "log"),
              control = cp4pr_ctrl())
  r <- residuals(fit, type = "pearson")
  rs <- residuals(fit, type = "pearson", scaled = TRUE)
  expect_equal(rs, r / sigma(fit))
  expect_false(isTRUE(all.equal(rs, r)))
})
