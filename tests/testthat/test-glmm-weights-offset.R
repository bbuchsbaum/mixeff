# GLMM prior weights, offsets, and grouped-binomial responses. The engine
# supports all three (new_with_weights / new_with_offset / Family::Binomial);
# these tests cover the wrapper wiring through the FFI.

ctl <- mm_control(verbose = -1)

test_that("glmm() accepts an offset and applies it (Poisson rate model)", {
  testthat::skip_if_not_installed("lme4")
  set.seed(42)
  n <- 300
  d <- data.frame(x = rnorm(n), g = factor(rep(1:30, each = 10)),
                  expo = runif(n, 0.5, 2))
  b <- rnorm(30, sd = 0.4)
  d$y <- rpois(n, exp(0.3 + 0.5 * d$x + b[as.integer(d$g)] + log(d$expo)))

  fit <- glmm(y ~ x + (1 | g), d, family = poisson,
              offset = log(d$expo), control = ctl)
  ref <- suppressMessages(
    lme4::glmer(y ~ x + (1 | g) + offset(log(expo)), d, family = poisson))
  # Slope agrees closely; intercept within the documented fast-PIRLS gap.
  expect_equal(unname(fixef(fit)[["x"]]), unname(lme4::fixef(ref)[["x"]]),
               tolerance = 5e-2)
  expect_equal(fit$offset, log(d$expo))
})

test_that("glmm() fits grouped binomial via cbind(successes, failures)", {
  testthat::skip_if_not_installed("lme4")
  data(cbpp, package = "lme4")
  fit <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
              cbpp, family = binomial, control = ctl)
  ref <- suppressMessages(
    lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                cbpp, family = binomial))
  expect_s3_class(fit, "mm_glmm")
  # Grouped-binomial fixed effects within the documented PIRLS-vs-Laplace gap.
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)), tolerance = 0.1)
  # Trial counts are stored as the prior weights.
  expect_equal(fit$weights, cbpp$size)
})

test_that("cbind and proportion+weights spellings give the same fit", {
  data(cbpp, package = "lme4")
  cb <- cbpp
  cb$prop <- cb$incidence / cb$size
  by_cbind <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   cbpp, family = binomial, control = ctl)
  by_weights <- glmm(prop ~ period + (1 | herd), cb, family = binomial,
                     weights = cb$size, control = ctl)
  expect_equal(unname(fixef(by_cbind)), unname(fixef(by_weights)),
               tolerance = 1e-8)
})

test_that("cbind response and weights= cannot be combined", {
  data(cbpp, package = "lme4")
  expect_error(
    glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
         cbpp, family = binomial, weights = cbpp$size, control = ctl),
    class = "mm_arg_error"
  )
})

test_that("plain 0/1 binomial still routes to Bernoulli (no regression)", {
  testthat::skip_if_not_installed("lme4")
  set.seed(1)
  d <- data.frame(y = rbinom(200, 1, 0.4), x = rnorm(200),
                  g = factor(rep(1:20, each = 10)))
  fit <- glmm(y ~ x + (1 | g), d, family = binomial, control = ctl)
  ref <- suppressMessages(lme4::glmer(y ~ x + (1 | g), d, family = binomial))
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)), tolerance = 5e-2)
  expect_null(fit$weights)
})

test_that("invalid weights/offset are rejected with a typed condition", {
  set.seed(1)
  d <- data.frame(y = rpois(50, 2), x = rnorm(50), g = factor(rep(1:5, each = 10)))
  expect_error(glmm(y ~ x + (1 | g), d, family = poisson,
                    weights = rep(-1, 50), control = ctl), class = "mm_arg_error")
  expect_error(glmm(y ~ x + (1 | g), d, family = poisson,
                    weights = rep(1, 3), control = ctl), class = "mm_arg_error")
  expect_error(glmm(y ~ x + (1 | g), d, family = poisson,
                    offset = rep(Inf, 50), control = ctl), class = "mm_arg_error")
})
