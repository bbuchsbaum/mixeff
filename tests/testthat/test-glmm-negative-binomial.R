# Negative-binomial GLMM parity (bd-01KT2N1GGD61RHZ68KZFM9KD49).
#
# Engine surface at pin 3b6ec69: NB2 with log link on the default profiled
# (PIRLS) path; fixed theta (MASS::negative.binomial style) and estimated
# theta (glmer.nb-style outer loop). The joint_laplace route is guarded off
# for NB upstream, so glmm() refuses it with a typed error.
#
# Tolerances: the profiled path is a documented estimator difference from
# glmer's Laplace fit (see the mm_estimator_notice), so fixef/logLik bounds
# are estimator-gap sized (observed ~8e-3 / ~0.24 on this fixture), not
# optimizer noise. Theta estimation additionally differs from glmer.nb's
# theta.ml alternation, so theta gets a relative bound.

nb_test_data <- function(seed = 42L) {
  set.seed(seed)
  n_g <- 30L
  n_per <- 20L
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n_g * n_per)
  b <- rnorm(n_g, sd = 0.5)
  mu <- exp(0.5 + 0.4 * x + b[as.integer(g)])
  y <- rnbinom(length(mu), size = 1.8, mu = mu)
  data.frame(y = y, x = x, g = g)
}

test_that("fixed-theta NB matches glmer(family = MASS::negative.binomial)", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")
  d <- nb_test_data()

  fit <- glmm(y ~ x + (1 | g), data = d,
              family = MASS::negative.binomial(1.8),
              control = mm_control(verbose = -1))
  ref <- suppressMessages(lme4::glmer(
    y ~ x + (1 | g), data = d, family = MASS::negative.binomial(1.8)
  ))

  expect_identical(names(fixef(fit)), names(lme4::fixef(ref)))
  expect_equal(fixef(fit), lme4::fixef(ref), tolerance = 2e-2)
  expect_equal(as.numeric(logLik(fit)), as.numeric(stats::logLik(ref)),
               tolerance = 1e-4)
  expect_identical(fit$family$family, "negative_binomial")
  expect_identical(fit$family$link, "log")
  expect_equal(fit$family$nb_theta, 1.8)
  expect_false(isTRUE(fit$family$nb_theta_estimated))
  expect_identical(fit_status(fit), "converged_interior")
})

test_that("estimated-theta NB tracks glmer.nb", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")
  d <- nb_test_data()

  fit <- glmm(y ~ x + (1 | g), data = d, family = mm_negative_binomial(),
              control = mm_control(verbose = -1))
  ref <- suppressMessages(suppressWarnings(
    lme4::glmer.nb(y ~ x + (1 | g), data = d)
  ))
  theta_ref <- lme4::getME(ref, "glmer.nb.theta")

  expect_identical(names(fixef(fit)), names(lme4::fixef(ref)))
  expect_equal(fixef(fit), lme4::fixef(ref), tolerance = 2e-2)
  # Theta estimation routes differ (engine outer loop vs glmer.nb's theta.ml
  # alternation); observed 1.704 vs 1.574 on this fixture.
  expect_equal(fit$family$nb_theta, theta_ref, tolerance = 0.15)
  expect_true(isTRUE(fit$family$nb_theta_estimated))
  expect_equal(as.numeric(logLik(fit)), as.numeric(stats::logLik(ref)),
               tolerance = 3e-4)

  # summary() renders without error and carries the coefficient table
  s <- summary(fit)
  expect_identical(rownames(s$coefficients), names(fixef(fit)))
})

test_that("NB accepts a fixed theta through mm_negative_binomial(theta = )", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("MASS")
  d <- nb_test_data()
  fit_a <- glmm(y ~ x + (1 | g), data = d, family = mm_negative_binomial(1.8),
                control = mm_control(verbose = -1))
  fit_b <- glmm(y ~ x + (1 | g), data = d,
                family = MASS::negative.binomial(1.8),
                control = mm_control(verbose = -1))
  expect_equal(fixef(fit_a), fixef(fit_b), tolerance = 1e-10)
  expect_equal(as.numeric(logLik(fit_a)), as.numeric(logLik(fit_b)),
               tolerance = 1e-10)
})

test_that("NB family validation refuses bad specs with typed errors", {
  d <- nb_test_data()

  expect_error(mm_negative_binomial(theta = -1), class = "mm_arg_error")
  expect_error(mm_negative_binomial(theta = c(1, 2)), class = "mm_arg_error")

  # joint_laplace is guarded off for NB at this engine pin.
  expect_error(
    glmm(y ~ x + (1 | g), data = d, family = mm_negative_binomial(),
         method = "joint_laplace", control = mm_control(verbose = -1)),
    class = "mm_inference_unavailable"
  )

  # Non-log links are outside the certified NB surface.
  fam_sqrt <- structure(
    list(family = "negative_binomial", link = "sqrt", theta = NULL),
    class = "family"
  )
  expect_error(
    glmm(y ~ x + (1 | g), data = d, family = fam_sqrt,
         control = mm_control(verbose = -1)),
    class = "mm_inference_unavailable"
  )
})
