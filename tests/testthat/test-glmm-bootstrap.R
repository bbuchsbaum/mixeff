# GLMM parametric bootstrap (Audit1 WI-2.8): the universal inference path.
# Simulate-and-refit replicates of the EFFECTIVE estimator; percentile
# intervals; replicate accounting attached, failures counted not dropped.

mm_boot_binom_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(7)
      g <- factor(rep(1:12, each = 12))
      x <- rnorm(144)
      eta <- -0.3 + 0.8 * x + rep(rnorm(12, sd = 0.5), each = 12)
      d <- data.frame(y = rbinom(144, 1, plogis(eta)), x = x, g = g)
      cache <<- glmm(y ~ x + (1 | g), d, family = binomial(),
                     control = mm_control(verbose = -1))
    }
    cache
  }
})

test_that("bootstrap confint works on a default profiled binomial GLMM", {
  fit <- mm_boot_binom_fit()
  ci <- confint(fit, method = "bootstrap", nsim = 60L, seed = 42L)
  expect_true(is.matrix(ci))
  expect_identical(rownames(ci), names(fixef(fit)))
  expect_true(all(is.finite(ci)))
  expect_true(all(ci[, 1L] < ci[, 2L]))

  acct <- attr(ci, "mm_bootstrap")
  expect_identical(acct$requested, acct$successful + acct$failed)
  expect_identical(acct$requested, 60L)
  expect_true(acct$successful >= 2L)
  expect_true(is.finite(acct$mcse))
  expect_true(acct$reliability$reliability %in%
                c("high", "moderate", "low", "not_available"))
  # nsim=60 is below the 999-replicate bar for a moderate grade.
  expect_identical(acct$reliability$reliability, "low")
  expect_identical(attr(ci, "mm_method"),
                   "glmm_parametric_bootstrap_percentile")
  expect_identical(attr(ci, "mm_estimator"), "pirls_profiled")
})

test_that("bootstrap confint is deterministic under an explicit seed", {
  fit <- mm_boot_binom_fit()
  ci1 <- confint(fit, method = "bootstrap", nsim = 25L, seed = 11L)
  ci2 <- confint(fit, method = "bootstrap", nsim = 25L, seed = 11L)
  expect_identical(unclass(ci1)[, ], unclass(ci2)[, ])
  ci3 <- confint(fit, method = "bootstrap", nsim = 25L, seed = 12L)
  expect_false(identical(unclass(ci1)[, ], unclass(ci3)[, ]))
})

test_that("NULL seed is governed by R's RNG (set.seed reproducibility)", {
  fit <- mm_boot_binom_fit()
  set.seed(303)
  ci1 <- confint(fit, method = "bootstrap", nsim = 20L)
  set.seed(303)
  ci2 <- confint(fit, method = "bootstrap", nsim = 20L)
  expect_identical(unclass(ci1)[, ], unclass(ci2)[, ])
})

test_that("bootstrap confint covers the negative-binomial family", {
  # NB has no joint-Laplace route, so the bootstrap is its ONLY interval
  # path. This template is FIXED-theta, so replicates condition on that
  # theta; an estimated-theta template would re-estimate it per replicate.
  skip_if_not_installed("MASS")
  set.seed(31)
  g <- factor(rep(1:10, each = 12))
  x <- rnorm(120)
  mu <- exp(0.6 + 0.4 * x + rep(rnorm(10, sd = 0.3), each = 12))
  y <- MASS::rnegbin(120, mu = mu, theta = 1.8)
  d <- data.frame(y = y, x = x, g = g)
  fit <- glmm(y ~ x + (1 | g), d, family = mm_negative_binomial(theta = 1.8),
              control = mm_control(verbose = -1))
  ci <- confint(fit, method = "bootstrap", nsim = 40L, seed = 5L)
  expect_true(all(is.finite(ci)))
  expect_true(all(ci[, 1L] < ci[, 2L]))
})

test_that("bootstrap SEs are on the scale of certified joint Wald SEs", {
  # Sanity anchor, not calibration: the replicate-sd SEs from the profiled
  # bootstrap should land within a broad factor of the certified joint
  # Wald SEs on a well-behaved model. Guards against unit errors (wrong
  # scale, wrong pivot order), not against fine miscalibration -- the DEC-2
  # simulation study owns the latter.
  fit <- mm_boot_binom_fit()
  ci <- confint(fit, method = "bootstrap", nsim = 80L, seed = 99L)
  boot_se <- attr(ci, "mm_bootstrap")$std_errors

  joint <- glmm(y ~ x + (1 | g), fit$model_frame, family = binomial(),
                method = "joint_laplace", control = mm_control(verbose = -1))
  wald_se <- summary(joint, tests = "coefficients")$coefficients[, "Std. Error"]
  ratio <- unname(boot_se[names(fixef(joint))] / wald_se)
  expect_true(all(is.finite(ratio)))
  expect_true(all(ratio > 0.5 & ratio < 2))
})

test_that("bootstrap argument validation is typed", {
  fit <- mm_boot_binom_fit()
  expect_error(confint(fit, method = "bootstrap", nsim = 1L),
               class = "mm_arg_error")
  expect_error(confint(fit, method = "bootstrap", nsim = 20L, seed = -1),
               class = "mm_arg_error")
})