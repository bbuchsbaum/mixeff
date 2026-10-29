test_that("bootstrap LRT uses each supplied ML fit's likelihood", {
  set.seed(146)
  g <- factor(rep(seq_len(12), each = 8))
  x <- rnorm(length(g))
  z <- rnorm(length(g))
  y <- 1 + 0.5 * x + 0.7 * z + rnorm(12, sd = 0.9)[as.integer(g)] +
    rnorm(length(g), sd = 0.4)
  d <- data.frame(y, x, z, g)
  null <- lmm(y ~ x + (1 | g), d, REML = FALSE,
              control = mm_control(verbose = -1, max_feval = 1L))
  alternative <- lmm(y ~ x + z + (1 | g), d, REML = FALSE,
                     control = mm_control(verbose = -1))
  null_refit <- lmm(y ~ x + (1 | g), d, REML = FALSE,
                    control = mm_control(verbose = -1))
  # Ensure this fixture distinguishes the supplied fit from a cold refit
  # under the alternative's controls.
  expect_gt(abs(as.numeric(logLik(null)) - as.numeric(logLik(null_refit))), 1e-7)
  boot <- parametric_bootstrap(null, alternative, nsim = 2L, seed = 101L)
  expect_equal(boot$observed,
               2 * (as.numeric(logLik(alternative)) - as.numeric(logLik(null))),
               tolerance = 1e-10)
  restored <- unserialize(serialize(null, NULL))
  boot_restored <- parametric_bootstrap(restored, alternative, nsim = 2L, seed = 101L)
  expect_equal(boot_restored$observed, boot$observed, tolerance = 1e-10)
})
