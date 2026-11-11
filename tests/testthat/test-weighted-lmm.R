# Case-weighted LMMs (mixeff-rs issues #15 and #16).
#
# #15: with sum(log(weights)) > 0 the engine's optimizers compared trial
#      points on a different scale from the start value, so a weighted fit
#      with vector theta returned the starting theta.
# #16: test_random_effect() on the only random term compared the weighted
#      LMM likelihood with an unweighted linear model.

mk_weighted_lmm_data <- function(seed = 20261109L) {
  set.seed(seed)
  n_groups <- 15L
  n_per <- 16L
  g <- factor(rep(sprintf("g%02d", seq_len(n_groups)), each = n_per))
  x <- rep(seq(-1.5, 1.5, length.out = n_per), n_groups)
  b0 <- rnorm(n_groups, sd = 0.7)
  b1 <- rnorm(n_groups, sd = 0.3)
  w <- runif(length(x), 0.5, 2)
  y <- 1 + b0[as.integer(g)] + (0.5 + b1[as.integer(g)]) * x +
    rnorm(length(x), sd = 0.4) / sqrt(w)
  data.frame(y = y, x = x, g = g, w = w)
}

test_that("weighted fits do not depend on the scale of the weights", {
  df <- mk_weighted_lmm_data()
  expect_gt(sum(log(df$w)), 0)

  scales <- c(0.1, 1, 10)
  fits <- lapply(scales, function(scale) {
    d <- df
    d$w <- d$w * scale
    lmm(y ~ x + (1 + x | g), d, weights = w, REML = FALSE,
        control = mm_control(verbose = -1))
  })
  loglik <- vapply(fits, function(f) as.numeric(logLik(f)), numeric(1))
  vc <- lapply(fits, function(f) as.data.frame(VarCorr(f)))
  # Random-effect SDs and correlation are invariant; the residual SD is
  # sigma for unit weight, so it scales with sqrt(scale).
  ranef_sdcor <- lapply(vc, function(v) v$sdcor[v$grp != "Residual"])
  residual_sd <- vapply(vc, function(v) v$sdcor[v$grp == "Residual"], numeric(1))

  expect_equal(loglik[[2L]], loglik[[1L]], tolerance = 1e-6)
  expect_equal(loglik[[3L]], loglik[[1L]], tolerance = 1e-6)
  expect_equal(ranef_sdcor[[2L]], ranef_sdcor[[1L]], tolerance = 1e-4)
  expect_equal(ranef_sdcor[[3L]], ranef_sdcor[[1L]], tolerance = 1e-4)
  expect_equal(residual_sd / sqrt(scales), rep(residual_sd[[2L]], 3L),
               tolerance = 1e-4)

  # The unit-scale fit must have left the start theta = (1, 0, 1).
  pattern <- lmm(y ~ x + (1 + x | g), df, weights = w, REML = FALSE,
                 control = mm_control(verbose = -1,
                                      optimizer = "pattern_search"))
  expect_equal(loglik[[2L]], as.numeric(logLik(pattern)), tolerance = 1e-6)
})

test_that("weighted fits match lme4", {
  skip_if_not_installed("lme4")
  df <- mk_weighted_lmm_data()
  fit <- lmm(y ~ x + (1 + x | g), df, weights = w, REML = FALSE,
             control = mm_control(verbose = -1))
  ref <- lme4::lmer(y ~ x + (1 + x | g), df, weights = w, REML = FALSE)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)),
               tolerance = 1e-4)
})

test_that("test_random_effect() on a weighted model uses the weighted reference", {
  df <- mk_weighted_lmm_data()
  for (weighted in c(FALSE, TRUE)) {
    fit <- if (weighted) {
      lmm(y ~ x + (1 | g), df, weights = w, REML = FALSE,
          control = mm_control(verbose = -1))
    } else {
      lmm(y ~ x + (1 | g), df, REML = FALSE,
          control = mm_control(verbose = -1))
    }
    # logLik.lm includes 0.5 * sum(log(w)), the weighted-LMM scale.
    reduced <- if (weighted) {
      stats::lm(y ~ x, df, weights = w)
    } else {
      stats::lm(y ~ x, df)
    }
    expected <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(reduced)))

    row <- test_random_effect(fit, "g", method = "boundary_lrt")$table[1L, ]
    expect_identical(row$status, "available")
    expect_equal(row$statistic, expected, tolerance = 1e-6)
  }
})
