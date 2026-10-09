# Gamma / inverse-Gaussian GLMMs: the engine's log-likelihood is the density
# at phi = mean unit deviance; glmer's logLik() additionally carries the
# family aic()'s "+2", i.e. sits exactly 1 lower. mixeff reports lme4's
# number (and AIC / BIC / -2 logLik built from it).

test_that("Gamma and inverse-Gaussian logLik/AIC/BIC match glmer", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  set.seed(3)
  d <- data.frame(g = factor(rep(1:15, each = 10)), x = rnorm(150))
  d$y <- rgamma(150, shape = 4,
                rate = 4 / exp(0.5 + 0.3 * d$x + rnorm(15, sd = 0.4)[d$g]))
  for (fam in list(stats::Gamma(link = "log"),
                   stats::inverse.gaussian(link = "log"))) {
    fit <- glmm(y ~ x + (1 | g), d, family = fam,
                control = mm_control(verbose = -1))
    ref <- lme4::glmer(y ~ x + (1 | g), d, family = fam)
    lab <- fam$family
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
                 tolerance = 1e-5, label = lab)
    expect_equal(attr(logLik(fit), "df"), attr(logLik(ref), "df"))
    expect_equal(AIC(fit), AIC(ref), tolerance = 1e-5, label = lab)
    expect_equal(BIC(fit), BIC(ref), tolerance = 1e-5, label = lab)
    expect_equal(sigma(fit), sigma(ref), tolerance = 1e-3, label = lab)
    expect_equal(fixef(fit), fixef(ref), tolerance = 1e-3, label = lab)
    expect_equal(deviance(fit), deviance(ref), tolerance = 1e-3, label = lab)
  }
})

test_that("the scale-family logLik shift leaves likelihood-ratio tests unchanged", {
  expect_identical(mm_glmm_loglik_shift("gamma"), 1)
  expect_identical(mm_glmm_loglik_shift("inverse_gaussian"), 1)
  expect_identical(mm_glmm_loglik_shift("poisson"), 0)
  expect_identical(mm_glmm_loglik_shift("negative_binomial"), 0)
})
