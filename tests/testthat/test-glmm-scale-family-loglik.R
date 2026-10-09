# Gamma / inverse-Gaussian GLMMs: lme4 2.1-0 dropped the family aic()'s "+2"
# from glmer's logLik() and profiles phi (sigma() = sqrt(phi),
# phi = deviance / (n - rank([X, Z]))). The engine reports the same numbers
# and mixeff passes them through unshifted. The reference is lme4 2.1-0:
# live glmer() when installed, else the stored values
# (helper-lme4-dispersion.R).

test_that("Gamma and inverse-Gaussian logLik/AIC/BIC match glmer", {
  skip_on_cran()
  for (key in c("scale_gamma_log", "scale_ig_log")) {
    case <- mm_disp_case(key)
    ref <- mm_disp_ref(key)
    fit <- glmm(case$formula, case$data(), family = case$family,
                control = mm_control(verbose = -1))
    lab <- key
    expect_equal(as.numeric(logLik(fit)), ref$logLik,
                 tolerance = 1e-5, label = lab)
    expect_equal(attr(logLik(fit), "df"), ref$df)
    expect_equal(AIC(fit), ref$AIC, tolerance = 1e-5, label = lab)
    expect_equal(BIC(fit), ref$BIC, tolerance = 1e-5, label = lab)
    expect_equal(sigma(fit), ref$sigma, tolerance = 1e-3, label = lab)
    expect_equal(unname(fixef(fit)), ref$fixef, tolerance = 1e-3, label = lab)
    expect_equal(deviance(fit), ref$deviance, tolerance = 1e-3, label = lab)
  }
})

test_that("scale-family logLik, deviance and AIC are the engine's, unshifted", {
  case <- mm_disp_case("scale_gamma_log")
  fit <- glmm(case$formula, case$data(), family = case$family,
              control = mm_control(verbose = -1))
  expect_identical(fit$logLik, as.numeric(fit$fit$log_likelihood))
  expect_identical(fit$AIC, as.numeric(fit$fit$aic))
  expect_identical(fit$BIC, as.numeric(fit$fit$bic))
  expect_false(exists("mm_glmm_loglik_shift", envir = asNamespace("mixeff"),
                      inherits = FALSE))
})
