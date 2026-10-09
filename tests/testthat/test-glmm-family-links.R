# Audit 2026-10 §4.7: every family/link pair the engine fits is exposed by
# glmm(), with glmer(nAGQ = 1) parity on fixed effects, theta, sigma and the
# log-likelihood. These are free-dispersion families, whose glmer() handling
# changed in lme4 2.1-0 (which the engine follows), so the reference is
# lme4 2.1-0: live glmer() when it is installed, else the stored lme4 2.1-0
# values (mm_disp_ref(), helper-lme4-dispersion.R). Pairs the engine lacks
# are refused with a typed condition naming the supported set.

mm_fl_expect_parity <- function(key, fixef_tol = 1e-2, theta_tol = 2e-2,
                                ll_tol = 1e-4) {
  case <- mm_disp_case(key)
  ref <- mm_disp_ref(key)
  m <- glmm(case$formula, case$data(), family = case$family,
            control = mm_control(verbose = -1, max_feval = 50000L))
  expect_identical(m$method, "joint_laplace")
  expect_identical(m$family$link, case$family$link)
  expect_equal(unname(fixef(m)), ref$fixef, tolerance = fixef_tol)
  expect_equal(unname(m$theta), ref$theta, tolerance = theta_tol)
  expect_equal(as.numeric(sigma(m)), ref$sigma, tolerance = theta_tol)
  expect_equal(as.numeric(logLik(m)), ref$logLik, tolerance = ll_tol)
  invisible(m)
}

test_that("Gamma with its default inverse link matches glmer", {
  mm_fl_expect_parity("fl_gamma_inverse")
})

test_that("inverse.gaussian log and inverse links match glmer", {
  mm_fl_expect_parity("fl_ig_log")
  mm_fl_expect_parity("fl_ig_inverse")
})

test_that("gaussian non-identity links match glmer", {
  mm_fl_expect_parity("fl_gaussian_log", fixef_tol = 1e-4, theta_tol = 1e-4)
  mm_fl_expect_parity("fl_gaussian_sqrt", fixef_tol = 1e-4, theta_tol = 1e-4)
  mm_fl_expect_parity("fl_gaussian_inverse", fixef_tol = 1e-3,
                      theta_tol = 1e-2)
})

test_that("stored lme4 2.1-0 dispersion references match live glmer", {
  skip_on_cran()
  skip_if_not(mm_disp_lme4_is_2_1(), "needs lme4 >= 2.1-0")
  for (key in names(mm_disp_cases())) {
    live <- mm_disp_ref_live(key)
    stored <- mm_disp_ref_stored(key)
    expect_equal(live[names(stored)], stored, tolerance = 1e-4, label = key)
  }
})

test_that("family/link pairs the engine lacks are refused, not swapped", {
  d <- mm_fl_data(function(x, re) rpois(length(x), 2) + 1)
  for (fam in list(binomial("cauchit"), binomial("log"), binomial("identity"),
                   poisson("identity"), inverse.gaussian())) {
    err <- tryCatch(
      glmm(y ~ x + (1 | g), d, family = fam, control = mm_control(verbose = -1)),
      error = identity
    )
    expect_s3_class(err, "mm_inference_unavailable")
    expect_match(conditionMessage(err), "Supported families", fixed = TRUE)
  }
  err <- tryCatch(
    glmm(y ~ x + (1 | g), d, family = gaussian(),
         control = mm_control(verbose = -1)),
    error = identity
  )
  expect_s3_class(err, "mm_inference_unavailable")
  expect_match(conditionMessage(err), "lmm()", fixed = TRUE)
})

test_that("update() reconstructs the new families", {
  m <- glmm(y ~ x + (1 | g),
            mm_fl_data(function(x, re) exp(1 + 0.3 * x + re) +
                         rnorm(length(x), sd = 0.3)),
            family = gaussian("log"), control = mm_control(verbose = -1))
  fam <- mixeff:::mm_glmm_family_from_info(m$family)
  expect_identical(c(fam$family, fam$link), c("gaussian", "log"))
  fam2 <- mixeff:::mm_glmm_family_from_info(list(family = "inverse_gaussian",
                                                 link = "log"))
  expect_identical(fam2$family, "inverse.gaussian")
})
