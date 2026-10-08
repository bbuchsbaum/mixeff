# Audit 2026-10 §4.7: every family/link pair the engine fits is exposed by
# glmm(), with glmer(nAGQ = 1) parity on fixed effects and theta. The
# log-likelihood is not compared: for dispersion families mixeff's logLik
# sits a stable ~+1.0 above glmer's (a dispersion-convention difference; see
# the Gamma FINDING in test-glmm-parity-cells.R). Pairs the engine lacks are
# refused with a typed condition naming the supported set.

mm_fl_rig <- function(n, mu, lambda) {
  nu <- rnorm(n)^2
  y <- mu + mu^2 * nu / (2 * lambda) -
    mu / (2 * lambda) * sqrt(4 * mu * lambda * nu + mu^2 * nu^2)
  ifelse(runif(n) <= mu / (mu + y), y, mu^2 / y)
}

mm_fl_data <- function(gen, seed = 11) {
  set.seed(seed)
  ng <- 30L
  per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  re <- rnorm(ng, sd = 0.2)[g]
  data.frame(y = gen(x, re), x = x, g = g)
}

mm_fl_expect_parity <- function(family, gen, fixef_tol = 1e-2,
                                theta_tol = 2e-2) {
  skip_if_not_installed("lme4")
  d <- mm_fl_data(gen)
  ref <- suppressWarnings(lme4::glmer(y ~ x + (1 | g), d, family = family))
  m <- glmm(y ~ x + (1 | g), d, family = family,
            control = mm_control(verbose = -1, max_feval = 50000L))
  expect_identical(m$method, "joint_laplace")
  expect_identical(m$family$link, family$link)
  expect_equal(unname(fixef(m)), unname(lme4::fixef(ref)),
               tolerance = fixef_tol)
  expect_equal(unname(m$theta), unname(lme4::getME(ref, "theta")),
               tolerance = theta_tol)
  invisible(m)
}

test_that("Gamma with its default inverse link matches glmer", {
  mm_fl_expect_parity(Gamma(), function(x, re) {
    mu <- 1 / (2 + 0.3 * x + re)
    rgamma(length(x), shape = 15, rate = 15 / mu)
  })
})

test_that("inverse.gaussian log and inverse links match glmer", {
  mm_fl_expect_parity(inverse.gaussian("log"), function(x, re) {
    mm_fl_rig(length(x), exp(0.5 + 0.3 * x + re), 20)
  })
  mm_fl_expect_parity(inverse.gaussian("inverse"), function(x, re) {
    mm_fl_rig(length(x), 1 / (2 + 0.2 * x + re), 20)
  })
})

test_that("gaussian non-identity links match glmer", {
  mm_fl_expect_parity(gaussian("log"), function(x, re) {
    exp(1 + 0.3 * x + re) + rnorm(length(x), sd = 0.3)
  }, fixef_tol = 1e-4, theta_tol = 1e-4)
  mm_fl_expect_parity(gaussian("sqrt"), function(x, re) {
    (2 + 0.3 * x + re)^2 + rnorm(length(x), sd = 0.3)
  }, fixef_tol = 1e-4, theta_tol = 1e-4)
  mm_fl_expect_parity(gaussian("inverse"), function(x, re) {
    100 / (2 + 0.2 * x + re) + rnorm(length(x), sd = 2)
  }, fixef_tol = 1e-3, theta_tol = 1e-2)
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
