# Audit1 WI-9.3 (finding I3): glmer-referenced GLMM parity previously covered
# only binomial-logit, poisson-log, and NB-log. This file closes the
# family/link grid with small, seeded, well-conditioned synthetic cells:
#
#   joint_laplace vs lme4::glmer(nAGQ = 1) parity:
#     binomial-probit, binomial-cloglog, poisson-sqrt, Gamma-log,
#     binomial-logit with a random slope (1 + x | g).
#   pirls_profiled self-validation (simulation recovery, NOT glmer parity):
#     the same new cells, fitted with the default profiled estimator on data
#     generated from known parameters.
#
# Tolerances mirror test-glmm-joint-laplace-parity.R: fixef within 5e-3,
# logLik within 5e-2. Observed margins at authoring (mixeff 0.2.0 installed,
# lme4 2.0.1, macOS arm64, R 4.5):
#   probit        max|fixef diff| 7.4e-06, |logLik diff| 8.2e-06
#   cloglog       max|fixef diff| 2.2e-04, |logLik diff| 1.9e-04
#   poisson-sqrt  max|fixef diff| 5.7e-06, |logLik diff| 7.1e-06
#   Gamma-log     max|fixef diff| 5.8e-04, |logLik diff| ~1.0 (see FINDING
#                 in the Gamma test: the logLik gap is a dispersion-handling
#                 convention difference, and glmer itself is off the exact
#                 marginal likelihood for Gamma; fixef-slope agreement is
#                 tight, so only fixef parity is asserted there)
#   logit-slope   max|fixef diff| 2.0e-05, |logLik diff| 4.6e-06,
#                 RE sd diffs < 5e-05, corr diff < 1e-04
#
# Every glmer reference below converges cleanly (no convergence warnings) on
# the authoring platform; DGPs were tuned for that rather than tolerances
# being loosened.

mm_cells_control <- function() mm_control(verbose = -1, max_feval = 50000L)

# fixef vectors are ordered (Intercept), then formula terms, identically in
# glmer and mixeff, so compare by position (the coefficient NAME encodings
# differ). Mirrors mm_expect_glmer_joint_parity() in
# test-glmm-joint-laplace-parity.R.
mm_cells_expect_joint_parity <- function(m, g, fixef_tol = 5e-3,
                                         loglik_tol = 5e-2) {
  bg <- unname(lme4::fixef(g))
  bm <- unname(fixef(m))
  testthat::expect_identical(m$method, "joint_laplace")
  testthat::expect_equal(length(bm), length(bg))
  testthat::expect_lt(max(abs(bm - bg)), fixef_tol)
  testthat::expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))),
                      loglik_tol)
}

# Center and rescale a random-effect draw so the realized effects have
# exactly zero mean and the nominal sd. Used by the recovery cells only:
# it removes the "luck of 30 normal draws" component of single-seed
# recovery error, so the bounds below measure estimator behavior (plus
# residual-family noise), not the RE realization.
mm_cells_center_scale <- function(v, sd_target) {
  v <- v - mean(v)
  v * sd_target / stats::sd(v)
}

## ---------------------------------------------------------------------------
## joint_laplace vs glmer(nAGQ = 1)
## ---------------------------------------------------------------------------

test_that("joint_laplace matches glmer on binomial-probit (random intercept)", {
  mm_skip_if_no_lme4()
  set.seed(101)
  ng <- 30L
  per <- 15L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.5)[as.integer(g)]
  y <- rbinom(n, 1, pnorm(0.3 + 0.6 * x + re))
  d <- data.frame(y = y, x = x, g = g)
  gl <- lme4::glmer(y ~ x + (1 | g), d, family = binomial(link = "probit"),
                    nAGQ = 1)
  m <- glmm(y ~ x + (1 | g), d, family = binomial(link = "probit"),
            method = "joint_laplace", control = mm_cells_control())
  mm_cells_expect_joint_parity(m, gl)
})

test_that("joint_laplace matches glmer on binomial-cloglog (random intercept)", {
  mm_skip_if_no_lme4()
  set.seed(102)
  ng <- 30L
  per <- 15L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.4)[as.integer(g)]
  y <- rbinom(n, 1, 1 - exp(-exp(-0.5 + 0.4 * x + re)))
  d <- data.frame(y = y, x = x, g = g)
  gl <- lme4::glmer(y ~ x + (1 | g), d, family = binomial(link = "cloglog"),
                    nAGQ = 1)
  m <- glmm(y ~ x + (1 | g), d, family = binomial(link = "cloglog"),
            method = "joint_laplace", control = mm_cells_control())
  mm_cells_expect_joint_parity(m, gl)
})

test_that("joint_laplace matches glmer on poisson-sqrt (random intercept)", {
  mm_skip_if_no_lme4()
  set.seed(103)
  ng <- 30L
  per <- 15L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.2)[as.integer(g)]
  # sqrt link: mu = eta^2; keep eta well away from 0 so the DGP is smooth.
  y <- rpois(n, (2 + 0.3 * x + re)^2)
  d <- data.frame(y = y, x = x, g = g)
  gl <- lme4::glmer(y ~ x + (1 | g), d, family = poisson(link = "sqrt"),
                    nAGQ = 1)
  m <- glmm(y ~ x + (1 | g), d, family = poisson(link = "sqrt"),
            method = "joint_laplace", control = mm_cells_control())
  mm_cells_expect_joint_parity(m, gl)
})

test_that("joint_laplace matches glmer fixed effects on Gamma-log (random intercept)", {
  mm_skip_if_no_lme4()
  set.seed(205)
  ng <- 40L
  per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.2)[as.integer(g)]
  mu <- exp(1 + 0.5 * x + re)
  y <- rgamma(n, shape = 15, rate = 15 / mu)
  d <- data.frame(y = y, x = x, g = g)
  gl <- lme4::glmer(y ~ x + (1 | g), d, family = Gamma(link = "log"),
                    nAGQ = 1)
  m <- glmm(y ~ x + (1 | g), d, family = Gamma(link = "log"),
            method = "joint_laplace", control = mm_cells_control())

  bg <- unname(lme4::fixef(gl))
  bm <- unname(fixef(m))
  expect_identical(m$method, "joint_laplace")
  expect_equal(length(bm), length(bg))
  expect_lt(max(abs(bm - bg)), 5e-3)
  # Observed at authoring: 5.8e-4 on this DGP.

  # FINDING (WI-9.3, Gamma-log logLik): the standard 5e-2 logLik tolerance is
  # NOT asserted for the Gamma family. Across a sweep of well-conditioned
  # Gamma DGPs (shape 5..50, RE sd 0.1..0.3, clean glmer convergence) the
  # mixeff-vs-glmer logLik gap is a stable ~+1.0 (mixeff higher), independent
  # of the DGP -- a systematic dispersion-handling convention difference, not
  # noise. Independent verification on the shape=5/sd=0.3/seed=104 DGP with a
  # 60-point Gauss-Hermite exact marginal likelihood (cross-checked by direct
  # optimization of that exact likelihood):
  #   exact-ML optimum:    beta (0.987, 0.482), RE sd 0.368, ll -737.11
  #   glmmTMB:             beta (0.987, 0.482), RE sd 0.368, ll -737.13
  #   glmer  reported ll:  -723.25 (exact ll at its estimates: -740.7)
  #   mixeff reported ll:  -722.14 (exact ll at its estimates: -742.7)
  # i.e. BOTH glmer and mixeff report Gamma "logLik" values that are not on
  # the exact marginal-likelihood scale (a known lme4 Gamma caveat; the
  # upstream mixeff-rs joint-optimizer contract explicitly states Gamma
  # dispersion handling makes objectives non-comparable across conventions).
  # glmer is therefore not a certifiable logLik reference for Gamma, and a
  # naive parity assertion would enshrine one arbitrary convention. The
  # assertion that WOULD be made under the standard tolerance is preserved
  # here, commented out, per the audit remediation protocol:
  # expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(gl))), 5e-2)
  # Observed gap at authoring: ~1.01 on this DGP.
  expect_true(is.finite(as.numeric(logLik(m))))

  # Gamma RE-sd parity is likewise not asserted: lme4 reports Gamma VarCorr
  # on a dispersion-scaled convention (glmer sd 0.548 vs mixeff 0.242 on the
  # seed-104 DGP, where glmer sd * sigma = 0.255 ~ mixeff's value), so the
  # two tables are not directly comparable.
})

test_that("joint_laplace matches glmer on a binomial-logit random-slope model", {
  mm_skip_if_no_lme4()
  set.seed(105)
  ng <- 40L
  per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  b0 <- rnorm(ng, sd = 0.8)
  b1 <- rnorm(ng, sd = 0.5)
  eta <- 0.4 + 0.7 * x + b0[as.integer(g)] + b1[as.integer(g)] * x
  y <- rbinom(n, 1, plogis(eta))
  d <- data.frame(y = y, x = x, g = g)
  gl <- lme4::glmer(y ~ x + (1 + x | g), d, family = binomial(), nAGQ = 1)
  m <- glmm(y ~ x + (1 + x | g), d, family = binomial(),
            method = "joint_laplace", control = mm_cells_control())
  mm_cells_expect_joint_parity(m, gl)

  # Random-effect structure also tracks glmer: intercept sd, slope sd, and
  # the intercept-slope correlation. Matched by (var1, var2) keys in the
  # lme4-shaped as.data.frame(VarCorr()) output.
  vg <- as.data.frame(lme4::VarCorr(gl))
  vm <- as.data.frame(VarCorr(m))
  pick <- function(tab, v1, v2) {
    row <- tab[tab$grp == "g" & tab$var1 == v1 &
                 (if (is.na(v2)) is.na(tab$var2) else !is.na(tab$var2) & tab$var2 == v2), ,
               drop = FALSE]
    expect_equal(nrow(row), 1L)
    row$sdcor[[1L]]
  }
  expect_lt(abs(pick(vm, "(Intercept)", NA) - pick(vg, "(Intercept)", NA)), 2e-2)
  expect_lt(abs(pick(vm, "x", NA) - pick(vg, "x", NA)), 2e-2)
  expect_lt(abs(pick(vm, "(Intercept)", "x") - pick(vg, "(Intercept)", "x")), 5e-2)
})

## ---------------------------------------------------------------------------
## pirls_profiled self-validation: single-seed simulation recovery
## ---------------------------------------------------------------------------
# The default profiled estimator is NOT glmer's estimator and is not expected
# to match glmer (see test-glmm-joint-laplace-parity.R), so it is validated
# against the data-generating truth instead. Scope of these checks, stated
# honestly: ONE seeded realization per cell, with the RE draw centered and
# rescaled to its nominal sd, and bounds set at roughly 3x the fit's own Wald
# SE (widened where a documented small-sample/link attenuation applies).
# They establish that the profiled estimator lands near the generating
# coefficients for that realization -- which catches gross family/link
# miswiring, sign flips, and scale errors -- and nothing more: they are NOT
# evidence of unbiasedness, of SE calibration, or of distributional
# correctness (any of those would need a many-replicate Monte Carlo design).
# Observed single-seed deviations at authoring are quoted per cell.

test_that("pirls_profiled recovers binomial-probit generating coefficients", {
  set.seed(111)
  ng <- 30L
  per <- 15L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- mm_cells_center_scale(rnorm(ng), 0.5)[as.integer(g)]
  y <- rbinom(n, 1, pnorm(0.3 + 0.6 * x + re))
  d <- data.frame(y = y, x = x, g = g)
  m <- glmm(y ~ x + (1 | g), d, family = binomial(link = "probit"),
            control = mm_cells_control())
  expect_identical(m$method, "pirls_profiled")
  est <- unname(fixef(m))
  # Observed deviations at authoring: +0.072, +0.062 (Wald SEs 0.133, 0.081).
  expect_lt(abs(est[[1L]] - 0.3), 0.40)
  expect_lt(abs(est[[2L]] - 0.6), 0.25)
})

test_that("pirls_profiled recovers binomial-cloglog generating coefficients", {
  set.seed(112)
  ng <- 30L
  per <- 15L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- mm_cells_center_scale(rnorm(ng), 0.4)[as.integer(g)]
  y <- rbinom(n, 1, 1 - exp(-exp(-0.5 + 0.4 * x + re)))
  d <- data.frame(y = y, x = x, g = g)
  m <- glmm(y ~ x + (1 | g), d, family = binomial(link = "cloglog"),
            control = mm_cells_control())
  expect_identical(m$method, "pirls_profiled")
  est <- unname(fixef(m))
  # Observed deviations at authoring: +0.047, -0.036 (Wald SEs 0.090, 0.073).
  expect_lt(abs(est[[1L]] - (-0.5)), 0.30)
  expect_lt(abs(est[[2L]] - 0.4), 0.25)
})

test_that("pirls_profiled recovers poisson-sqrt generating coefficients", {
  set.seed(113)
  ng <- 30L
  per <- 15L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- mm_cells_center_scale(rnorm(ng), 0.2)[as.integer(g)]
  y <- rpois(n, (2 + 0.3 * x + re)^2)
  d <- data.frame(y = y, x = x, g = g)
  m <- glmm(y ~ x + (1 | g), d, family = poisson(link = "sqrt"),
            control = mm_cells_control())
  expect_identical(m$method, "pirls_profiled")
  est <- unname(fixef(m))
  # Observed deviations at authoring: +0.044, +0.005 (Wald SEs 0.040, 0.025).
  expect_lt(abs(est[[1L]] - 2.0), 0.15)
  expect_lt(abs(est[[2L]] - 0.3), 0.10)
})

test_that("pirls_profiled recovers Gamma-log generating coefficients", {
  set.seed(114)
  ng <- 40L
  per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- mm_cells_center_scale(rnorm(ng), 0.2)[as.integer(g)]
  mu <- exp(1 + 0.5 * x + re)
  y <- rgamma(n, shape = 15, rate = 15 / mu)
  d <- data.frame(y = y, x = x, g = g)
  m <- glmm(y ~ x + (1 | g), d, family = Gamma(link = "log"),
            control = mm_cells_control())
  expect_identical(m$method, "pirls_profiled")
  est <- unname(fixef(m))
  # Observed deviations at authoring: +0.032, -0.010 (Wald SEs 0.015, 0.015).
  # The intercept bound is wider than 3x SE because Laplace-type Gamma
  # estimators carry a small intercept attenuation (see the exact-ML evidence
  # in the Gamma parity test above); 0.10 still bounds the error at 10% of
  # the true intercept.
  expect_lt(abs(est[[1L]] - 1.0), 0.10)
  expect_lt(abs(est[[2L]] - 0.5), 0.08)
})

test_that("pirls_profiled recovers a binomial-logit random-slope model", {
  set.seed(115)
  ng <- 40L
  per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  b0 <- mm_cells_center_scale(rnorm(ng), 0.8)
  b1 <- mm_cells_center_scale(rnorm(ng), 0.5)
  eta <- 0.4 + 0.7 * x + b0[as.integer(g)] + b1[as.integer(g)] * x
  y <- rbinom(n, 1, plogis(eta))
  d <- data.frame(y = y, x = x, g = g)
  m <- glmm(y ~ x + (1 + x | g), d, family = binomial(),
            control = mm_cells_control())
  expect_identical(m$method, "pirls_profiled")
  est <- unname(fixef(m))
  # Observed deviations at authoring: -0.091, -0.106 (Wald SEs 0.128, 0.168).
  expect_lt(abs(est[[1L]] - 0.4), 0.40)
  expect_lt(abs(est[[2L]] - 0.7), 0.50)

  # RE sds land in the neighborhood of the generating values (0.8, 0.5).
  # Observed at authoring: 0.751 and 0.499. Generous bound, same caveats as
  # above: one seed, coarse location check only.
  vm <- as.data.frame(VarCorr(m))
  sd_int <- vm$sdcor[vm$grp == "g" & vm$var1 == "(Intercept)" & is.na(vm$var2)]
  sd_slp <- vm$sdcor[vm$grp == "g" & vm$var1 == "x" & is.na(vm$var2)]
  expect_length(sd_int, 1L)
  expect_length(sd_slp, 1L)
  expect_lt(abs(sd_int - 0.8), 0.35)
  expect_lt(abs(sd_slp - 0.5), 0.35)
})
