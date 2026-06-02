# Routine (non-env-gated) certification that the native joint_laplace GLMM
# estimator tracks lme4::glmer's joint-Laplace reference. This is the
# glmer-equivalent path; the default pirls_profiled is a different (faster,
# profiled) approximation and is NOT expected to match glmer here.
#
# Observed agreement at authoring (native joint-Laplace, mixeff-rs @8618e91+):
#   cbpp binomial : max|fixef diff| 5.9e-4, logLik diff 2.4e-4
#   poisson (sim) : max|fixef diff| 3.2e-4, logLik diff 3.1e-5
# Tolerances below sit comfortably above those margins.

mm_jl_control <- function() mm_control(verbose = -1, max_feval = 50000L)

# fixef vectors are ordered (Intercept), then formula terms, identically in
# glmer and mixeff, so compare by position (mixeff uses a different coefficient
# *name* encoding, e.g. "period: 2" vs "period2").
mm_expect_glmer_joint_parity <- function(m, g, fixef_tol = 5e-3,
                                         loglik_tol = 5e-2) {
  bg <- unname(lme4::fixef(g))
  bm <- unname(fixef(m))
  testthat::expect_equal(length(bm), length(bg))
  testthat::expect_lt(max(abs(bm - bg)), fixef_tol)
  testthat::expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))),
                      loglik_tol)
}

test_that("joint_laplace matches glmer on cbpp (binomial, grouped)", {
  skip_if_not_installed("lme4")
  data(cbpp, package = "lme4")
  g <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   cbpp, family = binomial)
  m <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cbpp,
            family = binomial(), method = "joint_laplace",
            control = mm_jl_control())
  expect_identical(m$method, "joint_laplace")
  mm_expect_glmer_joint_parity(m, g)
  # random-effect SD also tracks glmer
  sd_g <- sqrt(unname(unlist(lme4::VarCorr(g)$herd)))
  vc <- as.data.frame(VarCorr(m))
  sd_m <- vc$sdcor[vc$grp == "herd" & is.na(vc$var2)][1]
  expect_lt(abs(sd_m - sd_g), 2e-2)
})

test_that("joint_laplace matches glmer on a Poisson random-intercept model", {
  skip_if_not_installed("lme4")
  set.seed(11)
  ng <- 30L; per <- 10L
  g <- factor(rep(seq_len(ng), each = per)); n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.5)[as.integer(g)]
  y <- rpois(n, exp(0.5 + 0.3 * x + re))
  df <- data.frame(y = y, x = x, g = g)
  gl <- lme4::glmer(y ~ x + (1 | g), df, family = poisson)
  m <- glmm(y ~ x + (1 | g), df, family = poisson(),
            method = "joint_laplace", control = mm_jl_control())
  expect_identical(m$method, "joint_laplace")
  mm_expect_glmer_joint_parity(m, gl)
})

test_that("the default profiled estimator is labelled distinctly from joint_laplace", {
  data(cbpp, package = "lme4")
  m <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cbpp,
            family = binomial(), control = mm_control(verbose = -1))
  # default is the profiled path, a different estimator than glmer's joint Laplace
  expect_identical(m$method, "pirls_profiled")
})
