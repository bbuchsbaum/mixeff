# Nakagawa R2 / ICC and performance/insight integration (pre-CRAN audit 4.8).

perf_ok <- function() {
  testthat::skip_if_not_installed("lme4")
  testthat::skip_if_not_installed("performance")
  testthat::skip_if_not_installed("insight")
}

test_that("mm_r2/mm_icc match performance on an LMM with random slopes", {
  perf_ok()
  data("sleepstudy", package = "lme4", envir = environment())
  m <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))
  r <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  expect_equal(unlist(mm_r2(m)), unlist(performance::r2_nakagawa(r)),
               tolerance = 1e-4)
  expect_equal(mm_icc(m), as.data.frame(unclass(performance::icc(r))),
               tolerance = 1e-4, ignore_attr = TRUE)
  vm <- insight::get_variance(m)
  vr <- insight::get_variance(r)
  comps <- setdiff(names(vr), "cor.slope_intercept")
  expect_equal(vm[comps], vr[comps], tolerance = 1e-4)
  # a small correlation: compare on the absolute scale
  expect_lt(abs(vm$cor.slope_intercept - vr$cor.slope_intercept), 1e-3)
  expect_equal(mm_variance_components(m), vm)
  # performance's own generics dispatch to mixeff fits.
  expect_equal(unlist(performance::r2(m)), unlist(mm_r2(m)))
  expect_equal(unclass(performance::icc(m))$ICC_adjusted,
               mm_icc(m)$ICC_adjusted)
  expect_true(insight::is_mixed_model(m))
})

test_that("GLMM R2/ICC match performance for binomial and Poisson", {
  perf_ok()
  # insight 1.5 changed the cbind() binomial and the null-model mean for the
  # log-normal approximation; mixeff follows the current convention.
  skip_if(utils::packageVersion("insight") < "1.5.0",
          "GLMM distribution-specific variance follows insight >= 1.5")
  data("cbpp", package = "lme4", envir = environment())
  f <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cbpp,
            family = binomial(), method = "joint_laplace",
            control = mm_control(verbose = -1))
  g <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   cbpp, family = binomial())
  expect_equal(unlist(mm_r2(f)), unlist(performance::r2_nakagawa(g)),
               tolerance = 1e-3)
  expect_equal(mm_icc(f)$ICC_adjusted, performance::icc(g)$ICC_adjusted,
               tolerance = 1e-3)

  data("grouseticks", package = "lme4", envir = environment())
  grouseticks$obs <- factor(seq_len(nrow(grouseticks)))
  f <- glmm(TICKS ~ YEAR + (1 | LOCATION) + (1 | obs), grouseticks,
            family = poisson(), method = "joint_laplace",
            control = mm_control(verbose = -1))
  g <- lme4::glmer(TICKS ~ YEAR + (1 | LOCATION) + (1 | obs), grouseticks,
                   family = poisson())
  r2m <- suppressWarnings(mm_r2(f))
  r2g <- suppressWarnings(performance::r2_nakagawa(g))
  expect_equal(unlist(r2m), unlist(r2g), tolerance = 1e-3)
  v <- suppressWarnings(mm_variance_components(f))
  expect_gt(v$var.dispersion, 0)
  expect_warning(mm_r2(f), "too close to zero")
})

test_that("Gamma and negative-binomial decompositions are finite", {
  set.seed(9)
  g <- factor(rep(1:15, each = 8))
  x <- rnorm(120)
  mu <- exp(2 + 0.4 * x + rep(rnorm(15, sd = 0.4), each = 8))
  d <- data.frame(yg = rgamma(120, shape = 5, rate = 5 / mu),
                  yn = rnbinom(120, mu = mu, size = 3), x = x, g = g)
  fg <- glmm(yg ~ x + (1 | g), d, family = Gamma(link = "log"),
             control = mm_control(verbose = -1))
  vg <- mm_variance_components(fg)
  # insight >= 1.5: log-normal approximation log1p(sigma^2) for a log link.
  expect_equal(vg$var.distribution, log1p(sigma(fg)^2))
  fgi <- glmm(yg ~ x + (1 | g), d, family = Gamma(link = "inverse"),
              control = mm_control(verbose = -1))
  expect_equal(mm_variance_components(fgi)$var.distribution, sigma(fgi)^2)
  # Live comparison needs lme4 >= 2.1-0, whose sigma() matches the engine's
  # dispersion for Gamma GLMMs, and insight >= 1.5's Gamma conventions.
  if (requireNamespace("lme4", quietly = TRUE) &&
      requireNamespace("insight", quietly = TRUE) &&
      requireNamespace("performance", quietly = TRUE) &&
      utils::packageVersion("lme4") >= "2.1-0" &&
      utils::packageVersion("insight") >= "1.5.0") {
    lg <- lme4::glmer(yg ~ x + (1 | g), d, family = Gamma(link = "log"))
    vl <- insight::get_variance(lg)
    comps <- intersect(names(vl), names(vg))
    expect_equal(vg[comps], vl[comps], tolerance = 1e-4)
    expect_equal(unlist(mm_r2(fg)), unlist(performance::r2_nakagawa(lg)),
                 tolerance = 1e-4)
  }
  fn <- glmm(yn ~ x + (1 | g), d, family = mm_negative_binomial(),
             control = mm_control(verbose = -1))
  r2 <- mm_r2(fn)
  expect_true(all(is.finite(unlist(r2))))
  expect_true(r2$R2_marginal < r2$R2_conditional)
})

test_that("singular fits and bad input are announced", {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = gl(10, 6))
  fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
  # A tolerance above every variance marks the fit singular.
  expect_warning(icc <- mm_icc(fit, tolerance = 1e6), "singular")
  expect_true(is.na(icc$ICC_adjusted))
  expect_warning(r2 <- mm_r2(fit, tolerance = 1e6), "singular")
  expect_true(is.na(r2$R2_conditional))
  expect_true(is.finite(r2$R2_marginal))
  expect_error(mm_r2(lm(y ~ x, d)), class = "mm_arg_error")
})

mm_perf_live_2_1 <- function() {
  requireNamespace("lme4", quietly = TRUE) &&
    requireNamespace("insight", quietly = TRUE) &&
    requireNamespace("performance", quietly = TRUE) &&
    utils::packageVersion("lme4") >= "2.1-0" &&
    utils::packageVersion("insight") >= "1.5.0"
}

test_that("Gaussian non-identity-link GLMMs use sigma^2, as insight 1.5 does", {
  for (key in c("fl_gaussian_log", "fl_gaussian_sqrt")) {
    case <- mm_disp_case(key)
    d <- case$data()
    m <- glmm(case$formula, d, family = case$family,
              control = mm_control(verbose = -1))
    v <- mm_variance_components(m)
    expect_equal(v$var.distribution, sigma(m)^2)
    expect_identical(v$var.dispersion, 0)
    expect_equal(v$var.residual, sigma(m)^2)
    r2 <- mm_r2(m)
    expect_true(all(is.finite(unlist(r2))))
    if (mm_perf_live_2_1()) {
      g <- suppressWarnings(lme4::glmer(case$formula, d, family = case$family))
      vl <- insight::get_variance(g)
      comps <- intersect(names(vl), names(v))
      expect_equal(v[comps], vl[comps], tolerance = 1e-4)
      expect_equal(unlist(r2), unlist(performance::r2_nakagawa(g)),
                   tolerance = 1e-4)
      expect_equal(mm_icc(m)$ICC_adjusted, performance::icc(g)$ICC_adjusted,
                   tolerance = 1e-4)
      # The registered methods give the same numbers on the mixeff fit.
      expect_equal(unlist(performance::r2(m)), unlist(r2))
    }
  }
})

test_that("inverse-Gaussian GLMMs are refused: insight's value is not a variance", {
  case <- mm_disp_case("fl_ig_log")
  d <- case$data()
  m <- glmm(case$formula, d, family = case$family,
            control = mm_control(verbose = -1))
  for (f in list(mm_variance_components, mm_r2, mm_icc)) {
    err <- expect_error(f(m), class = "mm_inference_unavailable")
    expect_identical(err$reason_code, "r2_distribution_variance_undefined")
  }
  if (requireNamespace("insight", quietly = TRUE)) {
    expect_warning(out <- insight::get_variance(m), "not defined")
    expect_true(is.na(out))
  }
  skip_if_not(mm_perf_live_2_1(), "needs lme4 >= 2.1-0 and insight >= 1.5.0")
  # What insight 1.5 does instead: its catch-all `resid.variance <- sig`,
  # i.e. sigma() = sqrt(phi) unsquared ...
  g <- suppressWarnings(lme4::glmer(case$formula, d, family = case$family))
  expect_equal(insight::get_variance(g)$var.distribution, stats::sigma(g))
  # ... which changes with the units of y although the log-link fixed and
  # random variances do not, so R2 / ICC would too.
  d100 <- transform(d, y = 100 * y)
  g100 <- suppressWarnings(lme4::glmer(case$formula, d100, family = case$family))
  v1 <- insight::get_variance(g)
  v100 <- insight::get_variance(g100)
  expect_equal(v100$var.fixed, v1$var.fixed, tolerance = 1e-3)
  expect_equal(v100$var.random, v1$var.random, tolerance = 1e-2)
  expect_equal(v100$var.distribution, v1$var.distribution / 10, tolerance = 1e-2)
})
