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
  expect_equal(vg$var.distribution, sigma(fg)^2)
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
