# simulate() for GLMMs with a free dispersion parameter (Gamma, inverse
# Gaussian, Gaussian with a non-identity link) against lme4 2.1-0's
# simulate.merMod(). mixeff consumes the RNG stream in lme4's order (all
# random-effect draws, then the response draws, with the same samplers), so
# with matching estimates the Gamma and Gaussian draws agree draw for draw:
# the check is exact equality of seeded draws up to the estimates' relative
# agreement (1e-4), against live glmer() under lme4 >= 2.1-0 and otherwise
# against stored lme4 2.1-0 draws (fixtures/lme4-2.1-dispersion-refs.json).
# The inverse Gaussian differs on purpose (lme4 uses shape w / sigma); its
# check is distributional.

sim_fit <- function(key, data = NULL) {
  case <- mm_disp_case(key)
  d <- data %||% case$data()
  glmm(case$formula, d, family = case$family,
       control = mm_control(verbose = -1))
}

sim_head <- function(s) {
  rows <- seq_len(mm_disp_sim_rows)
  c(s[[1]][rows], s[[2]][rows])
}

test_that("seeded Gamma and Gaussian-link draws equal lme4 2.1-0's", {
  for (key in c("cells_gamma_log", "fl_gamma_inverse", "fl_gaussian_log")) {
    m <- sim_fit(key)
    ref <- mm_disp_ref(key)
    expect_equal(sim_head(simulate(m, nsim = 2, seed = 7)), ref$sim_new,
                 tolerance = 1e-4, info = key)
    expect_equal(sim_head(simulate(m, nsim = 2, seed = 7, use.u = TRUE)),
                 ref$sim_cond, tolerance = 1e-4, info = key)
  }
})

test_that("full seeded draws match live lme4 2.1-0 simulate()", {
  skip_if_not(mm_disp_lme4_is_2_1(), "needs lme4 >= 2.1-0")
  for (key in c("cells_gamma_log", "fl_gaussian_log")) {
    case <- mm_disp_case(key)
    d <- case$data()
    m <- glmm(case$formula, d, family = case$family,
              control = mm_control(verbose = -1))
    g <- suppressWarnings(lme4::glmer(case$formula, d, family = case$family))
    for (re.form in list(NA, NULL, ~ (1 | g))) {
      a <- as.matrix(simulate(m, nsim = 3, seed = 11, re.form = re.form))
      b <- as.matrix(stats::simulate(g, nsim = 3, seed = 11, re.form = re.form))
      expect_equal(unname(a), unname(b), tolerance = 1e-4, info = key)
    }
  }
})

test_that("simulated moments follow the fitted variance functions", {
  # Conditional on the BLUPs the draws have mean fitted() and variance
  # phi * V(mu) / w, with phi = sigma()^2: V(mu) = mu^2 (Gamma), mu^3
  # (inverse Gaussian), 1 (Gaussian).
  skip_if_not_installed("statmod")
  vfun <- list(cells_gamma_log = function(mu) mu^2,
               fl_ig_log = function(mu) mu^3,
               fl_gaussian_log = function(mu) rep(1, length(mu)))
  for (key in names(vfun)) {
    m <- sim_fit(key)
    mu <- as.numeric(fitted(m))
    s <- as.matrix(simulate(m, nsim = 2000, seed = 3, use.u = TRUE))
    expect_equal(mean(rowMeans(s) / mu), 1, tolerance = 0.01, info = key)
    ratio <- apply(s, 1, stats::var) / vfun[[key]](mu)
    expect_equal(mean(ratio), sigma(m)^2, tolerance = 0.03, info = key)
  }
  # Prior weights divide the variance (lme4 ignores them for the Gaussian).
  case <- mm_disp_case("fl_gaussian_log")
  d <- case$data()
  d$w <- rep(c(1, 4), length.out = nrow(d))
  mw <- glmm(case$formula, d, family = case$family, weights = w,
             control = mm_control(verbose = -1))
  s <- as.matrix(simulate(mw, nsim = 2000, seed = 4, use.u = TRUE))
  v <- apply(s, 1, stats::var) * d$w
  expect_equal(mean(v), sigma(mw)^2, tolerance = 0.03)
})

test_that("inverse-Gaussian draws use shape w / phi (lme4 2.1-0 uses w / sigma)", {
  skip_if_not_installed("statmod")
  m <- sim_fit("fl_ig_log")
  mu <- as.numeric(fitted(m))
  # Exactly statmod::rinvgauss() with shape w / phi after the same seed.
  set.seed(5)
  expected <- statmod::rinvgauss(length(mu), mean = mu, shape = 1 / sigma(m)^2)
  expect_equal(simulate(m, seed = 5, use.u = TRUE)[[1]], expected)
  skip_if_not(mm_disp_lme4_is_2_1(), "needs lme4 >= 2.1-0")
  case <- mm_disp_case("fl_ig_log")
  g <- suppressWarnings(lme4::glmer(case$formula, case$data(),
                                    family = case$family))
  # lme4's draws are statmod::rinvgauss() with shape 1 / sigma() instead ...
  set.seed(5)
  lme4_like <- statmod::rinvgauss(length(mu), mean = mu, shape = 1 / sigma(m))
  expect_equal(stats::simulate(g, seed = 5, use.u = TRUE)[[1]], lme4_like,
               tolerance = 1e-4)
  # ... so their variance is sigma * mu^3, not the fitted phi * mu^3.
  s <- as.matrix(stats::simulate(g, nsim = 2000, seed = 3, use.u = TRUE))
  expect_equal(mean(apply(s, 1, stats::var) / mu^3), sigma(m), tolerance = 0.03)
})
