# lme4 2.1-0 references for GLMMs with a free dispersion parameter (Gamma,
# inverse Gaussian, Gaussian with a non-identity link).
#
# lme4 2.1-0 changed how glmer() handles these families (GH #557, #643, #936:
# PIRLS weights carry 1/phi, phi is profiled, sigma() = sqrt(phi) with
# phi = deviance / (n - rank([X, Z])), theta is the absolute random-effect
# SD, logLik drops the family aic()'s +2). mixeff's engine follows 2.1-0, so
# glmer() from older lme4 (CI installs lme4 2.0-6 from CRAN) is no longer a
# valid reference for them. Tests therefore ask mm_disp_ref() for the
# reference: with lme4 >= 2.1-0 installed it is computed live from glmer(),
# otherwise it is read from the stored lme4 2.1-0 values in
# fixtures/lme4-2.1-dispersion-refs.json. Regenerate that file with
# fixtures/make-lme4-2.1-dispersion-refs.R (needs lme4 >= 2.1-0).

mm_disp_lme4_is_2_1 <- function() {
  requireNamespace("lme4", quietly = TRUE) &&
    utils::packageVersion("lme4") >= "2.1-0"
}

# ---- data generators (shared by the tests and the fixture script) ----------

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

mm_disp_data_cells_gamma <- function() {
  set.seed(205)
  ng <- 40L
  per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.2)[as.integer(g)]
  mu <- exp(1 + 0.5 * x + re)
  y <- rgamma(n, shape = 15, rate = 15 / mu)
  data.frame(y = y, x = x, g = g)
}

mm_disp_data_scale <- function() {
  set.seed(3)
  d <- data.frame(g = factor(rep(1:15, each = 10)), x = rnorm(150))
  d$y <- rgamma(150, shape = 4,
                rate = 4 / exp(0.5 + 0.3 * d$x + rnorm(15, sd = 0.4)[d$g]))
  d
}

mm_disp_data_influence <- function() {
  set.seed(481)
  g <- gl(8, 8)
  x <- rnorm(64)
  mu <- exp(1 + 0.3 * x + rnorm(8, sd = 0.4)[g])
  data.frame(x = x, g = g, yg = rgamma(64, shape = 4, rate = 4 / mu))
}

# ---- the cases ---------------------------------------------------------------

mm_disp_cases <- function() {
  fl <- function(family, gen) {
    list(data = function() mm_fl_data(gen), formula = y ~ x + (1 | g),
         family = family)
  }
  list(
    fl_gamma_inverse = fl(stats::Gamma(), function(x, re) {
      mu <- 1 / (2 + 0.3 * x + re)
      rgamma(length(x), shape = 15, rate = 15 / mu)
    }),
    fl_ig_log = fl(stats::inverse.gaussian("log"), function(x, re) {
      mm_fl_rig(length(x), exp(0.5 + 0.3 * x + re), 20)
    }),
    fl_ig_inverse = fl(stats::inverse.gaussian("inverse"), function(x, re) {
      mm_fl_rig(length(x), 1 / (2 + 0.2 * x + re), 20)
    }),
    fl_gaussian_log = fl(stats::gaussian("log"), function(x, re) {
      exp(1 + 0.3 * x + re) + rnorm(length(x), sd = 0.3)
    }),
    fl_gaussian_sqrt = fl(stats::gaussian("sqrt"), function(x, re) {
      (2 + 0.3 * x + re)^2 + rnorm(length(x), sd = 0.3)
    }),
    fl_gaussian_inverse = fl(stats::gaussian("inverse"), function(x, re) {
      100 / (2 + 0.2 * x + re) + rnorm(length(x), sd = 2)
    }),
    cells_gamma_log = list(data = mm_disp_data_cells_gamma,
                           formula = y ~ x + (1 | g),
                           family = stats::Gamma(link = "log")),
    scale_gamma_log = list(data = mm_disp_data_scale,
                           formula = y ~ x + (1 | g),
                           family = stats::Gamma(link = "log")),
    scale_ig_log = list(data = mm_disp_data_scale,
                        formula = y ~ x + (1 | g),
                        family = stats::inverse.gaussian(link = "log")),
    influence_gamma_log = list(data = mm_disp_data_influence,
                               formula = yg ~ x + (1 | g),
                               family = stats::Gamma(link = "log"),
                               hat = TRUE)
  )
}

mm_disp_case <- function(key) {
  case <- mm_disp_cases()[[key]]
  if (is.null(case)) stop(sprintf("unknown dispersion reference case `%s`", key))
  case
}

# ---- reference extraction (lme4 >= 2.1-0) ------------------------------------

# Hat matrix diagonal with the working weights inside and outside the
# factorization (mixeff's documented GLMM convention), from lme4's own
# factors; and the matching Cook's distances.
mm_disp_working_hat <- function(fit) {
  sqrtW <- Matrix::Diagonal(x = sqrt(stats::weights(fit, type = "working")))
  v <- lme4::getME(fit, c("L", "Lambdat", "Zt", "RX", "X", "RZX"))
  CL <- Matrix::solve(v$L, Matrix::solve(v$L, v$Lambdat %*% v$Zt %*% sqrtW,
                                         system = "P"), system = "L")
  CR <- Matrix::solve(t(v$RX), t(v$X) %*% sqrtW - Matrix::crossprod(v$RZX, CL))
  as.numeric(Matrix::colSums(CR^2) + Matrix::colSums(CL^2))
}

mm_disp_ref_live <- function(key) {
  case <- mm_disp_case(key)
  d <- case$data()
  fit <- suppressWarnings(suppressMessages(
    lme4::glmer(case$formula, d, family = case$family)
  ))
  ll <- stats::logLik(fit)
  out <- list(
    fixef = unname(lme4::fixef(fit)),
    theta = unname(lme4::getME(fit, "theta")),
    logLik = as.numeric(ll),
    df = as.integer(attr(ll, "df")),
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit),
    sigma = stats::sigma(fit),
    deviance = stats::deviance(fit)
  )
  if (isTRUE(case$hat)) {
    h <- mm_disp_working_hat(fit)
    r <- stats::residuals(fit, type = "pearson")
    p <- ncol(lme4::getME(fit, "X"))
    out$hat_working <- h
    out$cooks_working <- unname((r / (1 - h))^2 * h / (stats::sigma(fit)^2 * p))
  }
  out
}

mm_disp_fixture_path <- function() {
  candidates <- c(
    testthat::test_path("fixtures", "lme4-2.1-dispersion-refs.json"),
    file.path("tests", "testthat", "fixtures", "lme4-2.1-dispersion-refs.json")
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) stop("fixtures/lme4-2.1-dispersion-refs.json is missing")
  hit
}

mm_disp_ref_stored <- function(key) {
  refs <- jsonlite::fromJSON(mm_disp_fixture_path(), simplifyVector = TRUE)
  ref <- refs$cases[[key]]
  if (is.null(ref)) stop(sprintf("no stored lme4 2.1 reference for `%s`", key))
  ref
}

# The lme4 2.1-0 reference for `key`: live glmer() when lme4 >= 2.1-0 is
# installed, else the stored values. Requires lme4 only for the live path.
mm_disp_ref <- function(key) {
  if (mm_disp_lme4_is_2_1()) mm_disp_ref_live(key) else mm_disp_ref_stored(key)
}
