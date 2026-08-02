# WI-9.4 (Audit1): equivalent-representation invariants. Each test asserts
# one algebraic identity between two spellings/presentations of the same
# model on the same data, fitted by the same estimator, so tolerances are
# tight (1e-8 unless a second optimizer run is involved).
#
# Inventory of pre-existing overlapping coverage (not duplicated here):
#   * cbind(s, f) vs proportion + weights fixef equality:
#     test-glmm-weights-offset.R ("cbind and proportion+weights spellings
#     give the same fit") -- only the missing logLik assertion is added
#     below.
#   * saveRDS/readRDS + revive() in-sample predict equality:
#     test-phase2-revive.R ("saveRDS/readRDS revival keeps durable
#     extractors usable", "serialized fits restart in a fresh R process")
#     and artifact-row preservation in test-inference.R ("saved fits
#     preserve Rust artifact inference rows") -- the newdata prediction
#     path after a plain readRDS (no explicit revive()) is added below.
#   * (1 | a/b) canonicalization to (1 | a) + (1 | a:b) at the
#     compile/diagnostic level: test-explain-model.R (nested canonical
#     text) and test-explain-model-snapshots.R (syntax_expansion
#     diagnostic); test-pw2-crossed-nested.R fits the explicit
#     (1 | team) + (1 | team:player_id) spelling against lme4 but never
#     compares the two mixeff spellings to each other -- the fit-level
#     identity is added below.
#   * `||` semantics at the compile level, including the factor-`||`
#     family divergence: test-compile-audit.R ("`(1 + x || g)` compiles
#     independent intercept and slope kernels", "a factor inside || emits
#     the double_bar_factor_term diagnostic"); the aphantasia suite uses
#     the explicit expansion formula but never asserts the numeric-`||`
#     fit-level identity -- added below.
#   * Row-shuffled PREDICTION consistency (not fit invariance):
#     test-glmm-predict.R ("predict.mm_glmm conditional newdata reproduces
#     in-sample fitted").
# No existing test covers row-permutation fit invariance, predictor
# centering, response scaling, or factor releveling.

mm_inv_ctl <- function() mm_control(verbose = -1)

# Shared LMM fixture: 60 rows, 10 groups, one covariate, one 2-level factor.
mm_inv_data <- function() {
  set.seed(101)
  n <- 60L
  df <- data.frame(
    g = factor(rep_len(sprintf("g%02d", 1:10), n)),
    x = round(rnorm(n), 3),
    A = factor(rep_len(c("a1", "a2", "a2"), n))
  )
  u <- rnorm(10, sd = 0.8)
  df$y <- round(
    2 + 0.5 * df$x + 0.9 * (df$A == "a2") + u[as.integer(df$g)] +
      rnorm(n, sd = 0.7),
    4
  )
  df
}

mm_inv_fit <- function(formula, data) {
  lmm(formula, data, REML = FALSE, control = mm_inv_ctl())
}

test_that("row-permutation invariance: shuffled rows give the identical fit", {
  # Identity: the marginal likelihood is a product over observations
  # (given the group structure), so permuting the rows of `data` must not
  # change the optimum: logLik, fixef, theta, and sigma are invariant.
  df <- mm_inv_data()
  fit <- mm_inv_fit(y ~ x + A + (1 | g), df)
  set.seed(7)
  perm <- df[sample(nrow(df)), , drop = FALSE]
  fit_p <- mm_inv_fit(y ~ x + A + (1 | g), perm)

  expect_equal(as.numeric(logLik(fit_p)), as.numeric(logLik(fit)),
               tolerance = 1e-8)
  expect_equal(fixef(fit_p), fixef(fit), tolerance = 1e-8)
  expect_equal(getME(fit_p, "theta"), getME(fit, "theta"), tolerance = 1e-8)
  expect_equal(sigma(fit_p), sigma(fit), tolerance = 1e-8)
})

test_that("predictor centering: slope invariant, intercept shifts by beta1 * mean(x)", {
  # Identity: for xc = x - mean(x) the model is a reparameterization of the
  # same column space, so logLik, sigma, and theta are unchanged, the slope
  # is identical, and intercept_c = intercept + slope * mean(x).
  df <- mm_inv_data()
  fit <- mm_inv_fit(y ~ x + A + (1 | g), df)
  dc <- df
  dc$xc <- dc$x - mean(dc$x)
  fit_c <- mm_inv_fit(y ~ xc + A + (1 | g), dc)

  expect_equal(as.numeric(logLik(fit_c)), as.numeric(logLik(fit)),
               tolerance = 1e-8)
  expect_equal(fixef(fit_c)[["xc"]], fixef(fit)[["x"]], tolerance = 1e-8)
  expect_equal(
    fixef(fit_c)[["(Intercept)"]],
    fixef(fit)[["(Intercept)"]] + fixef(fit)[["x"]] * mean(df$x),
    tolerance = 1e-8
  )
  expect_equal(fixef(fit_c)[["Aa2"]], fixef(fit)[["Aa2"]], tolerance = 1e-8)
  expect_equal(sigma(fit_c), sigma(fit), tolerance = 1e-8)
  expect_equal(getME(fit_c, "theta"), getME(fit, "theta"), tolerance = 1e-8)
})

test_that("response scaling: y -> c*y scales fixef/sigma and shifts logLik by -n*log|c|", {
  # Identity: if ys = c * y then the ML fit transforms as
  #   fixef_s = c * fixef,  sigma_s = |c| * sigma,  theta unchanged
  #   (theta is the relative covariance factor, scale-free),
  #   logLik_s = logLik - n * log|c|   (change of variables: the density of
  #   ys is f_y(ys / c) / |c|^n), equivalently
  #   deviance_s = deviance + 2 * n * log|c|.
  df <- mm_inv_data()
  fit <- mm_inv_fit(y ~ x + A + (1 | g), df)
  cc <- -2.5
  ds <- df
  ds$ys <- cc * ds$y
  fit_s <- mm_inv_fit(ys ~ x + A + (1 | g), ds)
  n <- nrow(df)

  expect_equal(unname(fixef(fit_s)), unname(cc * fixef(fit)),
               tolerance = 1e-8)
  expect_equal(sigma(fit_s), abs(cc) * sigma(fit), tolerance = 1e-8)
  expect_equal(getME(fit_s, "theta"), getME(fit, "theta"), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit_s)),
               as.numeric(logLik(fit)) - n * log(abs(cc)),
               tolerance = 1e-8)
  expect_equal(deviance(fit_s), deviance(fit) + 2 * n * log(abs(cc)),
               tolerance = 1e-8)
})

test_that("factor releveling: transformed coefficients reproduce the same cell means", {
  # Identity: releveling changes the treatment-coded basis by an invertible
  # linear map, so the fit (logLik, sigma) is unchanged and the implied cell
  # means mu(a, b) = X(a, b) beta -- the coefficient map evaluated on the
  # full factorial grid -- are identical across codings.
  set.seed(23)
  df <- mm_inv_data()
  df$B <- factor(rep_len(c("b1", "b2", "b3", "b3", "b2"), nrow(df)))
  df$y2 <- round(df$y + 0.6 * (df$B == "b2") -
                   0.8 * (df$A == "a2" & df$B == "b3"), 4)
  fit <- mm_inv_fit(y2 ~ A * B + (1 | g), df)

  dr <- df
  dr$A <- stats::relevel(dr$A, "a2")
  dr$B <- factor(as.character(dr$B), levels = c("b3", "b1", "b2"))
  fit_r <- mm_inv_fit(y2 ~ A * B + (1 | g), dr)

  expect_equal(as.numeric(logLik(fit_r)), as.numeric(logLik(fit)),
               tolerance = 1e-8)
  expect_equal(sigma(fit_r), sigma(fit), tolerance = 1e-8)

  grid <- expand.grid(A = levels(df$A), B = levels(df$B),
                      KEEP.OUT.ATTRS = FALSE)
  grid_r <- grid
  grid_r$A <- factor(as.character(grid$A), levels = levels(dr$A))
  grid_r$B <- factor(as.character(grid$B), levels = levels(dr$B))
  mu <- predict(fit, newdata = grid, re.form = NA)
  mu_r <- predict(fit_r, newdata = grid_r, re.form = NA)
  expect_equal(unname(mu_r), unname(mu), tolerance = 1e-8)
})

test_that("grouped binomial: cbind and proportion+weights spellings share the logLik", {
  # Identity: cbind(s, n - s) and (s/n with weights = n) encode the same
  # grouped-binomial likelihood, so the maximized logLik must be equal.
  # Extends test-glmm-weights-offset.R ("cbind and proportion+weights
  # spellings give the same fit"), which pins fixef equality only.
  mm_skip_if_no_lme4()
  data(cbpp, package = "lme4")
  cb <- cbpp
  cb$prop <- cb$incidence / cb$size
  by_cbind <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   cbpp, family = binomial, control = mm_inv_ctl())
  by_weights <- glmm(prop ~ period + (1 | herd), cb, family = binomial,
                     weights = cb$size, control = mm_inv_ctl())

  expect_equal(as.numeric(logLik(by_cbind)), as.numeric(logLik(by_weights)),
               tolerance = 1e-8)
  expect_equal(unname(fitted(by_cbind)), unname(fitted(by_weights)),
               tolerance = 1e-8)
})

test_that("nested notation (1 | a/b) equals explicit (1 | a) + (1 | a:b)", {
  # Identity: the `/` nesting operator is defined as sugar for the explicit
  # two-kernel spelling, so both formulas must reach the identical fit.
  # (test-explain-model.R checks the canonicalization text; this pins the
  # numbers.)
  set.seed(11)
  nn <- expand.grid(a = sprintf("s%02d", 1:6), b = c("u", "v", "w"),
                    rep = 1:4, KEEP.OUT.ATTRS = FALSE)
  nn$a <- factor(nn$a)
  nn$b <- factor(nn$b) # levels shared across `a`: nesting must be explicit
  ua <- rnorm(6, sd = 1.0)
  uab <- rnorm(18, sd = 0.6)
  nn$yy <- round(3 + ua[as.integer(nn$a)] +
                   uab[(as.integer(nn$a) - 1L) * 3L + as.integer(nn$b)] +
                   rnorm(nrow(nn), sd = 0.5), 4)

  fit_slash <- mm_inv_fit(yy ~ 1 + (1 | a / b), nn)
  fit_expl <- mm_inv_fit(yy ~ 1 + (1 | a) + (1 | a:b), nn)

  expect_equal(as.numeric(logLik(fit_slash)), as.numeric(logLik(fit_expl)),
               tolerance = 1e-8)
  expect_equal(fixef(fit_slash), fixef(fit_expl), tolerance = 1e-8)
  expect_equal(sort(unname(getME(fit_slash, "theta"))),
               sort(unname(getME(fit_expl, "theta"))),
               tolerance = 1e-8)
  expect_equal(sigma(fit_slash), sigma(fit_expl), tolerance = 1e-8)
})

test_that("numeric || equals its documented explicit expansion (1|g) + (0+x|g)", {
  # Identity: for a NUMERIC slope, (1 + x || g) is documented as exactly
  # (1 | g) + (0 + x | g) (test-compile-audit.R checks this at the
  # compile/kernel level; the factor-|| case is a documented family
  # divergence and is out of scope here). Both spellings parameterize the
  # same diagonal covariance family, so the optima must coincide. Two
  # independent optimizer runs, so 1e-6 instead of 1e-8.
  set.seed(13)
  dd <- data.frame(g = factor(rep(sprintf("h%02d", 1:12), each = 8)))
  dd$x <- round(rnorm(nrow(dd)), 3)
  ui <- rnorm(12, sd = 0.9)
  us <- rnorm(12, sd = 0.5)
  dd$y <- round(1 + 0.8 * dd$x + ui[as.integer(dd$g)] +
                  us[as.integer(dd$g)] * dd$x + rnorm(nrow(dd), sd = 0.4), 4)

  fit_bar <- mm_inv_fit(y ~ x + (1 + x || g), dd)
  fit_expl <- mm_inv_fit(y ~ x + (1 | g) + (0 + x | g), dd)

  expect_equal(as.numeric(logLik(fit_bar)), as.numeric(logLik(fit_expl)),
               tolerance = 1e-6)
  expect_equal(fixef(fit_bar), fixef(fit_expl), tolerance = 1e-6)
  expect_equal(sort(unname(getME(fit_bar, "theta"))),
               sort(unname(getME(fit_expl, "theta"))),
               tolerance = 1e-6)
  expect_equal(sigma(fit_bar), sigma(fit_expl), tolerance = 1e-6)
})

test_that("save/reload leaves newdata predictions byte-identical without revive()", {
  # Identity: JSON artifacts are the source of truth; a plain
  # saveRDS/readRDS round trip (no explicit revive()) must reproduce
  # newdata predictions exactly, both population-level and conditional.
  # (In-sample predict() equality after revive() is already covered by
  # test-phase2-revive.R.)
  df <- mm_inv_data()
  fit <- mm_inv_fit(y ~ x + A + (1 | g), df)
  nd <- df[c(3L, 17L, 41L), c("x", "A", "g")]
  p_pop <- predict(fit, newdata = nd, re.form = NA)
  p_cond <- predict(fit, newdata = nd)

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(fit, tf)
  restored <- readRDS(tf)

  expect_identical(predict(restored, newdata = nd, re.form = NA), p_pop)
  expect_identical(predict(restored, newdata = nd), p_cond)
  expect_identical(predict(restored), predict(fit))
})
