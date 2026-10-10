# PR-4 coverage: lme4-compatibility helpers (family(), weights(), REMLcrit(),
# the theta-map random-effect structure, getME() extras, VarCorr formatting,
# multi-model information tables) and the S3 hook registration in zzz.R.
# Small local fixtures; values are checked against lme4 where natural.

cp4_data <- function(seed = 401L) {
  set.seed(seed)
  g <- gl(10, 6)
  h <- factor(rep(letters[1:6], 10))
  x <- rnorm(60)
  f <- factor(rep(c("p", "q", "r"), 20))
  y <- 1 + 0.5 * x + rnorm(10, sd = 0.8)[g] + rnorm(6, sd = 0.5)[h] +
    rnorm(60, sd = 0.6)
  data.frame(y = y, x = x, g = g, h = h, f = f)
}

cp4_ctrl <- function() mm_control(verbose = -1)

test_that("zzz: mm_register_external_s3() installs package hooks and methods", {
  events <- c("lme4", "lmerTest", "emmeans", "generics", "lattice",
              "insight", "performance")
  # Load the optional namespaces first so the "already loaded" branches run.
  for (p in c("lme4", "lattice", "insight", "performance", "lmerTest",
              "emmeans", "generics")) {
    suppressMessages(requireNamespace(p, quietly = TRUE))
  }
  hook_names <- vapply(events, function(p) packageEvent(p, "onLoad"),
                       character(1))
  saved <- lapply(hook_names, getHook)
  on.exit({
    for (i in seq_along(hook_names)) {
      setHook(hook_names[[i]], NULL, "replace")
      for (h in saved[[i]]) setHook(hook_names[[i]], h, "append")
    }
  }, add = TRUE)
  mixeff:::mm_register_external_s3()
  after <- lapply(hook_names, getHook)
  expect_true(all(lengths(after) > lengths(saved)))
  # Each event gained a hook that re-runs mixeff's registration.
  ours <- lapply(after, function(hooks) {
    Filter(function(h) {
      is.function(h) && any(grepl("mm_register_", deparse(body(h))))
    }, hooks)
  })
  expect_true(all(lengths(ours) >= 1L))
  if (isNamespaceLoaded("lme4")) {
    expect_identical(
      getS3method("fixef", "mm_lmm", envir = asNamespace("lme4")),
      mixeff:::fixef.mm_lmm
    )
  }
  if (isNamespaceLoaded("lattice")) {
    expect_true(is.function(
      getS3method("dotplot", "mm_ranef", envir = asNamespace("lattice"))
    ))
  }
  # Firing the new hooks directly is harmless (idempotent re-registration).
  for (i in seq_along(ours)) {
    if (isNamespaceLoaded(events[[i]])) {
      expect_no_error(ours[[i]][[length(ours[[i]])]]())
    }
  }
})

test_that("family() maps every engine family to an R family object", {
  skip_if_not_installed("MASS")
  set.seed(402)
  g <- gl(8, 8)
  x <- rnorm(64)
  mu <- exp(0.3 + 0.2 * x + rnorm(8, sd = 0.3)[g])
  d <- data.frame(x = x, g = g, yg = rgamma(64, shape = 4, rate = 4 / mu),
                  yc = rnbinom(64, mu = mu, size = 2))
  gam <- glmm(yg ~ x + (1 | g), d, family = Gamma(link = "log"),
              control = cp4_ctrl())
  fam <- family(gam)
  expect_s3_class(fam, "family")
  expect_identical(fam$family, "Gamma")
  expect_identical(fam$link, "log")
  # lme4's capitalised family name is accepted as well as the engine's.
  gm <- gam
  gm$family$family <- "Gamma"
  expect_identical(family(gm)$family, "Gamma")

  nb <- glmm(yc ~ x + (1 | g), d, family = mm_negative_binomial(2),
             control = cp4_ctrl())
  fnb <- family(nb)
  expect_match(fnb$family, "Negative Binomial")
  ref <- MASS::negative.binomial(theta = 2)
  expect_equal(fnb$variance(c(1, 3)), ref$variance(c(1, 3)))

  # A Bernoulli-coded fit is still a binomial family in R.
  b <- nb
  b$family <- list(family = "bernoulli", link = "logit")
  expect_identical(family(b)$family, "binomial")
  # Unknown family names and missing NB theta are typed refusals.
  bad <- nb
  bad$family <- list(family = "tweedie", link = "log")
  expect_error(family(bad), class = "mm_inference_unavailable")
  expect_error(mixeff:::mm_nb_family(NA_real_),
               class = "mm_inference_unavailable")
  expect_error(mixeff:::mm_nb_family(numeric()),
               class = "mm_inference_unavailable")
})

test_that("the MASS-free negative-binomial family matches MASS", {
  skip_if_not_installed("MASS")
  # Re-evaluate mm_nb_family() where requireNamespace() reports MASS absent.
  fn <- mixeff:::mm_nb_family
  environment(fn) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = asNamespace("mixeff")
  )
  own <- fn(1.7)
  ref <- MASS::negative.binomial(theta = 1.7)
  expect_s3_class(own, "family")
  expect_identical(own$link, "log")
  expect_identical(own$theta, 1.7)
  expect_match(own$family, "Negative Binomial\\(1\\.7\\)")
  y <- c(0, 1, 4, 9)
  mu <- c(0.5, 1.2, 3.5, 7)
  wt <- c(1, 2, 1, 1)
  expect_equal(own$variance(mu), ref$variance(mu))
  expect_equal(own$dev.resids(y, mu, wt), ref$dev.resids(y, mu, wt))
  expect_equal(own$linkinv(own$linkfun(mu)), mu)
  expect_equal(own$mu.eta(log(mu)), ref$mu.eta(log(mu)))
  expect_true(own$validmu(mu))
  expect_false(own$validmu(c(1, 0)))
  expect_true(own$valideta(log(mu)))
  # MASS's aic() is -2 loglik (up to a constant in theta); ours is the full
  # -2 log-likelihood of the NB2 density.
  expect_equal(own$aic(y, length(y), mu, wt, 0),
               -2 * sum(wt * dnbinom(y, size = 1.7, mu = mu, log = TRUE)))
})

test_that("weights(), REMLcrit() and offsets follow lme4", {
  skip_if_not_installed("lme4")
  d <- cp4_data()
  ml <- lmm(y ~ x + (1 | g), d, REML = FALSE, control = cp4_ctrl())
  reml <- lmm(y ~ x + (1 | g), d, REML = TRUE, control = cp4_ctrl())
  ref <- lme4::lmer(y ~ x + (1 | g), d, REML = FALSE)
  expect_equal(REMLcrit(ml), lme4::REMLcrit(ref), tolerance = 1e-5)
  expect_equal(REMLcrit(reml), -2 * as.numeric(logLik(reml)))
  # merMod objects are forwarded to lme4.
  expect_equal(REMLcrit(ref), lme4::REMLcrit(ref))
  expect_error(REMLcrit(list()), class = "mm_arg_error")
  expect_equal(weights(ml), rep(1, 60))
  expect_error(weights(ml, type = "working"), class = "mm_arg_error")
  expect_identical(getME(ml, "REML"), 0L)
  expect_identical(getME(reml, "REML"), 2L)
  expect_identical(getME(ml, "REML"), lme4::getME(ref, "REML"))
  expect_equal(getME(ml, "offset"), rep(0, 60))
  expect_equal(getME(ml, "weights"), rep(1, 60))

  set.seed(403)
  gd <- data.frame(g = gl(8, 8), x = rnorm(64), off = log(runif(64, 1, 3)))
  gd$y <- rpois(64, exp(0.2 + 0.3 * gd$x + gd$off))
  gf <- glmm(y ~ x + offset(off) + (1 | g), gd, family = poisson(),
             control = cp4_ctrl())
  expect_error(REMLcrit(gf), class = "mm_inference_unavailable")
  expect_equal(getME(gf, "offset"), gd$off)
  expect_identical(getME(gf, "REML"), 0L)
  expect_equal(getME(gf, "weights"), rep(1, 64))
  gref <- lme4::glmer(y ~ x + offset(off) + (1 | g), gd, family = poisson())
  expect_equal(getME(gf, "offset"), lme4::getME(gref, "offset"))
  ww <- weights(gf, type = "working")
  expect_equal(ww, as.numeric(fitted(gf)), tolerance = 1e-8)
  expect_equal(ww, as.numeric(weights(gref, type = "working")),
               tolerance = 1e-3)
  # Without an offset term the GLMM offset is zero.
  gf0 <- glmm(y ~ x + (1 | g), gd, family = poisson(), control = cp4_ctrl())
  expect_equal(getME(gf0, "offset"), rep(0, 64))
})

test_that("the theta-map structure rejects inconsistent maps", {
  d <- cp4_data()
  fit <- lmm(y ~ x + (1 + x || g), d, control = cp4_ctrl())
  st <- mixeff:::mm_compute_re_structure(fit)
  expect_identical(st$perm, 1:2)
  expect_identical(vapply(st$terms, `[[`, "", "group"), c("g", "g"))

  none <- fit
  none$artifact$theta_maps <- list()
  expect_error(mixeff:::mm_compute_re_structure(none),
               class = "mm_inference_unavailable")
  # A singular-check on such a fit falls back to the fit status.
  expect_identical(mixeff:::mm_is_singular_theta(none, 1e-4), NA)
  expect_type(is_singular(none), "logical")

  far <- fit
  far$artifact$theta_maps[[1]]$map$theta_slots[[1]]$lambda_row <- 5L
  expect_error(mixeff:::mm_compute_re_structure(far),
               class = "mm_schema_error")

  out_of_range <- fit
  out_of_range$artifact$theta_maps[[2]]$map$theta_slots[[1]]$global_index <- 9L
  expect_error(mixeff:::mm_compute_re_structure(out_of_range),
               class = "mm_schema_error")

  dup <- fit
  dup$artifact$theta_maps[[2]]$map$theta_slots[[1]]$global_index <- 0L
  expect_error(mixeff:::mm_compute_re_structure(dup),
               class = "mm_schema_error")

  # A map without any basis is read as an intercept term.
  bare <- lmm(y ~ x + (1 | g), d, control = cp4_ctrl())
  nb <- bare
  nb$artifact$theta_maps[[1]]$map$user_basis <- NULL
  nb$artifact$theta_maps[[1]]$map$optimizer_basis <- NULL
  st2 <- mixeff:::mm_compute_re_structure(nb)
  expect_identical(st2$terms[[1]]$cnames, "(Intercept)")

  # Conditional modes that do not match the map are a schema error.
  wrong <- bare
  wrong$lazy_cache <- new.env(parent = emptyenv())
  names(wrong$random_effects$g) <- "zzz"
  expect_error(mixeff:::mm_re_bu(wrong), class = "mm_schema_error")
})

test_that("factor random slopes rebuild lme4's Z", {
  skip_if_not_installed("lme4")
  d <- cp4_data()
  fit <- lmm(y ~ x + (0 + f | g), d, control = cp4_ctrl())
  ref <- suppressMessages(lme4::lmer(y ~ x + (0 + f | g), d))
  expect_equal(as.matrix(getME(fit, "Z")), as.matrix(lme4::getME(ref, "Z")),
               ignore_attr = TRUE)
  expect_identical(getME(fit, "cnms"), lme4::getME(ref, "cnms"))
  frame <- d
  expect_null(mixeff:::mm_re_factor_basis_matrix("intercept", TRUE, frame))
  expect_null(mixeff:::mm_re_factor_basis_matrix("a b (", FALSE, frame))
  expect_null(mixeff:::mm_re_factor_basis_matrix("nope", FALSE, frame))
  mm <- mixeff:::mm_re_factor_basis_matrix(c("intercept", "f"),
                                           c(TRUE, FALSE), frame)
  expect_identical(colnames(mm), c("(Intercept)", "fq", "fr"))
})

test_that("random-effect basis columns are rebuilt from labels", {
  frame <- data.frame(a = c(1, 2, 3), b = c(2, 2, 0.5),
                      f = factor(c("u", "v", "u")), s = c("m", "n", "m"))
  bv <- mixeff:::mm_re_basis_values
  expect_equal(bv("(Intercept)", frame), c(1, 1, 1))
  expect_equal(bv("a", frame), c(1, 2, 3))
  expect_equal(bv("f: v", frame), c(0, 1, 0))
  expect_equal(bv("a:b", frame), c(2, 4, 1.5))
  expect_error(bv("s:a", frame), class = "mm_inference_unavailable")
  expect_error(bv("missing", frame), class = "mm_inference_unavailable")
})

test_that("GLMM condVar merges a split `||` term per grouping factor", {
  set.seed(404)
  g <- gl(12, 8)
  x <- rnorm(96)
  y <- rpois(96, exp(0.3 + 0.3 * x + rnorm(12, sd = 0.4)[g] +
                     rnorm(12, sd = 0.4)[g] * x))
  d <- data.frame(y = y, x = x, g = g)
  fit <- glmm(y ~ x + (1 + x || g), d, family = poisson(),
              control = cp4_ctrl())
  re <- ranef(fit, condVar = TRUE)
  pv <- attr(re$g, "postVar")
  expect_identical(dim(pv), c(2L, 2L, 12L))
  # Independent pieces: the merged block is diagonal and positive.
  expect_true(all(pv[1, 2, ] == 0))
  expect_true(all(pv[1, 1, ] > 0 & pv[2, 2, ] > 0))
  skip_if_not_installed("lme4")
  ref <- lme4::glmer(y ~ x + (1 + x || g), d, family = poisson())
  rpv <- attr(lme4::ranef(ref, condVar = TRUE)$g, "postVar")
  # lme4 1.1 keeps one array per `||` piece; lme4 2.x merges them.
  r11 <- if (is.list(rpv)) rpv[[1]][1, 1, ] else rpv[1, 1, ]
  r22 <- if (is.list(rpv)) rpv[[2]][1, 1, ] else rpv[2, 2, ]
  expect_equal(unname(pv[1, 1, ]), r11, tolerance = 0.05)
  expect_equal(unname(pv[2, 2, ]), r22, tolerance = 0.05)
})

test_that("VarCorr formatting pads mixed-size terms and handles none", {
  skip_if_not_installed("lme4")
  d <- cp4_data()
  fit <- lmm(y ~ x + (1 + x | g) + (1 | h), d, control = cp4_ctrl())
  vc <- VarCorr(fit)
  fm <- mixeff:::mm_format_varcorr(vc)
  expect_identical(colnames(fm)[1:3], c("Groups", "Name", "Std.Dev."))
  expect_identical(nrow(fm), 4L)
  expect_identical(unname(fm[, "Groups"]), c("g", "", "h", "Residual"))
  expect_identical(unname(fm[3:4, "Corr"]), c("", ""))
  expect_output(print(vc), "Corr")
  empty <- structure(list(), useSc = FALSE, class = c("mm_varcorr", "VarCorr.merMod"))
  fe <- mixeff:::mm_format_varcorr(empty)
  expect_identical(nrow(fe), 0L)
  expect_output(print(empty), "Random effects: none")
  expect_null(mixeff:::mm_varcorr_lme4_terms(list(), character(), fit))
  # Engine blocks that do not match the lme4 terms keep the engine layout.
  odd <- list(matrix(1, 1, 1, dimnames = list("w", "w")))
  expect_null(mixeff:::mm_varcorr_lme4_terms(odd, "g", fit))
  # ... as do blocks with engine columns no lme4 term accounts for.
  extra <- c(unclass(vc)[1:2], odd)
  expect_null(mixeff:::mm_varcorr_lme4_terms(extra, c("g", "h", "g"), fit))
  expect_length(mixeff:::mm_varcorr_lme4_terms(unclass(vc)[1:2], c("g", "h"), fit), 2L)
  expect_null(vc$not_a_component)
  expect_s3_class(vc$table, "data.frame")
})

test_that("AIC()/BIC() tables warn when the models use different rows", {
  d <- cp4_data()
  a <- lmm(y ~ x + (1 | g), d, REML = FALSE, control = cp4_ctrl())
  b <- lmm(y ~ x + (1 | g), d[-(1:6), ], REML = FALSE, control = cp4_ctrl())
  expect_warning(tab <- AIC(a, b), "same number of observations")
  expect_identical(names(tab), c("df", "AIC"))
  expect_equal(tab$AIC[[1]], AIC(a))
  expect_warning(tb <- BIC(a, b), "same number of observations")
  expect_equal(tb$BIC[[2]], BIC(b))
})
