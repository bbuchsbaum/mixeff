# PR-4 coverage: Nakagawa variance decomposition edge cases (link-specific
# distribution variances, no-intercept designs, slopes absent from X) and the
# insight-facing helpers.

cp4p_ctrl <- function() mm_control(verbose = -1)

cp4p_binomial <- function(seed = 431L) {
  set.seed(seed)
  g <- gl(12, 6)
  x <- rnorm(72)
  size <- rep(c(4L, 6L, 8L), 24)
  k <- rbinom(72, size, plogis(-0.2 + 0.5 * x + rnorm(12, sd = 0.5)[g]))
  data.frame(k = k, n = size, x = x, g = g)
}

test_that("binomial distribution variance follows the link and trials", {
  d <- cp4p_binomial()
  fit <- glmm(cbind(k, n - k) ~ x + (1 | g), d, family = binomial(),
              control = cp4p_ctrl())
  dv <- mixeff:::mm_distribution_variance
  # cbind() responses divide by the mean number of trials.
  expect_equal(dv(fit), (pi^2 / 3) / mean(d$n))
  v <- mm_variance_components(fit)
  expect_equal(v$var.distribution, (pi^2 / 3) / mean(d$n))
  probit <- fit
  probit$family$link <- "probit"
  expect_equal(dv(probit), 1 / mean(d$n))
  cll <- fit
  cll$family$link <- "cloglog"
  expect_equal(dv(cll), (pi^2 / 6) / mean(d$n))
  cauchit <- fit
  cauchit$family$link <- "cauchit"
  expect_error(dv(cauchit), class = "mm_inference_unavailable")

  # A 0/1 response has one trial per row.
  d$b <- as.integer(d$k > d$n / 2)
  fb <- glmm(b ~ x + (1 | g), d, family = binomial(), control = cp4p_ctrl())
  expect_identical(mixeff:::mm_binomial_trial_factor(fb), 1)
  expect_equal(dv(fb), pi^2 / 3)
})

test_that("Poisson distribution variance handles sqrt, bad links and low mu", {
  set.seed(432)
  g <- gl(10, 8)
  x <- rnorm(80)
  d <- data.frame(y = rpois(80, exp(0.1 + 0.2 * x + rnorm(10, sd = 0.3)[g])),
                  x = x, g = g)
  fit <- glmm(y ~ x + (1 | g), d, family = poisson(), control = cp4p_ctrl())
  dv <- mixeff:::mm_distribution_variance
  sq <- fit
  sq$family$link <- "sqrt"
  expect_identical(dv(sq), 0.25)
  idl <- fit
  idl$family$link <- "identity"
  expect_error(dv(idl), class = "mm_inference_unavailable")
  ig <- fit
  ig$family$family <- "inverse_gaussian"
  expect_error(dv(ig), class = "mm_inference_unavailable")
  # Counts this small trigger the lognormal-approximation warning ...
  expect_warning(v <- dv(fit), "too close to zero")
  expect_gt(v, 0)
  # ... which verbose = FALSE silences.
  expect_no_warning(dv(fit, warn = FALSE))
})

test_that("variance parts handle no-intercept designs and slopes outside X", {
  set.seed(433)
  g <- gl(10, 6)
  x <- rnorm(60)
  w <- rnorm(60)
  y <- 2 + x + rnorm(10)[g] + rnorm(10, sd = 0.5)[g] * w + rnorm(60, sd = 0.5)
  d <- data.frame(y = y, x = x, w = w, g = g)
  f0 <- lmm(y ~ 0 + x + (1 | g), d, control = cp4p_ctrl())
  v0 <- mm_variance_components(f0)
  vc <- VarCorr(f0)
  expect_equal(v0$var.random, as.numeric(attr(vc$g, "stddev"))^2,
               tolerance = 1e-8)
  expect_equal(v0$var.fixed, var(d$x * fixef(f0)[["x"]]))
  # A random slope on a variable that is not a fixed effect contributes no
  # variance through X (as insight does).
  fs <- lmm(y ~ x + (0 + w | g), d, control = cp4p_ctrl())
  vs <- mm_variance_components(fs)
  expect_identical(vs$var.random, 0)
  expect_identical(names(vs$var.slope), "g.w")
  skip_if_not_installed("lme4")
  skip_if_not_installed("insight")
  ref <- suppressWarnings(insight::get_variance(lme4::lmer(y ~ 0 + x + (1 | g), d)))
  expect_equal(v0$var.random, ref$var.random, tolerance = 1e-4)
})

test_that("an empty random-term basis is read as an intercept", {
  set.seed(434)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = gl(10, 6))
  d$y <- d$y + rnorm(10)[d$g]
  fit <- lmm(y ~ x + (1 | g), d, control = cp4p_ctrl())
  ref <- mm_variance_components(fit)
  nb <- fit
  nb$artifact$semantic_model$random_terms[[1]]$basis <- list()
  expect_equal(mm_variance_components(nb), ref)
})

test_that("insight helpers: get_variance(), find_random_slopes()", {
  set.seed(435)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = gl(10, 6))
  d$y <- d$y + rnorm(10)[d$g] + rnorm(10, sd = 0.4)[d$g] * d$x
  fit <- lmm(y ~ x + (x | g), d, control = cp4p_ctrl())
  gv <- mixeff:::get_variance.mm_fit
  all <- gv(fit)
  expect_equal(gv(fit, "random"), all$var.random)
  expect_equal(gv(fit, "fixed"), all$var.fixed)
  expect_equal(gv(fit, "slope"), all$var.slope)
  broken <- fit
  broken$artifact$semantic_model$random_terms[[1]]$basis <- list(list(bad = TRUE))
  broken$model_frame <- NULL
  expect_warning(out <- gv(broken), ".")
  expect_identical(out, NA)
  expect_silent(out2 <- gv(broken, verbose = FALSE))
  expect_identical(out2, NA)
  frs <- mixeff:::find_random_slopes.mm_fit
  expect_identical(frs(fit), list(random = "x"))
  ri <- lmm(y ~ x + (1 | g), d, control = cp4p_ctrl())
  expect_null(frs(ri))
  expect_true(mixeff:::is_mixed_model.mm_fit(ri))
})
