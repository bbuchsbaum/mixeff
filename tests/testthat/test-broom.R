# Tests for broom / broom.mixed tidy()/glance()/augment() support.

mm_broom_lmm_data <- function() {
  set.seed(303)
  ng <- 12L; per <- 8L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rep(seq(-1.5, 1.5, length.out = per), ng)
  re_int <- rnorm(ng, sd = 1.1)[as.integer(g)]
  re_slope <- rnorm(ng, sd = 0.5)[as.integer(g)]
  y <- 1 + 0.5 * x + re_int + re_slope * x + rnorm(n, sd = 0.7)
  data.frame(y = y, x = x, g = g)
}

mm_broom_binom_data <- function() {
  set.seed(404)
  ng <- 15L; per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.8)[as.integer(g)]
  y <- rbinom(n, 1, plogis(-0.2 + 0.7 * x + re))
  data.frame(y = y, x = x, g = g)
}

test_that("tidy.mm_lmm returns fixed and ran_pars rows in broom.mixed shape", {
  skip_if_not_installed("generics")
  fit <- lmm(y ~ x + (1 + x | g), mm_broom_lmm_data(),
             control = mm_control(verbose = -1))
  td <- generics::tidy(fit)

  expect_true(all(c("effect", "group", "term", "estimate") %in% names(td)))
  expect_setequal(unique(td$effect), c("fixed", "ran_pars"))
  # fixed effects
  fixed <- td[td$effect == "fixed", ]
  expect_setequal(fixed$term, c("(Intercept)", "x"))
  expect_equal(fixed$estimate, unname(fixef(fit)), tolerance = 1e-8)
  expect_true(all(is.finite(fixed$std.error)))
  # ran_pars: sd__ terms, a cor__ term, and the residual sd__Observation
  rp <- td[td$effect == "ran_pars", ]
  expect_true(any(grepl("^sd__", rp$term)))
  expect_true(any(grepl("^cor__", rp$term)))
  expect_true("sd__Observation" %in% rp$term)
})

test_that("tidy dispatches through broom and broom.mixed too", {
  skip_if_not_installed("broom")
  skip_if_not_installed("broom.mixed")
  fit <- lmm(y ~ x + (1 | g), mm_broom_lmm_data(),
             control = mm_control(verbose = -1))
  a <- broom::tidy(fit)
  b <- broom.mixed::tidy(fit)
  expect_true(is.data.frame(a))
  expect_identical(a$term, b$term)
  expect_equal(a$estimate, b$estimate)
})

test_that("tidy conf.int adds Wald intervals for fixed effects", {
  skip_if_not_installed("generics")
  fit <- lmm(y ~ x + (1 | g), mm_broom_lmm_data(),
             control = mm_control(verbose = -1))
  td <- generics::tidy(fit, effects = "fixed", conf.int = TRUE)
  expect_true(all(c("conf.low", "conf.high") %in% names(td)))
  expect_true(all(td$conf.low < td$estimate & td$estimate < td$conf.high))
})

test_that("glance.mm_lmm returns a one-row model summary", {
  skip_if_not_installed("generics")
  fit <- lmm(y ~ x + (1 | g), mm_broom_lmm_data(),
             control = mm_control(verbose = -1))
  gl <- generics::glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("nobs", "sigma", "logLik", "AIC", "BIC", "deviance",
                    "df.residual") %in% names(gl)))
  expect_equal(gl$nobs, nobs(fit))
  expect_equal(gl$AIC, AIC(fit), tolerance = 1e-8)
})

test_that("augment.mm_lmm adds .fitted and .resid", {
  skip_if_not_installed("generics")
  df <- mm_broom_lmm_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  ag <- generics::augment(fit)
  expect_true(all(c(".fitted", ".resid") %in% names(ag)))
  expect_equal(nrow(ag), nobs(fit))
  expect_equal(ag$.fitted, unname(fitted(fit)), tolerance = 1e-8)
})

test_that("tidy.mm_glmm returns Wald p-values for fixed effects", {
  skip_if_not_installed("generics")
  fit <- glmm(y ~ x + (1 | g), mm_broom_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  td <- generics::tidy(fit, effects = "fixed")
  expect_true("p.value" %in% names(td))
  expect_true(all(td$p.value >= 0 & td$p.value <= 1))
  expect_equal(td$estimate, unname(fixef(fit)), tolerance = 1e-8)
})

test_that("glance.mm_glmm works", {
  skip_if_not_installed("generics")
  fit <- glmm(y ~ x + (1 | g), mm_broom_binom_data(), family = binomial(),
              control = mm_control(verbose = -1))
  gl <- generics::glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$nobs, nobs(fit))
})
