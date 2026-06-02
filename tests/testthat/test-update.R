# Tests for update.mm_lmm / update.mm_glmm (lme4-style model re-fitting).

mm_update_lmm_data <- function() {
  set.seed(101)
  ng <- 12L
  per <- 8L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rep(seq(-1.5, 1.5, length.out = per), ng)
  z <- rnorm(n)
  re_int <- rnorm(ng, sd = 1.2)[as.integer(g)]
  re_slope <- rnorm(ng, sd = 0.6)[as.integer(g)]
  y <- 1 + 0.5 * x - 0.3 * z + re_int + re_slope * x + rnorm(n, sd = 0.7)
  data.frame(y = y, x = x, z = z, g = g)
}

mm_update_binom_data <- function() {
  set.seed(202)
  ng <- 15L
  per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  z <- rnorm(n)
  re <- rnorm(ng, sd = 0.8)[as.integer(g)]
  eta <- -0.2 + 0.7 * x - 0.4 * z + re
  y <- rbinom(n, size = 1, prob = plogis(eta))
  data.frame(y = y, x = x, z = z, g = g)
}

test_that("update.mm_lmm edits the formula and preserves random terms", {
  df <- mm_update_lmm_data()
  fit <- lmm(y ~ x + z + (1 + x | g), df, control = mm_control(verbose = -1))
  fit2 <- update(fit, . ~ . - z)

  expect_s3_class(fit2, "mm_lmm")
  expect_true("x" %in% names(fixef(fit2)))
  expect_false("z" %in% names(fixef(fit2)))
  # the random slope (1 + x | g) survives a fixed-effect-only edit
  expect_named(ranef(fit2), "g")
  expect_true("x" %in% colnames(ranef(fit2)$g))
})

test_that("update.mm_lmm toggles REML and changes the likelihood", {
  df <- mm_update_lmm_data()
  fit <- lmm(y ~ x + (1 | g), df, REML = TRUE, control = mm_control(verbose = -1))
  fit_ml <- update(fit, REML = FALSE)

  expect_false(isTRUE(fit_ml$REML))
  expect_false(isTRUE(all.equal(as.numeric(logLik(fit)),
                                as.numeric(logLik(fit_ml)))))
})

test_that("update.mm_lmm with no changes reproduces the fit", {
  df <- mm_update_lmm_data()
  fit <- lmm(y ~ x + z + (1 + x | g), df, control = mm_control(verbose = -1))
  again <- update(fit)
  expect_equal(unclass(fixef(again)), unclass(fixef(fit)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(again)), as.numeric(logLik(fit)),
               tolerance = 1e-8)
})

test_that("update.mm_lmm accepts new data", {
  df <- mm_update_lmm_data()
  fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  df2 <- df[seq_len(64), , drop = FALSE]
  df2$g <- droplevels(df2$g)
  fit2 <- update(fit, data = df2)
  expect_equal(nobs(fit2), 64L)
})

test_that("update.mm_lmm evaluate = FALSE returns an inspectable call", {
  df <- mm_update_lmm_data()
  fit <- lmm(y ~ x + z + (1 | g), df, control = mm_control(verbose = -1))
  cl <- update(fit, . ~ . - z, evaluate = FALSE)
  expect_true(is.call(cl))
  expect_identical(as.character(cl[[1L]]), "lmm")
  refit <- eval(cl)
  expect_s3_class(refit, "mm_lmm")
  expect_false("z" %in% names(fixef(refit)))
})

test_that("update.mm_glmm refits, reconstructs the family, and edits the formula", {
  df <- mm_update_binom_data()
  fit <- glmm(y ~ x + z + (1 | g), df, family = binomial(),
              control = mm_control(verbose = -1))
  fit2 <- update(fit, . ~ . - z)

  expect_s3_class(fit2, "mm_glmm")
  expect_true("x" %in% names(fixef(fit2)))
  expect_false("z" %in% names(fixef(fit2)))
  # family is reconstructed from stored info when not re-supplied
  expect_identical(fit2$family$family, "binomial")
  expect_identical(fit2$family$link, "logit")
})

test_that("update.mm_glmm with no changes reproduces the fit", {
  df <- mm_update_binom_data()
  fit <- glmm(y ~ x + z + (1 | g), df, family = binomial(),
              control = mm_control(verbose = -1))
  again <- update(fit)
  expect_equal(unclass(fixef(again)), unclass(fixef(fit)), tolerance = 1e-6)
})

test_that("update.mm_lmm matches an lme4 update() formula edit", {
  skip_if_not_installed("lme4")
  df <- mm_update_lmm_data()
  fit <- lmm(y ~ x + z + (1 + x | g), df, REML = TRUE,
             control = mm_control(verbose = -1))
  fit2 <- update(fit, . ~ . - z)

  l1 <- lme4::lmer(y ~ x + z + (1 + x | g), df, REML = TRUE)
  l2 <- update(l1, . ~ . - z)

  expect_equal(sort(names(fixef(fit2))), sort(names(lme4::fixef(l2))))
  expect_equal(unname(fixef(fit2)[names(lme4::fixef(l2))]),
               unname(lme4::fixef(l2)),
               tolerance = 1e-3)
})
