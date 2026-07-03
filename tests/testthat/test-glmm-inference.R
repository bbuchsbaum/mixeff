# Tests for GLMM fixed-effect inference: contrast / drop1 / anova.

mm_glmm_inf_data <- function() {
  set.seed(707)
  ng <- 18L
  per <- 14L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  z <- rnorm(n)
  re <- rnorm(ng, sd = 0.7)[as.integer(g)]
  y <- rbinom(n, 1, plogis(-0.3 + 0.8 * x - 0.5 * z + re))
  data.frame(y = y, x = x, z = z, g = g)
}

test_that("contrast.mm_glmm gives Wald z inference matching the coefficient", {
  # joint_laplace is the glmer-equivalent estimator whose fixed-effect Wald
  # inference is certified by the engine; contrast consumes that certified
  # covariance, so SE/z/p are available and match fit$std_errors.
  fit <- glmm(y ~ x + z + (1 | g), mm_glmm_inf_data(), family = binomial(),
              method = "joint_laplace", control = mm_control(verbose = -1))
  # contrast selecting the 'x' coefficient
  L <- c(0, 1, 0)
  ct <- contrast(fit, L)
  expect_s3_class(ct, "mm_contrast")
  expect_equal(ct$table$estimate, unname(fixef(fit)["x"]), tolerance = 1e-8)
  expect_equal(ct$table$std_error, unname(fit$std_errors["x"]), tolerance = 1e-8)
  expect_identical(ct$table$statistic_name, "z")
  expect_equal(ct$table$statistic, unname(fixef(fit)["x"] / fit$std_errors["x"]),
               tolerance = 1e-8)
  expect_true(ct$table$p_value >= 0 && ct$table$p_value <= 1)
})

test_that("contrast.mm_glmm withholds Wald inference for uncertified PIRLS fits", {
  # The default fast-PIRLS path is not certified for fixed-effect Wald
  # inference, so contrast reports the estimate but withholds SE/z/p with a
  # reason rather than fabricating numbers from the uncertified working Hessian
  # (no fake certainty), matching summary()/confint()/tidy().
  fit <- glmm(y ~ x + z + (1 | g), mm_glmm_inf_data(), family = binomial(),
              method = "pirls_profiled", control = mm_control(verbose = -1))
  ct <- contrast(fit, c(0, 1, 0))
  expect_equal(ct$table$estimate, unname(fixef(fit)["x"]), tolerance = 1e-8)
  expect_true(is.na(ct$table$std_error))
  expect_true(is.na(ct$table$statistic))
  expect_true(is.na(ct$table$p_value))
  expect_identical(ct$table$method, "not_computed")
  expect_true(is.character(ct$table$reason) && nzchar(ct$table$reason))
})

test_that("contrast.mm_glmm honours rhs and the wald alias", {
  fit <- glmm(y ~ x + z + (1 | g), mm_glmm_inf_data(), family = binomial(),
              method = "joint_laplace", control = mm_control(verbose = -1))
  a <- contrast(fit, c(0, 1, 0), rhs = 0.5, method = "wald")
  expect_equal(a$table$estimate, unname(fixef(fit)["x"]) - 0.5, tolerance = 1e-8)
  expect_identical(a$table$method, "asymptotic")
})

test_that("drop1.mm_glmm runs likelihood-ratio term deletions", {
  fit <- glmm(y ~ x + z + (1 | g), mm_glmm_inf_data(), family = binomial(),
              control = mm_control(verbose = -1))
  d <- drop1(fit, test = "Chisq")
  expect_s3_class(d, "mm_drop1")
  expect_setequal(d$table$dropped, c("x", "z"))
  expect_true(all(d$table$df == 1))
  expect_true(all(d$table$LRT >= 0))
  expect_true(all(d$table$p_value >= 0 & d$table$p_value <= 1))
  expect_true(all(d$table$method == "asymptotic_lrt"))
})

test_that("drop1.mm_glmm test='none' reports criteria without LRT", {
  fit <- glmm(y ~ x + z + (1 | g), mm_glmm_inf_data(), family = binomial(),
              control = mm_control(verbose = -1))
  d <- drop1(fit)
  expect_true(all(is.na(d$table$LRT)))
  expect_true(all(is.finite(d$table$AIC)))
})

test_that("anova.mm_glmm compares nested models by sequential LRT", {
  df <- mm_glmm_inf_data()
  full <- glmm(y ~ x + z + (1 | g), df, family = binomial(),
               control = mm_control(verbose = -1))
  reduced <- update(full, . ~ . - z)
  cmp <- anova(reduced, full)
  expect_s3_class(cmp, "mm_glmm_comparison")
  expect_equal(nrow(cmp$table), 2L)
  # rows ordered by npar ascending: reduced first, full second
  expect_equal(cmp$table$npar, sort(cmp$table$npar))
  expect_true(is.na(cmp$table$Chisq[1]))
  expect_true(is.finite(cmp$table$Chisq[2]) && cmp$table$Chisq[2] >= 0)
  expect_equal(cmp$table$Df[2], 1L)
  expect_true(cmp$table$p_value[2] >= 0 && cmp$table$p_value[2] <= 1)
})

test_that("anova.mm_glmm refuses single-model analysis with guidance", {
  fit <- glmm(y ~ x + z + (1 | g), mm_glmm_inf_data(), family = binomial(),
              control = mm_control(verbose = -1))
  expect_error(anova(fit), class = "mm_inference_unavailable")
})
