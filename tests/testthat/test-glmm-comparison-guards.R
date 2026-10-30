# These are public fits; the pair guard is also used by the sequential LRT.
test_that("GLMM LRT compares only nested fixed effects on the same likelihood", {
  set.seed(912)
  n <- 120L
  d <- data.frame(g = factor(rep(seq_len(12), each = 10)),
                  x = rnorm(n), z = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.2 + 0.5 * d$x))
  ctl <- mm_control(verbose = -1)
  fit <- function(formula, data = d, family = binomial(), ...) {
    glmm(formula, data, family = family, method = "pirls_profiled",
         control = ctl, ...)
  }
  null <- fit(y ~ 1 + (1 | g))
  x <- fit(y ~ x + (1 | g))
  z <- fit(y ~ z + (1 | g))
  xz <- fit(y ~ x + z + (1 | g))
  zx <- fit(y ~ z + x + (1 | g))
  guard <- mixeff:::mm_assert_glmm_lrt_pair
  expect_true(guard(null, x))
  expect_s3_class(anova(null, x), "anova")
  expect_true(guard(x, zx))
  expect_error(guard(x, z), "nested", class = "mm_inference_unavailable")
  expect_error(guard(xz, zx), "increasing rank", class = "mm_inference_unavailable")

  changed <- d
  changed$y <- 1 - changed$y
  other_y <- fit(y ~ x + (1 | g), changed)
  expect_error(anova(null, other_y), "response values", class = "mm_inference_unavailable")
  expect_error(guard(null, other_y), "response values", class = "mm_inference_unavailable")
  weighted <- fit(y ~ x + (1 | g), weights = rep(2, n))
  expect_error(guard(null, weighted), "prior weights", class = "mm_inference_unavailable")
  shifted <- fit(y ~ x + (1 | g), offset = rep(0.1, n))
  expect_error(guard(null, shifted), "offsets", class = "mm_inference_unavailable")
  poisson_fit <- fit(y ~ x + (1 | g), family = poisson())
  expect_error(guard(null, poisson_fit), "family and link", class = "mm_inference_unavailable")
})

test_that("GLMM LRT checks grouped trials and random-effect designs", {
  set.seed(913)
  n <- 100L
  d <- data.frame(g = factor(rep(seq_len(10), each = 10)), x = rnorm(n))
  d$trials <- rep(5, n)
  d$success <- rbinom(n, d$trials, plogis(0.2 + 0.4 * d$x))
  d$failure <- d$trials - d$success
  ctl <- mm_control(verbose = -1)
  fit <- function(formula, data) glmm(formula, data, family = binomial(),
                                     method = "pirls_profiled", control = ctl)
  null <- fit(cbind(success, failure) ~ 1 + (1 | g), d)
  alternative <- fit(cbind(success, failure) ~ x + (1 | g), d)
  guard <- mixeff:::mm_assert_glmm_lrt_pair
  expect_true(guard(null, alternative))
  changed <- d
  # Keep the observed proportion unchanged while changing trial counts.
  changed$success <- 2 * d$success
  changed$failure <- 2 * d$failure
  trials <- fit(cbind(success, failure) ~ x + (1 | g), changed)
  expect_error(guard(null, trials), "trial counts", class = "mm_inference_unavailable")
  changed <- d
  changed$g <- factor(rep(seq_len(10), times = 10))
  groups <- fit(cbind(success, failure) ~ x + (1 | g), changed)
  expect_error(guard(null, groups), "random-effects design", class = "mm_inference_unavailable")
  random_slope <- fit(cbind(success, failure) ~ x + (1 + x | g), d)
  expect_error(guard(null, random_slope), "random-effects structure", class = "mm_inference_unavailable")
})


test_that("free-dispersion LRTs require the same recorded estimation policy", {
  set.seed(777)
  n <- 160L
  d <- data.frame(g = factor(rep(seq_len(16), each = 10)), x = rnorm(n))
  group_effect <- rep(rnorm(16, sd = 0.3), each = 10)
  d$y <- rgamma(n, shape = 2, scale = exp(0.2 + 0.3 * d$x + group_effect) / 2)
  fit <- function(formula, ...) {
    glmm(formula, d, family = Gamma(link = "log"), method = "pirls_profiled",
         control = mm_control(verbose = -1, ...))
  }
  null <- fit(y ~ 1 + (1 | g))
  alternative <- fit(y ~ x + (1 | g), disp_method = "moment")
  legacy <- fit(y ~ x + (1 | g), disp_method = "old/buggy")
  uncorrected <- fit(y ~ x + (1 | g), disp_dof_correction = FALSE)
  expect_s3_class(anova(null, alternative), "anova")
  expect_error(anova(null, legacy), "dispersion policy", class = "mm_inference_unavailable")
  expect_error(anova(null, uncorrected), "dispersion policy", class = "mm_inference_unavailable")
  # Editing controls cannot relabel the estimator of a recorded solution.
  legacy$control$disp_method <- "moment"
  expect_error(anova(null, legacy), "dispersion policy", class = "mm_inference_unavailable")
})


test_that("negative-binomial LRTs preserve fixed or estimated theta semantics", {
  set.seed(914)
  d <- data.frame(g = factor(rep(seq_len(12), each = 12)), x = rnorm(144))
  b <- rep(rnorm(12, sd = 0.4), each = 12)
  d$y <- rnbinom(nrow(d), size = 1.8, mu = exp(0.5 + 0.5 * d$x + b))
  fit <- function(formula, theta = NULL) {
    glmm(formula, d, family = mm_negative_binomial(theta),
         method = "pirls_profiled", control = mm_control(verbose = -1))
  }
  fixed_null <- fit(y ~ 1 + (1 | g), 1.8)
  fixed_alt <- fit(y ~ x + (1 | g), 1.8)
  changed_theta <- fit(y ~ x + (1 | g), 3)
  estimated_null <- fit(y ~ 1 + (1 | g))
  estimated_alt <- fit(y ~ x + (1 | g))

  for (pair in list(list(fixed_null, fixed_alt),
                    list(estimated_null, estimated_alt))) {
    comparison <- anova(pair[[1]], pair[[2]])
    expect_s3_class(comparison, "anova")
    expect_equal(comparison$Df[2], 1)
    expected <- max(0, 2 * (as.numeric(logLik(pair[[2]])) -
                            as.numeric(logLik(pair[[1]]))))
    expect_equal(comparison$Chisq[2], expected)
    expect_equal(comparison[["Pr(>Chisq)"]][2],
                 pchisq(expected, df = 1, lower.tail = FALSE))
  }
  expect_true(estimated_null$family$nb_theta_estimated)
  expect_true(estimated_alt$family$nb_theta_estimated)
  # Estimated theta is a nuisance parameter refit in each model; its fitted
  # value need not agree. A different fixed theta changes the likelihood.
  expect_gt(abs(estimated_null$family$nb_theta - estimated_alt$family$nb_theta),
            1e-3)
  for (alternative in list(changed_theta, estimated_alt)) {
    err <- tryCatch(anova(fixed_null, alternative), error = identity)
    expect_s3_class(err, "mm_inference_unavailable")
    expect_identical(err$reason_code, "glmm_lrt_nb_theta_differs")
  }
})

test_that("GLMM LRT refuses mixing effective estimators", {
  set.seed(915)
  d <- data.frame(g = factor(rep(seq_len(12), each = 12)), x = rnorm(144))
  b <- rep(rnorm(12, sd = 0.6), each = 12)
  d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.5 * d$x + b))
  null <- glmm(y ~ 1 + (1 | g), d, family = binomial(),
               method = "pirls_profiled", control = mm_control(verbose = -1))
  alternative <- glmm(y ~ x + (1 | g), d, family = binomial(),
                      method = "joint_laplace", control = mm_control(verbose = -1))
  err <- tryCatch(anova(null, alternative), error = identity)
  expect_s3_class(err, "mm_inference_unavailable")
  expect_identical(err$reason_code, "glmm_lrt_estimator_differs")
})
