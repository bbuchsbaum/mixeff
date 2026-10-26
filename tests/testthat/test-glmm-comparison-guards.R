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
  expect_true(guard(x, zx))
  expect_error(guard(x, z), "nested", class = "mm_inference_unavailable")
  expect_error(guard(xz, zx), "increasing rank", class = "mm_inference_unavailable")

  changed <- d
  changed$y <- 1 - changed$y
  other_y <- fit(y ~ x + (1 | g), changed)
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
