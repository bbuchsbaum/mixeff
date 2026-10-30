# Tests for the caller optimizer-control escape hatch
# (mm_control optimizer / start / tolerance overrides), backed by mixeff-rs
# OptimizerControl. Default behaviour is the driver's automatic selection.

mm_oc_data <- function() {
  set.seed(1)
  ng <- 20L
  per <- 4L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  data.frame(y = rnorm(n), x = rnorm(n), g = g)
}

test_that("mm_control validates the optimizer-control fields", {
  expect_error(mm_control(optimizer = 42), class = "mm_arg_error")
  expect_error(mm_control(optimizer = c("a", "b")), class = "mm_arg_error")
  expect_error(mm_control(start = "a"), class = "mm_arg_error")
  expect_error(mm_control(ftol_rel = -1), class = "mm_arg_error")
  expect_error(mm_control(ftol_rel = c(1, 2)), class = "mm_arg_error")
  ctl <- mm_control(optimizer = "cobyla", start = 0.5, ftol_rel = 1e-8)
  expect_s3_class(ctl, "mm_control")
  expect_identical(ctl$optimizer, "cobyla")
  expect_equal(ctl$start, 0.5)
})

test_that("an optimizer override is honored and recorded in the certificate", {
  df <- mm_oc_data()
  f0 <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  f1 <- lmm(y ~ x + (1 | g), df,
            control = mm_control(verbose = -1, optimizer = "cobyla"))
  # the chosen optimizer is recorded (provenance / auditability)
  expect_identical(optimizer_certificate(f1)$raw$optimizer_name, "cobyla")
  # same model, different optimizer -> same answer within tolerance
  expect_equal(unname(fixef(f1)), unname(fixef(f0)), tolerance = 1e-3)
})

test_that("warm start (start = theta) reproduces the fit", {
  df <- mm_oc_data()
  f0 <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
  f2 <- lmm(y ~ x + (1 | g), df,
            control = mm_control(verbose = -1, start = f0$theta))
  expect_equal(unname(fixef(f2)), unname(fixef(f0)), tolerance = 1e-6)
  expect_equal(f2$theta, f0$theta, tolerance = 1e-4)
})

test_that("an unknown optimizer raises a typed error", {
  df <- mm_oc_data()
  expect_error(
    lmm(y ~ x + (1 | g), df,
        control = mm_control(verbose = -1, optimizer = "nope")),
    class = "mm_arg_error"
  )
})

test_that("a wrong-length warm start is rejected", {
  df <- mm_oc_data()
  # y ~ x + (1 | g) has a single theta; a length-5 start must be refused.
  expect_error(
    lmm(y ~ x + (1 | g), df,
        control = mm_control(verbose = -1, start = c(0.1, 0.2, 0.3, 0.4, 0.5))),
    class = "mm_condition"
  )
})

test_that("tolerance overrides are accepted and produce a fit", {
  df <- mm_oc_data()
  f <- lmm(y ~ x + (1 | g), df,
           control = mm_control(verbose = -1, ftol_rel = 1e-10,
                                ftol_abs = 1e-10, xtol_rel = 1e-9))
  expect_s3_class(f, "mm_lmm")
})

test_that("GLMM still honors max_feval through the optimizer-control path", {
  set.seed(2)
  ng <- 15L
  per <- 12L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  df <- data.frame(y = rbinom(n, 1, 0.5), x = rnorm(n), g = g)
  f <- glmm(y ~ x + (1 | g), df, family = binomial(), method = "joint_laplace",
            control = mm_control(verbose = -1, max_feval = 50000L))
  expect_s3_class(f, "mm_glmm")
})

test_that("scalar theta controls equal explicit vectors in crossed fits and refits", {
  set.seed(20261027)
  d <- expand.grid(subject = factor(seq_len(24)), item = factor(seq_len(20)))
  d$x <- runif(nrow(d), -1, 1)
  u <- rnorm(24, 0, 0.6)
  v <- rnorm(20, 0, 0.4)
  eta <- 0.2 + 0.4 * d$x + u[d$subject] + v[d$item]
  d$y <- eta + rnorm(nrow(d), 0, 0.7)
  d$binary <- rbinom(nrow(d), 1, plogis(eta))
  make_control <- function(n) mm_control(verbose = -1, optCtrl = list(
    xtol_abs = rep(1e-8, n), rhobeg = rep(0.1, n)))
  f <- y ~ x + (1 | subject) + (1 | item)
  scalar <- lmm(f, d, control = make_control(1))
  vector <- lmm(f, d, control = make_control(2))
  expect_equal(fixef(scalar), fixef(vector), tolerance = 1e-12)
  expect_equal(scalar$theta, vector$theta, tolerance = 1e-12)
  expect_equal(as.numeric(logLik(scalar)), as.numeric(logLik(vector)), tolerance = 1e-12)
  # Retain the scalar request so new/refitted models expand against their own
  # theta dimension rather than recycling an old expanded vector.
  expect_length(scalar$control$xtol_abs, 1L)
  reduced <- update(scalar, . ~ x + (1 | subject))
  expect_length(reduced$theta, 1L)
  restored <- unserialize(serialize(scalar, NULL))
  ml <- update(restored, REML = FALSE)
  ml_vector <- lmm(f, d, REML = FALSE, control = make_control(2))
  expect_equal(fixef(ml), fixef(ml_vector), tolerance = 1e-12)
  for (method in c("pirls_profiled", "joint_laplace")) {
    fit <- function(n) glmm(binary ~ x + (1 | subject) + (1 | item), d,
                            family = binomial(), method = method,
                            control = make_control(n))
    a <- fit(1)
    b <- fit(2)
    expect_equal(fixef(a), fixef(b), tolerance = 1e-12)
    expect_equal(a$theta, b$theta, tolerance = 1e-12)
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-12)
  }
  for (name in c("xtol_abs", "rhobeg")) {
    bad <- setNames(list(rep(0.1, 3)), name)
    expect_error(lmm(f, d, control = mm_control(verbose = -1, optCtrl = bad)),
                 "length 3 does not match theta length 2", class = "mm_condition")
    for (value in c(0, -1, Inf, NA_real_)) {
      expect_error(mm_control(optCtrl = setNames(list(value), name)),
                   class = "mm_arg_error")
    }
  }
})

test_that("scalar absolute tolerance broadcasts across ten covariance parameters", {
  set.seed(20261028)
  d <- data.frame(g = factor(rep(seq_len(20), each = 10)),
                  x = rnorm(200), z = rnorm(200), w = rnorm(200))
  X <- model.matrix(~ x + z + w, d)
  b <- matrix(rnorm(80, sd = 0.3), 20, 4)
  d$y <- drop(X %*% c(0.2, 0.4, -0.1, 0.3)) +
    rowSums(X * b[d$g, ]) + rnorm(200, sd = 0.7)
  fit <- function(tol) lmm(y ~ x + z + w + (1 + x + z + w | g), d,
                           control = mm_control(verbose = -1,
                                                optCtrl = list(xtol_abs = tol)))
  scalar <- fit(1e-8)
  vector <- fit(rep(1e-8, 10))
  expect_length(scalar$theta, 10L)
  expect_equal(scalar$theta, vector$theta, tolerance = 1e-12)
  expect_equal(fixef(scalar), fixef(vector), tolerance = 1e-12)
})
