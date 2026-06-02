# Tests for the caller optimizer-control escape hatch
# (mm_control optimizer / start / tolerance overrides), backed by mixeff-rs
# OptimizerControl. Default behaviour is the driver's automatic selection.

mm_oc_data <- function() {
  set.seed(1)
  ng <- 20L; per <- 4L
  g <- factor(rep(seq_len(ng), each = per)); n <- ng * per
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
  ng <- 15L; per <- 12L
  g <- factor(rep(seq_len(ng), each = per)); n <- ng * per
  df <- data.frame(y = rbinom(n, 1, 0.5), x = rnorm(n), g = g)
  f <- glmm(y ~ x + (1 | g), df, family = binomial(), method = "joint_laplace",
            control = mm_control(verbose = -1, max_feval = 50000L))
  expect_s3_class(f, "mm_glmm")
})
