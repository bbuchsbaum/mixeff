# No silent surgery (PRD §8.1): recognized lme4 arguments that the wrapper
# cannot honor must raise a typed condition, never be swallowed by `...` and
# return a plausible-but-wrong value. Cheap arguments that CAN be honored
# (REML matching the fit, pearson/scaled residuals, vcov correlation) must be
# honored correctly.

mk_ss_fit <- function() {
  suppressMessages(
    lmm(Reaction ~ Days + (Days | Subject), lme4::sleepstudy,
        control = mm_control(verbose = -1))
  )
}

test_that("logLik/deviance refuse a mismatched REML request", {
  m <- mk_ss_fit()  # REML = TRUE
  expect_error(logLik(m, REML = FALSE), class = "mm_inference_unavailable")
  expect_error(deviance(m, REML = FALSE), class = "mm_inference_unavailable")
})

test_that("logLik honors a REML request that matches the fit", {
  m <- mk_ss_fit()
  expect_equal(as.numeric(logLik(m, REML = TRUE)), as.numeric(logLik(m)))
})

test_that("predict refuses random.only and points to re.form", {
  m <- mk_ss_fit()
  expect_error(predict(m, random.only = TRUE), class = "mm_arg_error")
})

test_that("simulate refuses newparams / newdata / use.u", {
  m <- mk_ss_fit()
  expect_error(simulate(m, newparams = list(beta = 1)), class = "mm_arg_error")
  expect_error(simulate(m, newdata = lme4::sleepstudy[1:5, ]),
               class = "mm_arg_error")
  expect_error(simulate(m, use.u = TRUE), class = "mm_arg_error")
})

test_that("refit refuses newweights but still accepts control via ...", {
  m <- mk_ss_fit()
  expect_error(
    refit(m, lme4::sleepstudy$Reaction, newweights = rep(1, 180)),
    class = "mm_arg_error"
  )
  # control in ... must still work (not rejected)
  expect_s3_class(
    refit(m, lme4::sleepstudy$Reaction, control = mm_control(verbose = -1)),
    "mm_lmm"
  )
})

test_that("residuals honor type and scaled", {
  m <- mk_ss_fit()
  raw <- residuals(m)
  expect_equal(unname(residuals(m, type = "pearson")),
               unname(raw / sigma(m)))
  expect_equal(unname(residuals(m, scaled = TRUE)),
               unname(raw / sigma(m)))
  # working/deviance equal response for a Gaussian LMM
  expect_equal(unname(residuals(m, type = "working")), unname(raw))
})

test_that("vcov honors correlation = TRUE per the lme4 contract", {
  m <- mk_ss_fit()
  V <- vcov(m, correlation = TRUE)
  expect_equal(attr(V, "correlation"), stats::cov2cor(vcov(m)))
})
