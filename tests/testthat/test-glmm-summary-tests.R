## summary(mm_glmm, tests = "coefficients") -- Wald-z coefficient table
## once the upstream fixed-effect covariance payload is available.

mk_summary_glmm <- function(seed = 41L) {
  data <- mk_summary_glmm_data()
  glmm(
    y ~ x + (1 | subject),
    data,
    family = Gamma(link = "log"),
    method = "joint_laplace",
    nAGQ = 1L,
    control = mm_control(verbose = -1)
  )
}

mk_summary_glmm_data <- function() {
  y <- x <- numeric()
  subject <- character()
  for (group in 0:3) {
    for (obs in 0:4) {
      xv <- obs - 2
      eta <- 0.5 + 0.2 * xv + (group - 1.5) * 0.08
      y <- c(y, exp(eta) * (0.95 + 0.02 * ((group + obs) %% 3)))
      x <- c(x, xv)
      subject <- c(subject, paste0("g", group + 1L))
    }
  }
  data.frame(y = y, x = x, subject = factor(subject))
}

test_that("summary(mm_glmm, tests='coefficients') no longer raises mm_inference_unavailable", {
  fit <- mk_summary_glmm()
  expect_no_error(summary(fit, tests = "coefficients"))
  expect_s3_class(summary(fit, tests = "coefficients"), "summary.mm_glmm")
})

test_that("summary(mm_glmm, tests='none') still works", {
  fit <- mk_summary_glmm()
  s_none <- summary(fit, tests = "none")
  expect_s3_class(s_none, "summary.mm_glmm")
  expect_null(s_none$inference)
  ## The 'none' path should give a coefficient table with NA stats
  expect_true("Estimate" %in% colnames(s_none$coefficients))
  expect_true(all(is.na(s_none$coefficients$df) |
                  is.numeric(s_none$coefficients$df)))
})

test_that("Wald z = beta / SE and p = 2*pnorm(|z|) match mm_lincomb()", {
  fit <- mk_summary_glmm()
  s <- summary(fit, tests = "coefficients")
  coef <- s$coefficients

  ## Compute the Wald z for each row via mm_lincomb and compare
  bnms <- names(fixef(fit))
  for (nm in bnms) {
    out <- mm_lincomb(fit, setNames(1, nm))
    row <- coef[nm, , drop = FALSE]
    expect_equal(row[["Estimate"]],   out$estimate,  tolerance = 1e-12,
                 info = sprintf("Estimate mismatch on '%s'", nm))
    expect_equal(row[["Std. Error"]], out$std_error, tolerance = 1e-12,
                 info = sprintf("SE mismatch on '%s'", nm))
    expect_equal(row[["z value"]],    out$statistic, tolerance = 1e-10,
                 info = sprintf("z mismatch on '%s'", nm))
    expect_equal(row[["Pr(>|z|)"]],   out$p_value,   tolerance = 1e-10,
                 info = sprintf("p mismatch on '%s'", nm))
  }
})

test_that("summary surfaces mm_vcov_status reflecting the upstream payload reliability", {
  fit <- mk_summary_glmm()
  s <- summary(fit, tests = "coefficients")
  expect_true(is.list(s$vcov_status))
  expect_true(all(c("status", "method", "reliability", "reason") %in%
                  names(s$vcov_status)))
  expect_identical(s$vcov_status$status, "available")
  expect_identical(s$vcov_status$method, "asymptotic_wald_z")
})

test_that("print(summary(mm_glmm, tests='coefficients')) writes the reliability notice", {
  fit <- mk_summary_glmm()
  s <- summary(fit, tests = "coefficients")
  out <- capture.output(print(s))
  expect_true(any(grepl("Fixed effects:", out, fixed = TRUE)))
  expect_true(any(grepl("Wald-z reliability:", out, fixed = TRUE)))
})

test_that("Coefficient table uses lme4-comparable column labels", {
  fit <- mk_summary_glmm()
  s <- summary(fit, tests = "coefficients")
  cols <- colnames(s$coefficients)
  expect_true(all(c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in% cols))
})

test_that("summary keeps fast-PIRLS GLMM Wald rows explicitly unavailable", {
  fit <- glmm(
    y ~ x + (1 | subject),
    mk_summary_glmm_data(),
    family = Gamma(link = "log"),
    method = "pirls_profiled",
    nAGQ = 1L,
    control = mm_control(verbose = -1)
  )
  coef <- summary(fit, tests = "coefficients")$coefficients
  expect_true(all(coef$method == "not_computed"))
  expect_true(all(is.na(coef[["Pr(>|z|)"]])))
})
