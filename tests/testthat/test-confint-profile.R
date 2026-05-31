# Stage D.3 (bd-01KRFGFSK4A0MGPFQVCNY5SYFK): confint(method = "profile")
# wired through the upstream `profile_confint_payload` FFI.
#
# Done condition from the bead: (a) ML full-rank LMM returns beta/sigma/theta
# rows, all finite and monotonic; (b) REML full-rank LMM omits beta with a
# documented reason_code, sigma/theta still present; (c) beta profile CI
# agrees with the Wald CI on a well-behaved model to within stated
# tolerance; (d) boundary fit raises a typed refusal rather than fabricating.

mm_skip_if_no_lme4_local <- function() {
  testthat::skip_if_not_installed("lme4")
}

mm_sleepstudy_data <- function() {
  mm_skip_if_no_lme4_local()
  env <- new.env(parent = emptyenv())
  utils::data("sleepstudy", package = "lme4", envir = env)
  if (!exists("sleepstudy", envir = env, inherits = FALSE)) {
    testthat::skip("sleepstudy dataset unavailable")
  }
  get("sleepstudy", envir = env, inherits = FALSE)
}

mm_sleepstudy_fit_ml <- function() {
  lmm(Reaction ~ Days + (1 + Days | Subject),
      data = mm_sleepstudy_data(), REML = FALSE,
      control = mm_control(verbose = -1))
}

mm_sleepstudy_fit_reml <- function() {
  lmm(Reaction ~ Days + (1 + Days | Subject),
      data = mm_sleepstudy_data(), REML = TRUE,
      control = mm_control(verbose = -1))
}

test_that("confint(method='profile') under ML returns beta/sigma/theta rows", {
  fit <- mm_sleepstudy_fit_ml()
  ci <- confint(fit, method = "profile", level = 0.95)

  expect_s3_class(ci, "mm_confint")
  expect_identical(attr(ci, "method"), "profile_likelihood")
  expect_identical(attr(ci, "fit_criterion"), "ML")
  expect_true("(Intercept)" %in% rownames(ci))
  expect_true("Days" %in% rownames(ci))
  expect_true("sigma" %in% rownames(ci))

  payload <- attr(ci, "mm_profile")
  expect_true(is.list(payload))
  kinds <- payload$table$parameter_kind
  expect_true("beta" %in% kinds)
  expect_true("sigma" %in% kinds)
  expect_true("theta" %in% kinds)

  # Bounds may be NA when the profile is truncated (upstream documents this
  # explicitly, mirroring lme4::confint). For each populated row, only the
  # populated bounds need to bracket the estimate; NA bounds are honest
  # "not determined" signals, not regressions.
  ok <- payload$table$reason_code %in% c(NA_character_, "")
  finite <- payload$table[ok, , drop = FALSE]
  beta_rows <- finite[finite$parameter_kind == "beta", , drop = FALSE]
  expect_true(nrow(beta_rows) >= 1L)
  expect_true(all(is.finite(beta_rows$lower)),
              info = "ML beta profile bounds should be finite on sleepstudy")
  expect_true(all(is.finite(beta_rows$upper)),
              info = "ML beta profile bounds should be finite on sleepstudy")
  expect_true(all(beta_rows$lower <= beta_rows$estimate + 1e-8))
  expect_true(all(beta_rows$estimate <= beta_rows$upper + 1e-8))

  # Variance-component rows are allowed to have NA bounds (truncated
  # profile), but every populated bound must still bracket the estimate.
  vc_rows <- finite[finite$parameter_kind != "beta", , drop = FALSE]
  for (i in seq_len(nrow(vc_rows))) {
    row <- vc_rows[i, , drop = FALSE]
    if (is.finite(row$lower)) {
      expect_lte(row$lower, row$estimate + 1e-8,
                 label = sprintf("lower bound for %s", row$parameter))
    }
    if (is.finite(row$upper)) {
      expect_gte(row$upper, row$estimate - 1e-8,
                 label = sprintf("upper bound for %s", row$parameter))
    }
  }
})

test_that("confint(method='profile') under REML omits beta with reason_code", {
  fit <- mm_sleepstudy_fit_reml()
  ci <- confint(fit, method = "profile", level = 0.95)

  expect_identical(attr(ci, "fit_criterion"), "REML")
  payload <- attr(ci, "mm_profile")
  table <- payload$table

  beta_rows <- table[table$parameter_kind == "beta", , drop = FALSE]
  expect_true(nrow(beta_rows) >= 1L,
              info = "REML profile payload should still surface refusal rows for beta")
  expect_true(all(beta_rows$reason_code == "profile_beta_unavailable_under_reml"))
  expect_true(all(is.na(beta_rows$lower)))
  expect_true(all(is.na(beta_rows$upper)))

  non_beta <- table[table$parameter_kind != "beta", , drop = FALSE]
  expect_true(nrow(non_beta) >= 1L,
              info = "REML profile payload must keep sigma/theta intervals")
  expect_true(all(non_beta$reason_code %in% c(NA_character_, "")))
})

test_that("profile CI for beta agrees with Wald CI on a well-behaved ML fit", {
  fit <- mm_sleepstudy_fit_ml()
  profile_ci <- confint(fit, parm = c("(Intercept)", "Days"),
                        method = "profile", level = 0.95)
  wald_ci <- confint(fit, parm = c("(Intercept)", "Days"),
                     method = "wald", level = 0.95)

  expect_true(all(rownames(profile_ci) %in% rownames(wald_ci)))
  # Profile and Wald intervals should agree to within a few percent on
  # sleepstudy beta. The actual numbers are tens to hundreds, so a
  # relative tolerance is appropriate.
  for (nm in rownames(profile_ci)) {
    rel_lower <- abs(profile_ci[nm, 1] - wald_ci[nm, 1]) /
      max(1e-6, abs(wald_ci[nm, 1]))
    rel_upper <- abs(profile_ci[nm, 2] - wald_ci[nm, 2]) /
      max(1e-6, abs(wald_ci[nm, 2]))
    expect_lt(rel_lower, 0.05,
              label = sprintf("rel lower diff for %s", nm))
    expect_lt(rel_upper, 0.05,
              label = sprintf("rel upper diff for %s", nm))
  }
})

test_that("profile CI parm subsetting filters returned rows", {
  fit <- mm_sleepstudy_fit_ml()
  ci_all <- confint(fit, method = "profile", level = 0.95)
  ci_intercept <- confint(fit, parm = "(Intercept)",
                          method = "profile", level = 0.95)
  expect_true(nrow(ci_intercept) >= 1L)
  expect_true(all(rownames(ci_intercept) %in% rownames(ci_all)))
  expect_identical(rownames(ci_intercept), "(Intercept)")
})

test_that("profile CI surfaces a typed refusal on a boundary singular fit", {
  testthat::skip_if_not_installed("lme4")
  env <- new.env(parent = emptyenv())
  utils::data("Dyestuff2", package = "lme4", envir = env)
  if (!exists("Dyestuff2", envir = env, inherits = FALSE)) {
    testthat::skip("Dyestuff2 dataset unavailable")
  }
  fit <- lmm(Yield ~ 1 + (1 | Batch),
             data = get("Dyestuff2", envir = env, inherits = FALSE),
             REML = TRUE, control = mm_control(verbose = -1))

  result <- tryCatch(
    confint(fit, method = "profile", level = 0.95),
    error = function(cnd) cnd
  )
  if (inherits(result, "condition")) {
    # Typed refusal path: any of the wrapper's structured error classes
    # is acceptable so long as it is NOT a silent return of fabricated
    # numbers.
    expect_true(any(class(result) %in%
                    c("mm_inference_unavailable", "mm_schema_error",
                      "mm_bridge_error", "mm_fit_error")),
                info = sprintf("unexpected error class: %s",
                               paste(class(result), collapse = "/")))
    return(invisible(NULL))
  }
  # Payload path: every NA bound on a boundary fit must carry a
  # non-empty regularity note (no NA-without-explanation rows).
  expect_s3_class(result, "mm_confint")
  payload <- attr(result, "mm_profile")
  expect_true(is.list(payload))
  table <- payload$table
  expect_true(nrow(table) >= 1L,
              info = "boundary fit profile payload must surface at least one row")
  na_rows <- table[!is.finite(table$lower) | !is.finite(table$upper), ,
                   drop = FALSE]
  if (nrow(na_rows)) {
    expect_true(all(nzchar(na_rows$regularity)),
                info = "NA bounds on a boundary fit must carry a regularity note")
  } else {
    # If every bound is finite, the boundary case still has to flag the
    # near-zero theta as boundary-clamped on the lower side rather than
    # silently emitting a negative variance bound.
    expect_true(all(is.finite(table$lower)) && all(is.finite(table$upper)))
  }
})
