# Convergence-honesty invariant (Audit1 remediation WI-1.4).
#
# A `converged_*` fit status is a certificate: it must never coexist with a
# non-finite objective, log-likelihood, coefficient vector, theta vector, or
# fitted values. The engine enforces this for the primary optimum since
# 9b3b4ad (non_finite_objective -> not_optimized); this helper is the
# wrapper-side guard that a regression -- or a new path that skips the engine
# check (e.g. a fallback or revival route) -- cannot reach users as a false
# certificate.
#
# `mm_check_convergence_honesty()` works on plain field access (no S3), so
# tests can also feed it deliberately corrupted fake fits to prove it fires.

mm_check_convergence_honesty <- function(fit) {
  status <- as.character(fit$fit_status %||% "")
  # %||% only catches NULL; a corrupted fake can carry NA_character_, which
  # must fail cleanly (out of scope), not crash startsWith()/if().
  if (length(status) != 1L || is.na(status)) {
    status <- ""
  }
  if (!startsWith(status, "converged")) {
    return(character(0))
  }
  problems <- character(0)
  finite_all <- function(x) length(x) > 0L && all(is.finite(x))

  if (!finite_all(as.numeric(fit$logLik))) {
    problems <- c(problems, "logLik")
  }
  if (!finite_all(as.numeric(fit$beta))) {
    problems <- c(problems, "beta")
  }
  # theta may legitimately be length-0 only for fixed-effect-only models,
  # which this package refuses; require finite values when present.
  theta <- as.numeric(fit$theta)
  if (length(theta) && !all(is.finite(theta))) {
    problems <- c(problems, "theta")
  }
  fitted_vals <- as.numeric(fit$fitted)
  if (length(fitted_vals) && !all(is.finite(fitted_vals))) {
    problems <- c(problems, "fitted")
  }
  cert_obj <- fit$artifact$optimizer_certificate$objective_value
  if (!is.null(cert_obj) && !all(is.finite(as.numeric(cert_obj)))) {
    problems <- c(problems, "certificate_objective")
  }
  problems
}

mm_expect_convergence_honesty <- function(fit, label = "fit") {
  problems <- mm_check_convergence_honesty(fit)
  testthat::expect_identical(
    problems,
    character(0),
    info = sprintf(
      "%s: fit_status `%s` coexists with non-finite quantities: %s",
      label,
      as.character(fit$fit_status %||% ""),
      paste(problems, collapse = ", ")
    )
  )
}
