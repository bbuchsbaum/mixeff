#' @method summary mm_lmm
#' @export
summary.mm_lmm <- function(object, tests = c("coefficients", "none"),
                           method = c("auto", "satterthwaite",
                                      "kenward_roger", "bootstrap",
                                      "asymptotic", "none"), ...) {
  tests <- match.arg(tests)
  method <- match.arg(method)
  inference <- if (identical(tests, "coefficients")) {
    inf_method <- if (method %in% c("auto", "satterthwaite", "kenward_roger",
                                    "asymptotic", "none")) {
      method
    } else {
      "auto"
    }
    inference_table(object, method = inf_method)
  } else {
    NULL
  }
  coef <- mm_summary_coefficients(object, inference)

  out <- list(
    call = object$call,
    formula = object$formula,
    REML = object$REML,
    coefficients = coef,
    sigma = object$sigma,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    nobs = object$nobs,
    df_residual = object$df_residual,
    fit_status = object$fit_status,
    varcorr = object$varcorr,
    tests = tests,
    inference = inference,
    requested_method = method
  )
  class(out) <- "summary.mm_lmm"
  out
}

#' @method print summary.mm_lmm
#' @export
print.summary.mm_lmm <- function(x, ...) {
  cat(sprintf("Linear mixed model fit by %s\n", if (x$REML) "REML" else "ML"))
  cat(sprintf("Formula: %s\n", deparse1(x$formula)))
  cat(sprintf("Fit status: %s\n\n", x$fit_status))
  print(x$varcorr)
  cat("\nFixed effects:\n")
  print(x$coefficients)
  if (!is.null(x$inference)) {
    inf <- x$inference$table
    cat("\nInference status:\n")
    cols <- intersect(
      c("term", "method", "status", "reliability", "reliability_reason"),
      names(inf)
    )
    print(inf[, cols, drop = FALSE], row.names = FALSE)
    has_reason <- "reason" %in% names(inf) &&
      any(!is.na(inf$reason) & nzchar(inf$reason))
    if (has_reason) {
      cat("\nReasons:\n")
      with_reason <- inf[!is.na(inf$reason) & nzchar(inf$reason), , drop = FALSE]
      for (i in seq_len(nrow(with_reason))) {
        cat(sprintf("  %s: %s\n",
                    with_reason$term[[i]],
                    with_reason$reason[[i]]))
      }
    }
  }
  if (mm_summary_verbose(...)) {
    cat("\nInference rows are supplied by Rust; `reliability_reason` is a closed-enum warrant for the reliability grade.\n")
  }
  invisible(x)
}

#' @method summary mm_glmm
#' @export
summary.mm_glmm <- function(object, tests = c("none", "coefficients"), ...) {
  tests <- match.arg(tests)
  inference <- NULL
  vcov_status <- NULL
  if (identical(tests, "coefficients")) {
    inference <- mm_glmm_wald_z_inference(object)
    vcov_status <- attr(inference, "mm_vcov_status")
  }
  coef <- mm_summary_coefficients(object, inference)
  out <- list(
    call = object$call,
    formula = object$formula,
    family = object$family,
    method = object$method,
    nAGQ = object$nAGQ,
    coefficients = coef,
    dispersion = object$dispersion,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    nobs = object$nobs,
    df_residual = object$df_residual,
    fit_status = object$fit_status,
    varcorr = object$varcorr,
    tests = tests,
    inference = inference,
    vcov_status = vcov_status
  )
  class(out) <- "summary.mm_glmm"
  out
}

## Build a Wald-z coefficient inference table for an mm_glmm fit from
## the stored fixed-effect covariance payload. Falls back to NA stats
## when the payload reports itself unavailable (e.g. revived fits that
## lost the matrix); the diagonal-from-std_errors flavor counts as
## available for univariate Wald z, but the mm_vcov_status attribute
## still flags the moderate / unavailable reliability for callers.
mm_glmm_wald_z_inference <- function(object) {
  beta <- as.numeric(object$beta)
  names(beta) <- names(object$beta)
  Vfull <- stats::vcov(object)
  V <- as.matrix(unclass(Vfull))
  status <- list(
    status      = attr(Vfull, "mm_status")      %||% "available",
    method      = attr(Vfull, "mm_method")      %||% NA_character_,
    reliability = attr(Vfull, "mm_reliability") %||% NA_character_,
    reason      = attr(Vfull, "mm_reason")      %||%
                  attr(Vfull, "mm_unavailable_reason") %||% NA_character_
  )
  se <- if (!is.null(object$std_errors) &&
            length(object$std_errors) == length(beta)) {
    as.numeric(object$std_errors)
  } else {
    suppressWarnings(sqrt(diag(V)))
  }
  z <- beta / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  method_used <- if (identical(status$status, "available")) {
    "asymptotic"
  } else {
    "not_computed"
  }
  if (identical(method_used, "not_computed")) {
    z[] <- NA_real_; p[] <- NA_real_
  }
  rows <- data.frame(
    label             = names(beta),
    kind              = "coefficient",
    estimate          = unname(beta),
    std_error         = unname(se),
    denominator_df    = NA_real_,
    statistic         = unname(z),
    statistic_name    = ifelse(is.na(z), NA_character_, "z"),
    p_value           = unname(p),
    method            = method_used,
    status            = status$status,
    reliability       = status$reliability,
    reliability_reason = status$reason,
    reason            = status$reason,
    stringsAsFactors  = FALSE
  )
  out <- list(table = rows, raw = NULL)
  class(out) <- c("mm_inference_table", "list")
  attr(out, "mm_vcov_status") <- status
  out
}

#' @method print summary.mm_glmm
#' @export
print.summary.mm_glmm <- function(x, ...) {
  cat("Generalized linear mixed model fit\n")
  cat(sprintf("Formula: %s\n", deparse1(x$formula)))
  cat(sprintf("Family/link: %s/%s\n", x$family$family, x$family$link))
  cat(sprintf("Method: %s (nAGQ = %d)\n", x$method, x$nAGQ))
  cat(sprintf("Fit status: %s\n\n", x$fit_status))
  print(x$varcorr)
  cat("\nFixed effects:\n")
  print(x$coefficients)
  if (!is.null(x$vcov_status) && !is.null(x$inference)) {
    rel <- x$vcov_status$reliability
    if (!is.na(rel) && nzchar(rel) && !identical(rel, "available")) {
      cat(sprintf(
        "\nWald-z reliability: %s (%s).",
        rel,
        x$vcov_status$method %||% "no method tag"
      ))
      if (!is.na(x$vcov_status$reason) && nzchar(x$vcov_status$reason)) {
        cat(sprintf(" Reason: %s.", x$vcov_status$reason))
      }
      cat("\n")
    }
  }
  invisible(x)
}

# Gate the long inference-rows footer behind explicit verbosity. Default
# summary() / print(summary(.)) output is now compact; users opt in via
# print(summary(fit), verbose = TRUE) or by setting
# options(mixeff.verbose = 1L). Keeps the audit-first contract intact
# (the reliability_reason field is still on every row in `x$inference`)
# while shaving teaching-debt off the REPL surface.
mm_summary_verbose <- function(...) {
  args <- list(...)
  if (isTRUE(args$verbose)) return(TRUE)
  isTRUE(getOption("mixeff.verbose", 0L) >= 1L)
}

mm_summary_coefficients <- function(object, inference) {
  rows <- if (!is.null(inference)) inference$table else NULL
  if (!is.null(rows) && nrow(rows) && "kind" %in% names(rows)) {
    rows <- rows[rows$kind == "coefficient", , drop = FALSE]
  }
  if (is.null(rows) || !nrow(rows)) {
    rows <- data.frame(
      label = names(object$beta),
      estimate = unname(object$beta),
      std_error = unname(object$std_errors),
      denominator_df = NA_real_,
      statistic = NA_real_,
      statistic_name = NA_character_,
      p_value = NA_real_,
      method = "not_computed",
      stringsAsFactors = FALSE
    )
  } else {
    rows <- rows[match(names(object$beta), rows$label), , drop = FALSE]
    missing <- is.na(rows$label)
    if (any(missing)) {
      rows$label[missing] <- names(object$beta)[missing]
      rows$estimate[missing] <- unname(object$beta)[missing]
      rows$std_error[missing] <- unname(object$std_errors)[missing]
      rows$denominator_df[missing] <- NA_real_
      rows$statistic[missing] <- NA_real_
      rows$statistic_name[missing] <- NA_character_
      rows$p_value[missing] <- NA_real_
      rows$method[missing] <- "not_computed"
    }
  }

  stat_col <- mm_summary_statistic_column(rows$statistic_name)
  p_col <- mm_summary_p_value_column(rows$statistic_name)
  coef <- data.frame(
    Estimate = rows$estimate,
    `Std. Error` = rows$std_error,
    df = rows$denominator_df,
    method = rows$method,
    check.names = FALSE
  )
  coef[[stat_col]] <- rows$statistic
  coef[[p_col]] <- rows$p_value
  coef <- coef[, c("Estimate", "Std. Error", "df", stat_col, p_col, "method"),
               drop = FALSE]
  rownames(coef) <- rows$label
  coef
}

mm_summary_statistic_column <- function(statistic_name) {
  values <- unique(na.omit(as.character(statistic_name)))
  if (length(values) == 1L) {
    switch(values,
           z = "z value",
           t = "t value",
           f = "F value",
           chi_square = "Chisq",
           "statistic")
  } else {
    "statistic"
  }
}

mm_summary_p_value_column <- function(statistic_name) {
  values <- unique(na.omit(as.character(statistic_name)))
  if (length(values) == 1L) {
    switch(values,
           z = "Pr(>|z|)",
           t = "Pr(>|t|)",
           f = "Pr(>F)",
           chi_square = "Pr(>Chisq)",
           "p.value")
  } else {
    "p.value"
  }
}
