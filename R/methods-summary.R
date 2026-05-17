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
