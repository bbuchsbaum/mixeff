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
    if (identical(inf_method, "auto") &&
        !mm_boundary_df_method_unavailable(object, "satterthwaite")) {
      # "auto" resolves to the finite-sample Satterthwaite route whenever it
      # is feasible; the engine's fit-time cached table is the asymptotic
      # Wald-z fallback and stays the honest answer only where Satterthwaite
      # df are refused (singular / boundary optimum).
      inf_method <- "satterthwaite"
    }
    tbl <- inference_table(object, method = inf_method)
    if (identical(inf_method, "satterthwaite") && identical(method, "auto") &&
        !any(tbl$table$status == "available")) {
      # Satterthwaite was refused for every coefficient (the boundary can be
      # reached without is_singular() flagging it, e.g. a variance pinned at
      # zero on a not_optimized fit). Only an unrequested route may be
      # swapped: the cached rows are labeled asymptotic_wald_z and carry
      # their own engine fallback note, so nothing is hidden.
      tbl <- inference_table(object, method = "auto")
    }
    tbl
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
    varcorr = VarCorr(object),
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
  print(mm_summary_format_coef(x$coefficients))
  notes <- mm_fit_status_note(x$fit_status)
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
    # Engine-authored row notes are the prose warrant behind the method and
    # reliability columns (e.g. why a fallback or finite-difference route was
    # used); surface them instead of leaving them buried in the table.
    row_notes <- unique(unlist(inf$notes %||% list(), use.names = FALSE))
    row_notes <- row_notes[!is.na(row_notes) & nzchar(row_notes)]
    notes <- c(notes, row_notes)
  }
  mm_summary_print_notes(notes)
  if (mm_summary_verbose(...)) {
    cat("\nInference rows are supplied by Rust; `reliability_reason` is a closed-enum warrant for the reliability grade.\n")
  }
  invisible(x)
}

#' @method summary mm_glmm
#' @export
summary.mm_glmm <- function(object, tests = c("coefficients", "none"), ...) {
  # Coefficient tests are the default (matching lme4::glmer's summary): the
  # inference rows carry their own status/reliability labels, so when the fit
  # method cannot certify SE/z/p the columns are withheld WITH the reason
  # printed, rather than silently absent under a tests = "none" default.
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
    varcorr = VarCorr(object),
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
  parsed <- mm_json_parse_fixed_effect_inference_table(
    object$artifact$fixed_effect_inference_table %||% NULL
  )
  if (!is.null(parsed)) {
    rows <- parsed$table
    if ("kind" %in% names(rows)) {
      rows <- rows[rows$kind == "coefficient", , drop = FALSE]
    }
    if (nrow(rows)) {
      # Artifact rows carry engine-encoded labels; fit$beta is lme4-named.
      rows$label <- mm_coef_engine_to_lme4(rows$label, object$coef_map)
      rows <- rows[match(names(object$beta), rows$label), , drop = FALSE]
      out <- list(table = rows, raw = parsed$raw)
      class(out) <- c("mm_inference_table", "list")
      attr(out, "mm_vcov_status") <- mm_glmm_inference_rows_status(rows)
      return(out)
    }
  }

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
    z[] <- NA_real_
    p[] <- NA_real_
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

mm_glmm_inference_rows_status <- function(rows) {
  compact <- function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) paste(x, collapse = ", ") else NA_character_
  }
  reason <- compact(rows$reason)
  list(
    status = if (length(rows$status) &&
                 all(!is.na(rows$status) & rows$status == "available")) {
      "available"
    } else {
      compact(rows$status)
    },
    method = compact(rows$method),
    reliability = compact(rows$reliability),
    reason = reason
  )
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
  print(mm_summary_format_coef(x$coefficients))
  reason_printed <- FALSE
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
        reason_printed <- TRUE
      }
      cat("\n")
    }
  }
  notes <- c(
    mm_fit_status_note(x$fit_status),
    mm_glmm_withheld_inference_note(x, include_reason = !reason_printed)
  )
  mm_summary_print_notes(notes)
  invisible(x)
}

# One plain-language sentence when every test statistic in the table is
# withheld, so a summary() full of NA columns is never left unexplained.
# Points at the certified estimator when the fit used the uncertified
# default -- an available option reported as fact, not a model prescription.
mm_glmm_withheld_inference_note <- function(x, include_reason = TRUE) {
  coef <- x$coefficients
  stat_cols <- intersect(c("z value", "t value", "statistic"), names(coef))
  if (!length(stat_cols)) return(character())
  stats <- coef[[stat_cols[[1L]]]]
  if (!length(stats) || !all(is.na(stats))) return(character())
  note <- paste0(
    "test statistics and p-values are withheld: the fit's covariance payload ",
    "does not certify fixed-effect inference"
  )
  # The engine reason is skipped when the Wald-z reliability line above the
  # notes already printed it verbatim; repeating a paragraph-length warrant
  # twice reads as noise, not honesty.
  reason <- x$vcov_status$reason %||% NA_character_
  if (include_reason && !is.na(reason) && nzchar(reason)) {
    note <- sprintf("%s (engine reason: %s)", note, reason)
  }
  if (identical(x$method, "pirls_profiled")) {
    note <- paste0(
      note,
      ". Engine-certified Wald inference is available from a fit with ",
      'method = "joint_laplace".'
    )
  } else {
    note <- paste0(note, ".")
  }
  note
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

# A non-converged optimum is model state the user must not read past: repeat
# it as a plain-language note next to the tests, not only in the header line.
# Reports what happened; prescribes nothing (PRD R9).
mm_fit_status_note <- function(fit_status) {
  status <- as.character(fit_status %||% "")
  if (!nzchar(status) || startsWith(status, "converged")) return(character())
  sprintf(
    paste0(
      "fit status `%s`: the optimizer stopped without certifying an optimum; ",
      "estimates, variance components, and any tests above are reported from ",
      "the last accepted iterate."
    ),
    status
  )
}

mm_summary_print_notes <- function(notes) {
  notes <- unique(notes[!is.na(notes) & nzchar(notes)])
  if (!length(notes)) return(invisible(NULL))
  cat("\nNotes:\n")
  for (n in notes) {
    cat(sprintf("  %s\n", n))
  }
  invisible(NULL)
}

# Display copy of the coefficient table: p-values render through
# format.pval() so an underflowed value prints "< 1e-16" instead of the
# fabricated-certainty "0.000000e+00". The stored table stays numeric.
mm_summary_format_coef <- function(coef) {
  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "Pr(>F)", "Pr(>Chisq)", "p.value"),
                     names(coef))
  if (!length(p_col)) return(coef)
  out <- coef
  for (col in p_col) {
    out[[col]] <- format.pval(out[[col]], digits = 4, eps = 1e-16)
  }
  out
}

mm_summary_coefficients <- function(object, inference) {
  beta_names <- names(object$beta)
  beta <- unname(object$beta)
  se <- if (!is.null(object$std_errors) &&
            length(object$std_errors) == length(object$beta)) {
    unname(object$std_errors)
  } else {
    rep(NA_real_, length(object$beta))
  }
  rows <- if (!is.null(inference)) inference$table else NULL
  if (!is.null(rows) && nrow(rows) && "kind" %in% names(rows)) {
    rows <- rows[rows$kind == "coefficient", , drop = FALSE]
  }
  if (is.null(rows) || !nrow(rows)) {
    rows <- data.frame(
      label = beta_names,
      estimate = beta,
      std_error = se,
      denominator_df = NA_real_,
      statistic = NA_real_,
      statistic_name = NA_character_,
      p_value = NA_real_,
      method = "not_computed",
      stringsAsFactors = FALSE
    )
  } else {
    rows <- rows[match(beta_names, rows$label), , drop = FALSE]
    missing <- is.na(rows$label)
    if (any(missing)) {
      rows$label[missing] <- beta_names[missing]
      rows$estimate[missing] <- beta[missing]
      rows$std_error[missing] <- se[missing]
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
  cols <- c("Estimate", "Std. Error", "df", stat_col, p_col, "method")
  if (all(is.na(coef$df))) {
    # df is undefined for every row (e.g. asymptotic Wald z); an all-NA
    # column reads as breakage, so omit it rather than print NA.
    cols <- setdiff(cols, "df")
  }
  coef <- coef[, cols, drop = FALSE]
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
