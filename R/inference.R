#' Contrast fixed effects
#'
#' Note: this is not `emmeans::contrast`. `contrast()` is mixeff's fixed-effect
#' contrast front door. R validates the contrast matrix shape, then asks Rust
#' to evaluate estimability, method prerequisites, standard errors, degrees of
#' freedom, statistics, p-values, reliability, and unavailable reasons.
#'
#' @param fit A fitted `mm_lmm`.
#' @param L A numeric contrast vector or matrix with one column per fixed
#'   effect.
#' @param rhs Numeric right-hand side, recycled to the number of contrasts.
#' @param method Requested inference method.
#' @param bootstrap Optional [bootstrap_control()] object for
#'   `method = "bootstrap"`.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_contrast` object with a data-frame `table`. The `estimate`
#'   column is the tested difference, `L beta_hat - rhs`.
#'
#' @export
contrast <- function(fit, L, rhs = 0, method = c("auto", "satterthwaite",
                                                 "kenward_roger", "bootstrap",
                                                 "asymptotic", "none"),
                     bootstrap = NULL, ...) {
  UseMethod("contrast")
}

#' @rdname contrast
#' @export
contrast.mm_lmm <- function(fit, L, rhs = 0, method = c("auto", "satterthwaite",
                                                        "kenward_roger", "bootstrap",
                                                        "asymptotic", "none"),
                            bootstrap = NULL, ...) {
  method <- match.arg(method)
  L <- mm_contrast_matrix(L, fit)
  rhs <- rep(as.numeric(rhs), length.out = nrow(L))
  if (identical(method, "none")) {
    estimate <- as.numeric(L %*% fit$beta) - rhs
    table <- data.frame(
      contrast = rownames(L),
      estimate = estimate,
      rhs = rhs,
      std_error = NA_real_,
      df = NA_real_,
      statistic = NA_real_,
      statistic_name = NA_character_,
      p_value = NA_real_,
      method = "not_computed",
      requested_method = method,
      status = "not_assessed",
      reliability = "not_available",
      estimability = "not_assessed",
      reason = mm_inference_unavailable_reason(method),
      notes = I(rep(list(character()), nrow(L))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    obj <- list(table = table, L = L, rhs = rhs, requested_method = method, raw = NULL)
    class(obj) <- "mm_contrast"
    return(obj)
  }

  parsed <- mm_rust_contrast_table(fit, L, rhs, method, bootstrap = bootstrap)
  table <- parsed$table
  table$contrast <- table$label
  table$rhs <- rhs
  table$requested_method <- method
  table <- table[, c("contrast", "estimate", "rhs", "std_error", "df",
                     "statistic", "statistic_name", "p_value", "method",
                     "requested_method", "status", "reliability",
                     "estimability", "reason", "details", "notes"),
                 drop = FALSE]
  obj <- list(
    table = table,
    L = L,
    rhs = rhs,
    requested_method = method,
    raw = parsed$raw
  )
  class(obj) <- "mm_contrast"
  obj
}

#' Fixed-effect bootstrap control
#'
#' @param nsim Requested bootstrap replicate count.
#' @param seed Optional integer seed. `NULL` leaves the Rust RNG seed
#'   unspecified and records that state in row details.
#' @param failed_refit_policy How failed refits are accounted for. Stable Rust
#'   wire labels are `"exclude"`, `"count_extreme"`, and `"abort"`.
#'
#' @return A list used by `contrast(..., method = "bootstrap")`.
#'
#' @export
bootstrap_control <- function(nsim = 999L,
                              seed = NULL,
                              failed_refit_policy = c("exclude", "count_extreme", "abort")) {
  failed_refit_policy <- match.arg(failed_refit_policy)
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    mm_abort(
      message = "`nsim` must be a positive integer.",
      class = "mm_inference_unavailable",
      input = nsim
    )
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0)) {
    mm_abort(
      message = "`seed` must be `NULL` or a non-negative integer.",
      class = "mm_inference_unavailable",
      input = seed
    )
  }
  structure(
    list(
      requested_replicates = as.integer(nsim),
      seed = if (is.null(seed)) NULL else as.integer(seed),
      failed_refit_policy = failed_refit_policy
    ),
    class = "mm_bootstrap_control"
  )
}

#' @method print mm_contrast
#' @export
print.mm_contrast <- function(x, ...) {
  cat("Fixed-effect contrasts:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Test a fixed-effect term
#'
#' `test_effect()` asks Rust to construct fixed-effect term hypotheses and
#' returns the corresponding fixed-effect inference rows.
#'
#' @param fit A fitted `mm_lmm`.
#' @param term A fixed-effect term label.
#' @param method Requested inference method.
#' @param bootstrap Optional [bootstrap_control()] object for bootstrap-backed
#'   methods.
#' @param group Optional grouping factor for `method = "cluster_bootstrap"`.
#'   Required for crossed or multi-grouping-factor models. In schema 1.0.0,
#'   cluster resampling is an estimator-distribution target and term-level
#'   p-values return `not_assessed` with a stable reason code.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_effect_test` object.
#'
#' @export
test_effect <- function(fit, term, method = c("auto", "satterthwaite",
                                             "kenward_roger", "bootstrap",
                                             "bootstrap_lrt",
                                             "cluster_bootstrap",
                                             "asymptotic", "none"),
                        bootstrap = NULL, group = NULL, ...) {
  UseMethod("test_effect")
}

#' @rdname test_effect
#' @export
test_effect.mm_lmm <- function(fit, term, method = c("auto", "satterthwaite",
                                                     "kenward_roger", "bootstrap",
                                                     "bootstrap_lrt",
                                                     "cluster_bootstrap",
                                                     "asymptotic", "none"),
                               bootstrap = NULL, group = NULL, ...) {
  method <- match.arg(method)
  if (!is.character(term) || !length(term)) {
    mm_abort(
      message = "`term` must be a non-empty character vector.",
      class = "mm_inference_unavailable",
      input = term
    )
  }
  fixed_terms <- mm_fixed_effect_terms(fit)
  unknown <- setdiff(term, fixed_terms)
  if (length(unknown)) {
    mm_abort(
      message = sprintf("Unknown fixed-effect term(s): %s.",
                        paste(unknown, collapse = ", ")),
      class = "mm_inference_unavailable",
      input = unknown
    )
  }
  if (identical(method, "none")) {
    table <- mm_unavailable_effect_table(term, method)
  } else if (identical(method, "bootstrap")) {
    rows <- lapply(term, function(t) {
      mm_rust_term_bootstrap_row(fit, t, bootstrap)
    })
    table <- do.call(rbind, rows)
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reliability_reason", "reason",
                       "reason_code", "reason_detail", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
  } else if (identical(method, "bootstrap_lrt")) {
    rows <- lapply(term, function(t) {
      mm_rust_term_bootstrap_lrt_row(fit, t, bootstrap)
    })
    table <- do.call(rbind, rows)
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reliability_reason", "reason",
                       "reason_code", "reason_detail", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
  } else if (identical(method, "cluster_bootstrap")) {
    table <- mm_cluster_bootstrap_unavailable_effect_table(fit, term, group)
  } else {
    parsed <- mm_rust_term_table(fit, method)
    table <- parsed$table[parsed$table$term %in% term, , drop = FALSE]
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reason", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
  }
  obj <- list(table = table, requested_method = method)
  class(obj) <- "mm_effect_test"
  obj
}

#' @method print mm_effect_test
#' @export
print.mm_effect_test <- function(x, ...) {
  cat("Effect tests:\n")
  cols <- intersect(
    c("term", "num_df", "den_df", "statistic", "statistic_name", "p_value",
      "method", "status", "reason_code"),
    names(x$table)
  )
  show <- x$table[, cols, drop = FALSE]
  show <- show[, vapply(show, function(col) {
    if (is.numeric(col)) return(any(!is.na(col)))
    values <- as.character(col)
    any(!is.na(values) & nzchar(values))
  }, logical(1)), drop = FALSE]
  print(show, row.names = FALSE)
  hidden <- setdiff(names(x$table), names(show))
  if (length(hidden)) {
    cat(sprintf("Full audit columns available in `x$table` (%d hidden).\n",
                length(hidden)))
  }
  invisible(x)
}

#' Assess contrast estimability
#'
#' Routes each requested contrast row through the Rust fixed-effect inference
#' bridge and reports the upstream estimability assessment verbatim.
#' Returned rows carry `status` (the closed enum from upstream:
#' `estimable`, `not_estimable`, `aliased`, ...), a boolean `estimable`
#' convenience flag, the contrast `rank` and `requested_rank`, and a
#' stable `reason` populated only when the engine refuses the contrast.
#'
#' @param fit A fitted `mm_lmm`.
#' @param L Optional contrast vector or matrix. Defaults to the fixed-effect
#'   coefficient basis.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_estimability` object.
#'
#' @export
estimability <- function(fit, L = NULL, ...) {
  UseMethod("estimability")
}

#' @rdname estimability
#' @export
estimability.mm_lmm <- function(fit, L = NULL, ...) {
  if (is.null(L)) {
    L <- diag(length(fit$beta))
    colnames(L) <- names(fit$beta)
    rownames(L) <- names(fit$beta)
  } else {
    L <- mm_contrast_matrix(L, fit)
  }

  parsed <- tryCatch(
    mm_rust_contrast_table(fit, L, rep(0, nrow(L)), "auto", bootstrap = NULL),
    error = function(cnd) cnd
  )

  if (inherits(parsed, "condition")) {
    table <- data.frame(
      contrast = rownames(L),
      estimable = NA,
      status = "not_assessed",
      rank = NA_integer_,
      requested_rank = NA_integer_,
      reason = conditionMessage(parsed),
      stringsAsFactors = FALSE
    )
  } else {
    table <- mm_estimability_table_from_inference(parsed$table, rownames(L))
  }

  obj <- list(table = table, L = L)
  class(obj) <- "mm_estimability"
  obj
}

mm_estimability_table_from_inference <- function(inference_table, contrast_names) {
  rows <- lapply(seq_len(nrow(inference_table)), function(i) {
    payload <- inference_table$estimability[[i]]
    assessment <- payload$assessment %||% list()
    status <- as.character(assessment$status %||% "not_assessed")
    diagnostics <- assessment$diagnostics %||% list()
    reason <- if (identical(status, "estimable")) {
      NA_character_
    } else if (length(diagnostics)) {
      paste(vapply(diagnostics, function(d) as.character(d$message %||% d), character(1)),
            collapse = "; ")
    } else if (identical(status, "not_assessed")) {
      "rust_estimability_certificate_unavailable"
    } else {
      sprintf("contrast %s under the fitted fixed-effect design", status)
    }
    data.frame(
      contrast = contrast_names[[i]],
      estimable = if (identical(status, "estimable")) TRUE else
                  if (identical(status, "not_assessed")) NA else FALSE,
      status = status,
      rank = as.integer(assessment$rank %||% NA_integer_),
      requested_rank = as.integer(assessment$requested_rank %||% NA_integer_),
      reason = reason,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @method print mm_estimability
#' @export
print.mm_estimability <- function(x, ...) {
  cat("Estimability:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Degrees of freedom for a contrast
#'
#' @param fit A fitted `mm_lmm`.
#' @param L A contrast vector or matrix.
#' @param method Requested degrees-of-freedom method.
#' @param ... Reserved for future methods.
#'
#' @return A numeric vector of `NA` degrees of freedom with an
#'   `mm_unavailable_reason` attribute.
#'
#' @export
df_for_contrast <- function(fit, L, method = c("auto", "satterthwaite",
                                              "kenward_roger", "bootstrap",
                                              "asymptotic", "none"), ...) {
  UseMethod("df_for_contrast")
}

#' @rdname df_for_contrast
#' @export
df_for_contrast.mm_lmm <- function(fit, L, method = c("auto", "satterthwaite",
                                                      "kenward_roger", "bootstrap",
                                                      "asymptotic", "none"), ...) {
  method <- match.arg(method)
  L <- mm_contrast_matrix(L, fit)
  n <- nrow(L)
  out <- rep(NA_real_, n)
  names(out) <- rownames(L)

  if (identical(method, "none")) {
    attr(out, "method") <- "not_requested"
    attr(out, "requested_method") <- method
    attr(out, "mm_unavailable_reason") <- mm_inference_unavailable_reason(method)
    class(out) <- c("mm_df_for_contrast", "numeric")
    return(out)
  }

  parsed <- tryCatch(
    mm_rust_contrast_table(fit, L, rep(0, n), method, bootstrap = NULL),
    error = function(cnd) cnd
  )
  if (inherits(parsed, "condition")) {
    attr(out, "method") <- "unavailable"
    attr(out, "requested_method") <- method
    attr(out, "mm_unavailable_reason") <- conditionMessage(parsed)
    class(out) <- c("mm_df_for_contrast", "numeric")
    return(out)
  }

  table <- parsed$table
  df_values <- if ("df" %in% names(table)) {
    as.numeric(table$df)
  } else if ("denominator_df" %in% names(table)) {
    as.numeric(table$denominator_df)
  } else {
    rep(NA_real_, nrow(table))
  }
  out[] <- df_values[seq_along(out)]

  resolved_methods <- unique(stats::na.omit(as.character(table$method)))
  method_label <- if (length(resolved_methods) == 1L) resolved_methods else method
  reasons <- stats::na.omit(as.character(table$reason))
  reason_label <- if (all(is.na(out))) {
    if (length(reasons)) reasons[[1L]] else mm_inference_unavailable_reason(method)
  } else {
    NA_character_
  }

  attr(out, "method") <- method_label
  attr(out, "requested_method") <- method
  if (!is.na(reason_label)) {
    attr(out, "mm_unavailable_reason") <- reason_label
  }
  class(out) <- c("mm_df_for_contrast", "numeric")
  out
}

#' @method print mm_df_for_contrast
#' @export
print.mm_df_for_contrast <- function(x, ...) {
  print(unclass(x))
  cat(sprintf("method: %s\n", attr(x, "method") %||% "unavailable"))
  cat(sprintf("reason: %s\n", attr(x, "mm_unavailable_reason") %||% "not_recorded"))
  invisible(x)
}

#' @method confint mm_lmm
#' @importFrom stats confint
#' @export
confint.mm_lmm <- function(object, parm, level = 0.95,
                           method = c("wald", "bootstrap"),
                           bootstrap = NULL,
                           interval = c("percentile", "basic"), ...) {
  method <- match.arg(method)
  interval <- match.arg(interval)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    mm_abort(
      message = "`level` must be a single number between 0 and 1.",
      class = "mm_inference_unavailable",
      input = level
    )
  }
  terms <- names(object$beta)
  if (missing(parm)) {
    parm <- terms
  }
  if (is.numeric(parm)) {
    parm <- terms[parm]
  }
  unknown <- setdiff(parm, terms)
  if (length(unknown)) {
    mm_abort(
      message = sprintf("Unknown fixed-effect parameter(s): %s.",
                        paste(unknown, collapse = ", ")),
      class = "mm_inference_unavailable",
      input = unknown
    )
  }
  if (identical(method, "bootstrap")) {
    return(mm_bootstrap_confint(object, parm, level, bootstrap, interval))
  }
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)
  est <- object$beta[parm]
  se <- object$std_errors[parm]
  out <- cbind(est - crit * se, est + crit * se)
  colnames(out) <- c(sprintf("%.1f %%", 100 * alpha / 2),
                     sprintf("%.1f %%", 100 * (1 - alpha / 2)))
  attr(out, "method") <- "wald_asymptotic_from_stored_standard_errors"
  attr(out, "status") <- "not_certified_by_rust_inference_contract"
  mm_new_confint(out)
}

#' @method print mm_confint
#' @export
print.mm_confint <- function(x, ...) {
  cat("Confidence intervals:\n")
  show <- matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )
  print(show, ...)

  method <- attr(x, "method") %||% "not_recorded"
  status <- attr(x, "status") %||% "not_recorded"
  cat(sprintf("method: %s\n", method))
  interval <- attr(x, "interval")
  if (!is.null(interval)) {
    cat(sprintf("interval: %s\n", interval))
  }
  cat(sprintf("status: %s\n", status))

  bootstrap <- attr(x, "bootstrap")
  if (length(bootstrap)) {
    summary <- mm_bootstrap_confint_summary(x, bootstrap)
    cat("\nBootstrap run:\n")
    print(summary, row.names = FALSE)
    notes <- unique(unlist(lapply(bootstrap, function(payload) {
      payload$metadata$notes %||% character()
    }), use.names = FALSE))
    if (length(notes)) {
      cat("notes:\n")
      for (note in notes) {
        wrapped <- strwrap(note, width = 78, exdent = 4)
        cat(sprintf("  - %s\n", paste(wrapped, collapse = "\n    ")))
      }
    }
    cat("Full bootstrap payload available in `attr(x, \"bootstrap\")`.\n")
  }
  invisible(x)
}

mm_contrast_matrix <- function(L, fit) {
  if (is.null(L)) {
    mm_abort(
      message = "`L` must be a numeric contrast vector or matrix.",
      class = "mm_inference_unavailable",
      input = L
    )
  }
  if (is.vector(L)) {
    L <- matrix(as.numeric(L), nrow = 1L)
  }
  if (!is.matrix(L) || !is.numeric(L)) {
    mm_abort(
      message = "`L` must be a numeric contrast vector or matrix.",
      class = "mm_inference_unavailable",
      input = L
    )
  }
  if (ncol(L) != length(fit$beta)) {
    mm_abort(
      message = sprintf("`L` must have %d column(s), one for each fixed effect.",
                        length(fit$beta)),
      class = "mm_inference_unavailable",
      input = L
    )
  }
  colnames(L) <- colnames(L) %||% names(fit$beta)
  if (is.null(rownames(L))) {
    rownames(L) <- paste0("c", seq_len(nrow(L)))
  }
  L
}

mm_rust_contrast_table <- function(fit, L, rhs, method, bootstrap = NULL) {
  bridge <- mm_rust_fit_bridge_payload(fit)
  if (identical(method, "bootstrap") && !is.null(bootstrap)) {
    if (!inherits(bootstrap, "mm_bootstrap_control")) {
      bootstrap <- do.call(bootstrap_control, as.list(bootstrap))
    }
    bootstrap_json <- jsonlite::toJSON(
      unclass(bootstrap),
      auto_unbox = TRUE,
      null = "null"
    )
    json <- tryCatch(
      mm_fixed_effect_bootstrap_contrast_json(
        bridge$formula_string,
        isTRUE(fit$REML),
        bridge$spec_data$column_order,
        bridge$spec_data$numeric_columns,
        bridge$spec_data$categorical_values,
        bridge$spec_data$categorical_levels,
        bridge$weights,
        bridge$control_json,
        as.numeric(t(L)),
        as.integer(nrow(L)),
        as.integer(ncol(L)),
        as.character(rownames(L)),
        as.numeric(rhs),
        as.character(bootstrap_json)
      ),
      error = function(cnd) cnd
    )
  } else {
    json <- tryCatch(
      mm_fixed_effect_contrast_json(
        bridge$formula_string,
        isTRUE(fit$REML),
        bridge$spec_data$column_order,
        bridge$spec_data$numeric_columns,
        bridge$spec_data$categorical_values,
        bridge$spec_data$categorical_levels,
        bridge$weights,
        bridge$control_json,
        as.numeric(t(L)),
        as.integer(nrow(L)),
        as.integer(ncol(L)),
        as.character(rownames(L)),
        as.numeric(rhs),
        method
      ),
      error = function(cnd) cnd
    )
  }
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = method)
  }
  table <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse fixed-effect contrast JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  mm_json_parse_fixed_effect_inference_table(table)
}

mm_rust_term_table <- function(fit, method) {
  bridge <- mm_rust_fit_bridge_payload(fit)
  json <- tryCatch(
    mm_fixed_effect_term_json(
      bridge$formula_string,
      isTRUE(fit$REML),
      bridge$spec_data$column_order,
      bridge$spec_data$numeric_columns,
      bridge$spec_data$categorical_values,
      bridge$spec_data$categorical_levels,
      bridge$weights,
      bridge$control_json,
      method
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = method)
  }
  table <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse fixed-effect term JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  mm_json_parse_fixed_effect_inference_table(table)
}

mm_rust_fit_bridge_payload <- function(fit) {
  spec_data <- mm_translate_data(fit$model_frame)
  formula_string <- mm_coerce_formula_string(fit$formula)
  control_json <- jsonlite::toJSON(
    unclass(fit$control %||% mm_control()),
    auto_unbox = TRUE,
    null = "null"
  )
  list(
    spec_data = spec_data,
    formula_string = formula_string,
    weights = mm_bridge_weights(fit$weights),
    control_json = as.character(control_json)
  )
}

mm_unavailable_effect_table <- function(term, method) {
  data.frame(
    term = term,
    num_df = NA_real_,
    den_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = "not_computed",
    requested_method = method,
    status = "not_assessed",
    reliability = "not_available",
    reason = mm_inference_unavailable_reason(method),
    details = I(rep(list(NULL), length(term))),
    notes = I(rep(list(character()), length(term))),
    stringsAsFactors = FALSE
  )
}

mm_cluster_bootstrap_unavailable_effect_table <- function(fit, term, group = NULL) {
  groups <- names(ranef(fit))
  if (is.null(group)) {
    if (length(groups) == 1L) {
      group <- groups[[1L]]
      reason_code <- "bootstrap_cluster_resample_p_value_unavailable"
      reason <- paste(
        "cluster_resample is an estimator-distribution target in schema 1.0.0;",
        "it does not certify fixed-effect hypothesis-test p-values"
      )
    } else {
      group <- NA_character_
      reason_code <- "cluster_bootstrap_multifactor_ambiguous"
      reason <- paste(
        "cluster bootstrap requires an explicit `group` when a model has",
        "multiple random-effect grouping factors"
      )
    }
  } else {
    if (!is.character(group) || length(group) != 1L || is.na(group) || !nzchar(group)) {
      mm_abort(
        message = "`group` must be a single grouping-factor name.",
        class = "mm_inference_unavailable",
        input = group
      )
    }
    if (!group %in% groups) {
      mm_abort(
        message = sprintf("Unknown random-effect grouping factor `%s`. Known groups: %s.",
                          group, paste(groups, collapse = ", ")),
        class = "mm_inference_unavailable",
        input = group
      )
    }
    reason_code <- "bootstrap_cluster_resample_p_value_unavailable"
    reason <- paste(
      "cluster_resample is an estimator-distribution target in schema 1.0.0;",
      "it does not certify fixed-effect hypothesis-test p-values"
    )
  }

  data.frame(
    term = term,
    num_df = NA_real_,
    den_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = "cluster_bootstrap",
    requested_method = "cluster_bootstrap",
    status = "not_assessed",
    reliability = "not_available",
    reliability_reason = "not_available",
    reason = reason,
    reason_code = reason_code,
    reason_detail = reason,
    details = I(rep(list(list(
      bootstrap = list(
        target_kind = "cluster_resample",
        group = group,
        p_value_certified = FALSE
      )
    )), length(term))),
    notes = I(rep(list(character()), length(term))),
    stringsAsFactors = FALSE
  )
}

mm_fixed_effect_terms <- function(fit) {
  terms <- as.character(unlist(
    fit$artifact$semantic_model$fixed_terms %||% list(),
    use.names = FALSE
  ))
  if (!length(terms)) {
    terms <- names(fit$beta)
  }
  terms
}

mm_inference_unavailable_reason <- function(method) {
  if (identical(method, "none")) {
    return("inference_not_requested")
  }
  "rust_inference_certificate_unavailable"
}

mm_bootstrap_confint <- function(fit, parm, level, bootstrap, interval) {
  if (is.null(bootstrap)) bootstrap <- bootstrap_control()
  if (!inherits(bootstrap, "mm_bootstrap_control")) {
    bootstrap <- do.call(bootstrap_control, as.list(bootstrap))
  }
  rows <- lapply(parm, function(parameter) {
    payload <- mm_full_model_bootstrap_payload(fit, parameter, level, bootstrap)
    interval_row <- mm_select_bootstrap_interval(payload, level, interval)
    list(
      parameter = parameter,
      lower = interval_row$lower,
      upper = interval_row$upper,
      payload = payload
    )
  })

  out <- cbind(
    vapply(rows, `[[`, numeric(1), "lower"),
    vapply(rows, `[[`, numeric(1), "upper")
  )
  rownames(out) <- vapply(rows, `[[`, character(1), "parameter")
  alpha <- 1 - level
  colnames(out) <- c(sprintf("%.1f %%", 100 * alpha / 2),
                     sprintf("%.1f %%", 100 * (1 - alpha / 2)))
  attr(out, "method") <- "bootstrap_full_model_distribution"
  attr(out, "interval") <- interval
  attr(out, "status") <- "available"
  attr(out, "bootstrap") <- lapply(rows, `[[`, "payload")
  mm_new_confint(out)
}

mm_new_confint <- function(x) {
  class(x) <- c("mm_confint", "matrix")
  x
}

mm_bootstrap_confint_summary <- function(x, bootstrap) {
  parameters <- rownames(x)
  if (is.null(parameters)) {
    parameters <- rep(NA_character_, length(bootstrap))
  }
  rows <- lapply(seq_along(bootstrap), function(i) {
    metadata <- bootstrap[[i]]$metadata %||% list()
    seed <- metadata$seed_record$seed %||% NA
    data.frame(
      parameter = parameters[[i]],
      requested = as.integer(metadata$requested_replicates %||% NA_integer_),
      successful = as.integer(metadata$successful_replicates %||% NA_integer_),
      failed_refits = as.integer(metadata$failed_refits %||% NA_integer_),
      boundary_rate = round(as.numeric(metadata$boundary_rate %||% NA_real_), 3),
      seed = as.integer(seed),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mm_full_model_bootstrap_payload <- function(fit, parameter, level, bootstrap) {
  L <- matrix(0, nrow = 1L, ncol = length(fit$beta))
  colnames(L) <- names(fit$beta)
  rownames(L) <- parameter
  L[1L, parameter] <- 1
  bridge <- mm_rust_fit_bridge_payload(fit)
  bootstrap_json <- jsonlite::toJSON(
    unclass(bootstrap),
    auto_unbox = TRUE,
    null = "null"
  )
  json <- tryCatch(
    mm_full_model_bootstrap_contrast_json(
      bridge$formula_string,
      isTRUE(fit$REML),
      bridge$spec_data$column_order,
      bridge$spec_data$numeric_columns,
      bridge$spec_data$categorical_values,
      bridge$spec_data$categorical_levels,
      bridge$weights,
      bridge$control_json,
      as.numeric(t(L)),
      as.integer(nrow(L)),
      as.integer(ncol(L)),
      as.character(rownames(L)),
      0,
      as.character(bootstrap_json),
      as.numeric(level)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = "bootstrap")
  }
  tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse full-model bootstrap JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
}

mm_select_bootstrap_interval <- function(payload, level, interval) {
  intervals <- payload$intervals %||% list()
  hits <- vapply(intervals, function(x) {
    identical(mm_scalar_text(x$method), interval) &&
      isTRUE(all.equal(as.numeric(x$level), level, tolerance = 1e-12))
  }, logical(1))
  if (!any(hits)) {
    mm_abort(
      message = sprintf("Bootstrap payload did not contain a %s interval at level %.3f.",
                        interval, level),
      class = "mm_schema_error",
      input = payload
    )
  }
  row <- intervals[[which(hits)[[1L]]]]
  list(
    lower = as.numeric(row$lower),
    upper = as.numeric(row$upper)
  )
}

# Map a term name to a row-permutation L matrix that picks the coefficients
# belonging to the term. The intercept term `1` matches `(Intercept)` exactly;
# any other term `t` matches coefficient names of the form `t` or `t: <level>`
# (mixeff's factor-level labelling) and `t:other` (interaction prefix).
mm_term_to_l_matrix <- function(fit, term) {
  coef_names <- names(fit$beta)
  k <- length(coef_names)
  hits <- if (identical(term, "1")) {
    coef_names == "(Intercept)"
  } else {
    pattern <- paste0("^", regex_escape(term), "(: |$|:)")
    grepl(pattern, coef_names)
  }
  if (!any(hits)) {
    mm_abort(
      message = sprintf(
        "Term '%s' did not match any fixed-effect coefficient. Coefficients: %s",
        term, paste(coef_names, collapse = ", ")),
      class = "mm_inference_unavailable",
      input = term
    )
  }
  idx <- which(hits)
  L <- matrix(0, nrow = length(idx), ncol = k)
  for (i in seq_along(idx)) L[i, idx[[i]]] <- 1
  rownames(L) <- coef_names[idx]
  colnames(L) <- coef_names
  L
}

regex_escape <- function(s) {
  gsub("([.\\\\+*?\\[\\^\\]$(){}=!<>|:#-])", "\\\\\\1", s, perl = TRUE)
}

# Dispatch a term-level bootstrap test through the new Rust bridge entrypoint.
# Single-df terms produce a t-form row; multi-df terms produce an F-form row
# with `numerator_df` set to the effective restriction rank.
mm_rust_term_bootstrap_row <- function(fit, term, bootstrap) {
  if (is.null(bootstrap)) bootstrap <- bootstrap_control()
  if (!inherits(bootstrap, "mm_bootstrap_control")) {
    bootstrap <- do.call(bootstrap_control, as.list(bootstrap))
  }
  L <- mm_term_to_l_matrix(fit, term)
  rhs <- rep(0, nrow(L))
  bridge <- mm_rust_fit_bridge_payload(fit)
  bootstrap_json <- jsonlite::toJSON(
    unclass(bootstrap),
    auto_unbox = TRUE,
    null = "null"
  )
  json <- tryCatch(
    mm_fixed_effect_bootstrap_term_json(
      bridge$formula_string,
      isTRUE(fit$REML),
      bridge$spec_data$column_order,
      bridge$spec_data$numeric_columns,
      bridge$spec_data$categorical_values,
      bridge$spec_data$categorical_levels,
      bridge$weights,
      bridge$control_json,
      as.numeric(t(L)),
      as.integer(nrow(L)),
      as.integer(ncol(L)),
      as.character(term),
      as.numeric(rhs),
      as.character(bootstrap_json)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = "bootstrap")
  }
  table <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  parsed <- mm_json_parse_fixed_effect_inference_table(table)
  row <- parsed$table[1L, , drop = FALSE]
  row$term <- term
  row
}

# Dispatch a bootstrap likelihood-ratio test for `term`. Builds the reduced
# formula (drops the term from the fixed effects), refuses REML with a stable
# reason, and routes through the new `mm_bootstrap_lrt_json` bridge entrypoint.
# Returns a single inference-table-shaped row.
mm_rust_term_bootstrap_lrt_row <- function(fit, term, bootstrap) {
  if (isTRUE(fit$REML)) {
    return(mm_inference_row_unavailable(
      term = term,
      method = "bootstrap_lrt",
      reason = paste(
        "bootstrap likelihood-ratio test requires ML-fitted models;",
        "this fit is REML -- refit with `lmm(..., reml = FALSE)` and retry"
      ),
      reason_code = "bootstrap_lrt_requires_ml"
    ))
  }
  if (is.null(bootstrap)) bootstrap <- bootstrap_control()
  if (!inherits(bootstrap, "mm_bootstrap_control")) {
    bootstrap <- do.call(bootstrap_control, as.list(bootstrap))
  }
  reduced_formula <- tryCatch(
    update(fit$formula, stats::as.formula(sprintf(". ~ . - %s", term))),
    error = function(e) e
  )
  if (inherits(reduced_formula, "condition")) {
    return(mm_inference_row_unavailable(
      term = term,
      method = "bootstrap_lrt",
      reason = sprintf("could not drop term '%s' from formula: %s",
                       term, conditionMessage(reduced_formula)),
      reason_code = "bootstrap_lrt_reduced_formula_failed"
    ))
  }
  bridge <- mm_rust_fit_bridge_payload(fit)
  bootstrap_json <- jsonlite::toJSON(
    unclass(bootstrap),
    auto_unbox = TRUE,
    null = "null"
  )
  json <- tryCatch(
    mm_bootstrap_lrt_json(
      deparse1(reduced_formula),
      bridge$formula_string,
      bridge$spec_data$column_order,
      bridge$spec_data$numeric_columns,
      bridge$spec_data$categorical_values,
      bridge$spec_data$categorical_levels,
      bridge$weights,
      bridge$control_json,
      as.character(bootstrap_json)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    return(mm_inference_row_unavailable(
      term = term,
      method = "bootstrap_lrt",
      reason = conditionMessage(json),
      reason_code = "bootstrap_lrt_engine_refused"
    ))
  }
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  payload <- parsed$payload %||% list()
  data.frame(
    term = term,
    label = term,
    kind = "coefficient",
    estimate = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    numerator_df = NA_real_,
    denominator_df = NA_real_,
    statistic = parsed$observed_statistic %||% NA_real_,
    statistic_name = "chi_square",
    p_value = parsed$p_value %||% NA_real_,
    method = "bootstrap_lrt",
    status = if (is.null(parsed$p_value)) "not_assessed" else "available",
    reliability = if (is.null(parsed$p_value)) "not_available" else "moderate",
    reliability_reason = "bootstrap_monte_carlo_replicates",
    reason = NA_character_,
    reason_code = NA_character_,
    reason_detail = NA_character_,
    estimability = I(list(NULL)),
    details = I(list(list(
      bootstrap = list(
        target_kind = payload$metadata$target$kind %||%
          payload$target_kind %||% "likelihood_ratio",
        observed_statistic = parsed$observed_statistic %||% NA_real_,
        mcse = parsed$mcse %||% NA_real_,
        successful_replicates = payload$metadata$successful_replicates %||% NA_integer_,
        completed_replicates = payload$metadata$completed_replicates %||% NA_integer_,
        replicate_statistics = payload$replicate_statistics %||% list(),
        boundary_count = payload$metadata$boundary_count %||% NA_integer_,
        seed = payload$metadata$seed_record$seed %||% NA_integer_
      ),
      reduced_formula = deparse1(reduced_formula),
      alternative_formula = bridge$formula_string
    ))),
    notes = I(list(as.character(unlist(parsed$notes %||% list(),
                                       use.names = FALSE)))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# Construct a one-row "unavailable" inference-table-shaped data frame.
mm_inference_row_unavailable <- function(term, method, reason, reason_code = NA_character_) {
  data.frame(
    term = term,
    label = term,
    kind = "coefficient",
    estimate = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    numerator_df = NA_real_,
    denominator_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = method,
    status = "not_assessed",
    reliability = "not_available",
    reliability_reason = "inference_unavailable_by_policy",
    reason = reason,
    reason_code = reason_code,
    reason_detail = NA_character_,
    estimability = I(list(NULL)),
    details = I(list(NULL)),
    notes = I(list(character())),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
