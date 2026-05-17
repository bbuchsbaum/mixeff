#' Compare fitted mixeff models
#'
#' `compare()` is the namespace-qualified model-comparison front door. For LMMs
#' it reports likelihood, information criteria, and asymptotic likelihood-ratio
#' comparisons. REML fits are refit by ML when `refit_for_comparison = "auto"`
#' or `"ml"`; `"error"` refuses that comparison.
#'
#' @param object A fitted `mm_lmm`.
#' @param ... Additional fitted `mm_lmm` objects.
#' @param target Comparison target label.
#' @param method `"auto"` / `"lrt"` for asymptotic likelihood-ratio rows,
#'   `"aic"` for information criteria only, or `"bootstrap"` for a small
#'   parametric-bootstrap LRT when `nsim > 0`.
#' @param refit_for_comparison How to handle REML fits.
#' @param nsim Number of bootstrap simulations for `method = "bootstrap"`.
#' @param seed Optional bootstrap seed.
#'
#' @return An `mm_model_comparison` object with a data-frame `table`.
#'
#' @importFrom stats anova drop1
#' @export
compare <- function(object, ...) {
  UseMethod("compare")
}

#' @rdname compare
#' @export
compare.mm_lmm <- function(object,
                           ...,
                           target = c("fixed_effects", "random_effects", "prediction"),
                           method = c("auto", "lrt", "bootstrap", "aic"),
                           refit_for_comparison = c("auto", "error", "ml"),
                           nsim = 0L,
                           seed = NULL) {
  target <- match.arg(target)
  method <- match.arg(method)
  refit_for_comparison <- match.arg(refit_for_comparison)
  if (identical(method, "bootstrap")) {
    if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 0) {
      mm_abort(
        message = "`nsim` must be a non-negative integer for `compare(method = \"bootstrap\")`.",
        class = "mm_arg_error",
        input = nsim
      )
    }
    nsim <- as.integer(nsim)
  }
  fits <- c(list(object), list(...))
  if (!all(vapply(fits, inherits, logical(1), what = "mm_lmm"))) {
    mm_abort(
      message = "`compare()` requires fitted `mm_lmm` objects.",
      class = "mm_arg_error",
      input = fits
    )
  }
  mm_assert_comparable_lmm(fits)
  prepared <- mm_prepare_comparison_fits(fits, refit_for_comparison)
  fits <- prepared$fits
  table <- mm_compare_table(fits, method, prepared$refit)
  bootstrap <- NULL
  if (identical(method, "bootstrap") && nsim > 0L && length(fits) == 2L) {
    bootstrap <- parametric_bootstrap(
      fits[[1L]],
      fits[[2L]],
      nsim = nsim,
      seed = seed
    )
    last <- nrow(table)
    table$p_value[last] <- bootstrap$p_value
    table$method[last] <- "parametric_bootstrap_lrt"
    table$status[last] <- bootstrap$status
    succ <- bootstrap$successful_replicates
    table$reason[last] <- if (identical(bootstrap$status, "available")) {
      sprintf("parametric bootstrap LRT (%s/%d replicates%s)",
              if (is.na(succ)) "?" else as.character(succ),
              nsim,
              if (is.na(bootstrap$mcse)) "" else sprintf(", MCSE=%.4g", bootstrap$mcse))
    } else {
      bootstrap$reason %||% "parametric bootstrap LRT did not certify a p-value"
    }
  } else if (identical(method, "bootstrap")) {
    # Bootstrap was requested but cannot run: do not let the asymptotic LRT
    # p-value masquerade as the requested bootstrap result.
    table$p_value <- NA_real_
    table$method <- "bootstrap_not_run"
    table$status <- "bootstrap_not_run"
    table$reason <- "set nsim > 0 and compare exactly two models to run bootstrap"
  }
  obj <- list(
    table = table,
    fits = fits,
    target = target,
    method = method,
    refit_for_comparison = refit_for_comparison,
    bootstrap = bootstrap
  )
  class(obj) <- "mm_model_comparison"
  obj
}

#' @method print mm_model_comparison
#' @export
print.mm_model_comparison <- function(x, ...) {
  cat("Model comparison:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Parametric bootstrap likelihood-ratio comparison
#'
#' Runs the engine-certified parametric-bootstrap likelihood-ratio test
#' between two nested ML-fitted LMMs through the Rust
#' `mm_bootstrap_lrt_json` entry point. The smaller model (fewer estimated
#' parameters) is the reduced model; the larger is the alternative. The
#' returned object carries the engine's replicate accounting (successful and
#' completed replicates, boundary count, Monte-Carlo standard error, seed)
#' rather than a bare `mean()` p-value, so every reported number traces back
#' to a versioned Rust payload.
#'
#' The engine refuses REML fits: refit with `lmm(..., REML = FALSE)` before
#' calling. (`compare(method = "bootstrap")` refits REML to ML automatically.)
#'
#' @param null,alternative Fitted `mm_lmm` objects. Order is irrelevant; the
#'   model with fewer parameters is treated as the reduced model.
#' @param nsim Number of bootstrap replicates.
#' @param seed Optional bootstrap seed.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_parametric_bootstrap` object.
#'
#' @export
parametric_bootstrap <- function(null, alternative, nsim = 100L, seed = NULL, ...) {
  if (!inherits(null, "mm_lmm") || !inherits(alternative, "mm_lmm")) {
    mm_abort(
      message = "`parametric_bootstrap()` requires two fitted `mm_lmm` objects.",
      class = "mm_arg_error",
      input = list(null = null, alternative = alternative)
    )
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    mm_abort(
      message = "`nsim` must be a positive integer.",
      class = "mm_arg_error",
      input = nsim
    )
  }
  nsim <- as.integer(nsim)
  # Reduced = fewer estimated parameters; the engine LRT is direction-aware
  # so we order explicitly rather than trust the call order.
  if (alternative$dof < null$dof) {
    tmp <- null
    null <- alternative
    alternative <- tmp
  }
  if (isTRUE(null$REML) || isTRUE(alternative$REML)) {
    mm_abort(
      message = paste(
        "parametric bootstrap likelihood-ratio test requires ML-fitted",
        "models; refit with `lmm(..., REML = FALSE)` and retry"
      ),
      class = "mm_inference_unavailable",
      reason_code = "bootstrap_lrt_requires_ml",
      input = list(null_reml = isTRUE(null$REML),
                   alternative_reml = isTRUE(alternative$REML))
    )
  }
  mm_assert_bootstrap_lrt_pair(null, alternative)
  bootstrap <- bootstrap_control(nsim = nsim, seed = seed)
  bridge <- mm_rust_fit_bridge_payload(alternative)
  bootstrap_json <- jsonlite::toJSON(
    unclass(bootstrap),
    auto_unbox = TRUE,
    null = "null"
  )
  json <- tryCatch(
    mm_bootstrap_lrt_json(
      deparse1(null$formula),
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
    mm_abort(
      message = conditionMessage(json),
      class = "mm_inference_unavailable",
      reason_code = "bootstrap_lrt_engine_refused",
      parent = json
    )
  }
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  payload <- parsed$payload %||% list()
  meta <- payload$metadata %||% list()
  simulated <- as.numeric(unlist(payload$replicate_statistics %||% list(),
                                 use.names = FALSE))
  certified <- !is.null(parsed$p_value)
  out <- list(
    observed = parsed$observed_statistic %||% mm_lrt_stat(null, alternative),
    simulated = simulated,
    p_value = parsed$p_value %||% NA_real_,
    nsim = nsim,
    successful_replicates = meta$successful_replicates %||% NA_integer_,
    completed_replicates = meta$completed_replicates %||% NA_integer_,
    boundary_count = meta$boundary_count %||% NA_integer_,
    mcse = parsed$mcse %||% NA_real_,
    seed = meta$seed_record$seed %||% seed,
    status = if (certified) "available" else "not_assessed",
    reason = if (certified) {
      NA_character_
    } else {
      "parametric bootstrap likelihood-ratio test did not certify a p-value"
    },
    notes = as.character(unlist(parsed$notes %||% list(), use.names = FALSE)),
    reduced_formula = deparse1(null$formula),
    alternative_formula = bridge$formula_string
  )
  class(out) <- "mm_parametric_bootstrap"
  out
}

#' @method print mm_parametric_bootstrap
#' @export
print.mm_parametric_bootstrap <- function(x, ...) {
  fmt_int <- function(v) if (is.null(v) || is.na(v)) "NA" else format(as.integer(v))
  fmt_num <- function(v) if (is.null(v) || !is.finite(v)) "NA" else sprintf("%.6g", v)
  cat("Parametric bootstrap likelihood-ratio test:\n")
  cat(sprintf("  status:    %s\n", x$status %||% "unknown"))
  cat(sprintf("  observed:  %s\n", fmt_num(x$observed)))
  cat(sprintf("  requested replicates: %s\n", fmt_int(x$nsim)))
  cat(sprintf("  successful / completed: %s / %s\n",
              fmt_int(x$successful_replicates),
              fmt_int(x$completed_replicates)))
  if (!is.null(x$boundary_count) && !is.na(x$boundary_count) &&
      x$boundary_count > 0L) {
    cat(sprintf("  boundary replicates: %s\n", fmt_int(x$boundary_count)))
  }
  cat(sprintf("  MCSE:      %s\n", fmt_num(x$mcse)))
  cat(sprintf("  seed:      %s\n", fmt_int(x$seed)))
  if (identical(x$status, "available")) {
    cat(sprintf("  p.value:   %s\n", fmt_num(x$p_value)))
  } else {
    cat(sprintf("  p.value:   not certified -- %s\n",
                x$reason %||% "no reason recorded"))
  }
  notes <- x$notes %||% character()
  if (length(notes)) {
    cat("  notes:\n")
    for (n in notes) cat(sprintf("    - %s\n", n))
  }
  invisible(x)
}

#' @method anova mm_lmm
#' @export
anova.mm_lmm <- function(object, ..., type = c("III", "II", "I"),
                         method = c("auto", "satterthwaite", "kenward_roger",
                                    "bootstrap", "asymptotic", "none"),
                         refit_for_comparison = c("auto", "error", "ml")) {
  dots <- list(...)
  if (length(dots)) {
    cmp_method <- if (identical(match.arg(method), "bootstrap")) "bootstrap" else "auto"
    return(compare(
      object,
      ...,
      method = cmp_method,
      refit_for_comparison = match.arg(refit_for_comparison)
    ))
  }
  type <- match.arg(type)
  method <- match.arg(method)
  refit_for_comparison <- match.arg(refit_for_comparison)
  terms <- setdiff(mm_fixed_effect_terms(object), "1")
  if (identical(method, "none")) {
    table <- mm_unavailable_effect_table(terms, method)
  } else {
    parsed <- mm_rust_term_table(object, method)
    table <- parsed$table[parsed$table$term %in% terms, , drop = FALSE]
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reason", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
  }
  table$type <- type
  table <- table[, c("term", "type", setdiff(names(table), c("term", "type"))),
                 drop = FALSE]
  obj <- list(
    table = table,
    type = type,
    requested_method = method,
    refit_for_comparison = refit_for_comparison
  )
  class(obj) <- "mm_anova"
  obj
}

#' @method print mm_anova
#' @export
print.mm_anova <- function(x, ...) {
  cat(sprintf("Type %s analysis of fixed effects (method: %s):\n",
              x$type %||% "III", x$requested_method %||% "auto"))
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Drop one fixed-effect term at a time
#'
#' `drop1.mm_lmm()` refits reduced fixed-effect models and compares them to the
#' original fit. It is conservative: random-effect terms are preserved exactly,
#' and the reduced formulas are reported in the result table.
#'
#' @param object A fitted `mm_lmm`.
#' @param scope Optional character vector of fixed-effect terms to drop.
#' @param test Comparison test label. `"Chisq"` reports asymptotic LRT rows;
#'   `"none"` reports information criteria only.
#' @param refit_for_comparison How to handle REML fits.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_drop1` object.
#'
#' @method drop1 mm_lmm
#' @export
drop1.mm_lmm <- function(object,
                         scope = NULL,
                         test = c("none", "Chisq"),
                         refit_for_comparison = c("auto", "error", "ml"),
                         ...) {
  test <- match.arg(test)
  refit_for_comparison <- match.arg(refit_for_comparison)
  terms <- setdiff(mm_fixed_effect_terms(object), "1")
  if (!is.null(scope)) {
    terms <- intersect(terms, as.character(scope))
  }
  full <- mm_prepare_comparison_fits(list(object), refit_for_comparison)$fits[[1L]]
  rows <- lapply(terms, function(term) {
    reduced_formula <- mm_drop_fixed_term_formula(full, term)
    reduced <- lmm(reduced_formula, full$model_frame, REML = isTRUE(full$REML),
                   weights = full$weights,
                   control = mm_control(verbose = -1))
    stat <- mm_lrt_stat(reduced, full)
    df <- full$dof - reduced$dof
    data.frame(
      dropped = term,
      formula = deparse1(reduced_formula),
      df = df,
      logLik = as.numeric(logLik(reduced)),
      AIC = AIC(reduced),
      BIC = BIC(reduced),
      LRT = if (identical(test, "Chisq")) stat else NA_real_,
      p_value = if (identical(test, "Chisq") && df > 0) {
        stats::pchisq(stat, df = df, lower.tail = FALSE)
      } else {
        NA_real_
      },
      method = if (identical(test, "Chisq")) "asymptotic_lrt" else "none",
      stringsAsFactors = FALSE
    )
  })
  table <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      dropped = character(),
      formula = character(),
      df = numeric(),
      logLik = numeric(),
      AIC = numeric(),
      BIC = numeric(),
      LRT = numeric(),
      p_value = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    )
  }
  rownames(table) <- NULL
  obj <- list(table = table, full = full)
  class(obj) <- "mm_drop1"
  obj
}

#' @method print mm_drop1
#' @export
print.mm_drop1 <- function(x, ...) {
  cat("Single-term deletion table:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

mm_assert_comparable_lmm <- function(fits) {
  n <- vapply(fits, nobs, integer(1))
  if (length(unique(n)) != 1L) {
    mm_abort(
      message = "Compared models must have the same number of observations.",
      class = "mm_arg_error",
      input = n
    )
  }
  responses <- vapply(fits, mm_response_name, character(1))
  if (length(unique(responses)) != 1L) {
    mm_abort(
      message = "Compared models must use the same response variable.",
      class = "mm_arg_error",
      input = responses
    )
  }
}

mm_assert_bootstrap_lrt_pair <- function(null, alternative) {
  if (!isTRUE(alternative$dof > null$dof)) {
    mm_abort(
      message = paste(
        "Parametric bootstrap LRT requires nested models with the",
        "alternative estimating more parameters than the reduced model."
      ),
      class = "mm_arg_error",
      reason_code = "bootstrap_lrt_requires_nested_models",
      input = list(null_df = null$dof, alternative_df = alternative$dof)
    )
  }

  null_vars <- all.vars(null$formula)
  alt_names <- names(alternative$model_frame %||% data.frame())
  missing <- setdiff(null_vars, alt_names)
  if (length(missing)) {
    mm_abort(
      message = sprintf(
        "Parametric bootstrap LRT requires the reduced model variables to be present in the alternative model frame; missing: %s.",
        paste(missing, collapse = ", ")
      ),
      class = "mm_arg_error",
      reason_code = "bootstrap_lrt_requires_nested_model_frames",
      input = missing
    )
  }

  shared <- intersect(names(null$model_frame %||% data.frame()), alt_names)
  mismatched <- shared[!vapply(shared, function(nm) {
    identical(null$model_frame[[nm]], alternative$model_frame[[nm]])
  }, logical(1))]
  if (length(mismatched)) {
    mm_abort(
      message = sprintf(
        "Parametric bootstrap LRT requires compared fits to share identical model-frame values; mismatched column(s): %s.",
        paste(mismatched, collapse = ", ")
      ),
      class = "mm_arg_error",
      reason_code = "bootstrap_lrt_requires_same_observations",
      input = mismatched
    )
  }

  if (!identical(null$weights, alternative$weights)) {
    mm_abort(
      message = "Parametric bootstrap LRT requires compared fits to use identical case weights.",
      class = "mm_arg_error",
      reason_code = "bootstrap_lrt_requires_same_weights",
      input = list(null = null$weights, alternative = alternative$weights)
    )
  }
  invisible(TRUE)
}

mm_prepare_comparison_fits <- function(fits, refit_for_comparison) {
  has_reml <- vapply(fits, function(x) isTRUE(x$REML), logical(1))
  refit <- rep(FALSE, length(fits))
  if (any(has_reml)) {
    if (identical(refit_for_comparison, "error")) {
      mm_abort(
        message = "REML fits require `refit_for_comparison = \"auto\"` or `\"ml\"` for likelihood comparison.",
        class = "mm_inference_unavailable",
        input = has_reml
      )
    }
    fits <- lapply(fits, function(fit) {
      if (!isTRUE(fit$REML)) return(fit)
      lmm(fit$formula, fit$model_frame, REML = FALSE,
          weights = fit$weights,
          control = mm_control(verbose = -1))
    })
    refit <- has_reml
  }
  list(fits = fits, refit = refit)
}

mm_compare_table <- function(fits, method, refit) {
  ord <- order(vapply(fits, function(x) x$dof, numeric(1)))
  fits <- fits[ord]
  refit <- refit[ord]
  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    data.frame(
      model = paste0("m", i),
      formula = deparse1(fit$formula),
      nobs = nobs(fit),
      df = fit$dof,
      logLik = as.numeric(logLik(fit)),
      deviance = deviance(fit),
      AIC = AIC(fit),
      BIC = BIC(fit),
      REML = isTRUE(fit$REML),
      refit = isTRUE(refit[[i]]),
      delta_df = NA_real_,
      LRT = NA_real_,
      p_value = NA_real_,
      method = if (identical(method, "aic")) "none" else "asymptotic_lrt",
      status = if (identical(method, "aic")) "information_criteria" else "asymptotic_uncertified",
      reason = if (identical(method, "aic")) "" else "ordinary LRT is not a finite-sample mixed-model certificate",
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  if (!identical(method, "aic") && nrow(table) > 1L) {
    for (i in 2:nrow(table)) {
      table$delta_df[[i]] <- table$df[[i]] - table$df[[i - 1L]]
      table$LRT[[i]] <- pmax(0, table$deviance[[i - 1L]] - table$deviance[[i]])
      if (table$delta_df[[i]] > 0) {
        table$p_value[[i]] <- stats::pchisq(table$LRT[[i]],
                                            df = table$delta_df[[i]],
                                            lower.tail = FALSE)
      }
    }
  }
  rownames(table) <- NULL
  table
}

mm_lrt_stat <- function(null, alternative) {
  pmax(0, deviance(null) - deviance(alternative))
}

mm_drop_fixed_term_formula <- function(fit, term) {
  response <- mm_response_name(fit)
  fixed <- setdiff(mm_fixed_effect_terms(fit), c("1", term))
  fixed_rhs <- if (length(fixed)) paste(fixed, collapse = " + ") else "1"
  random <- vapply(
    fit$artifact$semantic_model$random_terms %||% list(),
    function(x) x$source_syntax$text %||% "",
    character(1)
  )
  random <- random[nzchar(random)]
  rhs <- paste(c(fixed_rhs, random), collapse = " + ")
  stats::as.formula(paste(response, "~", rhs), env = environment(fit$formula))
}
