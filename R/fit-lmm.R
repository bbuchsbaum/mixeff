#' Fit a linear mixed-effects model
#'
#' `lmm()` is mixeff's Phase 1 linear mixed-model fit driver. It compiles the
#' requested lme4-style formula, prints the same [explain_model()] view that
#' pre-fit audit users see, then delegates the numerical fit to the upstream
#' Rust `LinearMixedModel`.
#'
#' The returned object is deliberately serializable: fixed effects, theta,
#' sigma, likelihood summaries, fitted values, residuals, random effects, and
#' the post-fit compiler artifact are all stored directly on the R object. The
#' native Rust handle is treated as a rebuildable cache, not as the source of
#' truth.
#'
#' @param formula A two-sided lme4-style formula, e.g.
#'   `y ~ x + (1 + x | subject)`.
#' @param data A `data.frame` containing all variables in `formula`.
#' @param REML Logical; fit by restricted maximum likelihood when `TRUE`.
#' @param weights Optional positive numeric case weights, either a vector with
#'   one value per row or an expression evaluated in `data`.
#' @param control A list from [mm_control()].
#'
#' @return An object of class `mm_lmm`, also inheriting from `mm_fit` and
#'   `mm_compiled`.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   y = rnorm(40),
#'   x = rep(0:3, 10),
#'   subject = factor(rep(seq_len(10), each = 4))
#' )
#' fit <- lmm(y ~ x + (1 + x | subject), df, control = mm_control(verbose = -1))
#' fixef(fit)
#' VarCorr(fit)
#' }
#'
#' @export
lmm <- function(formula, data, REML = TRUE, weights = NULL,
                control = mm_control()) {
  call <- match.call()
  control <- mm_validate_control(control)
  weights <- mm_lmm_weights(substitute(weights), data, parent.frame())
  if (!is.logical(REML) || length(REML) != 1L || is.na(REML)) {
    mm_abort(
      message = "`REML` must be TRUE or FALSE.",
      class = "mm_fit_error",
      input = REML
    )
  }

  spec <- compile_model(formula, data)
  mm_validate_fit_structure(spec)
  if (control$verbose >= 0L) {
    print(explain_model(spec))
  }

  spec_data <- mm_translate_data(spec$model_frame)
  formula_string <- mm_coerce_formula_string(formula)
  control_json <- jsonlite::toJSON(unclass(control), auto_unbox = TRUE, null = "null")

  json <- tryCatch(
    .Call(
      wrap__mm_fit_lmm_json,
      formula_string,
      isTRUE(REML),
      spec_data$column_order,
      spec_data$numeric_columns,
      spec_data$categorical_values,
      spec_data$categorical_levels,
      mm_bridge_weights(weights),
      as.character(control_json)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, formula = formula_string)
  }

  fit_result <- mm_json_parse_lmm_fit(json)
  fit_summary <- mm_json_parse_fit_summary(fit_result$fit_summary)
  artifact <- mm_json_parse_artifact(fit_result$artifact_json)
  beta <- mm_named_numeric(fit_result$beta, fit_result$beta_names)
  std_errors <- mm_named_numeric(fit_result$std_errors, fit_result$beta_names)
  fixed_effect_vcov <- mm_fixed_effect_vcov_from_payload(
    artifact$fixed_effect_covariance_matrix,
    beta,
    std_errors
  )

  fit <- list(
    call           = call,
    formula        = formula,
    REML           = isTRUE(REML),
    control        = control,
    vars           = spec$vars,
    model_frame    = spec$model_frame,
    weights        = weights,
    artifact       = artifact,
    fit            = fit_result,
    fit_summary    = fit_summary,
    schema         = mm_object_schema(artifact),
    rust_handle    = NULL,
    lazy_cache     = mm_empty_lazy_cache(),
    beta           = beta,
    theta          = as.numeric(unlist(fit_result$theta, use.names = FALSE)),
    sigma          = as.numeric(fit_result$sigma),
    logLik         = as.numeric(fit_result$log_likelihood),
    deviance       = as.numeric(fit_result$deviance),
    AIC            = as.numeric(fit_result$aic),
    BIC            = as.numeric(fit_result$bic),
    nobs           = as.integer(fit_result$nobs),
    dof            = as.integer(fit_result$dof),
    df_residual    = as.integer(fit_result$df_residual),
    fit_status     = as.character(fit_result$fit_status),
    std_errors     = std_errors,
    fixed_effect_vcov = fixed_effect_vcov,
    fixed_fitted   = as.numeric(unlist(fit_result$fixed_fitted, use.names = FALSE)),
    fitted         = as.numeric(unlist(fit_result$fitted, use.names = FALSE)),
    residuals      = as.numeric(unlist(fit_result$residuals, use.names = FALSE)),
    random_effects = mm_ranef_from_terms(fit_result$ranef),
    varcorr        = mm_varcorr_from_result(fit_summary$varcorr %||% fit_result$varcorr)
  )
  class(fit) <- c("mm_lmm", "mm_fit", "mm_compiled")
  fit
}

# Pre-fit structural guards (audit-first; PRD §8.1). The Rust fit driver can
# panic ("Matrix index out of bounds") on empty data, or silently return a
# degenerate fit when a random effect is not identifiable (e.g. one
# observation per group). We refuse with a typed condition *before* crossing
# the bridge rather than letting a panic or a misleading fit escape.
mm_validate_fit_structure <- function(spec, lmm = TRUE) {
  n <- nrow(spec$model_frame)
  if (n == 0L) {
    mm_abort(
      message = "`data` has no rows (zero observations after handling missing values); cannot fit a model.",
      class = "mm_data_error"
    )
  }
  if (n == 1L) {
    mm_abort(
      message = "`data` has a single observation; a mixed model needs at least two rows to estimate any variance.",
      class = "mm_not_identifiable",
      nobs = n
    )
  }
  # The "levels < observations" rule guards against the random effect being
  # confounded with the residual scale. It applies to LMMs only: GLMMs have
  # no separate residual variance, so an observation-level random effect
  # (one level per row) is a valid, common overdispersion device that glmer
  # also accepts. Skip the per-group check for GLMMs.
  if (!isTRUE(lmm)) {
    return(invisible(spec))
  }
  for (g in mm_spec_grouping_columns(spec)) {
    cols <- g$columns
    if (!all(cols %in% names(spec$model_frame))) next
    n_levels <- nrow(unique(spec$model_frame[cols]))
    if (n_levels >= n) {
      mm_abort(
        message = sprintf(
          paste0("Grouping factor `%s` has %d level(s) for %d observation(s); ",
                 "with (at least) one level per observation the random effect ",
                 "and the residual variance are not separately identifiable. ",
                 "(lme4 errors here: number of levels of each grouping factor ",
                 "must be < number of observations.)"),
          g$label, n_levels, n
        ),
        class = "mm_not_identifiable",
        group = g$label,
        n_levels = n_levels,
        nobs = n
      )
    }
  }
  invisible(spec)
}

# Extract the grouping-factor column(s) for each random term from the compiled
# semantic model. The engine encodes the group as one of `single$name`,
# `cell$names`, or `interaction$names`.
mm_spec_grouping_columns <- function(spec) {
  terms <- spec$artifact$semantic_model$random_terms %||% list()
  out <- list()
  for (t in terms) {
    g <- t$group
    cols <- NULL
    if (is.character(g) && length(g) == 1L) {
      cols <- g
    } else if (is.list(g)) {
      if (!is.null(g$single$name)) {
        cols <- g$single$name
      } else if (!is.null(g$cell$names)) {
        cols <- unlist(g$cell$names, use.names = FALSE)
      } else if (!is.null(g$interaction$names)) {
        cols <- unlist(g$interaction$names, use.names = FALSE)
      }
    }
    if (!is.null(cols)) {
      cols <- as.character(cols)
      out[[length(out) + 1L]] <- list(columns = cols,
                                      label = paste(cols, collapse = ":"))
    }
  }
  out
}

mm_lmm_weights <- function(expr, data, enclos) {
  if (identical(expr, quote(NULL))) {
    return(NULL)
  }
  weights <- eval(expr, data, enclos)
  if (is.null(weights)) {
    return(NULL)
  }
  if (!is.numeric(weights) || length(weights) != nrow(data) ||
      anyNA(weights) || any(!is.finite(weights)) || any(weights <= 0)) {
    mm_abort(
      message = "`weights` must be a finite positive numeric vector with one value per row in `data`.",
      class = "mm_data_error",
      input = weights
    )
  }
  as.numeric(weights)
}

mm_bridge_weights <- function(weights) {
  if (is.null(weights)) numeric() else as.numeric(weights)
}

mm_json_parse_lmm_fit <- function(json) {
  if (!is.character(json) || length(json) != 1L || is.na(json) || !nzchar(json)) {
    mm_abort(
      message = "`json` must be a single non-empty character string.",
      class = "mm_schema_error",
      input = json
    )
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse LMM fit JSON: %s", conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  schema <- parsed$schema
  if (!is.list(schema) ||
      !identical(as.character(schema$schema_name), "mixeff.lmm_fit_result") ||
      !identical(as.character(schema$schema_version), "1")) {
    mm_abort(
      message = "LMM fit JSON has an unknown schema header.",
      class = "mm_schema_error",
      input = parsed
    )
  }
  parsed
}

mm_json_parse_fit_summary <- function(fit_summary) {
  if (!is.list(fit_summary)) {
    mm_abort(
      message = "LMM fit JSON is missing the fit-summary payload.",
      class = "mm_schema_error",
      input = fit_summary
    )
  }
  if (!identical(as.character(fit_summary$schema_name), "mixedmodels.fit_summary") ||
      !identical(as.character(fit_summary$schema_version), "1.0.0")) {
    mm_abort(
      message = "LMM fit JSON has an unknown fit-summary schema.",
      class = "mm_schema_error",
      input = fit_summary
    )
  }
  if (!is.list(fit_summary$varcorr)) {
    mm_abort(
      message = "LMM fit-summary payload is missing its VarCorr table.",
      class = "mm_schema_error",
      input = fit_summary
    )
  }
  fit_summary
}

mm_named_numeric <- function(values, names) {
  out <- as.numeric(unlist(values, use.names = FALSE))
  nm <- as.character(unlist(names, use.names = FALSE))
  if (length(out) == length(nm)) {
    names(out) <- nm
  }
  out
}

mm_abort_from_bridge <- function(cnd, ...) {
  parts <- mm_split_tagged_error(conditionMessage(cnd))
  cls <- if (!is.na(parts$tag)) parts$tag else "mm_bridge_error"
  msg <- if (!is.na(parts$tag)) parts$message else conditionMessage(cnd)
  mm_abort(
    message = msg,
    class = cls,
    ...,
    parent = cnd
  )
}
