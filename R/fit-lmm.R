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
#' @param subset Optional expression selecting rows of `data`, evaluated in
#'   `data` (as in [stats::lm()]).
#' @param na.action Optional function controlling missing-value handling,
#'   applied to the model variables before fitting (e.g. [stats::na.omit]).
#'   The default (`NULL`) refuses any `NA` in a model variable with a typed
#'   `mm_data_error` (audit-first: missing-data dropping must be opt-in). Pass
#'   `na.action = na.omit` for lme4's complete-case behaviour.
#' @param contrasts Optional named list of factor contrasts. The engine codes
#'   all factors with treatment contrasts; a request for any other coding is
#'   refused (recode the factor instead).
#' @param control A list from [mm_control()].
#'
#' @return An object of class `mm_lmm`, also inheriting from `mm_fit` and
#'   `mm_compiled`.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rnorm(80),
#'   x = rnorm(80),
#'   subject = factor(rep(seq_len(20), each = 4))
#' )
#' fit <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))
#' fixef(fit)
#' VarCorr(fit)
#' summary(fit)
#'
#' @export
lmm <- function(formula, data, REML = TRUE, weights = NULL,
                subset = NULL, na.action = NULL, contrasts = NULL,
                control = mm_control()) {
  call <- match.call()
  control <- mm_validate_control(control)
  if (!is.null(contrasts)) mm_reject_nontreatment_contrasts(contrasts)
  weights <- mm_lmm_weights(substitute(weights), data, parent.frame())
  if (!is.logical(REML) || length(REML) != 1L || is.na(REML)) {
    mm_abort(
      message = "`REML` must be TRUE or FALSE.",
      class = "mm_fit_error",
      input = REML
    )
  }

  # lme4-style data preparation: row subset, then NA handling. Both run before
  # the design is compiled so the engine sees a clean, complete frame.
  prep <- mm_prepare_fit_data(formula, data, substitute(subset), na.action,
                              weights, parent.frame())
  data <- prep$data
  weights <- prep$weights

  # lme4 parity: grouping variables must be categorical. Coerce non-factor /
  # non-character grouping columns to factors (announced, not silent) so an
  # integer subject/item ID does not hit the native "grouping factor not
  # categorical" refusal.
  data <- mm_apply_grouping_coercion(formula, data, control$verbose)

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
    varcorr        = mm_varcorr_from_result(
      fit_summary$varcorr %||% fit_result$varcorr,
      artifact = artifact
    )
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

# Recursively collect the random-effect bar terms (`(... | g)`, including the
# `||` uncorrelated form) from a formula's right-hand side. A pared-down
# equivalent of lme4::findbars() that takes no lme4 dependency. Returns a list
# of `call` objects, each a `|`/`||` expression whose third element is the
# grouping expression.
mm_find_bars <- function(term) {
  if (!is.call(term)) {
    return(list())
  }
  op <- term[[1L]]
  if (identical(op, quote(`|`)) || identical(op, quote(`||`))) {
    return(list(term))
  }
  if (identical(op, quote(`(`))) {
    return(mm_find_bars(term[[2L]]))
  }
  # Recurse into the operands of +, *, :, etc.
  do.call(c, lapply(as.list(term)[-1L], mm_find_bars))
}

# Names of the variables that appear as grouping factors (right of a `|`) in an
# lme4-style formula. `(1 + x | a:b)` and `(1 | a/b)` both contribute `a` and
# `b`. Used to decide which columns must be categorical before the fit.
mm_formula_grouping_vars <- function(formula) {
  if (!inherits(formula, "formula")) {
    return(character(0))
  }
  rhs <- formula[[length(formula)]]
  bars <- mm_find_bars(rhs)
  if (!length(bars)) {
    return(character(0))
  }
  unique(unlist(lapply(bars, function(b) all.vars(b[[3L]])), use.names = FALSE))
}

# Coerce any grouping variable that is not already categorical (factor or
# character) to a factor, matching lme4 / nlme / glmmTMB, which all silently
# factor() their grouping variables. The engine's data bridge sends
# numeric/integer/logical columns as numeric, so an integer ID column would
# otherwise be rejected by the native fit constructor ("grouping factor not
# categorical"). This is *not* silent surgery: the caller surfaces a typed,
# suppressible notice naming every coerced column (see
# `mm_grouping_coercion_notice()`).
#
# Character grouping columns are left untouched — they are already categorical
# in the bridge, and re-leveling them would only churn random-effect labels.
#
# @return list(data, coerced, classes): the (possibly modified) data frame, the
#   names of coerced columns, and a named character vector of their original
#   storage classes (parallel to `coerced`).
mm_coerce_grouping_factors <- function(formula, data) {
  gvars <- intersect(mm_formula_grouping_vars(formula), names(data))
  coerced <- character(0)
  classes <- character(0)
  for (v in gvars) {
    col <- data[[v]]
    if (is.factor(col) || is.character(col)) {
      next
    }
    classes[[v]] <- paste(class(col), collapse = "/")
    data[[v]] <- factor(col)
    coerced <- c(coerced, v)
  }
  list(data = data, coerced = coerced, classes = classes)
}

# Human-readable, suppressible notice describing a grouping-factor coercion.
mm_grouping_coercion_notice <- function(coerced, classes) {
  items <- vapply(
    coerced,
    function(v) sprintf("`%s` (%s)", v, classes[[v]]),
    character(1)
  )
  sprintf(
    paste0(
      "Coerced grouping variable%s %s to a factor for the random-effects ",
      "structure (lme4 does the same). The fit is unaffected; wrap the ",
      "variable in `factor()` yourself to control the level order. Silence ",
      "this notice with mm_control(verbose = -1)."
    ),
    if (length(coerced) > 1L) "s" else "",
    paste(items, collapse = ", ")
  )
}

# Coerce grouping variables to factors and (unless silenced) emit the notice.
# Shared by lmm() and glmm(); call before compile_model() so the whole audit /
# explain / fit pipeline sees consistent categorical grouping.
mm_apply_grouping_coercion <- function(formula, data, verbose) {
  res <- mm_coerce_grouping_factors(formula, data)
  if (length(res$coerced) && isTRUE(verbose >= 0L)) {
    mm_inform(
      mm_grouping_coercion_notice(res$coerced, res$classes),
      class = "mm_grouping_coercion_notice"
    )
  }
  res$data
}

# lme4-style data preparation shared by the fit drivers: apply a `subset`
# expression and an `na.action` to `data` (and the already-evaluated `weights`
# vector, kept aligned). Returns the processed `data` and `weights`. The NA
# policy is applied by *calling* the supplied na.action on the model variables,
# so na.omit / na.exclude / na.fail / na.pass all behave as their authors
# intend; the default (NULL) leaves NAs for compile_model() to refuse.
mm_prepare_fit_data <- function(formula, data, subset_expr, na.action,
                                weights, enclos) {
  if (!is.data.frame(data)) {
    mm_abort(message = "`data` must be a data.frame.", class = "mm_data_error",
             input = data)
  }
  n0 <- nrow(data)

  # ---- subset ----
  if (!identical(subset_expr, quote(NULL)) && !is.null(subset_expr)) {
    keep <- eval(subset_expr, data, enclos)
    if (is.logical(keep)) {
      if (length(keep) != n0) {
        mm_abort(
          message = sprintf("`subset` must have one value per row of `data` (%d).", n0),
          class = "mm_arg_error", input = keep
        )
      }
      keep[is.na(keep)] <- FALSE
      idx <- which(keep)
    } else if (is.numeric(keep)) {
      idx <- as.integer(keep)
      idx <- idx[!is.na(idx)]
    } else {
      mm_abort(message = "`subset` must evaluate to a logical or integer index.",
               class = "mm_arg_error", input = keep)
    }
    data <- data[idx, , drop = FALSE]
    if (!is.null(weights)) weights <- weights[idx]
  }

  # ---- na.action ----
  if (!is.null(na.action)) {
    if (!is.function(na.action)) na.action <- match.fun(na.action)
    vars <- intersect(all.vars(formula), names(data))
    sub <- data[, vars, drop = FALSE]
    cleaned <- na.action(sub)                 # may error (na.fail) or drop rows
    omitted <- attr(cleaned, "na.action")
    if (!is.null(omitted)) {
      drop_idx <- as.integer(omitted)
      data <- data[-drop_idx, , drop = FALSE]
      if (!is.null(weights)) weights <- weights[-drop_idx]
    }
  }

  list(data = data, weights = weights)
}

# The Rust engine codes every factor with treatment (dummy) contrasts. Accept a
# `contrasts` request only if it asks for that coding; otherwise refuse rather
# than silently fit a differently-coded design than the user asked for.
mm_reject_nontreatment_contrasts <- function(contrasts) {
  is_treatment <- function(x) {
    is.character(x) && length(x) == 1L &&
      x %in% c("contr.treatment", "contr.SAS")
  }
  ok <- is.list(contrasts) && length(contrasts) &&
    all(vapply(contrasts, is_treatment, logical(1)))
  if (!ok) {
    mm_abort(
      message = paste(
        "Custom `contrasts` are not supported: the engine codes all factors",
        "with treatment contrasts. Recode the factor (e.g. relevel(), or",
        "construct the desired numeric columns) to obtain a different coding."
      ),
      class = "mm_arg_error",
      input = contrasts
    )
  }
  invisible(TRUE)
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
