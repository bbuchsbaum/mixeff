# Native fitted-model handles and durable fitted-state restoration.
# Live handles are process-local caches. After serialization, native queries
# restore the engine-owned fitted snapshot without rerunning optimization.
# A different estimand (for example REML -> ML) still needs an explicit refit.
# Optimization controls are not part of the identity of a recorded solution.

# Whether lmm()/glmm() keep the fitted engine model. Documented in NEWS:
# the handle holds a copy of the model data for as long as the fit object
# lives.
mm_keep_handle <- function() {
  isTRUE(getOption("mixeff.keep_handle", TRUE))
}

mm_handle_key_lmm <- function(formula_string, reml, control_json) {
  list(kind = "lmm",
       formula = as.character(formula_string),
       reml = isTRUE(reml),
       control = as.character(control_json))
}

mm_handle_key_glmm <- function(formula_string, family, link, method, n_agq,
                               control_json) {
  list(kind = "glmm",
       formula = as.character(formula_string),
       family = as.character(family),
       link = as.character(link),
       method = as.character(method),
       n_agq = as.integer(n_agq),
       control = as.character(control_json))
}

# Label a freshly returned handle with the key it was fitted under. Returns
# NULL when the bridge returned no handle.
mm_keyed_handle <- function(handle, key) {
  if (is.null(handle) || !identical(typeof(handle), "externalptr")) {
    return(NULL)
  }
  attr(handle, "mm_handle_key") <- key
  handle
}

# Reuse a matching cache or restore the recorded fitted solution. Controls
# cannot change the baseline of a query; only an explicit refit does that.
mm_live_handle <- function(fit, key) {
  snapshot <- mm_fitted_state(fit)
  handle <- fit$rust_handle
  cached_key <- attr(handle, "mm_handle_key", exact = TRUE)
  if (!is.null(cached_key)) {
    cached_key$control <- key$control <- NULL
    # A substituted GLMM is queried under its effective estimator.
    if (identical(key$kind, "glmm")) {
      cached_key$method <- key$method <- NULL
    }
  }
  if (identical(cached_key, key) &&
      identical(attr(handle, "mm_fitted_state", exact = TRUE), snapshot) &&
      isTRUE(mm_handle_alive(handle))) {
    return(handle)
  }
  restored <- tryCatch(
    mm_restore_fit_handle(snapshot, key$kind),
    error = function(cnd) cnd
  )
  if (inherits(restored, "condition")) mm_abort_from_bridge(restored)
  attr(restored, "mm_fitted_state") <- snapshot
  mm_keyed_handle(restored, key)
}

# Live handle of an LMM for an engine refit under `reml` with the fit's own
# formula and control.
mm_lmm_live_handle <- function(fit, reml = isTRUE(fit$REML),
                               formula_string = NULL, control_json = NULL) {
  if (!inherits(fit, "mm_lmm") || !identical(isTRUE(reml), isTRUE(fit$REML))) return(NULL)
  formula_string <- formula_string %||%
    mm_coerce_formula_string(mm_engine_formula(fit))
  control_json <- control_json %||% mm_refit_control_json(fit)
  mm_live_handle(fit, mm_handle_key_lmm(formula_string, reml, control_json))
}

# Live handle of a GLMM for an engine refit with these estimator settings.
mm_glmm_live_handle <- function(fit, family, link, method, n_agq,
                                formula_string = NULL, control_json = NULL) {
  if (!inherits(fit, "mm_glmm")) return(NULL)
  if (!method %in% c(mm_glmm_requested_method(fit), mm_glmm_effective_method(fit)) ||
      !identical(as.integer(n_agq), as.integer(fit$nAGQ))) return(NULL)
  formula_string <- formula_string %||%
    mm_coerce_formula_string(mm_engine_formula(fit))
  control_json <- control_json %||% mm_refit_control_json(fit)
  mm_live_handle(fit, mm_handle_key_glmm(formula_string, family, link, method,
                                         n_agq, control_json))
}

# Placeholder data columns for a bridge call served from a live handle: the
# engine never reads them, so the model frame is not translated.
mm_empty_spec_data <- function() {
  list(
    column_order = character(),
    numeric_columns = list(),
    categorical_values = list(),
    categorical_levels = list(),
    categorical_ordered = character()
  )
}

# Training data for a bridge call: translated only when there is no live
# handle (or the call needs the data anyway, e.g. to fit a reduced model).
mm_bridge_spec_data <- function(fit, handle, need_data = FALSE) {
  if (is.null(handle) || isTRUE(need_data)) {
    mm_translate_data(mm_engine_frame(fit))
  } else {
    mm_empty_spec_data()
  }
}
