# Native fitted-model handles (pre-CRAN checklist §5.2(c)).
#
# `lmm()` / `glmm()` keep the fitted engine model alive behind an external
# pointer stored as `fit$rust_handle`. Follow-on computations (contrasts,
# summary() term tables, predictions and their variance, ranef(condVar),
# profiles, bootstraps, model comparison, convergence verification) run on
# that live model instead of rebuilding and refitting it cold from the stored
# model frame. The handle is a process-local cache, never the source of
# truth: the fit object still carries every durable result, and when the
# handle is absent or dead (after saveRDS()/readRDS(), in another process,
# or with `options(mixeff.keep_handle = FALSE)`) every verb takes the cold
# refit path, which produces identical results.
#
# Each handle carries a key (attribute `mm_handle_key`): the exact engine
# arguments it was fitted with (formula string, REML or GLMM family / link /
# method / nAGQ, and control JSON). A verb uses the handle only when the
# refit it would otherwise run has the same key -- e.g. a profile under ML of
# a REML fit, or a GLMM bootstrap of the fallback estimator of a substituted
# fit, still refits cold. Copies of a fit share the handle; that is safe
# because the engine only reads the cached model (mutating verbs work on a
# clone).

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

# The fit's live handle when it was fitted under `key`, else NULL.
mm_live_handle <- function(fit, key) {
  handle <- fit$rust_handle
  if (is.null(handle) || !identical(typeof(handle), "externalptr")) {
    return(NULL)
  }
  if (!identical(attr(handle, "mm_handle_key", exact = TRUE), key)) {
    return(NULL)
  }
  if (!isTRUE(mm_handle_alive(handle))) {
    return(NULL)
  }
  handle
}

# Live handle of an LMM for an engine refit under `reml` with the fit's own
# formula and control.
mm_lmm_live_handle <- function(fit, reml = isTRUE(fit$REML),
                               formula_string = NULL, control_json = NULL) {
  if (!inherits(fit, "mm_lmm") || is.null(fit$rust_handle)) return(NULL)
  formula_string <- formula_string %||%
    mm_coerce_formula_string(mm_engine_formula(fit))
  control_json <- control_json %||% mm_refit_control_json(fit)
  mm_live_handle(fit, mm_handle_key_lmm(formula_string, reml, control_json))
}

# Live handle of a GLMM for an engine refit with these estimator settings.
mm_glmm_live_handle <- function(fit, family, link, method, n_agq,
                                formula_string = NULL, control_json = NULL) {
  if (!inherits(fit, "mm_glmm") || is.null(fit$rust_handle)) return(NULL)
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
