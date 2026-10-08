# Internal refits of a stored fit (REML -> ML for comparison, reduced models
# for drop1 / random-term LRTs, refit()).
#
# lmm()/glmm() evaluate `weights`/`offset` lme4-style: the argument
# expression is looked up in `data` first. An internal call written as
# `lmm(..., data = fit$model_frame, weights = fit$weights)` therefore resolves
# `fit` to a DATA COLUMN when the user's data has one named `fit` (or `full`,
# `object`, ...). These helpers pass the stored vectors as literal values in
# the constructed call, so nothing is looked up in the data mask, and then
# replace the literal in the recorded call with a symbol so `fit$call` stays
# printable.

mm_internal_lmm <- function(formula, data, REML, weights = NULL,
                            offset = NULL,
                            control = mm_control(verbose = -1)) {
  cl <- as.call(list(
    quote(lmm),
    formula = formula,
    data = quote(data),
    REML = isTRUE(REML),
    weights = weights,
    offset = offset,
    control = quote(control)
  ))
  fit <- eval(cl)
  if (!is.null(weights)) fit$call$weights <- quote(weights)
  if (!is.null(offset)) fit$call$offset <- quote(offset)
  fit
}

mm_internal_glmm <- function(formula, data, family, weights = NULL,
                             offset = NULL, method, nAGQ,
                             control = mm_control(verbose = -1)) {
  cl <- as.call(list(
    quote(glmm),
    formula = formula,
    data = quote(data),
    family = quote(family),
    weights = weights,
    offset = offset,
    method = method,
    nAGQ = nAGQ,
    control = quote(control)
  ))
  fit <- eval(cl)
  if (!is.null(weights)) fit$call$weights <- quote(weights)
  if (!is.null(offset)) fit$call$offset <- quote(offset)
  fit
}

# The user's mm_control() for an internal refit, with the pre-fit explanation
# and notices silenced. `keep_start = FALSE` drops a warm start, for refits
# of a DIFFERENT model (reduced or null models) whose theta dimension can
# differ from the stored fit's.
mm_internal_control <- function(fit, keep_start = TRUE) {
  control <- fit$control %||% mm_control()
  control <- unclass(control)
  control$verbose <- -1L
  if (!isTRUE(keep_start)) control$start <- NULL
  mm_validate_control(control)
}
