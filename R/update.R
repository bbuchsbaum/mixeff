# update() methods for mixeff fits.
#
# Mirrors stats::update() / lme4's update() for the common model-building
# idioms: editing the formula (`. ~ . - x`), toggling REML, swapping
# weights/offset/family/control, or supplying new data. The re-fit reuses the
# stored model frame as the default data source so formula edits that *remove*
# terms or change estimation options work without re-supplying data. Formula
# edits that introduce a variable not present in the original model frame
# require an explicit `data =` argument (the model frame only carries variables
# the original model used).

#' Update and re-fit a mixeff model
#'
#' `update()` re-fits an [`mm_lmm`][lmm] or [`mm_glmm`][glmm] with a modified
#' model specification, mirroring [stats::update()] and `lme4`'s `update()` for
#' the common cases: changing the formula (`. ~ . - x`), toggling `REML`,
#' swapping `weights`/`offset`/`family`/`control`, or supplying new `data`.
#'
#' The re-fit reuses the fitted model frame ([model.frame()]) as the default
#' data source, so formula edits that *remove* terms or change estimation
#' options work without re-supplying data. A formula edit that introduces a
#' **new** variable absent from the original model frame requires an explicit
#' `data =` argument.
#'
#' @param object A fitted `mm_lmm` or `mm_glmm`.
#' @param formula. A formula-change applied with [stats::update.formula()];
#'   omit to keep the current formula. Random-effect terms (`(x | g)`,
#'   `(x || g)`) are preserved across `. ~ .` edits.
#' @param ... Arguments to override on the re-fit. For `mm_lmm`: `data`,
#'   `REML`, `weights`, `control`. For `mm_glmm`: additionally `family`,
#'   `offset`, `method`, `nAGQ`, `inference`.
#' @param evaluate If `TRUE` (default) re-fit and return the new model; if
#'   `FALSE` return the unevaluated call.
#'
#' @return A new fitted model of the same class as `object`, or an unevaluated
#'   call when `evaluate = FALSE`.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rnorm(80), x = rnorm(80), z = rnorm(80),
#'   g = factor(rep(seq_len(10), each = 8))
#' )
#' fit <- lmm(y ~ x + z + (1 | g), df, control = mm_control(verbose = -1))
#' # drop a fixed term
#' fit2 <- update(fit, . ~ . - z)
#' # refit by ML for a likelihood-ratio comparison
#' fit_ml <- update(fit, REML = FALSE)
#' fixef(fit2)
#'
#' @name update.mm
#' @importFrom stats update update.formula formula model.frame
#' @method update mm_lmm
#' @export
update.mm_lmm <- function(object, formula., ..., evaluate = TRUE) {
  overrides <- list(...)
  new_formula <- if (missing(formula.)) {
    stats::formula(object)
  } else {
    stats::update.formula(stats::formula(object), formula.)
  }
  args <- list(
    formula = new_formula,
    data    = mm_update_arg(overrides, "data", stats::model.frame(object)),
    REML    = mm_update_arg(overrides, "REML", isTRUE(object$REML)),
    weights = mm_update_arg(overrides, "weights", object$weights),
    control = mm_update_arg(overrides, "control",
                            object$control %||% mm_control(verbose = -1))
  )
  mm_refit_call("lmm", args, evaluate, parent.frame())
}

#' @rdname update.mm
#' @method update mm_glmm
#' @export
update.mm_glmm <- function(object, formula., ..., evaluate = TRUE) {
  overrides <- list(...)
  new_formula <- if (missing(formula.)) {
    stats::formula(object)
  } else {
    stats::update.formula(stats::formula(object), formula.)
  }
  args <- list(
    formula   = new_formula,
    data      = mm_update_arg(overrides, "data", stats::model.frame(object)),
    family    = mm_update_arg(overrides, "family",
                              mm_glmm_family_from_info(object$family)),
    weights   = mm_update_arg(overrides, "weights", object$weights),
    offset    = mm_update_arg(overrides, "offset", object$offset),
    method    = mm_update_arg(overrides, "method", object$method),
    nAGQ      = mm_update_arg(overrides, "nAGQ", object$nAGQ),
    inference = mm_update_arg(overrides, "inference", object$inference_request),
    control   = mm_update_arg(overrides, "control",
                              object$control %||% mm_control(verbose = -1))
  )
  mm_refit_call("glmm", args, evaluate, parent.frame())
}

# Pick an override value if the user supplied the name (even as NULL), else the
# stored default. Distinguishes "weights = NULL" (drop weights) from "weights
# not mentioned" (keep stored weights).
mm_update_arg <- function(overrides, name, default) {
  if (name %in% names(overrides)) overrides[[name]] else default
}

# Reconstruct an R family object from the stored engine family info so a GLMM
# can be re-fit without the caller re-supplying `family`. Inverse of
# mm_glmm_family_info(); family_info$family is the lowercased R-facing name.
mm_glmm_family_from_info <- function(info) {
  fam <- info$family
  link <- info$link
  switch(
    fam,
    binomial = stats::binomial(link = link),
    poisson  = stats::poisson(link = link),
    gamma    = stats::Gamma(link = link),
    mm_abort(
      message = sprintf(
        "Cannot reconstruct a family object for engine family `%s`; pass `family =` explicitly to update().",
        fam
      ),
      class = "mm_arg_error",
      input = info
    )
  )
}

# Build and (optionally) evaluate the re-fit call. Argument values are inlined
# into the call literally, so it is self-contained (the stored model frame is
# carried as data). evaluate = FALSE returns the call for inspection/deferred
# evaluation, matching stats::update(evaluate = FALSE).
mm_refit_call <- function(fitter_name, args, evaluate, envir) {
  cl <- as.call(c(list(as.name(fitter_name)), args))
  if (isTRUE(evaluate)) eval(cl, envir = envir) else cl
}
