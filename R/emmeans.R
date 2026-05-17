#' Optional emmeans support for mixeff LMMs
#'
#' These methods let `emmeans` build reference grids for `mm_lmm` objects when
#' the optional `emmeans` package is installed. They expose the same fixed-effect
#' design surface used by [mm_grid()] and [mm_means()].
#'
#' The current bridge is intentionally narrow: Gaussian LMMs only and
#' population fixed-effect means only. When the fitted artifact carries an
#' available `mixedmodels.fixed_effect_covariance_matrix` payload, `emmeans`
#' receives that full fixed-effect covariance matrix. Native [mm_predictions()],
#' [mm_means()], and [mm_comparisons()] remain the contract-preserving mixeff
#' surface because they preserve row-level status and reason fields.
#'
#' @param object A fitted `mm_lmm`.
#' @param trms,xlev,grid,... Arguments supplied by `emmeans`.
#' @param data Optional data override supplied by `emmeans`.
#'
#' @return Objects expected by `emmeans::recover_data()` and
#'   `emmeans::emm_basis()`.
#'
#' @keywords internal
#' @name emmeans-support
NULL

#' @rdname emmeans-support
#' @export
recover_data.mm_lmm <- function(object, data = NULL, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    mm_abort(
      message = "`recover_data.mm_lmm()` requires the optional emmeans package.",
      class = "mm_inference_unavailable"
    )
  }
  trms <- stats::delete.response(stats::terms(mm_fixed_formula(object)))
  frame <- data %||% object$model_frame
  emmeans::recover_data(
    object$call,
    trms,
    attr(object$model_frame, "na.action"),
    frame = frame,
    ...
  )
}

#' @rdname emmeans-support
#' @export
emm_basis.mm_lmm <- function(object, trms, xlev, grid,
                             method = c("auto", "satterthwaite",
                                        "kenward_roger", "asymptotic",
                                        "none"),
                             ...) {
  if (!requireNamespace("emmeans", quietly = TRUE) ||
      !requireNamespace("estimability", quietly = TRUE)) {
    mm_abort(
      message = "`emm_basis.mm_lmm()` requires the optional emmeans and estimability packages.",
      class = "mm_inference_unavailable"
    )
  }
  method <- match.arg(method)
  m <- stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X_train <- stats::model.matrix(object, type = "fixed")
  X <- stats::model.matrix(trms, m, contrasts.arg = attr(X_train, "contrasts"))
  if (ncol(X) != length(fixef(object))) {
    mm_abort(
      message = "The emmeans reference grid design does not match the fitted fixed effects.",
      class = "mm_inference_unavailable",
      expected = names(fixef(object)),
      observed = colnames(X)
    )
  }
  colnames(X) <- names(fixef(object))
  X <- X[, names(fixef(object)), drop = FALSE]

  bhat <- as.numeric(fixef(object))
  names(bhat) <- names(fixef(object))
  V <- mm_emmeans_vcov(object)
  dfargs <- list(object = object, method = method)
  dffun <- function(k, dfargs) {
    df <- getExportedValue("mixeff", "df_for_contrast")(
      dfargs$object, k, method = dfargs$method
    )
    value <- as.numeric(df[[1L]])
    if (is.finite(value)) value else Inf
  }
  attr(dffun, "mesg") <- sprintf("mixeff %s", method)
  misc <- list(initMesg = mm_emmeans_init_messages(V))

  list(
    X = X,
    bhat = bhat,
    nbasis = estimability::all.estble,
    V = V,
    dffun = dffun,
    dfargs = dfargs,
    misc = misc
  )
}

mm_emmeans_vcov <- function(object) {
  stats::vcov(object, type = "fixed")
}

mm_emmeans_init_messages <- function(V) {
  status <- attr(V, "mm_status") %||% "unknown"
  method <- attr(V, "mm_method") %||% "unknown"
  if (identical(status, "available")) {
    return(sprintf(
      "mixeff emmeans bridge: fixed-effect covariance from mixedmodels.fixed_effect_covariance_matrix (%s); prefer mm_means()/mm_comparisons() when row-level status and reasons are needed.",
      method
    ))
  }
  reason <- attr(V, "mm_unavailable_reason") %||%
    attr(V, "mm_reason") %||%
    "fixed_effect_covariance_matrix_unavailable"
  sprintf(
    "mixeff emmeans bridge: fixed-effect covariance unavailable (%s); standard errors may be unavailable. Prefer mm_means()/mm_comparisons() for contract-preserving status and reasons.",
    reason
  )
}
