#' Optional emmeans support for mixeff models
#'
#' These methods let `emmeans` build reference grids for `mm_lmm` and
#' `mm_glmm` objects when the optional `emmeans` package is installed. They
#' expose the same fixed-effect design surface used by [mm_grid()] and
#' [mm_means()].
#'
#' The bridge covers population fixed-effect means only. For GLMMs the whole
#' emmeans surface is gated by the package's inference-capability contract:
#' a default profiled fit refuses with a typed error (refit with
#' `method = "joint_laplace"` for certified Wald inference, or opt in to the
#' labelled working-Hessian approximation at fit time with
#' `inference = "working_hessian"`). When the fitted artifact carries a
#' usable `mixedmodels.fixed_effect_covariance_matrix` payload, `emmeans`
#' receives that full fixed-effect covariance matrix. Native
#' [mm_predictions()], [mm_means()], and [mm_comparisons()] remain the
#' contract-preserving mixeff surface for LMMs because they preserve
#' row-level status and reason fields.
#'
#' @param object A fitted `mm_lmm` or `mm_glmm`.
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
  # Build the basis in the engine's coefficient basis and aligned to
  # names(beta). emmeans supplies `trms`/`xlev`/`grid`, but R's default
  # contrasts (contr.poly for ordered factors) and interaction column order
  # disagree with the engine, so we reconstruct from the fit's own fixed
  # formula via the shared helper. See mm_engine_fixed_matrix() in predict.R.
  X <- mm_engine_fixed_matrix(object, grid)

  bhat <- as.numeric(fixef(object))
  names(bhat) <- names(fixef(object))
  V <- mm_emmeans_vcov(object)
  dfargs <- list(object = object, method = method)
  dffun <- function(k, dfargs) {
    df <- getExportedValue("mixeff", "df_for_contrast")(
      dfargs$object, k, method = dfargs$method
    )
    # df_for_contrast() returns an mm_* object; $df is the named numeric
    # vector (one entry per contrast row of k).
    value <- as.numeric(df$df)
    ifelse(is.finite(value), value, Inf)
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

#' @rdname emmeans-support
#' @export
recover_data.mm_glmm <- function(object, data = NULL, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    mm_abort(
      message = "`recover_data.mm_glmm()` requires the optional emmeans package.",
      class = "mm_inference_unavailable"
    )
  }
  # The emmeans surface is inference (SEs, z tests, CIs on the grid), so it
  # is gated by the same capability arbiter as every other Wald route
  # (WI-2.2). Refusing here, at grid recovery, stops the whole emmeans
  # pipeline with one actionable message instead of letting emm_basis()
  # fabricate inference from an uncertified covariance.
  cap <- mm_glmm_inference_capability(object)
  if (!isTRUE(cap$wald)) {
    mm_abort(
      message = mm_glmm_wald_refusal_message("emmeans grid summaries"),
      class = "mm_inference_unavailable",
      reason_code = cap$reason_code,
      input = object
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
emm_basis.mm_glmm <- function(object, trms, xlev, grid, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE) ||
      !requireNamespace("estimability", quietly = TRUE)) {
    mm_abort(
      message = "`emm_basis.mm_glmm()` requires the optional emmeans and estimability packages.",
      class = "mm_inference_unavailable"
    )
  }
  # Same arbiter gate as recover_data.mm_glmm: emm_basis() can also be
  # reached directly (emmeans passes fitted objects around), so both entry
  # points enforce the contract independently.
  cap <- mm_glmm_inference_capability(object)
  if (!isTRUE(cap$wald)) {
    mm_abort(
      message = mm_glmm_wald_refusal_message("emmeans grid summaries"),
      class = "mm_inference_unavailable",
      reason_code = cap$reason_code,
      input = object
    )
  }
  m <- stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X_train <- stats::model.matrix(object, type = "fixed")
  X <- stats::model.matrix(trms, m, contrasts.arg = attr(X_train, "contrasts"))
  # fixef() carries lme4/model.matrix-style names, so the reference-grid
  # design aligns by NAME; a positional rename would silently misassign
  # columns whenever the grid's column order differs from the fit's.
  missing_cols <- setdiff(names(fixef(object)), colnames(X))
  if (length(missing_cols)) {
    mm_abort(
      message = paste0(
        "The emmeans reference grid design does not match the fitted fixed ",
        "effects; missing column(s): ", paste(missing_cols, collapse = ", "), "."
      ),
      class = "mm_inference_unavailable",
      expected = names(fixef(object)),
      observed = colnames(X)
    )
  }
  X <- X[, names(fixef(object)), drop = FALSE]

  bhat <- as.numeric(fixef(object))
  names(bhat) <- names(fixef(object))
  V <- mm_emmeans_vcov(object)
  ## GLMM Wald is asymptotic z; df is +Inf so emmeans uses the
  ## standard-normal z reference distribution.
  dffun <- function(k, dfargs) Inf
  dfargs <- list()
  misc <- list(initMesg = mm_emmeans_init_messages(V))
  ## Hand emmeans the family link so type = "response" applies linkinv
  ## naturally. The family object on object$family is the standard R
  ## family() shape (linkfun, linkinv, mu.eta, valideta).
  fam <- object$family
  if (!is.null(fam) && is.list(fam) && !is.null(fam$link)) {
    misc <- emmeans::.std.link.labels(fam, misc)
  }

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
      paste0(
        "mixeff emmeans bridge: fixed-effect covariance from ",
        "mixedmodels.fixed_effect_covariance_matrix (%s); prefer ",
        "mm_means()/mm_comparisons() when row-level status and reasons are needed."
      ),
      method
    ))
  }
  if (identical(status, "available_noninferential")) {
    # Reachable only through the explicit fit-level opt-in
    # (`inference = "working_hessian"`): the capability gate refuses
    # profiled fits before emmeans gets this far otherwise.
    return(sprintf(
      paste0(
        "mixeff emmeans bridge: standard errors and tests use the ",
        "UNCERTIFIED working-Hessian approximation (%s; reliability ",
        "moderate), as requested via inference = \"working_hessian\". ",
        "Refit with method = \"joint_laplace\" for certified inference."
      ),
      method
    ))
  }
  reason <- attr(V, "mm_unavailable_reason") %||%
    attr(V, "mm_reason") %||%
    "fixed_effect_covariance_matrix_unavailable"
  sprintf(
    paste0(
      "mixeff emmeans bridge: fixed-effect covariance unavailable (%s); ",
      "standard errors may be unavailable. Prefer mm_means()/mm_comparisons() ",
      "for contract-preserving status and reasons."
    ),
    reason
  )
}
