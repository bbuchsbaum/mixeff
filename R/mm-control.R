#' Control mixeff fitting behavior
#'
#' `mm_control()` collects small R-side controls for [lmm()] and [glmm()].
#' `verbose = -1` suppresses the pre-fit [explain_model()] printout;
#' non-negative values print it once before optimization. `max_feval` caps the
#' number of objective evaluations the optimizer may use.
#'
#' @param verbose Integer verbosity level. Use `-1` to suppress the automatic
#'   model explanation.
#' @param max_feval Optional positive integer capping the optimizer's objective
#'   evaluations. This is most useful for [glmm()] with
#'   `method = "joint_laplace"`, whose native joint optimizer otherwise runs to
#'   an engine-chosen budget that can be expensive on large high-baseline
#'   models. Lower values trade some accuracy for speed. `NULL` (default) leaves
#'   the engine default in place.
#'
#' @return A list of class `mm_control`.
#'
#' @export
mm_control <- function(verbose = 0L, max_feval = NULL) {
  if (!is.numeric(verbose) || length(verbose) != 1L || is.na(verbose)) {
    mm_abort(
      message = "`verbose` must be a single numeric value.",
      class = "mm_arg_error",
      input = verbose
    )
  }
  out <- list(verbose = as.integer(verbose))
  if (!is.null(max_feval)) {
    if (!is.numeric(max_feval) || length(max_feval) != 1L ||
        is.na(max_feval) || max_feval < 1) {
      mm_abort(
        message = "`max_feval` must be `NULL` or a single positive integer.",
        class = "mm_arg_error",
        input = max_feval
      )
    }
    out$max_feval <- as.integer(max_feval)
  }
  class(out) <- "mm_control"
  out
}

mm_validate_control <- function(control) {
  if (missing(control) || is.null(control)) {
    return(mm_control())
  }
  if (!is.list(control)) {
    mm_abort(
      message = "`control` must be a list created by mm_control().",
      class = "mm_arg_error",
      input = control
    )
  }
  verbose <- control$verbose %||% 0L
  mm_control(verbose = verbose, max_feval = control$max_feval)
}
