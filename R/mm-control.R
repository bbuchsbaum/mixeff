#' Control mixeff fitting behavior
#'
#' `mm_control()` collects small R-side controls for [lmm()]. Phase 1.E uses
#' `verbose` only: `verbose = -1` suppresses the pre-fit [explain_model()]
#' printout; non-negative values print it once before optimization.
#'
#' @param verbose Integer verbosity level. Use `-1` to suppress the automatic
#'   model explanation.
#'
#' @return A list of class `mm_control`.
#'
#' @export
mm_control <- function(verbose = 0L) {
  if (!is.numeric(verbose) || length(verbose) != 1L || is.na(verbose)) {
    mm_abort(
      message = "`verbose` must be a single numeric value.",
      class = "mm_fit_error",
      input = verbose
    )
  }
  out <- list(verbose = as.integer(verbose))
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
      class = "mm_fit_error",
      input = control
    )
  }
  verbose <- control$verbose %||% 0L
  mm_control(verbose = verbose)
}
