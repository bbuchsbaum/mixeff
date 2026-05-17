#' Show requested, effective, and fitted model-state changes
#'
#' `changes()` summarizes the transitions recorded in the compiler artifact:
#' requested formula to effective formula, design-time reductions or covariance
#' transitions, and fitted covariance rank/status from the optimizer
#' certificate pass.
#'
#' @param fit A compiled `mm_spec` or fitted `mm_fit`.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_change_log` object with a data-frame `table` and the raw
#'   artifact fragments used to build it.
#'
#' @export
changes <- function(fit, ...) {
  UseMethod("changes")
}

#' @rdname changes
#' @export
changes.mm_compiled <- function(fit, ...) {
  artifact <- mm_compiled_artifact(fit)
  rows <- c(
    mm_change_formula_rows(artifact),
    mm_change_reduction_rows(artifact),
    mm_change_covariance_transition_rows(artifact),
    mm_change_effective_covariance_rows(artifact)
  )
  table <- if (length(rows)) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  } else {
    mm_change_empty_table()
  }

  obj <- list(
    table = table,
    reductions = artifact$reductions %||% list(),
    covariance_transitions = artifact$covariance_transitions %||% list(),
    effective_covariance = artifact$effective_covariance %||% list(),
    fit_status = artifact$optimizer_certificate$status %||% "not_assessed"
  )
  class(obj) <- "mm_change_log"
  obj
}

#' @method print mm_change_log
#' @export
print.mm_change_log <- function(x, ...) {
  cat("Model changes:\n")
  if (!nrow(x$table)) {
    cat("  none recorded\n")
    return(invisible(x))
  }
  print(x$table, row.names = FALSE)
  invisible(x)
}

mm_change_formula_rows <- function(artifact) {
  requested <- artifact$requested_formula %||% ""
  effective <- artifact$effective_formula %||% requested
  list(data.frame(
    stage = "semantic_ir",
    term_id = "",
    group = "",
    requested = requested,
    effective = effective,
    fitted = artifact$optimizer_certificate$status %||% "not_assessed",
    status = if (identical(requested, effective)) "unchanged" else "canonicalized",
    detail = "formula display",
    stringsAsFactors = FALSE
  ))
}

mm_change_reduction_rows <- function(artifact) {
  lapply(artifact$reductions %||% list(), function(reduction) {
    data.frame(
      stage = "design_time_reduction",
      term_id = mm_scalar_text(reduction$term_id),
      group = mm_scalar_text(reduction$group),
      requested = mm_list_text(reduction$from %||% reduction$requested),
      effective = mm_list_text(reduction$to %||% reduction$effective),
      fitted = artifact$optimizer_certificate$status %||% "not_assessed",
      status = mm_scalar_text(reduction$trigger %||% reduction$status, "reduction"),
      detail = mm_scalar_text(reduction$reason),
      stringsAsFactors = FALSE
    )
  })
}

mm_change_covariance_transition_rows <- function(artifact) {
  lapply(artifact$covariance_transitions %||% list(), function(transition) {
    data.frame(
      stage = "covariance_transition",
      term_id = mm_scalar_text(transition$term_id),
      group = mm_scalar_text(transition$group),
      requested = mm_scalar_text(transition$from %||% transition$requested_family),
      effective = mm_scalar_text(transition$to %||% transition$effective_family),
      fitted = artifact$optimizer_certificate$status %||% "not_assessed",
      status = mm_scalar_text(transition$reason %||% transition$status, "transition"),
      detail = mm_scalar_text(transition$detail),
      stringsAsFactors = FALSE
    )
  })
}

mm_change_effective_covariance_rows <- function(artifact) {
  lapply(artifact$effective_covariance %||% list(), function(summary) {
    requested_rank <- summary$requested_rank %||% NA_integer_
    supported_rank <- summary$supported_rank %||% NA_integer_
    detail <- sprintf(
      "requested rank %s; fitted rank %s",
      requested_rank,
      supported_rank
    )
    data.frame(
      stage = "certificate_time",
      term_id = mm_scalar_text(summary$term_id),
      group = "",
      requested = mm_list_text(summary$requested_basis),
      effective = mm_list_text(summary$requested_basis),
      fitted = mm_scalar_text(summary$status, "not_assessed"),
      status = mm_scalar_text(summary$status, "not_assessed"),
      detail = detail,
      stringsAsFactors = FALSE
    )
  })
}

mm_change_empty_table <- function() {
  data.frame(
    stage = character(),
    term_id = character(),
    group = character(),
    requested = character(),
    effective = character(),
    fitted = character(),
    status = character(),
    detail = character(),
    stringsAsFactors = FALSE
  )
}
