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
    fit_status = artifact$optimizer_certificate$status %||% "not_assessed",
    requested_formula = artifact$requested_formula %||% "",
    effective_formula = artifact$effective_formula %||%
      artifact$requested_formula %||% ""
  )
  class(obj) <- "mm_change_log"
  obj
}

#' @method print mm_change_log
#' @export
print.mm_change_log <- function(x, ...) {
  cat("Model changes:\n")
  lines <- mm_change_sentence_lines(x)
  status <- as.character(x$fit_status %||% "not_assessed")
  converged <- startsWith(status, "converged")
  assessed <- !identical(status, "not_assessed")
  if (length(lines)) {
    cat(sprintf("  %s\n", lines), sep = "")
    if (assessed && !converged) {
      cat(sprintf(
        "  Note: the optimizer stopped early (fit status `%s`); fitted-rank statements reflect the last accepted iterate.\n",
        status
      ))
    }
  } else if (!assessed) {
    cat("  none recorded at compile time.\n")
  } else if (converged) {
    cat("  none: the model was fitted as requested.\n")
  } else {
    cat(sprintf(
      "  none: no structural change was made; the optimizer stopped early (fit status `%s`).\n",
      status
    ))
  }
  cat("Stage-by-stage records available via $table.\n")
  invisible(x)
}

# One plain-language sentence per recorded change. The certificate-time rank
# summary is the canonical statement of a boundary event; reduction and
# transition entries that restate the same event (trigger
# "certificate_time_boundary" on a term whose rank summary already printed)
# are skipped here -- they remain in $table, so nothing is hidden.
mm_change_sentence_lines <- function(x) {
  labels <- mm_change_term_labels(x$effective_covariance)
  lines <- character()

  requested <- x$requested_formula %||% ""
  effective <- x$effective_formula %||% requested
  if (nzchar(requested) && !identical(requested, effective)) {
    lines <- c(lines, sprintf(
      "Formula canonicalized: wrote `%s`; effective form is `%s`.",
      requested, effective
    ))
  }

  rank_changed <- character()
  for (summary in x$effective_covariance %||% list()) {
    requested_rank <- summary$requested_rank %||% NA_integer_
    supported_rank <- summary$supported_rank %||% NA_integer_
    if (is.na(requested_rank) || is.na(supported_rank) ||
        requested_rank == supported_rank) {
      next
    }
    term <- mm_scalar_text(summary$term_id)
    rank_changed <- c(rank_changed, term)
    lines <- c(lines, sprintf(
      "Fitted covariance for %s: requested rank %s, fitted rank %s [%s].",
      mm_change_term_label(term, labels, summary$source_syntax),
      requested_rank, supported_rank,
      mm_scalar_text(summary$status, "reduced_rank")
    ))
  }

  for (reduction in x$reductions %||% list()) {
    term <- mm_scalar_text(reduction$affected_term %||% reduction$term_id)
    trigger <- mm_scalar_text(reduction$trigger %||% reduction$status,
                              "reduction")
    if (mm_change_is_rank_restatement(trigger, term, rank_changed)) next
    lines <- c(lines, sprintf(
      "%s: %s [%s].",
      mm_change_term_label(term, labels),
      mm_scalar_text(reduction$reason, "reduced"),
      trigger
    ))
  }

  for (transition in x$covariance_transitions %||% list()) {
    term <- mm_scalar_text(transition$affected_term %||% transition$term_id)
    trigger <- mm_scalar_text(transition$trigger %||% transition$reason,
                              "transition")
    if (mm_change_is_rank_restatement(trigger, term, rank_changed)) next
    lines <- c(lines, sprintf(
      "%s: covariance changed from %s to %s [%s].",
      mm_change_term_label(term, labels),
      mm_change_covariance_text(transition$from %||%
                                  transition$requested_family),
      mm_change_covariance_text(transition$to %||%
                                  transition$effective_family),
      trigger
    ))
  }

  lines
}

mm_change_is_rank_restatement <- function(trigger, term, rank_changed) {
  identical(trigger, "certificate_time_boundary") && term %in% rank_changed
}

# Map opaque term ids ("r0") to the user's own syntax ("(1 | s)") using the
# certificate-time summaries, which carry both.
mm_change_term_labels <- function(summaries) {
  labels <- list()
  for (summary in summaries %||% list()) {
    id <- mm_scalar_text(summary$term_id)
    syntax <- mm_scalar_text(summary$source_syntax)
    if (nzchar(id) && nzchar(syntax)) labels[[id]] <- syntax
  }
  labels
}

mm_change_term_label <- function(term, labels, source_syntax = NULL) {
  syntax <- mm_scalar_text(source_syntax)
  if (nzchar(syntax)) return(syntax)
  mapped <- labels[[term]]
  if (!is.null(mapped)) return(mapped)
  if (nzchar(term)) term else "(term)"
}

# Covariance family fields cross the boundary either as plain strings
# ("scalar") or as a single-key object with parameters
# (list(reduced_rank = list(rank = 0))); render the latter as
# "reduced_rank(rank = 0)".
mm_change_covariance_text <- function(value) {
  if (is.null(value)) return("(unknown)")
  if (is.character(value)) return(paste(value, collapse = ":"))
  if (is.list(value) && length(value) == 1L && !is.null(names(value)) &&
      nzchar(names(value)[[1L]])) {
    family <- names(value)[[1L]]
    inner <- value[[1L]]
    if (is.list(inner) && length(inner) && !is.null(names(inner))) {
      fields <- vapply(names(inner), function(nm) {
        sprintf("%s = %s", nm, mm_scalar_text(inner[[nm]]))
      }, character(1))
      return(sprintf("%s(%s)", family, paste(fields, collapse = ", ")))
    }
    return(family)
  }
  mm_scalar_text(value, "(unknown)")
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
