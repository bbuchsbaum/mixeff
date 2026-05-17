#' Inspect mixeff diagnostics and fit status
#'
#' `diagnostics()` returns the structured diagnostics carried by a compiled
#' spec or fitted model artifact. `fit_status()` is the compact status string
#' recorded by the optimizer certificate for fitted models.
#'
#' @param fit A compiled `mm_spec` or fitted `mm_fit`.
#' @param severity Optional character vector used to filter diagnostics by
#'   severity.
#' @param stage Optional character vector used to filter diagnostics by stage.
#' @param ... Reserved for future methods.
#'
#' @return `diagnostics()` returns an `mm_diagnostics` object containing the
#'   raw diagnostic list and a data-frame view. `fit_status()` returns a
#'   length-one character string.
#'
#' @export
diagnostics <- function(fit, severity = NULL, stage = NULL, ...) {
  UseMethod("diagnostics")
}

#' @rdname diagnostics
#' @export
diagnostics.mm_compiled <- function(fit, severity = NULL, stage = NULL, ...) {
  artifact <- mm_compiled_artifact(fit)
  raw <- mm_artifact_diagnostics(artifact)
  table <- mm_diagnostics_table(raw)
  keep <- seq_along(raw)
  if (!is.null(severity)) {
    keep <- keep[table$severity[keep] %in% severity]
  }
  if (!is.null(stage)) {
    keep <- keep[table$stage[keep] %in% stage]
  }
  table <- table[keep, , drop = FALSE]
  raw <- raw[keep]
  rownames(table) <- NULL

  out <- list(
    diagnostics = raw,
    table = table,
    severity = severity,
    stage = stage
  )
  class(out) <- "mm_diagnostics"
  out
}

#' @rdname diagnostics
#' @export
fit_status <- function(fit, ...) {
  UseMethod("fit_status")
}

#' @rdname diagnostics
#' @export
fit_status.mm_fit <- function(fit, ...) {
  fit$fit_status %||%
    fit$artifact$optimizer_certificate$status %||%
    "not_assessed"
}

#' @rdname diagnostics
#' @export
fit_status.mm_compiled <- function(fit, ...) {
  fit$artifact$optimizer_certificate$status %||% "not_assessed"
}

#' @method print mm_diagnostics
#' @export
print.mm_diagnostics <- function(x, ...) {
  cat("Diagnostics:\n")
  if (!nrow(x$table)) {
    cat("  none\n")
    return(invisible(x))
  }
  show <- x$table[, intersect(c("code", "severity", "stage", "affected_terms"),
                              names(x$table)), drop = FALSE]
  show <- unique(show)
  print(show, row.names = FALSE)
  if ("message" %in% names(x$table)) {
    cat("\nMessages:\n")
    messages <- unique(x$table[, intersect(c("code", "message"), names(x$table)),
                               drop = FALSE])
    for (i in seq_len(nrow(messages))) {
      msg <- messages$message[[i]]
      code <- messages$code[[i]]
      wrapped <- strwrap(msg, width = 78, exdent = 4)
      cat(sprintf("  %s: %s\n", code, paste(wrapped, collapse = "\n    ")))
    }
  }
  invisible(x)
}

mm_compiled_artifact <- function(x) {
  if (!inherits(x, "mm_compiled") || !is.list(x$artifact)) {
    mm_abort(
      message = "`x` must be an `mm_spec` or `mm_fit` carrying a parsed artifact.",
      class = "mm_schema_error",
      input = x
    )
  }
  x$artifact
}

mm_artifact_diagnostics <- function(artifact) {
  c(
    artifact$diagnostics %||% list(),
    artifact$design_audit$diagnostics %||% list(),
    artifact$optimizer_certificate$diagnostics %||% list()
  )
}

mm_diagnostics_table <- function(raw) {
  if (!length(raw)) {
    return(data.frame(
      code = character(),
      severity = character(),
      stage = character(),
      message = character(),
      affected_terms = character(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(raw, function(d) {
    data.frame(
      code = mm_scalar_text(d$code),
      severity = mm_scalar_text(d$severity),
      stage = mm_scalar_text(d$stage),
      message = mm_scalar_text(d$message),
      affected_terms = mm_list_text(d$affected_terms),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mm_scalar_text <- function(x, default = "") {
  if (is.null(x)) return(default)
  if (is.list(x)) {
    flat <- unlist(x, use.names = FALSE)
    if (!length(flat)) return(default)
    return(paste(as.character(flat), collapse = ":"))
  }
  if (!length(x)) return(default)
  as.character(x[[1L]])
}

mm_list_text <- function(x, default = "") {
  if (is.null(x)) return(default)
  flat <- unlist(x, use.names = FALSE)
  if (!length(flat)) default else paste(as.character(flat), collapse = ", ")
}
