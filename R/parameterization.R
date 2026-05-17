#' Inspect covariance parameterization
#'
#' `parameterization()` exposes the fitted theta/Lambda mapping recorded in the
#' compiler artifact. It is the R table view of the upstream theta-map and
#' covariance-parameter trace records.
#'
#' @param fit A compiled `mm_spec` or fitted `mm_fit`.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_theta_map` object with a data-frame `table` and raw trace
#'   records.
#'
#' @export
parameterization <- function(fit, ...) {
  UseMethod("parameterization")
}

#' @rdname parameterization
#' @export
parameterization.mm_compiled <- function(fit, ...) {
  artifact <- mm_compiled_artifact(fit)
  traces <- artifact$covariance_parameter_traces %||% list()
  rows <- lapply(traces, mm_parameterization_trace_row)
  table <- if (length(rows)) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  } else {
    mm_parameterization_empty_table()
  }
  obj <- list(table = table, traces = traces, theta_maps = artifact$theta_maps %||% list())
  class(obj) <- "mm_theta_map"
  obj
}

#' @method print mm_theta_map
#' @export
print.mm_theta_map <- function(x, ...) {
  cat("Covariance parameterization:\n")
  if (!nrow(x$table)) {
    cat("  no theta parameters\n")
    return(invisible(x))
  }
  cols <- intersect(
    c("term_id", "group", "source_syntax", "covariance_family",
      "theta_name", "theta_value", "theta_status", "varcorr_entries"),
    names(x$table)
  )
  print(x$table[, cols, drop = FALSE], row.names = FALSE)
  hidden <- setdiff(names(x$table), cols)
  if (length(hidden)) {
    cat(sprintf("Full theta/Lambda columns available in `x$table` (%d hidden).\n",
                length(hidden)))
  }
  invisible(x)
}

mm_parameterization_trace_row <- function(trace) {
  theta <- trace$theta %||% list()
  lambda <- trace$lambda %||% list()
  data.frame(
    term_id = mm_scalar_text(trace$term_id),
    group = mm_scalar_text(trace$group),
    source_syntax = mm_scalar_text(trace$source_syntax),
    covariance_family = mm_scalar_text(trace$covariance_family),
    user_basis = mm_list_text(trace$user_basis),
    optimizer_basis = mm_list_text(trace$optimizer_basis),
    theta_index = as.integer(theta$global_index %||% NA_integer_),
    theta_name = mm_scalar_text(theta$name),
    theta_status = mm_scalar_text(theta$status),
    constraint = mm_scalar_text(theta$constraint),
    theta_value = as.numeric(theta$value %||% NA_real_),
    lambda_row = as.integer(lambda$row %||% NA_integer_),
    lambda_col = as.integer(lambda$col %||% NA_integer_),
    lambda_row_basis = mm_scalar_text(lambda$row_basis),
    lambda_col_basis = mm_scalar_text(lambda$col_basis),
    lambda_value = as.numeric(lambda$value %||% NA_real_),
    varcorr_entries = mm_varcorr_entry_text(trace$varcorr_entries),
    stringsAsFactors = FALSE
  )
}

mm_varcorr_entry_text <- function(entries) {
  entries <- entries %||% list()
  if (!length(entries)) return("")
  parts <- vapply(entries, function(entry) {
    basis <- mm_list_text(entry$basis %||% entry$terms)
    value <- entry$value %||% NA_real_
    sprintf("%s[%s]=%s", mm_scalar_text(entry$kind), basis, signif(value, 6))
  }, character(1))
  paste(parts, collapse = "; ")
}

mm_parameterization_empty_table <- function() {
  data.frame(
    term_id = character(),
    group = character(),
    source_syntax = character(),
    covariance_family = character(),
    user_basis = character(),
    optimizer_basis = character(),
    theta_index = integer(),
    theta_name = character(),
    theta_status = character(),
    constraint = character(),
    theta_value = numeric(),
    lambda_row = integer(),
    lambda_col = integer(),
    lambda_row_basis = character(),
    lambda_col_basis = character(),
    lambda_value = numeric(),
    varcorr_entries = character(),
    stringsAsFactors = FALSE
  )
}
