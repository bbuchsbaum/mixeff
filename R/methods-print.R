#' @method print mm_lmm
#' @export
print.mm_lmm <- function(x, ...) {
  schema <- x$schema %||% mm_object_schema(x$artifact)
  cat(sprintf("Linear mixed model fit by %s\n", if (x$REML) "REML" else "ML"))
  cat(sprintf("Formula: %s\n", deparse1(x$formula)))
  cat(sprintf("Fit status: %s\n", x$fit_status))
  cat(mm_print_optimizer_line(x))
  cat(sprintf(
    "Artifact: %s v%s; crate: %s\n",
    schema$schema_name %||% "not_recorded",
    schema$schema_version %||% "not_recorded",
    schema$crate_version %||% "not_recorded"
  ))
  cat(sprintf(
    "nobs: %d, sigma: %.6g, logLik: %.6g\n",
    x$nobs, x$sigma, x$logLik
  ))
  cat("Fixed effects:\n")
  print(signif(fixef(x), 6))
  singular_lines <- mm_singular_render_lines(x)
  if (length(singular_lines)) {
    cat("\nFitted covariance state:\n")
    cat(paste(singular_lines, collapse = "\n"))
    cat("\n")
  }
  cat("Audit verbs: audit(), diagnostics(), inference_table(), model_report()\n")
  invisible(x)
}

# Render the singular-fit summary per PRD §9.5.6: describe effective rank,
# point to changes() and random_options(), never recommend a different
# model spelling. Returns character() when the fit is not singular.
mm_singular_render_lines <- function(fit) {
  if (!isTRUE(is_singular(fit))) return(character())
  summaries <- fit$artifact$effective_covariance %||% list()
  lines <- c("The fitted covariance matrix is rank-deficient.")
  for (summary in summaries) {
    requested <- summary$requested_rank %||% NA
    supported <- summary$supported_rank %||% NA
    term_id   <- summary$term_id %||% "term"
    if (!is.na(requested) && !is.na(supported) && requested != supported) {
      lines <- c(lines, sprintf(
        "  %s: requested rank %s; fitted effective rank %s.",
        term_id, requested, supported
      ))
    }
  }
  group <- mm_singular_first_group(fit)
  lines <- c(
    lines,
    "Use changes(fit) to see which dimension was unsupported.",
    sprintf("Use random_options(spec, group = %s) to inspect lower-dimensional covariance choices.",
            group)
  )
  lines
}

# Find the first random-effect grouping factor recorded on the fit, so the
# singular-render pointer can name a concrete group instead of a placeholder.
# Falls back to "<group>" if no card is recorded.
mm_singular_first_group <- function(fit) {
  cards <- fit$artifact$design_audit$random_term_cards %||% list()
  for (card in cards) {
    g <- card$group
    if (is.character(g) && length(g) == 1L) return(g)
    if (is.list(g)) {
      if (!is.null(g$single$name))     return(g$single$name)
      if (!is.null(g$cell$names))      return(paste(unlist(g$cell$names), collapse = ":"))
      if (!is.null(g$interaction$names)) return(paste(unlist(g$interaction$names), collapse = ":"))
    }
  }
  ranef_names <- names(fit$random_effects %||% list())
  if (length(ranef_names)) return(ranef_names[[1L]])
  "<group>"
}

mm_print_optimizer_line <- function(x) {
  cert <- x$artifact$optimizer_certificate %||% list()
  if (!length(cert)) return("")
  optimizer  <- mm_scalar_text(cert$optimizer_name, "not_recorded")
  iterations <- mm_scalar_text(cert$iterations,     "not_recorded")
  objective  <- cert$objective_value
  obj_text <- if (is.numeric(objective) && length(objective) == 1L && is.finite(objective)) {
    sprintf("%.6g", objective)
  } else {
    "not_recorded"
  }
  sprintf("Optimizer: %s; iterations: %s; objective: %s\n",
          optimizer, iterations, obj_text)
}

#' @method print mm_varcorr
#' @export
print.mm_varcorr <- function(x, ...) {
  cat("Variance components:\n")
  if (nrow(x$table)) {
    out <- x$table
    out$variance <- signif(out$variance, 6)
    out$std_dev <- signif(out$std_dev, 6)
    boundary <- if (!is.null(out$boundary)) isTRUE(any(out$boundary)) else FALSE
    if (!is.null(out$boundary)) {
      out$note <- ifelse(out$boundary, "[boundary]", "")
      out$boundary <- NULL
    }
    print(out, row.names = FALSE)
    if (boundary) {
      cat("[boundary]: variance component is at the boundary of the parameter space.\n")
    }
  } else {
    cat("  none\n")
  }
  if (is.finite(x$residual_sd)) {
    cat(sprintf("Residual std. dev.: %.6g\n", x$residual_sd))
  }
  invisible(x)
}

#' @method print mm_ranef
#' @export
print.mm_ranef <- function(x, ...) {
  cat("Random effects:\n")
  if (!length(x)) {
    cat("  none\n")
    return(invisible(x))
  }
  for (nm in names(x)) {
    cat(sprintf("$%s\n", nm))
    print(x[[nm]])
  }
  invisible(x)
}

#' @method print mm_coef
#' @export
print.mm_coef <- function(x, ...) {
  cat("Conditional coefficients:\n")
  if (!length(x)) {
    cat("  none\n")
    return(invisible(x))
  }
  for (nm in names(x)) {
    cat(sprintf("$%s\n", nm))
    print(x[[nm]])
  }
  invisible(x)
}
