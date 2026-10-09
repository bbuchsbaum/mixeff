#' @method print mm_lmm
#' @export
print.mm_lmm <- function(x, ...) {
  cat(sprintf("Linear mixed model fit by %s\n", if (x$REML) "REML" else "ML"))
  cat(sprintf("Formula: %s\n", deparse1(x$formula)))
  cat(sprintf("Fit status: %s\n", x$fit_status))
  cat(mm_print_optimizer_line(x))
  # Artifact/crate provenance is developer metadata, not model output;
  # reproducibility(fit) reports it (UX bar: print() shows nothing lme4's
  # print would not need explained).
  cat(sprintf(
    "nobs: %d, sigma: %.6g, logLik: %.6g\n",
    x$nobs, x$sigma, x$logLik
  ))
  cat("Fixed effects:\n")
  beta_shown <- fixef(x)
  aliased <- mm_aliased_coefficients(x)
  aliased <- intersect(aliased, names(beta_shown))
  if (length(aliased)) beta_shown[aliased] <- NA_real_
  print(signif(beta_shown, 6))
  if (length(aliased)) {
    cat(sprintf(
      "Aliased (not separately estimable, shown as NA): %s. See the design note in summary().\n",
      paste(aliased, collapse = ", ")
    ))
  }
  singular_lines <- mm_singular_render_lines(x)
  if (length(singular_lines)) {
    cat("\nFitted covariance state:\n")
    cat(paste(singular_lines, collapse = "\n"))
    cat("\n")
  }
  cat("Audit verbs: audit(), diagnostics(), inference_table(), model_report()\n")
  invisible(x)
}

#' @method print mm_glmm
#' @export
print.mm_glmm <- function(x, ...) {
  cat("Generalized linear mixed model fit\n")
  cat(sprintf("Formula: %s\n", deparse1(x$formula)))
  cat(sprintf("Family/link: %s/%s\n", x$family$family, x$family$link))
  cat(sprintf("Method: %s (nAGQ = %d)\n", x$method, x$nAGQ))
  cat(sprintf("Fit status: %s\n", x$fit_status))
  cat(mm_print_optimizer_line(x))
  # Artifact/crate provenance is developer metadata, not model output;
  # reproducibility(fit) reports it (UX bar: print() shows nothing lme4's
  # print would not need explained).
  cat(sprintf(
    "nobs: %d, dispersion: %.6g, logLik: %.6g\n",
    x$nobs, x$dispersion, x$logLik
  ))
  cat("Fixed effects:\n")
  beta_shown <- fixef(x)
  aliased <- mm_aliased_coefficients(x)
  aliased <- intersect(aliased, names(beta_shown))
  if (length(aliased)) beta_shown[aliased] <- NA_real_
  print(signif(beta_shown, 6))
  if (length(aliased)) {
    cat(sprintf(
      "Aliased (not separately estimable, shown as NA): %s. See the design note in summary().\n",
      paste(aliased, collapse = ", ")
    ))
  }
  cat("Audit verbs: audit(), diagnostics(), model_report()\n")
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
  lines <- c(lines, "Use changes(fit) to see which dimension was unsupported.")
  # Only advertise random_options() when it can actually run: it enumerates
  # slope-bearing spellings, so a model with no slope candidate anywhere
  # (intercept-only random terms and no non-intercept fixed effect) makes it
  # refuse -- a printed pointer must not error on the very fit that printed it.
  if (mm_singular_has_slope_candidate(fit)) {
    lines <- c(lines, sprintf(
      "Use random_options(spec, group = %s) to inspect lower-dimensional covariance choices.",
      group
    ))
  }
  lines
}

# Mirror of mm_default_slope()'s candidate search, from fit-artifact facts
# alone: a slope on any random-term card, else any non-intercept fixed term.
mm_singular_has_slope_candidate <- function(fit) {
  cards <- fit$artifact$design_audit$random_term_cards %||% list()
  for (card in cards) {
    slopes <- unlist(
      lapply(card$blocks %||% list(), function(block) block$slopes %||% list()),
      use.names = FALSE
    )
    if (length(slopes)) return(TRUE)
  }
  fixed <- unlist(fit$artifact$semantic_model$fixed_terms %||% list(),
                  use.names = FALSE)
  length(setdiff(fixed, "1")) > 0L
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
print.mm_varcorr <- function(x, digits = max(3, getOption("digits") - 2),
                             ...) {
  # lme4's layout (Groups / Name / Std.Dev. / Corr); mixeff's own markers
  # ([boundary], [design_weak_identifiability]) follow as footnotes. The
  # full-precision long table stays in `x$table`.
  if (mm_varcorr_is_internal(x)) {
    tbl <- x$table
    cat("Variance components:\n")
    if (nrow(tbl)) print(tbl, row.names = FALSE) else cat("  none\n")
    if (is.finite(x$residual_sd %||% NA_real_)) {
      cat(sprintf("Residual std. dev.: %.6g\n", x$residual_sd))
    }
    return(invisible(x))
  }
  weak_groups <- attr(x, "mm_design_weak_identifiability_groups") %||% character()
  formatted <- mm_format_varcorr(x, digits = digits)
  if (!nrow(formatted)) {
    cat("Random effects: none\n")
    return(invisible(x))
  }
  tbl <- x$table
  marks <- character(nrow(formatted))
  boundary <- FALSE
  lme4_boundary <- attr(x, "mm_boundary")
  if (!is.null(lme4_boundary)) {
    hit <- which(lme4_boundary %in% TRUE)
    if (length(hit) && length(lme4_boundary) <= length(marks)) {
      marks[hit] <- "[boundary]"
      boundary <- TRUE
    }
  } else if (is.data.frame(tbl) && nrow(tbl) && !is.null(tbl$boundary)) {
    n_re <- nrow(tbl)
    hit <- which(tbl$boundary %in% TRUE)
    if (length(hit) && n_re <= length(marks)) {
      marks[hit] <- "[boundary]"
      boundary <- TRUE
    }
  }
  if (any(nzchar(marks))) {
    formatted <- cbind(formatted, " " = marks)
  }
  print(formatted, quote = FALSE, ...)
  if (boundary) {
    cat("[boundary]: variance component is at the boundary of the parameter space.\n")
  }
  if (length(weak_groups)) {
    groups <- paste(sprintf("`%s`", weak_groups), collapse = ", ")
    cat(sprintf(
      paste0(
        "[design_weak_identifiability]: random intercept variance for %s ",
        "is weakly interpretable because its grouping indicators are ",
        "aliased with the fixed-effect design.\n"
      ),
      groups
    ))
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
