#' Explain the random-effects structure of a compiled model
#'
#' `explain_model()` renders the random-effects guidance surface for an
#' `mm_spec` returned by [compile_model()] or, in later phases, an `mm_fit`.
#' It formats the upstream `RandomTermCard` and diagnostic payloads; Rust
#' remains the source of truth for per-block English wording and design facts.
#'
#' @param spec An `mm_spec` produced by [compile_model()] or an `mm_fit`.
#'
#' @return An object of class `mm_explanation` carrying:
#' \describe{
#'   \item{`text`}{the rendered explanation text}
#'   \item{`cards`}{the upstream random-term cards}
#'   \item{`cross_card_constraints`}{report-level constraints between cards}
#'   \item{`diagnostics`}{the upstream diagnostics used for design notes}
#'   \item{`report`}{the parsed upstream `ModelAuditReport`}
#' }
#'
#' @section Errors:
#' Raises an `mm_schema_error` if `spec` is not an `mm_spec`/`mm_fit` or does
#' not carry a valid compiled artifact.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   y = rnorm(20),
#'   t = rep(0:3, 5),
#'   s = factor(rep(1:5, each = 4))
#' )
#' explain_model(compile_model(y ~ t + (1 | s), df))
#' }
#'
#' @seealso [compile_model()], [audit_design()].
#'
#' @export
explain_model <- function(spec) {
  if (!inherits(spec, c("mm_spec", "mm_fit"))) {
    mm_abort(
      message = "`spec` must be an `mm_spec` (from compile_model) or an `mm_fit` (from lmm).",
      class = "mm_schema_error",
      input = spec
    )
  }

  audit <- audit_design(spec)
  out <- list(
    text                   = mm_explanation_text(spec, audit),
    cards                  = audit$random_term_cards,
    cross_card_constraints = audit$cross_card_constraints,
    diagnostics            = audit$diagnostics,
    report                 = audit$report
  )
  class(out) <- "mm_explanation"
  out
}

#' @method format mm_explanation
#' @export
format.mm_explanation <- function(x, ...) {
  x$text
}

#' @method print mm_explanation
#' @export
print.mm_explanation <- function(x, ...) {
  cat(x$text)
  if (!grepl("\n$", x$text)) cat("\n")
  invisible(x)
}

# Emit the pre-fit explanation on the message stream (rlang::inform), not
# stdout: suppressMessages() and knitr's message=FALSE must be able to quiet
# the automatic printout, while an explicit print(explain_model(spec)) still
# writes to stdout like any other print method.
mm_inform_explanation <- function(spec) {
  mm_inform(format(explain_model(spec)), class = "mm_explanation_notice")
}

mm_explanation_text <- function(spec, audit) {
  cards <- audit$random_term_cards %||% list()
  constraints <- audit$cross_card_constraints %||% list()
  diagnostics <- audit$diagnostics %||% list()

  lines <- c(
    "Random effects explanation:",
    sprintf("  formula: %s", mm_spec_formula_text(spec)),
    ""
  )

  if (!length(cards)) {
    lines <- c(lines, "Random effects:", "  none")
  } else {
    lines <- c(lines, "Random effects:")
    lines <- c(lines, mm_group_block_summaries(cards))
    for (card in cards) {
      lines <- c(lines, mm_card_lines(card))
    }
  }

  if (length(constraints)) {
    lines <- c(lines, "", "Relationship between blocks:")
    for (constraint in constraints) {
      lines <- c(lines, mm_constraint_lines(constraint, cards))
    }
  }

  design_notes <- mm_design_note_lines(diagnostics)
  if (length(design_notes)) {
    lines <- c(lines, "", "Design notes:", design_notes)
  }

  repairs <- mm_repair_lines(diagnostics)
  if (length(repairs)) {
    lines <- c(lines, "", "Possible repairs, not applied automatically:", repairs)
  }

  fit_notes <- mm_fit_note_lines(diagnostics)
  if (length(fit_notes)) {
    lines <- c(lines, "", "Fit notes:", fit_notes)
  }

  singularity <- mm_singularity_lines(audit$report)
  if (length(singularity)) {
    lines <- c(lines, "", "Fitted covariance state:", singularity)
  }

  paste(lines, collapse = "\n")
}

mm_spec_formula_text <- function(spec) {
  artifact <- spec$artifact %||% list()
  artifact$effective_formula %||% artifact$requested_formula %||% deparse1(spec$formula)
}

mm_group_block_summaries <- function(cards) {
  groups <- vapply(cards, function(card) mm_group_label(card$group), character(1))
  out <- character()
  for (group in unique(groups)) {
    n <- sum(groups == group)
    if (n > 1L) {
      out <- c(out, sprintf("  %s has %d separate random-effect blocks.", group, n))
    }
  }
  out
}

mm_card_lines <- function(card) {
  group <- mm_group_label(card$group)
  lines <- c(
    sprintf("  %s:", card$term_id),
    sprintf("    wrote:      %s", card$original_fragment),
    sprintf("    canonical:  %s", card$canonical_fragment)
  )

  for (i in seq_along(card$blocks)) {
    block <- card$blocks[[i]]
    prefix <- if (length(card$blocks) > 1L) sprintf("block %d ", i) else ""
    lines <- c(
      lines,
      sprintf(
        "    %snamed form: re(group = %s, intercept = %s, slopes = %s, cov = \"%s\")",
        prefix,
        group,
        if (isTRUE(block$intercept)) "TRUE" else "FALSE",
        mm_slopes_label(block$slopes),
        mm_covariance_label(block$covariance)
      ),
      sprintf("    %sscope:      %s", prefix, block$english),
      sprintf(
        "    %scovariance: %s; theta parameters: %s",
        prefix,
        mm_covariance_label(block$covariance),
        as.character(block$theta_parameters)
      )
    )
  }

  lines <- c(
    lines,
    sprintf("    support:    %s", mm_support_label(card$design_support))
  )

  variation <- mm_variation_label(card$design_support$within_group_variation)
  if (nzchar(variation)) {
    lines <- c(lines, sprintf("    variation:  %s", variation))
  }
  lines
}

mm_group_label <- function(group) {
  if (is.null(group)) return("(unknown)")
  if (is.character(group)) return(group[[1L]])
  if (!is.null(group$single$name)) return(group$single$name)
  if (!is.null(group$cell$names)) {
    return(paste(unlist(group$cell$names, use.names = FALSE), collapse = ":"))
  }
  paste(unlist(group, use.names = FALSE), collapse = ":")
}

mm_slopes_label <- function(slopes) {
  slopes <- unlist(slopes %||% list(), use.names = FALSE)
  if (!length(slopes)) "NULL" else paste(slopes, collapse = ", ")
}

mm_covariance_label <- function(covariance) {
  label <- paste(unlist(covariance %||% "", use.names = FALSE), collapse = ":")
  if (identical(label, "")) "(unknown)" else label
}

mm_support_label <- function(support) {
  if (!is.list(support)) return("(not assessed)")
  parts <- c(support$status %||% "not_assessed")
  if (!is.null(support$group_levels)) {
    parts <- c(parts, sprintf("group levels: %s", support$group_levels))
  }
  if (!is.null(support$min_rows_per_group)) {
    parts <- c(parts, sprintf("min rows/group: %s", support$min_rows_per_group))
  }
  if (!is.null(support$median_rows_per_group)) {
    parts <- c(parts, sprintf("median rows/group: %s", support$median_rows_per_group))
  }
  paste(parts, collapse = "; ")
}

mm_variation_label <- function(variation) {
  if (!length(variation)) return("")
  entries <- vapply(names(variation), function(nm) {
    sprintf("%s=%s", nm, variation[[nm]])
  }, character(1))
  paste(entries, collapse = "; ")
}

mm_constraint_lines <- function(constraint, cards) {
  cards_label <- paste(unlist(constraint$between_cards %||% list(), use.names = FALSE),
                       collapse = " <-> ")
  basis_label <- paste(unlist(constraint$between_basis %||% list(), use.names = FALSE),
                       collapse = " <-> ")
  sprintf("  %s (%s): %s", cards_label, basis_label, constraint$reason)
}

mm_design_note_lines <- function(diagnostics) {
  mm_bucket_advice_lines(diagnostics, "design_note")
}

mm_fit_note_lines <- function(diagnostics) {
  mm_bucket_advice_lines(diagnostics, "fit_note")
}

mm_bucket_advice_lines <- function(diagnostics, bucket) {
  if (!length(diagnostics)) return(character())
  keep <- vapply(diagnostics, function(d) {
    identical(mm_diagnostic_bucket(d$code), bucket)
  }, logical(1))
  diagnostics <- diagnostics[keep]
  if (!length(diagnostics)) return(character())
  unique(vapply(diagnostics, function(d) {
    sprintf("  %s: %s", mm_scalar_text(d$code), d$message)
  }, character(1)))
}

mm_repair_lines <- function(diagnostics) {
  if (!length(diagnostics)) return(character())
  keep <- vapply(diagnostics, function(d) {
    identical(mm_diagnostic_bucket(d$code), "repair")
  }, logical(1))
  refused <- diagnostics[keep]
  if (!length(refused)) return(character())
  pairs <- unlist(
    lapply(refused, function(d) {
      actions <- unlist(d$suggested_actions %||% list(), use.names = FALSE)
      actions <- actions[nzchar(actions)]
      if (!length(actions)) actions <- d$message %||% ""
      actions <- actions[nzchar(actions)]
      lapply(actions, function(a) list(code = mm_scalar_text(d$code), action = a))
    }),
    recursive = FALSE
  )
  if (!length(pairs)) return(character())
  keys <- vapply(pairs, function(p) paste(p$code, p$action, sep = "\x1f"),
                 character(1))
  pairs <- pairs[!duplicated(keys)]
  vapply(seq_along(pairs), function(i) {
    sprintf("  %d. %s: %s", i, pairs[[i]]$code, pairs[[i]]$action)
  }, character(1))
}

mm_singularity_lines <- function(report) {
  sections <- report$sections %||% list()
  effective <- sections[vapply(sections, function(section) {
    identical(section$title, "Effective Covariance")
  }, logical(1))]
  if (!length(effective)) return(character())
  lines <- effective[[1L]]$lines %||% list()
  deficient <- lines[vapply(lines, function(line) {
    grepl("rank-deficient|rank deficient|effective rank", line$detail %||% "",
          ignore.case = TRUE)
  }, logical(1))]
  if (!length(deficient)) return(character())
  vapply(deficient, function(line) {
    sprintf("  %s: %s", line$label, line$detail)
  }, character(1))
}
