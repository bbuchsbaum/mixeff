#' Inspect nearby random-effect spellings for one grouping factor
#'
#' `random_options()` is an opt-in map over nearby random-effect structures
#' for a grouping factor. It recompiles each displayed spelling through the
#' same upstream audit path as [compile_model()], so support facts and block
#' meanings come from Rust-authored `RandomTermCard` records.
#'
#' @param spec An `mm_spec` from [compile_model()] or, in later phases, an
#'   `mm_fit`.
#' @param group Grouping factor to inspect. May be supplied bare
#'   (`group = subject`) or as a string.
#' @param slope Optional slope variable to use for nearby slope-bearing
#'   spellings. When omitted, the function uses the first current random
#'   slope for `group`, then any scope-note fixed effect for `group`, then
#'   the first non-intercept fixed effect.
#'
#' @return An object of class `mm_random_options` with an `options` data frame,
#'   the upstream candidate `cards`, and the candidate audit reports.
#'
#' @export
random_options <- function(spec, group, slope = NULL) {
  mm_assert_compiled_spec(spec)
  group_label <- if (missing(group)) {
    NULL
  } else {
    mm_expr_label(substitute(group), allow_null = FALSE)
  }
  slope_label <- if (missing(slope) || is.null(slope)) {
    NULL
  } else {
    mm_expr_label(substitute(slope), allow_null = TRUE)
  }

  if (is.null(group_label) || !nzchar(group_label)) {
    mm_abort(
      message = "`group` must name a grouping factor.",
      class = "mm_schema_error",
      input = group
    )
  }

  data <- mm_spec_model_frame(spec)
  audit <- audit_design(spec)
  current_cards <- mm_cards_for_group(audit$random_term_cards, group_label)
  if (!length(current_cards)) {
    mm_abort(
      message = sprintf("No random-effect cards found for group `%s`.", group_label),
      class = "mm_schema_error",
      input = group_label
    )
  }

  slope_label <- slope_label %||% mm_default_slope(spec, audit, group_label)
  if (is.null(slope_label) || !nzchar(slope_label)) {
    mm_abort(
      message = sprintf("No slope variable is available for group `%s`.", group_label),
      class = "mm_schema_error",
      input = group_label
    )
  }
  if (!slope_label %in% names(data)) {
    mm_abort(
      message = sprintf("Slope variable `%s` is not available in `spec$model_frame`.", slope_label),
      class = "mm_data_error",
      input = slope_label
    )
  }

  candidates <- mm_random_option_candidates(group_label, slope_label)
  current_fragment <- mm_current_random_fragment(current_cards)
  rows <- vector("list", length(candidates))
  cards <- vector("list", length(candidates))
  reports <- vector("list", length(candidates))

  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    candidate_spec <- compile_model(
      mm_candidate_formula(spec, candidate$fragment),
      data
    )
    candidate_audit <- audit_design(candidate_spec)
    candidate_cards <- mm_cards_for_group(candidate_audit$random_term_cards, group_label)

    rows[[i]] <- mm_random_option_row(
      candidate,
      candidate_cards,
      candidate_audit$cross_card_constraints,
      current = identical(candidate$fragment, current_fragment)
    )
    cards[[candidate$key]] <- candidate_cards
    reports[[candidate$key]] <- candidate_audit$report
  }

  options <- do.call(rbind, rows)
  rownames(options) <- NULL

  out <- list(
    group = group_label,
    slope = slope_label,
    options = options,
    cards = cards,
    constraints = lapply(reports, function(report) report$cross_card_constraints %||% list()),
    reports = reports
  )
  class(out) <- "mm_random_options"
  out
}

#' @method print mm_random_options
#' @export
print.mm_random_options <- function(x, ...) {
  cat(sprintf("Random-effect options for group: %s\n", x$group))
  cat("Current model:\n")
  current <- x$options[x$options$current, , drop = FALSE]
  if (nrow(current)) {
    for (i in seq_len(nrow(current))) {
      cat(sprintf("  %s <- this is what you wrote\n", current$formula[[i]]))
      cat(sprintf("  %s\n", current$plain_meaning[[i]]))
    }
  } else {
    cat("  (current spelling not in nearby map)\n")
  }

  cat("Nearby options:\n")
  for (i in seq_len(nrow(x$options))) {
    row <- x$options[i, , drop = FALSE]
    marker <- if (isTRUE(row$current[[1L]])) " <- this is what you wrote" else ""
    note <- if (nzchar(row$note[[1L]])) sprintf(" (%s)", row$note[[1L]]) else ""
    cat(sprintf("  %s%s%s\n", row$formula[[1L]], marker, note))
    cat(sprintf("    varying coefficients: %s\n", row$varying_coefficients[[1L]]))
    cat(sprintf("    covariance family:    %s\n", row$covariance_family[[1L]]))
    cat(sprintf("    theta parameters:     %s\n", row$theta_parameters[[1L]]))
    cat(sprintf("    design status:        %s\n", row$design_status[[1L]]))
    cat(sprintf("    plain meaning:        %s\n", row$plain_meaning[[1L]]))
  }
  invisible(x)
}

mm_random_option_candidates <- function(group, slope) {
  list(
    list(key = "punt",
         fragment = sprintf("(1 | %s)", group),
         note = ""),
    list(key = "slope_only",
         fragment = sprintf("(0 + %s | %s)", slope, group),
         note = ""),
    list(key = "split_uncorrelated",
         fragment = sprintf("(1 | %s) + (0 + %s | %s)", group, slope, group),
         note = ""),
    list(key = "double_bar_synonym",
         fragment = sprintf("(1 + %s || %s)", slope, group),
         note = ""),
    list(key = "full",
         fragment = sprintf("(1 + %s | %s)", slope, group),
         note = "")
  )
}

mm_random_option_row <- function(candidate, cards, constraints, current) {
  blocks <- unlist(lapply(cards, function(card) card$blocks), recursive = FALSE)
  basis <- unique(unlist(lapply(blocks, function(block) {
    block$basis %||% list()
  }), use.names = FALSE))
  varying <- basis
  varying[tolower(varying) == "intercept"] <- "intercept"
  covariances <- unique(vapply(blocks, function(block) {
    mm_covariance_label(block$covariance)
  }, character(1)))
  theta <- sum(as.integer(vapply(blocks, function(block) {
    block$theta_parameters %||% 0L
  }, integer(1))))
  status <- unique(vapply(cards, function(card) {
    card$design_support$status %||% "not_assessed"
  }, character(1)))
  meaning <- mm_option_plain_meaning(cards, constraints)

  data.frame(
    formula = candidate$fragment,
    varying_coefficients = if (length(varying)) paste(varying, collapse = ", ") else "(none)",
    covariance_family = mm_option_covariance_family(candidate$key, covariances),
    theta_parameters = theta,
    design_status = paste(status, collapse = "; "),
    plain_meaning = meaning,
    note = candidate$note,
    current = isTRUE(current),
    stringsAsFactors = FALSE
  )
}

mm_option_plain_meaning <- function(cards, constraints) {
  block_text <- unlist(lapply(cards, function(card) {
    vapply(card$blocks, function(block) block$english, character(1))
  }), use.names = FALSE)
  constraint_text <- vapply(constraints %||% list(), function(constraint) {
    constraint$reason %||% ""
  }, character(1))
  text <- unique(c(block_text, constraint_text[nzchar(constraint_text)]))
  if (!length(text)) "(not assessed)" else paste(text, collapse = " ")
}

mm_option_covariance_family <- function(key, covariances) {
  if (key %in% c("split_uncorrelated", "double_bar_synonym")) {
    return("diagonal via separate blocks")
  }
  paste(covariances, collapse = "; ")
}

mm_current_random_fragment <- function(cards) {
  originals <- vapply(cards, function(card) card$original_fragment, character(1))
  originals <- unique(originals)
  if (length(originals) == 1L) return(originals[[1L]])
  paste(originals, collapse = " + ")
}

mm_default_slope <- function(spec, audit, group) {
  cards <- mm_cards_for_group(audit$random_term_cards, group)
  card_slopes <- unique(unlist(lapply(cards, function(card) {
    unlist(lapply(card$blocks, function(block) block$slopes %||% list()),
           use.names = FALSE)
  }), use.names = FALSE))
  if (length(card_slopes)) return(card_slopes[[1L]])

  scope_notes <- audit$diagnostics[vapply(audit$diagnostics, function(d) {
    identical(d$code, "scope_note") && identical(d$payload$group, group)
  }, logical(1))]
  fixed_effects <- unique(vapply(scope_notes, function(d) {
    d$payload$fixed_effect %||% ""
  }, character(1)))
  fixed_effects <- fixed_effects[nzchar(fixed_effects)]
  if (length(fixed_effects)) return(fixed_effects[[1L]])

  fixed <- unlist(spec$artifact$semantic_model$fixed_terms %||% list(), use.names = FALSE)
  fixed <- setdiff(fixed, "1")
  fixed <- fixed[fixed %in% names(mm_spec_model_frame(spec))]
  if (length(fixed)) fixed[[1L]] else NULL
}

mm_candidate_formula <- function(spec, random_fragment) {
  semantic <- spec$artifact$semantic_model
  response <- semantic$response %||% all.vars(spec$formula)[[1L]]
  fixed <- unlist(semantic$fixed_terms %||% list("1"), use.names = FALSE)
  fixed <- fixed[nzchar(fixed)]
  fixed_rhs <- if (length(fixed)) paste(fixed, collapse = " + ") else "1"
  stats::as.formula(sprintf("%s ~ %s + %s", response, fixed_rhs, random_fragment))
}

mm_expr_label <- function(expr, allow_null = FALSE) {
  if (is.null(expr)) return(if (allow_null) NULL else "")
  if (is.character(expr) && length(expr) == 1L) return(expr)
  deparse1(expr)
}

mm_assert_compiled_spec <- function(spec) {
  if (!inherits(spec, c("mm_spec", "mm_fit"))) {
    mm_abort(
      message = "`spec` must be an `mm_spec` (from compile_model) or an `mm_fit` (from lmm).",
      class = "mm_schema_error",
      input = spec
    )
  }
  invisible(TRUE)
}

mm_spec_model_frame <- function(spec) {
  data <- spec$model_frame
  if (!is.data.frame(data)) {
    mm_abort(
      message = "`spec` does not carry a model frame; re-run compile_model() with this package version.",
      class = "mm_schema_error",
      input = spec
    )
  }
  data
}

mm_cards_for_group <- function(cards, group) {
  cards[vapply(cards, function(card) identical(mm_group_label(card$group), group), logical(1))]
}
