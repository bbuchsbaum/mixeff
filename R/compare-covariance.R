#' Compare covariance parameterizations for current random terms
#'
#' `compare_covariance()` is a compact alternate view of the same upstream
#' random-term cards used by [explain_model()] and [random_options()]. For each
#' current random-term card, it lays out the full, diagonal, and scalar
#' covariance families without ranking them.
#'
#' @param spec An `mm_spec` from [compile_model()] or, in later phases, an
#'   `mm_fit`.
#'
#' @return An object of class `mm_compare_covariance` with a `table` data
#'   frame and the upstream cards it was derived from.
#'
#' @export
compare_covariance <- function(spec) {
  mm_assert_compiled_spec(spec)
  audit <- audit_design(spec)
  cards <- audit$random_term_cards %||% list()
  constraints <- audit$cross_card_constraints %||% list()
  rows <- unlist(lapply(cards, mm_compare_covariance_card_rows), recursive = FALSE)
  table <- if (length(rows)) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(
      term_id = character(),
      group = character(),
      basis = character(),
      covariance_family = character(),
      theta_parameters = integer(),
      assumes_zero = character(),
      design_status = character(),
      current = logical(),
      stringsAsFactors = FALSE
    )
  }

  obj <- list(table = table, cards = cards, cross_card_constraints = constraints)
  class(obj) <- "mm_compare_covariance"
  obj
}

#' @method print mm_compare_covariance
#' @export
print.mm_compare_covariance <- function(x, ...) {
  cat("Covariance comparison:\n")
  if (!nrow(x$table)) {
    cat("  none\n")
    return(invisible(x))
  }

  for (i in seq_len(nrow(x$table))) {
    row <- x$table[i, , drop = FALSE]
    marker <- if (isTRUE(row$current[[1L]])) " <- current" else ""
    cat(sprintf(
      "  %s / %s / %s%s\n",
      row$term_id[[1L]],
      row$group[[1L]],
      row$covariance_family[[1L]],
      marker
    ))
    cat(sprintf("    basis:            %s\n", row$basis[[1L]]))
    cat(sprintf("    theta parameters: %s\n", row$theta_parameters[[1L]]))
    cat(sprintf("    assumes zero:     %s\n", row$assumes_zero[[1L]]))
    cat(sprintf("    design status:    %s\n", row$design_status[[1L]]))
  }
  if (length(x$cross_card_constraints %||% list())) {
    cat("Cross-card constraints:\n")
    for (constraint in x$cross_card_constraints) {
      cat(sprintf("  %s\n", constraint$reason %||% "not recorded"))
    }
  }
  invisible(x)
}

mm_compare_covariance_card_rows <- function(card) {
  block <- card$blocks[[1L]]
  basis <- unlist(block$basis %||% list(), use.names = FALSE)
  p <- length(basis)
  current <- mm_covariance_label(block$covariance)
  families <- c("full", "diagonal", "scalar")

  lapply(families, function(family) {
    data.frame(
      term_id = card$term_id,
      group = mm_group_label(card$group),
      basis = if (length(basis)) paste(basis, collapse = ", ") else "(none)",
      covariance_family = family,
      theta_parameters = mm_covariance_theta_count(family, p),
      assumes_zero = mm_covariance_zero_assumption(family, basis),
      design_status = card$design_support$status %||% "not_assessed",
      current = identical(family, current) ||
        (identical(family, "diagonal") && identical(current, "diagonal")),
      stringsAsFactors = FALSE
    )
  })
}

mm_covariance_theta_count <- function(family, p) {
  switch(
    family,
    full = p * (p + 1L) / 2L,
    diagonal = p,
    scalar = if (p > 0L) 1L else 0L,
    NA_integer_
  )
}

mm_covariance_zero_assumption <- function(family, basis) {
  p <- length(basis)
  if (p <= 1L) return("none")
  switch(
    family,
    full = "none",
    diagonal = "off-diagonal covariances",
    scalar = "off-diagonal covariances",
    "unknown"
  )
}
