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
  raw <- mm_model_diagnostics(fit)
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

mm_model_diagnostics <- function(fit) {
  artifact <- mm_compiled_artifact(fit)
  c(mm_artifact_diagnostics(artifact), mm_r_design_diagnostics(fit))
}

mm_r_design_diagnostics <- function(fit) {
  c(mm_design_weak_identifiability_diagnostics(fit))
}

mm_design_weak_identifiability_diagnostics <- function(
    fit, tolerance = sqrt(.Machine$double.eps)) {
  if (!inherits(fit, "mm_compiled") || !is.data.frame(fit$model_frame)) {
    return(list())
  }
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  if (!length(terms)) return(list())

  X <- tryCatch(
    stats::model.matrix(mm_fixed_formula(fit), data = fit$model_frame),
    error = function(cnd) NULL
  )
  if (is.null(X) || !nrow(X) || !ncol(X)) return(list())
  qr_x <- qr(X, tol = tolerance)
  fixed_rank <- qr_x$rank

  out <- list()
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    if (!mm_random_term_is_intercept_only(term)) next

    group_label <- tryCatch(
      mm_random_term_group_label(fit, term, i),
      error = function(cnd) NA_character_
    )
    if (is.na(group_label) || !nzchar(group_label)) next

    group <- tryCatch(
      mm_group_factor(fit$model_frame, group_label),
      error = function(cnd) NULL
    )
    if (is.null(group)) next

    G <- stats::model.matrix(~ 0 + group)
    fitted <- tryCatch(qr.fitted(qr_x, G), error = function(cnd) NULL)
    if (is.null(fitted)) next

    residual <- G - fitted
    rel <- norm(residual, type = "F") / max(norm(G, type = "F"), .Machine$double.eps)
    if (!is.finite(rel) || rel > tolerance) next

    term_id <- mm_scalar_text(term$id, sprintf("r%d", i - 1L))
    source <- term$source_syntax$text %||% sprintf("(1 | %s)", group_label)
    n_levels <- length(levels(group))
    group_rank <- qr(G, tol = tolerance)$rank
    message <- sprintf(
      paste0(
        "random intercept `%s` is aliased with the fixed-effect design ",
        "(projection residual %.2e); VarCorr(%s) is weakly interpretable. ",
        "Fixed effects and fitted values can remain stable, but this ",
        "variance component should not be treated as a certified lme4-equivalent target."
      ),
      group_label, rel, group_label
    )
    out[[length(out) + 1L]] <- list(
      code = "design_weak_identifiability",
      severity = "warning",
      stage = "design_audit",
      message = message,
      affected_terms = list(term_id),
      suggested_actions = list(
        sprintf("fit a reduced comparison model without `%s`", source),
        sprintf("use a less saturated fixed-effect structure if %s-level variance is substantively important", group_label)
      ),
      payload = list(
        group = group_label,
        term_id = term_id,
        source = source,
        n_levels = n_levels,
        fixed_rank = fixed_rank,
        group_rank = group_rank,
        projection_residual_rel = rel,
        tolerance = tolerance,
        diagnostic_class = "design_weak_identifiability"
      )
    )
  }
  out
}

mm_random_term_is_intercept_only <- function(term) {
  basis <- term$basis %||% list()
  if (!length(basis)) return(TRUE)
  length(basis) == 1L && identical(mm_basis_label(basis[[1L]]), "(Intercept)")
}

mm_design_weak_identifiability_groups <- function(fit) {
  diagnostics <- mm_design_weak_identifiability_diagnostics(fit)
  groups <- vapply(diagnostics, function(d) {
    d$payload$group %||% ""
  }, character(1))
  unique(groups[nzchar(groups)])
}

mm_append_r_design_diagnostic_text <- function(text, diagnostics) {
  design <- diagnostics[vapply(diagnostics, function(d) {
    identical(mm_diagnostic_bucket(d$code), "design_note")
  }, logical(1))]
  if (!length(design)) return(text)
  lines <- unique(vapply(design, function(d) {
    sprintf("  %s: %s", mm_scalar_text(d$code), d$message)
  }, character(1)))
  paste(c(text, "", "Wrapper design notes:", lines), collapse = "\n")
}

mm_diagnostics_table <- function(raw) {
  if (!length(raw)) {
    out <- data.frame(
      code = character(),
      severity = character(),
      stage = character(),
      message = character(),
      affected_terms = character(),
      stringsAsFactors = FALSE
    )
    return(mm_diagnostics_guard(out))
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
  mm_diagnostics_guard(do.call(rbind, rows))
}

# Stable registry pinning every Rust `DiagnosticCode` to an R-side advice
# bucket so explain_model() / audit_design() never silently drop a code.
# Mirrors src/rust/upstream/mixeff-rs/src/compiler/diagnostics.rs; the
# coverage gate in tests/testthat/test-diagnostic-formatters.R fails if
# Rust gains a variant the registry does not classify.
#
# Buckets:
#   "design_note" -> rendered under "Design notes:" in explain_model()
#                    (pre-fit / structural advisories)
#   "repair"      -> rendered under "Possible repairs, not applied
#                    automatically:" (refusals or blocker-class issues that
#                    require user intervention before a fit succeeds)
#   "fit_note"    -> rendered under "Fit notes:" in explain_model()
#                    (post-fit optimizer / inference state)
#   "raw_only"    -> intentionally absent from advice surfaces; the row is
#                    still visible via diagnostics(fit)$table. Use only when
#                    the variant is duplicative of another formatted code or
#                    when the variant is structurally pre-fit cosmetic.
mm_diagnostic_code_registry <- list(
  # --- Pedagogical taxonomy (Rust appended these to model the
  # advice-surface intentionally; the wording is user-facing).
  scope_note                           = list(bucket = "design_note"),
  support_note                         = list(bucket = "design_note"),
  syntax_expansion                     = list(bucket = "design_note"),
  covariance_assumption                = list(bucket = "design_note"),
  structural_refusal                   = list(bucket = "repair"),

  # --- Formula parsing / canonicalisation
  formula_canonicalized                = list(
    bucket    = "raw_only",
    rationale = "Duplicative of syntax_expansion / covariance_assumption which carry user-facing wording for the same canonicalisation event."
  ),
  formula_canonicalization_unsupported = list(bucket = "repair"),

  # --- Semantic IR / random-term audit
  duplicate_random_term                = list(bucket = "design_note"),
  conflicting_covariance               = list(bucket = "design_note"),
  crossing_likely_unintended           = list(bucket = "design_note"),
  random_slope_without_intercept       = list(bucket = "design_note"),
  fixed_random_redundant               = list(bucket = "design_note"),
  repeated_unit_unmodeled              = list(bucket = "design_note"),
  random_slope_unsupported             = list(bucket = "repair"),
  random_effect_few_levels             = list(bucket = "design_note"),
  covariance_too_rich                  = list(bucket = "design_note"),
  covariance_reduced                   = list(bucket = "design_note"),

  # --- Design audit / fixed effects
  fixed_effect_column_missing          = list(bucket = "repair"),
  fixed_effect_rank_deficient          = list(bucket = "design_note"),
  fixed_effect_empty_cell              = list(bucket = "design_note"),

  # --- Numerics / boundary / identifiability
  boundary_parameter                   = list(bucket = "fit_note"),
  near_unit_random_effect_correlation  = list(bucket = "fit_note"),
  binomial_separation                  = list(bucket = "repair"),
  not_identifiable                     = list(bucket = "repair"),
  invalid_agq_request                  = list(bucket = "repair"),

  # --- Optimizer / PIRLS state
  optimizer_not_assessed               = list(
    bucket    = "raw_only",
    rationale = "Pre-fit cosmetic; explain_model() of an unfitted spec already documents the absence of a fit via print() / Audit Summary."
  ),
  optimizer_nonconvergence             = list(bucket = "fit_note"),
  optimizer_recovery                   = list(bucket = "fit_note"),
  pirls_failure                        = list(bucket = "fit_note"),

  # --- Inference / serialization
  inference_unavailable                = list(bucket = "fit_note"),
  serialization_not_assessed           = list(
    bucket    = "raw_only",
    rationale = "Pre-fit serialization state; only relevant when inspecting the diagnostics table directly."
  ),

  # --- Catch-all
  unsupported                          = list(bucket = "repair")
)

mm_r_diagnostic_code_registry <- list(
  design_weak_identifiability = list(bucket = "design_note")
)

# Stable registry for `ResponseDiagnosticReason` (mixeff-rs
# src/model/batch.rs). No R-side surface currently exposes the batch
# engine's per-column response diagnostics; this registry is the
# forward-compat slot so the GLMM batch path (bd-01KRCKCYZ51H7H2WN8C5D7FNGT)
# cannot silently introduce a new variant. The coverage test compares
# names(registry) against the Rust enum.
mm_response_diagnostic_reason_registry <- list(
  non_finite_response = list(
    bucket    = "raw_only",
    rationale = "Batch-engine response diagnostic; no R surface currently exposes ResponseColumnDiagnostic."
  ),
  constant_response = list(
    bucket    = "raw_only",
    rationale = "Batch-engine response diagnostic; no R surface currently exposes ResponseColumnDiagnostic."
  ),
  boundary_theta = list(
    bucket    = "raw_only",
    rationale = "Batch-engine response diagnostic; no R surface currently exposes ResponseColumnDiagnostic."
  ),
  optimizer_failed = list(
    bucket    = "raw_only",
    rationale = "Batch-engine response diagnostic; no R surface currently exposes ResponseColumnDiagnostic."
  ),
  unsupported_mode = list(
    bucket    = "raw_only",
    rationale = "Batch-engine response diagnostic; no R surface currently exposes ResponseColumnDiagnostic."
  ),
  # A streaming sink requested stop before this column was fitted
  # (engine 1.0.0-rc.1, SinkFlow::Stop).
  sink_stopped = list(
    bucket    = "raw_only",
    rationale = "Batch-engine response diagnostic; no R surface currently exposes ResponseColumnDiagnostic."
  ),
  # Adaptive grouping refit this column individually because the group
  # theta missed its probe objective by more than the refinement tolerance
  # (engine 1.0.0-rc.1).
  adaptive_refinement = list(
    bucket    = "raw_only",
    rationale = "Batch-engine response diagnostic; no R surface currently exposes ResponseColumnDiagnostic."
  )
)

mm_diagnostic_bucket <- function(code) {
  if (!length(code)) return(NA_character_)
  code <- as.character(code)
  if (!nzchar(code)) return(NA_character_)
  spec <- mm_diagnostic_code_registry[[code]] %||%
    mm_r_diagnostic_code_registry[[code]]
  if (is.null(spec)) return(NA_character_)
  spec$bucket
}

# Session-scoped state for the warn-once forward-compat gate.
mm_unknown_diag_state <- new.env(parent = emptyenv())
mm_unknown_diag_state$seen <- character()

mm_diagnostics_guard <- function(table) {
  codes <- table$code
  if (!length(codes)) return(table)
  codes <- unique(codes[nzchar(codes)])
  known <- c(names(mm_diagnostic_code_registry), names(mm_r_diagnostic_code_registry))
  unknown <- setdiff(codes, known)
  if (!length(unknown)) return(table)
  attr(table, "mm_unrecognized_diagnostic_code") <- unknown
  new_unknown <- setdiff(unknown, mm_unknown_diag_state$seen)
  if (length(new_unknown)) {
    mm_unknown_diag_state$seen <- c(mm_unknown_diag_state$seen, new_unknown)
    warning(
      sprintf(
        paste0(
          "mixeff: unrecognized DiagnosticCode(s) from the Rust engine: %s. ",
          "Advice surfaces will not format these codes; please update ",
          "mm_diagnostic_code_registry in R/diagnostics.R."
        ),
        paste(new_unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  table
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
