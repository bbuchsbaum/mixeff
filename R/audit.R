#' Print the structured design audit for a compiled model spec
#'
#' `audit_design()` returns the user-facing audit report attached to an
#' `mm_spec`. The text is rendered by the upstream Rust crate (the
#' `mixedmodels.model_audit_report` schema's `Display` impl) — Rust authors
#' the wording, R formats nothing. Routing every printed audit line
#' through the upstream renderer is what enforces the R9 "no advice
#' creep" contract: drift in scope notes / tone is visible in one place
#' rather than scattered across R formatters.
#'
#' Phase 1.A scope: `audit_design()` accepts an `mm_spec` from
#' [compile_model()] and emits the report sections (Requested Model,
#' Model State, Fixed/Random Effects, Information Budget, Dependence
#' Paths, Parameterization Trace, Effective Covariance, Policy
#' Recommendations, Optimizer, Inference, Diagnostics). Sections that
#' depend on a fit (Optimizer / Inference) report `not applicable
#' before fitting` on a pre-fit spec.
#'
#' @param spec An `mm_spec` produced by [compile_model()] (Phase 1.A) or
#'   an `mm_fit` (post-Phase-1.E).
#'
#' @return An object of class `mm_audit` carrying:
#' \describe{
#'   \item{`text`}{the rendered report text (a single character string,
#'     newline-separated)}
#'   \item{`summary_text`}{the compact report rendered by the upstream
#'     `ModelAuditReport::render_summary` (Audit Summary plus the
#'     Requested Model section)}
#'   \item{`design_audit`}{the parsed `design_audit` field from the
#'     `CompiledModelArtifact` (random-term audits, fixed-effect rank,
#'     covariance kernel graph, ...) — `NULL` on uncompilable formulas}
#'   \item{`report`}{the parsed upstream `ModelAuditReport` v2, including
#'     Rust-authored `random_term_cards` for downstream explanation verbs}
#'   \item{`random_term_cards`}{the report's per-random-term cards, copied
#'     to the top level for convenient inspection}
#'   \item{`cross_card_constraints`}{report-level constraints between
#'     random-term cards}
#'   \item{`diagnostics`}{the parsed report diagnostics, falling back to
#'     artifact diagnostics when needed}
#' }
#' `print.mm_audit` defaults to the compact upstream-rendered summary in
#' `summary_text`. Use `print(x, full = TRUE)` for the complete upstream
#' report stored in `text`.
#'
#' @section Errors:
#' Raises an `mm_schema_error` if the supplied object does not carry a
#' parsed artifact with the expected schema header.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   y       = rnorm(20),
#'   x       = rnorm(20),
#'   subject = factor(rep(letters[1:5], each = 4))
#' )
#' audit_design(compile_model(y ~ x + (1 + x | subject), df))
#' }
#'
#' @seealso [compile_model()].
#'
#' @export
audit_design <- function(spec) {
  if (!inherits(spec, c("mm_spec", "mm_fit"))) {
    mm_abort(
      message = "`spec` must be an `mm_spec` (from compile_model) or an `mm_fit` (from lmm).",
      class = "mm_schema_error",
      input = spec
    )
  }
  artifact <- spec$artifact
  raw_json <- attr(artifact, "raw_json")
  if (is.null(raw_json) || !nzchar(raw_json)) {
    mm_abort(
      message = "Artifact is missing its raw JSON; cannot delegate to the upstream audit report renderer.",
      class = "mm_schema_error",
      input = artifact
    )
  }

  text <- mm_audit_bridge_call(mm_audit_report_text, raw_json)
  summary_text <- mm_audit_bridge_call(mm_audit_report_summary_text, raw_json)
  report_json <- mm_audit_bridge_call(mm_audit_report_json, raw_json)
  report <- mm_json_parse_audit_report(report_json)
  supplemental <- mm_r_design_diagnostics(spec)
  text <- mm_append_r_design_diagnostic_text(text, supplemental)
  summary_text <- mm_append_r_design_diagnostic_text(summary_text, supplemental)

  out <- list(
    text                   = text,
    summary_text           = summary_text,
    design_audit           = artifact$design_audit,
    report                 = report,
    random_term_cards      = report$random_term_cards %||% list(),
    cross_card_constraints = report$cross_card_constraints %||% list(),
    diagnostics            = c(report$diagnostics %||% artifact$diagnostics %||% list(),
                               supplemental)
  )
  class(out) <- "mm_audit"
  out
}

#' Audit a fitted mixeff model
#'
#' `audit()` is the post-fit alias for [audit_design()]. It renders the same
#' upstream-authored audit report, now backed by the fitted artifact carried by
#' an `mm_fit`.
#'
#' @param fit A fitted `mm_fit`, usually from [lmm()].
#' @param ... Reserved for future methods.
#'
#' @return An `mm_audit` object; see [audit_design()].
#'
#' @export
audit <- function(fit, ...) {
  UseMethod("audit")
}

#' @rdname audit
#' @export
audit.mm_fit <- function(fit, ...) {
  audit_design(fit)
}

#' @method print mm_audit
#' @export
print.mm_audit <- function(x, full = FALSE, ...) {
  text <- if (isTRUE(full)) x$text else x$summary_text %||% x$text
  cat(text)
  if (!grepl("\n$", text)) cat("\n")
  invisible(x)
}

mm_audit_bridge_call <- function(fn, raw_json) {
  out <- tryCatch(fn(raw_json), error = function(cnd) cnd)
  if (inherits(out, "condition")) {
    parts <- mm_split_tagged_error(conditionMessage(out))
    cls <- if (!is.na(parts$tag)) parts$tag else "mm_bridge_error"
    msg <- if (!is.na(parts$tag)) parts$message else conditionMessage(out)
    mm_abort(
      message = msg,
      class = cls,
      parent = out
    )
  }
  out
}
