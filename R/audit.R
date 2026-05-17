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
#' depend on a fit (Optimizer / Inference) report `not assessed` until
#' Phase 1.E lands `lmm()`.
#'
#' @param spec An `mm_spec` produced by [compile_model()] (Phase 1.A) or
#'   an `mm_fit` (post-Phase-1.E).
#'
#' @return An object of class `mm_audit` carrying:
#' \describe{
#'   \item{`text`}{the rendered report text (a single character string,
#'     newline-separated)}
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
#' `print.mm_audit` calls `cat()` on `text`, so calling
#' `audit_design(spec)` at the prompt prints the upstream-rendered
#' report once.
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

  text <- tryCatch(
    .Call(wrap__mm_audit_report_text, raw_json),
    error = function(cnd) cnd
  )
  if (inherits(text, "condition")) {
    parts <- mm_split_tagged_error(conditionMessage(text))
    cls <- if (!is.na(parts$tag)) parts$tag else "mm_bridge_error"
    msg <- if (!is.na(parts$tag)) parts$message else conditionMessage(text)
    mm_abort(
      message = msg,
      class = cls,
      parent = text
    )
  }

  report_json <- tryCatch(
    mm_audit_report_json(raw_json),
    error = function(cnd) cnd
  )
  if (inherits(report_json, "condition")) {
    parts <- mm_split_tagged_error(conditionMessage(report_json))
    cls <- if (!is.na(parts$tag)) parts$tag else "mm_bridge_error"
    msg <- if (!is.na(parts$tag)) parts$message else conditionMessage(report_json)
    mm_abort(
      message = msg,
      class = cls,
      parent = report_json
    )
  }
  report <- mm_json_parse_audit_report(report_json)

  out <- list(
    text                   = text,
    design_audit           = artifact$design_audit,
    report                 = report,
    random_term_cards      = report$random_term_cards %||% list(),
    cross_card_constraints = report$cross_card_constraints %||% list(),
    diagnostics            = report$diagnostics %||% artifact$diagnostics %||% list()
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
print.mm_audit <- function(x, ...) {
  cat(x$text)
  if (!grepl("\n$", x$text)) cat("\n")
  invisible(x)
}
