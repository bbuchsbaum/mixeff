#' Compile a mixed-effects model spec without fitting
#'
#' `compile_model()` parses the formula, runs the upstream semantic-IR /
#' design-audit pipeline against the supplied data, and returns an
#' `mm_spec` object — the audit-first analogue of the design-only step in
#' base `lm()`'s `model.frame()` / `model.matrix()` chain. Nothing is
#' optimized; nothing is fitted. `audit_design()`, `explain_model()`,
#' `random_options()`, and (in Phase 1.E) `lmm()` all consume the same
#' artifact.
#'
#' The compiled artifact is the structured truth: every print, summary,
#' and audit verb in mixeff reads back from it rather than re-deriving
#' meaning from formula text. R formats; Rust authors wording (PRD §9.6).
#'
#' Phase 1 compile scope: returns a populated `mm_spec` with the JSON
#' artifact attached. `explain_model()`, `random_options()`, and
#' `compare_covariance()` render random-effects guidance from upstream
#' random-term cards; the fit driver (`lmm()`) lands in 1.E.
#'
#' @param formula A two-sided lme4-style formula, e.g.
#'   `y ~ x + (1 + x | subject)`.
#' @param data A `data.frame` whose columns include every variable named
#'   in `formula`. Variables with missing values raise an
#'   `mm_data_error`; pass `na.omit(data)` explicitly if that is what you
#'   want.
#'
#' @return An object inheriting from `mm_spec` and containing:
#' \describe{
#'   \item{`call`}{the matched call}
#'   \item{`formula`}{the input formula}
#'   \item{`vars`}{character vector of variables read from `data`}
#'   \item{`model_frame`}{the data columns used to compile the artifact,
#'     retained so prefit audit views can evaluate nearby formula spellings}
#'   \item{`artifact`}{parsed JSON artifact (the
#'     `mixedmodels.compiled_model_artifact` v1 schema)}
#' }
#' The raw artifact JSON is attached as `attr(spec$artifact, "raw_json")`
#' so the post-compile FFI calls (e.g., the internal
#' `mm_audit_report_text` primitive) can round-trip without re-encoding.
#'
#' @section Errors:
#' Raises typed conditions (all inheriting from `mm_condition`):
#' \itemize{
#'   \item `mm_formula_error` — formula is not a two-sided R formula or
#'     fails parsing.
#'   \item `mm_data_error` — `data` is not a data.frame, refers to
#'     unknown variables, contains NAs in design columns, or has an
#'     unsupported column type.
#'   \item `mm_schema_error` — the artifact JSON returned by Rust does
#'     not match the wrapper's known schema set.
#' }
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   y       = rnorm(20),
#'   x       = rnorm(20),
#'   subject = factor(rep(letters[1:5], each = 4))
#' )
#' spec <- compile_model(y ~ x + (1 + x | subject), df)
#' audit_design(spec)
#' }
#'
#' @seealso [audit_design()] for the printed audit report.
#'
#' @export
compile_model <- function(formula, data) {
  call <- match.call()

  if (!inherits(formula, "formula")) {
    mm_abort(
      message = "`formula` must be a two-sided R formula (lhs ~ rhs).",
      class = "mm_formula_error",
      formula = formula
    )
  }
  if (length(formula) != 3L) {
    mm_abort(
      message = "`formula` must be two-sided (response on the left, predictors on the right).",
      class = "mm_formula_error",
      formula = formula
    )
  }
  if (!is.data.frame(data)) {
    mm_abort(
      message = "`data` must be a data.frame.",
      class = "mm_data_error",
      input = data
    )
  }

  vars <- all.vars(formula)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0L) {
    mm_abort(
      message = sprintf(
        "Variable(s) named in `formula` not found in `data`: %s.",
        paste(sprintf("`%s`", missing_vars), collapse = ", ")
      ),
      class = "mm_data_error",
      missing = missing_vars
    )
  }

  mm_check_no_na(data, vars)

  narrowed <- data[, vars, drop = FALSE]
  spec_data <- mm_translate_data(narrowed)
  formula_string <- mm_coerce_formula_string(formula)

  json <- tryCatch(
    .Call(
      wrap__mm_compile_model_json,
      formula_string,
      spec_data$column_order,
      spec_data$numeric_columns,
      spec_data$categorical_values,
      spec_data$categorical_levels,
      spec_data$categorical_ordered
    ),
    error = function(cnd) cnd
  )

  if (inherits(json, "condition")) {
    parts <- mm_split_tagged_error(conditionMessage(json))
    cls <- if (!is.na(parts$tag)) parts$tag else "mm_bridge_error"
    msg <- if (!is.na(parts$tag)) parts$message else conditionMessage(json)
    mm_abort(
      message = msg,
      class = cls,
      formula = formula_string,
      parent = json
    )
  }

  artifact <- mm_json_parse_artifact(json)

  spec <- list(
    call        = call,
    formula     = formula,
    vars        = vars,
    model_frame = narrowed,
    artifact    = artifact
  )
  class(spec) <- c("mm_spec", "mm_compiled")
  spec
}

#' @method print mm_spec
#' @export
print.mm_spec <- function(x, ...) {
  cat("<mm_spec>\n")
  cat(sprintf("  formula:           %s\n", deparse1(x$formula)))
  cat(sprintf("  effective formula: %s\n",
              x$artifact$effective_formula %||% x$artifact$requested_formula))
  fixed <- x$artifact$semantic_model$fixed_terms
  cat(sprintf("  fixed terms:       %s\n",
              if (length(fixed)) paste(unlist(fixed), collapse = ", ") else "(none)"))
  random <- x$artifact$semantic_model$random_terms
  cat(sprintf("  random terms:      %d\n", length(random)))
  cat(sprintf("  schema:            %s v%s\n",
              x$artifact$schema$schema_name,
              as.character(x$artifact$schema$schema_version)))
  cat("Use audit_design(spec) to view the structured design audit.\n")
  invisible(x)
}
