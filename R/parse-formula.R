#' Parse and canonicalize an lme4-style formula
#'
#' `mm_parse_formula()` parses a formula string through the Rust formula
#' parser and returns its canonical `Display` rendering. This is the Phase 0
#' round-trip primitive: equivalent formula spellings produce identical
#' canonical strings, so equivalence-class testing in R is just string
#' comparison on the canonical form.
#'
#' @param formula A single character string (length 1, non-NA, non-empty), or
#'   a one-sided / two-sided R `formula` object. R `formula` objects are
#'   coerced to character via `format()` before parsing.
#'
#' @return A single character string: the canonical rendering of the parsed
#'   formula.
#'
#' @section Errors:
#' Parse failures are signalled as a typed `mm_formula_error` condition (also
#' inheriting from `mm_condition` and `error`). The condition object carries
#' the original input string in its `formula` field. Catch with
#' `tryCatch(..., mm_formula_error = handler)`.
#'
#' @examples
#' mm_parse_formula("y ~ x + (1 | g)")
#' mm_parse_formula(y ~ x + (1 | g))
#'
#' @export
mm_parse_formula <- function(formula) {
  input <- mm_coerce_formula_string(formula)

  result <- tryCatch(
    .Call(wrap__mm_parse_formula, input),
    error = function(cnd) cnd
  )

  if (inherits(result, "condition")) {
    parts <- mm_split_tagged_error(conditionMessage(result))
    if (identical(parts$tag, "mm_formula_error")) {
      mm_abort(
        message = parts$message,
        class = "mm_formula_error",
        formula = input
      )
    }
    # Untagged error from the bridge — re-throw as a generic mm_condition so
    # callers can still catch it. This should be rare; an untagged bridge
    # error indicates either an extendr-level fault or a missing tag in
    # `src/rust/src/lib.rs`.
    mm_abort(
      message = conditionMessage(result),
      class = "mm_bridge_error",
      formula = input
    )
  }

  result
}

# Coerce supported formula inputs to a single character string suitable for
# the Rust parser. Accepts a length-1 character or an R `formula`.
mm_coerce_formula_string <- function(formula) {
  if (inherits(formula, "formula")) {
    s <- paste(trimws(format(formula)), collapse = " ")
    return(s)
  }
  if (!is.character(formula) || length(formula) != 1L || is.na(formula)) {
    mm_abort(
      message = "`formula` must be a single non-NA character string or an R formula object.",
      class = "mm_formula_error",
      formula = formula
    )
  }
  if (!nzchar(trimws(formula))) {
    mm_abort(
      message = "`formula` must be non-empty.",
      class = "mm_formula_error",
      formula = formula
    )
  }
  formula
}
