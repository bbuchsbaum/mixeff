#' Typed conditions raised by mixeff
#'
#' All structured errors and informational signals raised by the package
#' inherit from the `mm_condition` class so callers can catch them
#' generically (`tryCatch(..., mm_condition = handler)`), or by a more
#' specific subclass for finer-grained handling.
#'
#' The package ships `mm_formula_error`, `mm_schema_error`, and the
#' `mm_bridge_error` fallback for untagged Rust errors; `mm_data_error`
#' (raised by `compile_model()` on data shape / type problems);
#' `mm_fit_error` for fit construction / optimization failures; and
#' `mm_inference_unavailable` for inference, extractor, or prediction
#' requests the engine cannot certify on this fit. Further classes are
#' `mm_not_identifiable`, `mm_design_refusal`, and `mm_fit_not_optimized`
#' (see PRD §8.1).
#'
#' `mm_arg_error` is distinct from all of the above: it signals that the
#' *caller passed an invalid or malformed argument* (wrong type, value out
#' of range, unknown option name, missing required field), independent of
#' any model, engine, schema, or data semantics. Catching
#' `mm_inference_unavailable` must mean "the engine refuses to certify
#' this", not "you typed it wrong" — so plain argument validation raises
#' `mm_arg_error`, never the domain-refusal classes.
#'
#' @keywords internal
#' @name mm-conditions
NULL

# ---- internal helpers ------------------------------------------------------

# Reserved condition-field names — DO NOT use these as custom-data keys
# when raising typed conditions:
#
#   header, body, footer
#
# rlang's cnd_header() / cnd_body() / cnd_footer() formatters dispatch on
# these fields. If you set `header = <list>` on a condition (because you
# wanted to attach a parsed JSON header to an mm_schema_error, say),
# rlang's formatter will refuse with:
#
#     `header` field must be a character vector or a function.
#
# ...the moment any code path reads the condition message, including
# default print and many testthat expect_error() traces. We hit this in
# Phase 0 with mm_schema_error and renamed the field to `input`. Use
# `input`, `offending`, `value`, `subject`, etc. — anything but the three
# above. (See `R/schema.R::mm_json_negotiate()` for the canonical pattern.)

# Construct and signal a typed mixeff condition. All mixeff condition classes
# inherit from "mm_condition" (after the specific class) so callers can match
# either generically or specifically.
mm_abort <- function(message,
                     class,
                     ...,
                     call = rlang::caller_env(),
                     parent = NULL) {
  rlang::abort(
    message = message,
    class = c(class, "mm_condition"),
    ...,
    call = call,
    parent = parent
  )
}

# Strip the "mm_<name>: " prefix that the Rust bridge attaches to errors so
# they can be routed to the right typed condition on the R side. Returns a
# list with `tag` (character or NA) and `message` (the residual text).
mm_split_tagged_error <- function(message) {
  m <- regmatches(message, regexec("^(mm_[a-z_]+): (.*)$", message))[[1L]]
  if (length(m) == 3L) {
    list(tag = m[[2L]], message = m[[3L]])
  } else {
    list(tag = NA_character_, message = message)
  }
}
