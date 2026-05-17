#' Negotiate a JSON schema header against what `mixeff` supports
#'
#' Every artifact crossing the Rust-R bridge carries a header that names a
#' schema (e.g., `formula`, `artifact`, `audit`) and the version of that
#' schema (`v0`, `v1`, ...). `mm_json_negotiate()` validates a header
#' against the closed set of `(schema_name, schema_version)` pairs the
#' current wrapper build understands, and raises a typed
#' `mm_schema_error` on mismatch.
#'
#' This is the *fast-fail* primitive: any code path that consumes a Rust
#' artifact should call `mm_json_negotiate()` before parsing the body, so
#' a version skew between the Rust crate and the R wrapper produces a
#' single clean error rather than a confusing field-by-field decode failure.
#'
#' @param header A list with at least `schema_name` and `schema_version`
#'   as length-1 character strings. Additional fields (e.g.,
#'   `crate_version`, `package_version`) are accepted and ignored at this
#'   layer; downstream code records them on the `mm_fit` provenance.
#'
#' @return Invisibly returns `TRUE` on a clean match.
#'
#' @section Errors:
#' Any of the following raise an `mm_schema_error` (also inheriting from
#' `mm_condition` and `error`):
#' \itemize{
#'   \item `header` is not a list, or is missing `schema_name` /
#'     `schema_version`, or those fields are not length-1 character.
#'   \item `schema_name` is not in the wrapper's known set
#'     (see [mm_json_known_schemas()]).
#'   \item `schema_name` is known but `schema_version` does not match what
#'     the wrapper expects.
#' }
#' The condition object carries the original `header` in its `input` field.
#' (The field is *not* called `header` because rlang reserves that name on
#' conditions for `cnd_header()` formatting.)
#'
#' @examples
#' mm_json_negotiate(list(schema_name = "formula", schema_version = "v0"))
#' \dontrun{
#' # Raises mm_schema_error:
#' mm_json_negotiate(list(schema_name = "formula", schema_version = "v99"))
#' }
#'
#' @seealso [mm_json_known_schemas()] for the closed set,
#'   [mm_formula_manifest()] for the broader capability record.
#'
#' @export
mm_json_negotiate <- function(header) {
  if (!is.list(header)) {
    mm_abort(
      message = "`header` must be a list with `schema_name` and `schema_version`.",
      class = "mm_schema_error",
      input = header
    )
  }
  for (req in c("schema_name", "schema_version")) {
    val <- header[[req]]
    if (is.null(val) || !is.character(val) || length(val) != 1L || is.na(val)) {
      mm_abort(
        message = sprintf("`header$%s` must be a single non-NA character string.", req),
        class = "mm_schema_error",
        input = header
      )
    }
  }

  result <- tryCatch(
    .Call(wrap__mm_json_negotiate_one,
          header$schema_name,
          header$schema_version),
    error = function(cnd) cnd
  )

  if (inherits(result, "condition")) {
    parts <- mm_split_tagged_error(conditionMessage(result))
    msg <- if (identical(parts$tag, "mm_schema_error")) {
      parts$message
    } else {
      conditionMessage(result)
    }
    mm_abort(
      message = msg,
      class = "mm_schema_error",
      input = header
    )
  }

  invisible(TRUE)
}

#' Closed list of schema/version pairs the wrapper understands
#'
#' Returns the known-schema table that backs [mm_json_negotiate()]. New
#' schemas appear here as later phases add artifacts (compiled-model,
#' audit, theta_map, certificate, inference, reproducibility, prediction).
#'
#' @return A data frame with two character columns: `name` and `version`.
#'
#' @examples
#' mm_json_known_schemas()
#'
#' @seealso [mm_json_negotiate()].
#'
#' @export
mm_json_known_schemas <- function() {
  raw <- .Call(wrap__mm_json_known_schemas)
  data.frame(
    name    = as.character(raw$name),
    version = as.character(raw$version),
    stringsAsFactors = FALSE
  )
}
