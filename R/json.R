#' Parse a versioned JSON artifact emitted by the Rust bridge
#'
#' Internal: every artifact crossing the FFI carries a top-level `schema`
#' object — see [mm_json_negotiate()] for the contract. `mm_json_parse_artifact()`
#' is the typed parser used by Phase 1 verbs (`compile_model()`,
#' `audit_design()`, and later `lmm()`): it reads the JSON with
#' `jsonlite::fromJSON(simplifyVector = FALSE)` so deeply-nested mixed
#' structures stay as nested lists, validates the schema header against the
#' wrapper's known set, and returns the parsed list with the original JSON
#' string preserved as `attr(.,"raw_json")` so downstream code can hand the
#' bytes to other FFI calls (e.g., the internal `mm_audit_report_text`
#' primitive) without re-encoding.
#'
#' The artifact's `schema_version` is encoded as a JSON number upstream
#' (see `COMPILED_ARTIFACT_SCHEMA_VERSION: u32 = 1`); this helper coerces
#' it to a length-1 character before negotiating, so the wrapper-side
#' string contract holds.
#'
#' @keywords internal
#' @noRd
mm_json_parse_artifact <- function(json) {
  if (!is.character(json) || length(json) != 1L || is.na(json) || !nzchar(json)) {
    mm_abort(
      message = "`json` must be a single non-empty character string.",
      class = "mm_schema_error",
      input = json
    )
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse artifact JSON: %s", conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )

  schema <- parsed$schema
  if (!is.list(schema)) {
    mm_abort(
      message = "Artifact JSON is missing a `schema` header.",
      class = "mm_schema_error",
      input = parsed
    )
  }
  header <- list(
    schema_name    = as.character(schema$schema_name),
    schema_version = as.character(schema$schema_version)
  )
  mm_json_negotiate(header)

  attr(parsed, "raw_json") <- json
  parsed
}

#' Parse a structured audit report emitted by the Rust bridge
#'
#' Internal counterpart to `mm_json_parse_artifact()`. Audit reports carry
#' their schema fields at top level (`schema_name`, `schema_version`), and
#' may contain nested `random_term_cards`, each with its own schema header.
#'
#' @keywords internal
#' @noRd
mm_json_parse_audit_report <- function(json) {
  if (!is.character(json) || length(json) != 1L || is.na(json) || !nzchar(json)) {
    mm_abort(
      message = "`json` must be a single non-empty character string.",
      class = "mm_schema_error",
      input = json
    )
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse audit report JSON: %s", conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )

  mm_json_negotiate(list(
    schema_name    = as.character(parsed$schema_name),
    schema_version = as.character(parsed$schema_version)
  ))

  cards <- parsed$random_term_cards %||% list()
  for (card in cards) {
    mm_json_negotiate(list(
      schema_name    = as.character(card$schema_name),
      schema_version = as.character(card$schema_version)
    ))
  }

  attr(parsed, "raw_json") <- json
  parsed
}

#' Parse a Rust fixed-effect inference table
#'
#' @keywords internal
#' @noRd
mm_json_parse_fixed_effect_inference_table <- function(table) {
  if (is.null(table)) {
    return(NULL)
  }
  if (!is.list(table)) {
    mm_abort(
      message = "`fixed_effect_inference_table` must be a list.",
      class = "mm_schema_error",
      input = table
    )
  }

  mm_json_negotiate(list(
    schema_name    = as.character(table$schema_name),
    schema_version = as.character(table$schema_version)
  ))

  rows <- table$rows %||% list()
  parsed_rows <- lapply(rows, mm_fixed_effect_inference_row)
  data <- if (length(parsed_rows)) {
    out <- do.call(rbind, parsed_rows)
    rownames(out) <- NULL
    out
  } else {
    mm_fixed_effect_inference_empty_table()
  }

  obj <- list(
    table = data,
    raw = table,
    schema_name = as.character(table$schema_name),
    schema_version = as.character(table$schema_version),
    crate_version = as.character(table$crate_version %||% NA_character_)
  )
  class(obj) <- "mm_fixed_effect_inference_table"
  obj
}

mm_fixed_effect_inference_row <- function(row) {
  label <- mm_scalar_text(row$label)
  statistic_name <- mm_optional_text(row$statistic_name)
  data.frame(
    term = label,
    label = label,
    kind = mm_scalar_text(row$kind, "coefficient"),
    estimate = mm_optional_numeric(row$estimate),
    std_error = mm_optional_numeric(row$std_error),
    df = mm_optional_numeric(row$denominator_df),
    numerator_df = mm_optional_numeric(row$numerator_df),
    denominator_df = mm_optional_numeric(row$denominator_df),
    statistic = mm_optional_numeric(row$statistic),
    statistic_name = statistic_name,
    p_value = mm_optional_numeric(row$p_value),
    method = mm_scalar_text(row$method, "not_computed"),
    status = mm_scalar_text(row$status, "not_assessed"),
    reliability = mm_scalar_text(row$reliability, "not_available"),
    # No default here: the engine's closed-enum warrant is passed through
    # verbatim, and an absent field stays NA rather than being dressed up
    # as a reason string the engine never authored.
    reliability_reason = mm_optional_text(row$reliability_reason),
    reason = mm_optional_text(row$reason),
    reason_code = mm_optional_text(row$reason_code),
    reason_detail = mm_optional_text(row$reason_detail),
    estimability = I(list(row$estimability %||% NULL)),
    details = I(list(row$details %||% NULL)),
    notes = I(list(as.character(unlist(row$notes %||% list(), use.names = FALSE)))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mm_fixed_effect_inference_empty_table <- function() {
  data.frame(
    term = character(),
    label = character(),
    kind = character(),
    estimate = numeric(),
    std_error = numeric(),
    df = numeric(),
    numerator_df = numeric(),
    denominator_df = numeric(),
    statistic = numeric(),
    statistic_name = character(),
    p_value = numeric(),
    method = character(),
    status = character(),
    reliability = character(),
    reliability_reason = character(),
    reason = character(),
    reason_code = character(),
    reason_detail = character(),
    estimability = I(list()),
    details = I(list()),
    notes = I(list()),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mm_optional_numeric <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_real_)
  }
  out <- suppressWarnings(as.numeric(x[[1L]]))
  if (length(out) != 1L || is.na(out)) NA_real_ else out
}

mm_optional_text <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_character_)
  }
  out <- as.character(x[[1L]])
  if (length(out) != 1L || is.na(out)) NA_character_ else out
}
