# Test helpers wiring the parity ledger
# (`inst/extdata/expected-mismatches.json`) into per-field assertions.
#
# Contract:
#   - No ledger entry  -> strict expect_equal at the case tolerance.
#                         Failure means an unclassified parity regression.
#   - status pass      -> same as no entry (entry exists only for documentation).
#   - status expected_mismatch / upstream_bug
#                      -> assert observed_max_abs_diff <= entry$expected_max_abs_diff
#                         and observed_max_rel_diff <= entry$expected_max_rel_diff,
#                         then succeed() with the documented reason.
#                         Tightening (smaller diff) is good; growing past the
#                         recorded bound trips the test.
#   - status unsupported
#                      -> record but do not assert. The Rust contract does not
#                         certify this field for this case.
#
# Every call records a row in mm_scoreboard_state for the parity scoreboard
# output emitted by tests/testthat/test-parity-scoreboard.R.

mm_scoreboard_state <- new.env(parent = emptyenv())
mm_scoreboard_state$rows <- list()

mm_scoreboard_reset <- function() {
  mm_scoreboard_state$rows <- list()
}

mm_scoreboard_record <- function(case_id, field,
                                 observed_max_abs_diff,
                                 observed_max_rel_diff,
                                 case_tolerance,
                                 status,
                                 reason = NA_character_) {
  mm_scoreboard_state$rows <- c(
    mm_scoreboard_state$rows,
    list(list(
      case_id = case_id,
      field = field,
      case_tolerance = case_tolerance,
      observed_max_abs_diff = observed_max_abs_diff,
      observed_max_rel_diff = observed_max_rel_diff,
      status = status,
      reason = reason
    ))
  )
}

mm_scoreboard_table <- function() {
  rows <- mm_scoreboard_state$rows
  if (!length(rows)) {
    return(data.frame(
      case_id = character(),
      field = character(),
      case_tolerance = numeric(),
      observed_max_abs_diff = numeric(),
      observed_max_rel_diff = numeric(),
      status = character(),
      reason = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

mm_scoreboard_classify_status <- function(status) {
  if (identical(status, "unsupported")) {
    return("unsupported")
  }
  if (identical(status, "expected_mismatch") || identical(status, "upstream_bug")) {
    return("xfail")
  }
  "pass"
}

mm_scoreboard_compact <- function(cases, rows = mm_scoreboard_table()) {
  rows <- as.data.frame(rows, stringsAsFactors = FALSE)
  if (!nrow(rows)) {
    return(data.frame(
      case_id = character(),
      dataset = character(),
      formula = character(),
      reml = logical(),
      reference_versions = character(),
      compared_field = character(),
      tolerance = numeric(),
      observed_difference = numeric(),
      status = character(),
      reason = character(),
      stringsAsFactors = FALSE
    ))
  }

  case_map <- setNames(cases, vapply(cases, `[[`, character(1), "id"))
  refs <- mm_reference_versions()
  ref_label <- paste(names(refs), refs, sep = "=", collapse = ", ")

  case_lookup <- function(id) {
    entry <- case_map[[id]]
    if (is.null(entry)) return(list(dataset = NA_character_, formula = NA_character_, reml = NA))
    entry
  }

  data.frame(
    case_id = rows$case_id,
    dataset = vapply(rows$case_id, function(id) case_lookup(id)$dataset, character(1)),
    formula = vapply(rows$case_id, function(id) case_lookup(id)$formula, character(1)),
    reml = as.logical(vapply(rows$case_id, function(id) isTRUE(case_lookup(id)$reml), logical(1))),
    reference_versions = rep(ref_label, nrow(rows)),
    compared_field = rows$field,
    tolerance = rows$case_tolerance,
    observed_difference = ifelse(
      is.na(rows$observed_max_abs_diff) & is.na(rows$observed_max_rel_diff),
      NA_real_,
      pmax(rows$observed_max_abs_diff, rows$observed_max_rel_diff, na.rm = TRUE)
    ),
    status = vapply(rows$status, mm_scoreboard_classify_status, character(1)),
    reason = ifelse(is.na(rows$reason), "", rows$reason)
  )
}

mm_scoreboard_emit_json <- function(output_path,
                                   cases = NULL,
                                   rows = mm_scoreboard_table()) {
  if (is.null(cases)) {
    stop("cases is required to emit a compact R-facing parity scoreboard", call. = FALSE)
  }
  compact <- mm_scoreboard_compact(cases, rows)
  payload <- list(
    schema_version = 1L,
    schema_name = "mixeff.lme4_parity_scoreboard",
    generated_at = as.character(Sys.time()),
    reference_versions = mm_reference_versions(),
    rows = compact
  )
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, output_path, auto_unbox = TRUE, pretty = TRUE)
  invisible(compact)
}

mm_parity_diffs <- function(observed, expected) {
  observed_n <- mm_numeric_payload(observed)
  expected_n <- mm_numeric_payload(expected)
  if (length(observed_n) != length(expected_n)) {
    stop(sprintf(
      "Length mismatch comparing parity vectors: observed=%d expected=%d",
      length(observed_n), length(expected_n)
    ), call. = FALSE)
  }
  diffs <- abs(observed_n - expected_n)
  if (!length(diffs)) {
    return(list(abs = 0, rel = 0))
  }
  abs_max <- max(diffs, na.rm = TRUE)
  denom <- pmax(abs(expected_n), .Machine$double.eps)
  rel_max <- max(diffs / denom, na.rm = TRUE)
  list(abs = abs_max, rel = rel_max)
}

mm_assert_parity <- function(observed, expected, case_id, field,
                             tolerance, label,
                             mode = c("relative", "absolute")) {
  mode <- match.arg(mode)
  d <- mm_parity_diffs(observed, expected)
  entry <- tryCatch(mm_parity_lookup(case_id, field),
                    error = function(cnd) NULL)

  status <- entry$status %||% "pass"
  reason <- entry$reason %||% NA_character_
  mm_scoreboard_record(case_id, field, d$abs, d$rel, tolerance, status, reason)

  refs <- mm_reference_versions()
  ref_label <- paste(names(refs), refs, sep = "=", collapse = ", ")

  if (is.null(entry) || identical(status, "pass")) {
    if (identical(mode, "absolute")) {
      testthat::expect_lte(
        d$abs, tolerance,
        label = sprintf(
          "%s parity for case `%s` field `%s`: observed_max_abs_diff %.3e exceeds case absolute tolerance %.3e (%s)",
          label, case_id, field, d$abs, tolerance, ref_label
        )
      )
    } else {
      testthat::expect_equal(
        mm_numeric_payload(observed),
        mm_numeric_payload(expected),
        tolerance = tolerance,
        info = sprintf(
          "%s parity failed for case `%s` field `%s`; tolerance=%s; reference versions: %s",
          label, case_id, field,
          format(tolerance, scientific = TRUE), ref_label
        )
      )
    }
    return(invisible())
  }

  if (identical(status, "unsupported")) {
    testthat::succeed(message = sprintf(
      "documented `unsupported` for case `%s` field `%s`: %s",
      case_id, field, reason
    ))
    return(invisible())
  }

  bound_abs <- entry$expected_max_abs_diff
  bound_rel <- entry$expected_max_rel_diff
  if (!is.null(bound_abs)) {
    testthat::expect_lte(
      d$abs, bound_abs,
      label = sprintf(
        "documented `%s` for case `%s` field `%s`: observed_max_abs_diff %.3e exceeds recorded bound %.3e (%s); ledger reason: %s",
        status, case_id, field, d$abs, bound_abs, ref_label, reason
      )
    )
  }
  if (!is.null(bound_rel)) {
    testthat::expect_lte(
      d$rel, bound_rel,
      label = sprintf(
        "documented `%s` for case `%s` field `%s`: observed_max_rel_diff %.3e exceeds recorded bound %.3e (%s); ledger reason: %s",
        status, case_id, field, d$rel, bound_rel, ref_label, reason
      )
    )
  }
  if (is.null(bound_abs) && is.null(bound_rel)) {
    testthat::succeed(message = sprintf(
      "documented `%s` for case `%s` field `%s` (no numeric bound): %s",
      status, case_id, field, reason
    ))
  }
  invisible()
}
