# Coverage gate for Rust diagnostic enums versus R-side advice surfaces.
#
# Bead bd-01KS31H5CMZEN75QHN8G52GCJH: every emitted DiagnosticCode and
# ResponseDiagnosticReason variant must have a documented R-side outcome
# (formatter bucket or raw_only). This file pairs the upstream enum
# source with R/diagnostics.R's mm_diagnostic_code_registry to detect
# silent drift.

mm_extract_rust_enum <- function(path, enum_name) {
  if (!file.exists(path)) {
    testthat::skip(sprintf("vendored Rust source missing: %s", path))
  }
  text <- readLines(path, warn = FALSE)
  start <- grep(sprintf("pub enum %s\\b", enum_name), text)
  if (!length(start)) {
    testthat::skip(sprintf("could not locate enum %s in %s", enum_name, path))
  }
  start <- start[[1L]]
  open_brace <- start + which(grepl("\\{", text[seq(start, length(text))]))[[1L]] - 1L
  remainder <- text[seq(open_brace + 1L, length(text))]
  close_offset <- which(grepl("^\\}", remainder))[[1L]]
  body <- remainder[seq_len(close_offset - 1L)]
  body <- body[!grepl("^\\s*(//|$)", body)]
  variants <- regmatches(body, regexpr("[A-Z][A-Za-z0-9]+", body))
  variants <- variants[nzchar(variants)]
  unique(variants)
}

mm_camel_to_snake <- function(x) {
  x <- gsub("([A-Z])", "_\\1", x)
  x <- sub("^_", "", x)
  tolower(x)
}

mm_diagnostic_enum_path <- function() {
  file.path("..", "..", "src", "rust", "upstream", "mixeff-rs", "src",
            "compiler", "diagnostics.rs")
}

mm_batch_enum_path <- function() {
  file.path("..", "..", "src", "rust", "upstream", "mixeff-rs", "src",
            "model", "batch.rs")
}

test_that("mm_diagnostic_code_registry covers every Rust DiagnosticCode variant", {
  variants <- mm_extract_rust_enum(mm_diagnostic_enum_path(), "DiagnosticCode")
  snake <- mm_camel_to_snake(variants)
  registered <- names(mixeff:::mm_diagnostic_code_registry)
  missing <- setdiff(snake, registered)
  stale   <- setdiff(registered, snake)
  expect_identical(
    missing, character(),
    info = sprintf(
      "mm_diagnostic_code_registry is missing: %s",
      paste(missing, collapse = ", ")
    )
  )
  expect_identical(
    stale, character(),
    info = sprintf(
      "mm_diagnostic_code_registry has codes the Rust enum no longer emits: %s",
      paste(stale, collapse = ", ")
    )
  )
})

test_that("mm_response_diagnostic_reason_registry covers every Rust ResponseDiagnosticReason variant", {
  variants <- mm_extract_rust_enum(mm_batch_enum_path(),
                                   "ResponseDiagnosticReason")
  snake <- mm_camel_to_snake(variants)
  registered <- names(mixeff:::mm_response_diagnostic_reason_registry)
  missing <- setdiff(snake, registered)
  stale   <- setdiff(registered, snake)
  expect_identical(
    missing, character(),
    info = sprintf(
      "mm_response_diagnostic_reason_registry is missing: %s",
      paste(missing, collapse = ", ")
    )
  )
  expect_identical(
    stale, character(),
    info = sprintf(
      "mm_response_diagnostic_reason_registry has codes the Rust enum no longer emits: %s",
      paste(stale, collapse = ", ")
    )
  )
})

test_that("every registered diagnostic code is bound to a valid bucket", {
  buckets <- vapply(mixeff:::mm_diagnostic_code_registry,
                    function(spec) spec$bucket, character(1))
  expect_true(all(buckets %in% c("design_note", "repair", "fit_note", "raw_only")),
              info = sprintf("invalid buckets: %s",
                             paste(setdiff(unique(buckets),
                                           c("design_note", "repair",
                                             "fit_note", "raw_only")),
                                   collapse = ", ")))
  raw_only <- buckets == "raw_only"
  if (any(raw_only)) {
    rationales <- vapply(mixeff:::mm_diagnostic_code_registry[raw_only],
                         function(spec) spec$rationale %||% "", character(1))
    expect_true(all(nzchar(rationales)),
                info = sprintf(
                  "raw_only entries missing rationale: %s",
                  paste(names(rationales)[!nzchar(rationales)], collapse = ", ")
                ))
  }
})

test_that("each formattable bucket has a non-empty advice-line renderer", {
  diag <- function(code, message, severity = "info", stage = "design_audit",
                   suggested_actions = list()) {
    list(code = code, severity = severity, stage = stage, message = message,
         affected_terms = list(), suggested_actions = suggested_actions)
  }
  ds <- list(
    diag("scope_note", "scope note message"),
    diag("boundary_parameter", "theta on boundary", stage = "optimization"),
    diag("structural_refusal", "refusal message",
         suggested_actions = list("rewrite the term"))
  )
  expect_true(length(mixeff:::mm_design_note_lines(ds)) >= 1L)
  expect_true(length(mixeff:::mm_fit_note_lines(ds)) >= 1L)
  expect_true(length(mixeff:::mm_repair_lines(ds)) >= 1L)
})

test_that("mm_repair_lines falls back to message when suggested_actions is empty", {
  ds <- list(list(code = "not_identifiable", severity = "error",
                  stage = "estimability", message = "MLE does not exist",
                  affected_terms = list(),
                  suggested_actions = list()))
  lines <- mixeff:::mm_repair_lines(ds)
  expect_length(lines, 1L)
  expect_match(lines, "not_identifiable: MLE does not exist", fixed = TRUE)
})

test_that("mm_diagnostics_table warns once and attaches attribute on unknown code", {
  # Reset the warn-once state so the test is order-independent.
  state <- mixeff:::mm_unknown_diag_state
  state$seen <- character()
  on.exit({
    state$seen <- character()
  }, add = TRUE)

  fake <- list(
    list(code = "scope_note", severity = "info", stage = "design_audit",
         message = "known", affected_terms = list()),
    list(code = "fictional_future_code", severity = "warning",
         stage = "optimization", message = "from a future Rust release",
         affected_terms = list())
  )
  table <- withCallingHandlers(
    mixeff:::mm_diagnostics_table(fake),
    warning = function(w) {
      expect_match(conditionMessage(w), "unrecognized DiagnosticCode")
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(attr(table, "mm_unrecognized_diagnostic_code"),
                   "fictional_future_code")

  # Second call within the same session must not re-warn for the same code.
  table2 <- withCallingHandlers(
    mixeff:::mm_diagnostics_table(fake),
    warning = function(w) {
      fail(sprintf("re-warned within session for known unknown: %s",
                   conditionMessage(w)))
    }
  )
  expect_identical(attr(table2, "mm_unrecognized_diagnostic_code"),
                   "fictional_future_code")

  # Unknown codes are not classified into any advice bucket.
  bucket <- mixeff:::mm_diagnostic_bucket("fictional_future_code")
  expect_true(is.na(bucket))
})

test_that("formatted code keeps mm_diagnostic_bucket aligned with the formatter family", {
  reg <- mixeff:::mm_diagnostic_code_registry
  for (code in names(reg)) {
    expect_identical(
      mixeff:::mm_diagnostic_bucket(code), reg[[code]]$bucket,
      info = sprintf("registry/bucket helper drift for `%s`", code)
    )
  }
})
