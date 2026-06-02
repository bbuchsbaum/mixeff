test_that("lme4 parity mismatch ledger is shipped and well-formed", {
  path <- system.file("extdata", "expected-mismatches.json", package = "mixeff")
  expect_true(file.exists(path))
  ledger <- jsonlite::read_json(path, simplifyVector = FALSE)

  expect_identical(ledger$schema_version, 1L)
  expect_true(is.list(ledger$mismatches))
  required <- c("case_id", "field", "status", "reason")
  for (entry in ledger$mismatches) {
    for (key in required) {
      expect_true(key %in% names(entry),
                  info = sprintf("ledger entry missing `%s` key", key))
    }
    expect_true(entry$status %in% c("expected_mismatch", "upstream_bug",
                                    "unsupported",
                                    "design_weak_identifiability", "pass"),
                info = sprintf("ledger entry has unrecognized status `%s` for case `%s` field `%s`",
                               entry$status, entry$case_id, entry$field))
    if (entry$status %in% c("expected_mismatch", "upstream_bug")) {
      has_bound <- !is.null(entry$expected_max_abs_diff) ||
        !is.null(entry$expected_max_rel_diff)
      expect_true(has_bound,
                  info = sprintf("ledger entry for case `%s` field `%s` has status `%s` but no abs/rel bound",
                                 entry$case_id, entry$field, entry$status))
    }
  }

  keys <- vapply(ledger$mismatches, function(e) paste(e$case_id, e$field, sep = "::"),
                 character(1))
  expect_equal(length(unique(keys)), length(keys),
               info = "ledger contains duplicate (case_id, field) entries")
})

test_that("upstream mismatch report template is shipped", {
  path <- system.file("extdata", "upstream-mismatch-report-template.md",
                      package = "mixeff")
  expect_true(file.exists(path))
  text <- readLines(path, warn = FALSE)
  required <- c(
    "Dataset or generator",
    "Generator seed",
    "Formula",
    "REML",
    "lme4",
    "lmerTest",
    "pbkrtest",
    "mixeff commit",
    "mixedmodels commit",
    "Compared Fields",
    "tolerance",
    "reference value",
    "observed Rust value",
    "Minimal Reproducer"
  )
  for (pattern in required) {
    expect_true(any(grepl(pattern, text, fixed = TRUE)),
                info = sprintf("Template is missing `%s`", pattern))
  }
})

test_that("classic lme4 parity fixture manifest is valid", {
  manifest <- mm_lme4_parity_manifest()

  expect_identical(manifest$schema_version, 1L)
  expect_true(length(manifest$cases) >= 6L)
  ids <- vapply(manifest$cases, `[[`, character(1), "id")
  expect_equal(length(unique(ids)), length(ids))
  expect_true(all(vapply(manifest$cases, function(case) {
    all(c("dataset", "formula", "reml", "expected_status") %in% names(case))
  }, logical(1))))
})

test_that("classic lme4 LMM cases match core extractors within documented tolerances", {
  mm_skip_if_no_lme4()

  for (case in mm_lme4_parity_cases()) {
    mm_expect_core_lme4_parity(case)
  }
})

test_that("classic lme4 LMM cases match random-effect modes where labels align", {
  mm_skip_if_no_lme4()

  for (case in mm_lme4_parity_cases()) {
    mm_expect_ranef_lme4_parity(case)
  }
})

test_that("classic lme4 LMM cases match supported prediction semantics", {
  mm_skip_if_no_lme4()

  for (case in mm_lme4_parity_cases()) {
    mm_expect_prediction_lme4_parity(case)
  }
})
