test_that("classic lme4 parity checks emit a compact R-facing scoreboard artifact", {
  mm_skip_if_no_lme4()
  mm_scoreboard_reset()

  cases <- mm_lme4_parity_cases()
  for (case in cases) {
    mm_expect_core_lme4_parity(case)
    mm_expect_ranef_lme4_parity(case)
    mm_expect_prediction_lme4_parity(case)
  }

  rows <- mm_scoreboard_table()
  expect_true(nrow(rows) > 0L, info = "parity scoreboard captured at least one row")

  compact <- mm_scoreboard_compact(cases = cases, rows = rows)
  expect_true(nrow(compact) > 0L, info = "compact scoreboard has rows")
  expected <- c(
    "case_id", "dataset", "formula", "reml", "reference_versions",
    "compared_field", "tolerance", "observed_difference", "status", "reason"
  )
  expect_true(all(expected %in% names(compact)))

  output <- file.path(tempdir(), "mixeff-lme4-parity-scoreboard.json")
  mm_scoreboard_emit_json(output, cases = cases, rows = rows)
  expect_true(file.exists(output), info = sprintf("scoreboard artifact exists at `%s`", output))
})
