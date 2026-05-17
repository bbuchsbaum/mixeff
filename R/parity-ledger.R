# Internal helpers for the lme4 parity ledger
# (`inst/extdata/expected-mismatches.json`).
#
# These helpers are intentionally unexported. They classify per-(case_id, field)
# parity divergences between mixeff and lme4-family references so that test
# code can distinguish documented optimizer drift from unclassified regressions.
# See planning/PRD.md §3 (no bit-exactness) and §11 (test strategy) for the
# governing contract; bead bd-01KQF83XAN4CGAS176JRX8CR7E captures the rationale.

mm_parity_ledger_path <- function() {
  candidates <- c(
    system.file("extdata", "expected-mismatches.json", package = "mixeff"),
    file.path("inst", "extdata", "expected-mismatches.json")
  )
  hit <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(hit)) {
    stop(
      "expected-mismatches.json is not available; cannot consult parity ledger.",
      call. = FALSE
    )
  }
  hit[[1L]]
}

mm_parity_ledger <- function() {
  jsonlite::read_json(mm_parity_ledger_path(), simplifyVector = FALSE)
}

mm_parity_lookup <- function(case_id, field, ledger = mm_parity_ledger()) {
  hit <- Filter(
    function(entry) identical(entry$case_id, case_id) && identical(entry$field, field),
    ledger$mismatches
  )
  if (!length(hit)) return(NULL)
  hit[[1L]]
}
