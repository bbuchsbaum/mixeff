#!/usr/bin/env Rscript
# Download the OSF "Willingness to wait" (node ftexh) trial-level data and write
# the slimmed modeling fixtures used by the parity test.
#
# Writes:
#   tests/fixtures/osf_willingness_to_wait_study1a.csv  (1427 rows)
#   tests/fixtures/osf_willingness_to_wait_study1b.csv   (1600 rows)
#
# Provenance is documented in README.md. Needs network access (OSF).
# Run from the package root:
#   Rscript data-raw/osf-willingness-to-wait/reconstruct.R
#
# Tracked by mote bd-01KT3ZRCKWRZQFA4W7TXTGWAZ0.

sources <- list(
  study1a = "https://osf.io/download/u2zt4/",  # FullData_study1a.csv
  study1b = "https://osf.io/download/3bxkt/"   # FullData_replication.csv
)

# Only the columns the 9 published models touch. Keeps the fixture small and
# the modeling intent obvious; everything else (demographics, ARMQ items, raw
# response text) is dropped.
keep <- c("ID", "Title", "wait_choice", "Enjoyment", "arousal",
          "Q1_correct", "Q2_correct", "SVScore")

for (which in names(sources)) {
  tmp <- tempfile(fileext = ".csv")
  utils::download.file(sources[[which]], tmp, mode = "wb", quiet = TRUE)
  d <- utils::read.csv(tmp, stringsAsFactors = FALSE)
  missing <- setdiff(keep, names(d))
  if (length(missing)) {
    stop("OSF ", which, " file is missing expected columns: ",
         paste(missing, collapse = ", "))
  }
  d <- d[, keep]
  dest <- file.path("tests", "fixtures",
                    sprintf("osf_willingness_to_wait_%s.csv", which))
  utils::write.csv(d, dest, row.names = FALSE, na = "NA")
  message(sprintf("wrote %s (%d rows, %d cols)", dest, nrow(d), ncol(d)))
}
