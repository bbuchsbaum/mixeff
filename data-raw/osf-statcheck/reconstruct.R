#!/usr/bin/env Rscript
# Download + reconstruct the OSF Study-3 merged statcheck/badges data and write
# the Period-2 modeling fixture inst/extdata/osf-statcheck-t2.csv.
#
# Provenance and the lossy-re-save caveat are documented in README.md.
# Needs network access (OSF). Run from the package root:
#   Rscript data-raw/osf-statcheck/reconstruct.R
#
# Tracked by mote bd-01KT3ZB649HHX4ZYMKTQ5A18XX.

osf_url <- "https://osf.io/download/qn645/"   # 170329MergedDataStatcheckBadges.txt

tmp <- tempfile(fileext = ".txt")
utils::download.file(osf_url, tmp, mode = "wb", quiet = TRUE)

# The OSF file is UTF-16. Each physical line is the authors' original write.csv
# line wrapped in one pair of outer quotes with inner quotes doubled. For rows
# whose title (Source) contains commas the re-save dropped the title's
# field-quoting, so we parse the trailing 19 columns from the RIGHT and rejoin
# the leading fields as Source.
con <- file(tmp, encoding = "UTF-16LE")
raw <- readLines(con, warn = FALSE)
close(con)
raw <- raw[nzchar(raw)]
# strip BOM if present on first line
raw[1] <- sub("^﻿", "", raw[1])

inner <- gsub('""', '"', substr(raw, 2L, nchar(raw) - 1L), fixed = TRUE)
hdr <- scan(text = inner[1], what = "character", sep = ",", quote = "\"",
            quiet = TRUE)
ncol  <- length(hdr)          # 20
ntail <- ncol - 1L            # 19 well-formed trailing fields

rows <- lapply(inner[-1], function(ln) {
  f <- scan(text = ln, what = "character", sep = ",", quote = "\"",
            quiet = TRUE, strip.white = FALSE)
  n <- length(f)
  c(paste(f[seq_len(n - ntail)], collapse = ","), f[(n - ntail + 1L):n])
})
df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(df) <- hdr

stopifnot(nrow(df) == 20926L)
df$Error           <- as.integer(as.logical(df$Error))
df$DecisionError   <- as.integer(as.logical(df$DecisionError))
df$Year            <- as.integer(df$Year)
df$Period          <- as.integer(df$Period)
df$OpenData        <- as.integer(df$OpenData)
df$OpenMaterials   <- as.integer(df$OpenMaterials)
df$Preregistration <- as.integer(df$Preregistration)
stopifnot(!anyNA(df[, c("Error", "DecisionError", "Year", "Period",
                        "OpenData", "OpenMaterials", "Preregistration")]))

t2 <- df[df$Period == 2L, ]
t2$gid <- as.integer(factor(t2$Source))   # anonymized, deterministic group id
out <- t2[, c("gid", "Error", "DecisionError", "Year",
              "OpenData", "OpenMaterials", "Preregistration")]

dest <- file.path("inst", "extdata", "osf-statcheck-t2.csv")
utils::write.csv(out, dest, row.names = FALSE)
message(sprintf("wrote %s: %d rows, %d groups, %d Error events",
                dest, nrow(out), length(unique(out$gid)), sum(out$Error)))

# Fidelity check against the published number (requires lme4).
if (requireNamespace("lme4", quietly = TRUE)) {
  g0 <- lme4::glmer(Error ~ OpenPractice * Year + (1 | gid),
                    data = transform(out,
                      OpenPractice = OpenData == 1 | OpenMaterials == 1 |
                        Preregistration == 1),
                    family = binomial("logit"), nAGQ = 0)
  b <- summary(g0)$coefficients["OpenPracticeTRUE:Year", ]
  message(sprintf("fidelity (glmer nAGQ=0): OpenPractice:Year b=%.4f Z=%.3f p=%.4f  (paper: 0.7958 / 1.825 / .0679)",
                  b["Estimate"], b["z value"], b["Pr(>|z|)"]))
}
