#!/usr/bin/env Rscript
# Download the OSF gf3km ANT database and write the Cue-Onset fixtures.
#
# Writes:
#   tests/fixtures/osf_bl_alerting.csv.gz
#   tests/fixtures/osf_bl_orienting.csv.gz
#
# Executive-control and RT slices are not committed: their nobs do not
# recover Tables 5–6. Provenance is in README.md. Needs network (OSF).
# Run from the package root:
#   Rscript data-raw/osf-brannick-lacroix-2025/reconstruct.R
#
# Tracked by mote bd-01M0842MM295RFXXN4C1ZDE9CV.

root <- if (file.exists("DESCRIPTION")) {
  getwd()
} else {
  stop("Run reconstruct.R from the mixeff package root")
}

source(file.path(root, "data-raw", "osf-brannick-lacroix-2025",
                 "published-models.R"), local = TRUE)

cache_dir <- file.path(root, "tmp", "osf-brannick-lacroix-2025")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
raw_csv <- file.path(cache_dir, "8-25_Database.csv")

local_override <- Sys.getenv("OSF_BL_DATABASE", unset = "")
if (!nzchar(local_override) && length(commandArgs(TRUE))) {
  local_override <- commandArgs(TRUE)[[1L]]
}

download_osf <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 1e6) {
    message("using cached ", dest, " (", file.info(dest)$size, " bytes)")
    return(invisible(dest))
  }
  message("downloading ", url)
  utils::download.file(url, dest, mode = "wb", quiet = FALSE)
  dest
}

if (nzchar(local_override)) {
  if (!file.exists(local_override)) {
    stop("OSF_BL_DATABASE / argument does not exist: ", local_override)
  }
  message("using local ANT database ", local_override)
  raw_csv <- local_override
} else {
  download_osf(osf_bl_sources$database, raw_csv)
}
if (file.info(raw_csv)$size < 1e6) {
  stop("Downloaded ANT database looks too small: ", raw_csv)
}

header <- names(utils::read.csv(raw_csv, nrows = 1L, check.names = FALSE,
                                stringsAsFactors = FALSE,
                                fileEncoding = "UTF-8-BOM"))
message("source columns: ", paste(header, collapse = ", "))

need <- function(candidates, what) {
  hit <- candidates[candidates %in% header][1L]
  if (is.na(hit)) {
    stop("8-25_Database.csv is missing ", what, ". Available: ",
         paste(header, collapse = ", "))
  }
  hit
}

col_p     <- need(c("P", "Participant", "ID", "id"), "participant id")
col_blink <- need(c("BlinkLogistic", "Blink_Logistic", "Blink"),
                  "blink indicator")
col_ant   <- need(c("ANT_Cx", "ANT", "ANT_Admin"), "ANT administration")
col_music <- need(c("Music_Cx", "Music", "Listening"), "music condition")
col_cue   <- need(c("CUE_TYPE", "Cue", "CueType"), "cue type")
col_time  <- need(c("TimeStamp", "Timestamp", "TIME"), "timestamp")
col_alert <- if ("ALERTING_ORIENTING" %in% header) "ALERTING_ORIENTING" else NA_character_
col_orient <- if ("Orienting_ValidCenter" %in% header) {
  "Orienting_ValidCenter"
} else {
  NA_character_
}

keep <- unique(c(col_p, col_blink, col_ant, col_music, col_cue, col_time,
                 col_alert, col_orient))
keep <- keep[!is.na(keep)]

classes <- rep("NULL", length(header))
names(classes) <- header
classes[keep] <- "character"

message("reading ", length(keep), " columns from ", raw_csv)
d <- utils::read.csv(raw_csv, colClasses = classes, check.names = FALSE,
                     stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
                     na.strings = c("", "NA", "NaN"))
names(d)[names(d) == col_p]     <- "P"
names(d)[names(d) == col_blink] <- "BlinkRaw"
names(d)[names(d) == col_ant]   <- "ANT"
names(d)[names(d) == col_music] <- "Music"
names(d)[names(d) == col_cue]   <- "Cue"
names(d)[names(d) == col_time]  <- "TimeStamp"
if (!is.na(col_alert)) {
  names(d)[names(d) == col_alert] <- "ALERTING_ORIENTING"
}
if (!is.na(col_orient)) {
  names(d)[names(d) == col_orient] <- "Orienting_ValidCenter"
}

d$TimeStamp <- suppressWarnings(as.numeric(d$TimeStamp))
blink_chr <- toupper(trimws(as.character(d$BlinkRaw)))
d$Blink <- ifelse(blink_chr %in% c("1", "TRUE", "YES"), 1L,
                  ifelse(blink_chr %in% c("0", "FALSE", "NO"), 0L, NA_integer_))
d$P <- as.integer(factor(d$P))

in_cue_woi <- !is.na(d$TimeStamp) &
  d$TimeStamp >= osf_bl_cue_woi[[1L]] &
  d$TimeStamp <= osf_bl_cue_woi[[2L]]

alerting_rows <- in_cue_woi
if ("ALERTING_ORIENTING" %in% names(d)) {
  alerting_rows <- alerting_rows &
    toupper(as.character(d$ALERTING_ORIENTING)) == "ALERTING"
} else {
  cue_key <- tolower(as.character(d$Cue))
  alerting_rows <- alerting_rows &
    cue_key %in% c("double", "no", "nocue", "no cue")
}

orienting_rows <- in_cue_woi
if ("Orienting_ValidCenter" %in% names(d)) {
  orienting_rows <- orienting_rows &
    grepl("valid", tolower(as.character(d$Orienting_ValidCenter)), fixed = TRUE)
} else {
  cue_key <- tolower(as.character(d$Cue))
  orienting_rows <- orienting_rows & cue_key %in% c("center", "centre", "valid")
}

write_gz <- function(df, dest) {
  con <- gzfile(dest, open = "wt")
  utils::write.csv(df, con, row.names = FALSE, na = "NA")
  close(con)
  message(sprintf("wrote %s (%d rows, %d cols, %.1f KB)",
                  dest, nrow(df), ncol(df), file.info(dest)$size / 1024))
}

fx_dir <- file.path(root, "tests", "fixtures")
dir.create(fx_dir, showWarnings = FALSE)

alerting <- data.frame(
  P = d$P[alerting_rows],
  Blink = d$Blink[alerting_rows],
  ANT = d$ANT[alerting_rows],
  Cue = d$Cue[alerting_rows],
  Music = d$Music[alerting_rows],
  stringsAsFactors = FALSE
)
alerting <- alerting[!is.na(alerting$Blink), , drop = FALSE]
write_gz(alerting, file.path(fx_dir, "osf_bl_alerting.csv.gz"))

orienting <- data.frame(
  P = d$P[orienting_rows],
  Blink = d$Blink[orienting_rows],
  ANT = d$ANT[orienting_rows],
  Cue = d$Cue[orienting_rows],
  Music = d$Music[orienting_rows],
  stringsAsFactors = FALSE
)
orienting <- orienting[!is.na(orienting$Blink), , drop = FALSE]
write_gz(orienting, file.path(fx_dir, "osf_bl_orienting.csv.gz"))

message(sprintf("nobs alerting  %7d  (paper %d, delta %+d)",
                nrow(alerting), osf_bl_paper_nobs$alerting,
                nrow(alerting) - osf_bl_paper_nobs$alerting))
message(sprintf("nobs orienting %7d  (paper %d, delta %+d)",
                nrow(orienting), osf_bl_paper_nobs$orienting,
                nrow(orienting) - osf_bl_paper_nobs$orienting))

if (requireNamespace("lme4", quietly = TRUE)) {
  message("fidelity: glmer Table 3 High-ref (nAGQ=1). Printed table is not recovered.")
  a <- osf_bl_prepare_alerting(alerting)
  g0 <- suppressWarnings(lme4::glmer(
    osf_bl_forms$alerting, data = a, family = binomial("logit"),
    control = lme4::glmerControl(optimizer = "bobyqa")
  ))
  b <- lme4::fixef(g0)
  message(sprintf("  intercept=%.3f (paper -4.487)  CueNo=%.3f (paper 0.781)",
                  unname(b[["(Intercept)"]]), unname(b[["CueNo"]])))
}
