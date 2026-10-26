# Shared formulas, factor coding, and published-table anchors for
# Brannick & LaCroix (2025), Scientific Reports 15:27262.
# Auto-loaded by testthat; also sourced by
# data-raw/osf-brannick-lacroix-2025/{reconstruct,parity-harness,published-models}.R.
#
# Paper coding (Tables 3–4, High-dynamic / Happy as music reference):
#   Music HAPPY < SAD < SILENCE
#   ANT   ANT_1 < ANT_2
#   Cue   Double < No          (alerting)
#   Cue   Center < Valid       (orienting)
#
# Tracked by mote bd-01M0842MM295RFXXN4C1ZDE9CV.

osf_bl_view_only <- "0c0e8abeb85c4f97868a922eff278fa6"

osf_bl_sources <- list(
  project = "https://osf.io/gf3km/",
  analysis = "https://osf.io/7ryjh/",
  database = sprintf(
    "https://osf.io/download/9x5rz/?view_only=%s", osf_bl_view_only
  ),
  script_2024 = sprintf(
    "https://osf.io/download/35dfe/?view_only=%s", osf_bl_view_only
  ),
  supplemental_2025 = sprintf(
    "https://osf.io/download/685eda7ea310effdd5d1adee/?view_only=%s",
    osf_bl_view_only
  ),
  participant_xlsx = sprintf(
    "https://osf.io/download/tdsjr/?view_only=%s", osf_bl_view_only
  )
)

# Paper Cue-Onset WOI is written both as 2401–2901 and 2401–2902. Table 3/4
# nobs (34257 / 33449) are the reconstruction targets; the inclusive 2902
# bound is the one that can reach those counts.
osf_bl_cue_woi <- c(2401, 2902)

osf_bl_paper_nobs <- list(
  alerting = 34257L,
  orienting = 33449L
)

# Table 3 High-ref coefficients (Scientific Reports Table 3).
# The prose at p.7 quotes z = -2.657 / p = .008 for the three-way High x T2 x
# No-cue term; the table itself is -2.389 / .017. We treat the table as the
# numeric source of truth and record the prose mismatch in the README.
osf_bl_table3_high_ref <- c(
  `(Intercept)`                       = -4.487,
  `ANTT2`                             =  0.243,
  `CueNo`                             =  0.781,
  `MusicLow`                          = -0.206,
  `MusicSilence`                      = -0.011,
  `ANTT2:CueNo`                       = -0.146,
  `ANTT2:MusicLow`                    =  0.762,
  `ANTT2:MusicSilence`                =  0.109,
  `CueNo:MusicLow`                    =  0.532,
  `CueNo:MusicSilence`                = -0.163,
  `ANTT2:CueNo:MusicLow`              = -0.877,
  `ANTT2:CueNo:MusicSilence`          = -0.171
)

osf_bl_table4_high_ref <- c(
  `(Intercept)`                       = -4.084,
  `ANTT2`                             = -0.051,
  `CueValid`                          = -0.437,
  `MusicLow`                          =  0.160,
  `MusicSilence`                      = -0.042
)

osf_bl_forms <- list(
  alerting = Blink ~ Music * ANT * Cue + (1 + ANT | P),
  orienting = Blink ~ Music * ANT * Cue + (1 + ANT | P)
)

osf_bl_normalize_music <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("HIGH", "HIGH DYNAMIC", "HAPPY")] <- "HAPPY"
  x[x %in% c("LOW", "LOW DYNAMIC", "SAD")] <- "SAD"
  x[x %in% c("CONTROL", "NO MUSIC", "SILENCE")] <- "SILENCE"
  factor(x, levels = c("HAPPY", "SAD", "SILENCE"))
}

osf_bl_normalize_ant <- function(x) {
  x <- toupper(gsub("[^A-Z0-9]", "", as.character(x)))
  x[x %in% c("T1", "1", "ANT1")] <- "ANT_1"
  x[x %in% c("T2", "2", "ANT2")] <- "ANT_2"
  factor(x, levels = c("ANT_1", "ANT_2"))
}

osf_bl_normalize_cue <- function(x, which = c("alerting", "orienting")) {
  which <- match.arg(which)
  x <- trimws(as.character(x))
  key <- tolower(x)
  key[key %in% c("no", "nocue", "no_cue", "no cue")] <- "No"
  key[key %in% c("double", "doublecue", "double_cue")] <- "Double"
  key[key %in% c("center", "centre")] <- "Center"
  key[key %in% c("valid", "spatial", "validspatial")] <- "Valid"
  if (identical(which, "alerting")) {
    factor(key, levels = c("Double", "No"))
  } else {
    factor(key, levels = c("Center", "Valid"))
  }
}

osf_bl_prepare_alerting <- function(d) {
  out <- data.frame(
    P = as.character(d$P),
    Blink = as.integer(as.character(d$Blink)),
    ANT = osf_bl_normalize_ant(d$ANT),
    Cue = osf_bl_normalize_cue(d$Cue, "alerting"),
    Music = osf_bl_normalize_music(d$Music),
    stringsAsFactors = FALSE
  )
  out[stats::complete.cases(out), , drop = FALSE]
}

osf_bl_prepare_orienting <- function(d) {
  out <- data.frame(
    P = as.character(d$P),
    Blink = as.integer(as.character(d$Blink)),
    ANT = osf_bl_normalize_ant(d$ANT),
    Cue = osf_bl_normalize_cue(d$Cue, "orienting"),
    Music = osf_bl_normalize_music(d$Music),
    stringsAsFactors = FALSE
  )
  out[stats::complete.cases(out), , drop = FALSE]
}

osf_bl_map_fixef <- function(beta) {
  mapped <- names(beta)
  mapped <- gsub("ANTANT_2", "ANTT2", mapped, fixed = TRUE)
  mapped <- gsub("MusicSAD", "MusicLow", mapped, fixed = TRUE)
  mapped <- gsub("MusicSILENCE", "MusicSilence", mapped, fixed = TRUE)
  names(beta) <- mapped
  beta
}
