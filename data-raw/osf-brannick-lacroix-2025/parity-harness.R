#!/usr/bin/env Rscript
# Offline glmer vs mixeff parity harness on the committed Brannick &
# LaCroix 2025 Cue-Onset fixtures. Run from the package root:
#   Rscript data-raw/osf-brannick-lacroix-2025/parity-harness.R
#
# Tracked by mote bd-01M0842MM295RFXXN4C1ZDE9CV.

suppressMessages({
  library(lme4)
  library(mixeff)
})

root <- if (file.exists("DESCRIPTION")) getwd() else {
  stop("Run parity-harness.R from the mixeff package root")
}
source(file.path(root, "data-raw", "osf-brannick-lacroix-2025",
                 "published-models.R"), local = TRUE)

fx <- function(name) {
  f <- file.path(root, "tests", "fixtures", sprintf("osf_bl_%s.csv.gz", name))
  if (!file.exists(f)) {
    stop("fixture not found: ", f, " (run reconstruct.R first)")
  }
  utils::read.csv(gzfile(f), stringsAsFactors = FALSE)
}

ctl <- mm_control(verbose = -1, max_feval = 100000L)
bin <- binomial("logit")

compare_glmm <- function(id, d, fo) {
  g <- suppressWarnings(lme4::glmer(
    fo, data = d, family = bin,
    control = lme4::glmerControl(optimizer = "bobyqa")
  ))
  m <- glmm(fo, data = d, family = bin, method = "joint_laplace",
            control = ctl)
  bg <- unname(lme4::fixef(g))
  bm <- unname(fixef(m))
  cat(sprintf("%-12s n=%6d  dlogLik=%+8.1e  max|dFix|=%8.1e  status=%s\n",
              id, nrow(d),
              as.numeric(logLik(m)) - as.numeric(logLik(g)),
              max(abs(bm - bg)),
              m$fit_status))
  invisible(list(glmer = g, mixeff = m))
}

cat(sprintf("%-12s %s\n", "model", "parity"))
compare_glmm("alerting", osf_bl_prepare_alerting(fx("alerting")),
             osf_bl_forms$alerting)
compare_glmm("orienting", osf_bl_prepare_orienting(fx("orienting")),
             osf_bl_forms$orienting)
