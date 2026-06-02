#!/usr/bin/env Rscript
# Offline glmer-vs-mixeff parity harness on the committed OSF statcheck fixture.
#   Rscript data-raw/osf-statcheck/parity-harness.R
# Reproduces the in-the-wild findings recorded on mote bd-01KT3ZB649HHX4ZYMKTQ5A18XX.

suppressMessages({
  library(lme4)
  library(mixeff)
})

f <- system.file("extdata", "osf-statcheck-t2.csv", package = "mixeff")
if (!nzchar(f)) f <- "inst/extdata/osf-statcheck-t2.csv"
d <- utils::read.csv(f)
d$Source       <- factor(d$gid)
d$OpenPractice <- d$OpenData == 1 | d$OpenMaterials == 1 | d$Preregistration == 1
d$OtherBadge   <- d$OpenMaterials == 1 | d$Preregistration == 1
d$cYear        <- d$Year - 2015L           # centering -> well-conditioned
ix <- function(b) as.numeric(b[grep(":", names(b))])
ctl <- mm_control(verbose = -1, max_feval = 50000L)

cat("fixture:", nrow(d), "rows,", nlevels(d$Source), "groups,",
    sum(d$Error), "Error events\n\n")

## ---- 1. Estimation parity on well-conditioned (centered) models -------------
models <- list(
  "Error ~ OpenPractice + (1|Source)"          = Error ~ OpenPractice + (1 | Source),
  "Error ~ OpenPractice + cYear + (1|Source)"  = Error ~ OpenPractice + cYear + (1 | Source),
  "Error ~ OpenPractice * cYear + (1|Source)"  = Error ~ OpenPractice * cYear + (1 | Source),
  "Error ~ OpenData+OtherBadge+cYear+(1|Source)" =
    Error ~ OpenData + OtherBadge + cYear + (1 | Source)
)
cat("== Estimation parity (joint_laplace vs glmer, centered) ==\n")
for (nm in names(models)) {
  fo <- models[[nm]]
  g <- glmer(fo, data = d, family = binomial("logit"))
  m <- glmm(fo, data = d, family = binomial("logit"),
            method = "joint_laplace", control = ctl)
  bg <- unname(lme4::fixef(g)); bm <- unname(fixef(m))
  cat(sprintf("  %-46s max|dEst|=%.1e  dlogLik=%.1e\n", nm,
              max(abs(bm - bg)),
              as.numeric(logLik(m)) - as.numeric(logLik(g))))
}

## ---- 2. Gap E: raw-Year sub-optimal convergence reported as success ---------
cat("\n== Gap E: raw vs centered Year (interaction is offset-invariant) ==\n")
mr <- glmm(Error ~ OpenPractice * Year  + (1 | Source), d, binomial("logit"), method = "joint_laplace", control = ctl)
mc <- glmm(Error ~ OpenPractice * cYear + (1 | Source), d, binomial("logit"), method = "joint_laplace", control = ctl)
gc <- glmer(Error ~ OpenPractice * cYear + (1 | Source), d, family = binomial("logit"))
ll_raw <- suppressWarnings(as.numeric(logLik(mr)))
if (length(ll_raw) != 1L) ll_raw <- NA_real_
cat(sprintf("  raw Year : interaction=%.5f  status=%s  logLik=%s\n",
            ix(fixef(mr)), fit_status(mr), format(ll_raw)))
cat(sprintf("  centered : interaction=%.5f  status=%s  logLik=%.4f  (glmer=%.5f)\n",
            ix(fixef(mc)), fit_status(mc), as.numeric(logLik(mc)), ix(lme4::fixef(gc))))
cat("  -> raw != centered with status=converged_interior + non-finite logLik = upstream bd-01KT3Z64AY45NHA5144G2ZBMSY\n")

## ---- 3. Gap B: GLMM SEs uncertified ----------------------------------------
cat("\n== Gap B: GLMM inference is uncertified ==\n")
m1 <- glmm(Error ~ OpenPractice + (1 | Source), d, binomial("logit"), method = "joint_laplace", control = ctl)
print(summary(m1)$coefficients)
cat("  -> method=not_computed, z/p=NA; vcov exposed but uncertified = upstream bd-01KT3Z64YE5QN7626PQRJSJJVA\n")
