## Parity probe: lmm-dyestuff2-singular
## Dataset: Dyestuff2 (near-zero variance / boundary case)
## Formula: Yield ~ 1 + (1|Batch)
## Compares lme4/lmerTest vs mixeff on key quantities

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

# Resolve namespace conflicts: use explicit lme4:: for lme4 objects
# mixeff masks fixef, ranef, VarCorr, sigma — we wrap lme4 calls explicitly

cat("=== SESSION INFO ===\n")
cat("lme4:", as.character(packageVersion("lme4")), "\n")
cat("lmerTest:", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff:", as.character(packageVersion("mixeff")), "\n\n")

## ---- Dataset ----
data(Dyestuff2, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(Dyestuff2), "\n")
cat("ncol:", ncol(Dyestuff2), "\n")
print(head(Dyestuff2))
cat("\n")

## ============================================================
## 1. lme4 fit
## ============================================================
cat("=== lme4 FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(Yield ~ 1 + (1|Batch), data = Dyestuff2, REML = TRUE)
})
cat("lme4 wall time (sec):", t_lme4["elapsed"], "\n\n")

cat("--- lme4 summary ---\n")
print(summary(fit_lme4))
cat("\n")

cat("--- lme4 fixef ---\n")
print(lme4::fixef(fit_lme4))
cat("\n")

cat("--- lme4 VarCorr ---\n")
print(lme4::VarCorr(fit_lme4))
cat("\n")

cat("--- lme4 theta ---\n")
print(lme4::getME(fit_lme4, "theta"))
cat("\n")

cat("--- lme4 sigma ---\n")
print(stats::sigma(fit_lme4))
cat("\n")

cat("--- lme4 logLik ---\n")
print(logLik(fit_lme4))
cat("\n")

cat("--- lme4 AIC/BIC ---\n")
cat("AIC:", AIC(fit_lme4), "\n")
cat("BIC:", BIC(fit_lme4), "\n\n")

cat("--- lme4 vcov (fixed) ---\n")
print(vcov(fit_lme4))
cat("\n")

cat("--- lme4 ranef ---\n")
print(lme4::ranef(fit_lme4))
cat("\n")

cat("--- lme4 fitted (first 10) ---\n")
print(head(fitted(fit_lme4), 10))
cat("\n")

cat("--- lme4 isSingular ---\n")
cat("isSingular:", isSingular(fit_lme4), "\n\n")

cat("--- lme4 convergence warnings ---\n")
msgs <- tryCatch(fit_lme4@optinfo$conv$lme4$messages, error = function(e) NULL)
if (is.null(msgs) || length(msgs) == 0) {
  cat("No convergence warnings\n")
} else {
  cat(paste(msgs, collapse = "\n"), "\n")
}
cat("\n")

## ============================================================
## 2. mixeff fit
## ============================================================
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    lmm(Yield ~ 1 + (1|Batch), data = Dyestuff2,
        REML = TRUE, control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff wall time (sec):", t_mixeff["elapsed"], "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR !!!\n")
  cat("Class:", class(fit_mm), "\n")
  cat("Message:", conditionMessage(fit_mm), "\n\n")
  # Still try to continue with NAs for comparison
  mm_ok <- FALSE
} else {
  mm_ok <- TRUE
  cat("--- mixeff fit_status ---\n")
  cat("fit_status:", fit_mm$fit_status, "\n\n")

  cat("--- mixeff summary ---\n")
  print(tryCatch(summary(fit_mm), error = function(e) e))
  cat("\n")

  cat("--- mixeff fixef ---\n")
  print(tryCatch(fixef(fit_mm), error = function(e) e))
  cat("\n")

  cat("--- mixeff VarCorr ---\n")
  print(tryCatch(VarCorr(fit_mm), error = function(e) e))
  cat("\n")

  cat("--- mixeff theta ---\n")
  print(tryCatch(fit_mm$theta, error = function(e) e))
  cat("\n")

  cat("--- mixeff sigma ---\n")
  print(tryCatch(fit_mm$sigma, error = function(e) e))
  cat("\n")

  cat("--- mixeff logLik ---\n")
  print(tryCatch(logLik(fit_mm), error = function(e) e))
  cat("\n")

  cat("--- mixeff AIC/BIC ---\n")
  cat("AIC:", tryCatch(AIC(fit_mm), error = function(e) NA), "\n")
  cat("BIC:", tryCatch(BIC(fit_mm), error = function(e) NA), "\n\n")

  cat("--- mixeff vcov (fixed) ---\n")
  print(tryCatch(vcov(fit_mm), error = function(e) e))
  cat("\n")

  cat("--- mixeff ranef ---\n")
  print(tryCatch(ranef(fit_mm), error = function(e) e))
  cat("\n")

  cat("--- mixeff fitted (first 10) ---\n")
  print(tryCatch(head(fitted(fit_mm), 10), error = function(e) e))
  cat("\n")
}

## ============================================================
## 3. Numeric comparison
## ============================================================
cat("=== NUMERIC COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(name, lme4_val, mm_val, tol = 1e-6) {
  lme4_v <- as.numeric(lme4_val)
  mm_v   <- as.numeric(mm_val)
  if (any(is.na(mm_v))) {
    cat(sprintf("  %-20s  lme4=%-12s  mixeff=NA  STATUS=MISSING\n",
                name, paste(round(lme4_v, 6), collapse=",")))
    return(invisible(NA))
  }
  if (length(lme4_v) != length(mm_v)) {
    cat(sprintf("  %-20s  lme4 len=%d  mixeff len=%d  STATUS=LENGTH_MISMATCH\n",
                name, length(lme4_v), length(mm_v)))
    return(invisible(NA))
  }
  d <- max(abs(lme4_v - mm_v))
  status <- if (d <= tol) "WITHIN-TOL" else "DIVERGENT"
  cat(sprintf("  %-20s  lme4=%-12s  mixeff=%-12s  maxAbsDiff=%.2e  tol=%.0e  %s\n",
              name,
              paste(round(lme4_v, 6), collapse=","),
              paste(round(mm_v, 6), collapse=","),
              d, tol, status))
  invisible(d)
}

if (mm_ok) {
  cat("\n--- Fixed effects ---\n")
  compare("fixef (Intercept)",
          lme4::fixef(fit_lme4),
          fixef(fit_mm),
          tol = tols$fixef)

  cat("\n--- SE of fixed effects ---\n")
  lme4_se <- sqrt(diag(as.matrix(vcov(fit_lme4))))
  mm_se   <- fit_mm$std_errors
  compare("SE (Intercept)", lme4_se, mm_se, tol = tols$fixef)

  cat("\n--- vcov (fixed) ---\n")
  lme4_vcov <- as.numeric(vcov(fit_lme4))
  mm_vcov   <- as.numeric(vcov(fit_mm))
  compare("vcov[1,1]", lme4_vcov[1], mm_vcov[1], tol = tols$fixef)

  cat("\n--- theta ---\n")
  compare("theta",
          lme4::getME(fit_lme4, "theta"),
          fit_mm$theta,
          tol = tols$theta)

  cat("\n--- sigma ---\n")
  compare("sigma",
          stats::sigma(fit_lme4),
          fit_mm$sigma,
          tol = tols$sigma)

  cat("\n--- logLik ---\n")
  compare("logLik",
          as.numeric(logLik(fit_lme4)),
          as.numeric(logLik(fit_mm)),
          tol = tols$logLik)

  cat("\n--- AIC / BIC ---\n")
  compare("AIC", AIC(fit_lme4), AIC(fit_mm), tol = 1e-3)
  compare("BIC", BIC(fit_lme4), BIC(fit_mm), tol = 1e-3)

  cat("\n--- VarCorr (Batch var) ---\n")
  lme4_batchvar <- as.data.frame(lme4::VarCorr(fit_lme4))$vcov[1]
  mm_batchvar   <- tryCatch({
    vc <- VarCorr(fit_mm)
    # try to extract Batch variance
    if (is.list(vc)) {
      as.numeric(vc[[1]])
    } else {
      as.numeric(vc)
    }
  }, error = function(e) NA)
  compare("VarCorr Batch", lme4_batchvar, mm_batchvar, tol = tols$theta * max(stats::sigma(fit_lme4)^2, 1))

  cat("\n--- ranef correlation ---\n")
  lme4_ranef <- as.numeric(unlist(lme4::ranef(fit_lme4)))
  mm_ranef_raw <- tryCatch(ranef(fit_mm), error = function(e) NULL)
  if (!is.null(mm_ranef_raw)) {
    mm_ranef <- as.numeric(unlist(mm_ranef_raw))
    compare("ranef", lme4_ranef, mm_ranef, tol = 0.05)
  } else {
    cat("  ranef: mixeff returned NULL or error\n")
  }

  cat("\n--- fitted values (max abs diff) ---\n")
  lme4_fv <- as.numeric(fitted(fit_lme4))
  mm_fv   <- as.numeric(fitted(fit_mm))
  if (length(lme4_fv) == length(mm_fv)) {
    d_fit <- max(abs(lme4_fv - mm_fv))
    cat(sprintf("  fitted maxAbsDiff=%.6f\n", d_fit))
  } else {
    cat("  fitted length mismatch\n")
  }
} else {
  cat("  [skipped — mixeff returned an error]\n")
}

cat("\n=== BOUNDARY / SINGULAR BEHAVIOUR ===\n")
cat("lme4 isSingular:", isSingular(fit_lme4), "\n")
if (mm_ok) {
  cat("mixeff fit_status:", fit_mm$fit_status, "\n")
  # Does mixeff emit a diagnostic / warning about boundary?
  cat("mixeff theta (should be near 0):", fit_mm$theta, "\n")
  cat("lme4  theta (should be near 0):", lme4::getME(fit_lme4, "theta"), "\n")
}

cat("\n=== SPEED ===\n")
cat(sprintf("lme4  elapsed: %.4f s\n", t_lme4["elapsed"]))
cat(sprintf("mixeff elapsed: %.4f s\n", t_mixeff["elapsed"]))
cat(sprintf("ratio mixeff/lme4: %.2fx\n", t_mixeff["elapsed"] / t_lme4["elapsed"]))

cat("\n=== DONE ===\n")
