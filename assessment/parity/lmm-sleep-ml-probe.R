## Empirical parity probe: lmm-sleep-ml
## Dataset: sleepstudy  Formula: Reaction ~ Days + (Days|Subject)  REML = FALSE (ML)
## Compares lme4/lmerTest vs mixeff on all relevant quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

## ── 1. Data ──────────────────────────────────────────────────────────────────
data(sleepstudy, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(sleepstudy), "  ncol:", ncol(sleepstudy), "\n")
print(head(sleepstudy))
cat("\n")

## ── 2. Fit lme4 (ML) ─────────────────────────────────────────────────────────
cat("=== lme4 FIT (REML=FALSE) ===\n")
NREPS <- 10L
t_lme4 <- system.time({
  for (i in seq_len(NREPS)) {
    fit_lme4 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                     REML = FALSE)
  }
})
cat("lme4 wall-clock (seconds, mean of", NREPS, "reps):",
    t_lme4["elapsed"] / NREPS, "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4))
cat("-- SE (sqrt diag vcov) --\n"); print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov --\n"); print(as.matrix(vcov(fit_lme4)))
cat("-- theta (Cholesky factor) --\n"); print(lme4::getME(fit_lme4, "theta"))
cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))
cat("-- sigma --\n"); print(sigma(fit_lme4))
cat("-- logLik --\n"); print(logLik(fit_lme4))
cat("-- AIC --\n");    print(AIC(fit_lme4))
cat("-- BIC --\n");    print(BIC(fit_lme4))
cat("-- ranef (first 3 subjects) --\n"); print(head(lme4::ranef(fit_lme4)$Subject, 3))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
cat("-- convergence --\n"); cat("No convergence warning =",
  length(fit_lme4@optinfo$conv$lme4$messages) == 0, "\n\n")

## ── 3. Fit mixeff (ML) ───────────────────────────────────────────────────────
cat("=== mixeff FIT (REML=FALSE) ===\n")
t_mixeff <- system.time({
  for (i in seq_len(NREPS)) {
    fit_mm <- tryCatch(
      lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE,
          control = mm_control(verbose = -1L)),
      error = function(e) e
    )
  }
})
cat("mixeff wall-clock (seconds, mean of", NREPS, "reps):",
    t_mixeff["elapsed"] / NREPS, "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR !!!\n")
  cat(conditionMessage(fit_mm), "\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat("\n=== SPEED RATIO ===\n")
  cat("lme4 mean elapsed:", t_lme4["elapsed"] / NREPS, "s\n")
  cat("mixeff mean elapsed:", t_mixeff["elapsed"] / NREPS, "s\n")
  cat("ratio (lme4/mixeff):", (t_lme4["elapsed"] / NREPS) / (t_mixeff["elapsed"] / NREPS), "\n")
  quit(status = 0)
}

cat("fit_status:", fit_mm$fit_status, "\n\n")

cat("-- fixef --\n"); print(fixef(fit_mm))
cat("-- SE --\n"); print(fit_mm$std_errors)
cat("-- vcov --\n"); print(as.matrix(fit_mm$fixed_effect_vcov))
cat("-- theta --\n"); print(fit_mm$theta)
cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
cat("-- sigma --\n"); print(sigma(fit_mm))
cat("-- logLik --\n"); print(logLik(fit_mm))
cat("-- AIC --\n");    print(fit_mm$AIC)
cat("-- BIC --\n");    print(fit_mm$BIC)
cat("-- ranef (first 3 subjects) --\n"); print(head(ranef(fit_mm)$Subject, 3))
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("\n")

## ── 4. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  if (length(lme4_val) == 0 || length(mm_val) == 0) {
    cat(sprintf("  %-20s  SKIPPED (empty)\n", label))
    return(invisible(NULL))
  }
  if (length(lme4_val) != length(mm_val)) {
    cat(sprintf("  %-20s  LENGTH MISMATCH: lme4=%d mm=%d\n",
                label, length(lme4_val), length(mm_val)))
    return(invisible(NULL))
  }
  diffs <- abs(lme4_val - mm_val)
  max_d <- max(diffs)
  status <- if (max_d <= tol) "WITHIN-TOL" else "*** DIVERGED ***"
  cat(sprintf("  %-20s  max|diff|=%-14s  tol=%-10s  %s\n",
              label, sprintf(fmt, max_d), sprintf(fmt, tol), status))
  if (max_d > tol) {
    cat(sprintf("    lme4 : %s\n", paste(sprintf(fmt, lme4_val), collapse = "  ")))
    cat(sprintf("    mixeff: %s\n", paste(sprintf(fmt, mm_val),  collapse = "  ")))
  }
  invisible(max_d)
}

## fixef
compare("fixef[Intercept]", lme4::fixef(fit_lme4)["(Intercept)"],
        fixef(fit_mm)["(Intercept)"], tols$fixef)
compare("fixef[Days]",      lme4::fixef(fit_lme4)["Days"],
        fixef(fit_mm)["Days"], tols$fixef)

## SE
compare("SE[Intercept]",
        sqrt(diag(vcov(fit_lme4)))["(Intercept)"],
        fit_mm$std_errors["(Intercept)"], tols$fixef)
compare("SE[Days]",
        sqrt(diag(vcov(fit_lme4)))["Days"],
        fit_mm$std_errors["Days"], tols$fixef)

## vcov (flatten)
lme4_vcov_vec <- as.numeric(as.matrix(vcov(fit_lme4)))
mm_vcov_vec   <- as.numeric(as.matrix(fit_mm$fixed_effect_vcov))
compare("vcov (all elements)", lme4_vcov_vec, mm_vcov_vec, tols$fixef)

## theta (Cholesky factor elements)
lme4_theta <- lme4::getME(fit_lme4, "theta")
mm_theta   <- fit_mm$theta
compare("theta",    lme4_theta, mm_theta, tols$theta)

## sigma
compare("sigma",    sigma(fit_lme4), sigma(fit_mm), tols$sigma)

## logLik
compare("logLik",   as.numeric(logLik(fit_lme4)), fit_mm$logLik, tols$logLik)

## AIC / BIC
compare("AIC",      AIC(fit_lme4),  fit_mm$AIC, tols$logLik)
compare("BIC",      BIC(fit_lme4),  fit_mm$BIC, tols$logLik)

## ranef (vectorize all subjects)
lme4_ranef_subj <- as.matrix(lme4::ranef(fit_lme4)$Subject)
mm_ranef_subj   <- as.matrix(ranef(fit_mm)$Subject)
# align row order
common_subj <- intersect(rownames(lme4_ranef_subj), rownames(mm_ranef_subj))
if (length(common_subj) == nrow(lme4_ranef_subj)) {
  lme4_re_vec <- as.numeric(lme4_ranef_subj[common_subj, ])
  mm_re_vec   <- as.numeric(mm_ranef_subj[common_subj, ])
  compare("ranef (all subjects)", lme4_re_vec, mm_re_vec, tols$theta)
} else {
  cat("  ranef: subject sets differ — lme4:", nrow(lme4_ranef_subj),
      "  mixeff:", nrow(mm_ranef_subj), "\n")
}

## fitted values
compare("fitted (all obs)", as.numeric(fitted(fit_lme4)),
        as.numeric(fit_mm$fitted), tols$fixef)

## VarCorr variances  (Subject intercept variance, Subject slope variance, residual)
vc_lme4 <- lme4::VarCorr(fit_lme4)
vc_mm   <- VarCorr(fit_mm)

lme4_var_int  <- as.numeric(vc_lme4$Subject["(Intercept)", "(Intercept)"])
lme4_var_days <- as.numeric(vc_lme4$Subject["Days",        "Days"])
lme4_cor_id   <- as.numeric(attr(vc_lme4$Subject, "correlation")["(Intercept)", "Days"])
lme4_resid_sd <- as.numeric(attr(vc_lme4, "sc"))

mm_tbl        <- vc_mm$table
mm_var_int    <- mm_tbl$variance[mm_tbl$name == "(Intercept)"]
mm_var_days   <- mm_tbl$variance[mm_tbl$name == "Days"]
mm_resid_sd   <- vc_mm$residual_sd

compare("VarCorr var(Intercept)", sqrt(lme4_var_int),  mm_tbl$std_dev[mm_tbl$name == "(Intercept)"], tols$theta)
compare("VarCorr var(Days)",      sqrt(lme4_var_days), mm_tbl$std_dev[mm_tbl$name == "Days"],        tols$theta)
compare("VarCorr residual_sd",    lme4_resid_sd,       mm_resid_sd,                                  tols$sigma)

# correlation (Intercept, Days)
mm_cor_text <- mm_tbl$correlation[mm_tbl$name == "Days"]
cat(sprintf("  %-20s  lme4=%.4f  mixeff_text='%s'\n",
            "VarCorr cor(Int,Days)", lme4_cor_id, mm_cor_text))

cat("\n=== SPEED RATIO ===\n")
lme4_mean   <- t_lme4["elapsed"]   / NREPS
mixeff_mean <- t_mixeff["elapsed"] / NREPS
ratio <- lme4_mean / mixeff_mean
cat(sprintf("lme4  mean elapsed : %.4f s\n", lme4_mean))
cat(sprintf("mixeff mean elapsed: %.4f s\n", mixeff_mean))
cat(sprintf("ratio lme4/mixeff  : %.2fx  (%s)\n", ratio,
            ifelse(ratio > 1, "mixeff is FASTER", "mixeff is SLOWER")))

cat("\n=== SUMMARY ===\n")
cat("Tolerances: fixef/SE/vcov/fitted=1e-4, theta/ranef/VarCorr=1e-3, sigma=1e-4, logLik/AIC/BIC=1e-3\n")
