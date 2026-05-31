## Empirical parity probe: lmm-dyestuff-ri
## Dataset: Dyestuff  Formula: Yield ~ 1 + (1|Batch)  REML = TRUE
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
data(Dyestuff, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(Dyestuff), "  ncol:", ncol(Dyestuff), "\n")
print(head(Dyestuff))
cat("\n")

## ── 2. Fit lme4 ──────────────────────────────────────────────────────────────
cat("=== lme4 FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(Yield ~ 1 + (1 | Batch), data = Dyestuff, REML = TRUE)
})
cat("lme4 wall-clock (seconds):", t_lme4["elapsed"], "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov --\n"); print(as.matrix(vcov(fit_lme4)))
cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))
cat("-- sigma --\n"); print(sigma(fit_lme4))
cat("-- logLik --\n"); print(logLik(fit_lme4))
cat("-- AIC --\n");    print(AIC(fit_lme4))
cat("-- BIC --\n");    print(BIC(fit_lme4))
cat("-- ranef --\n"); print(lme4::ranef(fit_lme4))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
cat("-- convergence --\n"); cat("No convergence warning =",
  length(fit_lme4@optinfo$conv$lme4$messages) == 0, "\n\n")

## ── 3. Fit mixeff ─────────────────────────────────────────────────────────────
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    lmm(Yield ~ 1 + (1 | Batch), data = Dyestuff, REML = TRUE,
        control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff wall-clock (seconds):", t_mixeff["elapsed"], "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR !!!\n")
  cat(conditionMessage(fit_mm), "\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  quit(status = 0)
}

cat("fit_status:", fit_mm$fit_status, "\n\n")

cat("-- fixef --\n"); print(fixef(fit_mm))
cat("-- SE --\n"); print(fit_mm$std_errors)
cat("-- vcov --\n"); print(as.matrix(fit_mm$fixed_effect_vcov))
cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
cat("-- sigma --\n"); print(sigma(fit_mm))
cat("-- logLik --\n"); print(logLik(fit_mm))
cat("-- AIC --\n");    print(fit_mm$AIC)
cat("-- BIC --\n");    print(fit_mm$BIC)
cat("-- ranef --\n"); print(ranef(fit_mm))
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("\n")

## ── 4. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  diff     <- max(abs(lme4_val - mm_val))
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-30s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse=", "),
              paste(sprintf(fmt, mm_val),   collapse=", "),
              diff, tol, status))
  invisible(diff)
}

## fixef
compare("fixef (Intercept)",
        lme4::fixef(fit_lme4)[["(Intercept)"]],
        fixef(fit_mm)[["(Intercept)"]],
        tols$fixef)

## SE
compare("SE (Intercept)",
        sqrt(diag(vcov(fit_lme4)))[["(Intercept)"]],
        fit_mm$std_errors[["(Intercept)"]],
        tols$fixef)

## vcov diagonal
compare("vcov[1,1]",
        as.numeric(vcov(fit_lme4)[1,1]),
        as.numeric(fit_mm$fixed_effect_vcov[1,1]),
        tols$fixef^2)

## theta (Cholesky factor of RE covariance)
theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
cat(sprintf("%-30s lme4=%s  mm=%s\n", "theta",
            paste(sprintf("%.8f", theta_lme4), collapse=", "),
            paste(sprintf("%.8f", theta_mm),   collapse=", ")))
compare("theta", theta_lme4, theta_mm, tols$theta)

## sigma
compare("sigma",
        sigma(fit_lme4),
        sigma(fit_mm),
        tols$sigma)

## VarCorr: RE variance
vc_lme4 <- as.numeric(attr(lme4::VarCorr(fit_lme4)$Batch, "stddev"))^2
vc_mm_df <- VarCorr(fit_mm)
# extract numeric variance for Batch from mixeff VarCorr object
vc_mm <- tryCatch({
  # mm VarCorr may be a data.frame-like or a list
  if (is.data.frame(vc_mm_df)) {
    as.numeric(vc_mm_df[vc_mm_df$grp == "Batch", "vcov"])
  } else if (is.list(vc_mm_df)) {
    v <- vc_mm_df[["Batch"]]
    if (is.numeric(v)) v
    else if (is.matrix(v)) v[1,1]
    else NA_real_
  } else {
    NA_real_
  }
}, error = function(e) { cat("VarCorr extraction error:", conditionMessage(e), "\n"); NA_real_ })

cat(sprintf("%-30s lme4=%.8f  mm=%s\n", "VarCorr Batch var",
            vc_lme4,
            ifelse(is.na(vc_mm), "NA", sprintf("%.8f", vc_mm))))
if (!is.na(vc_mm)) compare("VarCorr Batch var", vc_lme4, vc_mm, tols$theta)

## logLik
compare("logLik",
        as.numeric(logLik(fit_lme4)),
        as.numeric(logLik(fit_mm)),
        tols$logLik)

## AIC
compare("AIC",
        AIC(fit_lme4),
        fit_mm$AIC,
        tols$logLik * 2)   # AIC = -2*logLik + 2k

## BIC
compare("BIC",
        BIC(fit_lme4),
        fit_mm$BIC,
        tols$logLik * 2)

## fitted values
fitted_lme4 <- fitted(fit_lme4)
fitted_mm   <- fit_mm$fitted
compare("fitted max abs diff",
        fitted_lme4,
        fitted_mm[seq_along(fitted_lme4)],
        tols$fixef)

## ranef
re_lme4 <- as.numeric(lme4::ranef(fit_lme4)$Batch[,1])
re_mm_df <- ranef(fit_mm)$Batch
re_mm    <- if (!is.null(re_mm_df)) as.numeric(re_mm_df[,1]) else NA_real_
cat(sprintf("\n%-30s\n", "ranef Batch (sorted):"))
cat("  lme4:", sprintf("%.6f", sort(re_lme4)), "\n")
cat("  mm:  ", sprintf("%.6f", sort(re_mm)),   "\n")
compare("ranef Batch max abs diff",
        sort(re_lme4), sort(re_mm), tols$fixef)

## speed ratio
cat(sprintf("\n%-30s lme4=%.4fs  mm=%.4fs  ratio(mm/lme4)=%.2fx\n",
            "wall-clock elapsed",
            t_lme4["elapsed"], t_mixeff["elapsed"],
            t_mixeff["elapsed"] / t_lme4["elapsed"]))

cat("\n=== SUMMARY TABLE ===\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
cat("All WITHIN-TOL = parity achieved.\n")
