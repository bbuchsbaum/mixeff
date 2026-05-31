## Empirical parity probe: lmm-sleep-ri
## Dataset: sleepstudy  Formula: Reaction ~ Days + (1|Subject)  REML = TRUE
## Compares lme4/lmerTest vs mixeff on all relevant quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:",    as.character(packageVersion("lme4")),    "\n")
cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff version:",  as.character(packageVersion("mixeff")),  "\n\n")

## ── 1. Data ──────────────────────────────────────────────────────────────────
data(sleepstudy, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(sleepstudy), "  ncol:", ncol(sleepstudy), "\n")
print(head(sleepstudy))
cat("\n")

## ── 2. Fit lme4 ──────────────────────────────────────────────────────────────
cat("=== lme4 FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(Reaction ~ Days + (1 | Subject),
                   data = sleepstudy, REML = TRUE)
})
cat("lme4 wall-clock (seconds):", t_lme4["elapsed"], "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov --\n"); print(as.matrix(vcov(fit_lme4)))
cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))
cat("-- sigma --\n"); print(stats::sigma(fit_lme4))
cat("-- logLik --\n"); print(logLik(fit_lme4))
cat("-- AIC --\n"); print(AIC(fit_lme4))
cat("-- BIC --\n"); print(BIC(fit_lme4))
cat("-- ranef (first 6) --\n"); print(head(lme4::ranef(fit_lme4)$Subject))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
cat("-- convergence --\n")
cat("No convergence warning =",
    length(fit_lme4@optinfo$conv$lme4$messages) == 0, "\n\n")

## ── 3. Fit mixeff ─────────────────────────────────────────────────────────────
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    lmm(Reaction ~ Days + (1 | Subject),
        data = sleepstudy, REML = TRUE,
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

cat("-- fixef --\n");         print(fixef(fit_mm))
cat("-- SE --\n");            print(fit_mm$std_errors)
cat("-- vcov --\n");          print(as.matrix(fit_mm$fixed_effect_vcov))
cat("-- VarCorr --\n");       print(VarCorr(fit_mm))
cat("-- sigma --\n");         print(sigma(fit_mm))
cat("-- logLik --\n");        print(logLik(fit_mm))
cat("-- AIC --\n");           print(fit_mm$AIC)
cat("-- BIC --\n");           print(fit_mm$BIC)
cat("-- ranef (first 6) --\n"); print(head(ranef(fit_mm)$Subject))
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
  cat(sprintf("%-34s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse = ", "),
              paste(sprintf(fmt, mm_val),   collapse = ", "),
              diff, tol, status))
  invisible(diff)
}

## fixef: (Intercept) and Days
fe_lme4 <- lme4::fixef(fit_lme4)
fe_mm   <- fixef(fit_mm)

compare("fixef (Intercept)",
        fe_lme4[["(Intercept)"]],
        fe_mm[["(Intercept)"]],
        tols$fixef)

compare("fixef Days",
        fe_lme4[["Days"]],
        fe_mm[["Days"]],
        tols$fixef)

## SE
se_lme4 <- sqrt(diag(vcov(fit_lme4)))  # vcov is stats:: generic, no conflict
se_mm   <- fit_mm$std_errors

compare("SE (Intercept)",
        se_lme4[["(Intercept)"]],
        se_mm[["(Intercept)"]],
        tols$fixef)

compare("SE Days",
        se_lme4[["Days"]],
        se_mm[["Days"]],
        tols$fixef)

## vcov diagonal
compare("vcov[1,1]",
        as.numeric(vcov(fit_lme4)[1, 1]),
        as.numeric(fit_mm$fixed_effect_vcov[1, 1]),
        tols$fixef^2)

compare("vcov[2,2]",
        as.numeric(vcov(fit_lme4)[2, 2]),
        as.numeric(fit_mm$fixed_effect_vcov[2, 2]),
        tols$fixef^2)

## theta (random-intercept only: one Cholesky factor)
theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
cat(sprintf("%-34s lme4=%s  mm=%s\n", "theta (raw)",
            paste(sprintf("%.8f", theta_lme4), collapse = ", "),
            paste(sprintf("%.8f", theta_mm),   collapse = ", ")))
compare("theta", theta_lme4, theta_mm, tols$theta)

## sigma
compare("sigma",
        stats::sigma(fit_lme4),
        sigma(fit_mm),
        tols$sigma)

## VarCorr: Subject RE variance
vc_lme4 <- as.numeric(attr(lme4::VarCorr(fit_lme4)$Subject, "stddev"))^2
vc_mm_obj <- VarCorr(fit_mm)
vc_mm <- tryCatch({
  # mm_varcorr is a list with $table (data.frame) and $residual_sd
  if (inherits(vc_mm_obj, "mm_varcorr") && is.data.frame(vc_mm_obj$table)) {
    tbl <- vc_mm_obj$table
    as.numeric(tbl[tbl$group == "Subject", "variance"])
  } else if (is.data.frame(vc_mm_obj)) {
    as.numeric(vc_mm_obj[vc_mm_obj$grp == "Subject", "vcov"])
  } else {
    NA_real_
  }
}, error = function(e) {
  cat("VarCorr extraction error:", conditionMessage(e), "\n")
  NA_real_
})

cat(sprintf("%-34s lme4=%.8f  mm=%s\n", "VarCorr Subject var",
            vc_lme4,
            ifelse(is.na(vc_mm), "NA", sprintf("%.8f", vc_mm))))
if (!is.na(vc_mm)) compare("VarCorr Subject var", vc_lme4, vc_mm, tols$theta)

## logLik
compare("logLik",
        as.numeric(stats::logLik(fit_lme4)),
        as.numeric(logLik(fit_mm)),
        tols$logLik)

## AIC
compare("AIC",
        stats::AIC(fit_lme4),
        fit_mm$AIC,
        tols$logLik * 2)

## BIC
compare("BIC",
        stats::BIC(fit_lme4),
        fit_mm$BIC,
        tols$logLik * 2)

## fitted values
fitted_lme4 <- as.numeric(stats::fitted(fit_lme4))
fitted_mm   <- as.numeric(fit_mm$fitted)
compare("fitted max abs diff",
        fitted_lme4,
        fitted_mm[seq_along(fitted_lme4)],
        tols$fixef)

## ranef (sorted by value for order-invariant comparison)
re_lme4 <- as.numeric(lme4::ranef(fit_lme4)$Subject[, 1])
re_mm_df <- ranef(fit_mm)$Subject
re_mm    <- if (!is.null(re_mm_df)) as.numeric(re_mm_df[, 1]) else NA_real_
cat(sprintf("\n%-34s\n", "ranef Subject (sorted):"))
cat("  lme4:", sprintf("%.5f", sort(re_lme4)), "\n")
cat("  mm:  ", sprintf("%.5f", sort(re_mm)),   "\n")
compare("ranef Subject max abs diff",
        sort(re_lme4), sort(re_mm), tols$fixef)

## ── 5. Speed comparison ───────────────────────────────────────────────────────
cat("\n=== SPEED COMPARISON ===\n")
NREPS <- 5
t_lme4_rep <- system.time(
  for (i in seq_len(NREPS))
    lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
)
t_mm_rep <- system.time(
  for (i in seq_len(NREPS))
    lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE,
        control = mm_control(verbose = -1L))
)
lme4_per  <- t_lme4_rep["elapsed"] / NREPS
mm_per    <- t_mm_rep["elapsed"]   / NREPS
cat(sprintf("lme4  mean/fit: %.4f s  (over %d reps)\n", lme4_per, NREPS))
cat(sprintf("mm    mean/fit: %.4f s  (over %d reps)\n", mm_per,   NREPS))
cat(sprintf("ratio (mm/lme4): %.2fx\n", mm_per / lme4_per))

cat("\n=== SUMMARY TABLE ===\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
cat("All WITHIN-TOL = parity achieved.\n")
