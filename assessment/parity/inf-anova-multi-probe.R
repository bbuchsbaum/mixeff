## Empirical parity probe: inf-anova-multi
## Dataset: sleepstudy
## Focus: anova(m1, m2) LRT — two nested LMMs compared via compare() / anova()
##        Compares ML-refit, Chisq statistic, Df, p-value between lme4 and mixeff.
## Also reports per-model fixef, SE, theta, sigma, logLik, AIC/BIC as context.

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

## ── 2. Define model pairs ────────────────────────────────────────────────────
## Pair A: random-intercept vs random-intercept+slope (fixed effects held equal)
## m1_A: Reaction ~ Days + (1 | Subject)        [smaller, RI only]
## m2_A: Reaction ~ Days + (Days | Subject)     [larger, RS model]
## This tests nested random-effect structures — classic LRT use-case.
##
## Pair B: same RE structure, differing fixed effects (intercept-only vs Days)
## m1_B: Reaction ~ 1 + (Days | Subject)        [smaller, intercept-only]
## m2_B: Reaction ~ Days + (Days | Subject)     [larger, Days fixed]

cat("=== MODEL PAIR A: RI vs RS (random effects differ) ===\n")
cat("m1_A: Reaction ~ Days + (1|Subject)  [random intercept only]\n")
cat("m2_A: Reaction ~ Days + (Days|Subject)  [random slope+intercept]\n\n")

## ── 3. Fit lme4 models ───────────────────────────────────────────────────────
cat("=== lme4 FITS ===\n")

## Pair A — lme4
fit_lme4_A1 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE)
fit_lme4_A2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)

cat("-- lme4 A1: Reaction ~ Days + (1|Subject) --\n")
cat("  logLik:", as.numeric(logLik(fit_lme4_A1)), "  df:", attr(logLik(fit_lme4_A1), "df"), "\n")
cat("  AIC:", AIC(fit_lme4_A1), "  BIC:", BIC(fit_lme4_A1), "\n")
cat("  fixef:", paste(round(lme4::fixef(fit_lme4_A1), 5), collapse = ", "), "\n")
cat("  theta:", paste(round(lme4::getME(fit_lme4_A1, "theta"), 5), collapse = ", "), "\n")
cat("  sigma:", round(sigma(fit_lme4_A1), 5), "\n\n")

cat("-- lme4 A2: Reaction ~ Days + (Days|Subject) --\n")
cat("  logLik:", as.numeric(logLik(fit_lme4_A2)), "  df:", attr(logLik(fit_lme4_A2), "df"), "\n")
cat("  AIC:", AIC(fit_lme4_A2), "  BIC:", BIC(fit_lme4_A2), "\n")
cat("  fixef:", paste(round(lme4::fixef(fit_lme4_A2), 5), collapse = ", "), "\n")
cat("  theta:", paste(round(lme4::getME(fit_lme4_A2, "theta"), 5), collapse = ", "), "\n")
cat("  sigma:", round(sigma(fit_lme4_A2), 5), "\n\n")

## lme4 anova (LRT)
cat("-- lme4 anova(m1_A, m2_A) [LRT] --\n")
lme4_anova_A <- anova(fit_lme4_A1, fit_lme4_A2, refit = FALSE)  # already ML
print(lme4_anova_A)
cat("\n")

## Extract lme4 LRT quantities for pair A
lme4_A_chisq <- lme4_anova_A$Chisq[2]
lme4_A_df    <- lme4_anova_A$`Chi Df`[2]
lme4_A_pval  <- lme4_anova_A$`Pr(>Chisq)`[2]
lme4_A_ll1   <- as.numeric(logLik(fit_lme4_A1))
lme4_A_ll2   <- as.numeric(logLik(fit_lme4_A2))

cat(sprintf("  lme4 LRT: Chisq=%.6f  Df=%d  p=%.8f\n",
            lme4_A_chisq, lme4_A_df, lme4_A_pval))
cat(sprintf("  lme4 logLik: m1=%.6f  m2=%.6f  diff=%.6f\n",
            lme4_A_ll1, lme4_A_ll2, lme4_A_ll2 - lme4_A_ll1))
cat("\n")

## Pair B — lme4
cat("=== MODEL PAIR B: fixed effects differ (intercept-only vs Days) ===\n")
cat("m1_B: Reaction ~ 1 + (Days|Subject)\n")
cat("m2_B: Reaction ~ Days + (Days|Subject)\n\n")

fit_lme4_B1 <- lmer(Reaction ~ 1 + (Days | Subject), data = sleepstudy, REML = FALSE)
fit_lme4_B2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)

cat("-- lme4 B1: Reaction ~ 1 + (Days|Subject) --\n")
cat("  logLik:", as.numeric(logLik(fit_lme4_B1)), "  df:", attr(logLik(fit_lme4_B1), "df"), "\n")
cat("  AIC:", AIC(fit_lme4_B1), "  BIC:", BIC(fit_lme4_B1), "\n")
cat("  fixef:", paste(round(lme4::fixef(fit_lme4_B1), 5), collapse = ", "), "\n\n")

cat("-- lme4 B2: Reaction ~ Days + (Days|Subject) [same as A2] --\n")
cat("  logLik:", as.numeric(logLik(fit_lme4_B2)), "  df:", attr(logLik(fit_lme4_B2), "df"), "\n\n")

cat("-- lme4 anova(m1_B, m2_B) [LRT] --\n")
lme4_anova_B <- anova(fit_lme4_B1, fit_lme4_B2, refit = FALSE)
print(lme4_anova_B)
cat("\n")

lme4_B_chisq <- lme4_anova_B$Chisq[2]
lme4_B_df    <- lme4_anova_B$`Chi Df`[2]
lme4_B_pval  <- lme4_anova_B$`Pr(>Chisq)`[2]

cat(sprintf("  lme4 LRT B: Chisq=%.6f  Df=%d  p=%.8f\n\n",
            lme4_B_chisq, lme4_B_df, lme4_B_pval))

## ── 4. Fit mixeff models ─────────────────────────────────────────────────────
cat("=== mixeff FITS ===\n")

fit_mm_A1 <- tryCatch(
  lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE,
      control = mm_control(verbose = -1L)),
  error = function(e) e
)
fit_mm_A2 <- tryCatch(
  lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE,
      control = mm_control(verbose = -1L)),
  error = function(e) e
)

if (inherits(fit_mm_A1, "condition")) {
  cat("!!! mixeff ERROR on A1 !!!\n")
  cat("class:", paste(class(fit_mm_A1), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm_A1), "\n")
} else {
  cat("-- mixeff A1: Reaction ~ Days + (1|Subject) --\n")
  cat("  fit_status:", fit_mm_A1$fit_status, "\n")
  cat("  logLik:", fit_mm_A1$logLik, "  df:", fit_mm_A1$dof, "\n")
  cat("  AIC:", fit_mm_A1$AIC, "  BIC:", fit_mm_A1$BIC, "\n")
  cat("  fixef:", paste(round(as.numeric(fixef(fit_mm_A1)), 5), collapse = ", "), "\n")
  cat("  theta:", paste(round(fit_mm_A1$theta, 5), collapse = ", "), "\n")
  cat("  sigma:", round(sigma(fit_mm_A1), 5), "\n\n")
}

if (inherits(fit_mm_A2, "condition")) {
  cat("!!! mixeff ERROR on A2 !!!\n")
  cat("class:", paste(class(fit_mm_A2), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm_A2), "\n")
} else {
  cat("-- mixeff A2: Reaction ~ Days + (Days|Subject) --\n")
  cat("  fit_status:", fit_mm_A2$fit_status, "\n")
  cat("  logLik:", fit_mm_A2$logLik, "  df:", fit_mm_A2$dof, "\n")
  cat("  AIC:", fit_mm_A2$AIC, "  BIC:", fit_mm_A2$BIC, "\n")
  cat("  fixef:", paste(round(as.numeric(fixef(fit_mm_A2)), 5), collapse = ", "), "\n")
  cat("  theta:", paste(round(fit_mm_A2$theta, 5), collapse = ", "), "\n")
  cat("  sigma:", round(sigma(fit_mm_A2), 5), "\n\n")
}

## mixeff compare / anova for pair A
cat("-- mixeff compare(A1, A2) [LRT] --\n")
mm_cmp_A <- tryCatch(
  compare(fit_mm_A1, fit_mm_A2),
  error = function(e) e
)
if (inherits(mm_cmp_A, "condition")) {
  cat("!!! mixeff compare() ERROR !!!\n")
  cat("class:", paste(class(mm_cmp_A), collapse = ", "), "\n")
  cat(conditionMessage(mm_cmp_A), "\n")
} else {
  print(mm_cmp_A)
  cat("\n")
}

cat("-- mixeff anova(A1, A2) [dispatches to compare()] --\n")
mm_anova_A <- tryCatch(
  anova(fit_mm_A1, fit_mm_A2),
  error = function(e) e
)
if (inherits(mm_anova_A, "condition")) {
  cat("!!! mixeff anova() ERROR !!!\n")
  cat("class:", paste(class(mm_anova_A), collapse = ", "), "\n")
  cat(conditionMessage(mm_anova_A), "\n\n")
} else {
  print(mm_anova_A)
  cat("\n")
}

## Pair B — mixeff
fit_mm_B1 <- tryCatch(
  lmm(Reaction ~ 1 + (Days | Subject), data = sleepstudy, REML = FALSE,
      control = mm_control(verbose = -1L)),
  error = function(e) e
)
fit_mm_B2 <- fit_mm_A2  # same model

cat("-- mixeff B1: Reaction ~ 1 + (Days|Subject) --\n")
if (inherits(fit_mm_B1, "condition")) {
  cat("!!! mixeff ERROR on B1 !!!\n")
  cat("class:", paste(class(fit_mm_B1), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm_B1), "\n\n")
} else {
  cat("  fit_status:", fit_mm_B1$fit_status, "\n")
  cat("  logLik:", fit_mm_B1$logLik, "  df:", fit_mm_B1$dof, "\n")
  cat("  AIC:", fit_mm_B1$AIC, "  BIC:", fit_mm_B1$BIC, "\n")
  cat("  fixef:", paste(round(as.numeric(fixef(fit_mm_B1)), 5), collapse = ", "), "\n\n")
}

cat("-- mixeff compare(B1, B2) [LRT] --\n")
if (!inherits(fit_mm_B1, "condition") && !inherits(fit_mm_B2, "condition")) {
  mm_cmp_B <- tryCatch(
    compare(fit_mm_B1, fit_mm_B2),
    error = function(e) e
  )
  if (inherits(mm_cmp_B, "condition")) {
    cat("!!! mixeff compare() B ERROR !!!\n")
    cat("class:", paste(class(mm_cmp_B), collapse = ", "), "\n")
    cat(conditionMessage(mm_cmp_B), "\n\n")
  } else {
    print(mm_cmp_B)
    cat("\n")
  }
} else {
  cat("  SKIPPED (upstream fit failed)\n\n")
  mm_cmp_B <- NULL
}

## ── 5. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4, chisq = 1e-3, pval = 1e-4)

cmp <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  if (length(lme4_val) == 0 || length(mm_val) == 0) {
    cat(sprintf("  %-30s  SKIPPED (empty)\n", label))
    return(invisible(NA))
  }
  if (length(lme4_val) != length(mm_val)) {
    cat(sprintf("  %-30s  LENGTH MISMATCH: lme4=%d mm=%d\n",
                label, length(lme4_val), length(mm_val)))
    return(invisible(NA))
  }
  diffs <- abs(lme4_val - mm_val)
  max_d <- max(diffs)
  status <- if (max_d <= tol) "WITHIN-TOL" else "*** DIVERGED ***"
  cat(sprintf("  %-30s  max|diff|=%-16s  tol=%-12s  %s\n",
              label, sprintf(fmt, max_d), sprintf(fmt, tol), status))
  if (max_d > tol) {
    cat(sprintf("    lme4 : %s\n", paste(sprintf(fmt, lme4_val), collapse = "  ")))
    cat(sprintf("    mixeff: %s\n", paste(sprintf(fmt, mm_val),  collapse = "  ")))
  }
  invisible(max_d)
}

cat("\n-- PAIR A: RI vs RS model --\n")

## per-model quantities for A1
if (!inherits(fit_mm_A1, "condition")) {
  cmp("A1 fixef[Intercept]", lme4::fixef(fit_lme4_A1)["(Intercept)"],
      fixef(fit_mm_A1)["(Intercept)"], tols$fixef)
  cmp("A1 fixef[Days]", lme4::fixef(fit_lme4_A1)["Days"],
      fixef(fit_mm_A1)["Days"], tols$fixef)
  cmp("A1 theta", lme4::getME(fit_lme4_A1, "theta"),
      fit_mm_A1$theta, tols$theta)
  cmp("A1 sigma", sigma(fit_lme4_A1), sigma(fit_mm_A1), tols$sigma)
  cmp("A1 logLik", as.numeric(logLik(fit_lme4_A1)), fit_mm_A1$logLik, tols$logLik)
  cmp("A1 AIC",    AIC(fit_lme4_A1),  fit_mm_A1$AIC, tols$logLik)
  cmp("A1 BIC",    BIC(fit_lme4_A1),  fit_mm_A1$BIC, tols$logLik)
}

if (!inherits(fit_mm_A2, "condition")) {
  cmp("A2 fixef[Intercept]", lme4::fixef(fit_lme4_A2)["(Intercept)"],
      fixef(fit_mm_A2)["(Intercept)"], tols$fixef)
  cmp("A2 fixef[Days]", lme4::fixef(fit_lme4_A2)["Days"],
      fixef(fit_mm_A2)["Days"], tols$fixef)
  cmp("A2 theta", lme4::getME(fit_lme4_A2, "theta"),
      fit_mm_A2$theta, tols$theta)
  cmp("A2 sigma", sigma(fit_lme4_A2), sigma(fit_mm_A2), tols$sigma)
  cmp("A2 logLik", as.numeric(logLik(fit_lme4_A2)), fit_mm_A2$logLik, tols$logLik)
  cmp("A2 AIC",    AIC(fit_lme4_A2),  fit_mm_A2$AIC, tols$logLik)
  cmp("A2 BIC",    BIC(fit_lme4_A2),  fit_mm_A2$BIC, tols$logLik)
}

## LRT comparison for pair A
cat("\n-- PAIR A LRT quantities --\n")
if (!inherits(mm_cmp_A, "condition")) {
  tbl_A <- mm_cmp_A$table
  # Second row is the comparison row (first is reference)
  mm_A_row <- tbl_A[tbl_A$lrt_available == TRUE, , drop = FALSE]
  if (nrow(mm_A_row) > 0) {
    mm_A_chisq <- mm_A_row$LRT[1]
    mm_A_df    <- mm_A_row$delta_df[1]
    mm_A_pval  <- mm_A_row$p_value[1]
    mm_A_ll1   <- tbl_A$logLik[1]
    mm_A_ll2   <- tbl_A$logLik[2]

    cat(sprintf("  lme4:   Chisq=%.6f  Df=%.0f  p=%.8f\n",
                lme4_A_chisq, lme4_A_df, lme4_A_pval))
    cat(sprintf("  mixeff: Chisq=%.6f  Df=%.0f  p=%.8f\n",
                mm_A_chisq, mm_A_df, mm_A_pval))

    cmp("Pair A LRT Chisq", lme4_A_chisq, mm_A_chisq, tols$chisq)
    cmp("Pair A LRT Df",    lme4_A_df,    mm_A_df,    0.5)   # integer; any diff is wrong
    cmp("Pair A LRT p-val", lme4_A_pval,  mm_A_pval,  tols$pval)
    cmp("Pair A logLik m1", lme4_A_ll1, mm_A_ll1, tols$logLik)
    cmp("Pair A logLik m2", lme4_A_ll2, mm_A_ll2, tols$logLik)
  } else {
    cat("  mixeff compare table has no lrt_available row for pair A\n")
    cat("  table:\n"); print(tbl_A)
  }
} else {
  cat("  mixeff compare() failed — LRT quantities unavailable\n")
}

## Pair B
cat("\n-- PAIR B: fixed effects differ --\n")
if (!inherits(fit_mm_B1, "condition")) {
  cmp("B1 fixef[Intercept]", lme4::fixef(fit_lme4_B1)["(Intercept)"],
      fixef(fit_mm_B1)["(Intercept)"], tols$fixef)
  cmp("B1 logLik", as.numeric(logLik(fit_lme4_B1)), fit_mm_B1$logLik, tols$logLik)
}

if (!is.null(mm_cmp_B) && !inherits(mm_cmp_B, "condition")) {
  tbl_B <- mm_cmp_B$table
  mm_B_row <- tbl_B[tbl_B$lrt_available == TRUE, , drop = FALSE]
  if (nrow(mm_B_row) > 0) {
    mm_B_chisq <- mm_B_row$LRT[1]
    mm_B_df    <- mm_B_row$delta_df[1]
    mm_B_pval  <- mm_B_row$p_value[1]

    cat(sprintf("  lme4:   Chisq=%.6f  Df=%.0f  p=%.8f\n",
                lme4_B_chisq, lme4_B_df, lme4_B_pval))
    cat(sprintf("  mixeff: Chisq=%.6f  Df=%.0f  p=%.8f\n",
                mm_B_chisq, mm_B_df, mm_B_pval))

    cmp("Pair B LRT Chisq", lme4_B_chisq, mm_B_chisq, tols$chisq)
    cmp("Pair B LRT Df",    lme4_B_df,    mm_B_df,    0.5)
    cmp("Pair B LRT p-val", lme4_B_pval,  mm_B_pval,  tols$pval)
  } else {
    cat("  mixeff compare table has no lrt_available row for pair B\n")
    cat("  table:\n"); print(tbl_B)
  }
} else if (!is.null(mm_cmp_B)) {
  cat("  mixeff compare() B failed\n")
}

## ── 6. REML-start compare (auto-refit path) ───────────────────────────────────
cat("\n=== PAIR A WITH REML FITS (auto-refit test) ===\n")
cat("Fitting REML models and using compare() with refit_for_comparison='auto'\n")

fit_lme4_A1_reml <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
fit_lme4_A2_reml <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)

cat("lme4 anova() with REML fits (lme4 auto-refits to ML internally):\n")
lme4_reml_anova <- tryCatch(
  anova(fit_lme4_A1_reml, fit_lme4_A2_reml),
  error = function(e) e
)
if (inherits(lme4_reml_anova, "condition")) {
  cat("  lme4 anova() error:", conditionMessage(lme4_reml_anova), "\n")
} else {
  print(lme4_reml_anova)
}

fit_mm_A1_reml <- tryCatch(
  lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE,
      control = mm_control(verbose = -1L)),
  error = function(e) e
)
fit_mm_A2_reml <- tryCatch(
  lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
      control = mm_control(verbose = -1L)),
  error = function(e) e
)

if (!inherits(fit_mm_A1_reml, "condition") && !inherits(fit_mm_A2_reml, "condition")) {
  cat("\nmixeff compare() with REML fits, refit_for_comparison='auto':\n")
  mm_reml_cmp <- tryCatch(
    compare(fit_mm_A1_reml, fit_mm_A2_reml, refit_for_comparison = "auto"),
    error = function(e) e
  )
  if (inherits(mm_reml_cmp, "condition")) {
    cat("!!! mixeff compare() REML-auto ERROR !!!\n")
    cat("class:", paste(class(mm_reml_cmp), collapse = ", "), "\n")
    cat(conditionMessage(mm_reml_cmp), "\n")
  } else {
    print(mm_reml_cmp)
    tbl_reml <- mm_reml_cmp$table
    mm_reml_row <- tbl_reml[tbl_reml$lrt_available == TRUE, , drop = FALSE]
    if (nrow(mm_reml_row) > 0) {
      lme4_reml_chisq <- lme4_reml_anova$Chisq[2]
      lme4_reml_pval  <- lme4_reml_anova$`Pr(>Chisq)`[2]
      cat(sprintf("\n  lme4 (post-refit): Chisq=%.6f  p=%.8f\n",
                  lme4_reml_chisq, lme4_reml_pval))
      cat(sprintf("  mixeff (auto-refit): Chisq=%.6f  p=%.8f\n",
                  mm_reml_row$LRT[1], mm_reml_row$p_value[1]))
      cmp("REML-auto Chisq", lme4_reml_chisq, mm_reml_row$LRT[1], tols$chisq)
      cmp("REML-auto p-val",  lme4_reml_pval,  mm_reml_row$p_value[1], tols$pval)
    }
  }
} else {
  cat("mixeff REML fits failed — skipping REML auto-refit path\n")
}

## ── 7. Speed ─────────────────────────────────────────────────────────────────
cat("\n=== SPEED COMPARISON ===\n")
NREPS <- 10L

t_lme4_compare <- system.time({
  for (i in seq_len(NREPS)) {
    f1 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE)
    f2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
    anova(f1, f2, refit = FALSE)
  }
})

t_mm_compare <- system.time({
  for (i in seq_len(NREPS)) {
    f1 <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE,
              control = mm_control(verbose = -1L))
    f2 <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE,
              control = mm_control(verbose = -1L))
    compare(f1, f2)
  }
})

lme4_mean   <- t_lme4_compare["elapsed"] / NREPS
mm_mean     <- t_mm_compare["elapsed"]   / NREPS
ratio       <- lme4_mean / mm_mean
cat(sprintf("lme4  mean elapsed (fit + anova, %d reps): %.4f s\n", NREPS, lme4_mean))
cat(sprintf("mixeff mean elapsed (fit + compare):       %.4f s\n", mm_mean))
cat(sprintf("ratio lme4/mixeff: %.2fx  (%s)\n", ratio,
            ifelse(ratio > 1, "mixeff is FASTER", "mixeff is SLOWER")))

cat("\n=== SUMMARY ===\n")
cat("Tolerances: fixef/SE/vcov/fitted=1e-4, theta/ranef=1e-3, sigma=1e-4, logLik/AIC/BIC=1e-3\n")
cat("LRT-specific: Chisq=1e-3, p-value=1e-4\n")
