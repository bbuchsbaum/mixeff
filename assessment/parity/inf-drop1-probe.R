## Empirical parity probe: inf-drop1
## Dataset: sleepstudy  Formula: Reaction ~ Days + (1|Subject)
## Focus: drop1 parity — lme4/lmerTest vs mixeff
##
## Key distinction: lme4/lmerTest drop1() uses Satterthwaite F-tests.
## mixeff drop1() uses asymptotic LRT (chi-sq). These are conceptually
## different but both test significance of fixed-effect terms.
##
## We compare:
##   - Base fit: fixef, SE/vcov, theta/VarCorr, sigma, logLik, AIC/BIC,
##               ranef, fitted, convergence
##   - drop1 inference: each engine's output, then explicit comparison
##     of what is comparable (reduced logLik, df, p-value scale)
##   - test_effect Satterthwaite vs lmerTest anova() for Days

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

## ── 2. Fit lme4 (REML=TRUE) ──────────────────────────────────────────────────
cat("=== lme4 FIT (REML=TRUE) ===\n")
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

## ── 3. lme4/lmerTest drop1 ───────────────────────────────────────────────────
cat("=== lme4/lmerTest DROP1 (via lmerTest; uses Satterthwaite F-test) ===\n")
cat("-- drop1(fit_lme4, test='Chisq') --\n")
d1_lme4_chisq_raw <- tryCatch(
  capture.output(print(drop1(fit_lme4, test = "Chisq"))),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(d1_lme4_chisq_raw)) cat(paste(d1_lme4_chisq_raw, collapse = "\n"), "\n")
d1_lme4_chisq <- tryCatch(drop1(fit_lme4, test = "Chisq"), error = function(e) NULL)
cat("\n")

cat("-- drop1(fit_lme4, test='none') --\n")
d1_lme4_none_raw <- tryCatch(
  capture.output(print(drop1(fit_lme4, test = "none"))),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(d1_lme4_none_raw)) cat(paste(d1_lme4_none_raw, collapse = "\n"), "\n")
d1_lme4_none <- tryCatch(drop1(fit_lme4, test = "none"), error = function(e) NULL)
cat("\n")

## lmerTest Satterthwaite anova
cat("-- lmerTest::anova(fit_lme4, type='III') --\n")
lmert_anova <- tryCatch(
  lmerTest::anova(fit_lme4, type = "III"),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(lmert_anova)) print(lmert_anova)
cat("\n")

## Also use stats::anova for LRT comparison
cat("=== lme4 REML=FALSE fit for LRT-based comparison ===\n")
fit_lme4_ml <- lmer(Reaction ~ Days + (1 | Subject),
                    data = sleepstudy, REML = FALSE)
fit_lme4_ml_reduced <- lmer(Reaction ~ 1 + (1 | Subject),
                             data = sleepstudy, REML = FALSE)
cat("-- stats::anova(full_ml, reduced_ml) -- LRT --\n")
lme4_lrt_raw <- tryCatch(
  anova(fit_lme4_ml_reduced, fit_lme4_ml),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(lme4_lrt_raw)) print(lme4_lrt_raw)
cat("\n")

cat("-- drop1(fit_lme4_ml, test='Chisq') via lmerTest (F-test) --\n")
d1_lme4_ml_raw <- tryCatch(
  capture.output(print(drop1(fit_lme4_ml, test = "Chisq"))),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(d1_lme4_ml_raw)) cat(paste(d1_lme4_ml_raw, collapse = "\n"), "\n")
d1_lme4_ml <- tryCatch(drop1(fit_lme4_ml, test = "Chisq"), error = function(e) NULL)
cat("\n")

## ── 4. Fit mixeff (REML=TRUE) ─────────────────────────────────────────────────
cat("=== mixeff FIT (REML=TRUE) ===\n")
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
  cat("!!! mixeff FIT ERROR !!!\n")
  cat(conditionMessage(fit_mm), "\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  quit(status = 0)
}

cat("fit_status:", fit_mm$fit_status, "\n\n")
cat("-- fixef --\n");         print(fixef(fit_mm))
cat("-- SE --\n");            print(fit_mm$std_errors)
cat("-- sigma --\n");         print(sigma(fit_mm))
cat("-- logLik --\n");        print(logLik(fit_mm))
cat("-- AIC --\n");           print(fit_mm$AIC)
cat("-- BIC --\n");           print(fit_mm$BIC)
cat("-- ranef (first 6) --\n"); print(head(ranef(fit_mm)$Subject))
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("\n")

## ── 5. mixeff drop1 ───────────────────────────────────────────────────────────
cat("=== mixeff DROP1 ===\n")
cat("-- drop1(fit_mm, test='Chisq') --\n")
t_drop1_mm <- system.time({
  d1_mm_chisq <- tryCatch(
    drop1(fit_mm, test = "Chisq"),
    error = function(e) e
  )
})
cat("drop1 wall-clock:", t_drop1_mm["elapsed"], "s\n\n")
if (inherits(d1_mm_chisq, "condition")) {
  cat("!!! mixeff drop1(Chisq) ERROR !!!\n")
  cat("class:", paste(class(d1_mm_chisq), collapse = ", "), "\n")
  cat(conditionMessage(d1_mm_chisq), "\n")
} else {
  cat("class of result:", paste(class(d1_mm_chisq), collapse = ", "), "\n")
  print(d1_mm_chisq)
  cat("\n-- Full table --\n")
  print(d1_mm_chisq$table, digits = 8)
}
cat("\n")

cat("-- drop1(fit_mm, test='none') --\n")
d1_mm_none <- tryCatch(
  drop1(fit_mm, test = "none"),
  error = function(e) e
)
if (inherits(d1_mm_none, "condition")) {
  cat("!!! mixeff drop1(none) ERROR !!!\n")
  cat(conditionMessage(d1_mm_none), "\n")
} else {
  print(d1_mm_none)
}
cat("\n")

## test_effect Satterthwaite
cat("-- test_effect(fit_mm, 'Days', method='satterthwaite') --\n")
te_sat <- tryCatch(
  test_effect(fit_mm, "Days", method = "satterthwaite"),
  error = function(e) e
)
if (inherits(te_sat, "condition")) {
  cat("!!! test_effect Satterthwaite ERROR !!!\n")
  cat(conditionMessage(te_sat), "\n")
} else {
  cat("-- test_effect table (relevant columns) --\n")
  cols_of_interest <- c("term", "num_df", "den_df", "statistic",
                         "statistic_name", "p_value", "method", "status")
  cols_present <- intersect(cols_of_interest, names(te_sat$table))
  print(te_sat$table[, cols_present, drop = FALSE], digits = 8)
}
cat("\n")

cat("-- test_effect(fit_mm, 'Days', method='kenward_roger') --\n")
te_kr <- tryCatch(
  test_effect(fit_mm, "Days", method = "kenward_roger"),
  error = function(e) e
)
if (inherits(te_kr, "condition")) {
  cat("!!! test_effect KR ERROR !!!\n")
  cat(conditionMessage(te_kr), "\n")
} else {
  cols_present <- intersect(cols_of_interest, names(te_kr$table))
  print(te_kr$table[, cols_present, drop = FALSE], digits = 8)
}
cat("\n")

## ── 6. Base fit numerical comparison ─────────────────────────────────────────
cat("=== BASE FIT NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  if (length(lme4_val) != length(mm_val)) {
    cat(sprintf("%-38s  LENGTH MISMATCH lme4=%d mm=%d\n",
                label, length(lme4_val), length(mm_val)))
    return(invisible(NA))
  }
  diff   <- max(abs(lme4_val - mm_val))
  status <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  if (length(lme4_val) <= 2) {
    cat(sprintf("%-38s lme4=%-24s mm=%-24s maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
                label,
                paste(sprintf(fmt, lme4_val), collapse=", "),
                paste(sprintf(fmt, mm_val),   collapse=", "),
                diff, tol, status))
  } else {
    cat(sprintf("%-38s maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
                label, diff, tol, status))
  }
  invisible(diff)
}

fe_lme4 <- lme4::fixef(fit_lme4)
fe_mm   <- fixef(fit_mm)
compare("fixef (Intercept)", fe_lme4[["(Intercept)"]], fe_mm[["(Intercept)"]], tols$fixef)
compare("fixef Days",        fe_lme4[["Days"]],        fe_mm[["Days"]],        tols$fixef)

se_lme4 <- sqrt(diag(vcov(fit_lme4)))
se_mm   <- fit_mm$std_errors
compare("SE (Intercept)", se_lme4[["(Intercept)"]], se_mm[["(Intercept)"]], tols$fixef)
compare("SE Days",        se_lme4[["Days"]],        se_mm[["Days"]],        tols$fixef)

theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
compare("theta", theta_lme4, theta_mm, tols$theta)

compare("sigma", stats::sigma(fit_lme4), sigma(fit_mm), tols$sigma)

compare("logLik (REML)",
        as.numeric(stats::logLik(fit_lme4)),
        as.numeric(logLik(fit_mm)),
        tols$logLik)

compare("AIC", stats::AIC(fit_lme4), fit_mm$AIC, tols$logLik * 2)
compare("BIC", stats::BIC(fit_lme4), fit_mm$BIC, tols$logLik * 2)

fitted_lme4 <- as.numeric(stats::fitted(fit_lme4))
fitted_mm   <- as.numeric(fit_mm$fitted)
compare("fitted max abs diff", fitted_lme4, fitted_mm[seq_along(fitted_lme4)], tols$fixef)

re_lme4 <- as.numeric(lme4::ranef(fit_lme4)$Subject[, 1])
re_mm_df <- ranef(fit_mm)$Subject
re_mm    <- if (!is.null(re_mm_df)) as.numeric(re_mm_df[, 1]) else NA_real_
compare("ranef Subject max abs diff", sort(re_lme4), sort(re_mm), tols$fixef)

cat("\n")

## ── 7. drop1 inference comparison ────────────────────────────────────────────
cat("=== DROP1 INFERENCE COMPARISON ===\n")

## 7a. mixeff asymptotic LRT vs lme4 explicit anova(ml_full, ml_reduced) LRT
cat("--- 7a. Asymptotic LRT: lme4 anova(ML models) vs mixeff drop1(Chisq) ---\n")
if (!is.null(lme4_lrt_raw) && !inherits(d1_mm_chisq, "condition")) {
  lme4_lrt_tbl <- as.data.frame(lme4_lrt_raw)
  ## Row 2 is the full model (alternative)
  lme4_ll_full    <- lme4_lrt_tbl[2, "logLik"]
  lme4_ll_reduced <- lme4_lrt_tbl[1, "logLik"]
  lme4_Chisq      <- lme4_lrt_tbl[2, "Chisq"]
  lme4_Df         <- lme4_lrt_tbl[2, "Df"]
  lme4_p          <- lme4_lrt_tbl[2, "Pr(>Chisq)"]

  cat(sprintf("lme4 anova LRT: ll_full=%.6f  ll_reduced=%.6f  Chisq=%.6f  Df=%g  p=%.3e\n",
              lme4_ll_full, lme4_ll_reduced, lme4_Chisq, lme4_Df, lme4_p))

  mm_row  <- d1_mm_chisq$table[d1_mm_chisq$table$dropped == "Days", , drop = FALSE]
  mm_lrt  <- mm_row$LRT
  mm_df   <- mm_row$df
  mm_p    <- mm_row$p_value
  mm_ll_r <- mm_row$logLik

  cat(sprintf("mixeff drop1:   ll_reduced=%.6f                LRT=%.6f   df=%g  p=%.3e\n",
              mm_ll_r, mm_lrt, mm_df, mm_p))
  cat("\n")

  ## Compare reduced-model logLik
  ll_diff <- abs(lme4_ll_reduced - mm_ll_r)
  ll_status <- if (ll_diff <= tols$logLik) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-40s lme4=%.8f  mm=%.8f  diff=%.3e  tol=%.0e  [%s]\n",
              "reduced logLik (drop Days)",
              lme4_ll_reduced, mm_ll_r, ll_diff, tols$logLik, ll_status))

  ## Compare LRT statistic
  lrt_diff <- abs(lme4_Chisq - mm_lrt)
  lrt_tol  <- 0.01
  lrt_status <- if (lrt_diff <= lrt_tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-40s lme4=%.8f  mm=%.8f  diff=%.3e  tol=%.2f  [%s]\n",
              "LRT statistic (drop Days)",
              lme4_Chisq, mm_lrt, lrt_diff, lrt_tol, lrt_status))

  ## Compare Df
  df_match <- isTRUE(lme4_Df == mm_df)
  cat(sprintf("%-40s lme4=%g  mm=%g  [%s]\n",
              "LRT Df (drop Days)",
              lme4_Df, mm_df, if (df_match) "MATCH" else "MISMATCH"))

  ## Compare p-values
  p_diff <- abs(lme4_p - mm_p)
  p_tol  <- 1e-20  # both should be astronomically small; compare on log scale
  cat(sprintf("%-40s lme4=%.3e  mm=%.3e  diff=%.3e\n",
              "p-value (drop Days)", lme4_p, mm_p, p_diff))
  cat(sprintf("  Note: lme4 p based on ML anova(); mixeff p based on asymptotic LRT;\n"))
  cat(sprintf("  LRT statistics differ because of REML vs ML refit handling.\n"))
  cat(sprintf("  lme4 LRT stat = 2*(ll_full-ll_reduced) = 2*(%.6f - %.6f) = %.6f\n",
              lme4_ll_full, lme4_ll_reduced, 2*(lme4_ll_full - lme4_ll_reduced)))
  cat(sprintf("  mixeff LRT stat = %.6f  (computed from ML refit)\n", mm_lrt))

} else {
  cat("  Cannot compare: one or both engines failed.\n")
  if (is.null(lme4_lrt_raw))          cat("  lme4 anova(ML): FAILED\n")
  if (inherits(d1_mm_chisq, "condition")) {
    cat("  mixeff drop1(Chisq): FAILED —", conditionMessage(d1_mm_chisq), "\n")
  }
}

cat("\n--- 7b. lme4/lmerTest drop1 (REML, Satterthwaite F) ---\n")
cat("lme4/lmerTest uses Satterthwaite F-test (not LRT) for REML drop1.\n")
cat("mixeff test_effect(method='satterthwaite') is the direct analog.\n\n")

if (!is.null(d1_lme4_chisq) && !inherits(te_sat, "condition")) {
  ## extract Days row from lmerTest drop1
  lmert_row <- d1_lme4_chisq["Days", , drop = FALSE]
  cat("lmerTest drop1 Days row:\n")
  print(lmert_row)
  cat("\nmixeff test_effect Satterthwaite:\n")
  cols_present2 <- intersect(c("term","num_df","den_df","statistic","statistic_name","p_value","method"),
                              names(te_sat$table))
  print(te_sat$table[, cols_present2, drop = FALSE], digits = 8)

  lmert_F   <- as.numeric(lmert_row[["F value"]])
  lmert_df2 <- as.numeric(lmert_row[["DenDF"]])
  lmert_p   <- as.numeric(lmert_row[["Pr(>F)"]])

  mm_F  <- as.numeric(te_sat$table[["statistic"]])
  mm_df2 <- as.numeric(te_sat$table[["den_df"]])
  mm_p  <- as.numeric(te_sat$table[["p_value"]])

  cat("\n")

  ## Note: lmerTest F vs mixeff t: for single df, F = t^2
  cat(sprintf("lmerTest F(1, %.4f) = %.6f\n", lmert_df2, lmert_F))
  cat(sprintf("mixeff t(%.4f)      = %.6f  => F equivalent = %.6f\n",
              mm_df2, mm_F, mm_F^2))

  F_equiv_diff <- abs(lmert_F - mm_F^2)
  cat(sprintf("%-40s lmerTest=%.6f  mm(t^2)=%.6f  diff=%.3e  tol=0.01  [%s]\n",
              "F-stat equiv (Days)",
              lmert_F, mm_F^2, F_equiv_diff,
              if (F_equiv_diff <= 0.01) "WITHIN-TOL" else "EXCEEDS-TOL"))

  df2_diff <- abs(lmert_df2 - mm_df2)
  cat(sprintf("%-40s lmerTest=%.4f  mm=%.4f  diff=%.3e  tol=0.10  [%s]\n",
              "Satterthwaite DenDF (Days)",
              lmert_df2, mm_df2, df2_diff,
              if (df2_diff <= 0.1) "WITHIN-TOL" else "EXCEEDS-TOL"))

  p_diff2 <- abs(lmert_p - mm_p)
  cat(sprintf("%-40s lmerTest=%.4e  mm=%.4e  diff=%.3e  tol=1e-10  [%s]\n",
              "Satterthwaite p-value (Days)",
              lmert_p, mm_p, p_diff2,
              if (p_diff2 <= 1e-10) "WITHIN-TOL" else "EXCEEDS-TOL"))
} else {
  if (is.null(d1_lme4_chisq)) cat("  lmerTest drop1: FAILED\n")
  if (inherits(te_sat, "condition")) {
    cat("  mixeff test_effect Satterthwaite: FAILED —", conditionMessage(te_sat), "\n")
  }
}

## ── 8. Speed comparison ───────────────────────────────────────────────────────
cat("\n=== SPEED COMPARISON ===\n")
NREPS <- 5

t_lme4_rep <- system.time(
  for (i in seq_len(NREPS)) {
    f <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE)
    drop1(f, test = "Chisq")
  }
)
t_mm_rep <- system.time(
  for (i in seq_len(NREPS)) {
    f <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE,
             control = mm_control(verbose = -1L))
    drop1(f, test = "Chisq")
  }
)
lme4_per <- t_lme4_rep["elapsed"] / NREPS
mm_per   <- t_mm_rep["elapsed"]   / NREPS
cat(sprintf("lme4  mean/fit+drop1: %.4f s  (over %d reps)\n", lme4_per, NREPS))
cat(sprintf("mm    mean/fit+drop1: %.4f s  (over %d reps)\n", mm_per,   NREPS))
cat(sprintf("ratio (mm/lme4): %.2fx\n\n", mm_per / lme4_per))

cat("=== END OF PROBE ===\n")
