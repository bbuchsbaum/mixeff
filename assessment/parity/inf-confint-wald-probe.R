## inf-confint-wald-probe.R
## Empirical parity probe: confint(method="wald") for sleepstudy
## Formula: Reaction ~ Days + (Days|Subject)
## Compares lme4 vs mixeff on Wald CIs plus supporting quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

## Explicit namespace wrappers to avoid masking conflicts
lme4_fixef  <- function(x) lme4::fixef(x)
lme4_ranef  <- function(x) lme4::ranef(x)
lme4_sigma  <- function(x) stats::sigma(x)
mm_fixef    <- function(x) mixeff::fixef(x)
mm_ranef    <- function(x) mixeff::ranef(x)
mm_sigma    <- function(x) stats::sigma(x)

cat("=== MODEL FITTING ===\n\n")

data(sleepstudy, package = "lme4")
form <- Reaction ~ Days + (Days | Subject)

## --- lme4 fit ---
cat("--- lme4 fit ---\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(form, data = sleepstudy, REML = TRUE)
})
cat(sprintf("lme4 fit time: %.4f s\n", t_lme4["elapsed"]))
print(summary(fit_lme4))

## --- mixeff fit ---
cat("\n--- mixeff fit ---\n")
t_mm <- system.time({
  fit_mm <- tryCatch(
    lmm(form, data = sleepstudy, REML = TRUE),
    error = function(e) { cat("ERROR fitting mixeff model:", conditionMessage(e), "\n"); NULL }
  )
})
cat(sprintf("mixeff fit time: %.4f s\n", t_mm["elapsed"]))
if (!is.null(fit_mm)) print(fit_mm)

cat("\n=== SUPPORTING QUANTITIES ===\n\n")

if (!is.null(fit_mm)) {

  ## fixed effects
  fe_lme4 <- lme4_fixef(fit_lme4)
  fe_mm   <- tryCatch(mm_fixef(fit_mm), error = function(e) { cat("fixef error:", conditionMessage(e), "\n"); NULL })
  cat("--- fixef ---\n")
  cat("lme4:   "); print(fe_lme4)
  cat("mixeff: "); print(fe_mm)
  if (!is.null(fe_mm)) {
    diff_fe <- abs(fe_lme4 - fe_mm)
    cat("abs diff: "); print(diff_fe)
    cat(sprintf("max abs diff fixef: %.2e  (tol 1e-4: %s)\n",
                max(diff_fe), if (max(diff_fe) < 1e-4) "PASS" else "FAIL"))
  }

  ## SE / vcov
  se_lme4 <- sqrt(diag(vcov(fit_lme4)))
  se_mm   <- tryCatch(sqrt(diag(vcov(fit_mm))), error = function(e) { cat("vcov error:", conditionMessage(e), "\n"); NULL })
  cat("\n--- SE (sqrt diag vcov) ---\n")
  cat("lme4:   "); print(se_lme4)
  cat("mixeff: "); print(se_mm)
  if (!is.null(se_mm)) {
    diff_se <- abs(se_lme4 - se_mm)
    cat("abs diff: "); print(diff_se)
    cat(sprintf("max abs diff SE: %.2e  (tol 1e-4: %s)\n",
                max(diff_se), if (max(diff_se) < 1e-4) "PASS" else "FAIL"))
  }

  ## theta / VarCorr
  th_lme4 <- lme4::getME(fit_lme4, "theta")
  th_mm   <- tryCatch(fit_mm$theta, error = function(e) { cat("theta error:", conditionMessage(e), "\n"); NULL })
  cat("\n--- theta ---\n")
  cat("lme4 theta:   "); print(th_lme4)
  cat("mixeff theta: "); print(th_mm)
  if (!is.null(th_mm) && length(th_lme4) == length(th_mm)) {
    diff_th <- abs(th_lme4 - th_mm)
    cat("abs diff: "); print(diff_th)
    cat(sprintf("max abs diff theta: %.2e  (tol 1e-3: %s)\n",
                max(diff_th), if (max(diff_th) < 1e-3) "PASS" else "FAIL"))
  }

  ## sigma
  sig_lme4 <- lme4_sigma(fit_lme4)
  sig_mm   <- tryCatch(mm_sigma(fit_mm), error = function(e) { cat("sigma error:", conditionMessage(e), "\n"); NULL })
  cat("\n--- sigma ---\n")
  cat(sprintf("lme4 sigma:   %.6f\n", sig_lme4))
  if (!is.null(sig_mm)) {
    cat(sprintf("mixeff sigma: %.6f\n", sig_mm))
    diff_sig <- abs(sig_lme4 - sig_mm)
    cat(sprintf("abs diff sigma: %.2e  (tol 1e-4: %s)\n",
                diff_sig, if (diff_sig < 1e-4) "PASS" else "FAIL"))
  }

  ## logLik
  ll_lme4 <- logLik(fit_lme4)
  ll_mm   <- tryCatch(logLik(fit_mm), error = function(e) { cat("logLik error:", conditionMessage(e), "\n"); NULL })
  cat("\n--- logLik ---\n")
  cat(sprintf("lme4 logLik:   %.6f\n", as.numeric(ll_lme4)))
  if (!is.null(ll_mm)) {
    cat(sprintf("mixeff logLik: %.6f\n", as.numeric(ll_mm)))
    diff_ll <- abs(as.numeric(ll_lme4) - as.numeric(ll_mm))
    cat(sprintf("abs diff logLik: %.2e  (tol 1e-3: %s)\n",
                diff_ll, if (diff_ll < 1e-3) "PASS" else "FAIL"))
  }

  ## AIC / BIC
  aic_lme4 <- AIC(fit_lme4); bic_lme4 <- BIC(fit_lme4)
  aic_mm   <- tryCatch(AIC(fit_mm), error = function(e) { cat("AIC error:", conditionMessage(e), "\n"); NULL })
  bic_mm   <- tryCatch(BIC(fit_mm), error = function(e) { cat("BIC error:", conditionMessage(e), "\n"); NULL })
  cat("\n--- AIC/BIC ---\n")
  cat(sprintf("lme4   AIC=%.4f  BIC=%.4f\n", aic_lme4, bic_lme4))
  if (!is.null(aic_mm) && !is.null(bic_mm)) {
    cat(sprintf("mixeff AIC=%.4f  BIC=%.4f\n", aic_mm, bic_mm))
    cat(sprintf("abs diff AIC: %.2e   BIC: %.2e\n",
                abs(aic_lme4 - aic_mm), abs(bic_lme4 - bic_mm)))
  }

  ## ranef
  re_lme4 <- lme4_ranef(fit_lme4)$Subject
  re_mm   <- tryCatch(mm_ranef(fit_mm)$Subject, error = function(e) { cat("ranef error:", conditionMessage(e), "\n"); NULL })
  cat("\n--- ranef (first 5 subjects) ---\n")
  cat("lme4:\n"); print(head(re_lme4, 5))
  if (!is.null(re_mm)) {
    cat("mixeff:\n"); print(head(re_mm, 5))
    # align by rownames
    common <- intersect(rownames(re_lme4), rownames(re_mm))
    diff_re <- abs(re_lme4[common, ] - re_mm[common, ])
    cat(sprintf("max abs diff ranef: %.2e\n", max(diff_re)))
  }

  ## fitted
  fv_lme4 <- fitted(fit_lme4)
  fv_mm   <- tryCatch(fitted(fit_mm), error = function(e) { cat("fitted error:", conditionMessage(e), "\n"); NULL })
  cat("\n--- fitted values ---\n")
  if (!is.null(fv_mm)) {
    diff_fv <- abs(fv_lme4 - fv_mm)
    cat(sprintf("max abs diff fitted: %.2e\n", max(diff_fv)))
  }

  ## convergence
  conv_lme4 <- length(fit_lme4@optinfo$conv$lme4) == 0
  conv_mm   <- tryCatch(isTRUE(fit_mm$converged), error = function(e) NA)
  cat("\n--- convergence ---\n")
  cat(sprintf("lme4 converged:   %s\n", conv_lme4))
  cat(sprintf("mixeff converged: %s\n", conv_mm))
}

cat("\n=== WALD CONFINT (primary cell target) ===\n\n")

## lme4 Wald confint
cat("--- lme4 confint(method='Wald') ---\n")
ci_lme4 <- tryCatch(
  confint(fit_lme4, method = "Wald"),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
print(ci_lme4)

## mixeff Wald confint
cat("\n--- mixeff confint(method='wald') ---\n")
ci_mm <- NULL
if (!is.null(fit_mm)) {
  ci_mm <- tryCatch(
    confint(fit_mm, method = "wald"),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  print(ci_mm)
}

## Comparison: only fixed-effect rows are directly comparable
## lme4 returns rows for .sig01/.sig02/.sigma (VarCorr params) + fixed effects
## mixeff Wald only covers fixed effects by design (per inference.R line ~887-895)
cat("\n--- Fixed-effect Wald CI comparison ---\n")
if (!is.null(ci_lme4) && !is.null(ci_mm)) {
  fe_names <- names(lme4_fixef(fit_lme4))
  # lme4 rows for fixed effects have same names
  ci_lme4_fe <- ci_lme4[fe_names, , drop = FALSE]
  # mixeff rows
  ci_mm_fe   <- as.matrix(ci_mm)[fe_names, , drop = FALSE]

  cat("lme4 fixed-effect CIs:\n")
  print(ci_lme4_fe)
  cat("\nmixeff fixed-effect CIs:\n")
  print(ci_mm_fe)

  # align column direction (both should be lower, upper)
  diff_ci <- abs(ci_lme4_fe - ci_mm_fe)
  cat("\nabs diff per bound:\n")
  print(diff_ci)
  cat(sprintf("max abs diff Wald CI (fixed effects): %.2e  (tol 1e-4: %s)\n",
              max(diff_ci), if (max(diff_ci) < 1e-4) "PASS" else "FAIL"))
}

## Does mixeff cover the VarCorr parameters in Wald CIs?
cat("\n--- Does mixeff Wald CI cover VarCorr/sigma parameters? ---\n")
if (!is.null(ci_mm)) {
  mm_rows <- rownames(as.matrix(ci_mm))
  lme4_vc_rows <- rownames(ci_lme4)[!rownames(ci_lme4) %in% names(lme4_fixef(fit_lme4))]
  cat(sprintf("lme4 Wald CI rows (non-fixef): %s\n", paste(lme4_vc_rows, collapse = ", ")))
  cat(sprintf("mixeff Wald CI rows:           %s\n", paste(mm_rows, collapse = ", ")))
  missing_vc <- setdiff(lme4_vc_rows, mm_rows)
  if (length(missing_vc)) {
    cat(sprintf("MISSING from mixeff Wald CI: %s\n", paste(missing_vc, collapse = ", ")))
  } else {
    cat("All lme4 Wald CI rows present in mixeff.\n")
  }
}

## parm subset selection
cat("\n--- parm subset selection (parm='Days') ---\n")
if (!is.null(fit_mm)) {
  ci_mm_days <- tryCatch(
    confint(fit_mm, parm = "Days", method = "wald"),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  print(ci_mm_days)
}

## level != 0.95
cat("\n--- level=0.90 ---\n")
ci_lme4_90 <- tryCatch(confint(fit_lme4, method = "Wald", level = 0.90), error = function(e) NULL)
ci_mm_90   <- if (!is.null(fit_mm)) tryCatch(confint(fit_mm, method = "wald", level = 0.90), error = function(e) NULL) else NULL
if (!is.null(ci_lme4_90) && !is.null(ci_mm_90)) {
  fe_names <- names(lme4_fixef(fit_lme4))
  diff_90 <- abs(ci_lme4_90[fe_names, ] - as.matrix(ci_mm_90)[fe_names, ])
  cat(sprintf("max abs diff Wald CI level=0.90: %.2e  (tol 1e-4: %s)\n",
              max(diff_90), if (max(diff_90) < 1e-4) "PASS" else "FAIL"))
}

## bad parm
cat("\n--- bad parm (should error gracefully) ---\n")
if (!is.null(fit_mm)) {
  tryCatch(
    confint(fit_mm, parm = "NONEXISTENT", method = "wald"),
    error = function(e) cat("Error (expected):", conditionMessage(e), "\n"),
    warning = function(w) cat("Warning:", conditionMessage(w), "\n")
  )
}

## bad level
cat("\n--- bad level (should error gracefully) ---\n")
if (!is.null(fit_mm)) {
  tryCatch(
    confint(fit_mm, level = 1.5, method = "wald"),
    error = function(e) cat("Error (expected):", conditionMessage(e), "\n"),
    warning = function(w) cat("Warning:", conditionMessage(w), "\n")
  )
}

## asymptotic synonym
cat("\n--- method='asymptotic' synonym ---\n")
if (!is.null(fit_mm)) {
  ci_asym <- tryCatch(
    confint(fit_mm, method = "asymptotic"),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ci_asym) && !is.null(ci_mm)) {
    diff_syn <- max(abs(as.matrix(ci_mm) - as.matrix(ci_asym)))
    cat(sprintf("wald vs asymptotic synonym max diff: %.2e (should be 0)\n", diff_syn))
  }
}

cat("\n=== SUMMARY TABLE ===\n")
cat("Quantity                   | lme4       | mixeff     | abs diff   | tol    | status\n")
cat("---------------------------|------------|------------|------------|--------|-------\n")

summarize_row <- function(label, v1, v2, tol) {
  if (is.null(v2) || any(is.na(v2))) {
    cat(sprintf("%-27s| %-10s | %-10s | %-10s | %-6s | %s\n",
                label, format(v1, digits=5), "ERROR/NA", "N/A", format(tol, scientific=TRUE), "FAIL"))
    return(invisible(NULL))
  }
  d <- max(abs(v1 - v2))
  status <- if (d < tol) "PASS" else "FAIL"
  cat(sprintf("%-27s| %-10s | %-10s | %-10s | %-6s | %s\n",
              label,
              format(round(v1[1], 5), nsmall=5),
              format(round(v2[1], 5), nsmall=5),
              format(d, scientific=TRUE, digits=3),
              format(tol, scientific=TRUE),
              status))
}

if (!is.null(fit_mm)) {
  fe_lme4 <- lme4_fixef(fit_lme4)
  fe_mm   <- tryCatch(mm_fixef(fit_mm), error = function(e) NULL)
  se_lme4 <- sqrt(diag(vcov(fit_lme4)))
  se_mm   <- tryCatch(sqrt(diag(vcov(fit_mm))), error = function(e) NULL)

  summarize_row("fixef[Intercept]",    fe_lme4["(Intercept)"],    if (!is.null(fe_mm)) fe_mm["(Intercept)"] else NULL, 1e-4)
  summarize_row("fixef[Days]",         fe_lme4["Days"],           if (!is.null(fe_mm)) fe_mm["Days"] else NULL, 1e-4)
  summarize_row("SE[Intercept]",       se_lme4["(Intercept)"],    if (!is.null(se_mm)) se_mm["(Intercept)"] else NULL, 1e-4)
  summarize_row("SE[Days]",            se_lme4["Days"],           if (!is.null(se_mm)) se_mm["Days"] else NULL, 1e-4)
  summarize_row("sigma",               lme4_sigma(fit_lme4),      tryCatch(mm_sigma(fit_mm), error=function(e) NULL), 1e-4)
  summarize_row("logLik",              as.numeric(logLik(fit_lme4)), tryCatch(as.numeric(logLik(fit_mm)), error=function(e) NULL), 1e-3)

  if (!is.null(ci_lme4) && !is.null(ci_mm)) {
    fe_names <- names(lme4_fixef(fit_lme4))
    ci_lme4_fe <- ci_lme4[fe_names, , drop = FALSE]
    ci_mm_fe   <- as.matrix(ci_mm)[fe_names, , drop = FALSE]
    summarize_row("Wald CI lower[Intercept]", ci_lme4_fe["(Intercept)", 1], ci_mm_fe["(Intercept)", 1], 1e-4)
    summarize_row("Wald CI upper[Intercept]", ci_lme4_fe["(Intercept)", 2], ci_mm_fe["(Intercept)", 2], 1e-4)
    summarize_row("Wald CI lower[Days]",      ci_lme4_fe["Days", 1],        ci_mm_fe["Days", 1], 1e-4)
    summarize_row("Wald CI upper[Days]",      ci_lme4_fe["Days", 2],        ci_mm_fe["Days", 2], 1e-4)
  }
}

cat("\n=== COVERAGE GAP NOTE ===\n")
cat("lme4 confint(method='Wald') returns rows for VarCorr parameters\n")
cat("(.sig01, .sig02, .sig03, .sigma) in addition to fixed effects.\n")
cat("mixeff confint(method='wald') covers ONLY fixed effects (by design).\n")
cat("This is a coverage gap vs lme4 Wald behaviour.\n")
cat("Attr status on mixeff output: 'not_certified_by_rust_inference_contract'\n")
