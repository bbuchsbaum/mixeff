## Parity probe: inf-ranef-condvar
## Dataset: sleepstudy
## Formula: Reaction ~ Days + (Days|Subject)
## Focus: ranef(condVar=TRUE) point estimates + conditional SDs vs lme4
## Written for cell "inf-ranef-condvar" in the parity assessment grid.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

# Helper: use %||% from mixeff namespace if not available
`%||%` <- function(a, b) if (!is.null(a)) a else b

cat("=== inf-ranef-condvar parity probe ===\n\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

data("sleepstudy", package = "lme4")

## ─── 1. Fit both models ────────────────────────────────────────────────────
cat("--- Fitting models ---\n")

t_lme4 <- system.time(
  fit_lme4 <- lmer(Reaction ~ Days + (Days | Subject),
                   data = sleepstudy, REML = TRUE)
)

t_mixeff <- system.time(
  fit_mm <- lmm(Reaction ~ Days + (Days | Subject),
                data = sleepstudy, REML = TRUE,
                control = mm_control(verbose = -1))
)

cat("lme4 wall time:   ", t_lme4["elapsed"], "s\n")
cat("mixeff wall time: ", t_mixeff["elapsed"], "s\n")
cat("Speed ratio (lme4/mixeff): ", round(t_lme4["elapsed"] / t_mixeff["elapsed"], 2), "\n\n")

## ─── 2. Convergence / refusal ─────────────────────────────────────────────
cat("--- Convergence / refusal status ---\n")
cat("lme4 converged: TRUE (no error)\n")
cat("mixeff converged:", !is.null(fit_mm), "(no error = TRUE)\n\n")

## ─── 3. Fixed effects comparison ──────────────────────────────────────────
cat("--- Fixed effects ---\n")
fe_lme4   <- lme4::fixef(fit_lme4)
fe_mm     <- fixef(fit_mm)           # mixeff::fixef
fe_diff   <- abs(fe_lme4 - fe_mm[names(fe_lme4)])
cat("lme4    fixef:", round(fe_lme4, 6), "\n")
cat("mixeff  fixef:", round(fe_mm[names(fe_lme4)], 6), "\n")
cat("abs diff:     ", round(fe_diff, 8), "\n")
cat("max abs diff (fixef): ", max(fe_diff), " [tol 1e-4]\n\n")

## ─── 4. SE / vcov comparison ──────────────────────────────────────────────
cat("--- SE / vcov ---\n")
se_lme4  <- sqrt(diag(as.matrix(lme4:::vcov.merMod(fit_lme4))))
se_mm    <- sqrt(diag(vcov(fit_mm)))
se_diff  <- abs(se_lme4 - se_mm[names(se_lme4)])
cat("lme4   SE:", round(se_lme4, 6), "\n")
cat("mixeff SE:", round(se_mm[names(se_lme4)], 6), "\n")
cat("abs diff: ", round(se_diff, 8), "\n")
cat("max abs diff (SE): ", max(se_diff), " [tol 1e-4]\n\n")

## ─── 5. theta / VarCorr ───────────────────────────────────────────────────
cat("--- VarCorr / theta ---\n")
vc_lme4 <- lme4::VarCorr(fit_lme4)
vc_mm   <- VarCorr(fit_mm)           # mixeff::VarCorr

# Extract theta (Cholesky factor elements) from lme4
theta_lme4 <- lme4::getME(fit_lme4, "theta")
cat("lme4 VarCorr (Subject):\n")
print(vc_lme4)
cat("\nmixeff VarCorr table:\n")
print(vc_mm$table)
cat("\nlme4 theta:", round(theta_lme4, 6), "\n")
sig_lme4 <- stats::sigma(fit_lme4)
sig_mm   <- sigma(fit_mm)            # mixeff::sigma
cat("lme4 sigma:", round(sig_lme4, 6), "\n")
cat("mixeff sigma:", round(sig_mm, 6), "\n")
cat("abs diff sigma:", abs(sig_lme4 - sig_mm), " [tol 1e-4]\n\n")

# Compare std devs from VarCorr
vc_lme4_sd <- attr(vc_lme4$Subject, "stddev")
mm_rows <- vc_mm$table[vc_mm$table$group == "Subject", ]

cat("lme4 Subject std devs  (Intercept, Days):", round(vc_lme4_sd, 6), "\n")
cat("mixeff Subject std devs (Intercept, Days):", round(mm_rows$std_dev, 6), "\n")
sd_diff <- abs(vc_lme4_sd - mm_rows$std_dev)
cat("abs diff std devs:", round(sd_diff, 8), "\n")
cat("max abs diff (theta/SD): ", max(sd_diff), " [tol 1e-3]\n\n")

## ─── 6. logLik / AIC / BIC ────────────────────────────────────────────────
cat("--- logLik / AIC / BIC ---\n")
ll_lme4  <- as.numeric(stats::logLik(fit_lme4))
ll_mm    <- as.numeric(logLik(fit_mm))
aic_lme4 <- stats::AIC(fit_lme4)
aic_mm   <- AIC(fit_mm)
bic_lme4 <- stats::BIC(fit_lme4)
bic_mm   <- BIC(fit_mm)

cat("lme4   logLik:", round(ll_lme4, 6), "\n")
cat("mixeff logLik:", round(ll_mm, 6), "\n")
cat("abs diff logLik:", abs(ll_lme4 - ll_mm), " [tol 1e-3]\n\n")
cat("lme4   AIC:", round(aic_lme4, 4), "\n")
cat("mixeff AIC:", round(aic_mm, 4), "\n")
cat("abs diff AIC:", abs(aic_lme4 - aic_mm), "\n\n")
cat("lme4   BIC:", round(bic_lme4, 4), "\n")
cat("mixeff BIC:", round(bic_mm, 4), "\n")
cat("abs diff BIC:", abs(bic_lme4 - bic_mm), "\n\n")

## ─── 7. ranef point estimates ─────────────────────────────────────────────
cat("--- ranef point estimates (condVar=FALSE) ---\n")
re_lme4 <- lme4::ranef(fit_lme4)$Subject
re_mm   <- ranef(fit_mm)$Subject     # mixeff::ranef

# Align rows by subject name
subj_order <- rownames(re_lme4)
re_mm_aligned <- re_mm[subj_order, , drop = FALSE]

re_diff_int  <- abs(re_lme4[["(Intercept)"]] - re_mm_aligned[["(Intercept)"]])
re_diff_days <- abs(re_lme4[["Days"]] - re_mm_aligned[["Days"]])

cat("lme4 ranef (first 6):\n")
print(round(head(re_lme4), 4))
cat("mixeff ranef (first 6):\n")
print(round(head(re_mm_aligned), 4))
cat("\nmax abs diff ranef (Intercept): ", max(re_diff_int), " [tol 1e-4]\n")
cat("max abs diff ranef (Days):       ", max(re_diff_days), " [tol 1e-4]\n\n")

## ─── 8. fitted values ─────────────────────────────────────────────────────
cat("--- Fitted values ---\n")
fv_lme4 <- stats::fitted(fit_lme4)
fv_mm   <- stats::fitted(fit_mm)
fv_diff <- abs(fv_lme4 - fv_mm)
cat("max abs diff fitted: ", max(fv_diff), " [tol 1e-4]\n\n")

## ─── 9. ranef(condVar=TRUE) — THE PRIMARY CELL FOCUS ─────────────────────
cat("=== PRIMARY FOCUS: ranef(condVar=TRUE) ===\n\n")

# lme4 conditional variances — use explicit namespace
re_lme4_cv <- lme4::ranef(fit_lme4, condVar = TRUE)
pv_lme4 <- attr(re_lme4_cv$Subject, "postVar")

cat("lme4 postVar array dim:", dim(pv_lme4), "\n")
# lme4 postVar has no dimnames on the 3rd dimension — use first slice
cat("lme4 postVar (first subject slice):\n")
print(round(pv_lme4[, , 1L], 6))
# Record subject ordering from ranef frame (lme4 keeps same order as ranef rows)
subj_names_lme4_order <- rownames(lme4::ranef(fit_lme4)$Subject)
cat("lme4 subject order:", subj_names_lme4_order, "\n\n")

# mixeff conditional variances
cat("Calling ranef(fit_mm, condVar=TRUE)...\n")
mm_condvar_result <- tryCatch(
  ranef(fit_mm, condVar = TRUE),   # mixeff::ranef
  error = function(e) e
)

common_subj <- NULL   # will be set in success branch

if (inherits(mm_condvar_result, "error")) {
  cat("ERROR from mixeff ranef(condVar=TRUE):\n")
  cat(conditionMessage(mm_condvar_result), "\n\n")
  cat("RESULT: mixeff ERRORS on condVar=TRUE\n")
} else {
  re_mm_cv <- mm_condvar_result

  # Check for unavailable-reason attribute (typed refusal)
  unavail_reason      <- attr(re_mm_cv, "mm_unavailable_reason")
  unavail_reason_subj <- attr(re_mm_cv$Subject, "mm_unavailable_reason")

  if (!is.null(unavail_reason) || !is.null(unavail_reason_subj)) {
    cat("mixeff TYPED REFUSAL for condVar=TRUE:\n")
    cat("  mm_unavailable_reason (list level):   ", unavail_reason %||% "NULL", "\n")
    cat("  mm_unavailable_reason (Subject level):", unavail_reason_subj %||% "NULL", "\n")
    cond_var_error <- attr(re_mm_cv$Subject, "mm_cond_var_error")
    cat("  mm_cond_var_error:", cond_var_error %||% "NULL", "\n\n")

    pv_mm <- attr(re_mm_cv$Subject, "postVar")
    cat("postVar from mixeff (typed refusal path):\n")
    if (is.null(pv_mm)) {
      cat("  postVar is NULL\n\n")
    } else {
      cat("  postVar dim:", dim(pv_mm), "\n")
      cat("  All NA?", all(is.na(pv_mm)), "\n")
      cat("  First slice:\n")
      print(pv_mm[, , 1])
      cat("\n")
    }
    cat("RESULT: mixeff gives typed refusal (NA postVar) rather than fabricating values.\n")
    cat("        This is a CAPABILITY GAP vs lme4 (lme4 returns real conditional variances).\n")
  } else {
    # Success path — compare arrays
    pv_mm <- attr(re_mm_cv$Subject, "postVar")

    cat("mixeff postVar array dim:", dim(pv_mm), "\n")

    if (is.null(pv_mm)) {
      cat("ERROR: postVar is NULL on success path — unexpected\n\n")
    } else {
      # lme4 postVar has no dimnames on dim 3; mixeff uses named subject levels.
      # The lme4 ranef rows and postVar slices share the same positional order.
      # Re-index lme4 postVar slices by the subject names from the ranef frame.
      subj_names_mm <- dimnames(pv_mm)[[3]]
      cat("mixeff subject order:", subj_names_mm, "\n")
      cat("lme4 subject order (from ranef frame):", subj_names_lme4_order, "\n\n")

      # Assign dimnames to lme4 postVar so we can align by name
      dimnames(pv_lme4) <- list(
        colnames(lme4::ranef(fit_lme4)$Subject),
        colnames(lme4::ranef(fit_lme4)$Subject),
        subj_names_lme4_order
      )

      # Subject 308 slice
      cat("lme4   postVar (Subject 308, slice):\n")
      print(round(pv_lme4[, , "308"], 6))
      cat("mixeff postVar (Subject 308, slice):\n")
      print(round(pv_mm[, , "308"], 6))
      cat("\n")

      # Check finiteness
      cat("All finite (mixeff postVar):", all(is.finite(pv_mm)), "\n")
      cat("Any NA:", any(is.na(pv_mm)), "\n\n")

      # Symmetry check
      sym_ok <- all(sapply(seq_len(dim(pv_mm)[[3]]), function(i) {
        sl <- pv_mm[, , i]
        max(abs(sl - t(sl))) < 1e-9
      }))
      cat("All slices symmetric:", sym_ok, "\n")

      # PSD check
      psd_ok <- all(sapply(seq_len(dim(pv_mm)[[3]]), function(i) {
        all(diag(pv_mm[, , i]) >= 0)
      }))
      cat("All diagonal entries >= 0:", psd_ok, "\n\n")

      # Align and compare using common subject names
      common_subj <- intersect(subj_names_lme4_order, subj_names_mm)
      pv_lme4_aligned <- pv_lme4[, , common_subj, drop = FALSE]
      pv_mm_aligned   <- pv_mm[, , common_subj, drop = FALSE]

      pv_diff <- abs(pv_lme4_aligned - pv_mm_aligned)
      cat("Comparison postVar arrays (", length(common_subj), "subjects):\n")
      cat("max abs diff (postVar): ", max(pv_diff), " [tol 1e-3]\n")
      cat("mean abs diff (postVar):", mean(pv_diff), "\n\n")

      # Per-element breakdown
      cat("Max diff by position [1,1] (Intercept var):",
          max(abs(pv_lme4_aligned[1,1,] - pv_mm_aligned[1,1,])), "\n")
      cat("Max diff by position [2,2] (Days var):     ",
          max(abs(pv_lme4_aligned[2,2,] - pv_mm_aligned[2,2,])), "\n")
      cat("Max diff by position [1,2] / [2,1] (cov): ",
          max(abs(pv_lme4_aligned[1,2,] - pv_mm_aligned[1,2,])), "\n\n")

      # Print comparison for all subjects
      cat("lme4   postVar [1,1] per subject:", round(pv_lme4_aligned[1,1,], 4), "\n")
      cat("mixeff postVar [1,1] per subject:", round(pv_mm_aligned[1,1,], 4), "\n\n")
      cat("lme4   postVar [2,2] per subject:", round(pv_lme4_aligned[2,2,], 4), "\n")
      cat("mixeff postVar [2,2] per subject:", round(pv_mm_aligned[2,2,], 4), "\n\n")

      # Conditional SDs (sqrt of diagonal postVar entries)
      csd_lme4_int  <- sqrt(pv_lme4_aligned[1,1,])
      csd_mm_int    <- sqrt(pv_mm_aligned[1,1,])
      csd_lme4_days <- sqrt(pv_lme4_aligned[2,2,])
      csd_mm_days   <- sqrt(pv_mm_aligned[2,2,])

      cat("--- Conditional SDs ---\n")
      cat("lme4   condSD (Intercept):", round(csd_lme4_int, 4), "\n")
      cat("mixeff condSD (Intercept):", round(csd_mm_int, 4), "\n")
      cat("max abs diff condSD (Intercept):", max(abs(csd_lme4_int - csd_mm_int)),
          " [tol 1e-3]\n\n")
      cat("lme4   condSD (Days):", round(csd_lme4_days, 4), "\n")
      cat("mixeff condSD (Days):", round(csd_mm_days, 4), "\n")
      cat("max abs diff condSD (Days):", max(abs(csd_lme4_days - csd_mm_days)),
          " [tol 1e-3]\n\n")

      within_tol <- max(pv_diff) <= 1e-3
      cat("RESULT: postVar", if (within_tol) "WITHIN tolerance 1e-3" else "BEYOND tolerance 1e-3",
          "\n")
    }
  }
}

## ─── 10. Summary table ────────────────────────────────────────────────────
cat("\n=== SUMMARY TABLE ===\n\n")
cat(sprintf("%-30s %15s %10s %10s\n", "Quantity", "max_abs_diff", "tol", "status"))
cat(strrep("-", 70), "\n")

quantities <- list(
  list(name = "fixef",           diff = max(fe_diff),                  tol = 1e-4),
  list(name = "SE(fixef)",       diff = max(se_diff),                  tol = 1e-4),
  list(name = "sigma",           diff = abs(sig_lme4 - sig_mm),        tol = 1e-4),
  list(name = "VarCorr SD",      diff = max(sd_diff),                  tol = 1e-3),
  list(name = "logLik",          diff = abs(ll_lme4 - ll_mm),          tol = 1e-3),
  list(name = "AIC",             diff = abs(aic_lme4 - aic_mm),        tol = 2e-3),
  list(name = "BIC",             diff = abs(bic_lme4 - bic_mm),        tol = 2e-3),
  list(name = "ranef(Intercept)",diff = max(re_diff_int),              tol = 1e-4),
  list(name = "ranef(Days)",     diff = max(re_diff_days),             tol = 1e-4),
  list(name = "fitted",          diff = max(fv_diff),                  tol = 1e-4)
)

for (q in quantities) {
  status <- if (is.finite(q$diff) && q$diff <= q$tol) "PASS" else "FAIL"
  cat(sprintf("%-30s %15.2e %10.2e %10s\n", q$name, q$diff, q$tol, status))
}

# condVar summary row
if (inherits(mm_condvar_result, "error")) {
  cat(sprintf("%-30s %15s %10.2e %10s\n", "postVar(condVar=TRUE)", "ERROR", 1e-3, "ERROR"))
} else {
  unavail <- !is.null(attr(mm_condvar_result, "mm_unavailable_reason")) ||
             !is.null(attr(mm_condvar_result$Subject, "mm_unavailable_reason"))
  if (unavail) {
    cat(sprintf("%-30s %15s %10.2e %10s\n",
                "postVar(condVar=TRUE)", "typed-refusal", 1e-3, "REFUSAL"))
  } else if (!is.null(common_subj)) {
    pv_mm2    <- attr(mm_condvar_result$Subject, "postVar")
    pv_lme4_2 <- pv_lme4[, , common_subj, drop = FALSE]
    pv_mm_2   <- pv_mm2[, , common_subj, drop = FALSE]
    pv_max    <- max(abs(pv_lme4_2 - pv_mm_2))
    status    <- if (pv_max <= 1e-3) "PASS" else "FAIL"
    cat(sprintf("%-30s %15.2e %10.2e %10s\n",
                "postVar(condVar=TRUE)", pv_max, 1e-3, status))
  }
}

cat("\nDone.\n")
