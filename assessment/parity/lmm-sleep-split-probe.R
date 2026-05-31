## ============================================================
## Parity probe: lmm-sleep-split
## Dataset: sleepstudy
## Formula: Reaction ~ Days + (1|Subject) + (0+Days|Subject)
## (split-block / || equivalent — uncorrelated random intercept + slope)
## ============================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

cat("=== Package versions ===\n")
cat("lme4:     ", as.character(packageVersion("lme4")), "\n")
cat("lmerTest: ", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff:   ", as.character(packageVersion("mixeff")), "\n\n")

data("sleepstudy", package = "lme4")

TOLERANCES <- list(
  fixef   = 1e-4,
  theta   = 1e-3,
  logLik  = 1e-3,
  sigma   = 1e-4
)

## ============================================================
## 1. FIT BOTH MODELS
## ============================================================

cat("=== Fitting lme4/lmerTest model ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmerTest::lmer(
    Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
    data = sleepstudy,
    REML = TRUE
  )
})
cat(sprintf("lme4 wall time: %.4f s\n\n", t_lme4["elapsed"]))

cat("=== Fitting mixeff model ===\n")
mixeff_error <- NULL
fit_mm <- NULL
t_mm <- system.time({
  fit_mm <- tryCatch(
    mixeff::lmm(
      Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
      data  = sleepstudy,
      REML  = TRUE,
      control = mixeff::mm_control(verbose = -1L)
    ),
    error = function(e) { mixeff_error <<- e; NULL }
  )
})
cat(sprintf("mixeff wall time: %.4f s\n\n", t_mm["elapsed"]))

if (!is.null(mixeff_error)) {
  cat("!!! mixeff ERROR (verbatim) !!!\n")
  cat(conditionMessage(mixeff_error), "\n\n")
} else {
  cat("mixeff fit succeeded.\n\n")
}

## ============================================================
## 2. EXTRACT QUANTITIES — fully namespace-qualified for lme4
## ============================================================

## ---- lme4 quantities (all namespace-qualified to avoid mixeff masking) -------
lme4_fixef  <- lme4::fixef(fit_lme4)
lme4_vcov   <- as.matrix(stats::vcov(fit_lme4))
lme4_se     <- sqrt(diag(lme4_vcov))
lme4_theta  <- lme4::getME(fit_lme4, "theta")
lme4_sigma  <- stats::sigma(fit_lme4)
lme4_loglik <- as.numeric(stats::logLik(fit_lme4))
lme4_aic    <- stats::AIC(fit_lme4)
lme4_bic    <- stats::BIC(fit_lme4)
lme4_ranef  <- lme4::ranef(fit_lme4)
lme4_fitted <- stats::fitted(fit_lme4)
lme4_vc     <- lme4::VarCorr(fit_lme4)

cat("=== lme4 quantities ===\n")
cat("fixef:\n"); print(lme4_fixef)
cat("SE:\n"); print(lme4_se)
cat("theta:\n"); print(lme4_theta)
cat("sigma:", lme4_sigma, "\n")
cat("logLik:", lme4_loglik, "\n")
cat("AIC:", lme4_aic, "  BIC:", lme4_bic, "\n")
cat("VarCorr:\n"); print(lme4_vc)
cat("ranef head:\n"); print(lapply(lme4_ranef, head, 3))
cat("\n")

## ---- mixeff quantities (if available) -----------------------
if (!is.null(fit_mm)) {
  mm_fixef  <- mixeff::fixef(fit_mm)
  mm_vcov   <- tryCatch(as.matrix(stats::vcov(fit_mm)), error = function(e) {
    cat("vcov() ERROR:", conditionMessage(e), "\n"); NULL })
  mm_se     <- if (!is.null(mm_vcov)) sqrt(diag(mm_vcov)) else {
    tryCatch(fit_mm$se, error = function(e) NULL) }
  mm_theta  <- fit_mm$theta
  mm_sigma  <- stats::sigma(fit_mm)
  mm_loglik <- as.numeric(stats::logLik(fit_mm))
  mm_aic    <- stats::AIC(fit_mm)
  mm_bic    <- stats::BIC(fit_mm)
  mm_ranef  <- tryCatch(mixeff::ranef(fit_mm), error = function(e) {
    cat("ranef() ERROR:", conditionMessage(e), "\n"); NULL })
  mm_fitted <- stats::fitted(fit_mm)
  mm_vc     <- tryCatch(mixeff::VarCorr(fit_mm), error = function(e) {
    cat("VarCorr() ERROR:", conditionMessage(e), "\n"); NULL })

  cat("=== mixeff quantities ===\n")
  cat("fixef:\n"); print(mm_fixef)
  cat("SE:\n"); print(mm_se)
  cat("theta:\n"); print(mm_theta)
  cat("sigma:", mm_sigma, "\n")
  cat("logLik:", mm_loglik, "\n")
  cat("AIC:", mm_aic, "  BIC:", mm_bic, "\n")
  cat("VarCorr:\n"); print(mm_vc)
  cat("ranef head:\n"); print(lapply(mm_ranef, head, 3))
  cat("\n")
}

## ============================================================
## 3. COMPARISONS
## ============================================================

tol_check <- function(label, lme4_val, mm_val, tol) {
  if (is.null(mm_val)) {
    cat(sprintf("  %-32s  MISSING (mixeff returned NULL)\n", label))
    return(invisible(NA))
  }
  # align by name if both are named
  if (!is.null(names(lme4_val)) && !is.null(names(mm_val))) {
    common <- intersect(names(lme4_val), names(mm_val))
    if (length(common) == 0) {
      cat(sprintf("  %-32s  NO COMMON NAMES lme4=[%s] mm=[%s]\n",
                  label,
                  paste(names(lme4_val), collapse=","),
                  paste(names(mm_val),   collapse=",")))
      return(invisible(NA))
    }
    lme4_val <- lme4_val[common]
    mm_val   <- mm_val[common]
  }
  lv <- as.numeric(lme4_val)
  mv <- as.numeric(mm_val)
  if (length(lv) != length(mv)) {
    cat(sprintf("  %-32s  LENGTH MISMATCH lme4=%d mm=%d\n",
                label, length(lv), length(mv)))
    return(invisible(NA))
  }
  diffs    <- abs(lv - mv)
  max_diff <- max(diffs)
  status   <- if (max_diff <= tol) "WITHIN-TOL" else "BEYOND-TOL"
  cat(sprintf("  %-32s  max|diff|=%.3e  tol=%.0e  [%s]\n",
              label, max_diff, tol, status))
  invisible(max_diff)
}

cat("=== Comparison summary ===\n")
if (is.null(fit_mm)) {
  cat("  mixeff FIT FAILED — no comparison possible.\n")
  cat("  Error class:", paste(class(mixeff_error), collapse=", "), "\n")
  cat("  Error message:", conditionMessage(mixeff_error), "\n")
} else {

  cat("--- Fixed effects ---\n")
  tol_check("fixef (Intercept)",  lme4_fixef["(Intercept)"],
            mm_fixef["(Intercept)"],  TOLERANCES$fixef)
  tol_check("fixef Days",          lme4_fixef["Days"],
            mm_fixef["Days"],          TOLERANCES$fixef)

  cat("--- Standard errors ---\n")
  if (!is.null(mm_se)) {
    tol_check("SE (Intercept)",    lme4_se["(Intercept)"],
              mm_se["(Intercept)"],    TOLERANCES$fixef)
    tol_check("SE Days",            lme4_se["Days"],
              mm_se["Days"],            TOLERANCES$fixef)
  } else {
    cat("  SE: not available from mixeff\n")
  }

  cat("--- Variance components (theta) ---\n")
  cat(sprintf("  lme4  theta (%d): %s\n", length(lme4_theta),
              paste(round(lme4_theta, 6), collapse=", ")))
  cat(sprintf("  mixeff theta (%d): %s\n", length(mm_theta),
              paste(round(mm_theta,   6), collapse=", ")))
  if (length(lme4_theta) == length(mm_theta)) {
    tol_check("theta (all)",       lme4_theta, mm_theta, TOLERANCES$theta)
  } else {
    cat(sprintf("  theta LENGTH MISMATCH: lme4=%d, mixeff=%d\n",
                length(lme4_theta), length(mm_theta)))
  }

  cat("--- sigma ---\n")
  tol_check("sigma",               lme4_sigma, mm_sigma, TOLERANCES$sigma)

  cat("--- Log-likelihood ---\n")
  tol_check("logLik",              lme4_loglik, mm_loglik, TOLERANCES$logLik)

  cat("--- Information criteria ---\n")
  tol_check("AIC",                 lme4_aic,   mm_aic,   TOLERANCES$logLik)
  tol_check("BIC",                 lme4_bic,   mm_bic,   TOLERANCES$logLik)

  cat("--- Random effects: intercept BLUPs ---\n")
  re_lme4_int <- lme4_ranef$Subject[["(Intercept)"]]
  re_mm_int   <- if (!is.null(mm_ranef)) {
    grp <- Filter(function(df) "(Intercept)" %in% names(df), mm_ranef)
    if (length(grp) > 0) grp[[1]][["(Intercept)"]] else NULL
  } else NULL
  tol_check("ranef intercept BLUPs", re_lme4_int, re_mm_int, TOLERANCES$fixef)

  cat("--- Random effects: slope BLUPs ---\n")
  re_lme4_slp <- lme4_ranef$Subject[["Days"]]
  re_mm_slp   <- if (!is.null(mm_ranef)) {
    grp <- Filter(function(df) "Days" %in% names(df), mm_ranef)
    if (length(grp) > 0) grp[[1]][["Days"]] else NULL
  } else NULL
  tol_check("ranef slope BLUPs", re_lme4_slp, re_mm_slp, TOLERANCES$fixef)

  cat("--- Fitted values ---\n")
  tol_check("fitted (all obs)",    lme4_fitted, mm_fitted, TOLERANCES$fixef)

  cat("\n--- Speed ---\n")
  cat(sprintf("  lme4   elapsed: %.4f s\n", t_lme4["elapsed"]))
  cat(sprintf("  mixeff elapsed: %.4f s\n", t_mm["elapsed"]))
  if (t_lme4["elapsed"] > 0) {
    ratio <- t_mm["elapsed"] / t_lme4["elapsed"]
    cat(sprintf("  mixeff/lme4 ratio: %.2fx  (%s)\n", ratio,
                if (ratio < 1) "mixeff faster" else "mixeff slower"))
  }
}

## ============================================================
## 4. CONVERGENCE / DIAGNOSTIC FLAGS
## ============================================================
cat("\n=== Convergence / diagnostic flags ===\n")
lme4_conv_msgs <- tryCatch({
  msgs <- fit_lme4@optinfo$conv$lme4$messages
  if (is.null(msgs)) "none" else paste(msgs, collapse="; ")
}, error = function(e) paste("(error reading optinfo:", conditionMessage(e), ")"))
cat("lme4 convergence messages:", lme4_conv_msgs, "\n")

if (!is.null(fit_mm)) {
  conv_val <- tryCatch(fit_mm$convergence, error = function(e) "(field absent)")
  sing_val <- tryCatch(fit_mm$singular,    error = function(e) "(field absent)")
  cat("mixeff convergence:", if (is.null(conv_val)) "(NULL/absent)" else conv_val, "\n")
  cat("mixeff singular:   ", if (is.null(sing_val)) "(NULL/absent)" else sing_val, "\n")
}

## ============================================================
## 5. VarCorr detail
## ============================================================
cat("\n=== VarCorr detail ===\n")
cat("--- lme4 VarCorr ---\n")
print(summary(lme4_vc))

if (!is.null(fit_mm) && !is.null(mm_vc)) {
  cat("--- mixeff VarCorr ---\n")
  print(mm_vc)
}

cat("\n=== Done ===\n")
