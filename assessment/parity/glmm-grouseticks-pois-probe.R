## Parity probe: glmm-grouseticks-pois
## Dataset: grouseticks (lme4)
## Formula: TICKS ~ YEAR + HEIGHT + (1|BROOD)
## Family: poisson(link="log")
## Compares lme4::glmer vs mixeff::glmm
## Tolerances: fixef 1e-4, theta 1e-3, logLik 1e-3, sigma 1e-4

suppressMessages({
  library(lme4)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
tryCatch({
  cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
}, error = function(e) cat("lmerTest: not loaded\n"))
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

## ---- Dataset ----
data(grouseticks)
cat("=== DATASET ===\n")
cat("nrow:", nrow(grouseticks), " ncol:", ncol(grouseticks), "\n")
print(head(grouseticks[, c("INDEX","TICKS","BROOD","HEIGHT","YEAR")], 6))
cat("YEAR levels:", paste(levels(grouseticks$YEAR), collapse=", "), "\n")
cat("BROOD n levels:", nlevels(grouseticks$BROOD), "\n\n")

FORMULA <- TICKS ~ YEAR + HEIGHT + (1 | BROOD)

## ---- lme4 fit ----
cat("=== lme4 FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- glmer(FORMULA, data = grouseticks, family = poisson(link = "log"))
})
cat("lme4 wall-clock (seconds):", round(t_lme4["elapsed"], 4), "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov --\n"); print(as.matrix(vcov(fit_lme4)))
cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))
cat("-- sigma --\n"); print(sigma(fit_lme4))
cat("-- logLik --\n"); print(logLik(fit_lme4))
cat("-- AIC --\n"); print(AIC(fit_lme4))
cat("-- BIC --\n"); print(BIC(fit_lme4))
cat("-- ranef (head) --\n")
re_lme4 <- lme4::ranef(fit_lme4)
print(head(re_lme4$BROOD, 8))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4), 6))
cat("-- convergence --\n")
conv_msgs <- fit_lme4@optinfo$conv$lme4
cat("Convergence warnings:", if (length(conv_msgs)) paste(conv_msgs, collapse="; ") else "none", "\n\n")

## ---- mixeff fit ----
cat("=== mixeff FIT ===\n")
t_mm <- system.time({
  fit_mm <- tryCatch(
    glmm(FORMULA, data = grouseticks, family = poisson(link = "log")),
    error = function(e) e
  )
})
cat("mixeff wall-clock (seconds):", round(t_mm["elapsed"], 4), "\n")

if (inherits(fit_mm, "condition")) {
  cat("** mixeff ERRORED **\n")
  cat("Error class:", paste(class(fit_mm), collapse=", "), "\n")
  cat("Error message:", conditionMessage(fit_mm), "\n\n")
} else {
  cat("fit_status:", fit_mm$fit_status, "\n\n")
  cat("-- fixef --\n"); print(fixef(fit_mm))
  cat("-- SE --\n"); print(fit_mm$std_errors)
  cat("-- vcov --\n"); print(fit_mm$fixed_effect_vcov)
  cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
  cat("-- sigma --\n"); print(sigma(fit_mm))
  cat("-- logLik --\n"); print(logLik(fit_mm))
  cat("-- AIC --\n"); print(AIC(fit_mm))
  cat("-- BIC --\n"); print(BIC(fit_mm))
  cat("-- ranef (head) --\n")
  re_mm <- ranef(fit_mm)
  print(head(re_mm$BROOD, 8))
  cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted, 6))
  cat("-- deviance --\n"); print(deviance(fit_mm))
}

## ---- Numerical comparison ----
cat("\n=== NUMERICAL COMPARISON ===\n")

fmt_cmp <- function(label, lme4_val, mm_val, tol, pad=40) {
  if (is.na(mm_val)) {
    cat(sprintf("%-*s lme4=%-20.10f  mm=NA            [MISSING]\n",
                pad, label, lme4_val))
    return(invisible(NA))
  }
  diff <- abs(lme4_val - mm_val)
  status <- if (diff <= tol) "[WITHIN-TOL]" else "[EXCEEDS-TOL*]"
  cat(sprintf("%-*s lme4=%-20.10f  mm=%-20.10f  maxAbsDiff=%.3e  tol=%.0e  %s\n",
              pad, label, lme4_val, mm_val, diff, tol, status))
  invisible(diff)
}

TOL_FIXEF  <- 1e-4
TOL_THETA  <- 1e-3
TOL_LOGLIK <- 1e-3
TOL_SIGMA  <- 1e-4
TOL_AIC    <- 2e-3
TOL_FITTED <- 1e-3

if (!inherits(fit_mm, "condition")) {
  # fixef — compare by position; note mixeff uses "YEAR: 96" lme4 uses "YEAR96"
  fe_lme4 <- lme4::fixef(fit_lme4)
  fe_mm   <- fixef(fit_mm)
  cat(sprintf("fixef names: lme4=[%s]  mm=[%s]\n",
              paste(names(fe_lme4), collapse=", "),
              paste(names(fe_mm),   collapse=", ")))
  n_fe <- min(length(fe_lme4), length(fe_mm))
  for (i in seq_len(n_fe)) {
    lbl <- paste0("fixef[", i, "] ", names(fe_lme4)[i])
    fmt_cmp(lbl, fe_lme4[[i]], fe_mm[[i]], TOL_FIXEF)
  }
  if (length(fe_lme4) != length(fe_mm)) {
    cat(sprintf("fixef length mismatch: lme4=%d  mm=%d\n", length(fe_lme4), length(fe_mm)))
  }

  # SE
  se_lme4 <- sqrt(diag(vcov(fit_lme4)))
  se_mm   <- fit_mm$std_errors
  n_se <- min(length(se_lme4), length(se_mm))
  for (i in seq_len(n_se)) {
    lbl <- paste0("SE[", i, "] ", names(se_lme4)[i])
    fmt_cmp(lbl, se_lme4[[i]], se_mm[[i]], TOL_FIXEF)
  }
  if (length(se_lme4) != length(se_mm)) {
    cat(sprintf("SE length mismatch: lme4=%d  mm=%d\n", length(se_lme4), length(se_mm)))
  }

  # vcov diagonal — by position
  vc_lme4 <- as.matrix(vcov(fit_lme4))
  vc_mm   <- fit_mm$fixed_effect_vcov
  n_vc <- min(nrow(vc_lme4), if (!is.null(vc_mm)) nrow(vc_mm) else 0L)
  for (i in seq_len(n_vc)) {
    nm <- rownames(vc_lme4)[i]
    mm_v <- vc_mm[i, i]
    fmt_cmp(paste0("vcov[", i, ",", i, "] ", nm), vc_lme4[i,i], mm_v, TOL_FIXEF^2 * 10)
  }

  # theta (Cholesky scale parameter)
  theta_lme4 <- lme4::getME(fit_lme4, "theta")
  theta_mm   <- fit_mm$theta
  if (length(theta_lme4) == length(theta_mm)) {
    for (i in seq_along(theta_lme4)) {
      fmt_cmp(paste0("theta[", i, "]"), theta_lme4[i], theta_mm[i], TOL_THETA)
    }
  } else {
    cat(sprintf("theta length mismatch: lme4=%d  mm=%d\n",
                length(theta_lme4), length(theta_mm)))
  }

  # VarCorr (BROOD variance)
  vc_lme4_obj <- lme4::VarCorr(fit_lme4)
  brood_var_lme4 <- as.numeric(vc_lme4_obj$BROOD[1,1])
  vc_mm_obj  <- VarCorr(fit_mm)
  brood_var_mm <- tryCatch({
    tbl <- vc_mm_obj$table
    tbl$variance[tbl$group == "BROOD"][1]
  }, error = function(e) NA_real_)
  fmt_cmp("VarCorr BROOD variance", brood_var_lme4, brood_var_mm, TOL_THETA)

  # sigma (dispersion; for Poisson fixed at 1)
  fmt_cmp("sigma (dispersion)", sigma(fit_lme4), sigma(fit_mm), TOL_SIGMA)

  # logLik
  fmt_cmp("logLik", as.numeric(logLik(fit_lme4)), as.numeric(logLik(fit_mm)), TOL_LOGLIK)

  # AIC / BIC
  fmt_cmp("AIC", AIC(fit_lme4), AIC(fit_mm), TOL_AIC)
  fmt_cmp("BIC", BIC(fit_lme4), BIC(fit_mm), TOL_AIC)

  # deviance
  fmt_cmp("deviance", deviance(fit_lme4), deviance(fit_mm), TOL_LOGLIK * 2)

  # ranef max abs diff (BROOD)
  re_brood_lme4 <- re_lme4$BROOD[["(Intercept)"]]
  re_brood_mm   <- re_mm$BROOD[["(Intercept)"]]
  # align by rownames
  common_broods <- intersect(rownames(re_lme4$BROOD), rownames(re_mm$BROOD))
  if (length(common_broods) > 0) {
    diff_re <- max(abs(re_lme4$BROOD[common_broods, "(Intercept)"] -
                       re_mm$BROOD[common_broods, "(Intercept)"]))
    cat(sprintf("%-40s maxAbsDiff=%.3e  tol=%.0e  %s\n",
                "ranef BROOD max abs diff",
                diff_re, TOL_THETA,
                if (diff_re <= TOL_THETA) "[WITHIN-TOL]" else "[EXCEEDS-TOL*]"))
  } else {
    cat("ranef BROOD: no common levels to compare\n")
  }

  # fitted max abs diff
  ft_lme4 <- fitted(fit_lme4)
  ft_mm   <- fit_mm$fitted
  if (length(ft_lme4) == length(ft_mm)) {
    diff_ft <- max(abs(ft_lme4 - ft_mm))
    cat(sprintf("%-40s maxAbsDiff=%.3e  tol=%.0e  %s\n",
                "fitted max abs diff",
                diff_ft, TOL_FITTED,
                if (diff_ft <= TOL_FITTED) "[WITHIN-TOL]" else "[EXCEEDS-TOL*]"))
  } else {
    cat(sprintf("fitted length mismatch: lme4=%d  mm=%d\n",
                length(ft_lme4), length(ft_mm)))
  }

  # speed
  cat(sprintf("\nwall-clock elapsed  lme4=%.4fs  mm=%.4fs  ratio(mm/lme4)=%.2fx",
              t_lme4["elapsed"], t_mm["elapsed"],
              t_mm["elapsed"] / t_lme4["elapsed"]))
  if (t_mm["elapsed"] < t_lme4["elapsed"]) {
    cat(sprintf("  (mixeff %.1fx faster)\n", t_lme4["elapsed"] / t_mm["elapsed"]))
  } else {
    cat(sprintf("  (lme4 %.1fx faster)\n", t_mm["elapsed"] / t_lme4["elapsed"]))
  }
} else {
  cat("** Cannot compare: mixeff errored **\n")
  cat("lme4 fixef:", paste(names(lme4::fixef(fit_lme4)),
                            round(lme4::fixef(fit_lme4), 6), sep="=",
                            collapse="  "), "\n")
  cat("lme4 logLik:", as.numeric(logLik(fit_lme4)), "\n")
}

cat("\n=== DONE ===\n")
