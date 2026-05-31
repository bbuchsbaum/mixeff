## inf-bootstrap probe — parametric bootstrap CI comparison
## Dataset: sleepstudy, Formula: Reaction ~ Days + (1|Subject)
## Compares lme4::bootMer vs mixeff confint(method="bootstrap")
##
## Run: Rscript assessment/parity/inf-bootstrap-probe.R

library(lme4)
library(mixeff)

set.seed(42)

cat("=================================================================\n")
cat("inf-bootstrap probe: sleepstudy ~ Reaction ~ Days + (1|Subject)\n")
cat("=================================================================\n\n")

data(sleepstudy, package = "lme4")

## ---------------------------------------------------------------
## 1. Fit both models
## ---------------------------------------------------------------
cat("--- Fitting lme4 model ---\n")
t_lme4_fit <- system.time({
  fit_lme4 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE)
})
cat(sprintf("lme4 fit time: %.3f sec\n\n", t_lme4_fit["elapsed"]))

cat("--- Fitting mixeff model ---\n")
t_mixeff_fit <- system.time({
  fit_mx <- tryCatch(
    lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy),
    error = function(e) e
  )
})
cat(sprintf("mixeff fit time: %.3f sec\n\n", t_mixeff_fit["elapsed"]))

if (inherits(fit_mx, "error")) {
  cat("MIXEFF FIT ERROR:\n")
  cat(conditionMessage(fit_mx), "\n")
  stop("Cannot proceed: mixeff fit failed.")
}

## ---------------------------------------------------------------
## 2. Basic fit quantities
## ---------------------------------------------------------------
cat("=================================================================\n")
cat("SECTION 1: Basic fit quantities\n")
cat("=================================================================\n\n")

# Fixed effects
cat("--- Fixed effects ---\n")
fe_lme4 <- lme4::fixef(fit_lme4)
fe_mx   <- fit_mx$beta
cat("lme4 fixef:   ", paste(round(fe_lme4, 6), collapse=", "), "\n")
cat("mixeff fixef: ", paste(round(fe_mx, 6), collapse=", "), "\n")
fe_diff <- abs(fe_lme4 - fe_mx)
cat("abs diff:     ", paste(round(fe_diff, 8), collapse=", "), "\n")
cat("tol 1e-4 pass:", all(fe_diff < 1e-4), "\n\n")

# SE / vcov
cat("--- Fixed-effect SEs ---\n")
se_lme4 <- sqrt(diag(vcov(fit_lme4)))
se_mx   <- fit_mx$std_errors
cat("lme4 SE:   ", paste(round(se_lme4, 6), collapse=", "), "\n")
cat("mixeff SE: ", paste(round(se_mx, 6), collapse=", "), "\n")
se_diff <- abs(se_lme4 - se_mx)
cat("abs diff:  ", paste(round(se_diff, 8), collapse=", "), "\n\n")

# Theta / VarCorr
cat("--- Random-effect theta ---\n")
theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mx   <- fit_mx$theta
cat("lme4 theta:   ", round(theta_lme4, 6), "\n")
cat("mixeff theta: ", round(theta_mx, 6), "\n")
theta_diff <- abs(theta_lme4 - theta_mx)
cat("abs diff:     ", round(theta_diff, 8), "\n")
cat("tol 1e-3 pass:", theta_diff < 1e-3, "\n\n")

# Sigma
cat("--- Residual sigma ---\n")
sigma_lme4 <- sigma(fit_lme4)
sigma_mx   <- fit_mx$sigma
cat("lme4 sigma:   ", round(sigma_lme4, 6), "\n")
cat("mixeff sigma: ", round(sigma_mx, 6), "\n")
sigma_diff <- abs(sigma_lme4 - sigma_mx)
cat("abs diff:     ", round(sigma_diff, 8), "\n")
cat("tol 1e-4 pass:", sigma_diff < 1e-4, "\n\n")

# logLik
cat("--- logLik ---\n")
ll_lme4 <- as.numeric(logLik(fit_lme4))
ll_mx   <- tryCatch(as.numeric(fit_mx$logLik), error=function(e) tryCatch(as.numeric(logLik(fit_mx)), error=function(e2) NA_real_))
cat("lme4 logLik:   ", round(ll_lme4, 6), "\n")
cat("mixeff logLik: ", round(ll_mx, 6), "\n")
ll_diff <- abs(ll_lme4 - ll_mx)
cat("abs diff:      ", round(ll_diff, 8), "\n")
cat("tol 1e-3 pass: ", ll_diff < 1e-3, "\n\n")

# AIC / BIC
cat("--- AIC / BIC ---\n")
aic_lme4 <- AIC(fit_lme4); bic_lme4 <- BIC(fit_lme4)
aic_mx   <- tryCatch(AIC(fit_mx), error=function(e) paste("ERROR:", conditionMessage(e)))
bic_mx   <- tryCatch(BIC(fit_mx), error=function(e) paste("ERROR:", conditionMessage(e)))
cat("lme4 AIC:", round(aic_lme4, 4), "  BIC:", round(bic_lme4, 4), "\n")
cat("mixeff AIC:", if(is.numeric(aic_mx)) round(aic_mx,4) else aic_mx,
    "  BIC:", if(is.numeric(bic_mx)) round(bic_mx,4) else bic_mx, "\n\n")

# ranef
cat("--- Random effects (ranef) ---\n")
ranef_lme4 <- lme4::ranef(fit_lme4)$Subject[,1]
ranef_mx   <- tryCatch(ranef(fit_mx)$Subject[,1],
                       error=function(e) paste("ERROR:", conditionMessage(e)))
if (is.numeric(ranef_mx)) {
  ranef_diff_max <- max(abs(sort(ranef_lme4) - sort(ranef_mx)))
  cat("max abs diff (sorted): ", round(ranef_diff_max, 8), "\n\n")
} else {
  cat("mixeff ranef:", ranef_mx, "\n\n")
}

# fitted
cat("--- Fitted values (first 6) ---\n")
fitted_lme4 <- fitted(fit_lme4)
fitted_mx   <- tryCatch(fitted(fit_mx), error=function(e) paste("ERROR:", conditionMessage(e)))
if (is.numeric(fitted_mx)) {
  fitted_diff_max <- max(abs(fitted_lme4 - fitted_mx))
  cat("max abs diff (fitted): ", round(fitted_diff_max, 8), "\n\n")
} else {
  cat("mixeff fitted:", fitted_mx, "\n\n")
}

## ---------------------------------------------------------------
## 3. Parametric bootstrap CI — lme4::bootMer
## ---------------------------------------------------------------
cat("=================================================================\n")
cat("SECTION 2: Parametric bootstrap CI — lme4::bootMer\n")
cat("=================================================================\n\n")

NSIM <- 499L

cat(sprintf("Running lme4::bootMer with nsim=%d (parametric, REML=FALSE)...\n", NSIM))
t_boot_lme4 <- system.time({
  boot_lme4 <- tryCatch(
    lme4::bootMer(fit_lme4, FUN = lme4::fixef, nsim = NSIM, type = "parametric",
            use.u = FALSE, re.form = NA, seed = 42),
    error = function(e) e
  )
})
cat(sprintf("lme4 bootMer time: %.3f sec\n\n", t_boot_lme4["elapsed"]))

if (inherits(boot_lme4, "error")) {
  cat("lme4 bootMer ERROR:\n")
  cat(conditionMessage(boot_lme4), "\n")
  lme4_boot_ci <- NULL
} else {
  lme4_boot_ci <- tryCatch(
    confint(boot_lme4, type = "perc"),
    error = function(e) { cat("lme4 confint error:", conditionMessage(e), "\n"); NULL }
  )
  cat("lme4 bootMer CI (percentile):\n")
  print(lme4_boot_ci)
  cat("\n")
}

## ---------------------------------------------------------------
## 4. Parametric bootstrap CI — mixeff confint(method="bootstrap")
## ---------------------------------------------------------------
cat("=================================================================\n")
cat("SECTION 3: Parametric bootstrap CI — mixeff\n")
cat("=================================================================\n\n")

cat(sprintf("Running mixeff confint(method='bootstrap') with nsim=%d...\n", NSIM))
t_boot_mx <- system.time({
  boot_mx <- tryCatch(
    confint(fit_mx, method = "bootstrap",
            bootstrap = bootstrap_control(nsim = NSIM, seed = 42L)),
    error = function(e) e,
    warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      withCallingHandlers(
        confint(fit_mx, method = "bootstrap",
                bootstrap = bootstrap_control(nsim = NSIM, seed = 42L)),
        warning = function(w2) invokeRestart("muffleWarning")
      )
    }
  )
})
cat(sprintf("mixeff bootstrap CI time: %.3f sec\n\n", t_boot_mx["elapsed"]))

cat("mixeff bootstrap CI output:\n")
if (inherits(boot_mx, "error")) {
  cat("ERROR:", conditionMessage(boot_mx), "\n\n")
  mx_boot_ci <- NULL
} else {
  print(boot_mx)
  mx_boot_ci <- boot_mx
  cat("\n")
}

## ---------------------------------------------------------------
## 5. Also test test_effect with bootstrap method
## ---------------------------------------------------------------
cat("=================================================================\n")
cat("SECTION 4: test_effect() with bootstrap method\n")
cat("=================================================================\n\n")

cat("Running test_effect(fit_mx, 'Days', method='bootstrap')...\n")
t_te_boot <- system.time({
  te_boot <- tryCatch(
    test_effect(fit_mx, "Days", method = "bootstrap",
                bootstrap = bootstrap_control(nsim = NSIM, seed = 42L)),
    error = function(e) e
  )
})
cat(sprintf("test_effect bootstrap time: %.3f sec\n\n", t_te_boot["elapsed"]))

if (inherits(te_boot, "error")) {
  cat("test_effect bootstrap ERROR:", conditionMessage(te_boot), "\n\n")
} else {
  cat("test_effect bootstrap result:\n")
  print(te_boot)
  cat("\n")
}

## Also test with bootstrap_lrt method
cat("Running test_effect(fit_mx, 'Days', method='bootstrap_lrt')...\n")
t_te_blrt <- system.time({
  te_blrt <- tryCatch(
    test_effect(fit_mx, "Days", method = "bootstrap_lrt",
                bootstrap = bootstrap_control(nsim = NSIM, seed = 42L)),
    error = function(e) e
  )
})
cat(sprintf("test_effect bootstrap_lrt time: %.3f sec\n\n", t_te_blrt["elapsed"]))

if (inherits(te_blrt, "error")) {
  cat("test_effect bootstrap_lrt ERROR:", conditionMessage(te_blrt), "\n\n")
} else {
  cat("test_effect bootstrap_lrt result:\n")
  print(te_blrt)
  cat("\n")
}

## ---------------------------------------------------------------
## 6. Distributional comparison of bootstrap distributions
## ---------------------------------------------------------------
cat("=================================================================\n")
cat("SECTION 5: Bootstrap distribution comparison\n")
cat("=================================================================\n\n")

if (!is.null(lme4_boot_ci) && !is.null(mx_boot_ci)) {
  # Extract CI values as numeric matrices for comparison
  lme4_ci_mat <- as.matrix(lme4_boot_ci)
  mx_ci_mat   <- as.matrix(mx_boot_ci)

  # Row-align: lme4 has (Intercept), Days; mx should too
  cat("lme4 CI:\n"); print(lme4_ci_mat)
  cat("mixeff CI:\n"); print(mx_ci_mat)

  # Match rows by name
  common_parms <- intersect(rownames(lme4_ci_mat), rownames(mx_ci_mat))
  if (length(common_parms) > 0) {
    diff_lower <- abs(lme4_ci_mat[common_parms, 1] - mx_ci_mat[common_parms, 1])
    diff_upper <- abs(lme4_ci_mat[common_parms, 2] - mx_ci_mat[common_parms, 2])
    cat("\nCI bound differences (lme4 vs mixeff):\n")
    cat("  lower bound abs diff:", paste(round(diff_lower, 4), collapse=", "), "\n")
    cat("  upper bound abs diff:", paste(round(diff_upper, 4), collapse=", "), "\n")
    max_ci_diff <- max(c(diff_lower, diff_upper))
    cat("  max abs diff across all CI bounds:", round(max_ci_diff, 6), "\n\n")
    cat("Note: bootstrap CI comparison is distributional (not exact) —\n")
    cat("differences arise from different RNG, refit strategy, etc.\n\n")
  } else {
    cat("No common parameter names between lme4 and mixeff CI output.\n")
    cat("lme4 rownames:", rownames(lme4_ci_mat), "\n")
    cat("mixeff rownames:", rownames(mx_ci_mat), "\n\n")
  }
} else {
  cat("Cannot compare: one or both bootstrap CI calls failed.\n\n")
}

## ---------------------------------------------------------------
## 7. Speed comparison
## ---------------------------------------------------------------
cat("=================================================================\n")
cat("SECTION 6: Speed summary\n")
cat("=================================================================\n\n")

cat(sprintf("Fit time:         lme4=%.3fs  mixeff=%.3fs  ratio=%.2fx\n",
            t_lme4_fit["elapsed"], t_mixeff_fit["elapsed"],
            t_mixeff_fit["elapsed"] / t_lme4_fit["elapsed"]))
cat(sprintf("Bootstrap CI:     lme4=%.3fs  mixeff=%.3fs  ratio=%.2fx\n",
            t_boot_lme4["elapsed"], t_boot_mx["elapsed"],
            t_boot_mx["elapsed"] / t_boot_lme4["elapsed"]))

cat("\n=================================================================\n")
cat("PROBE COMPLETE\n")
cat("=================================================================\n")
