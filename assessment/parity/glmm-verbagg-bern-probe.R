## glmm-verbagg-bern parity probe
## Dataset: VerbAgg (lme4)
## Formula: r2 ~ Anger + Gender + (1|id) + (1|item)
## Family:  binomial(logit)  — Bernoulli response
##
## Compares: lme4::glmer  vs  mixeff::glmm
## Quantities: fixef, SE, vcov, theta/VarCorr, sigma, logLik, AIC, BIC,
##             ranef, fitted, convergence/refusal status, wall-clock timing.

library(lme4)
library(mixeff)

# ── 0. Load data ──────────────────────────────────────────────────────────────
data("VerbAgg", package = "lme4")
# r2 is an ordered factor "N" / "Y"; need a 0/1 integer
VerbAgg$r2_bin <- as.integer(VerbAgg$r2 == "Y")
formula_str <- r2_bin ~ Anger + Gender + (1|id) + (1|item)

cat("=== Dataset summary ===\n")
cat("nrow:", nrow(VerbAgg), "\n")
cat("r2_bin distribution:\n")
print(table(VerbAgg$r2_bin))
cat("\n")

# ── 1. Fit lme4 model ─────────────────────────────────────────────────────────
cat("=== Fitting lme4::glmer ===\n")
t_lme4 <- system.time({
  fit_lme4 <- suppressMessages(suppressWarnings(
    glmer(formula_str, data = VerbAgg, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 2e5)))
  ))
})
cat("lme4 wall time (elapsed):", t_lme4["elapsed"], "s\n\n")

# ── 2. Fit mixeff model ───────────────────────────────────────────────────────
cat("=== Fitting mixeff::glmm ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    glmm(formula_str, data = VerbAgg,
         family = binomial(link = "logit"),
         control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff wall time (elapsed):", t_mixeff["elapsed"], "s\n\n")

if (inherits(fit_mm, "condition")) {
  cat("=== mixeff ERRORED ===\n")
  cat("Class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat("Message:", conditionMessage(fit_mm), "\n")
  cat("\nStopping — cannot compare quantities.\n")
  quit(status = 0)
}

cat("mixeff fit_status:", fit_mm$fit_status, "\n\n")

# ── 3. Fixed effects ──────────────────────────────────────────────────────────
cat("=== Fixed effects (fixef) ===\n")
fe_lme4  <- lme4::fixef(fit_lme4)
fe_mm    <- mixeff::fixef(fit_mm)

# Align names (mixeff may use different intercept label)
cat("lme4 names :", paste(names(fe_lme4), collapse = ", "), "\n")
cat("mixeff names:", paste(names(fe_mm),  collapse = ", "), "\n")

common_fe <- intersect(names(fe_lme4), names(fe_mm))
fe_diff   <- abs(fe_lme4[common_fe] - fe_mm[common_fe])

cat("\nFixed effect comparison (common terms):\n")
comp_fe <- data.frame(
  lme4   = fe_lme4[common_fe],
  mixeff = fe_mm[common_fe],
  abs_diff = fe_diff
)
print(comp_fe)
cat("Max |fixef diff|:", max(fe_diff), "\n")
cat("Tolerance (1e-4):", max(fe_diff) <= 1e-4, "\n\n")

# ── 4. Standard errors ────────────────────────────────────────────────────────
cat("=== Standard errors ===\n")
se_lme4 <- sqrt(diag(vcov(fit_lme4)))
se_mm   <- fit_mm$std_errors

cat("lme4 SE names :", paste(names(se_lme4), collapse = ", "), "\n")
cat("mixeff SE names:", paste(names(se_mm),  collapse = ", "), "\n")

common_se  <- intersect(names(se_lme4), names(se_mm))
se_diff    <- abs(se_lme4[common_se] - se_mm[common_se])

cat("\nSE comparison (common terms):\n")
comp_se <- data.frame(
  lme4   = se_lme4[common_se],
  mixeff = se_mm[common_se],
  abs_diff = se_diff
)
print(comp_se)
cat("Max |SE diff|:", max(se_diff), "\n\n")

# ── 5. vcov ───────────────────────────────────────────────────────────────────
cat("=== Variance-covariance matrix (vcov) ===\n")
vc_lme4 <- as.matrix(vcov(fit_lme4))
vc_mm   <- tryCatch(as.matrix(vcov(fit_mm)), error = function(e) {
  cat("vcov() error for mixeff:", conditionMessage(e), "\n"); NULL
})
if (!is.null(vc_mm)) {
  common_v <- intersect(rownames(vc_lme4), rownames(vc_mm))
  vc_diff  <- abs(vc_lme4[common_v, common_v] - vc_mm[common_v, common_v])
  cat("Max |vcov diff|:", max(vc_diff), "\n\n")
}

# ── 6. Random-effect variances / theta ────────────────────────────────────────
cat("=== Random-effect variances (VarCorr / theta) ===\n")
vc_lme4_re <- VarCorr(fit_lme4)
vc_mm_re   <- VarCorr(fit_mm)

cat("\nlme4 VarCorr:\n");  print(vc_lme4_re)
cat("\nmixeff VarCorr:\n"); print(vc_mm_re)

theta_lme4 <- getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
cat("\nlme4 theta: ", theta_lme4, "\n")
cat("mixeff theta:", theta_mm, "\n")

# theta has same length?
if (length(theta_lme4) == length(theta_mm)) {
  theta_diff <- abs(theta_lme4 - theta_mm)
  cat("Max |theta diff|:", max(theta_diff), "\n")
  cat("Tolerance (1e-3):", max(theta_diff) <= 1e-3, "\n\n")
} else {
  cat("theta length mismatch: lme4=", length(theta_lme4),
      "mixeff=", length(theta_mm), "\n\n")
}

# ── 7. sigma ─────────────────────────────────────────────────────────────────
cat("=== sigma (dispersion) ===\n")
sig_lme4 <- sigma(fit_lme4)
sig_mm   <- sigma(fit_mm)
cat("lme4 sigma :", sig_lme4, "\n")
cat("mixeff sigma:", sig_mm,  "\n")
cat("Max |sigma diff|:", abs(sig_lme4 - sig_mm), "\n")
cat("Tolerance (1e-4):", abs(sig_lme4 - sig_mm) <= 1e-4, "\n\n")

# ── 8. logLik ─────────────────────────────────────────────────────────────────
cat("=== logLik ===\n")
ll_lme4 <- as.numeric(logLik(fit_lme4))
ll_mm   <- as.numeric(logLik(fit_mm))
cat("lme4 logLik :", ll_lme4, "\n")
cat("mixeff logLik:", ll_mm,  "\n")
cat("Max |logLik diff|:", abs(ll_lme4 - ll_mm), "\n")
cat("Tolerance (1e-3):", abs(ll_lme4 - ll_mm) <= 1e-3, "\n\n")

# ── 9. AIC / BIC ──────────────────────────────────────────────────────────────
cat("=== AIC / BIC ===\n")
aic_lme4 <- AIC(fit_lme4); bic_lme4 <- BIC(fit_lme4)
aic_mm   <- AIC(fit_mm);   bic_mm   <- BIC(fit_mm)
cat("lme4  AIC:", aic_lme4, " BIC:", bic_lme4, "\n")
cat("mixeff AIC:", aic_mm,  " BIC:", bic_mm,  "\n")
cat("|AIC diff|:", abs(aic_lme4 - aic_mm), "\n")
cat("|BIC diff|:", abs(bic_lme4 - bic_mm), "\n\n")

# ── 10. ranef (first 5 rows each group) ───────────────────────────────────────
cat("=== Random effects (ranef, first 5 rows per group) ===\n")
re_lme4 <- ranef(fit_lme4)
re_mm   <- ranef(fit_mm)

cat("lme4 ranef groups  :", paste(names(re_lme4), collapse = ", "), "\n")
cat("mixeff ranef groups:", paste(names(re_mm),   collapse = ", "), "\n\n")

for (grp in intersect(names(re_lme4), names(re_mm))) {
  cat("Group:", grp, "\n")
  r_l <- re_lme4[[grp]]
  r_m <- re_mm[[grp]]
  # align by rownames
  common_ids <- intersect(rownames(r_l), rownames(r_m))
  col_l <- intersect(colnames(r_l), colnames(r_m))
  if (length(common_ids) > 0 && length(col_l) > 0) {
    diff_re <- abs(r_l[common_ids, col_l, drop = FALSE] -
                     r_m[common_ids, col_l, drop = FALSE])
    cat("  Max |ranef diff| (", grp, "):", max(diff_re), "\n")
    cat("  Head lme4:\n");  print(head(r_l[common_ids, col_l, drop = FALSE], 5))
    cat("  Head mixeff:\n"); print(head(r_m[common_ids, col_l, drop = FALSE], 5))
  } else {
    cat("  Cannot align ranef for group:", grp, "\n")
    cat("  lme4 rownames (first 5):", paste(head(rownames(r_l), 5), collapse=", "), "\n")
    cat("  mixeff rownames (first 5):", paste(head(rownames(r_m), 5), collapse=", "), "\n")
  }
  cat("\n")
}

# ── 11. Fitted values (first 10) ─────────────────────────────────────────────
cat("=== Fitted values (first 10) ===\n")
fv_lme4 <- fitted(fit_lme4)
fv_mm   <- fit_mm$fitted

cat("lme4 fitted (first 10):", head(fv_lme4, 10), "\n")
cat("mixeff fitted (first 10):", head(fv_mm, 10), "\n")
cat("Max |fitted diff| (all):", max(abs(fv_lme4 - fv_mm)), "\n\n")

# ── 12. Speed summary ─────────────────────────────────────────────────────────
cat("=== Speed ===\n")
cat("lme4  elapsed:", t_lme4["elapsed"],  "s\n")
cat("mixeff elapsed:", t_mixeff["elapsed"], "s\n")
cat("Ratio lme4/mixeff:", t_lme4["elapsed"] / t_mixeff["elapsed"], "\n\n")

# ── 13. Summary table ─────────────────────────────────────────────────────────
cat("=== SUMMARY TABLE ===\n")
tol_fixef <- 1e-4
tol_theta <- 1e-3
tol_ll    <- 1e-3
tol_sigma <- 1e-4

max_fixef_diff <- if (length(common_fe) > 0) max(abs(fe_lme4[common_fe] - fe_mm[common_fe])) else NA
max_se_diff    <- if (length(common_se) > 0) max(abs(se_lme4[common_se] - se_mm[common_se])) else NA
max_theta_diff <- if (length(theta_lme4) == length(theta_mm)) max(abs(theta_lme4 - theta_mm)) else NA
max_ll_diff    <- abs(ll_lme4 - ll_mm)
max_sig_diff   <- abs(sig_lme4 - sig_mm)
max_aic_diff   <- abs(aic_lme4 - aic_mm)
max_ranef_diff <- tryCatch({
  diffs <- sapply(intersect(names(re_lme4), names(re_mm)), function(grp) {
    r_l <- re_lme4[[grp]]; r_m <- re_mm[[grp]]
    cids <- intersect(rownames(r_l), rownames(r_m))
    ccol <- intersect(colnames(r_l), colnames(r_m))
    if (length(cids) > 0 && length(ccol) > 0)
      max(abs(r_l[cids, ccol, drop=FALSE] - r_m[cids, ccol, drop=FALSE]))
    else NA_real_
  })
  max(diffs, na.rm = TRUE)
}, error = function(e) NA_real_)
max_fitted_diff <- max(abs(fv_lme4 - fv_mm))

results <- data.frame(
  quantity    = c("fixef", "SE", "theta", "logLik", "sigma", "AIC", "ranef", "fitted"),
  lme4_val    = c(
    paste(round(fe_lme4, 5), collapse=" "),
    paste(round(se_lme4, 5), collapse=" "),
    paste(round(theta_lme4, 5), collapse=" "),
    round(ll_lme4, 4),
    round(sig_lme4, 6),
    round(aic_lme4, 4),
    "(see above)",
    "(see above)"
  ),
  max_abs_diff = c(
    round(max_fixef_diff, 8),
    round(max_se_diff, 8),
    round(max_theta_diff, 8),
    round(max_ll_diff, 8),
    round(max_sig_diff, 8),
    round(max_aic_diff, 8),
    round(max_ranef_diff, 8),
    round(max_fitted_diff, 8)
  ),
  within_tol  = c(
    isTRUE(max_fixef_diff <= tol_fixef),
    NA,
    isTRUE(max_theta_diff <= tol_theta),
    isTRUE(max_ll_diff    <= tol_ll),
    isTRUE(max_sig_diff   <= tol_sigma),
    NA, NA, NA
  ),
  stringsAsFactors = FALSE
)

print(results, row.names = FALSE)
cat("\n=== PROBE COMPLETE ===\n")
