## Parity probe: lmm-sleep-uncorr
## Dataset: sleepstudy
## Formula: Reaction ~ Days + (Days||Subject)
## Focus: zero-correlation || expansion

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n")
cat("\n")

data("sleepstudy", package = "lme4")

## ── 1. lme4/lmerTest fit ─────────────────────────────────────────────────────
cat("=== lme4/lmerTest FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmerTest::lmer(Reaction ~ Days + (Days || Subject),
                              data = sleepstudy, REML = TRUE)
})
cat("lme4 fit time:", t_lme4["elapsed"], "s\n\n")

# Use lme4:: explicitly to avoid masking by mixeff
fe_lme4  <- lme4::fixef(fit_lme4)
cat("--- fixef (lme4) ---\n")
print(fe_lme4)

cat("\n--- SE (lme4) ---\n")
se_lme4 <- sqrt(diag(as.matrix(vcov(fit_lme4))))
print(se_lme4)

cat("\n--- vcov (lme4) ---\n")
vcov_lme4 <- as.matrix(vcov(fit_lme4))
print(vcov_lme4)

cat("\n--- VarCorr (lme4) ---\n")
vc_lme4 <- lme4::VarCorr(fit_lme4)
print(vc_lme4)
vc_lme4_df <- as.data.frame(vc_lme4)
cat("\nVarCorr as data.frame:\n")
print(vc_lme4_df)

cat("\n--- theta (lme4) ---\n")
theta_lme4 <- lme4::getME(fit_lme4, "theta")
print(theta_lme4)

cat("\n--- sigma (lme4) ---\n")
sigma_lme4 <- sigma(fit_lme4)
cat(sigma_lme4, "\n")

cat("\n--- logLik (lme4) ---\n")
ll_lme4 <- as.numeric(logLik(fit_lme4))
cat(ll_lme4, "\n")

cat("\n--- AIC (lme4) ---\n")
aic_lme4 <- AIC(fit_lme4)
cat(aic_lme4, "\n")

cat("\n--- BIC (lme4) ---\n")
bic_lme4 <- BIC(fit_lme4)
cat(bic_lme4, "\n")

cat("\n--- ranef (lme4, first 5 rows) ---\n")
re_lme4 <- lme4::ranef(fit_lme4)
print(head(re_lme4$Subject, 5))

cat("\n--- fitted (lme4, first 10) ---\n")
fitted_lme4 <- fitted(fit_lme4)
cat(head(fitted_lme4, 10), "\n")

cat("\n--- convergence (lme4) ---\n")
conv_msgs <- fit_lme4@optinfo$conv$lme4$messages
if (is.null(conv_msgs) || length(conv_msgs) == 0) {
  cat("CONVERGED (no warnings)\n")
} else {
  cat("WARNINGS:", paste(conv_msgs, collapse = "; "), "\n")
}

## ── 2. mixeff fit ─────────────────────────────────────────────────────────────
cat("\n\n=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    mixeff::lmm(Reaction ~ Days + (Days || Subject),
                data = sleepstudy,
                REML = TRUE,
                control = mixeff::mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff fit time:", t_mixeff["elapsed"], "s\n\n")

if (inherits(fit_mm, "error")) {
  cat("ERROR from mixeff:\n")
  cat(conditionMessage(fit_mm), "\n")
  cat("Class:", paste(class(fit_mm), collapse = ", "), "\n")
  quit(save = "no")
}

cat("--- fit_status (mixeff) ---\n")
cat(fit_mm$fit_status, "\n")

cat("\n--- fixef (mixeff) ---\n")
fe_mm <- mixeff::fixef(fit_mm)
print(fe_mm)

cat("\n--- SE (mixeff) ---\n")
se_mm <- fit_mm$std_errors
print(se_mm)

cat("\n--- vcov (mixeff) ---\n")
print(fit_mm$fixed_effect_vcov)

cat("\n--- VarCorr (mixeff) ---\n")
vc_mm <- mixeff::VarCorr(fit_mm)
print(vc_mm$table)
cat("residual_sd:", vc_mm$residual_sd, "\n")

cat("\n--- theta (mixeff) ---\n")
theta_mm <- fit_mm$theta
print(theta_mm)

cat("\n--- sigma (mixeff) ---\n")
sigma_mm <- fit_mm$sigma
cat(sigma_mm, "\n")

cat("\n--- logLik (mixeff) ---\n")
ll_mm <- as.numeric(logLik(fit_mm))
cat(ll_mm, "\n")

cat("\n--- AIC (mixeff) ---\n")
aic_mm <- AIC(fit_mm)
cat(aic_mm, "\n")

cat("\n--- BIC (mixeff) ---\n")
bic_mm <- BIC(fit_mm)
cat(bic_mm, "\n")

cat("\n--- ranef (mixeff, first 5 rows) ---\n")
re_mm <- mixeff::ranef(fit_mm)
print(head(re_mm$Subject, 5))

cat("\n--- fitted (mixeff, first 10) ---\n")
fitted_mm <- fit_mm$fitted
cat(head(fitted_mm, 10), "\n")

## ── 3. Comparisons ────────────────────────────────────────────────────────────
cat("\n\n=== COMPARISONS ===\n")

TOL_FIXEF  <- 1e-4
TOL_THETA  <- 1e-3
TOL_LOGLIK <- 1e-3
TOL_SIGMA  <- 1e-4

# fixef
common_fe <- intersect(names(fe_lme4), names(fe_mm))
diff_fixef <- abs(fe_lme4[common_fe] - fe_mm[common_fe])
cat("\n--- fixef differences ---\n")
print(diff_fixef)
cat("Max |diff| fixef:", max(diff_fixef),
    " tol:", TOL_FIXEF,
    " PASS:", max(diff_fixef) <= TOL_FIXEF, "\n")

# SE
common_se <- intersect(names(se_lme4), names(se_mm))
if (length(common_se) > 0) {
  diff_se <- abs(se_lme4[common_se] - se_mm[common_se])
  cat("\n--- SE differences ---\n")
  print(diff_se)
  cat("Max |diff| SE:", max(diff_se),
      " tol:", TOL_FIXEF,
      " PASS:", max(diff_se) <= TOL_FIXEF, "\n")
} else {
  cat("\n--- SE: no common names to compare (lme4:", names(se_lme4),
      " vs mixeff:", names(se_mm), ") ---\n")
  # Compare by position if same length
  if (length(se_lme4) == length(se_mm)) {
    diff_se_pos <- abs(se_lme4 - se_mm)
    cat("By position max |diff|:", max(diff_se_pos), "\n")
  }
}

# sigma
diff_sigma <- abs(sigma_lme4 - sigma_mm)
cat("\n--- sigma ---\n")
cat("lme4:", sigma_lme4, "  mixeff:", as.numeric(sigma_mm),
    "  |diff|:", diff_sigma,
    "  PASS:", diff_sigma <= TOL_SIGMA, "\n")

# logLik
diff_ll <- abs(ll_lme4 - ll_mm)
cat("\n--- logLik ---\n")
cat("lme4:", ll_lme4, "  mixeff:", ll_mm,
    "  |diff|:", diff_ll,
    "  PASS:", diff_ll <= TOL_LOGLIK, "\n")

# AIC / BIC
diff_aic <- abs(aic_lme4 - aic_mm)
diff_bic <- abs(bic_lme4 - bic_mm)
cat("\n--- AIC ---\n")
cat("lme4:", aic_lme4, "  mixeff:", aic_mm,
    "  |diff|:", diff_aic, "\n")
cat("--- BIC ---\n")
cat("lme4:", bic_lme4, "  mixeff:", bic_mm,
    "  |diff|:", diff_bic, "\n")

# theta
cat("\n--- theta ---\n")
cat("lme4 theta:", theta_lme4, "\n")
cat("mixeff theta:", theta_mm, "\n")
theta_lme4_s <- sort(theta_lme4)
theta_mm_s   <- sort(theta_mm)
if (length(theta_lme4_s) == length(theta_mm_s)) {
  diff_theta <- abs(theta_lme4_s - theta_mm_s)
  cat("Sorted |diff| theta:", diff_theta, "\n")
  cat("Max |diff| theta:", max(diff_theta),
      " tol:", TOL_THETA,
      " PASS:", max(diff_theta) <= TOL_THETA, "\n")
} else {
  cat("theta LENGTH MISMATCH: lme4 =", length(theta_lme4),
      " mixeff =", length(theta_mm), "\n")
}

# VarCorr: compare variance components
cat("\n--- VarCorr variance components ---\n")
# lme4 gives us sdcor form; extract variances
vc_lme4_df2 <- as.data.frame(vc_lme4)
cat("lme4 VarCorr df:\n")
print(vc_lme4_df2)
cat("\nmixeff VarCorr table:\n")
print(vc_mm$table)

# ranef comparison
cat("\n--- ranef (all levels, column-by-column) ---\n")
re_lme4_sub <- re_lme4$Subject
re_mm_sub   <- re_mm$Subject
cat("lme4 columns:", names(re_lme4_sub), "\n")
cat("mixeff columns:", names(re_mm_sub), "\n")
common_levels <- intersect(rownames(re_lme4_sub), rownames(re_mm_sub))
re_lme4_a <- re_lme4_sub[common_levels, , drop = FALSE]
re_mm_a   <- re_mm_sub[common_levels, , drop = FALSE]
cat("n levels:", length(common_levels), "\n")
if (ncol(re_lme4_a) == ncol(re_mm_a)) {
  for (ci in seq_len(ncol(re_lme4_a))) {
    d <- abs(re_lme4_a[, ci] - re_mm_a[, ci])
    cat(sprintf("  col %d (%s vs %s): max |diff| = %.6f\n",
                ci, names(re_lme4_a)[ci], names(re_mm_a)[ci], max(d)))
  }
} else {
  cat("Column count differs! lme4:", ncol(re_lme4_a),
      " mixeff:", ncol(re_mm_a), "\n")
  # show full data
  cat("lme4 ranef:\n"); print(re_lme4_a)
  cat("mixeff ranef:\n"); print(re_mm_a)
}

# fitted comparison
cat("\n--- fitted values ---\n")
diff_fitted <- abs(fitted_lme4 - fitted_mm)
cat("Max |diff| fitted:", max(diff_fitted), "\n")
cat("Mean |diff| fitted:", mean(diff_fitted), "\n")

## ── 4. Speed ratio ────────────────────────────────────────────────────────────
cat("\n\n=== SPEED (5 reps each) ===\n")
N_REPS <- 5L
t_lme4_rep <- system.time(
  for (i in seq_len(N_REPS)) suppressMessages(
    lmerTest::lmer(Reaction ~ Days + (Days || Subject),
                   data = sleepstudy, REML = TRUE)
  )
)
t_mm_rep <- system.time(
  for (i in seq_len(N_REPS))
    mixeff::lmm(Reaction ~ Days + (Days || Subject),
                data = sleepstudy, REML = TRUE,
                control = mixeff::mm_control(verbose = -1L))
)

t_lme4_per <- t_lme4_rep["elapsed"] / N_REPS
t_mm_per   <- t_mm_rep["elapsed"]   / N_REPS
cat("lme4  mean per fit:", t_lme4_per, "s\n")
cat("mixeff mean per fit:", t_mm_per, "s\n")
cat("Speed ratio (mixeff/lme4):", t_mm_per / t_lme4_per, "\n")

cat("\n=== DONE ===\n")
