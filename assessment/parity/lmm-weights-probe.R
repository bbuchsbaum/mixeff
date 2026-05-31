## Empirical parity probe: lmm-weights
## Dataset: simulated (per-row weights)  Formula: y ~ x + (1|g)  REML = TRUE
## Focus: prior weights honored?
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

## ── 1. Simulate data with weights ────────────────────────────────────────────
set.seed(42)
n_groups <- 20
n_per    <- 10
N        <- n_groups * n_per

g <- factor(rep(paste0("g", seq_len(n_groups)), each = n_per))
x <- rnorm(N)

# True parameters
b0    <- 2.0
b1    <- 0.5
sigma_b <- 1.2   # random intercept SD
sigma_e <- 0.8   # residual SD

u <- rnorm(n_groups, 0, sigma_b)
e <- rnorm(N, 0, sigma_e)
y <- b0 + b1 * x + u[as.integer(g)] + e

# Per-row weights (simulating observation-level precision weights,
# e.g. inverse of known measurement variance)
# Use weights that vary meaningfully but not pathologically
wts <- runif(N, min = 0.5, max = 3.0)

dat <- data.frame(y = y, x = x, g = g, w = wts)

cat("=== DATASET ===\n")
cat("nrow:", nrow(dat), "  ncol:", ncol(dat), "\n")
cat("weights summary:\n"); print(summary(dat$w))
cat("\n")

## ── 2. Fit lme4 ──────────────────────────────────────────────────────────────
cat("=== lme4 FIT (with weights) ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(y ~ x + (1 | g),
                   data    = dat,
                   weights = dat$w,
                   REML    = TRUE)
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
cat("-- ranef (first 6) --\n"); print(head(lme4::ranef(fit_lme4)$g))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
cat("-- convergence --\n")
cat("No convergence warning =",
    length(fit_lme4@optinfo$conv$lme4$messages) == 0, "\n\n")

## Also fit WITHOUT weights for comparison reference
cat("=== lme4 FIT (no weights, reference) ===\n")
fit_lme4_nw <- lmer(y ~ x + (1 | g), data = dat, REML = TRUE)
cat("fixef no-weights:", paste(sprintf("%.6f", lme4::fixef(fit_lme4_nw)), collapse = ", "), "\n")
cat("fixef with-weights:", paste(sprintf("%.6f", lme4::fixef(fit_lme4)), collapse = ", "), "\n")
cat("(These should differ if weights are honored)\n\n")

## ── 3. Fit mixeff ─────────────────────────────────────────────────────────────
cat("=== mixeff FIT (with weights) ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    lmm(y ~ x + (1 | g),
        data    = dat,
        weights = dat$w,
        REML    = TRUE,
        control = mm_control(verbose = -1L)),
    error   = function(e) e,
    warning = function(w) w
  )
})
cat("mixeff wall-clock (seconds):", t_mixeff["elapsed"], "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR/WARNING !!!\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm), "\n")
  ## Still try to proceed in case it was a warning
  if (inherits(fit_mm, "error")) {
    cat("\n=== SUMMARY: mixeff failed to fit with weights ===\n")
    quit(status = 0)
  }
}

## If fit_mm is a warning, re-fit capturing result
if (inherits(fit_mm, "warning")) {
  cat("(Re-fitting to capture the actual fit object despite warning)\n")
  fit_mm <- withCallingHandlers(
    lmm(y ~ x + (1 | g), data = dat, weights = dat$w, REML = TRUE,
        control = mm_control(verbose = -1L)),
    warning = function(w) invokeRestart("muffleWarning")
  )
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
cat("-- ranef (first 6) --\n"); print(head(ranef(fit_mm)$g))
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("\n")

## ── 4. Weights-honored sanity check ──────────────────────────────────────────
cat("=== WEIGHTS-HONORED SANITY CHECK ===\n")
cat("Do lme4 and mixeff fixef differ from no-weights fit?\n")
fe_lme4    <- lme4::fixef(fit_lme4)
fe_lme4_nw <- lme4::fixef(fit_lme4_nw)
fe_mm      <- fixef(fit_mm)

diff_lme4_vs_nw <- max(abs(fe_lme4 - fe_lme4_nw))
diff_mm_vs_nw   <- max(abs(fe_mm   - fe_lme4_nw))
diff_lme4_vs_mm <- max(abs(fe_lme4 - fe_mm))

cat(sprintf("lme4(weighted) vs lme4(unweighted) max|diff|: %.6f\n", diff_lme4_vs_nw))
cat(sprintf("mm(weighted)   vs lme4(unweighted) max|diff|: %.6f\n", diff_mm_vs_nw))
cat(sprintf("lme4(weighted) vs mm(weighted)     max|diff|: %.6f  <- KEY\n\n", diff_lme4_vs_mm))
if (diff_lme4_vs_nw < 1e-6) {
  cat("WARNING: lme4 weighted and unweighted fixef are identical -- weights may be ignored or data is insensitive\n")
}
if (diff_lme4_vs_mm < 1e-6 && diff_lme4_vs_nw > 1e-6) {
  cat("GOOD: mixeff and lme4 weighted fixef agree and both differ from unweighted\n")
} else if (diff_lme4_vs_mm > 1e-4) {
  cat("DIVERGENCE: mixeff and lme4 weighted fixef disagree beyond 1e-4\n")
}
cat("\n")

## ── 5. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  diff     <- max(abs(lme4_val - mm_val))
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-36s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse = ", "),
              paste(sprintf(fmt, mm_val),   collapse = ", "),
              diff, tol, status))
  invisible(diff)
}

## fixef
fe_lme4 <- lme4::fixef(fit_lme4)
fe_mm   <- fixef(fit_mm)
compare("fixef (Intercept)", fe_lme4[["(Intercept)"]], fe_mm[["(Intercept)"]], tols$fixef)
compare("fixef x",           fe_lme4[["x"]],           fe_mm[["x"]],           tols$fixef)

## SE
se_lme4 <- sqrt(diag(vcov(fit_lme4)))
se_mm   <- fit_mm$std_errors
compare("SE (Intercept)", se_lme4[["(Intercept)"]], se_mm[["(Intercept)"]], tols$fixef)
compare("SE x",           se_lme4[["x"]],           se_mm[["x"]],           tols$fixef)

## vcov diagonal
compare("vcov[1,1]",
        as.numeric(vcov(fit_lme4)[1, 1]),
        as.numeric(fit_mm$fixed_effect_vcov[1, 1]),
        tols$fixef^2)
compare("vcov[2,2]",
        as.numeric(vcov(fit_lme4)[2, 2]),
        as.numeric(fit_mm$fixed_effect_vcov[2, 2]),
        tols$fixef^2)

## theta
theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
cat(sprintf("%-36s lme4=%s  mm=%s\n", "theta (raw)",
            paste(sprintf("%.8f", theta_lme4), collapse = ", "),
            paste(sprintf("%.8f", theta_mm),   collapse = ", ")))
compare("theta", theta_lme4, theta_mm, tols$theta)

## sigma
compare("sigma", stats::sigma(fit_lme4), sigma(fit_mm), tols$sigma)

## VarCorr: group RE variance
vc_lme4 <- as.numeric(attr(lme4::VarCorr(fit_lme4)$g, "stddev"))^2
vc_mm_df <- VarCorr(fit_mm)
vc_mm <- tryCatch({
  if (is.data.frame(vc_mm_df)) {
    as.numeric(vc_mm_df[vc_mm_df$grp == "g", "vcov"])
  } else if (is.list(vc_mm_df)) {
    v <- vc_mm_df[["g"]]
    if (is.numeric(v)) v
    else if (is.matrix(v)) v[1, 1]
    else NA_real_
  } else {
    NA_real_
  }
}, error = function(e) {
  cat("VarCorr extraction error:", conditionMessage(e), "\n")
  NA_real_
})
cat(sprintf("%-36s lme4=%.8f  mm=%s\n", "VarCorr g var",
            vc_lme4,
            ifelse(is.na(vc_mm), "NA", sprintf("%.8f", vc_mm))))
if (!is.na(vc_mm)) compare("VarCorr g var", vc_lme4, vc_mm, tols$theta)

## logLik
compare("logLik",
        as.numeric(stats::logLik(fit_lme4)),
        as.numeric(logLik(fit_mm)),
        tols$logLik)

## AIC
compare("AIC", stats::AIC(fit_lme4), fit_mm$AIC, tols$logLik * 2)

## BIC
compare("BIC", stats::BIC(fit_lme4), fit_mm$BIC, tols$logLik * 2)

## fitted values
fitted_lme4 <- as.numeric(stats::fitted(fit_lme4))
fitted_mm   <- as.numeric(fit_mm$fitted)
compare("fitted max abs diff",
        fitted_lme4,
        fitted_mm[seq_along(fitted_lme4)],
        tols$fixef)

## ranef (sorted for order-invariant comparison)
re_lme4 <- as.numeric(lme4::ranef(fit_lme4)$g[, 1])
re_mm_df <- ranef(fit_mm)$g
re_mm    <- if (!is.null(re_mm_df)) as.numeric(re_mm_df[, 1]) else NA_real_
cat(sprintf("\n%-36s\n", "ranef g (sorted):"))
cat("  lme4:", sprintf("%.5f", sort(re_lme4)), "\n")
cat("  mm:  ", sprintf("%.5f", sort(re_mm)),   "\n")
compare("ranef g max abs diff", sort(re_lme4), sort(re_mm), tols$fixef)

## ── 6. Speed comparison ───────────────────────────────────────────────────────
cat("\n=== SPEED COMPARISON ===\n")
NREPS <- 5
t_lme4_rep <- system.time(
  for (i in seq_len(NREPS))
    lmer(y ~ x + (1 | g), data = dat, weights = dat$w, REML = TRUE)
)
t_mm_rep <- system.time(
  for (i in seq_len(NREPS))
    lmm(y ~ x + (1 | g), data = dat, weights = dat$w, REML = TRUE,
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
cat("Key check: do weighted and unweighted fits differ (weights honored)?\n")
