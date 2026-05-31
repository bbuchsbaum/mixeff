## Empirical parity probe: speed-scaling
## Simulated data, increasing N (1e3, 1e4, 5e4) with 50-200 groups.
## Formula: y ~ x1 + x2 + (x1|g)
## Focus: wall-clock timing ratio at each size, plus full numerical parity check.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:",    as.character(packageVersion("lme4")),    "\n")
cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff version:",  as.character(packageVersion("mixeff")),  "\n\n")

## ── Helpers ──────────────────────────────────────────────────────────────────
set.seed(42)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  diff     <- max(abs(lme4_val - mm_val), na.rm = TRUE)
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-38s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse = ", "),
              paste(sprintf(fmt, mm_val),   collapse = ", "),
              diff, tol, status))
  invisible(diff)
}

sim_data <- function(N, n_groups) {
  g  <- factor(sample(seq_len(n_groups), N, replace = TRUE))
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  # True params: fixef = c(2, 0.5, -0.3), RE sd = 1.2 (intercept), 0.4 (x1), corr=0.3
  re_int  <- rnorm(n_groups, 0, 1.2)[as.integer(g)]
  re_x1   <- rnorm(n_groups, 0, 0.4)[as.integer(g)]
  y <- 2 + 0.5 * x1 + (-0.3) * x2 + re_int + re_x1 * x1 + rnorm(N, 0, 0.8)
  data.frame(y = y, x1 = x1, x2 = x2, g = g)
}

tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

## ── Size configurations ───────────────────────────────────────────────────────
configs <- list(
  list(N =  1000, n_groups =  50, label = "N=1e3,  G=50"),
  list(N = 10000, n_groups = 100, label = "N=1e4, G=100"),
  list(N = 50000, n_groups = 200, label = "N=5e4, G=200")
)

NREPS <- 3   # reps for timing (small to keep runtime reasonable)

all_results <- list()

for (cfg in configs) {
  N         <- cfg$N
  n_groups  <- cfg$n_groups
  label     <- cfg$label
  cat(sprintf("\n\n%s\n", strrep("=", 70)))
  cat(sprintf("=== SIZE: %s  (NREPS=%d for timing) ===\n", label, NREPS))
  cat(sprintf("%s\n", strrep("=", 70)))

  ## Simulate once for parity, then re-simulate inside loops for timing
  df <- sim_data(N, n_groups)
  cat(sprintf("Dataset: nrow=%d, groups=%d\n\n", nrow(df), nlevels(df$g)))

  ## ── lme4 single fit (for parity numbers) ─────────────────────────────────
  cat("--- lme4 single fit ---\n")
  t_lme4_single <- system.time({
    fit_lme4 <- tryCatch(
      lmer(y ~ x1 + x2 + (x1 | g), data = df, REML = TRUE),
      error = function(e) e
    )
  })

  if (inherits(fit_lme4, "condition")) {
    cat("!!! lme4 ERROR !!!\n")
    cat(conditionMessage(fit_lme4), "\n")
    next
  }
  cat(sprintf("lme4 single fit wall-clock: %.4f s\n", t_lme4_single["elapsed"]))
  conv_lme4 <- length(fit_lme4@optinfo$conv$lme4$messages) == 0
  cat(sprintf("lme4 convergence clean: %s\n", conv_lme4))

  ## ── mixeff single fit (for parity numbers) ────────────────────────────────
  cat("\n--- mixeff single fit ---\n")
  t_mm_single <- system.time({
    fit_mm <- tryCatch(
      lmm(y ~ x1 + x2 + (x1 | g), data = df, REML = TRUE,
          control = mm_control(verbose = -1L)),
      error = function(e) e
    )
  })

  mm_error <- inherits(fit_mm, "condition")
  if (mm_error) {
    cat("!!! mixeff ERROR !!!\n")
    cat(conditionMessage(fit_mm), "\n")
    cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  } else {
    cat(sprintf("mixeff single fit wall-clock: %.4f s\n", t_mm_single["elapsed"]))
    cat(sprintf("fit_status: %s\n", fit_mm$fit_status))
  }

  ## ── Numerical parity (only if both fits succeeded) ────────────────────────
  if (!mm_error) {
    cat(sprintf("\n--- Numerical parity (%s) ---\n", label))

    fe_lme4 <- lme4::fixef(fit_lme4)
    fe_mm   <- fixef(fit_mm)
    d_int  <- compare("fixef (Intercept)", fe_lme4["(Intercept)"], fe_mm["(Intercept)"], tols$fixef)
    d_x1   <- compare("fixef x1",          fe_lme4["x1"],           fe_mm["x1"],          tols$fixef)
    d_x2   <- compare("fixef x2",          fe_lme4["x2"],           fe_mm["x2"],          tols$fixef)

    se_lme4 <- sqrt(diag(vcov(fit_lme4)))
    se_mm   <- fit_mm$std_errors
    d_se_int <- compare("SE (Intercept)", se_lme4["(Intercept)"], se_mm["(Intercept)"], tols$fixef)
    d_se_x1  <- compare("SE x1",          se_lme4["x1"],           se_mm["x1"],          tols$fixef)
    d_se_x2  <- compare("SE x2",          se_lme4["x2"],           se_mm["x2"],          tols$fixef)

    theta_lme4 <- lme4::getME(fit_lme4, "theta")
    theta_mm   <- fit_mm$theta
    cat(sprintf("%-38s lme4=%s  mm=%s\n", "theta (raw)",
                paste(sprintf("%.8f", theta_lme4), collapse = ", "),
                paste(sprintf("%.8f", theta_mm),   collapse = ", ")))
    d_theta <- compare("theta", theta_lme4, theta_mm, tols$theta)

    d_sigma  <- compare("sigma",  stats::sigma(fit_lme4), sigma(fit_mm), tols$sigma)
    d_loglik <- compare("logLik", as.numeric(stats::logLik(fit_lme4)), as.numeric(logLik(fit_mm)), tols$logLik)
    d_aic    <- compare("AIC",    stats::AIC(fit_lme4), fit_mm$AIC, tols$logLik * 2)
    d_bic    <- compare("BIC",    stats::BIC(fit_lme4), fit_mm$BIC, tols$logLik * 2)

    fitted_lme4 <- as.numeric(stats::fitted(fit_lme4))
    fitted_mm   <- as.numeric(fit_mm$fitted)
    d_fitted <- compare("fitted max abs diff", fitted_lme4, fitted_mm[seq_along(fitted_lme4)], tols$fixef)

    ## ranef comparison — match by group level name
    re_lme4_df <- lme4::ranef(fit_lme4)$g
    re_mm_list <- ranef(fit_mm)
    re_mm_df   <- re_mm_list$g
    if (!is.null(re_lme4_df) && !is.null(re_mm_df)) {
      re_lme4_int <- re_lme4_df[, "(Intercept)"]
      re_lme4_x1  <- re_lme4_df[, "x1"]
      re_mm_int   <- re_mm_df[, "(Intercept)"]
      re_mm_x1    <- re_mm_df[, "x1"]
      d_re_int <- compare("ranef (Intercept) max abs diff", sort(re_lme4_int), sort(re_mm_int), tols$fixef)
      d_re_x1  <- compare("ranef x1 max abs diff",          sort(re_lme4_x1),  sort(re_mm_x1),  tols$fixef)
    } else {
      cat("ranef extraction: one or both returned NULL\n")
      d_re_int <- NA; d_re_x1 <- NA
    }

    max_diff <- max(d_int, d_x1, d_x2, d_sigma, d_loglik, na.rm = TRUE)
    cat(sprintf("\nMax abs diff across key quantities (fixef+sigma+logLik): %.3e\n", max_diff))
  }

  ## ── Timing benchmark (multiple reps) ─────────────────────────────────────
  cat(sprintf("\n--- Timing benchmark (%s, %d reps each) ---\n", label, NREPS))

  t_lme4_rep <- system.time(
    for (i in seq_len(NREPS)) {
      df_i <- sim_data(N, n_groups)
      lmer(y ~ x1 + x2 + (x1 | g), data = df_i, REML = TRUE)
    }
  )

  t_mm_rep <- system.time(
    for (i in seq_len(NREPS)) {
      df_i <- sim_data(N, n_groups)
      lmm(y ~ x1 + x2 + (x1 | g), data = df_i, REML = TRUE,
          control = mm_control(verbose = -1L))
    }
  )

  lme4_per <- t_lme4_rep["elapsed"] / NREPS
  mm_per   <- t_mm_rep["elapsed"]   / NREPS
  ratio    <- mm_per / lme4_per

  cat(sprintf("lme4  mean/fit: %.4f s  (over %d reps, includes data sim)\n", lme4_per, NREPS))
  cat(sprintf("mm    mean/fit: %.4f s  (over %d reps, includes data sim)\n", mm_per,   NREPS))
  cat(sprintf("ratio (mm/lme4): %.3fx  %s\n",
              ratio,
              if (ratio < 1) "(mixeff FASTER)" else if (ratio <= 1.5) "(roughly equal)" else "(mixeff SLOWER)"))

  all_results[[label]] <- list(
    N = N, n_groups = n_groups,
    lme4_per = lme4_per, mm_per = mm_per, ratio = ratio,
    mm_error = mm_error
  )
}

## ── Summary table ────────────────────────────────────────────────────────────
cat(sprintf("\n\n%s\n", strrep("=", 70)))
cat("=== SPEED SCALING SUMMARY TABLE ===\n")
cat(sprintf("%s\n", strrep("=", 70)))
cat(sprintf("%-20s  %8s  %8s  %8s  %10s\n",
            "Config", "lme4(s)", "mm(s)", "ratio", "verdict"))
cat(strrep("-", 60), "\n")
for (lbl in names(all_results)) {
  r <- all_results[[lbl]]
  if (r$mm_error) {
    verdict <- "mm-ERROR"
    cat(sprintf("%-20s  %8.4f  %8s  %8s  %10s\n",
                lbl, r$lme4_per, "ERROR", "N/A", verdict))
  } else {
    verdict <- if (r$ratio < 0.8) "FASTER" else if (r$ratio <= 1.5) "COMPARABLE" else "SLOWER"
    cat(sprintf("%-20s  %8.4f  %8.4f  %8.3fx  %10s\n",
                lbl, r$lme4_per, r$mm_per, r$ratio, verdict))
  }
}
cat(strrep("-", 60), "\n")
cat("\nNote: timing includes data simulation overhead (negligible vs fit time at large N).\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
cat("\nDone.\n")
