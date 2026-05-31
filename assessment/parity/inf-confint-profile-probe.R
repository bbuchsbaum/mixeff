#!/usr/bin/env Rscript
# Parity probe: confint(method = "profile") — cell "inf-confint-profile"
# Dataset: sleepstudy (lme4), Formula: Reaction ~ Days + (1|Subject)
# (simpler formula than (1+Days|Subject) as cell spec states)
# Compare lme4::confint(method="profile") vs mixeff confint(method="profile")

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

cat("=== inf-confint-profile parity probe ===\n\n")
cat("R version:", R.version$version.string, "\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

data("sleepstudy", package = "lme4")

formula_str <- "Reaction ~ Days + (1 | Subject)"
formula_obj <- Reaction ~ Days + (1 | Subject)

cat("Formula:", formula_str, "\n")
cat("Dataset: sleepstudy (", nrow(sleepstudy), "rows )\n\n")

# ─── 1. Fit models ───────────────────────────────────────────────────────────
cat("--- Fitting models ---\n")

t_lme4_fit <- system.time({
  fit_lme4_ml   <- lmer(formula_obj, data = sleepstudy, REML = FALSE)
  fit_lme4_reml <- lmer(formula_obj, data = sleepstudy, REML = TRUE)
})

t_mm_fit <- system.time({
  fit_mm_ml   <- lmm(formula_obj, data = sleepstudy, REML = FALSE,
                     control = mm_control(verbose = -1))
  fit_mm_reml <- lmm(formula_obj, data = sleepstudy, REML = TRUE,
                     control = mm_control(verbose = -1))
})

cat(sprintf("lme4   fit time: %.3f s\n", t_lme4_fit["elapsed"]))
cat(sprintf("mixeff fit time: %.3f s\n", t_mm_fit["elapsed"]))

# ─── 2. Core fit quantities comparison ────────────────────────────────────────
cat("\n--- Fixed effects (ML) ---\n")
lme4_fixef <- lme4::fixef(fit_lme4_ml)
mm_fixef   <- mixeff::fixef(fit_mm_ml)
cat("lme4  fixef:", paste(names(lme4_fixef), round(lme4_fixef, 6), sep = "=", collapse = ", "), "\n")
cat("mm    fixef:", paste(names(mm_fixef),   round(mm_fixef,   6), sep = "=", collapse = ", "), "\n")
fixef_diff <- abs(lme4_fixef - mm_fixef[names(lme4_fixef)])
cat("Max |diff| fixef:", max(fixef_diff, na.rm = TRUE), " (tol 1e-4)\n")

cat("\n--- sigma (ML) ---\n")
lme4_sigma <- sigma(fit_lme4_ml)
mm_sigma   <- fit_mm_ml$sigma
if (is.null(mm_sigma)) mm_sigma <- NA_real_
cat(sprintf("lme4  sigma: %.8f\n", lme4_sigma))
cat(sprintf("mm    sigma: %.8f\n", mm_sigma))
cat(sprintf("Max |diff| sigma: %g  (tol 1e-4)\n", abs(lme4_sigma - mm_sigma)))

cat("\n--- logLik (ML) ---\n")
lme4_ll <- as.numeric(stats::logLik(fit_lme4_ml))
mm_ll   <- as.numeric(stats::logLik(fit_mm_ml))
cat(sprintf("lme4  logLik: %.8f\n", lme4_ll))
cat(sprintf("mm    logLik: %.8f\n", mm_ll))
cat(sprintf("Max |diff| logLik: %g  (tol 1e-3)\n", abs(lme4_ll - mm_ll)))

cat("\n--- theta / VarCorr (ML) ---\n")
lme4_vc <- lme4::VarCorr(fit_lme4_ml)
lme4_theta <- lme4::getME(fit_lme4_ml, "theta")
mm_vc      <- mixeff::VarCorr(fit_mm_ml)
cat("lme4  theta:", round(lme4_theta, 8), "\n")
# mixeff theta stored in artifact or accessible via getME-like accessor
mm_theta_raw <- tryCatch(
  fit_mm_ml$theta,
  error = function(e) NA_real_
)
if (is.null(mm_theta_raw)) mm_theta_raw <- NA_real_
cat("mm    theta:", round(mm_theta_raw, 8), "\n")
if (length(lme4_theta) == length(mm_theta_raw) && !any(is.na(mm_theta_raw))) {
  theta_diff <- abs(lme4_theta - mm_theta_raw)
  cat("Max |diff| theta:", max(theta_diff), " (tol 1e-3)\n")
} else {
  cat("theta comparison: cannot compare (lengths differ or NA)\n")
}

# ─── 3. Profile CI — ML fit ───────────────────────────────────────────────────
cat("\n\n=== PROFILE CI: ML fit ===\n")

cat("\n--- lme4 profile CI (ML) ---\n")
t_lme4_prof_ml <- system.time({
  lme4_prof_ml <- tryCatch(
    suppressMessages(confint(fit_lme4_ml, method = "profile", level = 0.95)),
    error = function(e) { cat("lme4 profile ML ERROR:", conditionMessage(e), "\n"); NULL }
  )
})
cat(sprintf("lme4 profile CI (ML) time: %.3f s\n", t_lme4_prof_ml["elapsed"]))
if (!is.null(lme4_prof_ml)) {
  cat("lme4 profile CI (ML):\n")
  print(lme4_prof_ml)
}

cat("\n--- mixeff profile CI (ML) ---\n")
t_mm_prof_ml <- system.time({
  mm_prof_ml <- tryCatch(
    confint(fit_mm_ml, method = "profile", level = 0.95),
    error = function(e) { cat("mixeff profile ML ERROR:", conditionMessage(e), "\n"); e }
  )
})
cat(sprintf("mixeff profile CI (ML) time: %.3f s\n", t_mm_prof_ml["elapsed"]))
if (inherits(mm_prof_ml, "condition")) {
  cat("MIXEFF PROFILE CI ML: ERROR (shown above)\n")
  cat("Error class:", paste(class(mm_prof_ml), collapse = "/"), "\n")
} else {
  cat("mixeff profile CI (ML):\n")
  print(mm_prof_ml)
  cat("\nmm_profile payload structure:\n")
  payload <- attr(mm_prof_ml, "mm_profile")
  if (!is.null(payload)) {
    cat("  schema:", payload$schema$schema_name, payload$schema$schema_version, "\n")
    cat("  fit_criterion:", payload$fit_criterion, "\n")
    cat("  level:", payload$level, "\n")
    cat("  table:\n")
    print(payload$table)
  }
}

# ─── 4. Compare ML profile CI rows for beta ───────────────────────────────────
cat("\n\n=== COMPARISON: ML beta profile CI ===\n")
if (!is.null(lme4_prof_ml) && !inherits(mm_prof_ml, "condition")) {
  # lme4 uses ".beta." prefix and ".sig01", ".sigma" etc.
  # Map lme4 rows to mixeff parameter names
  lme4_rows <- rownames(lme4_prof_ml)
  cat("lme4 CI rownames:", paste(lme4_rows, collapse = ", "), "\n")
  cat("mm   CI rownames:", paste(rownames(mm_prof_ml), collapse = ", "), "\n\n")

  # Compare fixef rows
  beta_map <- list(
    "(Intercept)" = "(Intercept)",
    "Days"        = "Days"
  )
  for (mm_nm in names(beta_map)) {
    lme4_nm <- beta_map[[mm_nm]]
    # lme4 profile uses ".fixed." notation: check
    lme4_row <- if (lme4_nm %in% lme4_rows) lme4_nm else
      grep(paste0("^", lme4_nm, "$"), lme4_rows, value = TRUE)[1]
    if (!is.na(lme4_row) && nzchar(lme4_row) && mm_nm %in% rownames(mm_prof_ml)) {
      lme4_lo <- lme4_prof_ml[lme4_row, 1]
      lme4_hi <- lme4_prof_ml[lme4_row, 2]
      mm_lo   <- mm_prof_ml[mm_nm, 1]
      mm_hi   <- mm_prof_ml[mm_nm, 2]
      cat(sprintf("%-20s  lme4: [%.5f, %.5f]  mm: [%.5f, %.5f]  |diff_lo|=%g  |diff_hi|=%g\n",
                  mm_nm, lme4_lo, lme4_hi, mm_lo, mm_hi,
                  abs(lme4_lo - mm_lo), abs(lme4_hi - mm_hi)))
    } else {
      cat(sprintf("%-20s  lme4 row '%s' or mm row '%s' not found in CI matrix\n",
                  mm_nm, lme4_nm, mm_nm))
    }
  }

  # Compare sigma
  lme4_sigma_row <- grep("sigma", lme4_rows, ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(lme4_sigma_row) && "sigma" %in% rownames(mm_prof_ml)) {
    lme4_s_lo <- lme4_prof_ml[lme4_sigma_row, 1]
    lme4_s_hi <- lme4_prof_ml[lme4_sigma_row, 2]
    mm_s_lo   <- mm_prof_ml["sigma", 1]
    mm_s_hi   <- mm_prof_ml["sigma", 2]
    cat(sprintf("%-20s  lme4: [%.5f, %.5f]  mm: [%.5f, %.5f]  |diff_lo|=%g  |diff_hi|=%g\n",
                "sigma", lme4_s_lo, lme4_s_hi, mm_s_lo, mm_s_hi,
                abs(lme4_s_lo - mm_s_lo), abs(lme4_s_hi - mm_s_hi)))
  }

  # Compare theta
  lme4_theta_rows <- grep("theta|sd_|cor_|sig0", lme4_rows, ignore.case = TRUE, value = TRUE)
  mm_theta_rows   <- rownames(mm_prof_ml)[grepl("^theta", rownames(mm_prof_ml))]
  cat(sprintf("\nlme4 theta rows: %s\n", paste(lme4_theta_rows, collapse = ", ")))
  cat(sprintf("mm   theta rows: %s\n", paste(mm_theta_rows, collapse = ", ")))

  # Compute overall max abs diff for beta rows
  diffs_lo <- numeric()
  diffs_hi <- numeric()
  for (mm_nm in names(beta_map)) {
    lme4_nm <- beta_map[[mm_nm]]
    lme4_row <- if (lme4_nm %in% lme4_rows) lme4_nm else
      grep(paste0("^", lme4_nm, "$"), lme4_rows, value = TRUE)[1]
    if (!is.na(lme4_row) && nzchar(lme4_row) && mm_nm %in% rownames(mm_prof_ml)) {
      diffs_lo <- c(diffs_lo, abs(lme4_prof_ml[lme4_row, 1] - mm_prof_ml[mm_nm, 1]))
      diffs_hi <- c(diffs_hi, abs(lme4_prof_ml[lme4_row, 2] - mm_prof_ml[mm_nm, 2]))
    }
  }
  cat(sprintf("\nMax |diff| profile CI lower (beta): %g\n", max(diffs_lo, na.rm = TRUE)))
  cat(sprintf("Max |diff| profile CI upper (beta): %g\n", max(diffs_hi, na.rm = TRUE)))
  cat(sprintf("Max |diff| profile CI beta overall: %g\n",
              max(c(diffs_lo, diffs_hi), na.rm = TRUE)))
} else {
  cat("CANNOT COMPARE: one or both profile CIs failed.\n")
}

# ─── 5. Profile CI — REML fit ─────────────────────────────────────────────────
cat("\n\n=== PROFILE CI: REML fit ===\n")

cat("\n--- lme4 profile CI (REML) ---\n")
t_lme4_prof_reml <- system.time({
  lme4_prof_reml <- tryCatch(
    suppressMessages(confint(fit_lme4_reml, method = "profile", level = 0.95)),
    error = function(e) { cat("lme4 profile REML ERROR:", conditionMessage(e), "\n"); NULL }
  )
})
cat(sprintf("lme4 profile CI (REML) time: %.3f s\n", t_lme4_prof_reml["elapsed"]))
if (!is.null(lme4_prof_reml)) {
  cat("lme4 profile CI (REML):\n")
  print(lme4_prof_reml)
}

cat("\n--- mixeff profile CI (REML) ---\n")
t_mm_prof_reml <- system.time({
  mm_prof_reml <- tryCatch(
    confint(fit_mm_reml, method = "profile", level = 0.95),
    error = function(e) { cat("mixeff profile REML ERROR:", conditionMessage(e), "\n"); e }
  )
})
cat(sprintf("mixeff profile CI (REML) time: %.3f s\n", t_mm_prof_reml["elapsed"]))
if (inherits(mm_prof_reml, "condition")) {
  cat("MIXEFF PROFILE CI REML: ERROR (shown above)\n")
  cat("Error class:", paste(class(mm_prof_reml), collapse = "/"), "\n")
} else {
  cat("mixeff profile CI (REML):\n")
  print(mm_prof_reml)
  payload_r <- attr(mm_prof_reml, "mm_profile")
  if (!is.null(payload_r)) {
    cat("\nREML profile table:\n")
    print(payload_r$table)
  }
}

# ─── 6. Profile CI parm subsetting ───────────────────────────────────────────
cat("\n\n=== PROFILE CI: parm subsetting ===\n")
mm_prof_intercept <- tryCatch(
  confint(fit_mm_ml, parm = "(Intercept)", method = "profile", level = 0.95),
  error = function(e) { cat("mixeff parm subset ERROR:", conditionMessage(e), "\n"); e }
)
if (!inherits(mm_prof_intercept, "condition")) {
  cat("confint(parm='(Intercept)', method='profile') rownames:",
      paste(rownames(mm_prof_intercept), collapse = ", "), "\n")
  cat("Expected: only '(Intercept)'\n")
  cat("Pass:", identical(rownames(mm_prof_intercept), "(Intercept)"), "\n")
}

# ─── 7. Boundary fit (Dyestuff2) ─────────────────────────────────────────────
cat("\n\n=== PROFILE CI: boundary fit (Dyestuff2) ===\n")
data("Dyestuff2", package = "lme4")
fit_dy2_lme4 <- lmer(Yield ~ 1 + (1 | Batch), data = Dyestuff2, REML = TRUE)
fit_dy2_mm   <- tryCatch(
  lmm(Yield ~ 1 + (1 | Batch), data = Dyestuff2, REML = TRUE,
      control = mm_control(verbose = -1)),
  error = function(e) { cat("mixeff Dyestuff2 fit ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(fit_dy2_mm)) {
  mm_dy2_prof <- tryCatch(
    confint(fit_dy2_mm, method = "profile", level = 0.95),
    error = function(e) { cat("mixeff Dyestuff2 profile ERROR:", conditionMessage(e), "\n"); e }
  )
  if (inherits(mm_dy2_prof, "condition")) {
    cat("Boundary fit raises a condition (typed refusal or error):\n")
    cat("  class:", paste(class(mm_dy2_prof), collapse = "/"), "\n")
    cat("  message:", conditionMessage(mm_dy2_prof), "\n")
    typed <- any(class(mm_dy2_prof) %in%
                 c("mm_inference_unavailable", "mm_schema_error",
                   "mm_bridge_error", "mm_fit_error"))
    cat("  Is typed mixeff error:", typed, "\n")
  } else {
    cat("Boundary fit returned a CI (not an error):\n")
    print(mm_dy2_prof)
    payload_dy <- attr(mm_dy2_prof, "mm_profile")
    if (!is.null(payload_dy)) {
      cat("\nBoundary profile table:\n")
      print(payload_dy$table)
    }
  }
}

# ─── 8. Speed comparison ─────────────────────────────────────────────────────
cat("\n\n=== TIMING SUMMARY ===\n")
cat(sprintf("lme4   ML profile CI:   %.3f s\n", t_lme4_prof_ml["elapsed"]))
cat(sprintf("mixeff ML profile CI:   %.3f s\n", t_mm_prof_ml["elapsed"]))
if (t_lme4_prof_ml["elapsed"] > 0 && t_mm_prof_ml["elapsed"] > 0) {
  cat(sprintf("Speed ratio (lme4/mixeff): %.2fx\n",
              t_lme4_prof_ml["elapsed"] / t_mm_prof_ml["elapsed"]))
}
cat(sprintf("lme4   REML profile CI: %.3f s\n", t_lme4_prof_reml["elapsed"]))
cat(sprintf("mixeff REML profile CI: %.3f s\n", t_mm_prof_reml["elapsed"]))

cat("\n=== DONE ===\n")
