## Empirical parity probe: lmm-penicillin-crossed
## Dataset: Penicillin  Formula: diameter ~ 1 + (1|plate) + (1|sample)
## Focus: fully crossed random effects (two independent random intercepts)
## Compares lme4/lmerTest vs mixeff on all relevant quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

## ── 1. Data ──────────────────────────────────────────────────────────────────
data(Penicillin, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(Penicillin), "  ncol:", ncol(Penicillin), "\n")
cat("plate levels:", nlevels(Penicillin$plate), "\n")
cat("sample levels:", nlevels(Penicillin$sample), "\n")
print(head(Penicillin))
cat("\n")

## ── 2. Fit lme4 ──────────────────────────────────────────────────────────────
cat("=== lme4 FIT ===\n")
N_REPS <- 5L
t_lme4 <- system.time({
  for (i in seq_len(N_REPS)) {
    fit_lme4 <- lmer(diameter ~ 1 + (1 | plate) + (1 | sample),
                     data = Penicillin, REML = TRUE)
  }
})
cat("lme4 wall-clock (seconds, mean of", N_REPS, "reps):", t_lme4["elapsed"] / N_REPS, "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov --\n"); print(as.matrix(vcov(fit_lme4)))
cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))  # explicit namespace
cat("-- sigma --\n"); print(sigma(fit_lme4))
cat("-- logLik --\n"); print(logLik(fit_lme4))
cat("-- AIC --\n");    print(AIC(fit_lme4))
cat("-- BIC --\n");    print(BIC(fit_lme4))
cat("-- theta (getME) --\n"); print(lme4::getME(fit_lme4, "theta"))
cat("-- ranef plate (first 6) --\n"); print(head(lme4::ranef(fit_lme4)$plate))
cat("-- ranef sample (all) --\n"); print(lme4::ranef(fit_lme4)$sample)
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
conv_msgs <- fit_lme4@optinfo$conv$lme4$messages
cat("-- convergence --\n")
cat("No convergence warning =", length(conv_msgs) == 0, "\n")
if (length(conv_msgs) > 0) cat("Messages:", paste(conv_msgs, collapse="; "), "\n")
cat("\n")

## ── 3. Fit mixeff ─────────────────────────────────────────────────────────────
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  for (i in seq_len(N_REPS)) {
    fit_mm <- tryCatch(
      lmm(diameter ~ 1 + (1 | plate) + (1 | sample),
          data = Penicillin, REML = TRUE,
          control = mm_control(verbose = -1L)),
      error = function(e) e
    )
  }
})
cat("mixeff wall-clock (seconds, mean of", N_REPS, "reps):", t_mixeff["elapsed"] / N_REPS, "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR !!!\n")
  cat(conditionMessage(fit_mm), "\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  quit(status = 0)
}

cat("fit_status:", fit_mm$fit_status, "\n\n")

cat("-- fixef --\n"); print(fixef(fit_mm))
cat("-- SE --\n"); print(fit_mm$std_errors)
cat("-- vcov --\n"); print(as.matrix(fit_mm$fixed_effect_vcov))
cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
cat("-- sigma --\n"); print(sigma(fit_mm))
cat("-- logLik --\n"); print(logLik(fit_mm))
cat("-- AIC --\n");    print(fit_mm$AIC)
cat("-- BIC --\n");    print(fit_mm$BIC)
cat("-- theta --\n"); print(fit_mm$theta)
cat("-- ranef plate (first 6) --\n"); print(head(ranef(fit_mm)$plate))
cat("-- ranef sample (all) --\n"); print(ranef(fit_mm)$sample)
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("\n")

## ── 4. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  if (any(is.na(mm_val))) {
    cat(sprintf("%-35s lme4=%s  mm=NA  [MISSING]\n",
                label,
                paste(sprintf(fmt, lme4_val), collapse=", ")))
    return(invisible(NA_real_))
  }
  diff     <- max(abs(lme4_val - mm_val))
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-35s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse=", "),
              paste(sprintf(fmt, mm_val),   collapse=", "),
              diff, tol, status))
  invisible(diff)
}

## fixef
compare("fixef (Intercept)",
        lme4::fixef(fit_lme4)[["(Intercept)"]],
        fixef(fit_mm)[["(Intercept)"]],
        tols$fixef)

## SE
compare("SE (Intercept)",
        sqrt(diag(vcov(fit_lme4)))[["(Intercept)"]],
        fit_mm$std_errors[["(Intercept)"]],
        tols$fixef)

## vcov diagonal
compare("vcov[1,1]",
        as.numeric(vcov(fit_lme4)[1,1]),
        as.numeric(fit_mm$fixed_effect_vcov[1,1]),
        tols$fixef^2)

## theta: lme4 returns theta in order of grouping factors (alphabetical)
## For this model: (1|plate) and (1|sample), each theta is sd_re/sigma
theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
cat(sprintf("%-35s lme4=%s\n", "theta (lme4 named)",
            paste(sprintf("%s=%.8f", names(theta_lme4), theta_lme4), collapse=", ")))
cat(sprintf("%-35s mm=%s\n",   "theta (mm)",
            paste(sprintf("%.8f", theta_mm), collapse=", ")))

## Match theta by position (lme4 orders alphabetically: plate, sample)
if (length(theta_mm) == length(theta_lme4)) {
  compare("theta (all)", theta_lme4, theta_mm, tols$theta)
} else {
  cat(sprintf("theta length mismatch: lme4=%d mm=%d\n", length(theta_lme4), length(theta_mm)))
}

## sigma
compare("sigma",
        sigma(fit_lme4),
        sigma(fit_mm),
        tols$sigma)

## VarCorr: RE variances for plate and sample
vc_lme4_plate  <- as.numeric(attr(lme4::VarCorr(fit_lme4)$plate,  "stddev"))^2
vc_lme4_sample <- as.numeric(attr(lme4::VarCorr(fit_lme4)$sample, "stddev"))^2

vc_mm_obj <- VarCorr(fit_mm)
extract_vc_mm <- function(obj, grp_name) {
  tryCatch({
    # mixeff VarCorr is class mm_varcorr: list with $table (data.frame) and $residual_sd
    if (inherits(obj, "mm_varcorr")) {
      tbl <- obj$table
      v <- tbl[tbl$group == grp_name, "variance"]
      if (length(v) == 0) NA_real_ else as.numeric(v[1])
    } else if (is.data.frame(obj)) {
      # lme4-style data.frame with grp/vcov columns
      v <- obj[obj$grp == grp_name, "vcov"]
      if (length(v) == 0) NA_real_ else as.numeric(v[1])
    } else if (is.list(obj)) {
      v <- obj[[grp_name]]
      if (is.null(v)) NA_real_
      else if (is.matrix(v)) v[1,1]
      else as.numeric(v)
    } else NA_real_
  }, error = function(e) { cat("VarCorr extract error:", conditionMessage(e), "\n"); NA_real_ })
}

vc_mm_plate  <- extract_vc_mm(vc_mm_obj, "plate")
vc_mm_sample <- extract_vc_mm(vc_mm_obj, "sample")

cat(sprintf("%-35s lme4=%.8f  mm=%s\n", "VarCorr plate var",
            vc_lme4_plate,
            ifelse(is.na(vc_mm_plate), "NA", sprintf("%.8f", vc_mm_plate))))
cat(sprintf("%-35s lme4=%.8f  mm=%s\n", "VarCorr sample var",
            vc_lme4_sample,
            ifelse(is.na(vc_mm_sample), "NA", sprintf("%.8f", vc_mm_sample))))

if (!is.na(vc_mm_plate))  compare("VarCorr plate var",  vc_lme4_plate,  vc_mm_plate,  tols$theta)
if (!is.na(vc_mm_sample)) compare("VarCorr sample var", vc_lme4_sample, vc_mm_sample, tols$theta)

## logLik
compare("logLik",
        as.numeric(logLik(fit_lme4)),
        as.numeric(logLik(fit_mm)),
        tols$logLik)

## AIC / BIC
compare("AIC",
        AIC(fit_lme4),
        fit_mm$AIC,
        tols$logLik * 2)

compare("BIC",
        BIC(fit_lme4),
        fit_mm$BIC,
        tols$logLik * 2)

## fitted values (matched by row order)
fitted_lme4 <- fitted(fit_lme4)
fitted_mm   <- fit_mm$fitted
compare("fitted max abs diff",
        as.numeric(fitted_lme4),
        as.numeric(fitted_mm)[seq_along(fitted_lme4)],
        tols$fixef)

## ranef plate
re_lme4_plate <- as.numeric(lme4::ranef(fit_lme4)$plate[, 1])
re_mm_plate_df <- tryCatch(ranef(fit_mm)$plate, error = function(e) NULL)
re_mm_plate <- if (!is.null(re_mm_plate_df)) as.numeric(re_mm_plate_df[, 1]) else NA_real_

cat(sprintf("\n%-35s\n", "ranef plate (sorted, first 6):"))
cat("  lme4:", sprintf("%.6f", sort(re_lme4_plate)[1:6]), "\n")
if (!all(is.na(re_mm_plate))) cat("  mm:  ", sprintf("%.6f", sort(re_mm_plate)[1:6]), "\n")

if (!all(is.na(re_mm_plate)) && length(re_mm_plate) == length(re_lme4_plate)) {
  compare("ranef plate max abs diff",
          sort(re_lme4_plate), sort(re_mm_plate), tols$fixef)
} else {
  cat(sprintf("ranef plate: lme4 n=%d  mm n=%d\n", length(re_lme4_plate), length(re_mm_plate)))
}

## ranef sample
re_lme4_sample <- as.numeric(lme4::ranef(fit_lme4)$sample[, 1])
re_mm_sample_df <- tryCatch(ranef(fit_mm)$sample, error = function(e) NULL)
re_mm_sample <- if (!is.null(re_mm_sample_df)) as.numeric(re_mm_sample_df[, 1]) else NA_real_

cat(sprintf("\n%-35s\n", "ranef sample (sorted):"))
cat("  lme4:", sprintf("%.6f", sort(re_lme4_sample)), "\n")
if (!all(is.na(re_mm_sample))) cat("  mm:  ", sprintf("%.6f", sort(re_mm_sample)), "\n")

if (!all(is.na(re_mm_sample)) && length(re_mm_sample) == length(re_lme4_sample)) {
  compare("ranef sample max abs diff",
          sort(re_lme4_sample), sort(re_mm_sample), tols$fixef)
} else {
  cat(sprintf("ranef sample: lme4 n=%d  mm n=%d\n", length(re_lme4_sample), length(re_mm_sample)))
}

## speed ratio
t_lme4_mean   <- t_lme4["elapsed"] / N_REPS
t_mixeff_mean <- t_mixeff["elapsed"] / N_REPS
cat(sprintf("\n%-35s lme4=%.4fs  mm=%.4fs  ratio(mm/lme4)=%.2fx\n",
            "wall-clock (mean per rep)",
            t_lme4_mean, t_mixeff_mean,
            t_mixeff_mean / t_lme4_mean))

cat("\n=== SUMMARY TABLE ===\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
cat("All WITHIN-TOL = parity achieved.\n")
