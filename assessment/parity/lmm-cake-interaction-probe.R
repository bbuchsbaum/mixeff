## Empirical parity probe: lmm-cake-interaction
## Dataset: cake  Formula: angle ~ recipe * temperature + (1|recipe:replicate)
## Focus: fixed interaction + interaction grouping
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

## ── 1. Data ──────────────────────────────────────────────────────────────────
data(cake, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(cake), "  ncol:", ncol(cake), "\n")
cat("recipe levels:", levels(cake$recipe), "\n")
cat("temperature levels:", levels(cake$temperature), "\n")
cat("recipe:replicate groups:", nlevels(interaction(cake$recipe, cake$replicate)), "\n")
print(head(cake))
cat("\n")

## ── 2. Fit lme4 ──────────────────────────────────────────────────────────────
cat("=== lme4 FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
                   data = cake, REML = TRUE)
})
cat("lme4 wall-clock (seconds):", t_lme4["elapsed"], "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov (diagonal only) --\n"); print(diag(as.matrix(vcov(fit_lme4))))
cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))
cat("-- sigma --\n"); print(sigma(fit_lme4))
cat("-- logLik --\n"); print(logLik(fit_lme4))
cat("-- AIC --\n"); print(AIC(fit_lme4))
cat("-- BIC --\n"); print(BIC(fit_lme4))
cat("-- ranef (first 6) --\n"); print(head(lme4::ranef(fit_lme4)[["recipe:replicate"]]))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
cat("-- convergence --\n")
cat("No convergence warning =",
    length(fit_lme4@optinfo$conv$lme4$messages) == 0, "\n\n")

## ── 3. Fit mixeff ─────────────────────────────────────────────────────────────
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
        data = cake, REML = TRUE,
        control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff wall-clock (seconds):", t_mixeff["elapsed"], "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR !!!\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm), "\n")
  ## Print full condition for diagnosis
  cat("--- full condition ---\n")
  print(fit_mm)
  cat("\n=== SUMMARY ===\n")
  cat("mixeff: ERROR (see above)\n")
  cat("lme4:   succeeded\n")
  cat("Outcome: mixeff-error\n")
  quit(status = 0)
}

cat("fit_status:", fit_mm$fit_status, "\n\n")

cat("-- fixef --\n");         print(fixef(fit_mm))
cat("-- SE --\n");            print(fit_mm$std_errors)
cat("-- vcov (diagonal) --\n"); print(diag(as.matrix(fit_mm$fixed_effect_vcov)))
cat("-- VarCorr --\n");       print(VarCorr(fit_mm))
cat("-- sigma --\n");         print(sigma(fit_mm))
cat("-- logLik --\n");        print(logLik(fit_mm))
cat("-- AIC --\n");           print(fit_mm$AIC)
cat("-- BIC --\n");           print(fit_mm$BIC)
cat("-- ranef (first 6) --\n")
re_mm_raw <- tryCatch(ranef(fit_mm), error = function(e) {
  cat("ranef error:", conditionMessage(e), "\n"); NULL
})
if (!is.null(re_mm_raw)) print(head(re_mm_raw[[1]]))
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("\n")

## ── 4. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  if (any(is.na(mm_val))) {
    cat(sprintf("%-44s lme4=%s  mm=NA  [MISSING]\n",
                label,
                paste(sprintf(fmt, lme4_val), collapse = ", ")))
    return(invisible(NA_real_))
  }
  diff   <- max(abs(lme4_val - mm_val))
  status <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-44s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse = ", "),
              paste(sprintf(fmt, mm_val),   collapse = ", "),
              diff, tol, status))
  invisible(diff)
}

## ── fixef ──
fe_lme4 <- lme4::fixef(fit_lme4)
fe_mm   <- fixef(fit_mm)

## Align by name
fe_names <- names(fe_lme4)
cat("fixef names (lme4):", paste(fe_names, collapse = ", "), "\n")
cat("fixef names (mm):  ", paste(names(fe_mm), collapse = ", "), "\n\n")

for (nm in fe_names) {
  mm_val <- if (nm %in% names(fe_mm)) fe_mm[[nm]] else NA_real_
  compare(paste0("fixef[", nm, "]"),
          fe_lme4[[nm]], mm_val, tols$fixef)
}

## ── SE ──
cat("\n")
se_lme4 <- sqrt(diag(vcov(fit_lme4)))
se_mm   <- fit_mm$std_errors
for (nm in names(se_lme4)) {
  mm_val <- if (nm %in% names(se_mm)) se_mm[[nm]] else NA_real_
  compare(paste0("SE[", nm, "]"),
          se_lme4[[nm]], mm_val, tols$fixef)
}

## ── theta ──
cat("\n")
theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
cat(sprintf("%-44s lme4=%s  mm=%s\n", "theta (raw)",
            paste(sprintf("%.8f", theta_lme4), collapse = ", "),
            paste(sprintf("%.8f", as.numeric(theta_mm)), collapse = ", ")))
compare("theta", theta_lme4, theta_mm, tols$theta)

## ── sigma ──
cat("\n")
compare("sigma", sigma(fit_lme4), sigma(fit_mm), tols$sigma)

## ── VarCorr: RE variance ──
vc_lme4 <- as.numeric(attr(lme4::VarCorr(fit_lme4)[["recipe:replicate"]], "stddev"))^2
vc_mm_df <- VarCorr(fit_mm)
vc_mm <- tryCatch({
  grp_name <- "recipe:replicate"
  if (is.data.frame(vc_mm_df)) {
    rows <- vc_mm_df[vc_mm_df$grp == grp_name, , drop = FALSE]
    if (nrow(rows) > 0) as.numeric(rows$vcov[1]) else NA_real_
  } else if (is.list(vc_mm_df)) {
    v <- vc_mm_df[[grp_name]]
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
cat(sprintf("%-44s lme4=%.8f  mm=%s\n", "VarCorr recipe:replicate var",
            vc_lme4,
            ifelse(is.na(vc_mm), "NA", sprintf("%.8f", vc_mm))))
if (!is.na(vc_mm)) compare("VarCorr recipe:replicate var", vc_lme4, vc_mm, tols$theta)

## ── logLik ──
cat("\n")
compare("logLik",
        as.numeric(logLik(fit_lme4)),
        as.numeric(logLik(fit_mm)),
        tols$logLik)

## ── AIC/BIC ──
compare("AIC",
        AIC(fit_lme4),
        fit_mm$AIC,
        tols$logLik * 2)

compare("BIC",
        BIC(fit_lme4),
        fit_mm$BIC,
        tols$logLik * 2)

## ── fitted ──
fitted_lme4 <- as.numeric(fitted(fit_lme4))
fitted_mm   <- as.numeric(fit_mm$fitted)
compare("fitted max abs diff",
        fitted_lme4,
        fitted_mm[seq_along(fitted_lme4)],
        tols$fixef)

## ── ranef ──
cat("\n")
re_grp <- "recipe:replicate"
re_lme4 <- as.numeric(lme4::ranef(fit_lme4)[[re_grp]][, 1])
re_mm_df <- tryCatch(ranef(fit_mm), error = function(e) NULL)
re_mm <- if (!is.null(re_mm_df)) {
  nm <- re_grp
  if (nm %in% names(re_mm_df)) {
    as.numeric(re_mm_df[[nm]][, 1])
  } else if (length(re_mm_df) > 0) {
    as.numeric(re_mm_df[[1]][, 1])
  } else {
    NA_real_
  }
} else {
  NA_real_
}
cat(sprintf("%-44s\n", "ranef recipe:replicate (sorted, first 10):"))
cat("  lme4:", sprintf("%.5f", sort(re_lme4)[1:10]), "\n")
if (!anyNA(re_mm)) cat("  mm:  ", sprintf("%.5f", sort(re_mm)[1:10]), "\n")
compare("ranef recipe:replicate max abs diff",
        sort(re_lme4), if (anyNA(re_mm)) NA_real_ else sort(re_mm), tols$fixef)

## ── 5. Speed comparison ───────────────────────────────────────────────────────
cat("\n=== SPEED COMPARISON ===\n")
NREPS <- 5
t_lme4_rep <- system.time(
  for (i in seq_len(NREPS))
    lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
         data = cake, REML = TRUE)
)
t_mm_rep <- system.time(
  for (i in seq_len(NREPS))
    lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
        data = cake, REML = TRUE,
        control = mm_control(verbose = -1L))
)
lme4_per <- t_lme4_rep["elapsed"] / NREPS
mm_per   <- t_mm_rep["elapsed"]   / NREPS
cat(sprintf("lme4  mean/fit: %.4f s  (over %d reps)\n", lme4_per, NREPS))
cat(sprintf("mm    mean/fit: %.4f s  (over %d reps)\n", mm_per,   NREPS))
cat(sprintf("ratio (mm/lme4): %.2fx\n", mm_per / lme4_per))

cat("\n=== SUMMARY TABLE ===\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
cat("All WITHIN-TOL = parity achieved.\n")
