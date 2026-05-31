## Empirical parity probe: glmm-cbpp-agq
## Dataset: cbpp  Formula: cbind(incidence, size - incidence) ~ period + (1|herd)
## Family: binomial / logit
## Focus: nAGQ=1 (Laplace) AND nAGQ=2+ (AGQ) — note nlopt feature gate
## Compares lme4::glmer vs mixeff::glmm on all relevant quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

## ── 1. Data ──────────────────────────────────────────────────────────────────
data(cbpp, package = "lme4")
cat("=== DATASET (cbpp) ===\n")
cat("nrow:", nrow(cbpp), "  ncol:", ncol(cbpp), "\n")
print(head(cbpp))
cat("\n")

fmla <- cbind(incidence, size - incidence) ~ period + (1 | herd)

## ── 2. Fit lme4 nAGQ=1 (Laplace) ─────────────────────────────────────────────
cat("=== lme4 FIT (nAGQ=1, Laplace) ===\n")
t_lme4_1 <- system.time({
  fit_lme4_1 <- suppressMessages(suppressWarnings(
    lme4::glmer(fmla, data = cbpp, family = binomial(link = "logit"), nAGQ = 1L)
  ))
})
cat("lme4 nAGQ=1 wall-clock (seconds):", t_lme4_1["elapsed"], "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4_1))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4_1))))
cat("-- vcov --\n"); print(as.matrix(vcov(fit_lme4_1)))
cat("-- VarCorr (theta) --\n"); print(lme4::VarCorr(fit_lme4_1))
cat("-- theta --\n"); print(lme4::getME(fit_lme4_1, "theta"))
cat("-- sigma --\n"); print(sigma(fit_lme4_1))
cat("-- logLik --\n"); print(logLik(fit_lme4_1))
cat("-- AIC --\n"); print(AIC(fit_lme4_1))
cat("-- BIC --\n"); print(BIC(fit_lme4_1))
cat("-- ranef (first herd) --\n"); print(head(lme4::ranef(fit_lme4_1)$herd, 5))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4_1)))
cat("-- convergence --\n")
msgs <- fit_lme4_1@optinfo$conv$lme4$messages
cat("No convergence warning =", length(msgs) == 0, "\n\n")

## ── 3. Fit lme4 nAGQ=2 (2-point AGQ) ────────────────────────────────────────
cat("=== lme4 FIT (nAGQ=2, AGQ) ===\n")
t_lme4_2 <- system.time({
  fit_lme4_2 <- suppressMessages(suppressWarnings(
    lme4::glmer(fmla, data = cbpp, family = binomial(link = "logit"), nAGQ = 2L)
  ))
})
cat("lme4 nAGQ=2 wall-clock (seconds):", t_lme4_2["elapsed"], "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4_2))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4_2))))
cat("-- theta --\n"); print(lme4::getME(fit_lme4_2, "theta"))
cat("-- logLik --\n"); print(logLik(fit_lme4_2))
cat("-- AIC --\n"); print(AIC(fit_lme4_2))
cat("-- BIC --\n"); print(BIC(fit_lme4_2))

## ── 4. Fit lme4 nAGQ=10 (10-point AGQ) ───────────────────────────────────────
cat("=== lme4 FIT (nAGQ=10, AGQ) ===\n")
t_lme4_10 <- system.time({
  fit_lme4_10 <- suppressMessages(suppressWarnings(
    lme4::glmer(fmla, data = cbpp, family = binomial(link = "logit"), nAGQ = 10L)
  ))
})
cat("lme4 nAGQ=10 wall-clock (seconds):", t_lme4_10["elapsed"], "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4_10))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4_10))))
cat("-- theta --\n"); print(lme4::getME(fit_lme4_10, "theta"))
cat("-- logLik --\n"); print(logLik(fit_lme4_10))
cat("-- AIC --\n"); print(AIC(fit_lme4_10))
cat("-- BIC --\n"); print(BIC(fit_lme4_10))

## ── 5. Fit mixeff nAGQ=1 (pirls_profiled, Laplace) ───────────────────────────
cat("=== mixeff FIT (nAGQ=1, pirls_profiled) ===\n")
t_mm_1 <- system.time({
  fit_mm_1 <- tryCatch(
    glmm(fmla, data = cbpp, family = binomial(link = "logit"),
         method = "pirls_profiled", nAGQ = 1L,
         control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff nAGQ=1 wall-clock (seconds):", t_mm_1["elapsed"], "\n\n")

if (inherits(fit_mm_1, "condition")) {
  cat("!!! mixeff nAGQ=1 ERROR !!!\n")
  cat("class:", paste(class(fit_mm_1), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm_1), "\n\n")
  fit_mm_1_ok <- FALSE
} else {
  fit_mm_1_ok <- TRUE
  cat("fit_status:", fit_mm_1$fit_status, "\n")
  cat("nAGQ stored:", fit_mm_1$nAGQ, "\n\n")
  cat("-- fixef --\n"); print(fixef(fit_mm_1))
  cat("-- SE --\n"); print(fit_mm_1$std_errors)
  cat("-- vcov --\n"); print(as.matrix(fit_mm_1$fixed_effect_vcov))
  cat("-- VarCorr --\n"); print(VarCorr(fit_mm_1))
  cat("-- theta --\n"); print(fit_mm_1$theta)
  cat("-- sigma --\n"); print(sigma(fit_mm_1))
  cat("-- logLik --\n"); print(logLik(fit_mm_1))
  cat("-- AIC --\n"); print(fit_mm_1$AIC)
  cat("-- BIC --\n"); print(fit_mm_1$BIC)
  cat("-- ranef (first 5 herds) --\n"); print(head(ranef(fit_mm_1)$herd, 5))
  cat("-- fitted (first 6) --\n"); print(head(fit_mm_1$fitted))
  cat("\n")
}

## ── 6. Fit mixeff nAGQ=2 (pirls_profiled, AGQ=2) ─────────────────────────────
cat("=== mixeff FIT (nAGQ=2, pirls_profiled) ===\n")
t_mm_2 <- system.time({
  fit_mm_2 <- tryCatch(
    glmm(fmla, data = cbpp, family = binomial(link = "logit"),
         method = "pirls_profiled", nAGQ = 2L,
         control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff nAGQ=2 wall-clock (seconds):", t_mm_2["elapsed"], "\n\n")

if (inherits(fit_mm_2, "condition")) {
  cat("!!! mixeff nAGQ=2 ERROR/REFUSAL !!!\n")
  cat("class:", paste(class(fit_mm_2), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm_2), "\n\n")
  fit_mm_2_ok <- FALSE
} else {
  fit_mm_2_ok <- TRUE
  cat("fit_status:", fit_mm_2$fit_status, "\n")
  cat("nAGQ stored:", fit_mm_2$nAGQ, "\n\n")
  cat("-- fixef --\n"); print(fixef(fit_mm_2))
  cat("-- logLik --\n"); print(logLik(fit_mm_2))
  cat("-- AIC --\n"); print(fit_mm_2$AIC)
  cat("\n")
}

## ── 7. Fit mixeff joint_laplace (expected refusal — no nlopt) ────────────────
cat("=== mixeff FIT (joint_laplace — expected refusal) ===\n")
t_mm_jl <- system.time({
  fit_mm_jl <- tryCatch(
    glmm(fmla, data = cbpp, family = binomial(link = "logit"),
         method = "joint_laplace", nAGQ = 1L,
         control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff joint_laplace wall-clock (seconds):", t_mm_jl["elapsed"], "\n\n")

if (inherits(fit_mm_jl, "condition")) {
  cat("mixeff joint_laplace REFUSED (expected):\n")
  cat("class:", paste(class(fit_mm_jl), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm_jl), "\n\n")
} else {
  cat("mixeff joint_laplace succeeded (unexpected).\n\n")
}

## ── 8. Numerical comparison: nAGQ=1 ──────────────────────────────────────────
cat("=== NUMERICAL COMPARISON (lme4 nAGQ=1 vs mixeff nAGQ=1) ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  diff     <- max(abs(lme4_val - mm_val))
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-40s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse=", "),
              paste(sprintf(fmt, mm_val),   collapse=", "),
              diff, tol, status))
  invisible(diff)
}

if (fit_mm_1_ok) {
  fe_lme4 <- lme4::fixef(fit_lme4_1)
  fe_mm   <- fixef(fit_mm_1)
  # match names
  common_names <- intersect(names(fe_lme4), names(fe_mm))
  for (nm in common_names) {
    compare(paste0("fixef ", nm), fe_lme4[[nm]], fe_mm[[nm]], tols$fixef)
  }

  se_lme4 <- sqrt(diag(vcov(fit_lme4_1)))
  se_mm   <- fit_mm_1$std_errors
  for (nm in intersect(names(se_lme4), names(se_mm))) {
    compare(paste0("SE ", nm), se_lme4[[nm]], se_mm[[nm]], tols$fixef)
  }

  compare("theta",
          lme4::getME(fit_lme4_1, "theta"),
          fit_mm_1$theta,
          tols$theta)

  compare("logLik",
          as.numeric(logLik(fit_lme4_1)),
          as.numeric(logLik(fit_mm_1)),
          tols$logLik)

  compare("AIC",
          AIC(fit_lme4_1),
          fit_mm_1$AIC,
          tols$logLik * 2)

  compare("BIC",
          BIC(fit_lme4_1),
          fit_mm_1$BIC,
          tols$logLik * 2)

  # VarCorr / herd variance
  vc_lme4 <- as.numeric(lme4::VarCorr(fit_lme4_1)$herd[1,1])
  vc_mm <- tryCatch({
    vc_obj <- VarCorr(fit_mm_1)
    if (is.list(vc_obj) && !is.null(vc_obj$table)) {
      as.numeric(vc_obj$table$variance[1])
    } else if (is.data.frame(vc_obj)) {
      as.numeric(vc_obj[vc_obj$grp == "herd", "vcov"])
    } else {
      NA_real_
    }
  }, error = function(e) { cat("VarCorr extract error:", conditionMessage(e), "\n"); NA_real_ })

  cat(sprintf("%-40s lme4=%.8f  mm=%s\n", "VarCorr herd variance",
              vc_lme4, ifelse(is.na(vc_mm), "NA", sprintf("%.8f", vc_mm))))
  if (!is.na(vc_mm)) compare("VarCorr herd variance", vc_lme4, vc_mm, tols$theta)

  # fitted values
  fitted_lme4 <- fitted(fit_lme4_1)
  fitted_mm   <- fit_mm_1$fitted
  if (length(fitted_mm) == length(fitted_lme4)) {
    compare("fitted max abs diff", fitted_lme4, fitted_mm, tols$fixef)
  } else {
    cat(sprintf("%-40s lme4 n=%d  mm n=%d  LENGTHS DIFFER\n",
                "fitted lengths", length(fitted_lme4), length(fitted_mm)))
  }

  # ranef
  re_lme4 <- sort(as.numeric(lme4::ranef(fit_lme4_1)$herd[,1]))
  re_mm_df <- ranef(fit_mm_1)$herd
  if (!is.null(re_mm_df)) {
    re_mm <- sort(as.numeric(re_mm_df[,1]))
    if (length(re_mm) == length(re_lme4)) {
      compare("ranef herd max abs diff (sorted)", re_lme4, re_mm, tols$fixef)
    } else {
      cat(sprintf("%-40s lme4 n=%d  mm n=%d  LENGTHS DIFFER\n",
                  "ranef lengths", length(re_lme4), length(re_mm)))
    }
  } else {
    cat("ranef: mm returned NULL\n")
  }
}

## ── 9. Comparison: lme4 nAGQ=2 vs lme4 nAGQ=1 ───────────────────────────────
cat("\n=== lme4 AGQ SENSITIVITY: nAGQ=1 vs nAGQ=2 vs nAGQ=10 ===\n")
cat("(Shows how much AGQ matters for this dataset)\n\n")

lme4_fixef_compare <- function(fit_a, fit_b, label_a, label_b) {
  fe_a <- lme4::fixef(fit_a); fe_b <- lme4::fixef(fit_b)
  for (nm in intersect(names(fe_a), names(fe_b))) {
    diff <- abs(fe_a[[nm]] - fe_b[[nm]])
    cat(sprintf("fixef %-20s %s=%.8f  %s=%.8f  diff=%.3e\n",
                nm, label_a, fe_a[[nm]], label_b, fe_b[[nm]], diff))
  }
  ll_a <- as.numeric(logLik(fit_a)); ll_b <- as.numeric(logLik(fit_b))
  cat(sprintf("logLik                        %s=%.8f  %s=%.8f  diff=%.3e\n",
              label_a, ll_a, label_b, ll_b, abs(ll_a - ll_b)))
}

cat("--- lme4 nAGQ=1 vs nAGQ=2 ---\n")
lme4_fixef_compare(fit_lme4_1, fit_lme4_2, "nAGQ=1", "nAGQ=2")
cat("\n--- lme4 nAGQ=1 vs nAGQ=10 ---\n")
lme4_fixef_compare(fit_lme4_1, fit_lme4_10, "nAGQ=1", "nAGQ=10")
cat("\n--- lme4 nAGQ=2 vs nAGQ=10 ---\n")
lme4_fixef_compare(fit_lme4_2, fit_lme4_10, "nAGQ=2", "nAGQ=10")

## ── 10. Comparison: mixeff nAGQ=2 vs lme4 nAGQ=2 ─────────────────────────────
if (fit_mm_2_ok) {
  cat("\n=== NUMERICAL COMPARISON (lme4 nAGQ=2 vs mixeff nAGQ=2) ===\n")
  fe_lme4 <- lme4::fixef(fit_lme4_2)
  fe_mm   <- fixef(fit_mm_2)
  for (nm in intersect(names(fe_lme4), names(fe_mm))) {
    compare(paste0("fixef ", nm), fe_lme4[[nm]], fe_mm[[nm]], tols$fixef)
  }
  compare("theta",
          lme4::getME(fit_lme4_2, "theta"),
          fit_mm_2$theta,
          tols$theta)
  compare("logLik",
          as.numeric(logLik(fit_lme4_2)),
          as.numeric(logLik(fit_mm_2)),
          tols$logLik)
}

## ── 11. Speed summary ────────────────────────────────────────────────────────
cat("\n=== SPEED SUMMARY ===\n")
cat(sprintf("lme4 nAGQ=1:    %.4fs\n", t_lme4_1["elapsed"]))
cat(sprintf("lme4 nAGQ=2:    %.4fs\n", t_lme4_2["elapsed"]))
cat(sprintf("lme4 nAGQ=10:   %.4fs\n", t_lme4_10["elapsed"]))
if (fit_mm_1_ok) {
  cat(sprintf("mixeff nAGQ=1:  %.4fs  ratio(mm/lme4)=%.2fx\n",
              t_mm_1["elapsed"], t_mm_1["elapsed"] / t_lme4_1["elapsed"]))
}
if (fit_mm_2_ok) {
  cat(sprintf("mixeff nAGQ=2:  %.4fs  ratio(mm/lme4-nAGQ2)=%.2fx\n",
              t_mm_2["elapsed"], t_mm_2["elapsed"] / t_lme4_2["elapsed"]))
}

cat("\n=== DONE ===\n")
