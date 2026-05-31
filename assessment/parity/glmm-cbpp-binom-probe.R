## Empirical parity probe: glmm-cbpp-binom
## Dataset: cbpp  Formula: cbind(incidence, size - incidence) ~ period + (1|herd)
## Family: binomial, link: logit, method: Laplace (nAGQ=1)
## Compares lme4 vs mixeff on all relevant GLMM quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

## ── 1. Data ───────────────────────────────────────────────────────────────────
data(cbpp, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(cbpp), "  ncol:", ncol(cbpp), "\n")
print(head(cbpp))
cat("\n")

## ── 2. Fit lme4 ──────────────────────────────────────────────────────────────
cat("=== lme4 FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- glmer(
    cbind(incidence, size - incidence) ~ period + (1 | herd),
    data = cbpp,
    family = binomial,
    nAGQ = 1L   # Laplace
  )
})
cat("lme4 wall-clock (seconds):", t_lme4["elapsed"], "\n\n")

cat("-- fixef --\n"); print(lme4::fixef(fit_lme4))
cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov --\n"); print(as.matrix(vcov(fit_lme4)))
cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))
cat("-- theta (getME) --\n"); print(lme4::getME(fit_lme4, "theta"))
cat("-- sigma --\n"); print(sigma(fit_lme4))
cat("-- logLik --\n"); print(logLik(fit_lme4))
cat("-- AIC --\n");    print(AIC(fit_lme4))
cat("-- BIC --\n");    print(BIC(fit_lme4))
cat("-- deviance --\n"); print(deviance(fit_lme4))
cat("-- ranef --\n"); print(lme4::ranef(fit_lme4))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
cat("-- convergence --\n")
msgs <- fit_lme4@optinfo$conv$lme4$messages
cat("No convergence warning =", length(msgs) == 0, "\n")
if (length(msgs)) cat("Messages:", paste(msgs, collapse="; "), "\n")
cat("\n")

## ── 3. Fit mixeff ─────────────────────────────────────────────────────────────
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    glmm(
      cbind(incidence, size - incidence) ~ period + (1 | herd),
      data   = cbpp,
      family = binomial(link = "logit"),
      method = "pirls_profiled",
      nAGQ   = 1L,
      control = mm_control(verbose = -1L)
    ),
    error = function(e) e
  )
})
cat("mixeff wall-clock (seconds):", t_mixeff["elapsed"], "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR/REFUSAL !!!\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm), "\n")
  if (!is.null(fit_mm$reason_code)) cat("reason_code:", fit_mm$reason_code, "\n")
  cat("\n=== NUMERICAL COMPARISON ===\n")
  cat("mixeff failed — no numerical comparison possible.\n")
  quit(status = 0)
}

cat("fit_status:", fit_mm$fit_status, "\n\n")

cat("-- fixef --\n"); print(fixef(fit_mm))
cat("-- SE --\n"); print(fit_mm$std_errors)
cat("-- vcov --\n"); print(as.matrix(fit_mm$fixed_effect_vcov))
cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
cat("-- theta --\n"); print(fit_mm$theta)
cat("-- sigma --\n"); print(sigma(fit_mm))
cat("-- logLik --\n"); print(logLik(fit_mm))
cat("-- AIC --\n");    print(fit_mm$AIC)
cat("-- BIC --\n");    print(fit_mm$BIC)
cat("-- deviance --\n"); print(fit_mm$deviance)
cat("-- ranef --\n"); print(ranef(fit_mm))
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("\n")

## ── 4. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  diff     <- max(abs(lme4_val - mm_val))
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL*"
  cat(sprintf("%-35s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse=", "),
              paste(sprintf(fmt, mm_val),   collapse=", "),
              diff, tol, status))
  invisible(diff)
}

## fixef (4 coefficients: Intercept, period2, period3, period4)
fe_lme4 <- lme4::fixef(fit_lme4)
fe_mm   <- fixef(fit_mm)

# align by name
nm <- names(fe_lme4)
fe_mm_aligned <- fe_mm[nm]
for (n in nm) {
  compare(paste0("fixef ", n),
          fe_lme4[[n]],
          fe_mm_aligned[[n]],
          tols$fixef)
}
compare("fixef (all)",
        fe_lme4[nm],
        fe_mm_aligned[nm],
        tols$fixef)

## SE
se_lme4 <- sqrt(diag(vcov(fit_lme4)))
se_mm   <- fit_mm$std_errors[nm]
for (n in nm) {
  compare(paste0("SE ", n),
          se_lme4[[n]],
          se_mm[[n]],
          tols$fixef)
}
compare("SE (all)", se_lme4[nm], se_mm[nm], tols$fixef)

## theta
theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
cat(sprintf("%-35s lme4=%s  mm=%s\n", "theta",
            paste(sprintf("%.8f", theta_lme4), collapse=", "),
            paste(sprintf("%.8f", theta_mm),   collapse=", ")))
compare("theta", theta_lme4, theta_mm, tols$theta)

## sigma (dispersion)
compare("sigma (dispersion)",
        sigma(fit_lme4),
        sigma(fit_mm),
        tols$sigma)

## VarCorr: herd random-intercept variance
vc_lme4 <- as.numeric(lme4::VarCorr(fit_lme4)$herd)
vc_mm_obj <- VarCorr(fit_mm)
vc_mm <- tryCatch({
  if (is.list(vc_mm_obj) && !is.null(vc_mm_obj$table)) {
    tbl <- vc_mm_obj$table
    as.numeric(tbl[tbl$group == "herd", "variance"])
  } else if (is.data.frame(vc_mm_obj)) {
    as.numeric(vc_mm_obj[vc_mm_obj$group == "herd", "variance"])
  } else {
    NA_real_
  }
}, error = function(e) { cat("VarCorr extraction error:", conditionMessage(e), "\n"); NA_real_ })

cat(sprintf("%-35s lme4=%.8f  mm=%s\n", "VarCorr herd var",
            vc_lme4,
            ifelse(is.na(vc_mm), "NA", sprintf("%.8f", vc_mm))))
if (!is.na(vc_mm)) compare("VarCorr herd var", vc_lme4, vc_mm, tols$theta)

## logLik
compare("logLik",
        as.numeric(logLik(fit_lme4)),
        as.numeric(logLik(fit_mm)),
        tols$logLik)

## deviance (= -2 * logLik for GLMM)
compare("deviance",
        deviance(fit_lme4),
        fit_mm$deviance,
        tols$logLik * 2)

## AIC
compare("AIC",
        AIC(fit_lme4),
        fit_mm$AIC,
        tols$logLik * 2)

## BIC
compare("BIC",
        BIC(fit_lme4),
        fit_mm$BIC,
        tols$logLik * 2)

## fitted values (conditional mean on response scale)
fitted_lme4 <- fitted(fit_lme4)
fitted_mm   <- fit_mm$fitted
n_common <- min(length(fitted_lme4), length(fitted_mm))
compare("fitted max abs diff",
        fitted_lme4[seq_len(n_common)],
        fitted_mm[seq_len(n_common)],
        tols$fixef)

## ranef
re_lme4_df <- lme4::ranef(fit_lme4)$herd
re_mm_df   <- ranef(fit_mm)$herd
if (!is.null(re_mm_df) && !is.null(re_lme4_df)) {
  # align by herd level name
  lvl_common <- intersect(rownames(re_lme4_df), rownames(re_mm_df))
  re_lme4_v  <- as.numeric(re_lme4_df[lvl_common, 1])
  re_mm_v    <- as.numeric(re_mm_df[lvl_common, 1])
  cat(sprintf("%-35s\n", "ranef herd (Intercept), sorted:"))
  cat("  lme4:", sprintf("%.5f", sort(re_lme4_v)), "\n")
  cat("  mm:  ", sprintf("%.5f", sort(re_mm_v)),   "\n")
  compare("ranef herd max abs diff",
          sort(re_lme4_v), sort(re_mm_v), tols$fixef)
} else {
  cat("ranef: one or both fits returned NULL — cannot compare\n")
}

## speed ratio
cat(sprintf("\n%-35s lme4=%.4fs  mm=%.4fs  ratio(mm/lme4)=%.2fx\n",
            "wall-clock elapsed",
            t_lme4["elapsed"], t_mixeff["elapsed"],
            t_mixeff["elapsed"] / t_lme4["elapsed"]))

cat("\n=== SUMMARY TABLE ===\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
cat("WITHIN-TOL = parity achieved.  EXCEEDS-TOL* = divergence.\n")
