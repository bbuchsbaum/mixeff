## Empirical parity probe: glmm-gamma
## Dataset: simulated Gamma  Formula: y ~ x + (1|g)  Family: Gamma(link="log")
## Compares lme4 vs mixeff on all relevant quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

## ── 1. Simulate Gamma data ───────────────────────────────────────────────────
set.seed(42)
n_groups <- 20
n_per    <- 15
N        <- n_groups * n_per

g  <- factor(rep(seq_len(n_groups), each = n_per))
x  <- rnorm(N)
# true params: intercept=1.5, slope=0.4, random intercept sd=0.5, shape=4
u  <- rnorm(n_groups, 0, 0.5)[g]
mu <- exp(1.5 + 0.4 * x + u)
# shape parameter for Gamma: mean=mu, var=mu^2/shape => shape=4
y  <- rgamma(N, shape = 4, rate = 4 / mu)

dat <- data.frame(y = y, x = x, g = g)

cat("=== DATASET ===\n")
cat("nrow:", nrow(dat), "  ncol:", ncol(dat), "\n")
cat("y summary: min=", round(min(dat$y),4),
    " mean=", round(mean(dat$y),4),
    " max=", round(max(dat$y),4), "\n\n")

## ── 2. Fit lme4 ──────────────────────────────────────────────────────────────
cat("=== lme4 FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- tryCatch(
    glmer(y ~ x + (1 | g), data = dat, family = Gamma(link = "log"),
          control = glmerControl(optimizer = "bobyqa")),
    error = function(e) e
  )
})
cat("lme4 wall-clock (seconds):", t_lme4["elapsed"], "\n\n")

if (inherits(fit_lme4, "condition")) {
  cat("!!! lme4 ERROR !!!\n")
  cat(conditionMessage(fit_lme4), "\n")
} else {
  cat("-- fixef --\n"); print(lme4::fixef(fit_lme4))
  cat("-- SE --\n"); print(sqrt(diag(vcov(fit_lme4))))
  cat("-- vcov --\n"); print(as.matrix(vcov(fit_lme4)))
  cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))
  cat("-- sigma (dispersion) --\n"); print(sigma(fit_lme4))
  cat("-- logLik --\n"); print(logLik(fit_lme4))
  cat("-- AIC --\n");    print(AIC(fit_lme4))
  cat("-- BIC --\n");    print(BIC(fit_lme4))
  cat("-- ranef (first 6 groups) --\n")
  re <- lme4::ranef(fit_lme4)
  print(head(re$g))
  cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
  cat("-- convergence --\n")
  msgs <- fit_lme4@optinfo$conv$lme4$messages
  cat("No convergence warning =", is.null(msgs) || length(msgs) == 0, "\n\n")
}

## ── 3. Fit mixeff ─────────────────────────────────────────────────────────────
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    glmm(y ~ x + (1 | g), data = dat, family = Gamma(link = "log"),
         control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff wall-clock (seconds):", t_mixeff["elapsed"], "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR !!!\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm), "\n")
  if (!is.null(fit_mm$reason_code)) cat("reason_code:", fit_mm$reason_code, "\n")
  cat("\n=== VERDICT ===\n")
  cat("outcome: mixeff-error\n")
  cat("severity: major\n")
  quit(status = 0)
}

cat("fit_status:", fit_mm$fit_status, "\n\n")

cat("-- fixef --\n"); print(fixef(fit_mm))
cat("-- SE --\n"); print(fit_mm$std_errors)
cat("-- vcov --\n"); print(as.matrix(fit_mm$fixed_effect_vcov))
cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
cat("-- sigma (dispersion) --\n"); print(sigma(fit_mm))
cat("-- logLik --\n"); print(logLik(fit_mm))
cat("-- AIC --\n");    print(fit_mm$AIC)
cat("-- BIC --\n");    print(fit_mm$BIC)
cat("-- ranef (first 6 groups) --\n")
re_mm <- ranef(fit_mm)
if (!is.null(re_mm$g)) print(head(re_mm$g)) else cat("ranef$g is NULL\n")
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("\n")

## ── 4. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  diff     <- max(abs(lme4_val - mm_val))
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-35s lme4=%-20s  mm=%-20s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse=", "),
              paste(sprintf(fmt, mm_val),   collapse=", "),
              diff, tol, status))
  invisible(diff)
}

if (!inherits(fit_lme4, "condition")) {
  ## fixef
  compare("fixef (Intercept)",
          lme4::fixef(fit_lme4)[["(Intercept)"]],
          fixef(fit_mm)[["(Intercept)"]],
          tols$fixef)
  compare("fixef x",
          lme4::fixef(fit_lme4)[["x"]],
          fixef(fit_mm)[["x"]],
          tols$fixef)

  ## SE
  lme4_se <- sqrt(diag(vcov(fit_lme4)))
  mm_se   <- fit_mm$std_errors
  compare("SE (Intercept)",
          lme4_se[["(Intercept)"]],
          mm_se[["(Intercept)"]],
          tols$fixef)
  compare("SE x",
          lme4_se[["x"]],
          mm_se[["x"]],
          tols$fixef)

  ## vcov diagonal
  compare("vcov[1,1]",
          as.numeric(vcov(fit_lme4)[1,1]),
          as.numeric(fit_mm$fixed_effect_vcov[1,1]),
          tols$fixef * 10)
  compare("vcov[2,2]",
          as.numeric(vcov(fit_lme4)[2,2]),
          as.numeric(fit_mm$fixed_effect_vcov[2,2]),
          tols$fixef * 10)

  ## theta (Cholesky factor of RE covariance)
  theta_lme4 <- lme4::getME(fit_lme4, "theta")
  theta_mm   <- fit_mm$theta
  cat(sprintf("%-35s lme4=%s  mm=%s\n", "theta (raw)",
              paste(sprintf("%.8f", theta_lme4), collapse=", "),
              paste(sprintf("%.8f", theta_mm),   collapse=", ")))
  compare("theta", theta_lme4, theta_mm, tols$theta)

  ## sigma / dispersion
  compare("sigma (dispersion)",
          sigma(fit_lme4),
          sigma(fit_mm),
          tols$sigma)

  ## VarCorr: RE variance for group g
  vc_lme4 <- as.numeric(attr(lme4::VarCorr(fit_lme4)$g, "stddev"))^2
  vc_mm_raw <- VarCorr(fit_mm)
  vc_mm <- tryCatch({
    tbl <- vc_mm_raw$table
    if (!is.null(tbl) && "variance" %in% names(tbl)) {
      as.numeric(tbl$variance[1])
    } else if (is.data.frame(vc_mm_raw)) {
      as.numeric(vc_mm_raw[vc_mm_raw$grp == "g", "vcov"])
    } else {
      NA_real_
    }
  }, error = function(e) { cat("VarCorr extraction error:", conditionMessage(e), "\n"); NA_real_ })

  cat(sprintf("%-35s lme4=%.8f  mm=%s\n", "VarCorr g variance",
              vc_lme4,
              ifelse(is.na(vc_mm), "NA", sprintf("%.8f", vc_mm))))
  if (!is.na(vc_mm)) compare("VarCorr g variance", vc_lme4, vc_mm, tols$theta)

  ## logLik
  compare("logLik",
          as.numeric(logLik(fit_lme4)),
          as.numeric(logLik(fit_mm)),
          tols$logLik)

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

  ## fitted values
  fitted_lme4 <- fitted(fit_lme4)
  fitted_mm   <- fit_mm$fitted
  compare("fitted max abs diff",
          fitted_lme4,
          fitted_mm[seq_along(fitted_lme4)],
          tols$fixef)

  ## ranef
  re_lme4 <- as.numeric(lme4::ranef(fit_lme4)$g[,1])
  re_mm_df <- ranef(fit_mm)$g
  re_mm    <- if (!is.null(re_mm_df)) as.numeric(re_mm_df[,1]) else NA_real_
  cat(sprintf("\n%-35s\n", "ranef g (sorted first 6):"))
  cat("  lme4:", sprintf("%.6f", sort(re_lme4)[1:6]), "\n")
  if (!any(is.na(re_mm))) cat("  mm:  ", sprintf("%.6f", sort(re_mm)[1:6]),   "\n")
  compare("ranef g max abs diff",
          sort(re_lme4), sort(re_mm), tols$fixef)
}

## speed ratio
cat(sprintf("\n%-35s lme4=%.4fs  mm=%.4fs  ratio(mm/lme4)=%.2fx\n",
            "wall-clock elapsed",
            t_lme4["elapsed"], t_mixeff["elapsed"],
            t_mixeff["elapsed"] / t_lme4["elapsed"]))

cat("\n=== SUMMARY TABLE ===\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
