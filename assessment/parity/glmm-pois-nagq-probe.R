## glmm-pois-nagq-probe.R
## Empirical parity probe: GLMM Poisson, y ~ x + (1|g), nAGQ>1
## Uses explicit :: everywhere to avoid namespace masking.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

cat("=== Environment ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n")
cat("R version:", R.version$version.string, "\n\n")

## ---- 1. Simulate dataset ---------------------------------------------------
set.seed(42)
n_groups    <- 30
n_per_group <- 10
n_obs       <- n_groups * n_per_group

g          <- factor(rep(seq_len(n_groups), each = n_per_group))
x          <- rnorm(n_obs)
re_true_sd <- 0.4
b          <- rnorm(n_groups, 0, re_true_sd)
eta        <- 1 + 0.5 * x + b[as.integer(g)]
y          <- rpois(n_obs, exp(eta))
dat        <- data.frame(y = y, x = x, g = g)

cat("=== Dataset summary ===\n")
cat("N obs:", n_obs, " Groups:", n_groups, "\n")
cat("y range:", range(y), "  mean:", round(mean(y), 3), "\n\n")

## ---- Helper: tolerances ----------------------------------------------------
TOL_FIXEF  <- 1e-4
TOL_THETA  <- 1e-3
TOL_LOGLIK <- 1e-3
TOL_SIGMA  <- 1e-4

check_tol <- function(diff, tol, label) {
  within <- abs(diff) <= tol
  status <- if (within) "WITHIN-TOL" else "BEYOND-TOL"
  cat(sprintf("  %-36s abs_diff=%.2e  tol=%.0e  [%s]\n",
              label, abs(diff), tol, status))
  invisible(list(within = within, diff = abs(diff)))
}

## ============================================================
## SCENARIO A: nAGQ = 1 (Laplace)
## ============================================================
cat("================================================================\n")
cat("SCENARIO A: nAGQ = 1 (Laplace approximation)\n")
cat("================================================================\n\n")

cat("--- lme4 glmer (nAGQ=1) ---\n")
t_lme4_1 <- system.time({
  fit_lme4_1 <- tryCatch(
    lme4::glmer(y ~ x + (1|g), data = dat, family = poisson, nAGQ = 1),
    error = function(e) e
  )
})
if (inherits(fit_lme4_1, "condition")) {
  cat("lme4 ERROR:", conditionMessage(fit_lme4_1), "\n")
} else {
  conv_msgs <- tryCatch(
    unlist(lapply(fit_lme4_1@optinfo$conv$lme4, `[[`, "messages")),
    error = function(e) character(0))
  cat("lme4 converged:", !any(grepl("failed to converge", conv_msgs)), "\n")
  cat("fixef:", format(lme4::fixef(fit_lme4_1), digits=8), "\n")
  cat("theta:", format(lme4::getME(fit_lme4_1, "theta"), digits=8), "\n")
  cat("logLik:", format(as.numeric(stats::logLik(fit_lme4_1)), digits=10), "\n")
  cat("AIC:", format(stats::AIC(fit_lme4_1), digits=10), "\n")
  cat("BIC:", format(stats::BIC(fit_lme4_1), digits=10), "\n")
  cat("sigma:", format(stats::sigma(fit_lme4_1), digits=8), "\n")
  cat("nobs:", stats::nobs(fit_lme4_1), "  nAGQ:", fit_lme4_1@devcomp$dims["nAGQ"], "\n")
  cat("Wall time:", t_lme4_1["elapsed"], "s\n")
}

cat("\n--- mixeff glmm (nAGQ=1) ---\n")
t_mm_1 <- system.time({
  fit_mm_1 <- tryCatch(
    suppressMessages(mixeff::glmm(y ~ x + (1|g), data = dat,
                                  family = poisson, nAGQ = 1)),
    error = function(e) e
  )
})
if (inherits(fit_mm_1, "condition")) {
  cat("mixeff ERROR:", conditionMessage(fit_mm_1), "\n")
  cat("Class:", paste(class(fit_mm_1), collapse=", "), "\n")
} else {
  cat("fit_status:", fit_mm_1$fit_status, "\n")
  cat("method:", fit_mm_1$method, "  nAGQ:", fit_mm_1$nAGQ, "\n")
  cat("fixef:", format(mixeff::fixef(fit_mm_1), digits=8), "\n")
  cat("theta:", format(fit_mm_1$theta, digits=8), "\n")
  cat("logLik:", format(fit_mm_1$logLik, digits=10), "\n")
  cat("AIC:", format(fit_mm_1$AIC, digits=10), "\n")
  cat("BIC:", format(fit_mm_1$BIC, digits=10), "\n")
  cat("sigma:", format(fit_mm_1$sigma, digits=8), "\n")
  cat("nobs:", fit_mm_1$nobs, "\n")
  cat("Wall time:", t_mm_1["elapsed"], "s\n")
}

cat("\n--- Comparison nAGQ=1 ---\n")
if (!inherits(fit_lme4_1, "condition") && !inherits(fit_mm_1, "condition")) {
  fe_lme4 <- lme4::fixef(fit_lme4_1)
  fe_mm   <- mixeff::fixef(fit_mm_1)
  for (nm in intersect(names(fe_lme4), names(fe_mm))) {
    check_tol(fe_lme4[nm] - fe_mm[nm], TOL_FIXEF, paste0("fixef[", nm, "]"))
  }
  check_tol(lme4::getME(fit_lme4_1, "theta") - fit_mm_1$theta,
            TOL_THETA, "theta (RE sd)")
  check_tol(as.numeric(stats::logLik(fit_lme4_1)) - fit_mm_1$logLik,
            TOL_LOGLIK, "logLik")

  re_lme4  <- lme4::ranef(fit_lme4_1)$g[, 1]
  re_mm    <- mixeff::ranef(fit_mm_1)$g[, 1]
  lev_lme4 <- rownames(lme4::ranef(fit_lme4_1)$g)
  lev_mm   <- rownames(mixeff::ranef(fit_mm_1)$g)
  if (length(re_lme4) == length(re_mm) && setequal(lev_lme4, lev_mm)) {
    re_mm_al <- re_mm[match(lev_lme4, lev_mm)]
    cat(sprintf("  %-36s max_abs_diff=%.2e\n",
                "ranef (max across groups)", max(abs(re_lme4 - re_mm_al))))
  }

  fitted_lme4 <- stats::fitted(fit_lme4_1)
  fitted_mm   <- fit_mm_1$fitted
  if (length(fitted_lme4) == length(fitted_mm)) {
    cat(sprintf("  %-36s max_abs_diff=%.2e\n",
                "fitted values (max)", max(abs(fitted_lme4 - fitted_mm))))
  }

  vc_lme4 <- as.data.frame(lme4::VarCorr(fit_lme4_1))
  vc_mm   <- mixeff::VarCorr(fit_mm_1)$table
  if (nrow(vc_lme4) > 0 && nrow(vc_mm) > 0)
    check_tol(vc_lme4$sdcor[1] - vc_mm$std_dev[1], TOL_SIGMA, "VarCorr SD")
} else {
  cat("  (Cannot compare — one or both fits errored)\n")
}

## ============================================================
## SCENARIO B: nAGQ = 5 (adaptive Gauss-Hermite quadrature)
## ============================================================
cat("\n================================================================\n")
cat("SCENARIO B: nAGQ = 5 (adaptive Gauss-Hermite quadrature)\n")
cat("================================================================\n\n")

cat("--- lme4 glmer (nAGQ=5) ---\n")
t_lme4_5 <- system.time({
  fit_lme4_5 <- tryCatch(
    lme4::glmer(y ~ x + (1|g), data = dat, family = poisson, nAGQ = 5),
    error = function(e) e
  )
})
if (inherits(fit_lme4_5, "condition")) {
  cat("lme4 ERROR:", conditionMessage(fit_lme4_5), "\n")
} else {
  conv_msgs <- tryCatch(
    unlist(lapply(fit_lme4_5@optinfo$conv$lme4, `[[`, "messages")),
    error = function(e) character(0))
  cat("lme4 converged:", !any(grepl("failed to converge", conv_msgs)), "\n")
  cat("fixef:", format(lme4::fixef(fit_lme4_5), digits=8), "\n")
  cat("theta:", format(lme4::getME(fit_lme4_5, "theta"), digits=8), "\n")
  cat("logLik:", format(as.numeric(stats::logLik(fit_lme4_5)), digits=10), "\n")
  cat("AIC:", format(stats::AIC(fit_lme4_5), digits=10), "\n")
  cat("BIC:", format(stats::BIC(fit_lme4_5), digits=10), "\n")
  cat("sigma:", format(stats::sigma(fit_lme4_5), digits=8), "\n")
  cat("nobs:", stats::nobs(fit_lme4_5), "  nAGQ:", fit_lme4_5@devcomp$dims["nAGQ"], "\n")
  cat("Wall time:", t_lme4_5["elapsed"], "s\n")
}

cat("\n--- mixeff glmm (nAGQ=5, method='pirls_profiled') ---\n")
t_mm_5 <- system.time({
  fit_mm_5 <- tryCatch(
    suppressMessages(mixeff::glmm(y ~ x + (1|g), data = dat, family = poisson,
                                  nAGQ = 5, method = "pirls_profiled")),
    error = function(e) e
  )
})
if (inherits(fit_mm_5, "condition")) {
  cat("mixeff ERROR:", conditionMessage(fit_mm_5), "\n")
  cat("Class:", paste(class(fit_mm_5), collapse=", "), "\n")
  if (inherits(fit_mm_5, "mm_inference_unavailable") ||
      inherits(fit_mm_5, "mm_fit_error") ||
      inherits(fit_mm_5, "mm_arg_error")) {
    cat("==> TYPED REFUSAL\n")
  }
} else {
  cat("fit_status:", fit_mm_5$fit_status, "\n")
  cat("method:", fit_mm_5$method, "  nAGQ:", fit_mm_5$nAGQ, "\n")
  cat("fixef:", format(mixeff::fixef(fit_mm_5), digits=8), "\n")
  cat("theta:", format(fit_mm_5$theta, digits=8), "\n")
  cat("logLik:", format(fit_mm_5$logLik, digits=10), "\n")
  cat("AIC:", format(fit_mm_5$AIC, digits=10), "\n")
  cat("BIC:", format(fit_mm_5$BIC, digits=10), "\n")
  cat("sigma:", format(fit_mm_5$sigma, digits=8), "\n")
  cat("nobs:", fit_mm_5$nobs, "\n")
  cat("Wall time:", t_mm_5["elapsed"], "s\n")
}

cat("\n--- mixeff glmm (nAGQ=5, method='joint_laplace') — expect refusal ---\n")
fit_mm_5_jl <- tryCatch(
  suppressMessages(mixeff::glmm(y ~ x + (1|g), data = dat, family = poisson,
                                nAGQ = 5, method = "joint_laplace")),
  error = function(e) e
)
if (inherits(fit_mm_5_jl, "condition")) {
  cat("Result: ERROR (expected)\n")
  cat("Class:", paste(class(fit_mm_5_jl), collapse=", "), "\n")
  cat("Message:", conditionMessage(fit_mm_5_jl), "\n")
  typed <- inherits(fit_mm_5_jl, "mm_arg_error") ||
           inherits(fit_mm_5_jl, "mm_fit_error") ||
           inherits(fit_mm_5_jl, "mm_inference_unavailable")
  cat(if (typed) "==> TYPED REFUSAL: honest and labelled\n"
      else "==> UNTYPED error\n")
} else {
  cat("Result: FIT (unexpected — should have refused)\n")
}

## ============================================================
## SCENARIO C: Cross-nAGQ comparison
## ============================================================
cat("\n================================================================\n")
cat("SCENARIO C: Cross-nAGQ comparison\n")
cat("================================================================\n\n")

cat("lme4 nAGQ=1 vs nAGQ=5:\n")
if (!inherits(fit_lme4_1, "condition") && !inherits(fit_lme4_5, "condition")) {
  fe1 <- lme4::fixef(fit_lme4_1); fe5 <- lme4::fixef(fit_lme4_5)
  for (nm in intersect(names(fe1), names(fe5)))
    cat(sprintf("  fixef[%s]: nAGQ1=%.8f  nAGQ5=%.8f  diff=%.2e\n",
        nm, fe1[nm], fe5[nm], fe1[nm] - fe5[nm]))
  th1 <- lme4::getME(fit_lme4_1, "theta"); th5 <- lme4::getME(fit_lme4_5, "theta")
  cat(sprintf("  theta:  nAGQ1=%.8f  nAGQ5=%.8f  diff=%.2e\n", th1, th5, th1-th5))
  ll1 <- as.numeric(stats::logLik(fit_lme4_1))
  ll5 <- as.numeric(stats::logLik(fit_lme4_5))
  cat(sprintf("  logLik: nAGQ1=%.7f  nAGQ5=%.7f  diff=%.2e\n", ll1, ll5, ll1-ll5))
}

if (!inherits(fit_mm_1, "condition") && !inherits(fit_mm_5, "condition")) {
  cat("\nmixeff nAGQ=1 vs nAGQ=5 (pirls_profiled):\n")
  fe1 <- mixeff::fixef(fit_mm_1); fe5 <- mixeff::fixef(fit_mm_5)
  for (nm in intersect(names(fe1), names(fe5)))
    cat(sprintf("  fixef[%s]: nAGQ1=%.8f  nAGQ5=%.8f  diff=%.2e\n",
        nm, fe1[nm], fe5[nm], fe1[nm] - fe5[nm]))
  cat(sprintf("  theta:  nAGQ1=%.8f  nAGQ5=%.8f  diff=%.2e\n",
    fit_mm_1$theta, fit_mm_5$theta, fit_mm_1$theta - fit_mm_5$theta))
  cat(sprintf("  logLik: nAGQ1=%.7f  nAGQ5=%.7f  diff=%.2e\n",
    fit_mm_1$logLik, fit_mm_5$logLik, fit_mm_1$logLik - fit_mm_5$logLik))
}

## ============================================================
## SCENARIO D: mixeff(nAGQ=5) vs lme4(nAGQ=5) parity
## ============================================================
cat("\n================================================================\n")
cat("SCENARIO D: mixeff(nAGQ=5,pirls_profiled) vs lme4(nAGQ=5) parity\n")
cat("================================================================\n\n")

if (!inherits(fit_lme4_5, "condition") && !inherits(fit_mm_5, "condition")) {
  fe_lme4 <- lme4::fixef(fit_lme4_5); fe_mm <- mixeff::fixef(fit_mm_5)
  for (nm in intersect(names(fe_lme4), names(fe_mm)))
    check_tol(fe_lme4[nm] - fe_mm[nm], TOL_FIXEF, paste0("fixef[", nm, "]"))
  check_tol(lme4::getME(fit_lme4_5, "theta") - fit_mm_5$theta, TOL_THETA, "theta")
  check_tol(as.numeric(stats::logLik(fit_lme4_5)) - fit_mm_5$logLik,
            TOL_LOGLIK, "logLik")
  cat("\n  NOTE: lme4=joint-AGQ (response constants included);\n")
  cat("  mixeff=profiled-PIRLS-AGQ (response constants dropped).\n")
  cat("  logLik gap is by documented convention; fixef/theta parity is the\n")
  cat("  release gate. See glmm_support_contract.md §Approximation Semantics.\n")
} else if (!inherits(fit_lme4_5, "condition") && inherits(fit_mm_5, "condition")) {
  cat("mixeff nAGQ=5 unavailable; fallback: mixeff nAGQ=1 vs lme4 nAGQ=5:\n")
  if (!inherits(fit_mm_1, "condition")) {
    fe_lme4 <- lme4::fixef(fit_lme4_5); fe_mm <- mixeff::fixef(fit_mm_1)
    for (nm in intersect(names(fe_lme4), names(fe_mm)))
      check_tol(fe_lme4[nm] - fe_mm[nm], TOL_FIXEF,
                paste0("fixef[", nm, "] (mm nAGQ1 vs lme4 nAGQ5)"))
    check_tol(lme4::getME(fit_lme4_5, "theta") - fit_mm_1$theta,
              TOL_THETA, "theta (mm1 vs lme4_5)")
    check_tol(as.numeric(stats::logLik(fit_lme4_5)) - fit_mm_1$logLik,
              TOL_LOGLIK, "logLik (mm1 vs lme4_5)")
  }
} else {
  cat("  (Cannot compare — one or both fits failed)\n")
}

## ============================================================
## SCENARIO E: Timing
## ============================================================
cat("\n================================================================\n")
cat("SCENARIO E: Timing\n")
cat("================================================================\n\n")
cat(sprintf("lme4   nAGQ=1: %.3f s\n", t_lme4_1["elapsed"]))
cat(sprintf("mixeff nAGQ=1: %.3f s\n", t_mm_1["elapsed"]))
cat(sprintf("lme4   nAGQ=5: %.3f s\n", t_lme4_5["elapsed"]))
if (!inherits(fit_mm_5, "condition"))
  cat(sprintf("mixeff nAGQ=5: %.3f s\n", t_mm_5["elapsed"]))
if (t_lme4_1["elapsed"] > 0 && t_mm_1["elapsed"] > 0) {
  ratio <- t_lme4_1["elapsed"] / t_mm_1["elapsed"]
  cat(sprintf("Speed ratio nAGQ=1 (lme4/mixeff): %.2fx  (%s)\n",
    ratio, if (ratio > 1) "mixeff faster" else "lme4 faster"))
}

## ============================================================
## SUMMARY
## ============================================================
cat("\n================================================================\n")
cat("SUMMARY\n")
cat("================================================================\n\n")
cat("nAGQ=1 (Laplace):\n")
cat("  lme4 fit:", if (inherits(fit_lme4_1, "condition")) "ERROR" else "OK", "\n")
cat("  mixeff fit:", if (inherits(fit_mm_1, "condition")) "ERROR" else "OK", "\n")
if (!inherits(fit_mm_1, "condition"))
  cat("  mixeff method:", fit_mm_1$method, "nAGQ:", fit_mm_1$nAGQ, "\n")

cat("\nnAGQ=5 (AGQ):\n")
cat("  lme4 fit:", if (inherits(fit_lme4_5, "condition")) "ERROR" else "OK", "\n")
cat("  mixeff fit:", if (inherits(fit_mm_5, "condition"))
  paste("ERROR:", conditionMessage(fit_mm_5)) else "OK", "\n")
if (!inherits(fit_mm_5, "condition"))
  cat("  mixeff method:", fit_mm_5$method, "nAGQ:", fit_mm_5$nAGQ, "\n")

cat("\nnAGQ=5 + method=joint_laplace (expected refusal):\n")
cat("  mixeff:", if (inherits(fit_mm_5_jl, "condition"))
  "REFUSED (expected)" else "FITTED (unexpected)", "\n")
if (inherits(fit_mm_5_jl, "condition"))
  cat("  Class:", paste(class(fit_mm_5_jl), collapse=", "), "\n")

cat("\nDone.\n")
