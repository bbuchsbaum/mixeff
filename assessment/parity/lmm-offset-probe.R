## lmm-offset parity probe
## Cell: lmm-offset  Formula: y ~ x + offset(o) + (1|g)
## Tests whether mixeff supports offset() in lmm()

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

set.seed(42)
n_groups <- 20
n_per_group <- 5
n <- n_groups * n_per_group

g <- factor(rep(seq_len(n_groups), each = n_per_group))
x <- rnorm(n)
o <- rnorm(n, mean = 0.5, sd = 0.3)   # offset values (e.g. log exposure)
b_re <- rnorm(n_groups, sd = 0.8)     # true random intercepts
y <- 1.5 + 0.7 * x + o + b_re[as.integer(g)] + rnorm(n, sd = 0.5)

dat <- data.frame(y = y, x = x, o = o, g = g)

cat("=== Dataset summary ===\n")
cat(sprintf("n=%d, n_groups=%d\n", n, n_groups))
cat(sprintf("y: mean=%.3f, sd=%.3f\n", mean(dat$y), sd(dat$y)))
cat(sprintf("x: mean=%.3f, sd=%.3f\n", mean(dat$x), sd(dat$x)))
cat(sprintf("o: mean=%.3f, sd=%.3f\n", mean(dat$o), sd(dat$o)))
cat("\n")

## ─────────────────────────────────────────────
## 1. lme4/lmerTest fit
## ─────────────────────────────────────────────
cat("=== lme4 fit ===\n")
lme4_fit <- tryCatch(
  lmerTest::lmer(y ~ x + offset(o) + (1 | g), data = dat, REML = TRUE),
  error = function(e) { cat("lme4 ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(lme4_fit)) {
  cat("lme4 converged: TRUE\n")
  fe_lme4_print <- lme4::fixef(lme4_fit)
  cat("fixef:\n"); print(fe_lme4_print)
  cat("SE:\n"); print(sqrt(diag(as.matrix(vcov(lme4_fit)))))
  cat("theta:", lme4::getME(lme4_fit, "theta"), "\n")
  sigma_lme4 <- sigma(lme4_fit)
  cat("sigma:", sigma_lme4, "\n")
  ll_lme4 <- logLik(lme4_fit)
  cat("logLik:", as.numeric(ll_lme4), "\n")
  cat("AIC:", AIC(lme4_fit), "\n")
  cat("BIC:", BIC(lme4_fit), "\n")
  ranef_lme4 <- lme4::ranef(lme4_fit)$g[, 1]
  cat("ranef(g) head:", head(ranef_lme4, 5), "\n")
  fitted_lme4 <- fitted(lme4_fit)
  cat("fitted head:", head(fitted_lme4, 5), "\n")
} else {
  cat("lme4 fit FAILED\n")
}

cat("\n")

## ─────────────────────────────────────────────
## 2. mixeff fit
## ─────────────────────────────────────────────
cat("=== mixeff fit ===\n")
mm_fit <- tryCatch(
  mixeff::lmm(y ~ x + offset(o) + (1 | g), data = dat,
              REML = TRUE, control = mm_control(verbose = -1)),
  error = function(e) {
    cat("mixeff ERROR:", conditionMessage(e), "\n")
    cat("mixeff ERROR class:", paste(class(e), collapse = ", "), "\n")
    e
  }
)

mixeff_ok <- !inherits(mm_fit, "condition")

if (mixeff_ok) {
  cat("mixeff converged/fit_status:", mm_fit$fit_status, "\n")
  cat("fixef:\n"); print(mm_fit$beta)
  cat("SE:\n"); print(mm_fit$std_errors)
  cat("theta:", mm_fit$theta, "\n")
  cat("sigma:", mm_fit$sigma, "\n")
  cat("logLik:", mm_fit$logLik, "\n")
  cat("AIC:", mm_fit$AIC, "\n")
  cat("BIC:", mm_fit$BIC, "\n")
  ranef_mm <- unlist(mm_fit$random_effects)
  cat("ranef head:", head(ranef_mm, 5), "\n")
  cat("fitted head:", head(mm_fit$fitted, 5), "\n")
}

cat("\n")

## ─────────────────────────────────────────────
## 3. Comparison
## ─────────────────────────────────────────────
cat("=== Comparison ===\n")

TOL_FIXEF  <- 1e-4
TOL_THETA  <- 1e-3
TOL_LOGLIK <- 1e-3
TOL_SIGMA  <- 1e-4

if (!is.null(lme4_fit) && mixeff_ok) {
  ## fixef
  fe_lme4 <- lme4::fixef(lme4_fit)
  fe_mm   <- mm_fit$beta
  # align by name
  common_names <- intersect(names(fe_lme4), names(fe_mm))
  if (length(common_names) > 0) {
    diff_fixef <- abs(fe_lme4[common_names] - fe_mm[common_names])
    cat("fixef diffs (absolute):\n"); print(diff_fixef)
    cat(sprintf("  max |diff| fixef = %.2e  (tol %.0e) => %s\n",
                max(diff_fixef), TOL_FIXEF,
                ifelse(max(diff_fixef) <= TOL_FIXEF, "WITHIN-TOL", "DIVERGENT")))
  } else {
    cat("WARNING: no common fixef names to compare\n")
    cat("  lme4 names:", paste(names(fe_lme4), collapse=", "), "\n")
    cat("  mixeff names:", paste(names(fe_mm), collapse=", "), "\n")
  }

  ## SE
  se_lme4 <- sqrt(diag(as.matrix(vcov(lme4_fit))))
  se_mm   <- mm_fit$std_errors
  common_se <- intersect(names(se_lme4), names(se_mm))
  if (length(common_se) > 0) {
    diff_se <- abs(se_lme4[common_se] - se_mm[common_se])
    cat("SE diffs:\n"); print(diff_se)
    cat(sprintf("  max |diff| SE = %.2e\n", max(diff_se)))
  }

  ## theta
  th_lme4 <- lme4::getME(lme4_fit, "theta")
  th_mm   <- mm_fit$theta
  diff_theta <- abs(th_lme4 - th_mm[seq_along(th_lme4)])
  cat(sprintf("  max |diff| theta = %.2e  (tol %.0e) => %s\n",
              max(diff_theta), TOL_THETA,
              ifelse(max(diff_theta) <= TOL_THETA, "WITHIN-TOL", "DIVERGENT")))

  ## sigma
  diff_sigma <- abs(sigma_lme4 - mm_fit$sigma)
  cat(sprintf("  |diff| sigma = %.2e  (tol %.0e) => %s\n",
              diff_sigma, TOL_SIGMA,
              ifelse(diff_sigma <= TOL_SIGMA, "WITHIN-TOL", "DIVERGENT")))

  ## logLik
  diff_ll <- abs(as.numeric(ll_lme4) - mm_fit$logLik)
  cat(sprintf("  |diff| logLik = %.2e  (tol %.0e) => %s\n",
              diff_ll, TOL_LOGLIK,
              ifelse(diff_ll <= TOL_LOGLIK, "WITHIN-TOL", "DIVERGENT")))

  ## AIC / BIC
  diff_aic <- abs(AIC(lme4_fit) - mm_fit$AIC)
  diff_bic <- abs(BIC(lme4_fit) - mm_fit$BIC)
  cat(sprintf("  |diff| AIC = %.2e\n", diff_aic))
  cat(sprintf("  |diff| BIC = %.2e\n", diff_bic))

  ## ranef
  re_mm_named <- mm_fit$random_effects
  re_mm_vec <- if (is.list(re_mm_named)) unlist(re_mm_named) else re_mm_named
  if (length(re_mm_vec) == length(ranef_lme4)) {
    diff_re <- abs(sort(ranef_lme4) - sort(re_mm_vec))
    cat(sprintf("  max |diff| ranef = %.2e\n", max(diff_re)))
  } else {
    cat(sprintf("  ranef length mismatch: lme4=%d mixeff=%d\n",
                length(ranef_lme4), length(re_mm_vec)))
  }

  ## fitted
  if (length(fitted_lme4) == length(mm_fit$fitted)) {
    diff_fit <- abs(fitted_lme4 - mm_fit$fitted)
    cat(sprintf("  max |diff| fitted = %.2e\n", max(diff_fit)))
  }

  cat("\n=== VERDICT ===\n")
  all_within <- (length(common_names) > 0 && max(diff_fixef) <= TOL_FIXEF &&
                 max(diff_theta) <= TOL_THETA &&
                 diff_sigma <= TOL_SIGMA &&
                 diff_ll <= TOL_LOGLIK)
  cat(if (all_within) "within-tol" else "divergent", "\n")

} else if (!is.null(lme4_fit) && !mixeff_ok) {
  cat("lme4 succeeded, mixeff FAILED/REFUSED\n")
  cat("Error class:", paste(class(mm_fit), collapse=", "), "\n")
  cat("Error message:", conditionMessage(mm_fit), "\n")
  cat("\n=== VERDICT: mixeff-error ===\n")
} else if (is.null(lme4_fit) && !mixeff_ok) {
  cat("Both failed.\n=== VERDICT: both-error ===\n")
} else {
  cat("lme4 FAILED, mixeff succeeded — unexpected.\n")
}

## ─────────────────────────────────────────────
## 4. Also test: offset as column (not formula term)
## ─────────────────────────────────────────────
cat("\n=== Alternative: offset as argument (not in formula) ===\n")
mm_fit2 <- tryCatch(
  mixeff::lmm(y ~ x + (1 | g), data = dat,
              REML = TRUE, control = mm_control(verbose = -1)),
  error = function(e) {
    cat("mixeff no-offset ERROR:", conditionMessage(e), "\n"); NULL
  }
)
if (!is.null(mm_fit2)) {
  cat("mixeff no-offset fit_status:", mm_fit2$fit_status, "\n")
  cat("(This model excludes offset — shown for reference only)\n")
}

## ─────────────────────────────────────────────
## 5. Manifest / capability check
## ─────────────────────────────────────────────
cat("\n=== mixeff formula manifest ===\n")
mfest <- tryCatch(mm_formula_manifest(), error = function(e) NULL)
if (!is.null(mfest)) {
  cat("formula_features$operators:\n")
  print(mfest$formula_features$operators)
  cat("formula_features$transformations:\n")
  print(mfest$formula_features$transformations)
}

cat("\nDone.\n")
