## Empirical parity probe: lmm-contrasts
## Dataset: simulated factor with 3+ levels
## Formula: y ~ fac + (1|g)   REML = TRUE
## Focus: factor coding parity — treatment (default), sum/contr.sum, helmert,
##         and whether mixeff accepts a `contrasts` arg like lme4 does.
##
## Quantities compared: fixef, SE, vcov, theta, sigma, logLik, AIC, BIC,
##   ranef, fitted, convergence, and what happens when the lme4-style
##   `contrasts=` argument is passed to mixeff::lmm().

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

## ── 0. Session info ────────────────────────────────────────────────────────────
cat("=== SESSION INFO ===\n")
cat("lme4 version:",    as.character(packageVersion("lme4")),    "\n")
cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff version:",  as.character(packageVersion("mixeff")),  "\n\n")

## ── 1. Simulate data ───────────────────────────────────────────────────────────
set.seed(42)
n_groups  <- 12      # random-intercept groups
n_per_grp <- 10      # observations per group
n_levels  <- 4       # factor levels A/B/C/D

g   <- factor(rep(paste0("g", seq_len(n_groups)), each = n_per_grp))
fac <- factor(rep(rep(LETTERS[seq_len(n_levels)], length.out = n_per_grp),
                  n_groups))
b_int <- 5.0
b_fac <- c(0, 2, -1, 3)          # treatment contrasts vs reference "A"
re    <- rnorm(n_groups, 0, 1.5)  # random intercepts
eps   <- rnorm(length(g), 0, 1.0)
y     <- b_int +
         b_fac[as.integer(fac)] +
         re[as.integer(g)] +
         eps

df <- data.frame(y = y, fac = fac, g = g)

cat("=== DATASET ===\n")
cat("nrow:", nrow(df), "  ncol:", ncol(df), "\n")
cat("factor levels:", levels(df$fac), "\n")
cat("groups:", nlevels(df$g), "\n\n")

## ── helper: print + compare ────────────────────────────────────────────────────
compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  diff     <- max(abs(lme4_val - mm_val))
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-38s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse = ", "),
              paste(sprintf(fmt, mm_val),   collapse = ", "),
              diff, tol, status))
  invisible(diff)
}

tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

## ══════════════════════════════════════════════════════════════════════════════
## SCENARIO 1: Default treatment contrasts (R default)
## ══════════════════════════════════════════════════════════════════════════════
cat("=== SCENARIO 1: treatment contrasts (default) ===\n\n")

cat("--- lme4 fit ---\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(y ~ fac + (1 | g), data = df, REML = TRUE)
})
cat("lme4 wall-clock:", t_lme4["elapsed"], "s\n\n")

cat("-- fixef --\n");   print(lme4::fixef(fit_lme4))
cat("-- SE --\n");      print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov --\n");    print(as.matrix(vcov(fit_lme4)))
cat("-- VarCorr --\n"); print(lme4::VarCorr(fit_lme4))
cat("-- sigma --\n");   print(stats::sigma(fit_lme4))
cat("-- logLik --\n");  print(logLik(fit_lme4))
cat("-- AIC --\n");     print(AIC(fit_lme4))
cat("-- BIC --\n");     print(BIC(fit_lme4))
cat("-- ranef (first 6) --\n"); print(head(lme4::ranef(fit_lme4)$g))
cat("-- fitted (first 6) --\n"); print(head(fitted(fit_lme4)))
cat("-- convergence --\n")
cat("No conv warning:", length(fit_lme4@optinfo$conv$lme4$messages) == 0, "\n\n")

cat("--- mixeff fit ---\n")
t_mm <- system.time({
  fit_mm <- tryCatch(
    lmm(y ~ fac + (1 | g), data = df, REML = TRUE,
        control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff wall-clock:", t_mm["elapsed"], "s\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR !!!\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm), "\n\n")
  fit_mm_ok <- FALSE
} else {
  fit_mm_ok <- TRUE
  cat("fit_status:", fit_mm$fit_status, "\n\n")
  cat("-- fixef --\n");   print(fixef(fit_mm))
  cat("-- SE --\n");      print(fit_mm$std_errors)
  cat("-- vcov --\n");    print(as.matrix(fit_mm$fixed_effect_vcov))
  cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
  cat("-- sigma --\n");   print(sigma(fit_mm))
  cat("-- logLik --\n");  print(logLik(fit_mm))
  cat("-- AIC --\n");     print(fit_mm$AIC)
  cat("-- BIC --\n");     print(fit_mm$BIC)
  cat("-- ranef (first 6) --\n"); print(head(ranef(fit_mm)$g))
  cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
  cat("\n")
}

cat("--- numerical comparison (treatment contrasts) ---\n")
if (fit_mm_ok) {
  fe_lme4 <- lme4::fixef(fit_lme4)
  fe_mm   <- fixef(fit_mm)
  # mixeff uses "fac: B" while lme4 uses "facB" — compare by position
  cat("lme4 coef names:", paste(names(fe_lme4), collapse=", "), "\n")
  cat("mm   coef names:", paste(names(fe_mm),   collapse=", "), "\n\n")
  nms_lme4 <- names(fe_lme4)
  nms_mm   <- names(fe_mm)
  for (i in seq_along(nms_lme4)) {
    compare(paste0("fixef[", i, "] ", nms_lme4[i]),
            fe_lme4[[i]], fe_mm[[i]], tols$fixef)
  }

  se_lme4 <- sqrt(diag(vcov(fit_lme4)))
  se_mm   <- fit_mm$std_errors
  for (i in seq_along(nms_lme4)) {
    compare(paste0("SE[", i, "] ", nms_lme4[i]),
            se_lme4[[i]], se_mm[[i]], tols$fixef)
  }

  # vcov full matrix
  vc_lme4_mat <- as.matrix(vcov(fit_lme4))
  vc_mm_mat   <- as.matrix(fit_mm$fixed_effect_vcov)
  vcov_diff   <- max(abs(vc_lme4_mat - vc_mm_mat))
  cat(sprintf("%-38s maxAbsDiff=%.3e  tol=1e-04  [%s]\n",
              "vcov max abs diff",
              vcov_diff, if (vcov_diff <= 1e-4) "WITHIN-TOL" else "EXCEEDS-TOL"))

  # theta
  theta_lme4 <- lme4::getME(fit_lme4, "theta")
  theta_mm   <- fit_mm$theta
  compare("theta", theta_lme4, theta_mm, tols$theta)

  compare("sigma",
          stats::sigma(fit_lme4), sigma(fit_mm), tols$sigma)

  # VarCorr (g random intercept variance)
  vc_lme4_var <- as.numeric(attr(lme4::VarCorr(fit_lme4)$g, "stddev"))^2
  vc_mm_obj   <- VarCorr(fit_mm)
  vc_mm_var   <- tryCatch({
    if (inherits(vc_mm_obj, "mm_varcorr") && is.data.frame(vc_mm_obj$table)) {
      tbl <- vc_mm_obj$table
      as.numeric(tbl[tbl$group == "g", "variance"])
    } else NA_real_
  }, error = function(e) NA_real_)
  cat(sprintf("%-38s lme4=%.8f  mm=%s\n",
              "VarCorr g variance",
              vc_lme4_var,
              ifelse(is.na(vc_mm_var), "NA", sprintf("%.8f", vc_mm_var))))
  if (!is.na(vc_mm_var))
    compare("VarCorr g variance", vc_lme4_var, vc_mm_var, tols$theta)

  compare("logLik",
          as.numeric(logLik(fit_lme4)),
          as.numeric(logLik(fit_mm)),
          tols$logLik)
  compare("AIC",
          AIC(fit_lme4), fit_mm$AIC, tols$logLik * 2)
  compare("BIC",
          BIC(fit_lme4), fit_mm$BIC, tols$logLik * 2)

  # fitted
  fitted_lme4 <- as.numeric(fitted(fit_lme4))
  fitted_mm   <- as.numeric(fit_mm$fitted)
  compare("fitted max abs diff",
          fitted_lme4, fitted_mm[seq_along(fitted_lme4)], tols$fixef)

  # ranef (sorted, order-invariant)
  re_lme4 <- sort(as.numeric(lme4::ranef(fit_lme4)$g[, 1]))
  re_mm_df <- ranef(fit_mm)$g
  re_mm    <- if (!is.null(re_mm_df)) sort(as.numeric(re_mm_df[, 1])) else NA_real_
  compare("ranef g max abs diff",
          re_lme4, re_mm, tols$fixef)

  # speed
  cat(sprintf("\nwall-clock  lme4=%.4fs  mm=%.4fs  ratio=%.2fx\n",
              t_lme4["elapsed"], t_mm["elapsed"],
              t_mm["elapsed"] / t_lme4["elapsed"]))
} else {
  cat("(skipped — mixeff errored)\n")
}

## ══════════════════════════════════════════════════════════════════════════════
## SCENARIO 2: Sum contrasts (contr.sum) pre-set on data factor
## ══════════════════════════════════════════════════════════════════════════════
cat("\n=== SCENARIO 2: sum contrasts (contr.sum) pre-set on data factor ===\n\n")

df2 <- df
contrasts(df2$fac) <- contr.sum(n_levels)

cat("--- lme4 fit (contr.sum via contrasts attribute on factor) ---\n")
t_lme4_s <- system.time({
  fit_lme4_s <- lmer(y ~ fac + (1 | g), data = df2, REML = TRUE)
})
cat("lme4 wall-clock:", t_lme4_s["elapsed"], "s\n")
cat("-- fixef --\n"); print(lme4::fixef(fit_lme4_s))
cat("-- SE --\n");    print(sqrt(diag(vcov(fit_lme4_s))))
cat("-- logLik --\n"); print(logLik(fit_lme4_s))
cat("\n")

cat("--- lme4 fit (contr.sum via contrasts= argument) ---\n")
t_lme4_s2 <- system.time({
  fit_lme4_s2 <- lmer(y ~ fac + (1 | g), data = df,
                      contrasts = list(fac = "contr.sum"), REML = TRUE)
})
cat("lme4 wall-clock:", t_lme4_s2["elapsed"], "s\n")
cat("-- fixef --\n"); print(lme4::fixef(fit_lme4_s2))
cat("-- SE --\n");    print(sqrt(diag(vcov(fit_lme4_s2))))
cat("-- logLik --\n"); print(logLik(fit_lme4_s2))
cat("\n")

cat("--- mixeff fit (sum contrasts via pre-set attribute on factor in df2) ---\n")
t_mm_s <- system.time({
  fit_mm_s <- tryCatch(
    lmm(y ~ fac + (1 | g), data = df2, REML = TRUE,
        control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff wall-clock:", t_mm_s["elapsed"], "s\n\n")

if (inherits(fit_mm_s, "condition")) {
  cat("!!! mixeff ERROR (sum contrasts via factor attribute) !!!\n")
  cat("class:", paste(class(fit_mm_s), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm_s), "\n\n")
  fit_mm_s_ok <- FALSE
} else {
  fit_mm_s_ok <- TRUE
  cat("fit_status:", fit_mm_s$fit_status, "\n")
  cat("-- fixef --\n"); print(fixef(fit_mm_s))
  cat("-- SE --\n");    print(fit_mm_s$std_errors)
  cat("-- logLik --\n"); print(logLik(fit_mm_s))
  cat("\n")
}

cat("--- mixeff: contrasts= argument (like lme4) ---\n")
# Does mixeff::lmm() accept a contrasts= argument?
has_contrasts_arg <- "contrasts" %in% names(formals(mixeff::lmm))
cat("mixeff::lmm has 'contrasts' formal argument:", has_contrasts_arg, "\n")

fit_mm_carg <- tryCatch(
  lmm(y ~ fac + (1 | g), data = df,
      contrasts = list(fac = "contr.sum"),
      REML = TRUE,
      control = mm_control(verbose = -1L)),
  error = function(e) e
)
if (inherits(fit_mm_carg, "condition")) {
  cat("mixeff + contrasts= arg RESULT: ERROR\n")
  cat("class:", paste(class(fit_mm_carg), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm_carg), "\n")
} else {
  cat("mixeff + contrasts= arg RESULT: fit succeeded\n")
  cat("fit_status:", fit_mm_carg$fit_status, "\n")
  cat("-- fixef --\n"); print(fixef(fit_mm_carg))
}
cat("\n")

cat("--- numerical comparison (sum contrasts via factor attribute) ---\n")
if (fit_mm_s_ok) {
  fe_lme4_s <- lme4::fixef(fit_lme4_s)
  fe_mm_s   <- fixef(fit_mm_s)
  cat("lme4(sum) coef names:", paste(names(fe_lme4_s), collapse=", "), "\n")
  cat("mm(sum)   coef names:", paste(names(fe_mm_s),   collapse=", "), "\n\n")
  for (i in seq_along(fe_lme4_s)) {
    compare(paste0("fixef(sum)[", i, "] ", names(fe_lme4_s)[i]),
            fe_lme4_s[[i]], fe_mm_s[[i]], tols$fixef)
  }
  compare("SE(sum) max abs diff",
          sqrt(diag(vcov(fit_lme4_s))),
          as.numeric(fit_mm_s$std_errors),
          tols$fixef)
  compare("logLik(sum)",
          as.numeric(logLik(fit_lme4_s)),
          as.numeric(logLik(fit_mm_s)),
          tols$logLik)
  compare("sigma(sum)",
          stats::sigma(fit_lme4_s), sigma(fit_mm_s), tols$sigma)

  # Check: sum-contrast logLik == treatment-contrast logLik (invariant)
  ll_treat <- as.numeric(logLik(fit_lme4))
  ll_sum   <- as.numeric(logLik(fit_lme4_s))
  cat(sprintf("\nlme4 logLik: treatment=%.6f  sum=%.6f  diff=%.3e  [should be ~0]\n",
              ll_treat, ll_sum, abs(ll_treat - ll_sum)))
  if (fit_mm_ok && fit_mm_s_ok) {
    ll_mm_treat <- as.numeric(logLik(fit_mm))
    ll_mm_sum   <- as.numeric(logLik(fit_mm_s))
    cat(sprintf("mm  logLik: treatment=%.6f  sum=%.6f  diff=%.3e  [should be ~0]\n",
                ll_mm_treat, ll_mm_sum, abs(ll_mm_treat - ll_mm_sum)))
  }
} else {
  cat("(skipped — mixeff errored with sum contrasts)\n")
}

## ══════════════════════════════════════════════════════════════════════════════
## SCENARIO 3: Helmert contrasts (contr.helmert)
## ══════════════════════════════════════════════════════════════════════════════
cat("\n=== SCENARIO 3: Helmert contrasts (contr.helmert) ===\n\n")

df3 <- df
contrasts(df3$fac) <- contr.helmert(n_levels)

cat("--- lme4 fit (contr.helmert) ---\n")
fit_lme4_h <- lmer(y ~ fac + (1 | g), data = df3, REML = TRUE)
cat("-- fixef --\n"); print(lme4::fixef(fit_lme4_h))
cat("-- logLik --\n"); print(logLik(fit_lme4_h))
cat("\n")

cat("--- mixeff fit (contr.helmert via pre-set attribute) ---\n")
fit_mm_h <- tryCatch(
  lmm(y ~ fac + (1 | g), data = df3, REML = TRUE,
      control = mm_control(verbose = -1L)),
  error = function(e) e
)
if (inherits(fit_mm_h, "condition")) {
  cat("!!! mixeff ERROR (Helmert) !!!\n")
  cat("class:", paste(class(fit_mm_h), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm_h), "\n\n")
  fit_mm_h_ok <- FALSE
} else {
  fit_mm_h_ok <- TRUE
  cat("fit_status:", fit_mm_h$fit_status, "\n")
  cat("-- fixef --\n"); print(fixef(fit_mm_h))
  cat("-- logLik --\n"); print(logLik(fit_mm_h))
  cat("\n")
}

cat("--- numerical comparison (Helmert) ---\n")
if (fit_mm_h_ok) {
  fe_lme4_h <- lme4::fixef(fit_lme4_h)
  fe_mm_h   <- fixef(fit_mm_h)
  cat("lme4(helm) coef names:", paste(names(fe_lme4_h), collapse=", "), "\n")
  cat("mm(helm)   coef names:", paste(names(fe_mm_h),   collapse=", "), "\n\n")
  for (i in seq_along(fe_lme4_h)) {
    compare(paste0("fixef(helm)[", i, "] ", names(fe_lme4_h)[i]),
            fe_lme4_h[[i]], fe_mm_h[[i]], tols$fixef)
  }
  compare("logLik(helm)",
          as.numeric(logLik(fit_lme4_h)),
          as.numeric(logLik(fit_mm_h)),
          tols$logLik)
  compare("sigma(helm)",
          stats::sigma(fit_lme4_h), sigma(fit_mm_h), tols$sigma)
} else {
  cat("(skipped — mixeff errored with Helmert contrasts)\n")
}

## ══════════════════════════════════════════════════════════════════════════════
## SCENARIO 4: 5-level factor (more levels)
## ══════════════════════════════════════════════════════════════════════════════
cat("\n=== SCENARIO 4: 5-level factor (more levels, treatment coding) ===\n\n")

set.seed(99)
n5    <- 5
fac5  <- factor(rep(LETTERS[1:n5], length.out = nrow(df)))
b5    <- c(0, 1.5, -2, 3, -0.5)
y5    <- b_int + b5[as.integer(fac5)] + re[as.integer(g)] + rnorm(nrow(df), 0, 1)
df5   <- data.frame(y = y5, fac = fac5, g = g)

fit_lme4_5 <- lmer(y ~ fac + (1 | g), data = df5, REML = TRUE)
cat("lme4 fixef:", round(lme4::fixef(fit_lme4_5), 5), "\n")

fit_mm_5 <- tryCatch(
  lmm(y ~ fac + (1 | g), data = df5, REML = TRUE,
      control = mm_control(verbose = -1L)),
  error = function(e) e
)
if (inherits(fit_mm_5, "condition")) {
  cat("!!! mixeff ERROR (5-level factor) !!!\n")
  cat("class:", paste(class(fit_mm_5), collapse = ", "), "\n")
  cat(conditionMessage(fit_mm_5), "\n\n")
  fit_mm_5_ok <- FALSE
} else {
  fit_mm_5_ok <- TRUE
  cat("mm    fixef:", round(fixef(fit_mm_5), 5), "\n")
}

cat("--- numerical comparison (5 levels) ---\n")
if (fit_mm_5_ok) {
  fe_lme4_5 <- lme4::fixef(fit_lme4_5)
  fe_mm_5   <- fixef(fit_mm_5)
  cat("lme4(5lev) coef names:", paste(names(fe_lme4_5), collapse=", "), "\n")
  cat("mm(5lev)   coef names:", paste(names(fe_mm_5),   collapse=", "), "\n\n")
  for (i in seq_along(fe_lme4_5)) {
    compare(paste0("fixef5[", i, "] ", names(fe_lme4_5)[i]),
            fe_lme4_5[[i]], fe_mm_5[[i]], tols$fixef)
  }
  compare("logLik5",
          as.numeric(logLik(fit_lme4_5)),
          as.numeric(logLik(fit_mm_5)),
          tols$logLik)
  compare("sigma5",
          stats::sigma(fit_lme4_5), sigma(fit_mm_5), tols$sigma)
} else {
  cat("(skipped — mixeff errored)\n")
}

## ══════════════════════════════════════════════════════════════════════════════
## OVERALL SUMMARY
## ══════════════════════════════════════════════════════════════════════════════
cat("\n=== OVERALL SUMMARY ===\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
cat("Scenario 1 (treatment, 4 levels): fit_mm_ok =", fit_mm_ok, "\n")
cat("Scenario 2 (sum contrasts):       fit_mm_s_ok =",
    if (exists("fit_mm_s_ok")) fit_mm_s_ok else FALSE, "\n")
cat("Scenario 3 (Helmert):             fit_mm_h_ok =",
    if (exists("fit_mm_h_ok")) fit_mm_h_ok else FALSE, "\n")
cat("Scenario 4 (5-level):             fit_mm_5_ok =",
    if (exists("fit_mm_5_ok")) fit_mm_5_ok else FALSE, "\n")
cat("mixeff::lmm has contrasts= arg:", has_contrasts_arg, "\n")
