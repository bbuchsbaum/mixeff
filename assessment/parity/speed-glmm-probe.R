## Empirical parity probe: speed-glmm
## Dataset: simulated binomial, N=1e4. Formula: y ~ x + (1|g)
## Family: binomial, link: logit, method: Laplace (nAGQ=1)
## Focus: timing ratio AND numerical parity on all GLMM quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:   ", as.character(packageVersion("lme4")), "\n")
cat("mixeff version: ", as.character(packageVersion("mixeff")), "\n\n")

## ── 1. Simulate data (N=10000, 50 groups) ─────────────────────────────────────
set.seed(42)
N  <- 10000L
ng <- 50L
g  <- factor(sample(seq_len(ng), N, replace = TRUE))
x  <- rnorm(N)
u  <- rnorm(ng, sd = 0.8)[as.integer(g)]   # true random intercept sd = 0.8
eta <- 0.5 + 0.3 * x + u
p   <- plogis(eta)
y   <- rbinom(N, 1L, p)
dat <- data.frame(y = y, x = x, g = g)

cat("=== DATASET ===\n")
cat("N:", N, "  groups:", ng, "  y=1 proportion:", round(mean(y), 4), "\n\n")

## ── 2. Warm-up run (trigger JIT / load shared libs) ─────────────────────────
invisible(suppressWarnings(
  lme4::glmer(y ~ x + (1|g), data = dat, family = binomial, nAGQ = 1L,
              control = lme4::glmerControl(optimizer = "bobyqa"))
))
invisible(tryCatch(
  mixeff::glmm(y ~ x + (1|g), data = dat, family = binomial(link="logit"),
               method = "pirls_profiled", nAGQ = 1L,
               control = mm_control(verbose = -1L)),
  error = function(e) NULL
))

## ── 3. Timed fits (3 reps each) ───────────────────────────────────────────────
NREP <- 3L

cat("=== TIMING (", NREP, "reps each after warm-up) ===\n", sep = "")

## lme4 timing
lme4_times <- numeric(NREP)
for (i in seq_len(NREP)) {
  t0 <- proc.time()["elapsed"]
  suppressWarnings(
    lme4::glmer(y ~ x + (1|g), data = dat, family = binomial, nAGQ = 1L,
                control = lme4::glmerControl(optimizer = "bobyqa"))
  )
  lme4_times[i] <- proc.time()["elapsed"] - t0
}
cat("lme4  times (s):", sprintf("%.4f", lme4_times), "\n")
cat("lme4  median  :", sprintf("%.4f", median(lme4_times)), "s\n")
cat("lme4  mean    :", sprintf("%.4f", mean(lme4_times)), "s\n\n")

## mixeff timing
mm_times <- numeric(NREP)
mm_error <- NULL
for (i in seq_len(NREP)) {
  t0 <- proc.time()["elapsed"]
  fit_tmp <- tryCatch(
    mixeff::glmm(y ~ x + (1|g), data = dat, family = binomial(link="logit"),
                 method = "pirls_profiled", nAGQ = 1L,
                 control = mm_control(verbose = -1L)),
    error = function(e) e
  )
  mm_times[i] <- proc.time()["elapsed"] - t0
  if (inherits(fit_tmp, "condition")) {
    mm_error <- fit_tmp
    break
  }
}

if (!is.null(mm_error)) {
  cat("!!! mixeff ERROR/REFUSAL during timing !!!\n")
  cat("class  :", paste(class(mm_error), collapse = ", "), "\n")
  cat("message:", conditionMessage(mm_error), "\n")
  if (!is.null(mm_error$reason_code)) cat("reason_code:", mm_error$reason_code, "\n")
  cat("\nSpeed ratio: CANNOT COMPUTE (mixeff failed)\n")
  quit(status = 0)
}

cat("mixeff times (s):", sprintf("%.4f", mm_times), "\n")
cat("mixeff median   :", sprintf("%.4f", median(mm_times)), "s\n")
cat("mixeff mean     :", sprintf("%.4f", mean(mm_times)), "s\n\n")

ratio_median <- median(mm_times) / median(lme4_times)
ratio_mean   <- mean(mm_times)   / mean(lme4_times)
cat(sprintf("Speed ratio (mixeff/lme4): median=%.2fx  mean=%.2fx\n",
            ratio_median, ratio_mean))
cat(sprintf("  (ratio < 1 = mixeff FASTER; > 1 = mixeff SLOWER)\n\n"))

## ── 4. Reference fits (final single run for parity) ──────────────────────────
cat("=== REFERENCE FITS (for numerical parity) ===\n")
t_lme4 <- system.time(
  fit_lme4 <- suppressWarnings(
    lme4::glmer(y ~ x + (1|g), data = dat, family = binomial, nAGQ = 1L,
                control = lme4::glmerControl(optimizer = "bobyqa"))
  )
)
cat("lme4 elapsed:", t_lme4["elapsed"], "s\n\n")

t_mm <- system.time(
  fit_mm <- tryCatch(
    mixeff::glmm(y ~ x + (1|g), data = dat, family = binomial(link="logit"),
                 method = "pirls_profiled", nAGQ = 1L,
                 control = mm_control(verbose = -1L)),
    error = function(e) e
  )
)
cat("mixeff elapsed:", t_mm["elapsed"], "s\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff ERROR in reference fit !!!\n")
  cat("class  :", paste(class(fit_mm), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm), "\n")
  if (!is.null(fit_mm$reason_code)) cat("reason_code:", fit_mm$reason_code, "\n")
  quit(status = 0)
}

cat("fit_status:", fit_mm$fit_status, "\n\n")

## ── 5. Print key quantities ───────────────────────────────────────────────────
cat("=== lme4 KEY QUANTITIES ===\n")
cat("-- fixef --\n");         print(lme4::fixef(fit_lme4))
cat("-- SE --\n");            print(sqrt(diag(vcov(fit_lme4))))
cat("-- theta (getME) --\n"); print(lme4::getME(fit_lme4, "theta"))
cat("-- sigma --\n");         print(sigma(fit_lme4))
cat("-- logLik --\n");        print(as.numeric(logLik(fit_lme4)))
cat("-- AIC --\n");           print(AIC(fit_lme4))
cat("-- BIC --\n");           print(BIC(fit_lme4))
cat("-- deviance --\n");      print(deviance(fit_lme4))
cat("-- fitted head --\n");   print(head(fitted(fit_lme4)))
cat("-- ranef head --\n");    print(head(lme4::ranef(fit_lme4)$g))
msgs <- fit_lme4@optinfo$conv$lme4$messages
cat("-- convergence --\n");
cat("  No convergence warning =", length(msgs) == 0, "\n")
if (length(msgs)) cat("  Messages:", paste(msgs, collapse="; "), "\n")
cat("\n")

cat("=== mixeff KEY QUANTITIES ===\n")
cat("-- fixef --\n");   print(fixef(fit_mm))
cat("-- SE --\n");      print(fit_mm$std_errors)
cat("-- theta --\n");   print(fit_mm$theta)
cat("-- sigma --\n");   print(sigma(fit_mm))
cat("-- logLik --\n");  print(logLik(fit_mm))
cat("-- AIC --\n");     print(fit_mm$AIC)
cat("-- BIC --\n");     print(fit_mm$BIC)
cat("-- deviance --\n"); print(fit_mm$deviance)
cat("-- fitted head --\n"); print(head(fit_mm$fitted))
cat("-- ranef head --\n");
re_mm <- tryCatch(ranef(fit_mm), error = function(e) { cat("ranef error:", conditionMessage(e), "\n"); NULL })
if (!is.null(re_mm)) print(head(re_mm$g))
cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
cat("\n")

## ── 6. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  diff     <- max(abs(lme4_val - mm_val), na.rm = TRUE)
  status   <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL*"
  cat(sprintf("%-38s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label, diff, tol, status))
  invisible(diff)
}

## fixef
fe_lme4 <- lme4::fixef(fit_lme4)
fe_mm   <- fixef(fit_mm)
nm <- names(fe_lme4)
fe_mm_al <- fe_mm[nm]
for (n in nm) compare(paste0("fixef[", n, "]"), fe_lme4[[n]], fe_mm_al[[n]], tols$fixef)
d_fixef <- compare("fixef (all)", fe_lme4[nm], fe_mm_al[nm], tols$fixef)

## SE
se_lme4 <- sqrt(diag(vcov(fit_lme4)))
se_mm   <- fit_mm$std_errors[nm]
for (n in nm) compare(paste0("SE[", n, "]"), se_lme4[[n]], se_mm[[n]], tols$fixef)
d_se <- compare("SE (all)", se_lme4[nm], se_mm[nm], tols$fixef)

## theta
theta_lme4 <- lme4::getME(fit_lme4, "theta")
theta_mm   <- fit_mm$theta
d_theta <- compare("theta", theta_lme4, theta_mm, tols$theta)

## sigma
d_sigma <- compare("sigma (dispersion)", sigma(fit_lme4), sigma(fit_mm), tols$sigma)

## VarCorr: group random-intercept variance
vc_lme4 <- as.numeric(lme4::VarCorr(fit_lme4)$g)
vc_mm_obj <- VarCorr(fit_mm)
vc_mm <- tryCatch({
  if (is.list(vc_mm_obj) && !is.null(vc_mm_obj$table)) {
    tbl <- vc_mm_obj$table
    as.numeric(tbl[tbl$group == "g", "variance"])
  } else if (is.data.frame(vc_mm_obj)) {
    as.numeric(vc_mm_obj[vc_mm_obj$group == "g", "variance"])
  } else {
    NA_real_
  }
}, error = function(e) { cat("VarCorr extraction error:", conditionMessage(e), "\n"); NA_real_ })

cat(sprintf("%-38s  lme4=%.6f  mm=%s\n", "VarCorr g variance",
            vc_lme4,
            ifelse(is.na(vc_mm), "NA", sprintf("%.6f", vc_mm))))
d_varcorr <- if (!is.na(vc_mm)) compare("VarCorr g variance", vc_lme4, vc_mm, tols$theta) else NA_real_

## logLik
d_loglik <- compare("logLik", as.numeric(logLik(fit_lme4)),
                    as.numeric(logLik(fit_mm)), tols$logLik)

## deviance
d_deviance <- compare("deviance", deviance(fit_lme4), fit_mm$deviance, tols$logLik * 2)

## AIC
d_aic <- compare("AIC", AIC(fit_lme4), fit_mm$AIC, tols$logLik * 2)

## BIC
d_bic <- compare("BIC", BIC(fit_lme4), fit_mm$BIC, tols$logLik * 2)

## fitted
n_common <- min(length(fitted(fit_lme4)), length(fit_mm$fitted))
d_fitted <- compare("fitted max abs diff",
                    fitted(fit_lme4)[seq_len(n_common)],
                    fit_mm$fitted[seq_len(n_common)],
                    tols$fixef)

## ranef
re_lme4 <- lme4::ranef(fit_lme4)$g
if (!is.null(re_mm) && !is.null(re_lme4)) {
  lvl <- intersect(rownames(re_lme4), rownames(re_mm$g))
  if (length(lvl) > 0) {
    d_ranef <- compare("ranef g max abs diff",
                       as.numeric(re_lme4[lvl, 1]),
                       as.numeric(re_mm$g[lvl, 1]),
                       tols$fixef)
  } else {
    cat("ranef: no common level names found\n")
  }
}

## ── 7. Speed summary ──────────────────────────────────────────────────────────
cat("\n=== SPEED SUMMARY ===\n")
cat(sprintf("N=%d, ng=%d, formula=y~x+(1|g), family=binomial(logit), Laplace\n", N, ng))
cat(sprintf("lme4   median wall-clock : %.4f s\n", median(lme4_times)))
cat(sprintf("mixeff median wall-clock : %.4f s\n", median(mm_times)))
cat(sprintf("Speed ratio (mixeff/lme4): %.2fx\n", ratio_median))
if (ratio_median < 1) {
  cat(sprintf("  => mixeff is %.2fx FASTER than lme4\n", 1/ratio_median))
} else {
  cat(sprintf("  => mixeff is %.2fx SLOWER than lme4\n", ratio_median))
}

cat("\n=== OVERALL STATUS ===\n")
cat("Tolerances: fixef=1e-4, SE=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
all_diffs <- c(
  fixef  = d_fixef,
  SE     = d_se,
  theta  = d_theta,
  sigma  = d_sigma,
  logLik = d_loglik,
  deviance = d_deviance,
  AIC    = d_aic,
  BIC    = d_bic,
  fitted = d_fitted
)
max_overall <- max(all_diffs, na.rm = TRUE)
cat(sprintf("Max abs diff across all quantities: %.3e\n", max_overall))
cat("Outcome: ")
if (all(c(d_fixef, d_se) <= 1e-4) &&
    d_theta <= 1e-3 &&
    d_sigma <= 1e-4 &&
    d_loglik <= 1e-3) {
  cat("WITHIN-TOL (parity achieved)\n")
} else {
  cat("EXCEEDS-TOL* (divergence detected)\n")
}
