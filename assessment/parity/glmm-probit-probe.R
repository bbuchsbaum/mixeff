## glmm-probit parity probe
## Fits y ~ x + (1|g) with binomial(probit) in lme4 and mixeff,
## then compares fixef, SE, theta, logLik, AIC, BIC, ranef, fitted, convergence.
##
## NOTE: mixeff masks several lme4 generics (fixef, ranef, VarCorr, getME …).
## All lme4 accessor calls use the lme4:: prefix to avoid dispatch errors.

library(lme4)
library(mixeff)

set.seed(42)
N   <- 200
ng  <- 20
g   <- rep(seq_len(ng), each = N / ng)
x   <- rnorm(N)
b0  <- 0.3
b1  <- 0.8
u   <- rnorm(ng, sd = 0.6)
eta <- b0 + b1 * x + u[g]
y   <- rbinom(N, 1, pnorm(eta))

dat <- data.frame(y = y, x = x, g = factor(g))

# ── lme4 ─────────────────────────────────────────────────────────────────────
t_lme4 <- system.time({
  fit_lme4 <- lme4::glmer(y ~ x + (1|g), data = dat,
                          family = binomial(link = "probit"))
})

# ── mixeff ───────────────────────────────────────────────────────────────────
t_mixeff <- system.time({
  fit_mx <- tryCatch(
    mixeff::glmm(y ~ x + (1|g), data = dat,
                 family = binomial(link = "probit")),
    error = function(e) e
  )
})

cat("=== mixeff fit class ===\n")
print(class(fit_mx))

if (inherits(fit_mx, "condition")) {
  cat("=== mixeff ERROR ===\n")
  cat(conditionMessage(fit_mx), "\n")
  quit(save = "no")
}

# ── fixef ─────────────────────────────────────────────────────────────────────
fe_lme4 <- lme4::fixef(fit_lme4)
fe_mx   <- mixeff::fixef(fit_mx)

cat("\n=== fixef lme4 ===\n");  print(fe_lme4)
cat("=== fixef mixeff ===\n"); print(fe_mx)

fe_diff <- abs(fe_mx[names(fe_lme4)] - fe_lme4)
cat("=== fixef |diff| ===\n"); print(fe_diff)
cat("MAX fixef |diff|:", max(fe_diff), "\n")

# ── SE / vcov ────────────────────────────────────────────────────────────────
se_lme4 <- sqrt(diag(as.matrix(lme4::vcov.merMod(fit_lme4))))
se_mx   <- tryCatch(
  sqrt(diag(as.matrix(vcov(fit_mx)))),
  error = function(e) { cat("vcov(fit_mx) error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(se_mx)) {
  se_diff <- abs(se_mx - se_lme4)
  cat("\n=== SE lme4 ===\n");  print(se_lme4)
  cat("=== SE mixeff ===\n"); print(se_mx)
  cat("=== SE |diff| ===\n"); print(se_diff)
  cat("MAX SE |diff|:", max(se_diff), "\n")
}

# ── theta (random-effect Cholesky factor) ─────────────────────────────────────
th_lme4 <- lme4::getME(fit_lme4, "theta")
th_mx   <- fit_mx$theta

cat("\n=== theta lme4 ===\n");  print(th_lme4)
cat("=== theta mixeff ===\n"); print(th_mx)
th_diff <- abs(th_mx - th_lme4)
cat("=== theta |diff| ===\n"); print(th_diff)
cat("MAX theta |diff|:", max(th_diff), "\n")

# ── VarCorr ───────────────────────────────────────────────────────────────────
vc_lme4 <- lme4::VarCorr(fit_lme4)
cat("\n=== VarCorr lme4 ===\n"); print(vc_lme4)
cat("=== VarCorr mixeff ===\n")
vc_mx <- tryCatch(
  VarCorr(fit_mx),
  error = function(e) { cat("VarCorr(fit_mx) error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(vc_mx)) print(vc_mx)

# ── sigma (dispersion) ────────────────────────────────────────────────────────
sig_lme4 <- sigma(fit_lme4)
sig_mx   <- fit_mx$sigma

cat("\n=== sigma lme4:", sig_lme4, "===\n")
cat("=== sigma mixeff:", sig_mx, "===\n")
cat("sigma |diff|:", abs(sig_mx - sig_lme4), "\n")

# ── logLik ────────────────────────────────────────────────────────────────────
ll_lme4 <- as.numeric(logLik(fit_lme4))
ll_mx   <- fit_mx$logLik

cat("\n=== logLik lme4:", ll_lme4, "===\n")
cat("=== logLik mixeff:", ll_mx, "===\n")
cat("logLik |diff|:", abs(ll_mx - ll_lme4), "\n")

# ── AIC / BIC ─────────────────────────────────────────────────────────────────
aic_lme4 <- AIC(fit_lme4)
bic_lme4 <- BIC(fit_lme4)
aic_mx   <- fit_mx$AIC
bic_mx   <- fit_mx$BIC

cat("\n=== AIC lme4:", aic_lme4, " mixeff:", aic_mx, " |diff|:", abs(aic_mx - aic_lme4), "===\n")
cat("=== BIC lme4:", bic_lme4, " mixeff:", bic_mx, " |diff|:", abs(bic_mx - bic_lme4), "===\n")

# ── ranef ─────────────────────────────────────────────────────────────────────
re_lme4 <- lme4::ranef(fit_lme4)$g[, 1]
re_mx   <- tryCatch({
  r <- mixeff::ranef(fit_mx)
  if (is.list(r)) unlist(r[[1]], use.names = FALSE) else as.numeric(r)
}, error = function(e) {
  cat("ranef(fit_mx) error:", conditionMessage(e), "\n"); NULL
})

if (!is.null(re_mx) && length(re_mx) == length(re_lme4)) {
  re_diff <- abs(re_mx - re_lme4)
  cat("\n=== ranef (first 5) lme4 ===\n");  print(head(re_lme4))
  cat("=== ranef (first 5) mixeff ===\n"); print(head(re_mx))
  cat("MAX ranef |diff|:", max(re_diff), "\n")
} else {
  cat("\n=== ranef length mismatch or error ===\n")
  cat("lme4 len:", length(re_lme4), " mixeff len:", length(re_mx), "\n")
}

# ── fitted values ─────────────────────────────────────────────────────────────
fv_lme4 <- fitted(fit_lme4)
fv_mx   <- fit_mx$fitted

cat("\n=== fitted (first 5) lme4 ===\n");  print(head(fv_lme4))
cat("=== fitted (first 5) mixeff ===\n"); print(head(fv_mx))
if (!is.null(fv_mx) && length(fv_mx) == length(fv_lme4)) {
  fv_diff <- abs(fv_mx - fv_lme4)
  cat("MAX fitted |diff|:", max(fv_diff), "\n")
}

# ── convergence / status ──────────────────────────────────────────────────────
cat("\n=== lme4 convergence warnings ===\n")
msgs <- fit_lme4@optinfo$conv$lme4$messages
cat("n messages:", length(msgs), "\n")
if (length(msgs) > 0) cat(paste(msgs, collapse = "\n"), "\n")

cat("=== mixeff fit_status ===\n")
cat(fit_mx$fit_status, "\n")

# ── timing ────────────────────────────────────────────────────────────────────
cat("\n=== Timing ===\n")
cat("lme4   elapsed:", t_lme4["elapsed"], "s\n")
cat("mixeff elapsed:", t_mixeff["elapsed"], "s\n")
cat("speed ratio (lme4/mixeff):", t_lme4["elapsed"] / t_mixeff["elapsed"], "\n")
