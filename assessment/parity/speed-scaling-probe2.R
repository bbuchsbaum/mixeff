## speed-scaling-probe2.R — quiet version, only prints structured output
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

# Suppress all print output from fits
quiet <- function(expr) suppressMessages(suppressWarnings(capture.output(expr)))

set.seed(42)

sim_data <- function(N, n_groups) {
  g  <- factor(sample(seq_len(n_groups), N, replace = TRUE))
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  re_int <- rnorm(n_groups, 0, 1.2)[as.integer(g)]
  re_x1  <- rnorm(n_groups, 0, 0.4)[as.integer(g)]
  y <- 2 + 0.5*x1 + (-0.3)*x2 + re_int + re_x1*x1 + rnorm(N, 0, 0.8)
  data.frame(y=y, x1=x1, x2=x2, g=g)
}

compare <- function(label, v1, v2, tol) {
  d <- max(abs(as.numeric(v1) - as.numeric(v2)), na.rm=TRUE)
  status <- if (d <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("  %-36s diff=%.3e tol=%.0e [%s]\n", label, d, tol, status))
  invisible(d)
}

tols <- list(fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4)
NREPS <- 3

configs <- list(
  list(N=1000,  G=50,  label="N=1e3, G=50"),
  list(N=10000, G=100, label="N=1e4, G=100"),
  list(N=50000, G=200, label="N=5e4, G=200")
)

results <- list()

for (cfg in configs) {
  N <- cfg$N; G <- cfg$G; lbl <- cfg$label
  cat(sprintf("\n=== %s ===\n", lbl))
  df <- sim_data(N, G)

  ## single fit for parity
  fit_lme4 <- NULL; fit_mm <- NULL
  quiet(fit_lme4 <- lmer(y ~ x1 + x2 + (x1|g), data=df, REML=TRUE))
  fit_mm <- tryCatch(
    lmm(y ~ x1 + x2 + (x1|g), data=df, REML=TRUE, control=mm_control(verbose=-1L)),
    error=function(e) e
  )

  mm_err <- inherits(fit_mm, "condition")
  if (mm_err) {
    cat("  mixeff ERROR:", conditionMessage(fit_mm), "\n")
  } else {
    cat(sprintf("  fit_status: %s\n", fit_mm$fit_status))
    fe4  <- lme4::fixef(fit_lme4); fmm <- fixef(fit_mm)
    compare("fixef (Intercept)", fe4["(Intercept)"], fmm["(Intercept)"], tols$fixef)
    compare("fixef x1",          fe4["x1"],           fmm["x1"],          tols$fixef)
    compare("fixef x2",          fe4["x2"],           fmm["x2"],          tols$fixef)
    se4 <- sqrt(diag(vcov(fit_lme4))); semm <- fit_mm$std_errors
    compare("SE (Intercept)", se4["(Intercept)"], semm["(Intercept)"], tols$fixef)
    compare("SE x1",          se4["x1"],           semm["x1"],          tols$fixef)
    compare("SE x2",          se4["x2"],           semm["x2"],          tols$fixef)
    th4  <- lme4::getME(fit_lme4,"theta"); thmm <- fit_mm$theta
    cat(sprintf("  theta lme4: %s\n", paste(sprintf("%.6f",th4), collapse=", ")))
    cat(sprintf("  theta mm:   %s\n", paste(sprintf("%.6f",thmm), collapse=", ")))
    compare("theta",  th4,  thmm,  tols$theta)
    compare("sigma",  stats::sigma(fit_lme4), sigma(fit_mm), tols$sigma)
    compare("logLik", as.numeric(stats::logLik(fit_lme4)), as.numeric(logLik(fit_mm)), tols$logLik)
    compare("AIC",    stats::AIC(fit_lme4), fit_mm$AIC, tols$logLik*2)
    compare("BIC",    stats::BIC(fit_lme4), fit_mm$BIC, tols$logLik*2)
    compare("fitted", stats::fitted(fit_lme4), fit_mm$fitted[seq_len(N)], tols$fixef)
    re4  <- lme4::ranef(fit_lme4)$g
    remm <- ranef(fit_mm)$g
    if (!is.null(re4) && !is.null(remm)) {
      compare("ranef intercept", sort(re4[,"(Intercept)"]), sort(remm[,"(Intercept)"]), tols$fixef)
      compare("ranef x1",        sort(re4[,"x1"]),          sort(remm[,"x1"]),          tols$fixef)
    }
    max_d <- max(
      abs(as.numeric(fe4) - as.numeric(fmm)),
      abs(stats::sigma(fit_lme4) - sigma(fit_mm)),
      abs(as.numeric(stats::logLik(fit_lme4)) - as.numeric(logLik(fit_mm))),
      na.rm=TRUE
    )
    cat(sprintf("  MAX abs diff (fixef+sigma+logLik): %.3e\n", max_d))
  }

  ## timing
  t4 <- system.time(for(i in seq_len(NREPS)) {
    dfi <- sim_data(N, G)
    quiet(lmer(y ~ x1 + x2 + (x1|g), data=dfi, REML=TRUE))
  })
  tmm <- system.time(for(i in seq_len(NREPS)) {
    dfi <- sim_data(N, G)
    lmm(y ~ x1 + x2 + (x1|g), data=dfi, REML=TRUE, control=mm_control(verbose=-1L))
  })
  lme4_per <- t4["elapsed"] / NREPS
  mm_per   <- tmm["elapsed"] / NREPS
  ratio    <- mm_per / lme4_per
  verdict  <- if (ratio < 0.8) "FASTER" else if (ratio <= 1.5) "COMPARABLE" else "SLOWER"
  cat(sprintf("  lme4 mean/fit: %.4f s\n", lme4_per))
  cat(sprintf("  mm   mean/fit: %.4f s\n", mm_per))
  cat(sprintf("  ratio (mm/lme4): %.3fx  [%s]\n", ratio, verdict))
  results[[lbl]] <- list(lme4=lme4_per, mm=mm_per, ratio=ratio, verdict=verdict, mm_err=mm_err)
}

cat("\n=== SPEED SCALING SUMMARY ===\n")
cat(sprintf("%-20s  %8s  %8s  %8s  %10s\n","Config","lme4(s)","mm(s)","ratio","verdict"))
cat(strrep("-",58),"\n")
for (lbl in names(results)) {
  r <- results[[lbl]]
  if (r$mm_err) {
    cat(sprintf("%-20s  %8.4f  %8s  %8s  %10s\n", lbl, r$lme4, "ERROR","N/A","ERROR"))
  } else {
    cat(sprintf("%-20s  %8.4f  %8.4f  %8.3fx  %10s\n", lbl, r$lme4, r$mm, r$ratio, r$verdict))
  }
}
cat(strrep("-",58),"\n")
cat("Done.\n")
