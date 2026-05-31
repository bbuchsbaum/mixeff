## Empirical parity probe: lmm-pastes-nested
## Dataset: Pastes  Formula: strength ~ 1 + (1|batch/cask)  REML = TRUE
## Focus: nested grouping expansion  (batch/cask expands to batch + batch:cask)
## Compares lme4/lmerTest vs mixeff on all relevant quantities.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:", as.character(packageVersion("lme4")), "\n")
cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
cat("mixeff version:", as.character(packageVersion("mixeff")), "\n\n")

## ── 1. Data ──────────────────────────────────────────────────────────────────
data(Pastes, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(Pastes), "  ncol:", ncol(Pastes), "\n")
cat("columns:", paste(names(Pastes), collapse = ", "), "\n")
print(head(Pastes))
cat("Unique batch:", length(unique(Pastes$batch)), "\n")
cat("Unique cask:", length(unique(Pastes$cask)), "\n")
cat("\n")

## ── 2. Fit lme4 ──────────────────────────────────────────────────────────────
cat("=== lme4 FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(strength ~ 1 + (1 | batch/cask), data = Pastes, REML = TRUE)
})
cat("lme4 wall-clock (seconds):", t_lme4["elapsed"], "\n\n")

cat("-- fixef --\n");       print(lme4::fixef(fit_lme4))
cat("-- SE --\n");          print(sqrt(diag(vcov(fit_lme4))))
cat("-- vcov --\n");        print(as.matrix(vcov(fit_lme4)))
cat("-- VarCorr --\n");     print(lme4::VarCorr(fit_lme4))
cat("-- sigma --\n");       print(sigma(fit_lme4))
cat("-- logLik --\n");      print(logLik(fit_lme4))
cat("-- AIC --\n");         print(AIC(fit_lme4))
cat("-- BIC --\n");         print(BIC(fit_lme4))
cat("-- ranef names --\n"); print(names(lme4::ranef(fit_lme4)))
cat("-- ranef batch (head) --\n");      print(head(lme4::ranef(fit_lme4)$batch))
cat("-- ranef batch:cask (head) --\n"); print(head(lme4::ranef(fit_lme4)[["batch:cask"]]))
cat("-- fitted (first 6) --\n");        print(head(fitted(fit_lme4)))
cat("-- convergence --\n"); cat("No convergence warning =",
  length(fit_lme4@optinfo$conv$lme4$messages) == 0, "\n\n")

## Print what lme4 theta looks like for this nested model
theta_lme4 <- lme4::getME(fit_lme4, "theta")
cat("-- theta (length", length(theta_lme4), ") --\n"); print(theta_lme4)
cat("\n")

## ── 3. Fit mixeff ─────────────────────────────────────────────────────────────
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    lmm(strength ~ 1 + (1 | batch/cask), data = Pastes, REML = TRUE,
        control = mm_control(verbose = -1L)),
    error   = function(e) e,
    warning = function(w) w
  )
})
cat("mixeff wall-clock (seconds):", t_mixeff["elapsed"], "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff CONDITION !!!\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm), "\n")

  ## Try alternative: explicit expansion (1|batch) + (1|batch:cask)
  cat("\n--- Trying explicit expansion: (1|batch) + (1|batch:cask) ---\n")
  t_mixeff2 <- system.time({
    fit_mm2 <- tryCatch(
      lmm(strength ~ 1 + (1 | batch) + (1 | batch:cask), data = Pastes, REML = TRUE,
          control = mm_control(verbose = -1L)),
      error   = function(e) e,
      warning = function(w) w
    )
  })
  cat("mixeff (explicit) wall-clock (seconds):", t_mixeff2["elapsed"], "\n\n")
  if (inherits(fit_mm2, "condition")) {
    cat("!!! mixeff explicit expansion also failed !!!\n")
    cat("class:", paste(class(fit_mm2), collapse = ", "), "\n")
    cat("message:", conditionMessage(fit_mm2), "\n")
    fit_mm <- fit_mm2   # keep for summary
    quit(status = 0)
  } else {
    fit_mm <- fit_mm2
    t_mixeff <- t_mixeff2
    cat("Explicit expansion succeeded.\n\n")
  }
}

cat("fit_status:", fit_mm$fit_status, "\n\n")

cat("-- fixef --\n");  print(fixef(fit_mm))
cat("-- SE --\n");     print(fit_mm$std_errors)
cat("-- vcov --\n");   print(as.matrix(fit_mm$fixed_effect_vcov))
cat("-- VarCorr --\n"); print(VarCorr(fit_mm))
cat("-- sigma --\n");  print(sigma(fit_mm))
cat("-- logLik --\n"); print(logLik(fit_mm))
cat("-- AIC --\n");    print(AIC(fit_mm))
cat("-- BIC --\n");    print(BIC(fit_mm))
cat("-- ranef names --\n"); print(names(ranef(fit_mm)))
cat("-- ranef (all, head of each) --\n")
for (grp in names(ranef(fit_mm))) {
  cat("  group:", grp, "\n")
  print(head(ranef(fit_mm)[[grp]]))
}
cat("-- fitted (first 6) --\n"); print(head(fit_mm$fitted))
cat("-- theta --\n"); print(fit_mm$theta)
cat("\n")

## ── 4. Numerical comparison ───────────────────────────────────────────────────
cat("=== NUMERICAL COMPARISON ===\n")
tols <- list(fixef = 1e-4, theta = 1e-3, logLik = 1e-3, sigma = 1e-4)

compare_q <- function(label, lme4_val, mm_val, tol, fmt = "%.8f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  if (any(is.na(mm_val))) {
    cat(sprintf("%-35s lme4=%s  mm=NA  [MISSING]\n",
                label, paste(sprintf(fmt, lme4_val), collapse=", ")))
    return(invisible(NA_real_))
  }
  diff   <- max(abs(lme4_val - mm_val))
  status <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-35s lme4=%s  mm=%s  maxAbsDiff=%.3e  tol=%.0e  [%s]\n",
              label,
              paste(sprintf(fmt, lme4_val), collapse=", "),
              paste(sprintf(fmt, mm_val),   collapse=", "),
              diff, tol, status))
  invisible(diff)
}

## fixef
compare_q("fixef (Intercept)",
          lme4::fixef(fit_lme4)[["(Intercept)"]],
          fixef(fit_mm)[["(Intercept)"]],
          tols$fixef)

## SE
compare_q("SE (Intercept)",
          sqrt(diag(vcov(fit_lme4)))[["(Intercept)"]],
          fit_mm$std_errors[["(Intercept)"]],
          tols$fixef)

## vcov[1,1]
compare_q("vcov[1,1]",
          as.numeric(vcov(fit_lme4)[1,1]),
          as.numeric(fit_mm$fixed_effect_vcov[1,1]),
          tols$fixef^2)

## sigma
compare_q("sigma",
          sigma(fit_lme4),
          sigma(fit_mm),
          tols$sigma)

## logLik
compare_q("logLik",
          as.numeric(logLik(fit_lme4)),
          as.numeric(logLik(fit_mm)),
          tols$logLik)

## AIC
compare_q("AIC",
          AIC(fit_lme4),
          AIC(fit_mm),
          tols$logLik * 2)

## BIC
compare_q("BIC",
          BIC(fit_lme4),
          BIC(fit_mm),
          tols$logLik * 2)

## theta — nested model has 2 theta parameters: batch and batch:cask (or batch:cask only)
theta_mm <- fit_mm$theta
cat(sprintf("\n%-35s lme4=%s  mm=%s\n", "theta (all)",
            paste(sprintf("%.8f", theta_lme4), collapse=", "),
            paste(sprintf("%.8f", theta_mm),   collapse=", ")))
if (length(theta_lme4) == length(theta_mm)) {
  compare_q("theta (all)", theta_lme4, theta_mm, tols$theta)
} else {
  cat(sprintf("  theta length mismatch: lme4=%d, mm=%d — skipping element-wise\n",
              length(theta_lme4), length(theta_mm)))
}

## VarCorr: extract batch and batch:cask variances from lme4
vc_lme4    <- lme4::VarCorr(fit_lme4)
var_batch_lme4 <- as.numeric(vc_lme4$batch)[1]
var_casc_lme4  <- as.numeric(vc_lme4[["batch:cask"]])[1]
cat(sprintf("\n%-35s lme4_batch=%.8f  lme4_batch:cask=%.8f\n",
            "VarCorr variances", var_batch_lme4, var_casc_lme4))

## Extract from mixeff VarCorr
vc_mm <- VarCorr(fit_mm)
cat("mixeff VarCorr structure:\n"); print(vc_mm)

extract_vc_var <- function(vc, grpname) {
  tryCatch({
    if (is.data.frame(vc)) {
      rows <- vc[grepl(grpname, vc$grp, fixed=TRUE), ]
      if (nrow(rows) > 0) return(as.numeric(rows$vcov[1]))
    } else if (is.list(vc)) {
      v <- vc[[grpname]]
      if (is.null(v)) {
        # try partial match
        nm <- grep(grpname, names(vc), value=TRUE, fixed=TRUE)
        if (length(nm) > 0) v <- vc[[nm[1]]]
      }
      if (!is.null(v)) {
        if (is.numeric(v))   return(v[1])
        if (is.matrix(v))    return(v[1,1])
        if (is.data.frame(v)) return(as.numeric(v[1,1]))
      }
    }
    NA_real_
  }, error = function(e) NA_real_)
}

var_batch_mm <- extract_vc_var(vc_mm, "batch")
var_casc_mm  <- extract_vc_var(vc_mm, "batch:cask")
cat(sprintf("%-35s mm_batch=%s  mm_batch:cask=%s\n",
            "VarCorr variances (mm)",
            ifelse(is.na(var_batch_mm), "NA", sprintf("%.8f", var_batch_mm)),
            ifelse(is.na(var_casc_mm),  "NA", sprintf("%.8f", var_casc_mm))))

if (!is.na(var_batch_mm))
  compare_q("VarCorr batch var", var_batch_lme4, var_batch_mm, tols$theta)
if (!is.na(var_casc_mm))
  compare_q("VarCorr batch:cask var", var_casc_lme4, var_casc_mm, tols$theta)

## ranef comparison — batch
re_lme4_batch <- lme4::ranef(fit_lme4)$batch
re_mm_batch   <- ranef(fit_mm)$batch
if (!is.null(re_mm_batch) && !is.null(re_lme4_batch)) {
  # align by row names
  common <- intersect(rownames(re_lme4_batch), rownames(re_mm_batch))
  if (length(common) > 0) {
    compare_q("ranef batch (Intercept)",
              as.numeric(re_lme4_batch[common, 1]),
              as.numeric(re_mm_batch[common, 1]),
              tols$fixef)
  } else {
    # compare sorted
    compare_q("ranef batch sorted",
              sort(as.numeric(re_lme4_batch[,1])),
              sort(as.numeric(re_mm_batch[,1])),
              tols$fixef)
  }
} else {
  cat(sprintf("%-35s lme4 has ranef$batch; mm ranef$batch = %s\n",
              "ranef batch", ifelse(is.null(re_mm_batch), "NULL", "present")))
}

## ranef comparison — batch:cask
re_lme4_casc <- lme4::ranef(fit_lme4)[["batch:cask"]]
re_mm_casc   <- ranef(fit_mm)[["batch:cask"]]
if (!is.null(re_mm_casc) && !is.null(re_lme4_casc)) {
  common2 <- intersect(rownames(re_lme4_casc), rownames(re_mm_casc))
  if (length(common2) > 0) {
    compare_q("ranef batch:cask (Intercept)",
              as.numeric(re_lme4_casc[common2, 1]),
              as.numeric(re_mm_casc[common2, 1]),
              tols$fixef)
  } else {
    compare_q("ranef batch:cask sorted",
              sort(as.numeric(re_lme4_casc[,1])),
              sort(as.numeric(re_mm_casc[,1])),
              tols$fixef)
  }
} else {
  cat(sprintf("%-35s lme4 has ranef batch:cask; mm = %s\n",
              "ranef batch:cask", ifelse(is.null(re_mm_casc), "NULL", "present")))
}

## fitted values
fitted_lme4 <- fitted(fit_lme4)
fitted_mm   <- fit_mm$fitted
compare_q("fitted max abs diff",
          as.numeric(fitted_lme4),
          as.numeric(fitted_mm[seq_along(fitted_lme4)]),
          tols$fixef)

## speed ratio
cat(sprintf("\n%-35s lme4=%.4fs  mm=%.4fs  ratio(mm/lme4)=%.2fx\n",
            "wall-clock elapsed",
            t_lme4["elapsed"], t_mixeff["elapsed"],
            t_mixeff["elapsed"] / t_lme4["elapsed"]))

cat("\n=== SUMMARY TABLE ===\n")
cat("Tolerances: fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4\n")
cat("All WITHIN-TOL = parity achieved.\n")
