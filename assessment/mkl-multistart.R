# Bounded diagnostic; no production estimator or parity tolerance changes.
# Run after installing mixeff: Rscript assessment/mkl-multistart.R . OUTPUT
# The SD8.2 start is a constructed stress case, not the lost MKL fit.
.libPaths(c(Sys.getenv("MIXEFF_DIAGNOSTIC_LIBRARY", "check/mixeff.Rcheck"), .libPaths()))
diagnostic <- "assessment/mkl-diagnostic.R"
if (!file.exists(diagnostic)) {
  diagnostic <- "release-evidence/cancellation-2026-10-26/mkl-diagnostic.R"
}
source(diagnostic)
tight <- lme4::lmerControl(optimizer = "bobyqa", calc.derivs = TRUE,
  optCtrl = list(rhobeg = 0.01, rhoend = 1e-9, maxfun = 20000))
# lme4 2.1 also gates checkConv itself by sample/parameter count.
for (field in c("check.conv.nobsmax", "check.conv.nparmax")) {
  if (!is.null(tight$checkConv[[field]])) tight$checkConv[[field]] <- Inf
}
# Independent lme4 initial structure: unit diagonal and zero correlations.
# Nominate8.2 using the residual scale of an ordinary fixed-effects model,
# independent of either mixed-model fit. Sigma is profiled during fitting.
nominal_sigma <- sigma(lm(RT ~ modality, data = d))
stress_theta <- unlist(lapply(names(cnms), function(g) {
  if (g == "stim") c(1, 0, 8.2 / nominal_sigma)
  else if (g == "PID") c(1, 0, 1)
  else stop("Unexpected random-effect group")
}), use.names = FALSE)
names(stress_theta) <- names(theta_ref)
cat("\nPREDETERMINED STRESS START (relative Cholesky):\n")
print(stress_theta)
cat("Fixed-effects residual SD:", nominal_sigma,
    "; nominal stim slope SD8.2, not reconstruction of original failed fit.\n")
starts <- list(default_endpoint = theta_ref, stress_sd8_2 = stress_theta,
               engine = theta_mm)
for (label in names(starts)) {
  a <- list(formula = f, data = d, REML = TRUE, control = tight)
  if (!is.null(starts[[label]])) a$start <- list(theta = starts[[label]])
  run(paste0("multistart_", label), do.call(lme4::lmer, a))
}
selected <- fits[paste0("multistart_", names(starts))]
stopifnot(all(vapply(selected, inherits, logical(1), "merMod")))
result <- do.call(rbind, lapply(names(selected), function(label) {
  fit <- selected[[label]]; vc <- as.data.frame(lme4::VarCorr(fit))
  idx <- which(vc$grp == "stim" & vc$var1 == "modality" & is.na(vc$var2))
  stopifnot(length(idx) == 1L)
  data.frame(start = label, stim_slope_sd = vc$sdcor[idx],
    common_objective = devfun(theta_from_vc(fit)),
    optimizer_code = paste(fit@optinfo$conv$opt, collapse = ","),
    convergence_messages = paste(fit@optinfo$conv$lme4$messages, collapse = "; "),
    warnings = paste(records[[label]]$warnings, collapse = "; "))
}))
result$objective_excess <- result$common_objective - min(result$common_objective)
cat("\nBOUNDED MULTISTART SUMMARY:\n"); print(result, row.names = FALSE)
write.csv(result, file.path(outdir, "multistart-summary.csv"), row.names = FALSE)
saveRDS(list(starts = starts, nominal_sigma = nominal_sigma, control = tight,
  fits = selected, records = records, summary = result, session = sessionInfo(),
  fixture_md5 = tools::md5sum(fixture)), file.path(outdir, "multistart-all.rds"))
cat("Engine common objective excess above best bounded reference:",
    devfun(theta_mm) - min(result$common_objective), "\n")
cat("No assertion tolerance was changed. This is not a fixed-SD likelihood profile.\n")
