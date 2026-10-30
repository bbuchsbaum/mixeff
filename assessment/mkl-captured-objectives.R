# Analyze original check-process fits only after the full check has finished.
check_library <- Sys.getenv("MIXEFF_DIAGNOSTIC_LIBRARY", "check/mixeff.Rcheck")
.libPaths(c(normalizePath(check_library, mustWork = TRUE), .libPaths()))
source("assessment/mkl-diagnostic.R")
matched <- 0L
capture_dir <- Sys.getenv("MIXEFF_DIAG_OUT", "/tmp/mkl-pairs")
for (path in list.files(capture_dir, pattern = "\\.rds$", full.names = TRUE)) {
  record <- readRDS(path)
  if (!isTRUE(record$REML) ||
      !identical(paste(deparse(record$formula), collapse = " "),
                 paste(deparse(f), collapse = " "))) next
  matched <- matched + 1L
  cat("\nORIGINAL CHECK FIT:", path, "\n")
  print(record[c("formula", "REML", "control", "warnings", "session", "software",
                 "locale", "environment", "options", "library_paths", "package_paths",
                 "reference_optinfo", "covariance", "logLik")])
  original_theta <- lapply(record$pair, theta_from_vc)
  cat("ORIGINAL COVARIANCES IN COMMON OBJECTIVE:\n")
  print(vapply(original_theta, devfun, numeric(1)))
  run(paste0("tight_original_", basename(path)),
      lme4::lmer(f, d, REML = TRUE, start = list(theta = original_theta$lme4),
                 control = tight))
}
stopifnot(matched > 0L)
saveRDS(records, file.path(outdir, "records-with-captured-fits.rds"))
