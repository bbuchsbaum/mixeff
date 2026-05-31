.onLoad <- function(libname, pkgname) {
  mm_register_external_s3()
  invisible()
}

.onUnload <- function(libpath) {
  library.dynam.unload("mixeff", libpath)
}

# Tiny null-coalescing helper used by Phase 1+ verbs. Defined here so every
# R/ file in the package can rely on it without depending on rlang's
# `%||%` operator at top level (rlang is already in Imports for the
# typed-condition path; the bare operator is a one-liner not worth
# pulling).
`%||%` <- function(x, y) if (is.null(x)) y else x

mm_register_external_s3 <- function() {
  setHook(
    packageEvent("lme4", "onLoad"),
    function(...) mm_register_lme4_s3(),
    action = "append"
  )
  setHook(
    packageEvent("emmeans", "onLoad"),
    function(...) mm_register_emmeans_s3(),
    action = "append"
  )
  if ("lme4" %in% loadedNamespaces()) {
    mm_register_lme4_s3()
  }
  if ("emmeans" %in% loadedNamespaces()) {
    mm_register_emmeans_s3()
  }
}

mm_register_lme4_s3 <- function() {
  ns <- asNamespace("lme4")
  registerS3method("fixef", "mm_lmm", fixef.mm_lmm, envir = ns)
  registerS3method("fixef", "mm_glmm", fixef.mm_glmm, envir = ns)
  registerS3method("ranef", "mm_lmm", ranef.mm_lmm, envir = ns)
  registerS3method("ranef", "mm_glmm", ranef.mm_glmm, envir = ns)
  registerS3method("VarCorr", "mm_lmm", VarCorr.mm_lmm, envir = ns)
  registerS3method("VarCorr", "mm_glmm", VarCorr.mm_glmm, envir = ns)
  registerS3method("getME", "mm_lmm", getME.mm_lmm, envir = ns)
  registerS3method("refit", "mm_lmm", refit.mm_lmm, envir = ns)
  invisible(TRUE)
}

mm_register_emmeans_s3 <- function() {
  emmeans::.emm_register("mm_lmm", "mixeff")
  invisible(TRUE)
}
