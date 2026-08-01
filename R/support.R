#' The supported model envelope
#'
#' `supported_models()` returns the package's machine-readable support
#' registry: which GLMM family/link/estimator combinations are supported,
#' which features are supported, refused, or carry mixeff-specific
#' semantics, and the package's explicit non-goals. The registry ships as
#' `inst/support/model-support.json` and is the single source of truth for
#' the migration vignette's support table; the test suite asserts that the
#' fitting code enforces exactly this envelope (unsupported combinations
#' refuse with the declared condition classes and reason codes).
#'
#' mixeff aims for tested statistical agreement with lme4 on this envelope
#' — a documented, closed surface — not blanket equivalence with everything
#' `lmer()`/`glmer()` accept.
#'
#' @return A list of class `mm_model_support` with elements
#'   `glmm_families` (data frame: family, links, estimators, notes),
#'   `features` (data frame: feature, lmm, glmm, notes), `non_goals`
#'   (character), and `registry` (the parsed JSON, verbatim).
#'
#' @examples
#' supported_models()
#'
#' @export
supported_models <- function() {
  registry <- mm_model_support_registry()
  fam <- registry$glmm_families
  fam_df <- data.frame(
    family = vapply(fam, function(f) f$family, character(1)),
    links = vapply(fam, function(f)
      paste(unlist(f$links), collapse = ", "), character(1)),
    estimators = vapply(fam, function(f)
      paste(unlist(f$estimators), collapse = ", "), character(1)),
    notes = vapply(fam, function(f)
      as.character(f$notes %||% ""), character(1)),
    stringsAsFactors = FALSE
  )
  feats <- registry$features
  feat_df <- data.frame(
    feature = vapply(feats, function(f) f$feature, character(1)),
    lmm = vapply(feats, function(f) as.character(f$lmm), character(1)),
    glmm = vapply(feats, function(f) as.character(f$glmm), character(1)),
    notes = vapply(feats, function(f)
      as.character(f$notes %||% ""), character(1)),
    stringsAsFactors = FALSE
  )
  out <- list(
    glmm_families = fam_df,
    features = feat_df,
    non_goals = as.character(unlist(registry$non_goals)),
    registry = registry
  )
  class(out) <- "mm_model_support"
  out
}

mm_model_support_registry <- function() {
  path <- system.file("support", "model-support.json", package = "mixeff")
  if (!nzchar(path)) {
    # Source-tree fallback (devtools::load_all before install).
    path <- file.path("inst", "support", "model-support.json")
  }
  if (!file.exists(path)) {
    mm_abort(
      message = "The model-support registry (inst/support/model-support.json) is missing from this installation.",
      class = "mm_schema_error",
      input = path
    )
  }
  jsonlite::read_json(path)
}

#' @method print mm_model_support
#' @export
print.mm_model_support <- function(x, ...) {
  cat("mixeff supported model envelope",
      sprintf("(registry %s v%s)\n\n",
              x$registry$schema_name, x$registry$schema_version))
  cat("GLMM families:\n")
  print(x$glmm_families, row.names = FALSE, right = FALSE)
  cat("\nFeatures:\n")
  print(x$features[, c("feature", "lmm", "glmm")], row.names = FALSE,
        right = FALSE)
  cat("\nNon-goals:", paste(x$non_goals, collapse = "; "), "\n")
  cat("\nGaussian LMMs support the common random-intercept, random-slope,\n")
  cat("nested, and crossed structures; see vignette(\"lme4-migration\") for\n")
  cat("the argument-by-argument map and the feature notes column for\n")
  cat("mixeff-specific semantics.\n")
  invisible(x)
}
