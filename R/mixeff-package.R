#' mixeff: Audit-First Mixed-Effects Models via the 'mixedmodels' Rust Crate
#'
#' An R wrapper for the `mixedmodels` Rust crate. The package is audit-first:
#' every printed claim traces back to a versioned JSON artifact produced by
#' the Rust compiler, and the package refuses to fabricate inference results
#' the engine cannot certify. See `vignette("intro", package = "mixeff")` for
#' an overview and the demystification surface for random-effects syntax
#' (Phase 1+).
#'
#' @importFrom stats AIC BIC coef deviance df.residual fitted formula logLik model.frame nobs predict residuals setNames sigma update
#' @keywords internal
"_PACKAGE"
