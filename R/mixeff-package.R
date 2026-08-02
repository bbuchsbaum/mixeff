#' mixeff: Mixed-Effects Models via the 'mixeff-rs' Rust Crate
#'
#' An R wrapper for the `mixeff-rs` Rust crate. `lmm()` and `glmm()` fit
#' linear and generalized linear mixed-effects models from lme4-style
#' formulas. Fitted models store the compiled model, optimizer result, and
#' inference metadata as data, so a result can be inspected, saved, and
#' reloaded; where an inference method is unavailable, the result is either
#' withheld with a stable reason code or labeled with the method actually
#' used, never silently swapped. See
#' `vignette("mixeff", package = "mixeff")` for an overview and the
#' supported model contract, and `vignette("lmm", package = "mixeff")`
#' for what random-effects formulas mean.
#'
#' @importFrom stats AIC BIC coef deviance df.residual fitted formula logLik model.frame nobs predict residuals setNames sigma update
#' @keywords internal
"_PACKAGE"
