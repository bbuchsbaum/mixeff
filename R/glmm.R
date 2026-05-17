#' Fit a generalized linear mixed model
#'
#' `glmm()` is the Phase 4 GLMM boundary. It validates the R-side family/link
#' request, compiles the model formula, and records the requested objective
#' approximation metadata. The current Rust bridge does not yet expose a GLMM
#' fit primitive, so the function raises a typed `mm_fit_error` instead of
#' pretending to fit.
#'
#' @param formula A two-sided lme4-style formula.
#' @param data A `data.frame`.
#' @param family A family object or family constructor such as [binomial()] or
#'   [poisson()].
#' @param random Reserved for the native random-effect constructor path.
#' @param weights,subset,na.action,contrasts Reserved for future parity with
#'   [lmm()].
#' @param approximation Objective approximation: `"laplace"` or `"agq"`.
#' @param nAGQ Number of adaptive Gauss-Hermite quadrature points. `1` is
#'   Laplace; values above `1` require `approximation = "agq"`.
#' @param inference Requested inference method.
#' @param control A list from [mm_control()].
#' @param ... Reserved for future use.
#'
#' @return This function currently raises `mm_fit_error`.
#'
#' @importFrom stats na.omit
#' @export
glmm <- function(formula,
                 data,
                 family,
                 random = NULL,
                 weights = NULL,
                 subset = NULL,
                 na.action = na.omit,
                 contrasts = NULL,
                 approximation = c("laplace", "agq"),
                 nAGQ = 1L,
                 inference = c("auto", "none", "asymptotic", "bootstrap"),
                 control = mm_control(),
                 ...) {
  call <- match.call()
  approximation <- match.arg(approximation)
  inference <- match.arg(inference)
  control <- mm_validate_control(control)
  family_info <- mm_glmm_family_info(family)
  nAGQ <- mm_glmm_validate_nagq(nAGQ, approximation)

  if (!is.null(random) || !is.null(weights) || !is.null(subset) ||
      !identical(na.action, na.omit) || !is.null(contrasts)) {
    mm_abort(
      message = "`random`, `weights`, `subset`, custom `na.action`, and `contrasts` are reserved for the fitted GLMM bridge.",
      class = "mm_fit_error",
      input = call
    )
  }

  spec <- compile_model(formula, data)
  if (control$verbose >= 0L) {
    print(explain_model(spec))
  }

  metadata <- list(
    family = family_info,
    approximation = list(
      method = approximation,
      nAGQ = nAGQ,
      objective = if (nAGQ <= 1L) "laplace" else "adaptive_gauss_hermite"
    ),
    inference_request = inference
  )
  mm_abort(
    message = "GLMM fitting is not available in this mixeff build because the Rust bridge does not expose a GLMM fit primitive yet.",
    class = "mm_fit_error",
    formula = deparse1(formula),
    metadata = metadata,
    spec = spec
  )
}

mm_glmm_family_info <- function(family) {
  if (missing(family)) {
    mm_abort(
      message = "`family` is required for `glmm()`.",
      class = "mm_arg_error"
    )
  }
  fam <- if (is.function(family)) family() else family
  if (!inherits(fam, "family")) {
    mm_abort(
      message = "`family` must be an R family object or family constructor.",
      class = "mm_arg_error",
      input = family
    )
  }
  supported <- list(
    binomial = c("logit", "probit", "cloglog"),
    poisson = c("log", "sqrt")
  )
  family_name <- fam$family
  link_name <- fam$link
  if (!family_name %in% names(supported) ||
      !link_name %in% supported[[family_name]]) {
    mm_abort(
      message = sprintf("GLMM family/link `%s/%s` is not supported.",
                        family_name, link_name),
      class = "mm_fit_error",
      input = list(family = family_name, link = link_name)
    )
  }
  list(
    family = family_name,
    link = link_name
  )
}

mm_glmm_validate_nagq <- function(nAGQ, approximation) {
  if (!is.numeric(nAGQ) || length(nAGQ) != 1L || is.na(nAGQ) || nAGQ < 1) {
    mm_abort(
      message = "`nAGQ` must be a single positive integer.",
      class = "mm_arg_error",
      input = nAGQ
    )
  }
  nAGQ <- as.integer(nAGQ)
  if (identical(approximation, "laplace") && nAGQ != 1L) {
    mm_abort(
      message = "`approximation = \"laplace\"` requires `nAGQ = 1`.",
      class = "mm_arg_error",
      input = nAGQ
    )
  }
  if (identical(approximation, "agq") && nAGQ <= 1L) {
    mm_abort(
      message = "`approximation = \"agq\"` requires `nAGQ > 1`.",
      class = "mm_arg_error",
      input = nAGQ
    )
  }
  nAGQ
}
