#' Fit a generalized linear mixed model
#'
#' `glmm()` validates the R-side family/link request, compiles the model
#' formula, and delegates the numerical fit to the upstream Rust
#' `GeneralizedLinearMixedModel`. The default `method = "pirls_profiled"` is
#' the labelled fast-PIRLS path. `method = "joint_laplace"` uses the upstream
#' labelled joint Laplace route (`fast = FALSE`, `nAGQ = 1`) backed by the
#' native dependency-light optimizer in this vendored build.
#'
#' @param formula A two-sided lme4-style formula.
#' @param data A `data.frame`.
#' @param family A supported GLMM family object or family constructor. The
#'   certified 1.0 surface is: [binomial()] with `"logit"`, `"probit"`, or
#'   `"cloglog"` links; [poisson()] with `"log"` or `"sqrt"` links;
#'   [Gamma()] with `"log"` link; and negative binomial (NB2, `"log"` link)
#'   via [mm_negative_binomial()] (theta estimated, like `lme4::glmer.nb()`)
#'   or `MASS::negative.binomial(theta)` (fixed theta).
#' @param random Reserved for the native random-effect constructor path.
#' @param weights Optional prior weights. For binomial models these are trial
#'   counts for proportion responses; weights must be positive and finite.
#' @param offset Optional fixed linear-predictor offset; values must be finite.
#' @param subset,na.action,contrasts Reserved for future parity with [lmm()].
#' @param method GLMM estimation method. `"pirls_profiled"` is the default
#'   fast-PIRLS profiled path. `"joint_laplace"` requests the labelled joint
#'   Laplace route and requires `nAGQ <= 1`. The joint route tracks the lme4
#'   joint-Laplace reference far more closely than the profiled path on
#'   high-baseline models, at a higher optimizer cost; cap that cost with
#'   `mm_control(max_feval = )`. The default profiled path is **not** glmer's
#'   estimator and its coefficients do not match `glmer()` exactly; when
#'   `method` is left at its default, `glmm()` emits an informational notice
#'   to that effect (suppress with `mm_control(verbose = -1)`). Use
#'   `method = "joint_laplace"` for glmer-equivalent estimates.
#' @param nAGQ Number of adaptive Gauss-Hermite quadrature points. `1` is the
#'   Laplace setting. Values above `1` are allowed on the profiled path and
#'   are rejected for `method = "joint_laplace"` in the R wrapper.
#' @param inference Requested inference method.
#' @param control A list from [mm_control()].
#' @param ... Reserved for future use.
#'
#' @return An object of class `mm_glmm`, also inheriting from `mm_fit` and
#'   `mm_compiled`.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rbinom(120, 1, 0.5),
#'   x = rnorm(120),
#'   g = factor(rep(seq_len(12), each = 10))
#' )
#' fit <- glmm(y ~ x + (1 | g), df, family = binomial(),
#'             control = mm_control(verbose = -1))
#' fixef(fit)
#' # glmer-equivalent (joint Laplace) estimates:
#' fit_joint <- glmm(y ~ x + (1 | g), df, family = binomial(),
#'                   method = "joint_laplace",
#'                   control = mm_control(verbose = -1))
#' fixef(fit_joint)
#'
#' @importFrom stats na.omit
#' @export
glmm <- function(formula,
                 data,
                 family,
                 random = NULL,
                 weights = NULL,
                 offset = NULL,
                 subset = NULL,
                 na.action = na.omit,
                 contrasts = NULL,
                 method = c("pirls_profiled", "joint_laplace"),
                 nAGQ = 1L,
                 inference = c("auto", "none", "asymptotic", "bootstrap"),
                 control = mm_control(),
                 ...) {
  call <- match.call()
  method_explicit <- !missing(method)
  method <- match.arg(method)
  inference <- match.arg(inference)
  control <- mm_validate_control(control)
  family_info <- mm_glmm_family_info(family)
  nAGQ <- mm_glmm_validate_nagq(nAGQ, method)

  if (identical(family_info$family, "negative_binomial") &&
      identical(method, "joint_laplace")) {
    # The engine's NB theta loop routes through fast-PIRLS internally, so the
    # joint route is not wired for NB at this pin; refusing here gives a
    # clear message instead of the engine's internal optimizer-guard error.
    mm_abort(
      message = paste0(
        "Negative-binomial GLMMs support the default profiled method only; ",
        "method = \"joint_laplace\" is not yet available for this family. ",
        "Drop the method argument to use the profiled fit."
      ),
      class = "mm_inference_unavailable",
      reason_code = "nb_joint_laplace_unavailable",
      input = method
    )
  }

  if (!is.null(random) || !is.null(subset) ||
      !identical(na.action, na.omit) || !is.null(contrasts)) {
    mm_abort(
      message = "`random`, `subset`, custom `na.action`, and `contrasts` are reserved for the fitted GLMM bridge.",
      class = "mm_fit_error",
      input = call
    )
  }
  weights <- mm_glmm_validate_weights(weights, data, "weights")
  offset <- mm_glmm_validate_weights(offset, data, "offset", positive = FALSE)

  # Resolve binomial responses: translate a cbind(successes, failures) LHS into
  # a proportion response + trial-count weights, and pick the engine family
  # ("binomial" for grouped/weighted data, "bernoulli" for 0/1 responses).
  prep <- mm_glmm_binomial_prep(formula, data, family_info, weights)
  formula <- prep$formula
  data <- prep$data
  weights <- prep$weights
  engine_family <- prep$engine_family
  if (identical(family_info$family, "negative_binomial")) {
    # Theta mode rides the family string (see mm_glmm_nb_engine_spec).
    engine_family <- mm_glmm_nb_engine_spec(family_info)
  }

  # lme4 parity: grouping variables must be categorical. Coerce non-factor /
  # non-character grouping columns to factors (announced, not silent) so an
  # integer subject/item ID does not hit the native "grouping factor not
  # categorical" refusal.
  data <- mm_apply_grouping_coercion(formula, data, control$verbose)

  spec <- compile_model(formula, data)
  mm_validate_fit_structure(spec, lmm = FALSE)
  if (control$verbose >= 0L) {
    mm_inform_explanation(spec)
    # No silent surgery on the estimator choice: when the user did not pick a
    # method, surface that the default profiled path is NOT glmer's estimator
    # and point to the certified glmer-equivalent route. Suppressed by
    # mm_control(verbose = -1) (as used in loops/bootstrap).
    if (!method_explicit && identical(method, "pirls_profiled")) {
      mm_inform(mm_glmm_profiled_default_notice(), class = "mm_estimator_notice")
    }
  }

  spec_data <- mm_translate_data(spec$model_frame)
  formula_string <- mm_coerce_formula_string(formula)
  control_json <- jsonlite::toJSON(unclass(control), auto_unbox = TRUE,
                                   null = "null")

  json <- tryCatch(
    .Call(
      wrap__mm_fit_glmm_json,
      formula_string,
      engine_family,
      family_info$link,
      method,
      nAGQ,
      spec_data$column_order,
      spec_data$numeric_columns,
      spec_data$categorical_values,
      spec_data$categorical_levels,
      spec_data$categorical_ordered,
      mm_bridge_weights(weights),
      mm_bridge_weights(offset),
      as.character(control_json)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(
      json,
      formula = formula_string,
      metadata = list(
        family = family_info,
        method = method,
        nAGQ = nAGQ,
        inference_request = inference
      ),
      spec = spec
    )
  }

  fit_result <- mm_json_parse_glmm_fit(json)
  fit_summary <- mm_json_parse_fit_summary(fit_result$fit_summary)
  artifact <- mm_json_parse_artifact(fit_result$artifact_json)
  beta <- mm_named_numeric(fit_result$beta, fit_result$beta_names)
  std_errors <- mm_named_numeric(fit_result$std_errors, fit_result$beta_names)
  fixed_effect_vcov <- mm_fixed_effect_vcov_from_payload(
    artifact$fixed_effect_covariance_matrix,
    beta,
    std_errors
  )

  if (identical(family_info$family, "negative_binomial")) {
    # Record the fitted (or fixed) NB2 size parameter on the family so
    # summaries can report it; nb_theta_estimated stays TRUE for estimated
    # fits so downstream engine refits reproduce the original estimation.
    family_info$nb_theta <- as.numeric(fit_result$nb_theta %||%
                                         family_info$nb_theta)
  }

  fit <- list(
    call           = call,
    formula        = formula,
    family         = family_info,
    method         = as.character(fit_result$method %||% method),
    nAGQ           = as.integer(fit_result$n_agq %||% nAGQ),
    inference_request = inference,
    control        = control,
    vars           = spec$vars,
    model_frame    = spec$model_frame,
    weights        = weights,
    offset         = offset,
    artifact       = artifact,
    fit            = fit_result,
    fit_summary    = fit_summary,
    schema         = mm_object_schema(artifact),
    rust_handle    = NULL,
    lazy_cache     = mm_empty_lazy_cache(),
    beta           = beta,
    theta          = as.numeric(unlist(fit_result$theta, use.names = FALSE)),
    sigma          = as.numeric(fit_result$dispersion),
    dispersion     = as.numeric(fit_result$dispersion),
    logLik         = as.numeric(fit_result$log_likelihood),
    deviance       = as.numeric(fit_result$deviance),
    AIC            = as.numeric(fit_result$aic),
    BIC            = as.numeric(fit_result$bic),
    nobs           = as.integer(fit_result$nobs),
    dof            = as.integer(fit_result$dof),
    df_residual    = as.integer(fit_result$df_residual),
    fit_status     = as.character(fit_result$fit_status),
    std_errors     = std_errors,
    fixed_effect_vcov = fixed_effect_vcov,
    fixed_fitted   = NULL,
    fitted         = as.numeric(unlist(fit_result$fitted, use.names = FALSE)),
    residuals      = as.numeric(unlist(fit_result$residuals, use.names = FALSE)),
    random_effects = mm_ranef_from_terms(fit_result$ranef),
    varcorr        = mm_varcorr_from_result(
      fit_summary$varcorr %||% fit_result$varcorr,
      artifact = artifact
    )
  )
  fit <- mm_apply_lme4_coef_naming(fit)
  class(fit) <- c("mm_glmm", "mm_fit", "mm_compiled")
  fit
}

mm_glmm_profiled_default_notice <- function() {
  paste0(
    "glmm() is using the default method = \"pirls_profiled\": a fast profiled ",
    "(PIRLS) approximation whose coefficients do NOT exactly match ",
    "lme4::glmer(). For glmer-equivalent (joint Laplace) estimates, use ",
    "method = \"joint_laplace\". Silence this message with ",
    "mm_control(verbose = -1)."
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
  family_name <- fam$family
  link_name <- fam$link
  nb_theta <- NULL
  nb_theta_estimated <- FALSE
  # MASS::negative.binomial(theta) encodes theta in the family name
  # ("Negative Binomial(2.5)"): fixed-theta NB2, fit conditional on theta.
  # mm_negative_binomial() yields family "negative_binomial" with an optional
  # $theta; theta = NULL requests glmer.nb-style estimation.
  if (grepl("^Negative Binomial\\(", family_name)) {
    nb_theta <- suppressWarnings(
      as.numeric(sub("^Negative Binomial\\(([^)]*)\\)$", "\\1", family_name))
    )
    if (!length(nb_theta) || is.na(nb_theta) || !is.finite(nb_theta) ||
        nb_theta <= 0) {
      mm_abort(
        message = sprintf(
          "Could not read a positive finite theta from family `%s`.",
          family_name
        ),
        class = "mm_arg_error",
        input = family
      )
    }
    family_name <- "negative_binomial"
  } else if (identical(family_name, "negative_binomial")) {
    nb_theta <- fam$theta
    nb_theta_estimated <- is.null(nb_theta)
  }
  supported <- mm_glmm_supported_family_links()
  if (!family_name %in% names(supported) ||
      !link_name %in% supported[[family_name]]) {
    mm_abort_glmm_unsupported_family_link(family_name, link_name)
  }
  list(
    family = if (identical(family_name, "Gamma")) "gamma" else family_name,
    link = link_name,
    nb_theta = nb_theta,
    nb_theta_estimated = nb_theta_estimated
  )
}

#' Negative-binomial family for `glmm()`
#'
#' NB2 family (variance `mu + mu^2/theta`, log link). With `theta = NULL`
#' (the default) the size parameter is estimated alongside the model, matching
#' `lme4::glmer.nb()`. Supplying a positive `theta` fits conditional on that
#' value, matching `glmer(family = MASS::negative.binomial(theta))` — which
#' `glmm()` also accepts directly.
#'
#' @param theta Optional positive NB2 size (dispersion) parameter. `NULL`
#'   estimates theta from the data.
#' @return A `family` object accepted by [glmm()].
#' @examples
#' fam <- mm_negative_binomial()      # glmer.nb-style: theta estimated
#' fam_fixed <- mm_negative_binomial(theta = 2.5)
#' @export
mm_negative_binomial <- function(theta = NULL) {
  if (!is.null(theta)) {
    if (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta) ||
        theta <= 0) {
      mm_abort(
        message = "`theta` must be a single positive finite number (or NULL to estimate it).",
        class = "mm_arg_error",
        input = theta
      )
    }
    theta <- as.numeric(theta)
  }
  structure(
    list(family = "negative_binomial", link = "log", theta = theta),
    class = "family"
  )
}

# Engine wire spec for a family: negative-binomial rides its theta mode in
# the family string so no FFI signatures change ("negative_binomial" =
# estimate theta; "negative_binomial:theta=<v>" = fixed theta).
mm_glmm_nb_engine_spec <- function(family_info) {
  if (is.null(family_info$nb_theta) || isTRUE(family_info$nb_theta_estimated)) {
    "negative_binomial"
  } else {
    sprintf("negative_binomial:theta=%.17g", family_info$nb_theta)
  }
}

mm_glmm_supported_family_links <- function() {
  list(
    binomial = c("logit", "probit", "cloglog"),
    poisson = c("log", "sqrt"),
    Gamma = c("log"),
    negative_binomial = "log"
  )
}

mm_glmm_supported_family_link_table <- function() {
  supported <- mm_glmm_supported_family_links()
  do.call(rbind, lapply(names(supported), function(family) {
    data.frame(
      family = family,
      link = supported[[family]],
      stringsAsFactors = FALSE
    )
  }))
}

mm_abort_glmm_unsupported_family_link <- function(family, link) {
  reason_code <- "unsupported_glmm_family_link"
  mm_abort(
    message = sprintf(
      "GLMM family/link `%s/%s` is outside the certified upstream contract.",
      family, link
    ),
    class = "mm_inference_unavailable",
    reason_code = reason_code,
    family = family,
    link = link,
    supported = mm_glmm_supported_family_link_table(),
    input = list(family = family, link = link)
  )
}

# Resolve a binomial GLMM response. Handles two grouped-binomial spellings:
#   * cbind(successes, failures) ~ ...  -> proportion response + trial weights
#   * proportion ~ ..., weights = trials
# and chooses the engine family: "binomial" (grouped/weighted) vs "bernoulli"
# (0/1 responses). Non-binomial families pass through untouched.
mm_glmm_binomial_prep <- function(formula, data, family_info, weights) {
  engine_family <- family_info$family
  if (!identical(family_info$family, "binomial")) {
    return(list(formula = formula, data = data, weights = weights,
                engine_family = engine_family))
  }
  lhs <- formula[[2L]]
  is_cbind <- is.call(lhs) && identical(as.character(lhs[[1L]]), "cbind")
  if (is_cbind) {
    if (!is.null(weights)) {
      mm_abort(
        message = "Supply either a `cbind(successes, failures)` response or `weights =`, not both.",
        class = "mm_arg_error"
      )
    }
    env <- environment(formula) %||% parent.frame()
    succ <- as.numeric(eval(lhs[[2L]], data, env))
    fail <- as.numeric(eval(lhs[[3L]], data, env))
    if (length(succ) != nrow(data) || length(fail) != nrow(data)) {
      mm_abort(
        message = "`cbind()` response columns must each have one value per row of `data`.",
        class = "mm_data_error"
      )
    }
    n <- succ + fail
    if (any(!is.finite(n)) || any(n <= 0)) {
      mm_abort(
        message = "`cbind(successes, failures)` trial totals must be finite and positive.",
        class = "mm_data_error"
      )
    }
    respname <- ".mm_binomial_response"
    data[[respname]] <- succ / n
    weights <- as.numeric(n)
    formula[[2L]] <- as.name(respname)
  } else {
    response_name <- as.character(lhs)
    col <- data[[response_name]]
    if (is.factor(col)) {
      lvls <- levels(col)
      if (length(lvls) != 2L) {
        mm_abort(
          message = sprintf(
            "Binomial response '%s' is a factor with %d level(s); exactly 2 are required. Levels found: %s",
            response_name, length(lvls), paste(lvls, collapse = ", ")
          ),
          class = "mm_data_error"
        )
      }
      mm_inform(
        sprintf(
          "Coercing factor response '%s' to 0/1: treating '%s' as 1 (success).",
          response_name, lvls[[2L]]
        ),
        class = "mm_factor_coercion"
      )
      data[[response_name]] <- as.integer(col) - 1L
    } else if (is.logical(col)) {
      data[[response_name]] <- as.integer(col)
    }
  }
  engine_family <- if (is.null(weights)) "bernoulli" else "binomial"
  list(formula = formula, data = data, weights = weights,
       engine_family = engine_family)
}

# Validate a prior-weights or offset vector for glmm(). Returns NULL when not
# supplied, otherwise a numeric vector of length nrow(data). Weights must be
# positive; offsets need only be finite.
mm_glmm_validate_weights <- function(x, data, label, positive = TRUE) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != nrow(data) || anyNA(x) ||
      any(!is.finite(x)) || (positive && any(x <= 0))) {
    mm_abort(
      message = sprintf(
        "`%s` must be a %s numeric vector with one value per row of `data` (%d).",
        label, if (positive) "positive, finite" else "finite", nrow(data)
      ),
      class = "mm_arg_error",
      input = x
    )
  }
  as.numeric(x)
}

mm_glmm_validate_nagq <- function(nAGQ, method) {
  if (!is.numeric(nAGQ) || length(nAGQ) != 1L || is.na(nAGQ) || nAGQ < 1) {
    mm_abort(
      message = "`nAGQ` must be a single positive integer.",
      class = "mm_arg_error",
      input = nAGQ
    )
  }
  nAGQ <- as.integer(nAGQ)
  if (identical(method, "joint_laplace") && nAGQ > 1L) {
    mm_abort(
      message = "`method = \"joint_laplace\"` requires `nAGQ <= 1`.",
      class = "mm_arg_error",
      input = nAGQ
    )
  }
  nAGQ
}

mm_json_parse_glmm_fit <- function(json) {
  if (!is.character(json) || length(json) != 1L || is.na(json) || !nzchar(json)) {
    mm_abort(
      message = "`json` must be a single non-empty character string.",
      class = "mm_schema_error",
      input = json
    )
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse GLMM fit JSON: %s", conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  schema <- parsed$schema
  if (!is.list(schema) ||
      !identical(as.character(schema$schema_name), "mixeff.glmm_fit_result") ||
      !identical(as.character(schema$schema_version), "1")) {
    mm_abort(
      message = "GLMM fit JSON has an unknown schema header.",
      class = "mm_schema_error",
      input = parsed
    )
  }
  parsed
}
