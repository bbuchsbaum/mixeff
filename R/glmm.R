#' Fit a generalized linear mixed model
#'
#' `glmm()` validates the R-side family/link request, compiles the model
#' formula, and delegates the numerical fit to the upstream Rust
#' `GeneralizedLinearMixedModel`. The default `method = "joint_laplace"` is
#' the labeled joint Laplace route (`fast = FALSE`, `nAGQ = 1`), the estimator
#' `lme4::glmer()` uses by default, backed by the native dependency-light
#' optimizer in this vendored build. `method = "pirls_profiled"` is the
#' labeled fast-PIRLS profiled path (lme4's `nAGQ = 0`).
#'
#' @param formula A two-sided lme4-style formula.
#' @param data A `data.frame`.
#' @param family A supported GLMM family object or family constructor. The
#'   supported surface (every family/link pair the engine fits) is:
#'   [binomial()] with `"logit"`, `"probit"`, or `"cloglog"` links;
#'   [poisson()] with `"log"` or `"sqrt"` links; [Gamma()] with `"inverse"`
#'   (R's default) or `"log"` links; [inverse.gaussian()] with `"inverse"` or
#'   `"log"` links (R's default `"1/mu^2"` link is not available and is
#'   refused); [gaussian()] with non-identity `"log"`, `"inverse"`, or
#'   `"sqrt"` links (a Gaussian identity-link model is an LMM: use [lmm()]);
#'   and negative binomial (NB2, `"log"` link) via [mm_negative_binomial()]
#'   (theta estimated, like `lme4::glmer.nb()`) or
#'   `MASS::negative.binomial(theta)` (fixed theta). Binomial `"cauchit"`,
#'   `"log"`, and `"identity"` and Poisson `"identity"` links are not
#'   available in the engine and are refused with a typed
#'   `mm_inference_unavailable` error naming the supported set.
#' @param random Reserved for the native random-effect constructor path.
#' @param weights Optional prior weights: a column of `data` (as in lme4,
#'   `weights = trials`) or a numeric vector. For binomial models these are
#'   trial counts for proportion responses; weights must be positive and
#'   finite.
#' @param offset Optional fixed linear-predictor offset: a column of `data` or
#'   a numeric vector; values must be finite. `offset()` terms in the formula
#'   are also supported and are added to it (as in [stats::glm()]).
#' @param subset,na.action,contrasts As in [lmm()]: `subset` selects rows
#'   (evaluated in `data`), `na.action` controls missing values (the default
#'   `NULL` refuses `NA` in a model variable with a typed `mm_data_error`;
#'   pass `na.omit` for lme4's complete-case behaviour, or `na.exclude` to pad
#'   `fitted()`/`residuals()`/`predict()` back to the original rows), and
#'   `contrasts` is honoured only when it names the engine's coding
#'   (`contr.treatment` for unordered, `contr.poly` for ordered factors).
#' @param method GLMM estimation method. `"joint_laplace"` (the default)
#'   is the joint Laplace route, matching `glmer()`'s default estimator, and
#'   requires `nAGQ = 1`; it certifies Wald standard errors, tests, and
#'   intervals. Its optimizer cost is higher than the profiled path's; cap it
#'   with `mm_control(max_feval = )`. `"pirls_profiled"` is the fast profiled
#'   PIRLS path (equivalent to lme4's `nAGQ = 0` fast estimate); its
#'   coefficients do not match `glmer(nAGQ = 1)` exactly and Wald inference
#'   is withheld. When `method` is not supplied and the request needs the
#'   profiled path -- a negative-binomial family (the joint route is not
#'   available for it yet), `nAGQ > 1`, or `inference = "working_hessian"` --
#'   `glmm()` uses `"pirls_profiled"`
#'   and says so with an `mm_estimator_notice` message (silence with
#'   `mm_control(verbose = -1)`). An explicit `method = "joint_laplace"` with a
#'   negative-binomial family or `nAGQ != 1` is refused, never swapped.
#' @param nAGQ Number of adaptive Gauss-Hermite quadrature points. `1` (the
#'   default) is the Laplace approximation. `0` requests lme4's PIRLS-only
#'   fast estimate and selects `method = "pirls_profiled"` (refused together
#'   with an explicit `method = "joint_laplace"`). Values above `1` run on the
#'   profiled path only and, as in lme4, require a model with a single scalar
#'   random-effect term (e.g. `(1 | g)`); other models are refused.
#' @param inference Requested inference posture. The default `"auto"` keeps
#'   the certified contract: Wald standard errors, tests, and intervals are
#'   available only when the engine certifies them (currently
#'   `method = "joint_laplace"`, the default); the profiled estimator
#'   withholds them with a typed refusal. `"working_hessian"` is an explicit opt-in
#'   that unlocks the UNCERTIFIED profiled working-Hessian approximation on
#'   every inference route; each resulting row is labelled
#'   `wald_z_working_hessian` with reliability `moderate`. Its standard
#'   errors ran about 11% smaller than `glmer()`'s on the package's
#'   reference dataset (anti-conservative), so treat it as an exploration
#'   and screening tool, not a reporting route; see `inference_options()`.
#'   `"none"`, `"asymptotic"`, and `"bootstrap"` are accepted and recorded
#'   but currently equivalent to `"auto"`.
#' @param control A list from [mm_control()].
#' @param ... Reserved for future use.
#'
#'
#' @details
#' Optimization runs inside a single native call with no progress output: the
#' pre-fit explanation block (when `verbose >= 0`) is the last thing printed
#' before the fitted result returns. The fit checks for a user interrupt
#' (Ctrl-C / Esc) between optimizer evaluations and stops with a typed
#' `mm_interrupted` error (the check happens at evaluation boundaries, so a
#' single very expensive evaluation finishes first). Evaluation budgets are bounded (a bounded budget caps optimizer
#' iterations; it does not prove every native evaluation terminates); runtime on
#' large problems is governed by `mm_control(max_feval = )`.
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
#' # Default: glmer-equivalent joint Laplace estimates.
#' fit <- glmm(y ~ x + (1 | g), df, family = binomial(),
#'             control = mm_control(verbose = -1))
#' fixef(fit)
#' # Fast profiled PIRLS path (lme4's nAGQ = 0):
#' fit_fast <- glmm(y ~ x + (1 | g), df, family = binomial(),
#'                  method = "pirls_profiled",
#'                  control = mm_control(verbose = -1))
#' fixef(fit_fast)
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
                 na.action = NULL,
                 contrasts = NULL,
                 method = c("joint_laplace", "pirls_profiled"),
                 nAGQ = 1L,
                 inference = c("auto", "none", "asymptotic", "bootstrap",
                               "working_hessian"),
                 control = mm_control(),
                 ...) {
  call <- match.call()
  # lme4 semantics: `weights` and `offset` may name columns of `data`
  # (`weights = trials`), falling back to the calling environment.
  weights <- eval(substitute(weights), data, parent.frame())
  offset <- eval(substitute(offset), data, parent.frame())
  method_explicit <- !missing(method)
  method <- match.arg(method)
  inference <- match.arg(inference)
  control <- mm_validate_control(control)
  family_info <- mm_glmm_family_info(family)
  # Default estimator: joint Laplace (glmer nAGQ = 1). Requests the joint
  # route cannot serve resolve to the profiled path ANNOUNCED, never silently
  # (see mm_glmm_resolve_method); explicit requests are honoured or refused.
  resolved <- mm_glmm_resolve_method(method, method_explicit, family_info,
                                     nAGQ, inference)
  method <- resolved$method
  nAGQ <- mm_glmm_validate_nagq(resolved$nAGQ, method)

  if (identical(family_info$family, "negative_binomial") &&
      identical(method, "joint_laplace")) {
    # The engine's NB theta loop routes through fast-PIRLS internally, so the
    # joint route is not wired for NB at this pin; refusing here gives a
    # clear message instead of the engine's internal optimizer-guard error.
    mm_abort(
      message = paste0(
        "Negative-binomial GLMMs support the profiled method only; ",
        "method = \"joint_laplace\" is not yet available for this family. ",
        "Use method = \"pirls_profiled\" (or drop the method argument, which ",
        "selects it with a notice)."
      ),
      class = "mm_inference_unavailable",
      reason_code = "nb_joint_laplace_unavailable",
      input = method
    )
  }

  if (!is.null(random)) {
    mm_abort(
      message = "`random` is reserved for the native random-effect constructor path; write random terms in `formula`.",
      class = "mm_fit_error",
      input = call
    )
  }
  if (!is.null(contrasts)) mm_reject_nontreatment_contrasts(contrasts, data)
  if (is.data.frame(data)) {
    weights <- mm_glmm_validate_weights(weights, data, "weights")
    offset <- mm_glmm_validate_weights(offset, data, "offset", positive = FALSE)
  }

  # Same lme4-style data preparation as lmm() (R/model-data.R): transforms,
  # subset, na.action, unused levels, grouping coercion, formula expansion.
  # Formula offset() terms are added to `offset`.
  mprep <- mm_prepare_model_data(formula, data, substitute(subset), na.action,
                                 weights, offset, parent.frame(),
                                 control$verbose, lmm = FALSE)
  data <- mprep$data
  weights <- mprep$weights
  offset <- mprep$offset
  engine_formula <- mprep$formula_engine

  # Resolve binomial responses: translate a cbind(successes, failures) LHS into
  # a proportion response + trial-count weights, and pick the engine family
  # ("binomial" for grouped/weighted data, "bernoulli" for 0/1 responses).
  prep <- mm_glmm_binomial_prep(formula, data, family_info, weights)
  formula <- prep$formula
  data <- prep$data
  weights <- prep$weights
  engine_formula[[2L]] <- formula[[2L]]
  if (!is.null(mprep$na_action)) attr(data, "na.action") <- mprep$na_action
  if (!is.null(attr(mprep$data, "mm_predvars"))) {
    attr(data, "mm_predvars") <- attr(mprep$data, "mm_predvars")
  }
  engine_family <- prep$engine_family
  if (identical(family_info$family, "negative_binomial")) {
    # Theta mode rides the family string (see mm_glmm_nb_engine_spec).
    engine_family <- mm_glmm_nb_engine_spec(family_info)
  }

  spec <- compile_model(engine_formula, data)
  mm_validate_fit_structure(spec, lmm = FALSE)
  mm_scaling_advisory(spec, control$verbose)
  if (control$verbose >= 0L) {
    mm_inform_explanation(spec)
    # No silent surgery on the estimator choice: when the user did not pick a
    # method, surface that the default profiled path is NOT glmer's estimator
    # and point to the certified glmer-equivalent route. Suppressed by
    # mm_control(verbose = -1) (as used in loops/bootstrap).
    if (!is.null(resolved$notice)) {
      mm_inform(resolved$notice, class = "mm_estimator_notice")
    }
    if (identical(method, "joint_laplace")) {
      # The joint route runs to an engine-chosen evaluation budget inside a
      # single silent native call (~0.1-0.2s per evaluation on ~20k rows), so
      # large fits can take minutes with no output. Say so up front.
      mm_inform(
        paste0(
          "method = \"joint_laplace\" optimizes to an engine-chosen budget in ",
          "a single native call with no progress output; large fits can take ",
          "minutes. Cap the budget with mm_control(max_feval = ) if needed. ",
          "Silence this message with mm_control(verbose = -1)."
        ),
        class = "mm_runtime_notice"
      )
    }
  }

  spec_data <- mm_translate_data(spec$model_frame)
  formula_string <- mm_coerce_formula_string(engine_formula)
  control_json <- jsonlite::toJSON(unclass(control), auto_unbox = TRUE,
                                   null = "null", digits = NA)

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

  fit_result <- mm_bridge_fit_result(json, mm_json_parse_glmm_fit)
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
    # Exact engine family string used at fit time (bernoulli/binomial split,
    # NB theta-mode spelling). Stored so refit-based verbs (e.g.
    # verify_convergence) reconstruct the same model without re-deriving the
    # prep rules; older saved fits lack it and fall back to reconstruction.
    engine_family  = engine_family,
    method         = as.character(fit_result$method %||% method),
    nAGQ           = as.integer(fit_result$n_agq %||% nAGQ),
    # lme4's nAGQ = 0 maps to the profiled path with nAGQ = 1 on the wire;
    # keep the caller's request on record.
    nAGQ_requested = if (isTRUE(resolved$nagq0)) 0L else as.integer(nAGQ),
    inference_request = inference,
    control        = control,
    vars           = spec$vars,
    model_frame    = data,
    weights        = weights,
    offset         = offset,
    offset_arg     = mprep$offset_arg,
    engine_formula = if (!identical(engine_formula, formula)) engine_formula,
    expansion      = mprep$expansion,
    na.action      = mprep$na_action,
    artifact       = artifact,
    # Raw payload minus the n-length vectors, which are stored once below.
    fit            = fit_result[setdiff(names(fit_result),
                                        c("fitted", "fixed_fitted", "residuals"))],
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
  fit <- mm_apply_lme4_group_labels(fit)
  class(fit) <- c("mm_glmm", "mm_fit", "mm_compiled")
  # No silent surgery on the estimator that produced the numbers: when the
  # engine substituted a fallback for the requested method (typed
  # `estimator_substitution` record, engine f82c646+), say so at fit time.
  # The summary() note repeats this ungated, so verbose = -1 loops stay
  # quiet without hiding the substitution from readers of the result.
  substitution <- mm_estimator_substitution(fit)
  if (!is.null(substitution) && control$verbose >= 0L) {
    mm_inform(
      mm_glmm_substitution_notice(substitution),
      class = "mm_estimator_substitution_notice"
    )
  }
  fit
}

mm_glmm_substitution_notice <- function(sub) {
  sprintf(
    paste0(
      "glmm(method = \"%s\"): the requested estimator did not certify ",
      "(status `%s`, code `%s`); the result is a labelled `%s` fallback fit. ",
      "Its coefficients are the fallback estimator's, not %s's, and Wald ",
      "inference is withheld. Evidence: ",
      "optimizer_certificate(fit)$raw$estimator_substitution; check the ",
      "optimum with verify_convergence(fit)."
    ),
    sub$requested_method, sub$requested_fit_status,
    sub$requested_return_code, sub$effective_method, sub$requested_method
  )
}

# Resolve the estimator. `method` defaults to "joint_laplace" (glmer's
# nAGQ = 1 Laplace estimate). When `method` was NOT given and the request
# needs the profiled path, the profiled path is used and announced (typed
# `mm_estimator_notice`, silenced by verbose = -1):
#   * negative-binomial families (the joint route is not wired for NB);
#   * nAGQ > 1 (adaptive quadrature runs on the profiled path only);
#   * inference = "working_hessian" (an opt-in for the profiled estimator's
#     uncertified working-Hessian Wald approximation);
#   * nAGQ = 0 (lme4's PIRLS-only fast estimate == the profiled path).
# Explicit method requests are never swapped: joint_laplace with NB or with
# nAGQ != 1 is refused downstream.
mm_glmm_resolve_method <- function(method, method_explicit, family_info, nAGQ,
                                   inference = "auto") {
  nagq0 <- is.numeric(nAGQ) && length(nAGQ) == 1L && !is.na(nAGQ) &&
    nAGQ == 0
  if (nagq0) {
    if (method_explicit && identical(method, "joint_laplace")) {
      mm_abort(
        message = paste0(
          "`nAGQ = 0` requests lme4's PIRLS-only fast estimate, which is ",
          "method = \"pirls_profiled\"; it cannot be combined with ",
          "method = \"joint_laplace\"."
        ),
        class = "mm_arg_error",
        input = nAGQ
      )
    }
    # nAGQ = 0 is a request for the fast path itself, so no notice is needed
    # beyond the documentation: it is what the user asked for.
    return(list(method = "pirls_profiled", nAGQ = 1L, notice = NULL,
                nagq0 = TRUE))
  }
  out <- list(method = method, nAGQ = nAGQ, notice = NULL, nagq0 = FALSE)
  if (method_explicit) return(out)
  is_nb <- identical(family_info$family, "negative_binomial")
  agq <- is.numeric(nAGQ) && length(nAGQ) == 1L && !is.na(nAGQ) && nAGQ > 1
  wh <- identical(inference, "working_hessian")
  if (is_nb || agq || wh) {
    out$method <- "pirls_profiled"
    why <- if (is_nb) {
      "the joint Laplace route is not yet available for negative-binomial families"
    } else if (agq) {
      sprintf("adaptive quadrature (nAGQ = %d) runs on the profiled path only",
              as.integer(nAGQ))
    } else {
      paste0("inference = \"working_hessian\" is the profiled estimator's ",
             "opt-in approximation")
    }
    out$notice <- paste0(
      "glmm() is using method = \"pirls_profiled\" instead of the default ",
      "\"joint_laplace\": ", why, ". The profiled (PIRLS) estimate does not ",
      "exactly match lme4::glmer(). Pass method = \"pirls_profiled\" to make ",
      "this explicit, or silence this message with mm_control(verbose = -1)."
    )
  }
  out
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
    family = mm_glmm_engine_family_name(family_name),
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

# R family name -> engine/family_info name.
mm_glmm_engine_family_name <- function(family_name) {
  switch(family_name,
         Gamma = "gamma",
         inverse.gaussian = "inverse_gaussian",
         family_name)
}

# Every family/link pair the engine fits (mixeff-rs
# validate_supported_glmm_family_link). gaussian/identity is an LMM (lmm());
# the engine has no cauchit or 1/mu^2 link, and no binomial log/identity or
# poisson identity GLMM.
mm_glmm_supported_family_links <- function() {
  list(
    binomial = c("logit", "probit", "cloglog"),
    poisson = c("log", "sqrt"),
    Gamma = c("inverse", "log"),
    inverse.gaussian = c("inverse", "log"),
    gaussian = c("log", "inverse", "sqrt"),
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
  supported <- mm_glmm_supported_family_link_table()
  by_family <- vapply(split(supported$link, supported$family),
                      paste, character(1), collapse = "/")
  supported_text <- paste(sprintf("%s (%s)", names(by_family), by_family),
                          collapse = ", ")
  why <- if (identical(family, "gaussian") && identical(link, "identity")) {
    "A Gaussian identity-link mixed model is a linear mixed model; fit it with lmm(). "
  } else if (identical(link, "1/mu^2")) {
    "The engine has no 1/mu^2 link (inverse.gaussian's default); use link = \"inverse\" or \"log\". "
  } else if (identical(link, "cauchit")) {
    "The engine has no cauchit link. "
  } else {
    ""
  }
  mm_abort(
    message = sprintf(
      paste0(
        "The %s family with link `%s` is not supported by glmm(). %sSupported ",
        "families: %s. For other families, use lme4::glmer()."
      ),
      family, link, why, supported_text
    ),
    class = "mm_inference_unavailable",
    reason_code = reason_code,
    family = family,
    link = link,
    supported = supported,
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
  if (!is.numeric(nAGQ) || length(nAGQ) != 1L || is.na(nAGQ) || nAGQ < 1 ||
      nAGQ != round(nAGQ)) {
    mm_abort(
      message = "`nAGQ` must be 0 or a single positive integer.",
      class = "mm_arg_error",
      input = nAGQ
    )
  }
  nAGQ <- as.integer(nAGQ)
  if (identical(method, "joint_laplace") && nAGQ > 1L) {
    mm_abort(
      message = paste0(
        "`method = \"joint_laplace\"` requires `nAGQ = 1`; adaptive ",
        "quadrature (nAGQ > 1) is available with method = \"pirls_profiled\" ",
        "for a single scalar random effect."
      ),
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
