# GLMM fixed-effect inference: Wald contrasts, single-term LRT deletion, and
# model-comparison anova. The upstream contract certifies asymptotic Wald
# (z) inference for GLMMs (the working-Hessian covariance backs summary());
# we build on that with the same conventions as the LMM surface.

## ---- Single inference-capability arbiter (Audit1 WI-2.1) -------------------
##
## Every GLMM route that could emit a standard error, test statistic, p-value,
## or Wald interval must consult THIS function -- never the covariance
## payload's status (which reads available/available_noninferential even for
## the uncertified working Hessian) and never incidental storage like
## fit$std_errors. Sources of truth, in order:
##   1. the engine's certified fixed-effect inference table (payload 1);
##   2. the user's explicit fit-level opt-in `inference = "working_hessian"`,
##      which unlocks the labelled, uncertified working-Hessian approximation.
## An ABSENT inference table means NOT capable (an absent certificate is not
## a certificate).
mm_glmm_inference_capability <- function(fit) {
  parsed <- mm_json_parse_fixed_effect_inference_table(
    fit$artifact$fixed_effect_inference_table %||% NULL
  )
  st <- NULL
  if (!is.null(parsed)) {
    rows <- parsed$table
    if ("kind" %in% names(rows)) {
      rows <- rows[rows$kind == "coefficient", , drop = FALSE]
    }
    if (nrow(rows)) {
      st <- mm_glmm_inference_rows_status(rows)
    }
  }
  status <- as.character(st$status %||% "unavailable")
  if (identical(status, "available")) {
    return(list(
      wald = TRUE,
      source = "certified_inference_table",
      method = "asymptotic_wald_z",
      reliability = as.character(st$reliability %||% "moderate"),
      status = "available",
      reason_code = NA_character_,
      reason = NA_character_
    ))
  }
  if (identical(as.character(fit$inference_request %||% "auto"),
                "working_hessian")) {
    V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
    vstatus <- attr(V, "mm_status") %||% "unavailable"
    if (!is.null(V) &&
        vstatus %in% c("available", "available_noninferential") &&
        all(is.finite(as.matrix(unclass(V))))) {
      return(list(
        wald = TRUE,
        source = "working_hessian_opt_in",
        method = "wald_z_working_hessian",
        reliability = "moderate",
        status = "available_noninferential",
        reason_code = NA_character_,
        reason = paste0(
          "uncertified working-Hessian approximation requested via ",
          "inference = \"working_hessian\""
        )
      ))
    }
  }
  list(
    wald = FALSE,
    source = "certified_inference_table",
    method = "not_computed",
    reliability = "not_available",
    status = status,
    reason_code = "glmm_wald_uncertified",
    reason = as.character(
      st$reason %||%
        "certified GLMM Wald inference is unavailable for this estimator"
    )
  )
}

## One refusal message for every Wald-refusing GLMM route: the two certified
## paths first, then the single discovery pointer. Working-Hessian is
## deliberately NOT named here -- it is discoverable via inference_options()
## with its caveat attached (Audit1 DEC-1 progressive disclosure).
mm_glmm_wald_refusal_message <- function(what = "Standard errors and Wald tests") {
  paste0(
    what, " are withheld for this GLMM estimator: the engine does not ",
    "certify the profiled working-Hessian covariance for Wald inference. ",
    "Use confint(fit, method = \"bootstrap\") for parametric-bootstrap ",
    "intervals (all supported families), or refit with ",
    "method = \"joint_laplace\" for certified Wald inference (not available ",
    "for negative binomial). inference_options(fit) lists every route."
  )
}

## Status list in the shape routes attach as an attribute (and
## mm_lincomb_status_from_vcov historically produced).
mm_glmm_capability_status <- function(cap) {
  list(
    status = cap$status,
    method = cap$method,
    reliability = cap$reliability,
    reason = cap$reason
  )
}

#' @rdname contrast
#' @details
#' For `mm_glmm` fits, contrasts use an asymptotic Wald z-test built from the
#' stored fixed-effect covariance (the GLMM contract does not provide
#' finite-sample df), so `method` accepts only `"asymptotic"` (alias `"wald"`).
#' @method contrast mm_glmm
#' @export
contrast.mm_glmm <- function(fit, L, rhs = 0,
                             method = c("asymptotic", "wald"), ...) {
  method <- match.arg(method)
  if (identical(method, "wald")) method <- "asymptotic"
  L <- mm_contrast_matrix(L, fit)
  rhs <- rep(as.numeric(rhs), length.out = nrow(L))
  beta <- fit$beta
  Vfull <- stats::vcov(fit)
  V <- as.matrix(unclass(Vfull))
  dimnames(V) <- list(names(beta), names(beta))

  estimate <- as.numeric(L %*% beta) - rhs

  # Gate Wald inference through the single capability arbiter (WI-2.1): the
  # certified inference table, plus the explicit working-Hessian opt-in.
  # For an incapable fit we report the contrast estimate but withhold SE/z/p
  # rather than fabricate inference from an uncertified covariance ("no fake
  # certainty"), matching summary.mm_glmm / confint.mm_glmm / tidy.mm_glmm.
  cap <- mm_glmm_inference_capability(fit)
  if (isTRUE(cap$wald)) {
    se <- sqrt(diag(L %*% V %*% t(L)))
    statistic <- estimate / se
    p <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
    statistic_name <- "z"
    row_method <- cap$method
    reason <- cap$reason
  } else {
    se <- rep(NA_real_, nrow(L))
    statistic <- rep(NA_real_, nrow(L))
    p <- rep(NA_real_, nrow(L))
    statistic_name <- NA_character_
    row_method <- "not_computed"
    reason <- cap$reason %||% "certified GLMM Wald inference is unavailable"
  }

  table <- data.frame(
    contrast = rownames(L),
    estimate = estimate,
    rhs = rhs,
    std_error = se,
    df = NA_real_,
    statistic = statistic,
    statistic_name = statistic_name,
    p_value = p,
    method = row_method,
    requested_method = method,
    status = cap$status,
    reliability = cap$reliability,
    estimability = "not_assessed",
    reason = reason,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  obj <- list(table = table, L = L, rhs = rhs, requested_method = method,
              raw = NULL)
  class(obj) <- "mm_contrast"
  obj
}

#' Drop one fixed-effect term at a time from a GLMM
#'
#' Refits reduced fixed-effect GLMMs (random-effect terms preserved exactly)
#' and compares each to the full fit by asymptotic likelihood-ratio test,
#' mirroring `drop1(glmerMod, test = "Chisq")`.
#'
#' @param object A fitted `mm_glmm`.
#' @param scope Optional character vector of fixed-effect terms to drop.
#' @param test `"Chisq"` reports asymptotic LRT rows; `"none"` reports
#'   information criteria only.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_drop1` object.
#'
#' @method drop1 mm_glmm
#' @export
drop1.mm_glmm <- function(object, scope = NULL, test = c("none", "Chisq"),
                          ...) {
  test <- match.arg(test)
  terms <- setdiff(mm_fixed_effect_terms(object), "1")
  if (!is.null(scope)) {
    terms <- intersect(terms, as.character(scope))
  }
  family <- mm_glmm_family_from_info(object$family)
  rows <- lapply(terms, function(term) {
    reduced_formula <- mm_drop_fixed_term_formula(object, term)
    reduced <- glmm(reduced_formula, object$model_frame, family = family,
                    weights = object$weights, offset = object$offset,
                    method = object$method, nAGQ = object$nAGQ,
                    control = mm_control(verbose = -1))
    stat <- mm_lrt_stat(reduced, object)
    df <- object$dof - reduced$dof
    data.frame(
      dropped = term,
      formula = deparse1(reduced_formula),
      df = df,
      logLik = as.numeric(logLik(reduced)),
      AIC = AIC(reduced),
      BIC = BIC(reduced),
      LRT = if (identical(test, "Chisq")) stat else NA_real_,
      p_value = if (identical(test, "Chisq") && df > 0) {
        stats::pchisq(stat, df = df, lower.tail = FALSE)
      } else {
        NA_real_
      },
      method = if (identical(test, "Chisq")) "asymptotic_lrt" else "none",
      stringsAsFactors = FALSE
    )
  })
  table <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(dropped = character(), formula = character(), df = numeric(),
               logLik = numeric(), AIC = numeric(), BIC = numeric(),
               LRT = numeric(), p_value = numeric(), method = character(),
               stringsAsFactors = FALSE)
  }
  rownames(table) <- NULL
  obj <- list(table = table, full = object)
  class(obj) <- "mm_drop1"
  obj
}

#' Analysis of deviance for GLMMs
#'
#' With two or more fitted models, `anova()` performs a sequential
#' likelihood-ratio comparison (like `anova(glmer1, glmer2)`). For a single
#' model, fixed-effect tests are routed to [drop1()] (term LRTs),
#' [summary()] (Wald z), or [contrast()] (custom Wald contrasts), which the
#' GLMM contract supports directly.
#'
#' @param object A fitted `mm_glmm`.
#' @param ... Additional fitted models to compare.
#'
#' @return An `mm_model_comparison` object (multi-model case).
#'
#' @method anova mm_glmm
#' @export
anova.mm_glmm <- function(object, ...) {
  dots <- list(...)
  if (!length(dots)) {
    mm_abort(
      message = paste(
        "Single-model `anova()` is not provided for GLMMs. Use",
        "`drop1(fit, test = \"Chisq\")` for term likelihood-ratio tests,",
        "`summary(fit)` for Wald z-tests, or `contrast(fit, L)` for custom",
        "Wald contrasts."
      ),
      class = "mm_inference_unavailable",
      reason_code = "glmm_single_model_anova_unavailable"
    )
  }
  fits <- c(list(object), dots)
  if (!all(vapply(fits, inherits, logical(1), "mm_glmm"))) {
    mm_abort(
      message = "All models compared with `anova()` must be `mm_glmm` fits.",
      class = "mm_arg_error"
    )
  }
  mm_glmm_lrt_comparison(fits)
}

# Sequential likelihood-ratio comparison of nested GLMMs, ordered by parameter
# count. GLMM likelihoods are ML, so no REML refit is needed.
mm_glmm_lrt_comparison <- function(fits) {
  mm_assert_comparable_lmm(fits)
  npar <- vapply(fits, function(f) as.integer(f$dof), integer(1))
  ord <- order(npar)
  fits <- fits[ord]
  npar <- npar[ord]
  loglik <- vapply(fits, function(f) as.numeric(f$logLik), numeric(1))
  dev <- vapply(fits, function(f) as.numeric(f$deviance), numeric(1))
  aic <- vapply(fits, function(f) AIC(f), numeric(1))
  bic <- vapply(fits, function(f) BIC(f), numeric(1))
  chisq <- c(NA_real_, pmax(0, -diff(dev)))
  ddf <- c(NA_integer_, diff(npar))
  pval <- c(NA_real_,
            stats::pchisq(chisq[-1L], df = ddf[-1L], lower.tail = FALSE))
  table <- data.frame(
    model = vapply(fits, function(f) deparse1(f$formula), character(1)),
    npar = npar,
    AIC = aic,
    BIC = bic,
    logLik = loglik,
    deviance = dev,
    Chisq = chisq,
    Df = ddf,
    p_value = pval,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(table) <- NULL
  obj <- list(table = table, method = "asymptotic_lrt")
  class(obj) <- "mm_glmm_comparison"
  obj
}

#' @method print mm_glmm_comparison
#' @export
print.mm_glmm_comparison <- function(x, ...) {
  cat("GLMM model comparison (sequential likelihood-ratio test):\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

## ---- Typed refusals for GLMM verbs without a validated route (WI-2.4) -----
##
## These verbs exist for LMM fits only. Before this, a GLMM user got a bare
## S3 "no applicable method" dispatch error -- an undeclared gap rather than
## a typed refusal. Each refusal names the closest supported route.

mm_glmm_verb_refusal <- function(verb, alternative) {
  mm_abort(
    message = paste0(
      "`", verb, "()` is not implemented for GLMM fits in this version. ",
      alternative
    ),
    class = "mm_inference_unavailable",
    reason_code = "glmm_verb_unavailable"
  )
}

#' @method test_effect mm_glmm
#' @export
test_effect.mm_glmm <- function(fit, term, ...) {
  mm_glmm_verb_refusal(
    "test_effect",
    paste0(
      "For single-term tests on a GLMM, use drop1() (likelihood-ratio) or ",
      "contrast() on a certified fit (method = \"joint_laplace\")."
    )
  )
}

#' @method mm_means mm_glmm
#' @export
mm_means.mm_glmm <- function(fit, specs, ...) {
  mm_glmm_verb_refusal(
    "mm_means",
    paste0(
      "For GLMM marginal means, use the emmeans bridge on a certified fit ",
      "(method = \"joint_laplace\")."
    )
  )
}

#' @method mm_comparisons mm_glmm
#' @export
mm_comparisons.mm_glmm <- function(fit, specs, ...) {
  mm_glmm_verb_refusal(
    "mm_comparisons",
    paste0(
      "For GLMM marginal comparisons, use the emmeans bridge on a certified ",
      "fit (method = \"joint_laplace\")."
    )
  )
}

#' @method mm_predictions mm_glmm
#' @export
mm_predictions.mm_glmm <- function(fit, ...) {
  mm_glmm_verb_refusal(
    "mm_predictions",
    "For GLMM predictions with uncertainty, use predict(fit, se.fit = TRUE), which carries its own per-row engine certificate."
  )
}

#' @method mm_grid mm_glmm
#' @export
mm_grid.mm_glmm <- function(fit, specs, ...) {
  mm_glmm_verb_refusal(
    "mm_grid",
    paste0(
      "For GLMM reference grids, use the emmeans bridge on a certified fit ",
      "(method = \"joint_laplace\")."
    )
  )
}

## ---- GLMM parametric bootstrap (WI-2.8) ------------------------------------
##
## The universal inference path for GLMMs: works for every supported family
## cell, including negative binomial (which has no joint-Laplace route).
## Replicates are simulated from the fitted model under fresh random-effect
## draws and refitted with the EFFECTIVE estimator, so the intervals
## characterize the sampling distribution of the estimator that actually
## produced the reported numbers. NB replicate refits condition on the
## template's fitted theta. Failed replicates are counted, never silently
## dropped. Deterministic given `seed`; when `seed` is NULL one is drawn
## from R's RNG so set.seed() governs reproducibility.
mm_glmm_parametric_bootstrap_confint <- function(object, parm, level,
                                                 nsim = 999L, seed = NULL) {
  if (identical(mm_glmm_effective_method(object), "joint_laplace")) {
    # Engine limitation at pin f82c646: GeneralizedLinearMixedModel::refit
    # hardcodes the fast path (generalized/optimizer.rs:19-23) and the
    # dependency-light build then refuses the joint optimizer left in
    # optsum, so every replicate refit of a joint template fails
    # deterministically. Refuse up front with the honest reason instead of
    # burning nsim doomed refits (upstream bead filed; certified Wald is
    # available on joint fits anyway).
    mm_abort(
      message = paste0(
        "The parametric bootstrap is not available for joint_laplace fits ",
        "at this engine pin: replicate refits cannot re-run the joint ",
        "estimator. Use the certified Wald intervals ",
        "(confint(fit, method = \"asymptotic\")), or bootstrap the profiled ",
        "estimator by refitting with the default method."
      ),
      class = "mm_inference_unavailable",
      reason_code = "glmm_bootstrap_joint_laplace_unavailable",
      input = object$method
    )
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      nsim < 2 || nsim != as.integer(nsim)) {
    mm_abort(
      message = "`nsim` must be a single whole number >= 2.",
      class = "mm_arg_error",
      input = nsim
    )
  }
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0) {
    mm_abort(
      message = "`seed` must be a single non-negative integer.",
      class = "mm_arg_error",
      input = seed
    )
  }

  payload <- mm_rust_glmm_refit_payload(object, mm_glmm_effective_method(object))
  options_json <- jsonlite::toJSON(
    list(nsim = as.integer(nsim), seed = as.numeric(seed)),
    auto_unbox = TRUE
  )
  json <- tryCatch(
    mm_glmm_parametric_bootstrap_json(payload, as.character(options_json)),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json)
  }
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  # json_f64 encodes non-finite values as the strings "NaN"/"Infinity".
  num1 <- function(v) suppressWarnings(as.numeric(v %||% NA_real_))
  beta_names <- mm_coef_engine_to_lme4(
    as.character(unlist(parsed$beta_names)),
    object$coef_map
  )
  reps <- parsed$bootstrap$fits %||% list()
  p <- length(beta_names)
  objective <- vapply(reps, function(r) num1(r$objective), numeric(1))
  beta_mat <- t(vapply(
    reps,
    function(r) {
      b <- suppressWarnings(as.numeric(unlist(r$beta)))
      if (length(b) != p) rep(NA_real_, p) else b
    },
    numeric(p)
  ))
  colnames(beta_mat) <- beta_names
  ok <- is.finite(objective) &
    apply(beta_mat, 1L, function(b) all(is.finite(b)))
  successful <- sum(ok)
  failed <- length(reps) - successful
  if (successful < 2L) {
    mm_abort(
      message = sprintf(
        paste0("GLMM parametric bootstrap produced %d successful ",
               "replicate(s) out of %d; intervals need at least 2. The ",
               "replicate accounting is attached to this condition."),
        successful, length(reps)
      ),
      class = "mm_inference_unavailable",
      reason_code = "glmm_bootstrap_insufficient_replicates",
      requested = length(reps),
      successful = successful,
      failed = failed
    )
  }

  keep <- beta_mat[ok, parm, drop = FALSE]
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  ci <- t(apply(keep, 2L, stats::quantile, probs = probs, names = FALSE))
  dimnames(ci) <- list(parm, paste(format(100 * probs, trim = TRUE), "%"))
  se <- apply(keep, 2L, stats::sd)
  mcse <- max(se / sqrt(successful))

  attr(ci, "mm_method") <- "glmm_parametric_bootstrap_percentile"
  attr(ci, "mm_estimator") <- mm_glmm_effective_method(object)
  attr(ci, "mm_bootstrap") <- list(
    requested = length(reps),
    successful = successful,
    failed = failed,
    seed = as.numeric(seed),
    std_errors = se,
    mcse = mcse,
    reliability = mm_bootstrap_reliability(TRUE, successful, mcse)
  )
  ci
}
