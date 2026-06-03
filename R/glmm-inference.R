# GLMM fixed-effect inference: Wald contrasts, single-term LRT deletion, and
# model-comparison anova. The upstream contract certifies asymptotic Wald
# (z) inference for GLMMs (the working-Hessian covariance backs summary());
# we build on that with the same conventions as the LMM surface.

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

  # Gate Wald inference on the certified fixed-effect inference table -- the same
  # signal summary()/confint() consume. vcov()'s mm_status reads "available" even
  # for the uncertified fast-PIRLS working Hessian, so it cannot be the gate (see
  # the GLMM Wald uptake / upstream bd-01KT3Z64YE5QN7626PQRJSJJVA). For an
  # uncertified fit we report the contrast estimate but withhold SE/z/p rather
  # than fabricate inference from an uncertified covariance ("no fake certainty"),
  # matching summary.mm_glmm / confint.mm_glmm and tidy.mm_glmm.
  status <- attr(mm_glmm_wald_z_inference(fit), "mm_vcov_status")
  certified <- identical(status$status %||% "available", "available")
  if (certified) {
    se <- sqrt(diag(L %*% V %*% t(L)))
    statistic <- estimate / se
    p <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
    statistic_name <- "z"
    row_method <- "asymptotic"
    reason <- NA_character_
  } else {
    se <- rep(NA_real_, nrow(L))
    statistic <- rep(NA_real_, nrow(L))
    p <- rep(NA_real_, nrow(L))
    statistic_name <- NA_character_
    row_method <- "not_computed"
    reason <- status$reason %||% "certified GLMM Wald inference is unavailable"
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
    status = status,
    reliability = "asymptotic_wald",
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
