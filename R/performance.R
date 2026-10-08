# Nakagawa R2 and ICC (pre-CRAN audit 4.8): native mm_r2()/mm_icc(), plus
# methods registered into insight/performance (both in Suggests) so that
# performance::icc(), performance::r2() and insight::get_variance() work on
# mixeff fits. The variance decomposition follows insight's get_variance()
# for merMod fits (insight 0.19.x).

#' Variance decomposition, Nakagawa R2 and ICC for mixeff fits
#'
#' `mm_variance_components()` splits the variance of a fitted mixed model
#' into the components used by Nakagawa and Schielzeth's R2 and the
#' intraclass correlation: the variance of the fixed-effect linear predictor,
#' the random-effect variance, and the residual (distribution-specific plus
#' observation-level) variance. `mm_r2()` returns the marginal and
#' conditional R2 and `mm_icc()` the adjusted and unadjusted ICC. The
#' numbers follow `insight::get_variance()`, `performance::r2_nakagawa()` and
#' `performance::icc()` on the equivalent lme4 fit.
#'
#' When the insight and performance packages are installed, mixeff also
#' registers methods so that `performance::icc(fit)`, `performance::r2(fit)`
#' and `insight::get_variance(fit)` work directly on mixeff fits (not
#' `by_group = TRUE`, which needs lme4 internals).
#'
#' * Fixed variance: `var(X beta)`.
#' * Random variance: for each random term whose grouping factor has fewer
#'   levels than observations, the mean over observations of `z_i' Sigma z_i`,
#'   where `z_i` holds the term's coefficients' columns of the fixed-effect
#'   design (a random slope without a matching fixed column is dropped, as
#'   insight does).
#' * Distribution-specific variance: `sigma^2` for an LMM; `pi^2/3`, `1` and
#'   `pi^2/6` for binomial logit, probit and cloglog links; `sigma^2` for
#'   Gamma; and the log-normal approximation `log(1 + V(mu)/mu^2)` for Poisson
#'   (`V(mu) = mu`) and negative-binomial (`V(mu) = mu (1 + mu/theta)`) log-link
#'   models, with `mu = exp(b0)` from the intercept of the null model
#'   (`y ~ 1` with the same random effects), which is refitted.
#' * Observation-level random effects (one level per observation) count
#'   towards the residual variance of a GLMM.
#'
#' If a random-effect variance is below `tolerance` (a singular fit), the
#' random variance is not computed: `mm_r2()` and `mm_icc()` return `NA`
#' for the quantities that need it and signal a warning, as performance does.
#'
#' @param fit A fitted `mm_lmm` or `mm_glmm`.
#' @param tolerance Variances below this value mark the fit as singular.
#'
#' @return `mm_variance_components()`: a named list with `var.fixed`,
#'   `var.random`, `var.residual`, `var.distribution`, `var.dispersion`,
#'   `var.intercept`, `var.slope` and `cor.slope_intercept` (the layout of
#'   `insight::get_variance()`). `mm_r2()`: a list with `R2_conditional` and
#'   `R2_marginal`. `mm_icc()`: a one-row data frame with `ICC_adjusted`,
#'   `ICC_conditional` and `ICC_unadjusted`.
#'
#' @examples
#' set.seed(1)
#' g <- factor(rep(seq_len(10), each = 6))
#' df <- data.frame(x = rnorm(60), g = g)
#' df$y <- 1 + 0.5 * df$x + rnorm(10)[g] + rnorm(60)
#' fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
#' mm_r2(fit)
#' mm_icc(fit)
#'
#' @export
mm_variance_components <- function(fit, tolerance = 1e-5) {
  if (!inherits(fit, c("mm_lmm", "mm_glmm"))) {
    mm_abort(
      message = "`fit` must be a fitted mixeff model (`mm_lmm` or `mm_glmm`).",
      class = "mm_arg_error",
      input = class(fit)
    )
  }
  mm_variance_parts(fit, tolerance)
}

mm_variance_parts <- function(fit, tolerance = 1e-5, warn = TRUE) {
  X <- as.matrix(stats::model.matrix(fit, type = "fixed"))
  beta <- as.numeric(fixef(fit))
  n <- nrow(X)
  var_fixed <- stats::var(as.vector(X %*% beta))

  terms <- fit$artifact$semantic_model$random_terms %||% list()
  labels <- vapply(seq_along(terms), function(i) {
    mm_random_term_group_label(fit, terms[[i]], i)
  }, character(1))
  vc <- lapply(seq_along(terms), function(i) {
    term <- terms[[i]]
    cn <- vapply(term$basis %||% list(), mm_basis_label, character(1))
    if (!length(cn)) cn <- "(Intercept)"
    mm_random_term_covariance(fit, term$id %||% sprintf("r%d", i - 1L), cn,
                              labels[[i]])
  })
  names(vc) <- make.unique(labels)
  nlev <- vapply(labels, function(l) nlevels(mm_group_factor(fit$model_frame, l)),
                 numeric(1))
  obs_level <- nlev == n

  Xi <- X
  if (!"(Intercept)" %in% colnames(Xi)) {
    Xi <- cbind(`(Intercept)` = 1, Xi)
  }
  sigma_sum <- function(S) {
    rn <- intersect(rownames(S), colnames(Xi))
    if (!length(rn)) return(0)
    S <- S[rn, rn, drop = FALSE]
    Z <- Xi[, rn, drop = FALSE]
    sum((Z %*% S) * Z) / n
  }
  singular <- any(vapply(vc, function(S) any(abs(diag(S)) < tolerance),
                         logical(1)))
  var_random <- if (singular) {
    if (warn) {
      warning(
        "Some random-effect variances are (near) zero, so the random-effect ",
        "variance is not computed (singular fit; lower `tolerance` to force it).",
        call. = FALSE
      )
    }
    NULL
  } else {
    sum(vapply(vc[!obs_level], sigma_sum, numeric(1)))
  }

  var_dispersion <- 0
  if (inherits(fit, "mm_glmm")) {
    var_dispersion <- sum(vapply(vc[obs_level], sigma_sum, numeric(1)))
    var_distribution <- mm_distribution_variance(fit, warn = warn)
  } else {
    var_distribution <- fit$sigma^2
  }

  intercepts <- vapply(vc, function(S) {
    identical(rownames(S)[[1L]], "(Intercept)")
  }, logical(1))
  var_intercept <- vapply(vc[intercepts], function(S) S[1L, 1L], numeric(1))
  var_slope <- unlist(lapply(names(vc), function(g) {
    S <- vc[[g]]
    keep <- rownames(S) != "(Intercept)"
    if (!any(keep)) return(NULL)
    stats::setNames(diag(S)[keep], paste(g, rownames(S)[keep], sep = "."))
  }))

  cor_si <- unlist(lapply(names(vc)[intercepts], function(g) {
    S <- vc[[g]]
    if (nrow(S) < 2L || any(diag(S) <= 0)) return(NULL)
    r <- stats::cov2cor(S)[-1L, 1L]
    names(r) <- if (length(r) == 1L) g else paste(g, rownames(S)[-1L], sep = ".")
    r
  }))

  out <- list(
    var.fixed = var_fixed,
    var.random = var_random,
    var.residual = var_distribution + var_dispersion,
    var.distribution = var_distribution,
    var.dispersion = var_dispersion,
    var.intercept = if (length(var_intercept)) var_intercept,
    var.slope = var_slope,
    cor.slope_intercept = cor_si
  )
  out[!vapply(out, is.null, logical(1))]
}

mm_distribution_variance <- function(fit, warn = TRUE) {
  fam <- fit$family$family
  link <- fit$family$link
  bad_link <- function() {
    mm_abort(
      message = sprintf(
        "The distribution-specific variance is not defined for family `%s` with link `%s`.",
        fam, link
      ),
      class = "mm_inference_unavailable",
      reason_code = "r2_family_link_unavailable",
      input = c(fam, link)
    )
  }
  switch(
    fam,
    binomial = switch(link, logit = pi^2 / 3, probit = 1, cloglog = pi^2 / 6,
                      bad_link()),
    gamma = fit$dispersion^2,
    poisson = ,
    negative_binomial = {
      if (!identical(link, "log")) {
        if (identical(link, "sqrt")) return(0.25)
        bad_link()
      }
      mu <- exp(mm_null_model_intercept(fit))
      if (warn && mu < 6) {
        warning(sprintf(
          "mu of %0.1f is too close to zero; the distribution-specific variance may be unreliable.",
          mu
        ), call. = FALSE)
      }
      v <- if (identical(fam, "poisson")) {
        mu
      } else {
        mu * (1 + mu / fit$family$nb_theta)
      }
      log1p(v / mu^2)
    },
    bad_link()
  )
}

# Intercept of the null model y ~ 1 + <same random terms>, refitted with the
# original family, weights, offset and method.
mm_null_model_intercept <- function(fit) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  re <- unique(vapply(terms, function(t) {
    mm_scalar_text(t$source_syntax$written %||% t$source_syntax$text)
  }, character(1)))
  f <- stats::as.formula(
    paste(mm_response_name(fit), "~", paste(c("1", re), collapse = " + ")),
    env = environment(fit$formula)
  )
  control <- fit$control %||% mm_control(verbose = -1)
  control$verbose <- -1L
  null_fit <- do.call(glmm, list(
    formula = f, data = fit$model_frame,
    family = mm_glmm_family_from_info(fit$family),
    weights = fit$weights, offset = fit$offset,
    method = fit$method, nAGQ = fit$nAGQ,
    inference = "none", control = control
  ))
  unname(fixef(null_fit)[[1L]])
}

#' @rdname mm_variance_components
#' @export
mm_r2 <- function(fit, tolerance = 1e-5) {
  v <- mm_variance_components(fit, tolerance)
  if (is.null(v$var.random)) {
    r2_marginal <- v$var.fixed / (v$var.fixed + v$var.residual)
    r2_conditional <- NA_real_
  } else {
    total <- v$var.fixed + v$var.random + v$var.residual
    r2_marginal <- v$var.fixed / total
    r2_conditional <- (v$var.fixed + v$var.random) / total
  }
  list(R2_conditional = c(`Conditional R2` = r2_conditional),
       R2_marginal = c(`Marginal R2` = r2_marginal))
}

#' @rdname mm_variance_components
#' @export
mm_icc <- function(fit, tolerance = 1e-5) {
  v <- mm_variance_components(fit, tolerance)
  if (is.null(v$var.random)) {
    return(data.frame(ICC_adjusted = NA_real_, ICC_conditional = NA_real_,
                      ICC_unadjusted = NA_real_))
  }
  adjusted <- v$var.random / (v$var.random + v$var.residual)
  unadjusted <- v$var.random / (v$var.fixed + v$var.random + v$var.residual)
  data.frame(ICC_adjusted = adjusted, ICC_conditional = unadjusted,
             ICC_unadjusted = unadjusted)
}

# ---- insight / performance integration (registered lazily, see zzz.R) ----

get_variance.mm_fit <- function(x, component = c("all", "fixed", "random",
                                                 "residual", "distribution",
                                                 "dispersion", "intercept",
                                                 "slope"),
                                verbose = TRUE, tolerance = 1e-5, ...) {
  component <- match.arg(component)
  v <- tryCatch(mm_variance_parts(x, tolerance, warn = verbose),
                error = function(cnd) {
                  if (verbose) warning(conditionMessage(cnd), call. = FALSE)
                  NULL
                })
  if (is.null(v)) return(NA)
  if (identical(component, "all")) return(v)
  v[[paste0("var.", component)]]
}

is_mixed_model.mm_fit <- function(x) TRUE

find_random_slopes.mm_fit <- function(x) {
  terms <- x$artifact$semantic_model$random_terms %||% list()
  slopes <- unique(unlist(lapply(terms, function(t) {
    cn <- vapply(t$basis %||% list(), mm_basis_label, character(1))
    setdiff(cn, "(Intercept)")
  })))
  if (!length(slopes)) return(NULL)
  list(random = slopes)
}

r2.mm_fit <- function(model, ...) {
  performance::r2_nakagawa(model, ...)
}

mm_register_insight_s3 <- function() {
  ns <- asNamespace("insight")
  for (cls in c("mm_lmm", "mm_glmm")) {
    registerS3method("get_variance", cls, get_variance.mm_fit, envir = ns)
    registerS3method("is_mixed_model", cls, is_mixed_model.mm_fit, envir = ns)
    registerS3method("find_random_slopes", cls, find_random_slopes.mm_fit,
                     envir = ns)
  }
  invisible(TRUE)
}

mm_register_performance_s3 <- function() {
  ns <- asNamespace("performance")
  registerS3method("r2", "mm_lmm", r2.mm_fit, envir = ns)
  registerS3method("r2", "mm_glmm", r2.mm_fit, envir = ns)
  invisible(TRUE)
}
