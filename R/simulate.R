#' Refit a mixeff model with a new response
#'
#' `refit()` fits the same model to a new response, keeping every other
#' setting of the original fit, like `lme4::refit()`. For an `mm_lmm` it calls
#' [lmm()] with the stored model frame, `REML` setting and prior weights; for
#' an `mm_glmm` it calls [glmm()] with the stored family (negative-binomial
#' fits re-estimate theta when the original did), prior weights, offset,
#' `method`, `nAGQ` and control settings (verbosity silenced).
#'
#' `newresp` may be a numeric vector, or a one-column data frame such as the
#' output of [simulate()] with `nsim = 1`. For binomial GLMMs it may also be a
#' two-column `(successes, failures)` matrix (the shape `simulate()` returns
#' for a `cbind()` response; the trial counts become the new prior weights),
#' a two-level factor, or a logical vector.
#'
#' @param object A fitted `mm_lmm` or `mm_glmm`.
#' @param newresp The new response (see Details).
#' @param ... `control =` to override the stored control settings.
#'
#' @return A new fit of the same class as `object`.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rpois(60, 3), x = rnorm(60),
#'   g = factor(rep(seq_len(10), each = 6))
#' )
#' fit <- glmm(y ~ x + (1 | g), df, family = poisson(),
#'             control = mm_control(verbose = -1))
#' ystar <- simulate(fit, nsim = 1, seed = 2)
#' fixef(refit(fit, ystar))
#'
refit <- function(object, newresp, ...) {
  UseMethod("refit")
}

#' @rdname refit
#' @export
refit.default <- function(object, newresp, ...) {
  if (inherits(object, "merMod") && requireNamespace("lme4", quietly = TRUE)) {
    if (missing(newresp)) {
      return(lme4::refit(object, ...))
    }
    return(lme4::refit(object, newresp, ...))
  }
  mm_abort(
    message = "`refit()` has no method for this object.",
    class = "mm_arg_error",
    input = object
  )
}

#' @rdname refit
#' @export
refit.mm_lmm <- function(object, newresp, ...) {
  mm_reject_unsupported_dots(
    list(...), "refit",
    c(newweights = "changing prior weights on refit is not supported; refit with lmm(..., weights = ).")
  )
  newresp <- mm_refit_unwrap_resp(newresp)
  if (!is.numeric(newresp) || length(newresp) != nobs(object) ||
      anyNA(newresp)) {
    mm_abort(
      message = "`newresp` must be a numeric vector with one value per observation and no missing values.",
      class = "mm_arg_error",
      input = newresp
    )
  }
  data <- object$model_frame
  # An LMM offset is fitted as the adjusted response `.mm_offset_response`;
  # the new response replaces the user's response column and the offset is
  # re-applied by lmm().
  response <- if (!is.null(object$offset)) {
    all.vars(object$formula[[2L]])[[1L]]
  } else {
    mm_response_name(object)
  }
  data[[response]] <- as.numeric(newresp)
  control <- list(...)$control %||% mm_control(verbose = -1)
  fit <- lmm(object$formula, data, REML = isTRUE(object$REML),
             weights = object$weights, offset = mm_fit_offset_arg(object),
             control = control)
  fit$refit <- list(
    source = "refit",
    original_fit_status = fit_status(object)
  )
  fit
}

#' Simulate responses from a mixeff fit
#'
#' Draws responses from a fitted `mm_lmm` or `mm_glmm` following the
#' semantics of `lme4::simulate.merMod()`:
#'
#' * `re.form = NA` (the default) or `~0`, equivalently `use.u = FALSE`,
#'   draws **new** random effects from the fitted random-effect covariance
#'   for every simulation (unconditional, parametric-bootstrap simulation);
#' * `re.form = NULL`, equivalently `use.u = TRUE`, conditions on the fitted
#'   conditional modes (BLUPs) and simulates only the response noise around
#'   `fitted(object)`.
#'
#' Partial random-effect formulas are refused with a typed
#' `mm_inference_unavailable` error.
#'
#' LMM responses are Gaussian with standard deviation `sigma / sqrt(w)` for
#' prior weights `w` (lme4 ignores the prior weights here). GLMM responses are
#' family-aware: binomial draws use the prior weights as trial counts and
#' return proportions (`y / w`, 0/1 values for Bernoulli fits), or a
#' two-column `(successes, failures)` matrix per simulation when the model was
#' fitted with a `cbind()` response, as lme4 does; Poisson and
#' negative-binomial (fitted `theta`) draws return counts; Gamma draws use
#' shape `w / phi` with dispersion `phi = sigma(object)^2` and random effects
#' on their fitted (absolute) scale. (lme4 1.1 uses shape `sigma(object)`
#' for Gamma and draws Gamma random effects from the unscaled relative factor;
#' mixeff deliberately uses the model's own variance function instead.) A
#' binomial fit whose response was a factor is simulated as 0/1 numeric
#' (lme4 returns a factor).
#'
#' The random-number stream is consumed in lme4's order (all random-effect
#' draws first, terms in lme4's internal order, then the response draws), so
#' when the estimates agree a seeded `simulate()` reproduces lme4's draws.
#' Unlike lme4, a supplied `seed` does not leave the global RNG state changed.
#'
#' @param object A fitted `mm_lmm` or `mm_glmm`.
#' @param nsim Number of simulated responses.
#' @param seed Optional random seed.
#' @param use.u Logical; `TRUE` is the same as `re.form = NULL` and `FALSE`
#'   the same as `re.form = NA`. Specify at most one of `use.u` and `re.form`.
#' @param re.form `NA` or `~0` (default `NA`) to simulate new random effects,
#'   `NULL` to condition on the fitted random effects, or a formula such as
#'   `~ (1 | g)` to condition on those terms' fitted random effects and draw
#'   new ones for the others (as lme4).
#' @param ... Reserved; lme4 arguments that cannot be honoured (`newdata`,
#'   `newparams`, ...) are refused with a typed error.
#'
#' @return A data frame with one column per simulation (`sim_1`, ...), rows
#'   named like the model frame. For a binomial `cbind()` fit each column is a
#'   two-column matrix. Attributes `seed`, `mm_method` and `mm_re_form` record
#'   the seed, the simulation route, and the resolved conditioning
#'   (`"unconditional"` or `"conditional"`).
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rnorm(60), x = rnorm(60),
#'   g = factor(rep(seq_len(10), each = 6))
#' )
#' fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
#' head(simulate(fit, nsim = 2, seed = 1))
#' # condition on the fitted random effects
#' head(simulate(fit, nsim = 1, seed = 1, re.form = NULL))
#'
#' df$count <- rpois(60, 2)
#' gfit <- glmm(count ~ x + (1 | g), df, family = poisson(),
#'              control = mm_control(verbose = -1))
#' head(simulate(gfit, nsim = 2, seed = 1))
#'
#' @importFrom stats simulate
#' @method simulate mm_lmm
#' @export
simulate.mm_lmm <- function(object, nsim = 1, seed = NULL, use.u = FALSE,
                            re.form = NA, ...) {
  mm_reject_unsupported_dots(
    list(...), "simulate",
    c(newparams = "simulating from user-supplied parameters is not yet supported; refit the model you want to simulate from.",
      newdata = "simulating on new data is not yet supported; the simulation uses the fitted model frame.")
  )
  nsim <- mm_simulate_check_nsim(nsim)
  target <- mm_simulate_target(re.form, use.u, missing(re.form), missing(use.u))
  n <- nobs(object)
  # Prior weights scale the residual variance: sigma^2 / w.
  w <- object$weights %||% rep(1, n)
  resolved <- mm_simulate_resolve_partial(object, target, re.form)
  target <- resolved$target
  out <- mm_with_seed(seed, {
    eta <- if (identical(target, "conditional")) {
      matrix(as.numeric(fitted(object)), n, nsim)
    } else if (identical(target, "partial")) {
      mm_simulate_partial_eta(object, resolved, nsim)
    } else {
      mm_simulate_fixed_eta(object) + mm_simulate_new_re(object, nsim)
    }
    eta + matrix(stats::rnorm(n * nsim), n, nsim) * (object$sigma / sqrt(w))
  })
  out <- as.data.frame(out)
  names(out) <- paste0("sim_", seq_len(nsim))
  rownames(out) <- rownames(object$model_frame)
  attr(out, "seed") <- seed
  attr(out, "mm_method") <- "r_side_gaussian_parametric"
  attr(out, "mm_re_form") <- target
  out
}

mm_random_term_covariance <- function(fit, term_id, basis_labels,
                                      group_label = NULL) {
  p <- length(basis_labels)
  if (!p) return(matrix(0, 0, 0))
  labels <- ifelse(basis_labels == "(Intercept)", "intercept", basis_labels)
  Sigma <- diag(0, p)
  dimnames(Sigma) <- list(basis_labels, basis_labels)

  seen <- stats::setNames(rep(FALSE, p), basis_labels)
  traces <- fit$artifact$covariance_parameter_traces %||% list()
  traces <- traces[vapply(traces, function(x) identical(x$term_id, term_id), logical(1))]
  entries <- unlist(lapply(traces, function(x) x$varcorr_entries %||% list()),
                    recursive = FALSE)
  for (entry in entries) {
    kind <- mm_scalar_text(entry$kind)
    basis <- as.character(unlist(entry$basis %||% list(), use.names = FALSE))
    value <- as.numeric(entry$value %||% NA_real_)
    if (!length(basis) || !is.finite(value)) next
    basis <- ifelse(basis == "intercept", "(Intercept)", basis)
    if (identical(kind, "standard_deviation") && length(basis) == 1L &&
        basis %in% basis_labels) {
      Sigma[basis, basis] <- value^2
      seen[[basis]] <- TRUE
    }
  }
  for (entry in entries) {
    kind <- mm_scalar_text(entry$kind)
    basis <- as.character(unlist(entry$basis %||% list(), use.names = FALSE))
    value <- as.numeric(entry$value %||% NA_real_)
    if (!identical(kind, "correlation") || length(basis) != 2L ||
        !is.finite(value)) next
    basis <- ifelse(basis == "intercept", "(Intercept)", basis)
    if (all(basis %in% basis_labels)) {
      sd1 <- sqrt(Sigma[basis[[1L]], basis[[1L]]])
      sd2 <- sqrt(Sigma[basis[[2L]], basis[[2L]]])
      Sigma[basis[[1L]], basis[[2L]]] <- value * sd1 * sd2
      Sigma[basis[[2L]], basis[[1L]]] <- value * sd1 * sd2
    }
  }

  # Fall back to the VarCorr table only for a standard deviation the traces
  # did not report at all. A reported SD of exactly 0 (a singular fit) is a
  # real value, and the lookup must stay within this term's grouping factor:
  # matching on the coefficient name alone would borrow another term's
  # variance (e.g. the subject intercept for a zero-variance item intercept).
  missing_diag <- !seen
  if (any(missing_diag)) {
    vc <- fit$varcorr$table
    if (!is.null(group_label) && "group" %in% names(vc) &&
        any(vc$group == group_label)) {
      vc <- vc[vc$group == group_label, , drop = FALSE]
    } else if (length(unique(vc$group)) > 1L) {
      vc <- vc[0L, , drop = FALSE]
    }
    for (j in which(missing_diag)) {
      label <- labels[[j]]
      row <- vc[vc$name %in% c(label, basis_labels[[j]]) &
                  vc$variance >= 0, , drop = FALSE]
      if (nrow(row)) Sigma[j, j] <- row$variance[[1L]]
    }
  }
  Sigma
}

mm_rmvnorm <- function(n, Sigma) {
  p <- nrow(Sigma)
  if (!p) return(matrix(numeric(), nrow = n, ncol = 0L))
  Sigma[!is.finite(Sigma)] <- 0
  Sigma <- (Sigma + t(Sigma)) / 2
  root <- tryCatch(chol(Sigma), error = function(cnd) NULL)
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  if (!is.null(root)) {
    out <- z %*% root
  } else {
    eig <- eigen(Sigma, symmetric = TRUE)
    vals <- pmax(eig$values, 0)
    out <- z %*% (eig$vectors %*% diag(sqrt(vals), nrow = p))
  }
  colnames(out) <- colnames(Sigma)
  out
}

mm_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}
