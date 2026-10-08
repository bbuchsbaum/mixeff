# Family-aware simulation shared by simulate.mm_lmm() and simulate.mm_glmm(),
# and refit() for GLMMs.
#
# Semantics follow lme4::simulate.merMod():
#   * re.form = NA (default) or ~0, or use.u = FALSE: draw NEW random effects
#     from the fitted covariance (unconditional simulation);
#   * re.form = NULL, or use.u = TRUE: condition on the fitted conditional
#     modes (BLUPs), i.e. simulate only the response noise around fitted().
# The random-number stream is consumed in lme4's order (all spherical
# random-effect draws for all nsim, terms ordered as lme4 orders them, levels
# within terms, then the response draws), so with matching estimates a seeded
# mixeff simulation reproduces lme4's draws.

# Resolve re.form / use.u into "unconditional" or "conditional".
mm_simulate_target <- function(re.form, use.u, re_form_missing, use_u_missing) {
  if (!use_u_missing) {
    if (!re_form_missing) {
      mm_abort(
        message = "Specify only one of `use.u` and `re.form`.",
        class = "mm_arg_error",
        input = list(use.u = use.u, re.form = re.form)
      )
    }
    if (!is.logical(use.u) || length(use.u) != 1L || is.na(use.u)) {
      mm_abort(
        message = "`use.u` must be TRUE or FALSE.",
        class = "mm_arg_error",
        input = use.u
      )
    }
    return(if (use.u) "conditional" else "unconditional")
  }
  target <- mm_prediction_target(re.form)
  if (identical(target, "population")) return("unconditional")
  if (identical(target, "conditional")) return("conditional")
  mm_abort(
    message = paste(
      "`re.form` must be NA or ~0 (simulate new random effects, the lme4",
      "default) or NULL (condition on the fitted random effects); partial",
      "random-effect formulas are not available for simulation."
    ),
    class = "mm_inference_unavailable",
    reason_code = "simulate_partial_re_form_unavailable",
    input = re.form
  )
}

mm_simulate_check_nsim <- function(nsim) {
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    mm_abort(
      message = "`nsim` must be a positive integer.",
      class = "mm_arg_error",
      input = nsim
    )
  }
  as.integer(nsim)
}

# Lower-triangular factor L with L %*% t(L) == Sigma that tolerates positive
# semi-definite (singular) covariance blocks: a zero pivot yields a zero
# column, as lme4's relative factor Lambda does for a boundary fit.
mm_chol_lower_psd <- function(Sigma) {
  p <- nrow(Sigma)
  L <- matrix(0, p, p, dimnames = dimnames(Sigma))
  if (!p) return(L)
  Sigma[!is.finite(Sigma)] <- 0
  Sigma <- (Sigma + t(Sigma)) / 2
  tol <- max(abs(diag(Sigma)), 0) * 1e-12
  for (j in seq_len(p)) {
    s <- Sigma[j, j] - sum(L[j, seq_len(j - 1L)]^2)
    if (s <= tol) next
    L[j, j] <- sqrt(s)
    if (j < p) {
      for (i in seq.int(j + 1L, p)) {
        L[i, j] <- (Sigma[i, j] - sum(L[i, seq_len(j - 1L)] * L[j, seq_len(j - 1L)])) /
          L[j, j]
      }
    }
  }
  L
}

# Random-term layout for simulation: per term the grouping factor, basis
# values, the lower factor of its absolute covariance, in lme4's term order.
mm_simulate_term_layout <- function(fit) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  n <- nobs(fit)
  out <- lapply(seq_along(terms), function(i) {
    term <- terms[[i]]
    term_id <- term$id %||% sprintf("r%d", i - 1L)
    group_label <- mm_random_term_group_label(fit, term, i)
    group <- mm_group_factor(fit$model_frame, group_label)
    basis <- term$basis %||% list()
    basis_labels <- vapply(basis, mm_basis_label, character(1))
    basis_values <- lapply(basis, mm_basis_values, frame = fit$model_frame)
    if (!length(basis_values)) {
      basis_labels <- "(Intercept)"
      basis_values <- list(rep(1, n))
    }
    Sigma <- mm_random_term_covariance(fit, term_id, basis_labels, group_label)
    list(
      group = group,
      nlevels = nlevels(group),
      basis = basis_values,
      L = mm_chol_lower_psd(Sigma)
    )
  })
  nl <- vapply(out, `[[`, integer(1), "nlevels")
  # lme4::mkReTrms orders terms by decreasing number of levels, using
  # rev(order(nl)) only when the formula order is not already non-increasing.
  if (length(nl) > 1L && any(diff(nl) > 0)) {
    out <- out[rev(order(nl))]
  }
  out
}

# n x nsim matrix of simulated random-effect contributions (Z b) on the
# linear-predictor scale.
mm_simulate_new_re <- function(fit, nsim) {
  n <- nobs(fit)
  layout <- mm_simulate_term_layout(fit)
  if (!length(layout)) return(matrix(0, n, nsim))
  q_term <- vapply(layout, function(t) t$nlevels * length(t$basis), numeric(1))
  q <- sum(q_term)
  u <- matrix(stats::rnorm(q * nsim), nrow = q, ncol = nsim)
  out <- matrix(0, n, nsim)
  start <- 0
  for (t in layout) {
    k <- t$nlevels
    p <- length(t$basis)
    rows <- start + seq_len(k * p)
    start <- start + k * p
    idx <- as.integer(t$group)
    for (s in seq_len(nsim)) {
      # lme4 lays out u level-major within a term (level 1: coef 1..p, ...).
      U <- matrix(u[rows, s], nrow = k, ncol = p, byrow = TRUE)
      b <- U %*% t(t$L)
      for (j in seq_len(p)) {
        out[, s] <- out[, s] + t$basis[[j]] * b[idx, j]
      }
    }
  }
  out
}

# Linear predictor without random effects (fixed part + offset).
mm_simulate_fixed_eta <- function(fit) {
  if (inherits(fit, "mm_glmm")) {
    as.numeric(stats::predict(fit, re.form = NA, type = "link"))
  } else {
    as.numeric(fit$fixed_fitted)
  }
}

# Is the binomial response of a fit a cbind(successes, failures) matrix?
mm_glmm_cbind_colnames <- function(fit) {
  if (!identical(mm_response_name(fit), ".mm_binomial_response")) return(NULL)
  call_formula <- fit$call$formula
  if (is.name(call_formula)) {
    call_formula <- tryCatch(eval(call_formula, environment(fit$formula)),
                             error = function(cnd) NULL)
  }
  lhs <- if (is.call(call_formula) && length(call_formula) == 3L) {
    call_formula[[2L]]
  }
  if (is.call(lhs) && identical(as.character(lhs[[1L]]), "cbind") &&
      length(lhs) == 3L) {
    # model.response() of cbind(a, b - a) names only symbol arguments.
    vapply(as.list(lhs)[2:3], function(e) {
      if (is.name(e)) as.character(e) else ""
    }, character(1))
  } else {
    c("", "")
  }
}

#' @rdname simulate.mm_lmm
#' @method simulate mm_glmm
#' @export
simulate.mm_glmm <- function(object, nsim = 1, seed = NULL, use.u = FALSE,
                             re.form = NA, ...) {
  mm_reject_unsupported_dots(
    list(...), "simulate",
    c(newparams = "simulating from user-supplied parameters is not yet supported; refit the model you want to simulate from.",
      newdata = "simulating on new data is not yet supported; the simulation uses the fitted model frame.",
      family = "the simulation family is the fitted family; refit with glmm(family = ) to simulate from another.",
      weights = "the simulation uses the fitted prior weights.")
  )
  nsim <- mm_simulate_check_nsim(nsim)
  target <- mm_simulate_target(re.form, use.u, missing(re.form), missing(use.u))
  fam_info <- object$family
  family <- mm_glmm_family_from_info(fam_info)
  n <- nobs(object)
  wts <- object$weights %||% rep(1, n)
  cbind_names <- if (identical(fam_info$family, "binomial")) {
    mm_glmm_cbind_colnames(object)
  }
  if (identical(fam_info$family, "binomial") && any(abs(wts %% 1) > 1e-8)) {
    mm_abort(
      message = "Cannot simulate binomial responses from non-integer prior weights (trial counts).",
      class = "mm_inference_unavailable",
      reason_code = "simulate_noninteger_binomial_weights",
      input = wts
    )
  }

  out <- mm_with_seed(seed, {
    mu <- if (identical(target, "conditional")) {
      matrix(as.numeric(fitted(object)), n, nsim)
    } else {
      eta <- mm_simulate_fixed_eta(object) + mm_simulate_new_re(object, nsim)
      family$linkinv(eta)
    }
    ftd <- as.numeric(mu)
    ntot <- n * nsim
    sims <- switch(
      fam_info$family,
      binomial = if (!is.null(cbind_names)) {
        lapply(seq_len(nsim), function(s) {
          y <- stats::rbinom(n, size = wts, prob = mu[, s])
          yy <- cbind(y, wts - y)
          colnames(yy) <- cbind_names
          yy
        })
      } else {
        split(stats::rbinom(ntot, size = wts, prob = ftd) / wts, gl(nsim, n))
      },
      poisson = split(stats::rpois(ntot, ftd), gl(nsim, n)),
      negative_binomial = split(
        stats::rnbinom(ntot, mu = ftd, size = fam_info$nb_theta),
        gl(nsim, n)
      ),
      gamma = {
        # Var(y) = phi * mu^2 with phi = sigma^2; prior weights scale the
        # shape (phi / w), as glm's Gamma()$simulate does.
        shape <- wts / object$dispersion^2
        split(stats::rgamma(ntot, shape = shape, rate = shape / ftd),
              gl(nsim, n))
      },
      mm_abort(
        message = sprintf("Simulation is not implemented for GLMM family `%s`.",
                          fam_info$family),
        class = "mm_inference_unavailable",
        reason_code = "simulate_family_unavailable",
        input = fam_info$family
      )
    )
    sims <- lapply(sims, function(x) {
      if (!is.matrix(x)) x <- as.numeric(x)
      x
    })
    names(sims) <- paste0("sim_", seq_len(nsim))
    structure(sims, class = "data.frame",
              row.names = rownames(object$model_frame) %||% seq_len(n))
  })
  attr(out, "seed") <- seed
  attr(out, "mm_method") <- "r_side_glmm_parametric"
  attr(out, "mm_re_form") <- target
  out
}

#' @rdname refit
#' @method refit mm_glmm
#' @export
refit.mm_glmm <- function(object, newresp, ...) {
  mm_reject_unsupported_dots(
    list(...)[setdiff(names(list(...)), "control")], "refit",
    c(newweights = "changing prior weights on refit is not supported; refit with glmm(..., weights = ).")
  )
  if (missing(newresp)) {
    mm_abort(
      message = "`newresp` is required: supply the new response vector.",
      class = "mm_arg_error"
    )
  }
  n <- nobs(object)
  weights <- object$weights
  newresp <- mm_refit_unwrap_resp(newresp)
  if (identical(object$family$family, "binomial")) {
    if (is.matrix(newresp) && ncol(newresp) == 2L) {
      trials <- rowSums(newresp)
      if (nrow(newresp) != n || anyNA(trials) || any(trials <= 0)) {
        mm_abort(
          message = "A two-column `newresp` must hold (successes, failures) with positive totals, one row per observation.",
          class = "mm_arg_error",
          input = newresp
        )
      }
      newresp <- newresp[, 1L] / trials
      weights <- as.numeric(trials)
    } else if (is.factor(newresp)) {
      newresp <- as.integer(newresp != levels(newresp)[[1L]])
    } else if (is.logical(newresp)) {
      newresp <- as.integer(newresp)
    }
  }
  if (!is.numeric(newresp) || is.matrix(newresp) && ncol(newresp) != 1L ||
      length(newresp) != n || anyNA(newresp)) {
    mm_abort(
      message = "`newresp` must be a numeric vector with one value per observation and no missing values.",
      class = "mm_arg_error",
      input = newresp
    )
  }
  if (identical(object$family$family, "binomial")) {
    bernoulli <- is.null(weights)
    if (bernoulli && any(!newresp %in% c(0, 1))) {
      mm_abort(
        message = "This binomial fit has a 0/1 response; `newresp` must contain only 0 and 1 (or be a two-level factor, or a (successes, failures) matrix).",
        class = "mm_arg_error",
        input = newresp
      )
    }
    if (!bernoulli && any(newresp < 0 | newresp > 1)) {
      mm_abort(
        message = "This binomial fit has a proportion response; `newresp` must lie in [0, 1].",
        class = "mm_arg_error",
        input = newresp
      )
    }
  }
  data <- object$model_frame
  data[[mm_response_name(object)]] <- as.numeric(newresp)
  control <- list(...)$control
  if (is.null(control)) {
    control <- object$control %||% mm_control(verbose = -1)
    control$verbose <- -1L
  }
  # do.call() inlines the values: glmm() evaluates `weights`/`offset` in
  # `data` first, so passing symbols could pick up same-named columns.
  fit <- do.call(glmm, list(
    formula = object$formula, data = data,
    family = mm_glmm_family_from_info(object$family),
    weights = weights, offset = object$offset,
    method = object$method, nAGQ = object$nAGQ,
    inference = object$inference_request %||% "auto",
    control = control
  ))
  fit$refit <- list(
    source = "refit",
    original_fit_status = fit_status(object)
  )
  fit
}

# lme4's refit() accepts a one-column data frame / list, as returned by
# simulate(nsim = 1).
mm_refit_unwrap_resp <- function(newresp) {
  if (is.list(newresp)) {
    if (length(newresp) != 1L) {
      mm_abort(
        message = "`newresp` given as a list or data frame must have exactly one column.",
        class = "mm_arg_error",
        input = newresp
      )
    }
    newresp <- newresp[[1L]]
  }
  newresp
}
