# Leverage and influence measures for mixeff fits (pre-CRAN audit 4.5):
# hatvalues(), cooks.distance() and case/group-deletion influence().

# Response-scale working quantities shared by the hat matrix and Cook's
# distance: the observed response, the fitted mean, the prior weights, the
# working weights of the (penalized) weighted least-squares problem at the
# optimum, and the lme4-style Pearson residuals and dispersion.
mm_influence_parts <- function(model) {
  n <- nobs(model)
  y <- as.numeric(mm_response_vector(model))
  mu <- mm_insample(model, fitted(model))
  prior <- as.numeric(model$weights %||% rep(1, n))
  if (inherits(model, "mm_glmm")) {
    nb <- identical(model$family$family, "negative_binomial")
    # mm_negative_binomial() carries no link functions; take them from R.
    family <- if (nb) {
      stats::make.link(model$family$link %||% "log")
    } else {
      mm_glmm_family_from_info(model$family)
    }
    eta <- family$linkfun(mu)
    var_mu <- if (nb) {
      mu + mu^2 / model$family$nb_theta
    } else {
      family$variance(mu)
    }
    working <- prior * family$mu.eta(eta)^2 / var_mu
    pearson <- (y - mu) * sqrt(prior) / sqrt(var_mu)
    dispersion <- if (model$family$family %in% c("gamma", "inverse_gaussian")) {
      model$dispersion^2
    } else {
      1
    }
  } else {
    working <- prior
    pearson <- (y - mu) * sqrt(prior)
    dispersion <- model$sigma^2
  }
  list(y = y, mu = mu, prior = prior, working = working, pearson = pearson,
       dispersion = dispersion)
}

#' Leverage, Cook's distance and deletion influence for mixeff fits
#'
#' `hatvalues()` returns the diagonal of the hat matrix of the fitted mixed
#' model, `cooks.distance()` the Cook's distances computed from them, and
#' `influence()` refits the model with each case (or each level of a grouping
#' variable) deleted, as `lme4`'s methods of the same names do.
#'
#' The hat matrix is that of the penalized least-squares problem at the
#' optimum, `H = W^{1/2} C M^{-1} C' W^{1/2}` with `C = [X, Z Lambda]` and
#' `M = C' W C + diag(0_p, I_q)`, where `Lambda` is the relative covariance
#' factor and `W` the weights. For an LMM, `W` holds the prior weights and the
#' result equals `lme4::hatvalues()` on the same fit. For a GLMM, `W` holds the
#' working (IRLS) weights at convergence, as [stats::hatvalues()] does for a
#' `glm`; lme4 instead combines the working weights inside the factorization
#' with the prior weights outside it, so the two differ for GLMMs (lme4 warns
#' that its GLMM hat matrix "may not make sense").
#'
#' `cooks.distance()` follows `lme4`: `D_i = (r_i / (1 - h_i))^2 h_i /
#' (phi p)` with Pearson residuals `r_i` (`(y - mu) sqrt(w) / sqrt(V(mu))`),
#' dispersion `phi` (`sigma^2` for an LMM, 1 for binomial, Poisson and
#' negative-binomial GLMMs, `sigma^2` for Gamma) and `p` the rank of the
#' fixed-effect design.
#'
#' `influence()` deletes one case (default) or one level of `groups` at a
#' time and refits with the original settings (`REML`, family, weights,
#' offset, method, control). Deletion refits that fail are recorded as `NA`
#' rows with `converged = FALSE` rather than dropped. The result has class
#' `mm_influence` and the component names of lme4's `influence.merMod`
#' result; [stats::dfbeta()], [stats::dfbetas()] and [stats::cooks.distance()]
#' have methods for it (Cook's distance there is the deletion version based
#' on the change in fixed effects).
#'
#' @param model A fitted `mm_lmm` or `mm_glmm` (or, for the `dfbeta()`,
#'   `dfbetas()` and `cooks.distance()` methods, an `mm_influence` result).
#' @param fullHatMatrix If `TRUE`, return the full `n x n` hat matrix instead
#'   of its diagonal.
#' @param groups For `influence()`: `NULL` (default) to delete single cases,
#'   or the name(s) of grouping column(s) whose levels are deleted one at a
#'   time. Columns are looked up in `data` if given, else in the model frame.
#'   Several names delete the levels of their interaction.
#' @param data Optional data frame with one row per observation (in model
#'   frame order) used only to look up `groups`.
#' @param do.coef If `FALSE`, `influence()` returns only the hat values
#'   without refitting.
#' @param ncores Number of cores for the deletion refits (via
#'   `parallel::mclapply()` on non-Windows platforms); default 1.
#' @param which For `dfbeta()`: `"fixed"` (default) for the fixed effects or
#'   `"var.cov"` for the variance-covariance components.
#' @param ... Unused.
#'
#' @return `hatvalues()` and `cooks.distance()`: a named numeric vector (or a
#'   matrix for `fullHatMatrix = TRUE`). `influence()`: an `mm_influence` list.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rnorm(40), x = rnorm(40),
#'   g = factor(rep(seq_len(8), each = 5))
#' )
#' fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
#' head(hatvalues(fit))
#' head(cooks.distance(fit))
#' infl <- influence(fit, groups = "g")
#' dfbeta(infl)
#' cooks.distance(infl)
#'
#' @name mm_influence
#' @importFrom stats hatvalues cooks.distance influence dfbeta dfbetas
NULL

#' @rdname mm_influence
#' @method hatvalues mm_lmm
#' @export
hatvalues.mm_lmm <- function(model, fullHatMatrix = FALSE, ...) {
  parts <- mm_influence_parts(model)
  X <- stats::model.matrix(model, type = "fixed")
  # getME() returns Z and Lambda as a consistent pair in the engine's
  # (lme4's) term order.
  Z <- getME(model, "Z")
  Lambda <- getME(model, "Lambda")
  if (ncol(Z) != nrow(Lambda)) {
    mm_abort(
      message = "The random-effect design and covariance factor have incompatible sizes; hat values are unavailable for this fit.",
      class = "mm_inference_unavailable",
      reason_code = "hatvalues_design_mismatch",
      input = c(ncol(Z), nrow(Lambda))
    )
  }
  sqrt_w <- sqrt(parts$working)
  A <- cbind(Matrix::Matrix(X * sqrt_w, sparse = TRUE),
             Matrix::Diagonal(x = sqrt_w) %*% (Z %*% Lambda))
  p <- ncol(X)
  q <- ncol(Z)
  M <- Matrix::crossprod(A) +
    Matrix::Diagonal(x = c(rep(0, p), rep(1, q)))
  M <- Matrix::forceSymmetric(M)
  ch <- Matrix::Cholesky(M, LDL = FALSE, perm = TRUE)
  CL <- Matrix::solve(ch, Matrix::solve(ch, Matrix::t(A), system = "P"),
                      system = "L")
  nm <- rownames(model$model_frame)
  if (isTRUE(fullHatMatrix)) {
    out <- as.matrix(Matrix::crossprod(CL))
    dimnames(out) <- list(nm, nm)
    return(out)
  }
  out <- Matrix::colSums(CL^2)
  names(out) <- nm
  out
}

#' @rdname mm_influence
#' @method hatvalues mm_glmm
#' @export
hatvalues.mm_glmm <- hatvalues.mm_lmm

#' @rdname mm_influence
#' @method cooks.distance mm_lmm
#' @export
cooks.distance.mm_lmm <- function(model, ...) {
  parts <- mm_influence_parts(model)
  X <- stats::model.matrix(model, type = "fixed")
  p <- qr(X)$rank
  h <- hatvalues(model)
  res <- (parts$pearson / (1 - h))^2 * h / (parts$dispersion * p)
  res[is.infinite(res)] <- NaN
  names(res) <- rownames(model$model_frame)
  res
}

#' @rdname mm_influence
#' @method cooks.distance mm_glmm
#' @export
cooks.distance.mm_glmm <- cooks.distance.mm_lmm

#' @rdname mm_influence
#' @method influence mm_lmm
#' @export
influence.mm_lmm <- function(model, groups = NULL, data = NULL,
                             do.coef = TRUE, ncores = 1L, ...) {
  if (...length()) {
    mm_abort(
      message = sprintf("`influence()` does not support argument(s) %s.",
                        paste(sprintf("`%s`", names(list(...))), collapse = ", ")),
      class = "mm_arg_error"
    )
  }
  if (!isTRUE(do.coef)) {
    return(structure(list(hat = hatvalues(model)), class = "mm_influence"))
  }
  n <- nobs(model)
  frame <- model$model_frame
  if (is.null(groups)) {
    del <- rownames(frame) %||% as.character(seq_len(n))
    group_name <- "case"
  } else {
    src <- if (is.null(data)) frame else as.data.frame(data)
    if (nrow(src) != n) {
      mm_abort(
        message = sprintf("`data` must have one row per observation (%d).", n),
        class = "mm_arg_error",
        input = nrow(src)
      )
    }
    if (!is.character(groups) || !all(groups %in% names(src))) {
      mm_abort(
        message = "`groups` must name column(s) of `data` (or of the model frame).",
        class = "mm_arg_error",
        input = groups
      )
    }
    del <- if (length(groups) == 1L) {
      as.character(src[[groups]])
    } else {
      do.call(paste, c(unname(as.list(src[groups])), sep = "."))
    }
    group_name <- paste(groups, collapse = ".")
  }
  unique_del <- unique(del)

  fixed <- fixef(model)
  vc <- mm_influence_vc(model)
  one <- function(level) {
    keep <- del != level
    fit <- tryCatch(mm_influence_refit(model, keep), error = function(cnd) cnd)
    if (inherits(fit, "condition")) {
      return(list(fixed = rep(NA_real_, length(fixed)),
                  vc = rep(NA_real_, length(vc)),
                  vcov = matrix(NA_real_, length(fixed), length(fixed),
                                dimnames = list(names(fixed), names(fixed))),
                  converged = FALSE, status = "refit_error"))
    }
    status <- fit_status(fit)
    list(fixed = unname(fixef(fit)), vc = unname(mm_influence_vc(fit)),
         vcov = as.matrix(stats::vcov(fit)),
         converged = startsWith(status, "converged"), status = status)
  }
  ncores <- as.integer(ncores)
  res <- if (ncores > 1L && .Platform$OS.type != "windows" &&
             requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(unique_del, one, mc.cores = ncores)
  } else {
    lapply(unique_del, one)
  }
  fixed_1 <- do.call(rbind, lapply(res, `[[`, "fixed"))
  dimnames(fixed_1) <- list(unique_del, names(fixed))
  vc_1 <- do.call(rbind, lapply(res, `[[`, "vc"))
  dimnames(vc_1) <- list(unique_del, names(vc))
  vcov_1 <- lapply(res, `[[`, "vcov")
  converged <- vapply(res, `[[`, logical(1), "converged")
  status <- vapply(res, `[[`, character(1), "status")
  names(vcov_1) <- names(converged) <- names(status) <- unique_del

  out <- list(fixed, fixed_1, vc, vc_1, as.matrix(stats::vcov(model)), vcov_1,
              group_name, unique_del, converged, status)
  names(out) <- c(
    "fixed.effects", sprintf("fixed.effects[-%s]", group_name),
    "var.cov.comps", sprintf("var.cov.comps[-%s]", group_name),
    "vcov", sprintf("vcov[-%s]", group_name),
    "groups", "deleted", "converged", "fit.status"
  )
  class(out) <- "mm_influence"
  out
}

#' @rdname mm_influence
#' @method influence mm_glmm
#' @export
influence.mm_glmm <- influence.mm_lmm

# Variance components in lme4's influence() layout: sigma^2, then for each
# random term the lower triangle (by column) of its covariance matrix.
mm_influence_vc <- function(fit) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  vc <- if (inherits(fit, "mm_glmm") &&
            !fit$family$family %in% c("gamma", "inverse_gaussian")) 1 else fit$sigma^2
  names(vc) <- "sigma^2"
  labels <- make.unique(vapply(seq_along(terms), function(i) {
    mm_random_term_group_label(fit, terms[[i]], i)
  }, character(1)))
  prefix <- if (length(terms) > 1L) paste0(labels, ":") else rep("", length(terms))
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    term_id <- term$id %||% sprintf("r%d", i - 1L)
    basis <- term$basis %||% list()
    cn <- vapply(basis, mm_basis_label, character(1))
    if (!length(cn)) cn <- "(Intercept)"
    V <- mm_random_term_covariance(fit, term_id, cn,
                                   mm_random_term_group_label(fit, term, i))
    e <- outer(cn, cn, function(a, b) paste0("C[", a, ",", b, "]"))
    diag(e) <- paste0("v[", cn, "]")
    v <- V[lower.tri(V, diag = TRUE)]
    names(v) <- paste0(prefix[[i]], e[lower.tri(e, diag = TRUE)])
    vc <- c(vc, v)
  }
  vc
}

# Refit `model` on the rows in `keep`, with all original settings.
mm_influence_refit <- function(model, keep) {
  data <- model$model_frame[keep, , drop = FALSE]
  control <- model$control %||% mm_control(verbose = -1)
  control$verbose <- -1L
  w <- if (is.null(model$weights)) NULL else model$weights[keep]
  if (inherits(model, "mm_glmm")) {
    do.call(glmm, list(
      formula = model$formula, data = data,
      family = mm_glmm_family_from_info(model$family),
      weights = w,
      offset = if (is.null(model$offset)) NULL else model$offset[keep],
      method = model$method, nAGQ = model$nAGQ,
      inference = "none", control = control
    ))
  } else {
    do.call(lmm, list(
      formula = model$formula, data = data, REML = isTRUE(model$REML),
      weights = w, control = control
    ))
  }
}

#' @rdname mm_influence
#' @method dfbeta mm_influence
#' @export
dfbeta.mm_influence <- function(model, which = c("fixed", "var.cov"), ...) {
  which <- match.arg(which)
  key <- if (identical(which, "fixed")) "fixed.effects" else "var.cov.comps"
  b0 <- model[[key]]
  b <- if (length(model$groups)) {
    model[[sprintf("%s[-%s]", key, model$groups)]]
  }
  if (is.null(b)) {
    mm_abort(
      message = "This influence result has no deletion refits (`do.coef = FALSE`).",
      class = "mm_arg_error"
    )
  }
  b - matrix(b0, nrow = nrow(b), ncol = ncol(b), byrow = TRUE)
}

#' @rdname mm_influence
#' @method dfbetas mm_influence
#' @export
dfbetas.mm_influence <- function(model, ...) {
  vlist <- model[[sprintf("vcov[-%s]", model$groups)]]
  se <- t(vapply(vlist, function(v) sqrt(diag(v)), numeric(nrow(vlist[[1L]]))))
  if (nrow(vlist[[1L]]) == 1L) se <- t(se)
  dfbeta(model) / se
}

#' @rdname mm_influence
#' @method cooks.distance mm_influence
#' @export
cooks.distance.mm_influence <- function(model, ...) {
  db <- dfbeta(model)
  n <- nrow(db)
  p <- ncol(db)
  vinv <- (n - p) / (n * p) * solve(model$vcov)
  out <- vapply(seq_len(n), function(i) {
    as.numeric(db[i, ] %*% vinv %*% db[i, ])
  }, numeric(1))
  names(out) <- rownames(db)
  out
}

#' @export
print.mm_influence <- function(x, ...) {
  if (is.null(x$groups)) {
    cat("<mm_influence> hat values only\n")
    print(utils::head(x$hat))
    return(invisible(x))
  }
  cat(sprintf("<mm_influence> deletion of %d %s level(s); %d refit(s) not converged\n",
              length(x$deleted), x$groups, sum(!x$converged)))
  cat("Change in fixed effects (dfbeta), first rows:\n")
  print(utils::head(dfbeta(x)))
  invisible(x)
}
