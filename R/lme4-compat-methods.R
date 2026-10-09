# lme4-shaped extractor surface: residual types, family(), weights(),
# REMLcrit(), the penalized-least-squares pieces behind getME(), the
# lme4-compatible VarCorr() object, and multi-model AIC()/BIC().
#
# Everything here is recomputed R-side from quantities the fit stores
# (beta, theta, the theta map, conditional modes, fitted values, the model
# frame); nothing calls the engine.

# ---- families ---------------------------------------------------------------

#' Family objects and prior weights for mixeff fits
#'
#' `family()` returns the R [family()] object of a fit: [gaussian()] for
#' [lmm()] fits and the GLMM family (with its link) for [glmm()] fits. For
#' negative-binomial fits it is `MASS::negative.binomial(theta)` at the
#' fitted (or fixed) theta, like `lme4::glmer.nb()`; when MASS is not
#' installed an equivalent family object with the same variance, deviance and
#' link functions is returned.
#'
#' `weights()` mirrors `lme4:::weights.merMod()`: `type = "prior"` returns the
#' prior weights, a vector of ones when the model was fitted without weights;
#' `type = "working"` returns the final PIRLS working weights of a GLMM.
#'
#' @param object A fitted `mm_lmm` or `mm_glmm` object.
#' @param type `"prior"` or `"working"` (GLMMs only).
#' @param ... Unused.
#' @return A `family` object, or a numeric vector of weights.
#' @examples
#' set.seed(1)
#' df <- data.frame(y = rnorm(60), x = rnorm(60),
#'                  g = factor(rep(seq_len(10), each = 6)))
#' fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
#' family(fit)
#' head(weights(fit))
#' @name mm_family
NULL

#' @rdname mm_family
#' @importFrom stats family
#' @export
family.mm_lmm <- function(object, ...) {
  stats::gaussian()
}

#' @rdname mm_family
#' @export
family.mm_glmm <- function(object, ...) {
  mm_glmm_family_object(object)
}

mm_glmm_family_object <- function(fit) {
  info <- fit$family %||% list()
  fam <- as.character(info$family %||% NA_character_)
  link <- as.character(info$link %||% NA_character_)
  switch(
    fam,
    binomial = stats::binomial(link = link),
    bernoulli = stats::binomial(link = link),
    poisson = stats::poisson(link = link),
    gamma = stats::Gamma(link = link),
    Gamma = stats::Gamma(link = link),
    negative_binomial = mm_nb_family(as.numeric(info$nb_theta), link),
    mm_abort(
      message = sprintf("No R family object is known for GLMM family `%s`.",
                        fam),
      class = "mm_inference_unavailable",
      reason_code = "glmm_family_object_unavailable"
    )
  )
}

mm_nb_family <- function(theta, link = "log") {
  if (!length(theta) || !is.finite(theta) || theta <= 0) {
    mm_abort(
      message = "The negative-binomial fit does not carry a finite theta.",
      class = "mm_inference_unavailable",
      reason_code = "nb_theta_unavailable"
    )
  }
  if (requireNamespace("MASS", quietly = TRUE)) {
    return(MASS::negative.binomial(theta = theta, link = link))
  }
  lk <- stats::make.link(link)
  structure(
    list(
      family = sprintf("Negative Binomial(%s)", format(round(theta, 4))),
      link = link,
      linkfun = lk$linkfun,
      linkinv = lk$linkinv,
      variance = function(mu) mu + mu^2 / theta,
      dev.resids = function(y, mu, wt) {
        2 * wt * (y * log(pmax(1, y) / mu) -
                    (y + theta) * log((y + theta) / (mu + theta)))
      },
      aic = function(y, n, mu, wt, dev) {
        term <- (y + theta) * log(mu + theta) - y * log(mu) +
          lgamma(y + 1) - theta * log(theta) + lgamma(theta) -
          lgamma(theta + y)
        2 * sum(term * wt)
      },
      mu.eta = lk$mu.eta,
      validmu = function(mu) all(mu > 0),
      valideta = lk$valideta,
      theta = theta
    ),
    class = "family"
  )
}

# Prior weights as lme4 stores them: ones when the fit is unweighted.
mm_prior_weights <- function(fit) {
  w <- fit$weights
  if (is.null(w)) rep(1, nobs(fit)) else as.numeric(w)
}

mm_glmm_offset <- function(fit) {
  off <- fit$offset
  if (is.null(off)) rep(0, nobs(fit)) else as.numeric(off)
}

# Response, conditional mean, prior weights and family of a GLMM, on lme4's
# scale (binomial responses are proportions with trial-count weights).
mm_glmm_resp <- function(fit) {
  fam <- mm_glmm_family_object(fit)
  y <- as.numeric(mm_response_vector(fit))
  mu <- as.numeric(fit$fitted)
  list(y = y, mu = mu, wt = mm_prior_weights(fit), family = fam,
       eta = fam$linkfun(mu))
}

# PIRLS working weights at the converged conditional modes:
# prior weight * mu.eta(eta)^2 / variance(mu) (lme4's Xwts^2).
mm_glmm_working_weights <- function(fit) {
  r <- mm_glmm_resp(fit)
  r$wt * r$family$mu.eta(r$eta)^2 / r$family$variance(r$mu)
}

#' @rdname mm_family
#' @importFrom stats weights
#' @export
weights.mm_lmm <- function(object, type = c("prior", "working"), ...) {
  type <- match.arg(type)
  if (identical(type, "working")) {
    mm_abort(
      message = "Working weights are only available for GLMMs (as in lme4).",
      class = "mm_arg_error",
      input = type
    )
  }
  mm_prior_weights(object)
}

#' @rdname mm_family
#' @export
weights.mm_glmm <- function(object, type = c("prior", "working"), ...) {
  type <- match.arg(type)
  if (identical(type, "working")) {
    return(mm_glmm_working_weights(object))
  }
  mm_prior_weights(object)
}

# ---- residuals ----------------------------------------------------------------

mm_glmm_residuals <- function(fit, type) {
  r <- mm_glmm_resp(fit)
  switch(
    type,
    deviance = {
      d <- sqrt(pmax(r$family$dev.resids(r$y, r$mu, r$wt), 0))
      ifelse(r$y > r$mu, d, -d)
    },
    pearson = (r$y - r$mu) * sqrt(r$wt / r$family$variance(r$mu)),
    working = (r$y - r$mu) / r$family$mu.eta(r$eta),
    response = r$y - r$mu
  )
}

# ---- REMLcrit -------------------------------------------------------------

#' REML criterion of a linear mixed model
#'
#' `REMLcrit()` mirrors `lme4::REMLcrit()`: the REML criterion (`-2` times the
#' restricted log-likelihood) at the fitted parameters. For a REML fit this
#' is `-2 * logLik(fit)`. For an ML fit it is the REML criterion evaluated at
#' the ML estimates of the covariance parameters, computed from the
#' penalized least-squares decomposition (`getME(fit, "devcomp")`), exactly
#' as lme4 does. `lme4::REMLcrit()` is a plain function rather than a
#' generic, so mixeff exports its own; it forwards `merMod` objects to lme4.
#'
#' @param object A fitted `mm_lmm` (or an lme4 `merMod`).
#' @return A single number.
#' @examples
#' set.seed(1)
#' df <- data.frame(y = rnorm(60), x = rnorm(60),
#'                  g = factor(rep(seq_len(10), each = 6)))
#' fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
#' REMLcrit(fit)
#' @export
REMLcrit <- function(object) {
  if (inherits(object, "merMod") && requireNamespace("lme4", quietly = TRUE)) {
    return(lme4::REMLcrit(object))
  }
  if (inherits(object, "mm_glmm")) {
    mm_abort(
      message = "The REML criterion is only defined for linear mixed models (as in lme4).",
      class = "mm_inference_unavailable",
      reason_code = "remlcrit_glmm_unavailable"
    )
  }
  if (!inherits(object, "mm_lmm")) {
    mm_abort(
      message = "`REMLcrit()` expects a fitted mixeff LMM.",
      class = "mm_arg_error",
      input = object
    )
  }
  if (isTRUE(object$REML)) {
    return(-2 * as.numeric(object$logLik))
  }
  dc <- mm_devcomp(object)
  cmp <- dc$cmp
  n <- dc$dims[["n"]]
  nmp <- n - dc$dims[["p"]]
  ldW <- sum(log(mm_prior_weights(object)))
  unname(-ldW + cmp[["ldL2"]] + cmp[["ldRX2"]] +
           nmp * (1 + log(2 * pi * cmp[["pwrss"]]) - log(nmp)))
}

# ---- random-effect structure from the theta map ----------------------------

# One record per random-effect term in lme4's mkReTrms order: grouping
# label, lme4-named basis columns, the level order used by ranef(), the
# relative-covariance template (theta index per Lambda cell, indexing the
# lme4-ordered theta `$theta`) and theta lower bounds. Built from the stored
# theta map, so Lambda/Z/u/b are mutually consistent with fit$theta.
#
# The engine (like MixedModels.jl) orders its terms by decreasing number of
# random effects and keeps a `||` term as one diagonal block; lme4 splits
# `||` into scalar terms (the engine's theta map already has one entry per
# such piece, `term_id` = its formula position) and reorders the formula's
# terms by decreasing number of levels only when some later term has more
# levels than an earlier one (`rev(order(nl))`, so ties come out reversed).
# `$perm` maps lme4 theta positions to engine positions:
# `$theta == fit$theta[$perm]`.
mm_re_structure <- function(fit) {
  .mm_lazy(fit, "re_structure", mm_compute_re_structure)
}

mm_compute_re_structure <- function(fit) {
  maps <- fit$artifact$theta_maps %||% list()
  theta <- as.numeric(fit$theta)
  if (!length(maps)) {
    mm_abort(
      message = "The fit does not carry a theta map; the random-effect structure cannot be rebuilt.",
      class = "mm_inference_unavailable",
      reason_code = "theta_map_unavailable"
    )
  }
  re <- fit$random_effects %||% list()
  frame <- fit$model_frame
  terms <- lapply(seq_along(maps), function(k) {
    map <- maps[[k]]$map %||% maps[[k]]
    group <- mm_group_lme4_label(map$group, fallback = as.character(map$group))
    # Nested groups the R layer relabelled (`a/h` -> lme4's `h:a`).
    gm <- fit$group_map %||% list()
    if (length(gm)) {
      keys <- vapply(names(gm), mm_group_lme4_label, character(1),
                     fallback = "")
      hit <- match(group, c(keys, names(gm)))
      if (!is.na(hit)) {
        group <- as.character(gm[[(hit - 1L) %% length(gm) + 1L]]$label %||% group)
      }
    }
    raw_basis <- as.character(unlist(map$user_basis %||% map$optimizer_basis,
                                     use.names = FALSE))
    if (!length(raw_basis)) raw_basis <- "intercept"
    is_int <- raw_basis %in% c("intercept", "(Intercept)", "1")
    cnames <- raw_basis
    cnames[is_int] <- "(Intercept)"
    if (any(!is_int)) {
      cnames[!is_int] <- mm_re_colnames_lme4(raw_basis[!is_int], frame)
    }
    slots <- map$theta_slots %||% list()
    rows <- vapply(slots, function(sl) {
      as.integer(sl$lambda_row %||% sl$basis_row)
    }, integer(1))
    cols <- vapply(slots, function(sl) {
      as.integer(sl$lambda_col %||% sl$basis_col)
    }, integer(1))
    # A `||` term split into scalar pieces reports Lambda positions relative
    # to the original block; re-base them on this term.
    base <- if (length(rows)) min(c(rows, cols)) else 0L
    basis_matrix <- NULL
    if (length(rows) && max(c(rows, cols)) - base + 1L > length(cnames)) {
      # A factor in the basis spans several columns; the engine codes them
      # like R's model.matrix() (lme4's coding), so rebuild them that way.
      basis_matrix <- mm_re_factor_basis_matrix(raw_basis, is_int, frame)
      if (!is.null(basis_matrix)) cnames <- colnames(basis_matrix)
    }
    p <- length(cnames)
    Tidx <- matrix(0L, p, p)
    gidx <- integer(length(slots))
    lower <- numeric(length(slots))
    if (length(rows) && max(c(rows, cols)) - base >= p) {
      mm_abort(
        message = "The theta map's Lambda positions do not fit the term's basis.",
        class = "mm_schema_error"
      )
    }
    for (s in seq_along(slots)) {
      sl <- slots[[s]]
      i <- rows[[s]] - base + 1L
      j <- cols[[s]] - base + 1L
      g <- as.integer(sl$global_index) + 1L
      Tidx[i, j] <- g
      gidx[[s]] <- g
      cons <- sl$constraint
      lower[[s]] <- if (is.list(cons) && !is.null(cons$lower_bound)) {
        as.numeric(cons$lower_bound$lower %||% 0)
      } else {
        -Inf
      }
    }
    if (any(gidx > length(theta)) || any(gidx < 1L)) {
      mm_abort(
        message = "The theta map does not match the fitted theta vector.",
        class = "mm_schema_error"
      )
    }
    gf <- mm_group_factor(frame, group)
    levels <- if (!is.null(re[[group]])) rownames(re[[group]]) else levels(gf)
    term_id <- as.character(map$term_id %||% "")
    semantic <- suppressWarnings(as.integer(sub("^r", "", term_id)))
    list(group = group, raw_basis = raw_basis, cnames = cnames, p = p,
         basis_matrix = basis_matrix, Tidx = Tidx, theta_index = gidx, lower = lower,
         factor = factor(as.character(gf), levels = levels),
         levels = levels,
         semantic = if (is.na(semantic)) k - 1L else semantic)
  })
  # lme4 term order: formula order, then mkReTrms' reordering.
  terms <- terms[order(vapply(terms, `[[`, integer(1), "semantic"))]
  nl <- vapply(terms, function(tm) length(tm$levels), integer(1))
  if (length(nl) > 1L && any(diff(nl) > 0)) {
    terms <- terms[rev(order(nl))]
  }
  # Renumber theta in lme4 order: term by term, lower-triangle cells in
  # column-major order within a term (lme4's theta layout).
  perm <- integer()
  lower <- numeric()
  for (k in seq_along(terms)) {
    tm <- terms[[k]]
    cells <- which(tm$Tidx > 0L, arr.ind = TRUE)
    cells <- cells[order(cells[, 2L], cells[, 1L]), , drop = FALSE]
    engine_idx <- tm$Tidx[cells]
    new_idx <- length(perm) + seq_along(engine_idx)
    tm$lower <- tm$lower[match(engine_idx, tm$theta_index)]
    tm$Tidx[cells] <- new_idx
    tm$engine_theta_index <- engine_idx
    tm$theta_index <- new_idx
    perm <- c(perm, engine_idx)
    lower <- c(lower, tm$lower)
    terms[[k]] <- tm
  }
  if (length(perm) != length(theta) || anyDuplicated(perm)) {
    mm_abort(
      message = "The theta map does not cover the fitted theta vector exactly once.",
      class = "mm_schema_error"
    )
  }
  list(terms = terms, lower = lower, perm = perm, theta = theta[perm])
}

# R model.matrix() columns of a random-effect basis given as variable labels
# (`(0 + f | g)` -> `fp`, `fq`, `fr`; `(1 + f | g)` -> intercept + contrasts).
mm_re_factor_basis_matrix <- function(raw_basis, is_int, frame) {
  vars <- raw_basis[!is_int]
  if (!length(vars)) return(NULL)
  labels <- vapply(vars, function(v) {
    if (v %in% names(frame)) sprintf("`%s`", v) else v
  }, character(1))
  fo <- tryCatch(
    stats::as.formula(paste("~", paste(labels, collapse = " + "),
                            if (any(is_int)) "" else "- 1")),
    error = function(cnd) NULL
  )
  if (is.null(fo)) return(NULL)
  mm <- tryCatch(stats::model.matrix(fo, data = frame),
                 error = function(cnd) NULL)
  if (is.null(mm)) return(NULL)
  attr(mm, "assign") <- NULL
  attr(mm, "contrasts") <- NULL
  mm
}

# Columns of a term's basis as a list of numeric vectors.
mm_re_term_values <- function(tm, frame) {
  if (!is.null(tm$basis_matrix)) {
    return(lapply(seq_len(ncol(tm$basis_matrix)),
                  function(j) as.numeric(tm$basis_matrix[, j])))
  }
  lapply(tm$raw_basis, mm_re_basis_values, frame = frame)
}

mm_re_basis_values <- function(label, frame) {
  if (label %in% c("intercept", "(Intercept)", "1")) {
    return(rep(1, nrow(frame)))
  }
  if (label %in% names(frame)) {
    v <- frame[[label]]
    if (is.numeric(v) || is.logical(v)) return(as.numeric(v))
  }
  m <- regmatches(label, regexec("^(.*): (.*)$", label))[[1L]]
  if (length(m) == 3L && m[[2L]] %in% names(frame)) {
    return(as.numeric(as.character(frame[[m[[2L]]]]) == m[[3L]]))
  }
  parts <- strsplit(label, ":", fixed = TRUE)[[1L]]
  if (length(parts) > 1L && all(parts %in% names(frame)) &&
      all(vapply(frame[parts], function(v) is.numeric(v) || is.logical(v),
                 logical(1)))) {
    return(Reduce(`*`, lapply(frame[parts], as.numeric)))
  }
  mm_abort(
    message = sprintf(
      "Cannot rebuild random-effect basis column `%s` from the stored model frame.",
      label
    ),
    class = "mm_inference_unavailable",
    reason_code = "random_basis_rebuild_unavailable",
    input = label
  )
}

# Per-term transposed model matrices (lme4's Ztlist is one per term; within
# a term, rows are level-major then basis column, as in lme4's Zt).
mm_re_term_zt <- function(fit, tm) {
  frame <- fit$model_frame
  n <- nrow(frame)
  vals <- mm_re_term_values(tm, frame)
  lev <- as.integer(tm$factor)
  nl <- length(tm$levels)
  ok <- !is.na(lev)
  i <- unlist(lapply(seq_len(tm$p), function(j) (lev[ok] - 1L) * tm$p + j))
  jj <- rep(which(ok), tm$p)
  x <- unlist(lapply(vals, function(v) v[ok]))
  rn <- paste(rep(tm$levels, each = tm$p), rep(tm$cnames, nl), sep = ".")
  if (tm$p == 1L) rn <- tm$levels
  Matrix::sparseMatrix(i = i, j = jj, x = x, dims = c(nl * tm$p, n),
                       dimnames = list(rn, rownames(frame)))
}

mm_re_zt <- function(fit) {
  .mm_lazy(fit, "Zt_engine", function(fit) {
    st <- mm_re_structure(fit)
    do.call(rbind, lapply(st$terms, mm_re_term_zt, fit = fit))
  })
}

# Lambdat with lme4's sparsity pattern; `Lind` maps each stored entry to its
# theta index.
mm_re_lambdat <- function(fit) {
  .mm_lazy(fit, "Lambdat_engine", function(fit) {
    st <- mm_re_structure(fit)
    ii <- integer()
    jj <- integer()
    xx <- integer()
    offset <- 0L
    for (tm in st$terms) {
      nz <- which(tm$Tidx > 0L, arr.ind = TRUE)  # rows i, cols j of T
      for (l in seq_along(tm$levels)) {
        base <- offset + (l - 1L) * tm$p
        # Lambdat = kronecker(I, t(T)): entry (j, i) holds T[i, j]
        ii <- c(ii, base + nz[, 2L])
        jj <- c(jj, base + nz[, 1L])
        xx <- c(xx, tm$Tidx[nz])
      }
      offset <- offset + length(tm$levels) * tm$p
    }
    tmpl <- Matrix::sparseMatrix(i = ii, j = jj, x = as.numeric(xx),
                                 dims = c(offset, offset))
    Lind <- as.integer(tmpl@x)
    Lt <- tmpl
    Lt@x <- st$theta[Lind]
    list(Lambdat = Lt, Lind = Lind)
  })
}

mm_pinv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  keep <- s$d > tol * max(c(s$d, 0))
  if (!any(keep)) return(matrix(0, ncol(A), nrow(A)))
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# Conditional modes stacked in Z column order (b) and the spherical modes u
# (b = Lambda u; u is the minimum-norm solution, which is the penalized
# optimum when Lambda is singular).
mm_re_bu <- function(fit) {
  .mm_lazy(fit, "re_bu", function(fit) {
    st <- mm_re_structure(fit)
    theta <- st$theta
    re <- fit$random_effects
    b <- list()
    u <- list()
    for (tm in st$terms) {
      df <- re[[tm$group]]
      if (is.null(df) || !all(tm$cnames %in% names(df))) {
        mm_abort(
          message = sprintf(
            "Conditional modes for term `%s` do not match the theta map.",
            tm$group
          ),
          class = "mm_schema_error"
        )
      }
      B <- as.matrix(df[tm$levels, tm$cnames, drop = FALSE])  # levels x p
      Tm <- matrix(0, tm$p, tm$p)
      Tm[tm$Tidx > 0] <- theta[tm$Tidx[tm$Tidx > 0]]
      U <- B %*% t(mm_pinv(Tm))
      b[[length(b) + 1L]] <- as.vector(t(B))
      u[[length(u) + 1L]] <- as.vector(t(U))
    }
    list(b = unlist(b), u = unlist(u))
  })
}

# Penalized least-squares decomposition at the fitted parameters, with LMM
# prior weights or GLMM PIRLS working weights: the sparse Cholesky factor L
# of Lambda'Z'WZ Lambda + I (same fill-reducing permutation as lme4), RZX,
# RX, and lme4's devcomp.
mm_pls <- function(fit) {
  .mm_lazy(fit, "pls", mm_compute_pls)
}

mm_compute_pls <- function(fit) {
  glmm <- inherits(fit, "mm_glmm")
  X <- stats::model.matrix(fit, type = "fixed")
  X <- X[, names(fit$beta), drop = FALSE]
  Zt <- mm_re_zt(fit)
  Lt <- mm_re_lambdat(fit)$Lambdat
  w <- if (glmm) mm_glmm_working_weights(fit) else mm_prior_weights(fit)
  sw <- sqrt(w)
  LtZt <- Lt %*% Zt                       # q x n
  LtZtW <- LtZt %*% Matrix::Diagonal(x = sw)
  L <- Matrix::Cholesky(Matrix::tcrossprod(LtZtW), LDL = FALSE, Imult = 1)
  XW <- X * sw
  RZX <- Matrix::solve(L, Matrix::solve(L, LtZtW %*% XW, system = "P"),
                       system = "L")
  RZX <- as.matrix(RZX)
  XtX <- crossprod(XW)
  RX <- chol(XtX - crossprod(RZX))
  dimnames(RX) <- list(colnames(X), colnames(X))
  bu <- mm_re_bu(fit)
  ussq <- sum(bu$u^2)
  n <- nrow(X)
  p <- ncol(X)
  ldL2 <- 2 * as.numeric(Matrix::determinant(L, logarithm = TRUE,
                                             sqrt = TRUE)$modulus)
  ldRX2 <- 2 * sum(log(diag(RX)))
  if (glmm) {
    wrss <- sum(mm_glmm_residuals(fit, "pearson")^2)
    drsum <- sum(mm_glmm_residuals(fit, "deviance")^2)
  } else {
    wrss <- sum(mm_prior_weights(fit) * as.numeric(fit$residuals)^2)
    drsum <- NA_real_
  }
  pwrss <- wrss + ussq
  st <- mm_re_structure(fit)
  if (glmm) {
    use_sc <- identical(fit$family$family, "gamma")
    cmp <- c(ldL2 = ldL2, ldRX2 = ldRX2, wrss = wrss, ussq = ussq,
             pwrss = pwrss, drsum = drsum, REML = NA_real_,
             dev = -2 * as.numeric(fit$logLik),
             sigmaML = if (use_sc) as.numeric(sigma(fit)) else 1,
             sigmaREML = NA_real_)
    dims <- c(N = n, n = n, p = p, nmp = n - p, q = nrow(Zt),
              nth = length(fit$theta), nAGQ = as.integer(fit$nAGQ %||% 1L),
              compDev = 1L, useSc = as.integer(use_sc), reTrms = 1L,
              spFe = 0L, REML = 0L, GLMM = 1L, NLMM = 0L)
  } else {
    reml <- isTRUE(fit$REML)
    cmp <- c(ldL2 = ldL2, ldRX2 = ldRX2, wrss = wrss, ussq = ussq,
             pwrss = pwrss, drsum = NA_real_,
             REML = if (reml) -2 * as.numeric(fit$logLik) else NA_real_,
             dev = if (reml) NA_real_ else -2 * as.numeric(fit$logLik),
             sigmaML = sqrt(pwrss / n),
             sigmaREML = sqrt(pwrss / (n - p)))
    dims <- c(N = n, n = n, p = p, nmp = n - p, q = nrow(Zt),
              nth = length(fit$theta), useSc = 1L, reTrms = 1L, spFe = 0L,
              REML = if (reml) p else 0L, GLMM = 0L, NLMM = 0L)
  }
  storage.mode(dims) <- "integer"
  list(L = L, RZX = RZX, RX = RX, devcomp = list(cmp = cmp, dims = dims),
       LtZt = LtZt, weights = w, structure = st)
}

mm_devcomp <- function(fit) mm_pls(fit)$devcomp

# Conditional covariance of the random effects from the Laplace / PLS
# decomposition: sigma^2 * Lambda (Lambda'Z'WZ Lambda + I)^{-1} Lambda',
# returned per grouping factor as lme4's p x p x nlevels postVar arrays.
mm_r_cond_var_postvars <- function(fit) {
  pls <- mm_pls(fit)
  st <- pls$structure
  Lt <- mm_re_lambdat(fit)$Lambdat
  # M'M = Lambda A^{-1} Lambda' with M = L^{-1} P Lambda'
  M <- Matrix::solve(pls$L, Matrix::solve(pls$L, Lt, system = "P"),
                     system = "L")
  s2 <- as.numeric(sigma(fit))^2
  out <- list()
  offset <- 0L
  for (tm in st$terms) {
    nl <- length(tm$levels)
    arr <- array(0, dim = c(tm$p, tm$p, nl),
                 dimnames = list(tm$cnames, tm$cnames, tm$levels))
    for (l in seq_len(nl)) {
      cols <- offset + (l - 1L) * tm$p + seq_len(tm$p)
      arr[, , l] <- s2 * as.matrix(Matrix::crossprod(M[, cols, drop = FALSE]))
    }
    offset <- offset + nl * tm$p
    g <- tm$group
    if (is.null(out[[g]])) {
      out[[g]] <- arr
    } else {
      out[[g]] <- mm_merge_block_diag_postvar(out[[g]], arr, g)
    }
  }
  out
}

# ---- getME ----------------------------------------------------------------

mm_getme_names <- function() {
  c("X", "Z", "Zt", "Ztlist", "mmList", "y", "mu", "u", "b", "Gp", "Tp", "L",
    "Lambda", "Lambdat", "Lind", "Tlist", "A", "RX", "RZX", "sigma", "flist",
    "fixef", "beta", "theta", "ST", "REML", "is_REML", "n_rtrms", "n_rfacs",
    "N", "n", "p", "q", "p_i", "l_i", "q_i", "k", "m_i", "m", "cnms",
    "devcomp", "offset", "lower", "devfun", "glmer.nb.theta")
}

mm_getme_one <- function(object, one) {
  glmm <- inherits(object, "mm_glmm")
  st <- function() mm_re_structure(object)
  tnames <- function() {
    vapply(st()$terms, function(tm) tm$group, character(1))
  }
  theta_names <- function() {
    unlist(lapply(st()$terms, function(tm) {
      nm <- character(length(tm$theta_index))
      k <- 0L
      for (j in seq_len(tm$p)) for (i in seq_len(tm$p)) {
        if (tm$Tidx[i, j] > 0L) {
          k <- k + 1L
          nm[[k]] <- if (i == j) {
            paste(tm$group, tm$cnames[[i]], sep = ".")
          } else {
            paste(tm$group, tm$cnames[[i]], tm$cnames[[j]], sep = ".")
          }
        }
      }
      nm
    }))
  }
  refuse <- function(why) {
    mm_abort(
      message = sprintf("`getME()` component `%s` is not available: %s",
                        one, why),
      class = "mm_inference_unavailable",
      reason_code = "getme_component_unavailable",
      component = one,
      input = one
    )
  }
  switch(
    one,
    X = stats::model.matrix(object, type = "fixed"),
    Z = Matrix::t(mm_re_zt(object)),
    Zt = mm_re_zt(object),
    Ztlist = {
      out <- list()
      for (tm in st()$terms) {
        zt <- mm_re_term_zt(object, tm)
        for (j in seq_len(tm$p)) {
          rows <- seq.int(j, by = tm$p, length.out = length(tm$levels))
          piece <- zt[rows, , drop = FALSE]
          rownames(piece) <- tm$levels
          out[[paste(tm$group, tm$cnames[[j]], sep = ".")]] <- piece
        }
      }
      out
    },
    mmList = {
      out <- lapply(st()$terms, function(tm) {
        m <- do.call(cbind, mm_re_term_values(tm, object$model_frame))
        colnames(m) <- tm$cnames
        m
      })
      names(out) <- vapply(st()$terms, function(tm) {
        paste(paste(tm$cnames, collapse = " + "), "|", tm$group)
      }, character(1))
      out
    },
    y = {
      y <- as.numeric(mm_response_vector(object))
      names(y) <- rownames(object$model_frame)
      y
    },
    mu = as.numeric(object$fitted),
    u = mm_re_bu(object)$u,
    b = {
      b <- mm_re_bu(object)$b
      matrix(b, ncol = 1L)
    },
    Gp = as.integer(c(0L, cumsum(vapply(st()$terms, function(tm) {
      tm$p * length(tm$levels)
    }, numeric(1))))),
    Tp = {
      m_i <- vapply(st()$terms, function(tm) length(tm$theta_index), integer(1))
      stats::setNames(as.integer(c(0L, cumsum(m_i))), c("", tnames()))
    },
    L = mm_pls(object)$L,
    Lambda = Matrix::t(mm_re_lambdat(object)$Lambdat),
    Lambdat = mm_re_lambdat(object)$Lambdat,
    Lind = mm_re_lambdat(object)$Lind,
    Tlist = {
      out <- lapply(st()$terms, function(tm) {
        Tm <- matrix(0, tm$p, tm$p, dimnames = list(tm$cnames, tm$cnames))
        Tm[tm$Tidx > 0] <- st()$theta[tm$Tidx[tm$Tidx > 0]]
        Tm
      })
      names(out) <- tnames()
      out
    },
    ST = {
      out <- lapply(st()$terms, function(tm) {
        Tm <- matrix(0, tm$p, tm$p, dimnames = list(tm$cnames, tm$cnames))
        Tm[tm$Tidx > 0] <- st()$theta[tm$Tidx[tm$Tidx > 0]]
        S <- diag(Tm)
        ST <- Tm
        if (tm$p > 1L) {
          for (j in seq_len(tm$p)) {
            if (S[[j]] != 0) ST[, j] <- Tm[, j] / S[[j]]
          }
        }
        diag(ST) <- S
        ST
      })
      names(out) <- tnames()
      out
    },
    A = mm_re_lambdat(object)$Lambdat %*% mm_re_zt(object),
    RX = mm_pls(object)$RX,
    RZX = mm_pls(object)$RZX,
    sigma = as.numeric(sigma(object)),
    flist = {
      groups <- unique(tnames())
      out <- lapply(groups, function(g) {
        tm <- st()$terms[[match(g, tnames())]]
        tm$factor
      })
      names(out) <- groups
      attr(out, "assign") <- match(tnames(), groups)
      out
    },
    fixef = object$beta,
    beta = unname(object$beta),
    theta = stats::setNames(st()$theta, theta_names()),
    REML = if (glmm) 0L else if (isTRUE(object$REML)) length(object$beta) else 0L,
    is_REML = !glmm && isTRUE(object$REML),
    n_rtrms = length(st()$terms),
    n_rfacs = length(unique(tnames())),
    N = as.integer(nobs(object)),
    n = as.integer(nobs(object)),
    p = length(object$beta),
    q = nrow(mm_re_zt(object)),
    p_i = stats::setNames(vapply(st()$terms, function(tm) tm$p, integer(1)),
                          tnames()),
    l_i = stats::setNames(vapply(st()$terms, function(tm) length(tm$levels),
                                 integer(1)), tnames()),
    q_i = stats::setNames(vapply(st()$terms, function(tm) {
      tm$p * length(tm$levels)
    }, numeric(1)), tnames()),
    k = length(st()$terms),
    m_i = stats::setNames(vapply(st()$terms, function(tm) {
      length(tm$theta_index)
    }, integer(1)), tnames()),
    m = length(object$theta),
    cnms = {
      out <- lapply(st()$terms, function(tm) tm$cnames)
      names(out) <- tnames()
      out
    },
    devcomp = mm_devcomp(object),
    offset = if (glmm) mm_glmm_offset(object) else rep(0, nobs(object)),
    weights = mm_prior_weights(object),
    lower = stats::setNames(st()$lower, theta_names()),
    devfun = refuse(paste(
      "mixeff optimizes inside the Rust engine and does not expose an R",
      "deviance function; use profile() or verify_convergence() instead."
    )),
    glmer.nb.theta = {
      if (!glmm || !identical(object$family$family, "negative_binomial")) {
        refuse("the fit is not a negative-binomial GLMM.")
      }
      as.numeric(object$family$nb_theta)
    },
    mm_abort(
      message = sprintf(
        "`getME()` component `%s` is unknown; available components are %s.",
        one, paste(c(mm_getme_names(), "weights", "ALL"), collapse = ", ")
      ),
      class = "mm_arg_error",
      component = one,
      input = one
    )
  )
}

mm_getme <- function(object, name) {
  if (missing(name) || !is.character(name) || !length(name)) {
    mm_abort(
      message = "`name` must be a non-empty character vector.",
      class = "mm_arg_error",
      input = if (missing(name)) NULL else name
    )
  }
  if (identical(name, "ALL")) {
    nms <- setdiff(mm_getme_names(), "devfun")
    if (!(inherits(object, "mm_glmm") &&
          identical(object$family$family, "negative_binomial"))) {
      nms <- setdiff(nms, "glmer.nb.theta")
    }
    out <- lapply(nms, function(one) {
      tryCatch(mm_getme_one(object, one),
               error = function(cnd) NULL)
    })
    names(out) <- nms
    return(out[!vapply(out, is.null, logical(1))])
  }
  values <- lapply(name, function(one) mm_getme_one(object, one))
  names(values) <- name
  if (length(values) == 1L) values[[1L]] else values
}

# ---- singularity --------------------------------------------------------------

# lme4::isSingular(): any theta with a zero lower bound below `tol`.
mm_is_singular_theta <- function(x, tol) {
  st <- tryCatch(mm_re_structure(x), error = function(cnd) NULL)
  if (is.null(st)) return(NA)
  any(st$theta[st$lower == 0] < tol)
}

# ---- VarCorr --------------------------------------------------------------

# Build lme4's VarCorr.merMod shape from the stored (internal) variance
# component record: a named list with one covariance matrix per grouping
# factor carrying "stddev" and "correlation" attributes, and "sc"/"useSc"
# attributes on the list. mixeff's long table, residual SD and raw
# components stay reachable as `$table`, `$residual_sd`, `$components_raw`.
mm_varcorr_lme4 <- function(vc, fit) {
  comps <- vc$components_raw %||% list()
  frame <- fit$model_frame
  mats <- lapply(comps, function(comp) {
    nm <- mm_re_colnames_lme4(comp$names, frame)
    nm[nm %in% c("intercept", "1")] <- "(Intercept)"
    sd <- as.numeric(comp$std_dev)
    p <- length(sd)
    R <- diag(1, p)
    corr <- as.numeric(comp$correlations)
    if (p > 1L && length(corr)) {
      for (i in 2:p) {
        off <- (i - 1L) * (i - 2L) / 2L
        for (j in seq_len(i - 1L)) {
          R[i, j] <- R[j, i] <- corr[off + j]
        }
      }
    }
    dimnames(R) <- list(nm, nm)
    S <- R * outer(sd, sd)
    attr(S, "stddev") <- stats::setNames(sd, nm)
    attr(S, "correlation") <- R
    S
  })
  groups <- vapply(comps, function(comp) as.character(comp$group),
                   character(1))
  # lme4 layout: one entry per lme4 random-effect term, in mkReTrms order,
  # named make.unique(group) -- a `||` term gives `g`, `g.1`, ... -- cut
  # from the engine's per-block matrices (a diagonal block's off-diagonal
  # cells are structural zeros).
  lme4_mats <- mm_varcorr_lme4_terms(mats, groups, fit)
  if (!is.null(lme4_mats)) {
    mats <- lme4_mats
    groups <- attr(lme4_mats, "groups")
    attr(mats, "groups") <- NULL
  }
  names(mats) <- make.unique(groups, sep = ".")
  comps <- lapply(seq_along(mats), function(i) {
    S <- mats[[i]]
    R <- attr(S, "correlation")
    list(group = names(mats)[[i]], names = colnames(S),
         std_dev = unname(attr(S, "stddev")),
         correlations = if (ncol(S) > 1L) {
           unlist(lapply(2:ncol(S), function(r) R[r, seq_len(r - 1L)]),
                  use.names = FALSE)
         } else numeric())
  })
  glmm <- inherits(fit, "mm_glmm")
  use_sc <- if (glmm) identical(fit$family$family, "gamma") else TRUE
  sc <- if (use_sc) as.numeric(sigma(fit)) else 1
  attr(mats, "sc") <- sc
  attr(mats, "useSc") <- use_sc
  attr(mats, "mm_table") <- vc$table
  attr(mats, "mm_residual_sd") <- vc$residual_sd
  attr(mats, "mm_components_raw") <- comps
  attr(mats, "mm_boundary") <- mm_varcorr_boundary_flag(
    unlist(lapply(comps, `[[`, "std_dev"), use.names = FALSE),
    as.numeric(vc$residual_sd %||% NA_real_)
  )
  attr(mats, "mm_design_weak_identifiability_groups") <-
    attr(vc, "mm_design_weak_identifiability_groups")
  class(mats) <- c("mm_varcorr", "VarCorr.merMod")
  mats
}

# Re-cut per-block VarCorr matrices into lme4's terms (see mm_re_structure);
# NULL when the structure cannot be matched (the engine layout is kept).
mm_varcorr_lme4_terms <- function(mats, groups, fit) {
  st <- tryCatch(mm_re_structure(fit), error = function(cnd) NULL)
  if (is.null(st) || !length(mats)) return(NULL)
  used <- list()
  out <- list()
  out_groups <- character()
  for (tm in st$terms) {
    hit <- which(groups == tm$group & vapply(mats, function(S) {
      all(tm$cnames %in% colnames(S))
    }, logical(1)))
    if (!length(hit)) return(NULL)
    S <- mats[[hit[[1L]]]]
    keep <- match(tm$cnames, colnames(S))
    sub <- S[keep, keep, drop = FALSE]
    attributes(sub) <- list(dim = dim(sub), dimnames = dimnames(sub))
    attr(sub, "stddev") <- attr(S, "stddev")[keep]
    attr(sub, "correlation") <- attr(S, "correlation")[keep, keep, drop = FALSE]
    out[[length(out) + 1L]] <- sub
    out_groups <- c(out_groups, tm$group)
    used[[length(used) + 1L]] <- paste(hit[[1L]], keep)
  }
  # Every engine column must be accounted for exactly once.
  total <- sum(vapply(mats, ncol, integer(1)))
  if (length(unique(unlist(used))) != total ||
      length(unlist(used)) != total) {
    return(NULL)
  }
  attr(out, "groups") <- out_groups
  out
}

mm_varcorr_is_internal <- function(x) {
  "table" %in% names(unclass(x))
}

#' @export
`$.mm_varcorr` <- function(x, name) {
  ux <- unclass(x)
  if (name %in% names(ux)) return(ux[[name]])
  key <- c(table = "mm_table", residual_sd = "mm_residual_sd",
           components_raw = "mm_components_raw")
  if (name %in% names(key)) return(attr(x, key[[name]], exact = TRUE))
  NULL
}

mm_format_varcorr <- function(x, digits = max(3, getOption("digits") - 2)) {
  use_sc <- isTRUE(attr(x, "useSc"))
  sds <- c(lapply(unclass(x), attr, "stddev"),
           if (use_sc) list(Residual = unname(attr(x, "sc"))))
  lens <- lengths(sds)
  nr <- sum(lens)
  out <- array("", c(nr, 3L), list(rep.int("", nr),
                                   c("Groups", "Name", "Std.Dev.")))
  if (!nr) return(out)
  out[1 + cumsum(lens) - lens, "Groups"] <- names(lens)
  out[, "Name"] <- c(unlist(lapply(unclass(x), colnames)), if (use_sc) "")
  out[, "Std.Dev."] <- format(unlist(sds), digits = digits)
  if (any(lens > 1L)) {
    maxlen <- max(lens)
    co <- do.call(rbind, lapply(unclass(x), function(m) {
      cm <- attr(m, "correlation")
      dig <- max(2, digits - 2)
      cc <- format(round(cm, dig), nsmall = dig)
      cc[!lower.tri(cc)] <- ""
      if (nrow(cc) >= maxlen) return(cc)
      cbind(cc, matrix("", nrow(cc), maxlen - nrow(cc)))
    }))[, -maxlen, drop = FALSE]
    if (nrow(co) < nr) co <- rbind(co, matrix("", nr - nrow(co), ncol(co)))
    colnames(co) <- c("Corr", rep.int("", max(0L, ncol(co) - 1L)))
    out <- cbind(out, co, deparse.level = 0L)
  }
  out
}

# ---- multi-model information criteria ------------------------------------

mm_information_table <- function(objects, labels, k, criterion) {
  ll <- lapply(objects, stats::logLik)
  df <- vapply(ll, function(l) as.numeric(attr(l, "df")), numeric(1))
  val <- vapply(seq_along(ll), function(i) {
    l <- ll[[i]]
    if (identical(criterion, "AIC")) {
      -2 * as.numeric(l) + k * df[[i]]
    } else {
      -2 * as.numeric(l) + log(as.numeric(attr(l, "nobs"))) * df[[i]]
    }
  }, numeric(1))
  nobs_all <- vapply(ll, function(l) as.numeric(attr(l, "nobs") %||% NA),
                     numeric(1))
  if (length(unique(nobs_all)) > 1L) {
    warning("models are not all fitted to the same number of observations")
  }
  out <- data.frame(df = df, val)
  names(out)[[2L]] <- criterion
  row.names(out) <- make.unique(labels)
  out
}
