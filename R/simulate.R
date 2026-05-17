#' Refit a mixeff LMM with a new response
#'
#' `refit()` fits the same model formula to a new response by calling [lmm()]
#' with the stored model frame and `REML` setting.
#'
#' @param object A fitted `mm_lmm`.
#' @param newresp Numeric response for `refit()`.
#' @param ... Reserved for future methods.
#'
#' @return A new `mm_lmm`.
#'
#' @export
refit <- function(object, newresp, ...) {
  UseMethod("refit")
}

#' @rdname refit
#' @export
refit.mm_lmm <- function(object, newresp, ...) {
  if (!is.numeric(newresp) || length(newresp) != nobs(object) ||
      anyNA(newresp)) {
    mm_abort(
      message = "`newresp` must be a numeric vector with one value per observation and no missing values.",
      class = "mm_arg_error",
      input = newresp
    )
  }
  data <- object$model_frame
  data[[mm_response_name(object)]] <- as.numeric(newresp)
  control <- list(...)$control %||% mm_control(verbose = -1)
  fit <- lmm(object$formula, data, REML = isTRUE(object$REML),
             weights = object$weights, control = control)
  fit$refit <- list(
    source = "refit",
    original_fit_status = fit_status(object)
  )
  fit
}

#' Simulate from a mixeff LMM
#'
#' Draws Gaussian responses from the stored fixed effects, random-effect
#' covariance summaries, and residual scale.
#'
#' @param object A fitted `mm_lmm`.
#' @param nsim Number of simulated responses.
#' @param seed Optional random seed.
#' @param re.form Random-effects conditioning. `NULL` simulates new random
#'   effects; `NA` simulates from the population-level mean only.
#' @param ... Reserved for future methods.
#'
#' @return A data frame of simulated responses.
#'
#' @importFrom stats simulate
#' @method simulate mm_lmm
#' @export
simulate.mm_lmm <- function(object, nsim = 1, seed = NULL, re.form = NULL, ...) {
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    mm_abort(
      message = "`nsim` must be a positive integer.",
      class = "mm_arg_error",
      input = nsim
    )
  }
  nsim <- as.integer(nsim)
  target <- mm_prediction_target(re.form)
  if (!target %in% c("conditional", "population")) {
    mm_abort(
      message = "`re.form` requests beyond NULL and NA are not available for simulation.",
      class = "mm_inference_unavailable",
      input = re.form
    )
  }

  out <- mm_with_seed(seed, {
    sims <- replicate(nsim, mm_simulate_once(object, target), simplify = FALSE)
    as.data.frame(stats::setNames(sims, paste0("sim_", seq_len(nsim))),
                  check.names = FALSE)
  })
  rownames(out) <- rownames(object$model_frame)
  attr(out, "seed") <- seed
  attr(out, "mm_method") <- "r_side_gaussian_parametric"
  out
}

mm_simulate_once <- function(fit, target) {
  eta <- if (identical(target, "population")) {
    fit$fixed_fitted
  } else {
    mm_simulate_random_mean(fit)
  }
  as.numeric(eta + stats::rnorm(nobs(fit), sd = fit$sigma))
}

mm_simulate_random_mean <- function(fit) {
  eta <- as.numeric(fit$fixed_fitted)
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    term_id <- term$id %||% sprintf("r%d", i - 1L)
    group_label <- mm_random_term_group_label(fit, term, i)
    group <- mm_group_factor(fit$model_frame, group_label)
    levels <- levels(group)
    basis <- term$basis %||% list()
    basis_labels <- vapply(basis, mm_basis_label, character(1))
    basis_values <- lapply(basis, mm_basis_values, frame = fit$model_frame)
    if (!length(basis_values)) {
      basis_labels <- "(Intercept)"
      basis_values <- list(rep(1, nobs(fit)))
    }
    Sigma <- mm_random_term_covariance(fit, term_id, basis_labels)
    draws <- mm_rmvnorm(length(levels), Sigma)
    idx <- as.integer(group)
    for (j in seq_along(basis_values)) {
      eta <- eta + basis_values[[j]] * draws[idx, j]
    }
  }
  eta
}

mm_random_term_covariance <- function(fit, term_id, basis_labels) {
  p <- length(basis_labels)
  if (!p) return(matrix(0, 0, 0))
  labels <- ifelse(basis_labels == "(Intercept)", "intercept", basis_labels)
  Sigma <- diag(0, p)
  dimnames(Sigma) <- list(basis_labels, basis_labels)

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

  missing_diag <- diag(Sigma) <= 0
  if (any(missing_diag)) {
    vc <- fit$varcorr$table
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
