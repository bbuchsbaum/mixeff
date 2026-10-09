#' Parametric bootstrap of any statistic of a fitted LMM
#'
#' `mm_bootmer()` is mixeff's counterpart of `lme4::bootMer()`: it simulates
#' `nsim` responses from the fitted model, refits the model to each, and
#' evaluates `FUN` on every refit. The result has class `c("mm_bootmer",
#' "boot")` with the fields `boot::boot.ci()` reads (`t0`, `t`, `R`, `sim`,
#' `statistic`, ...), so `boot::boot.ci(b, type = "perc", index = 1)` works
#' as it does for a `bootMer` result.
#'
#' The function is named `mm_bootmer()` rather than `bootMer()` because
#' `lme4::bootMer()` is not an S3 generic: exporting a `bootMer()` would mask
#' lme4's for `merMod` fits.
#'
#' Replicates use the same model, `REML` setting, weights, and
#' [mm_control()] as `fit`; each refit warm-starts from the fitted
#' covariance parameters (as lme4's `refit()` does). A replicate whose refit
#' fails, or whose `FUN` errors, is recorded as an `NA` row and counted in
#' `attr(, "bootFail")` (lme4's convention) -- never silently dropped.
#'
#' @param fit A fitted `mm_lmm`. GLMM fits are refused (the GLMM parametric
#'   bootstrap is `confint(fit, method = "bootstrap")`).
#' @param FUN A function of a fitted `mm_lmm` returning a numeric vector
#'   (for example `fixef`, or `function(m) c(fixef(m), sigma = sigma(m))`).
#' @param nsim Number of bootstrap replicates.
#' @param seed Optional seed; `NULL` uses the current R RNG state, so
#'   `set.seed()` governs reproducibility.
#' @param type Only `"parametric"` is available. lme4's `"semiparametric"`
#'   residual bootstrap is refused.
#' @param use.u `FALSE` (default) draws new random effects for each
#'   replicate (lme4's default, the unconditional bootstrap); `TRUE`
#'   conditions on the fitted random effects and resamples only the
#'   residual noise.
#' @param verbose Print a progress dot per replicate.
#'
#' @return An object of class `c("mm_bootmer", "boot")`: a list with `t0`
#'   (`FUN(fit)`), `t` (an `nsim` x `length(t0)` matrix), `R`, `data`,
#'   `seed`, `statistic`, `sim = "parametric"`, `call`, and `mle`, plus
#'   attributes `bootFail` (failed replicates) and `boot.fail.msgs`.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:8, each = 5)), x = rep(0:4, 8))
#' d$y <- 1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(40)
#' fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
#' b <- mm_bootmer(fit, fixef, nsim = 20, seed = 2)
#' apply(b$t, 2, sd)
#' if (requireNamespace("boot", quietly = TRUE)) {
#'   boot::boot.ci(b, type = "perc", index = 2)
#' }
#'
#' @seealso [confint()] with `method = "bootstrap"` for engine-run
#'   fixed-effect bootstrap intervals; [simulate.mm_lmm()]; [refit()].
#' @export
mm_bootmer <- function(fit, FUN, nsim = 1L, seed = NULL,
                       type = c("parametric", "semiparametric"),
                       use.u = FALSE, verbose = FALSE) {
  call <- match.call()
  if (inherits(fit, "mm_glmm")) {
    mm_abort(
      message = paste0(
        "`mm_bootmer()` refits through `refit()`, which is not implemented ",
        "for GLMM fits. Use confint(fit, method = \"bootstrap\") for the ",
        "engine-run GLMM parametric bootstrap."
      ),
      class = "mm_inference_unavailable",
      reason_code = "glmm_bootmer_unimplemented",
      input = class(fit)[[1L]]
    )
  }
  if (!inherits(fit, "mm_lmm")) {
    mm_abort(message = "`fit` must be a fitted `mm_lmm`.",
             class = "mm_arg_error", input = fit)
  }
  type <- match.arg(type)
  if (!identical(type, "parametric")) {
    mm_abort(
      message = paste0(
        "`type = \"semiparametric\"` (lme4's residual bootstrap) is not ",
        "available; use type = \"parametric\"."
      ),
      class = "mm_inference_unavailable",
      reason_code = "semiparametric_bootstrap_unimplemented",
      input = type
    )
  }
  FUN <- match.fun(FUN)
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1 ||
      nsim != round(nsim)) {
    mm_abort(message = "`nsim` must be a single positive integer.",
             class = "mm_arg_error", input = nsim)
  }
  nsim <- as.integer(nsim)
  if (!is.logical(use.u) || length(use.u) != 1L || is.na(use.u)) {
    mm_abort(message = "`use.u` must be TRUE or FALSE.",
             class = "mm_arg_error", input = use.u)
  }

  t0 <- FUN(fit)
  if (!is.numeric(t0) || !length(t0)) {
    mm_abort(message = "`FUN(fit)` must return a non-empty numeric vector.",
             class = "mm_arg_error", input = t0)
  }

  rng_state <- mm_rng_state(seed)
  ysim <- mm_with_seed(seed, {
    if (isTRUE(use.u)) {
      w <- fit$weights %||% rep(1, nobs(fit))
      replicate(nsim, as.numeric(fit$fitted) +
                  stats::rnorm(nobs(fit), sd = fit$sigma / sqrt(w)),
                simplify = FALSE)
    } else {
      as.list(simulate(fit, nsim = nsim))
    }
  })

  control <- mm_internal_control(fit)
  if (is.null(control$start) && length(fit$theta) &&
      all(is.finite(fit$theta))) {
    control$start <- as.numeric(fit$theta)
  }
  t <- matrix(NA_real_, nrow = nsim, ncol = length(t0),
              dimnames = list(NULL, names(t0)))
  fail_msgs <- character()
  for (i in seq_len(nsim)) {
    res <- tryCatch({
      refit_i <- refit(fit, ysim[[i]], control = control)
      val <- FUN(refit_i)
      if (!is.numeric(val) || length(val) != length(t0)) {
        stop("`FUN` returned a value of a different length than `FUN(fit)`.")
      }
      val
    }, error = function(cnd) cnd)
    if (inherits(res, "condition")) {
      fail_msgs <- c(fail_msgs, conditionMessage(res))
    } else {
      t[i, ] <- res
    }
    if (isTRUE(verbose)) cat(".")
  }
  if (isTRUE(verbose)) cat("\n")

  out <- list(
    t0 = t0,
    t = t,
    R = nsim,
    data = fit$model_frame,
    seed = rng_state,
    statistic = FUN,
    sim = "parametric",
    call = call,
    ran.gen = "simulate(<mm_lmm>, 1, *)",
    mle = list(beta = fit$beta, theta = fit$theta, sigma = fit$sigma),
    use.u = use.u
  )
  n_fail <- length(fail_msgs)
  attr(out, "bootFail") <- n_fail
  attr(out, "boot.fail.msgs") <- table(fail_msgs)
  class(out) <- c("mm_bootmer", "boot")
  if (n_fail > 0L) {
    mm_inform(
      sprintf(paste0("mm_bootmer(): %d of %d replicate(s) failed and are ",
                     "recorded as NA rows; see attr(, \"boot.fail.msgs\")."),
              n_fail, nsim),
      class = "mm_bootstrap_failure_notice"
    )
  }
  out
}

#' @method print mm_bootmer
#' @export
print.mm_bootmer <- function(x, ...) {
  cat(sprintf("mixeff parametric bootstrap (mm_bootmer): %d replicates%s\n",
              x$R, if (isTRUE(x$use.u)) ", conditional on u" else ""))
  fail <- attr(x, "bootFail") %||% 0L
  if (fail > 0L) cat(sprintf("failed replicates: %d (NA rows)\n", fail))
  bias <- colMeans(x$t, na.rm = TRUE) - x$t0
  se <- apply(x$t, 2L, stats::sd, na.rm = TRUE)
  tab <- data.frame(original = x$t0, bias = bias, std.error = se,
                    check.names = FALSE)
  print(tab)
  invisible(x)
}

#' PCA of the random-effects covariance (lme4's `rePCA()`)
#'
#' For each grouping factor, the singular values and left singular vectors
#' of the relative covariance factor (the block-diagonal Cholesky factor of
#' the random-effect covariance divided by the residual variance, lme4's
#' `Lambda`), returned as `prcomp` objects. Near-zero standard deviations
#' flag a random-effect structure that is over-parameterised (singular).
#' `summary()` gives each component's proportion of variance.
#'
#' Singular values match `lme4::rePCA()`; the rotation columns are
#' eigenvectors of the same matrix and may differ from lme4's in sign.
#'
#' @param x A fitted `mm_lmm` or `mm_glmm` (for GLMM families without a
#'   dispersion parameter the covariance itself is used, as in lme4).
#'
#' @return An object of class `c("mm_prcomplist", "prcomplist")`: a named list (one entry per
#'   grouping factor) of `prcomp` objects.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:10, each = 6)), x = rep(0:5, 10))
#' d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(10, sd = 0.3)[d$g] * d$x +
#'   rnorm(60)
#' fit <- lmm(y ~ x + (x | g), d, control = mm_control(verbose = -1))
#' summary(rePCA(fit))
#' @export
rePCA <- function(x) {
  UseMethod("rePCA")
}

#' @rdname rePCA
#' @export
rePCA.default <- function(x) {
  fwd <- mm_forward_foreign_generic("rePCA", x)
  if (!is.null(fwd)) return(fwd$value)
  mm_abort(message = "`rePCA()` has no method for this object.",
           class = "mm_arg_error", input = x)
}

#' @rdname rePCA
#' @export
rePCA.mm_lmm <- function(x) {
  blocks <- mm_relative_covariance_blocks(x)
  groups <- unique(names(blocks))
  names(groups) <- groups
  out <- lapply(groups, function(g) {
    parts <- blocks[names(blocks) == g]
    S <- as.matrix(Matrix::bdiag(parts))
    labels <- unlist(lapply(parts, colnames), use.names = FALSE)
    S <- (S + t(S)) / 2
    eig <- eigen(S, symmetric = TRUE)
    pc <- list(sdev = sqrt(pmax(eig$values, 0)),
               rotation = eig$vectors,
               center = FALSE, scale = FALSE)
    dimnames(pc$rotation) <- list(labels, NULL)
    class(pc) <- "prcomp"
    pc
  })
  class(out) <- c("mm_prcomplist", "prcomplist")
  out
}

#' @rdname rePCA
#' @export
rePCA.mm_glmm <- rePCA.mm_lmm

#' @method summary mm_prcomplist
#' @export
summary.mm_prcomplist <- function(object, ...) {
  lapply(object, summary)
}

# Per random term, the covariance matrix divided by the residual variance
# (lme4's Lambda %*% t(Lambda)), named by grouping factor. GLMM covariance
# parameters are absolute (lme4 >= 2.1-0 for the free-dispersion families
# too), so GLMM blocks are not rescaled.
mm_relative_covariance_blocks <- function(fit) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  scale2 <- if (inherits(fit, "mm_glmm")) 1 else as.numeric(fit$sigma)^2
  if (!length(scale2) || !is.finite(scale2) || scale2 <= 0) scale2 <- 1
  blocks <- lapply(seq_along(terms), function(i) {
    term <- terms[[i]]
    term_id <- term$id %||% sprintf("r%d", i - 1L)
    basis <- term$basis %||% list()
    labels <- vapply(basis, mm_basis_label, character(1))
    if (!length(labels)) labels <- "(Intercept)"
    group <- mm_random_term_group_label(fit, term, i)
    S <- mm_random_term_covariance(fit, term_id, labels, group) / scale2
    list(group = group, S = S)
  })
  out <- lapply(blocks, `[[`, "S")
  names(out) <- vapply(blocks, `[[`, character(1), "group")
  out
}
