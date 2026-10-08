#' Per-group linear (or generalized linear) model fits
#'
#' `mm_lmlist()` is mixeff's counterpart of `lme4::lmList()`: it fits
#' `lm()` (or `glm()` when `family` is given) separately within each level
#' of a grouping factor, using lme4's `response ~ predictors | group`
#' formula. It is named `mm_lmlist()` so attaching mixeff does not mask
#' `lme4::lmList()` / `nlme::lmList()`.
#'
#' `coef()` returns one row per group (lme4's shape); `confint()` returns a
#' `group x (lower, upper) x coefficient` array. With `pool = TRUE` (the
#' default for `lm` fits, as in lme4) `confint()` uses the residual standard
#' deviation pooled across groups, with the pooled residual degrees of
#' freedom, and `sigma()` returns that pooled value; `confint(x, pool =
#' FALSE)` gives each group's own `lm` t interval (Wald z for GLMs), which
#' is what `confint()` on an `lme4::lmList()` fit returns. A group
#' whose fit fails is kept as `NULL` and reported, never silently dropped.
#'
#' @param formula A formula `y ~ x | g`.
#' @param data A data frame.
#' @param family Optional GLM family; `NULL` fits `lm()`.
#' @param subset,weights,na.action Passed to `lm()`/`glm()` within each
#'   group (`subset` is applied before splitting). `weights` may name a
#'   column of `data`.
#' @param pool Pool the residual standard deviation across groups for
#'   inference. Defaults to `TRUE` for `lm` fits and `FALSE` for GLMs.
#'
#' @return An object of class `mm_lmlist`: a named list of `lm`/`glm` fits
#'   (one per group) with attributes `call`, `pool`, `groups`, and
#'   `failed`.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:6, each = 8)), x = rep(0:7, 6))
#' d$y <- 2 + rnorm(6)[d$g] + (0.5 + rnorm(6, sd = 0.2)[d$g]) * d$x +
#'   rnorm(48)
#' ml <- mm_lmlist(y ~ x | g, d)
#' coef(ml)
#' confint(ml)[, , "x"]
#' @export
mm_lmlist <- function(formula, data, family = NULL, subset = NULL,
                      weights = NULL, na.action = stats::na.omit,
                      pool = is.null(family)) {
  call <- match.call()
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.call(formula[[3L]]) ||
      !identical(formula[[3L]][[1L]], as.name("|"))) {
    mm_abort(message = "`formula` must have the form `response ~ predictors | group`.",
             class = "mm_formula_error", input = formula)
  }
  if (!is.data.frame(data)) {
    mm_abort(message = "`data` must be a data.frame.", class = "mm_data_error",
             input = data)
  }
  rhs <- formula[[3L]]
  group_expr <- rhs[[3L]]
  model_formula <- formula
  model_formula[[3L]] <- rhs[[2L]]
  environment(model_formula) <- environment(formula)

  w <- eval(substitute(weights), data, parent.frame())
  if (!is.null(w)) data$.mm_lmlist_weights <- w
  sub <- substitute(subset)
  if (!is.null(sub)) {
    keep <- eval(sub, data, parent.frame())
    data <- data[keep & !is.na(keep), , drop = FALSE]
  }
  group <- eval(group_expr, data, environment(formula))
  group <- if (is.factor(group)) droplevels(group) else factor(group)
  pieces <- split(data, group)

  fam <- family
  if (is.character(fam)) fam <- get(fam, mode = "function")
  if (is.function(fam)) fam <- fam()
  fits <- lapply(pieces, function(piece) {
    args <- list(formula = model_formula, data = piece, na.action = na.action)
    if (!is.null(w)) args$weights <- piece$.mm_lmlist_weights
    tryCatch(
      if (is.null(fam)) do.call(stats::lm, args) else
        do.call(stats::glm, c(args, list(family = fam))),
      error = function(cnd) cnd
    )
  })
  failed <- vapply(fits, inherits, logical(1), what = "condition")
  if (any(failed)) {
    mm_inform(
      sprintf("mm_lmlist(): the fit failed for %d group(s): %s",
              sum(failed), paste(names(fits)[failed], collapse = ", ")),
      class = "mm_lmlist_failure_notice"
    )
    msgs <- vapply(fits[failed], conditionMessage, character(1))
    fits[failed] <- list(NULL)
  } else {
    msgs <- character()
  }
  structure(fits, class = "mm_lmlist", call = call, pool = isTRUE(pool),
            groups = group, failed = msgs, family = fam)
}

#' @method print mm_lmlist
#' @export
print.mm_lmlist <- function(x, ...) {
  cat("Per-group fits (mm_lmlist):", deparse1(attr(x, "call")), "\n")
  cat("Coefficients:\n")
  print(coef(x))
  if (isTRUE(attr(x, "pool"))) {
    cat(sprintf("\nDegrees of freedom: %d total; residual standard error (pooled): %.5g\n",
                sum(vapply(Filter(Negate(is.null), x), function(f) length(f$residuals), integer(1))),
                sigma(x)))
  }
  invisible(x)
}

#' @method coef mm_lmlist
#' @export
coef.mm_lmlist <- function(object, ...) {
  ok <- Filter(Negate(is.null), unclass(object))
  nms <- unique(unlist(lapply(ok, function(f) names(stats::coef(f)))))
  out <- matrix(NA_real_, nrow = length(object), ncol = length(nms),
                dimnames = list(names(object), nms))
  for (g in names(ok)) {
    cf <- stats::coef(ok[[g]])
    out[g, names(cf)] <- cf
  }
  as.data.frame(out, check.names = FALSE)
}

#' @method sigma mm_lmlist
#' @export
sigma.mm_lmlist <- function(object, ...) {
  ok <- Filter(Negate(is.null), unclass(object))
  rss <- sum(vapply(ok, function(f) sum(stats::weighted.residuals(f)^2), numeric(1)))
  dfr <- sum(vapply(ok, stats::df.residual, numeric(1)))
  sqrt(rss / dfr)
}

#' @method confint mm_lmlist
#' @export
confint.mm_lmlist <- function(object, parm, level = 0.95,
                              pool = attr(object, "pool"), ...) {
  cf <- coef(object)
  nms <- colnames(cf)
  if (!missing(parm)) nms <- if (is.numeric(parm)) nms[parm] else intersect(nms, parm)
  a <- (1 - level) / 2
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE,
                      digits = 3), "%")
  out <- array(NA_real_, dim = c(nrow(cf), 2L, length(nms)),
               dimnames = list(rownames(cf), pct, nms))
  pool <- isTRUE(pool) && is.null(attr(object, "family"))
  if (pool) {
    s <- sigma(object)
    ok <- Filter(Negate(is.null), unclass(object))
    dfp <- sum(vapply(ok, stats::df.residual, numeric(1)))
    q <- stats::qt(1 - a, dfp)
  }
  for (g in rownames(cf)) {
    f <- object[[g]]
    if (is.null(f)) next
    est <- stats::coef(f)
    if (pool) {
      V <- stats::vcov(f) / stats::sigma(f)^2 * s^2
      se <- sqrt(diag(V))
      lo <- est - q * se
      hi <- est + q * se
    } else {
      ci <- if (inherits(f, "glm")) {
        stats::confint.default(f, level = level)
      } else {
        stats::confint(f, level = level)
      }
      lo <- ci[, 1L]
      hi <- ci[, 2L]
    }
    keep <- intersect(nms, names(est))
    out[g, 1L, keep] <- lo[keep]
    out[g, 2L, keep] <- hi[keep]
  }
  class(out) <- c("mm_lmlist_confint", class(out))
  out
}

#' @method print mm_lmlist_confint
#' @export
print.mm_lmlist_confint <- function(x, ...) {
  print(unclass(x), ...)
  invisible(x)
}

#' @method fixef mm_lmlist
#' @export
fixef.mm_lmlist <- function(object, ...) {
  colMeans(coef(object), na.rm = TRUE)
}

#' Refit a model with every available optimizer
#'
#' `mm_allfit()` is mixeff's counterpart of `lme4::allFit()`: it refits
#' `fit` once per optimizer (through [update()] with
#' `mm_control(optimizer = )`, keeping the rest of the fit's control) and
#' summarises the fixed effects, covariance parameters, log-likelihood,
#' convergence status and run time side by side, so a fit's sensitivity to
#' the optimizer is visible. An optimizer that is not compiled into this
#' build, or whose refit fails, is reported as a failed row with its error
#' message -- never dropped.
#'
#' @param fit A fitted `mm_lmm` or `mm_glmm`.
#' @param optimizers Character vector of [mm_control()] optimizer names.
#'   The default tries every name `mm_control()` accepts.
#' @param verbose Print each optimizer as it runs.
#'
#' @return An object of class `mm_allfit`: a named list of fits (or
#'   conditions for failed refits). `summary()` returns a list with
#'   `which.OK`, `msgs`, `fixef`, `theta`, `llik`, `fit_status`, and
#'   `times`, mirroring `summary(lme4::allFit(...))`.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:8, each = 5)), x = rep(0:4, 8))
#' d$y <- 1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(40)
#' fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
#' af <- mm_allfit(fit, optimizers = c("auto", "pattern_search", "cobyla"))
#' summary(af)$llik
#' @export
mm_allfit <- function(fit, optimizers = mm_allfit_optimizers(),
                      verbose = FALSE) {
  if (!inherits(fit, c("mm_lmm", "mm_glmm"))) {
    mm_abort(message = "`fit` must be a fitted `mm_lmm` or `mm_glmm`.",
             class = "mm_arg_error", input = fit)
  }
  if (!is.character(optimizers) || !length(optimizers) || anyNA(optimizers)) {
    mm_abort(message = "`optimizers` must be a non-empty character vector.",
             class = "mm_arg_error", input = optimizers)
  }
  base <- unclass(mm_internal_control(fit))
  base$start <- NULL
  out <- list()
  times <- numeric()
  for (opt in optimizers) {
    if (isTRUE(verbose)) cat(opt, "\n")
    ctl <- base
    ctl$optimizer <- if (identical(opt, "auto")) NULL else opt
    ctl <- mm_validate_control(ctl)
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(stats::update(fit, control = ctl), error = function(cnd) cnd)
    times[[opt]] <- proc.time()[["elapsed"]] - t0
    out[[opt]] <- res
  }
  structure(out, class = "mm_allfit", times = times,
            original = fit)
}

mm_allfit_optimizers <- function() {
  c("auto", "bobyqa", "newuoa", "cobyla", "pattern_search", "trust_bq",
    "prima_bobyqa", "prima_cobyla", "prima_lincoa", "prima_newuoa")
}

#' @method summary mm_allfit
#' @export
summary.mm_allfit <- function(object, ...) {
  ok <- !vapply(object, inherits, logical(1), what = "condition")
  fits <- unclass(object)[ok]
  msgs <- lapply(unclass(object), function(x) {
    if (inherits(x, "condition")) conditionMessage(x) else NULL
  })
  stack <- function(f) {
    if (!length(fits)) return(NULL)
    do.call(rbind, lapply(fits, f))
  }
  list(
    which.OK = ok,
    msgs = msgs,
    fixef = stack(function(m) fixef(m)),
    theta = stack(function(m) m$theta),
    llik = vapply(fits, function(m) as.numeric(stats::logLik(m)), numeric(1)),
    fit_status = vapply(fits, fit_status, character(1)),
    times = attr(object, "times")[ok]
  )
}

#' @method print mm_allfit
#' @export
print.mm_allfit <- function(x, ...) {
  s <- summary(x)
  cat(sprintf("mm_allfit: %d of %d optimizer(s) succeeded\n",
              sum(s$which.OK), length(s$which.OK)))
  if (any(s$which.OK)) {
    tab <- data.frame(logLik = s$llik, fit_status = s$fit_status,
                      seconds = round(s$times, 3), check.names = FALSE)
    tab$logLik_diff <- tab$logLik - max(tab$logLik)
    print(tab)
  }
  failed <- names(s$which.OK)[!s$which.OK]
  if (length(failed)) {
    cat("Failed:\n")
    for (f in failed) cat(sprintf("  %s: %s\n", f, s$msgs[[f]]))
  }
  invisible(x)
}
