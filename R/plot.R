# Diagnostic plots for mixeff fits and random effects (pre-CRAN audit 4.5).
# Base graphics; lattice dotplot()/qqmath() methods for ranef() output are
# registered lazily when lattice is loaded (see R/zzz.R).

#' Diagnostic plots for mixeff fits
#'
#' `plot()` on a fitted model draws Pearson residuals against fitted values,
#' the default diagnostic of `lme4`'s `plot.merMod()`, using base graphics.
#' `which = 2` draws a normal quantile-quantile plot of the Pearson residuals
#' instead (also available as `qqnorm(fit)`), and `which = 3` the
#' scale-location plot (`sqrt(|r|)` against fitted values).
#'
#' Pearson residuals follow lme4's definition, `(y - mu) sqrt(w) / sqrt(V(mu))`
#' with prior weights `w` and variance function `V` (`V = 1` for an LMM), on
#' the response scale; they are not divided by the residual standard
#' deviation.
#'
#' @param x,y A fitted `mm_lmm` or `mm_glmm`.
#' @param which Integer vector selecting the plots: `1` residuals vs fitted
#'   (default), `2` normal Q-Q of the residuals, `3` scale-location.
#' @param smooth If `TRUE` (default), add a `lowess()` smooth to plots 1 and 3.
#' @param ask Ask before each new plot when several are drawn on an
#'   interactive device.
#' @param main,xlab,ylab Optional titles overriding the defaults.
#' @param ... Further graphical parameters passed to [graphics::plot()] or
#'   [stats::qqnorm()].
#'
#' @return `x`, invisibly. `qqnorm()` returns the list from
#'   [stats::qqnorm()] invisibly.
#'
#' @seealso [plot.mm_ranef()] for caterpillar and Q-Q plots of the random
#'   effects; [hatvalues.mm_lmm()] for leverage and Cook's distance.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rnorm(60), x = rnorm(60),
#'   g = factor(rep(seq_len(10), each = 6))
#' )
#' fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
#' plot(fit)
#' plot(fit, which = 2)
#' qqnorm(fit)
#'
#' @name plot.mm_fit
#' @importFrom stats qqnorm
#' @method plot mm_lmm
#' @export
plot.mm_lmm <- function(x, which = 1L, smooth = TRUE,
                        ask = length(which) > 1L && grDevices::dev.interactive(),
                        main = NULL, xlab = NULL, ylab = NULL, ...) {
  if (!is.numeric(which) || !length(which) || anyNA(which) ||
      !all(which %in% 1:3)) {
    mm_abort(
      message = "`which` must contain plot numbers 1 (residuals vs fitted), 2 (normal Q-Q) and/or 3 (scale-location).",
      class = "mm_arg_error",
      input = which
    )
  }
  if (isTRUE(ask)) {
    old <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(old), add = TRUE)
  }
  parts <- mm_influence_parts(x)
  res <- parts$pearson
  ftd <- parts$mu
  for (w in which) {
    if (w == 1L) {
      graphics::plot(ftd, res,
                     xlab = xlab %||% "Fitted values",
                     ylab = ylab %||% "Pearson residuals",
                     main = main %||% "Residuals vs fitted", ...)
      graphics::abline(h = 0, lty = 2, col = "grey50")
      if (isTRUE(smooth)) mm_plot_smooth(ftd, res)
    } else if (w == 2L) {
      mm_qqnorm_residuals(res, main = main, xlab = xlab, ylab = ylab, ...)
    } else {
      sres <- sqrt(abs(res))
      graphics::plot(ftd, sres,
                     xlab = xlab %||% "Fitted values",
                     ylab = ylab %||% expression(sqrt(abs("Pearson residuals"))),
                     main = main %||% "Scale-location", ...)
      if (isTRUE(smooth)) mm_plot_smooth(ftd, sres)
    }
  }
  invisible(x)
}

#' @rdname plot.mm_fit
#' @method plot mm_glmm
#' @export
plot.mm_glmm <- plot.mm_lmm

#' @rdname plot.mm_fit
#' @method qqnorm mm_lmm
#' @export
qqnorm.mm_lmm <- function(y, main = NULL, xlab = NULL, ylab = NULL, ...) {
  invisible(mm_qqnorm_residuals(mm_influence_parts(y)$pearson,
                                main = main, xlab = xlab, ylab = ylab, ...))
}

#' @rdname plot.mm_fit
#' @method qqnorm mm_glmm
#' @export
qqnorm.mm_glmm <- qqnorm.mm_lmm

mm_qqnorm_residuals <- function(res, main = NULL, xlab = NULL, ylab = NULL,
                                ...) {
  out <- stats::qqnorm(res,
                       main = main %||% "Normal Q-Q plot of Pearson residuals",
                       xlab = xlab %||% "Standard normal quantiles",
                       ylab = ylab %||% "Pearson residuals", ...)
  stats::qqline(res, lty = 2, col = "grey50")
  invisible(out)
}

mm_plot_smooth <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) > 3L) {
    graphics::lines(stats::lowess(x[ok], y[ok]), col = "red")
  }
}

#' Caterpillar and Q-Q plots of random effects
#'
#' `plot()` on the result of [ranef()] draws a caterpillar plot (dot chart)
#' of the conditional modes for each grouping factor and coefficient, with
#' levels sorted by their value and, when the object carries conditional
#' variances (`ranef(fit, condVar = TRUE)`), `level` prediction intervals,
#' mirroring `lattice::dotplot()` on `lme4`'s `ranef()` output. `qqnorm()`
#' draws normal Q-Q plots of the conditional modes, one panel per grouping
#' factor and coefficient.
#'
#' When the lattice package is loaded, `lattice::dotplot()` and
#' `lattice::qqmath()` also work on `ranef()` output and return one trellis
#' object per grouping factor, as they do for lme4.
#'
#' Intervals are omitted for coefficients whose conditional variances are not
#' available (for example GLMM fits where `ranef(condVar = TRUE)` returns
#' `NA` variances).
#'
#' @param x,y An `mm_ranef` object from [ranef()].
#' @param level Coverage of the intervals (default 0.95).
#' @param type `"caterpillar"` (default) or `"qq"`.
#' @param ... Further graphical parameters passed to [graphics::plot()].
#'
#' @return `x` (or `y`), invisibly.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rnorm(60), x = rnorm(60),
#'   g = factor(rep(seq_len(10), each = 6))
#' )
#' fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
#' re <- ranef(fit, condVar = TRUE)
#' plot(re)
#' qqnorm(re)
#'
#' @method plot mm_ranef
#' @export
plot.mm_ranef <- function(x, level = 0.95, type = c("caterpillar", "qq"),
                          ...) {
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    mm_abort(
      message = "`level` must be a single number between 0 and 1.",
      class = "mm_arg_error",
      input = level
    )
  }
  panels <- mm_ranef_panels(x)
  if (!length(panels)) return(invisible(x))
  old <- graphics::par(mfrow = grDevices::n2mfrow(length(panels)))
  on.exit(graphics::par(old), add = TRUE)
  z <- stats::qnorm((1 + level) / 2)
  for (pn in panels) {
    if (identical(type, "qq")) {
      stats::qqnorm(pn$values, main = pn$title,
                    xlab = "Standard normal quantiles",
                    ylab = "Conditional modes", ...)
      stats::qqline(pn$values, lty = 2, col = "grey50")
      next
    }
    ord <- order(pn$values)
    v <- pn$values[ord]
    se <- pn$se[ord]
    lo <- v - z * se
    hi <- v + z * se
    xlim <- range(c(v, lo, hi), finite = TRUE)
    k <- length(v)
    graphics::plot(v, seq_len(k), xlim = xlim, yaxt = "n", pch = 16,
                   xlab = pn$coef, ylab = "", main = pn$title, ...)
    graphics::axis(2, at = seq_len(k), labels = pn$levels[ord], las = 1,
                   cex.axis = if (k > 20L) 0.6 else 0.8)
    graphics::abline(v = 0, lty = 2, col = "grey50")
    ok <- is.finite(lo) & is.finite(hi)
    if (any(ok)) {
      graphics::segments(lo[ok], seq_len(k)[ok], hi[ok], seq_len(k)[ok])
    }
  }
  invisible(x)
}

#' @rdname plot.mm_ranef
#' @method qqnorm mm_ranef
#' @export
qqnorm.mm_ranef <- function(y, ...) {
  plot.mm_ranef(y, type = "qq", ...)
  invisible(y)
}

# One panel per (grouping factor, coefficient): values, conditional SDs
# (NA when unavailable), level labels and a title.
mm_ranef_panels <- function(x) {
  out <- list()
  for (g in names(x)) {
    df <- x[[g]]
    pv <- attr(df, "postVar") %||% attr(df, "condVar")
    for (j in seq_along(df)) {
      se <- if (is.array(pv) && length(dim(pv)) == 3L &&
                dim(pv)[[1L]] >= j && dim(pv)[[3L]] == nrow(df)) {
        sqrt(pmax(pv[j, j, ], 0))
      } else {
        rep(NA_real_, nrow(df))
      }
      out[[length(out) + 1L]] <- list(
        values = as.numeric(df[[j]]),
        se = as.numeric(se),
        levels = rownames(df) %||% as.character(seq_len(nrow(df))),
        coef = names(df)[[j]],
        title = sprintf("%s: %s", g, names(df)[[j]])
      )
    }
  }
  out
}

# Long-format data for the lattice methods (lme4's asDf0 layout).
mm_ranef_long <- function(x, g) {
  df <- x[[g]]
  pv <- attr(df, "postVar") %||% attr(df, "condVar")
  ss <- utils::stack(df)
  ss$ind <- factor(as.character(ss$ind), levels = names(df))
  lv <- rownames(df) %||% as.character(seq_len(nrow(df)))
  ss$.nn <- rep.int(stats::reorder(factor(lv), df[[1L]], FUN = mean),
                    ncol(df))
  se <- unlist(lapply(seq_along(df), function(j) {
    if (is.array(pv) && length(dim(pv)) == 3L) sqrt(pmax(pv[j, j, ], 0)) else
      rep(NA_real_, nrow(df))
  }))
  ss$se <- if (all(is.na(se))) NULL else se
  ss
}

# lattice::dotplot() method for ranef() output (registered lazily).
dotplot.mm_ranef <- function(x, data, main = TRUE, level = 0.95, ...) {
  rng <- stats::qnorm((1 + level) / 2)
  prepanel_ci <- function(x, y, se, subscripts, ...) {
    if (is.null(se)) return(list())
    x <- as.numeric(x)
    hw <- rng * as.numeric(se[subscripts])
    list(xlim = range(x - hw, x + hw, finite = TRUE))
  }
  panel_ci <- function(x, y, se, subscripts, pch = 16, ...) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    lattice::panel.abline(h = unique(y), col = "grey90")
    lattice::panel.abline(v = 0, col = "grey50", lty = 2)
    if (!is.null(se)) {
      se <- as.numeric(se[subscripts])
      ok <- is.finite(se)
      lattice::panel.segments(x[ok] - rng * se[ok], y[ok],
                              x[ok] + rng * se[ok], y[ok], col = "black")
    }
    lattice::panel.xyplot(x, y, pch = pch, ...)
  }
  f <- function(g) {
    ss <- mm_ranef_long(x, g)
    lattice::dotplot(.nn ~ values | ind, data = ss, se = ss$se,
                     prepanel = prepanel_ci, panel = panel_ci, xlab = NULL,
                     main = if (isTRUE(main)) g, ...)
  }
  stats::setNames(lapply(names(x), f), names(x))
}

# lattice::qqmath() method for ranef() output (registered lazily).
qqmath.mm_ranef <- function(x, data, main = TRUE, level = 0.95, ...) {
  rng <- stats::qnorm((1 + level) / 2)
  f <- function(g) {
    df <- x[[g]]
    ss <- mm_ranef_long(x, g)
    nr <- nrow(df)
    nc <- ncol(df)
    ord <- unlist(lapply(df, order)) + rep((0:(nc - 1L)) * nr, each = nr)
    ind <- gl(nc, nr, labels = names(df))
    q <- rep(stats::qnorm((seq_len(nr) - 0.5) / nr), nc)
    se <- if (is.null(ss$se)) rep(NA_real_, nr * nc) else ss$se
    panel_ci <- function(x, y, se, subscripts, pch = 16, ...) {
      lattice::panel.grid(h = -1, v = -1)
      lattice::panel.abline(v = 0)
      se <- as.numeric(se[subscripts])
      ok <- is.finite(se)
      lattice::panel.segments(x[ok] - rng * se[ok], y[ok],
                              x[ok] + rng * se[ok], y[ok], col = "black")
      lattice::panel.xyplot(x, y, pch = pch, ...)
    }
    lattice::xyplot(q ~ unlist(df)[ord] | ind[ord], se = se[ord],
                    panel = panel_ci,
                    scales = list(x = list(relation = "free")),
                    ylab = "Standard normal quantiles", xlab = NULL,
                    main = if (isTRUE(main)) g, ...)
  }
  stats::setNames(lapply(names(x), f), names(x))
}
