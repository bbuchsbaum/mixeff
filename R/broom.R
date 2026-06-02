# broom / broom.mixed support for mixeff fits.
#
# tidy(), glance(), and augment() are generics owned by the `generics` package
# and re-exported by `broom` / `broom.mixed`. We define methods here and
# register them into the `generics` namespace on load (see zzz.R), mirroring how
# this package registers fixef/ranef/VarCorr into lme4 -- so `broom.mixed::tidy`
# (which IS `generics::tidy`) dispatches to these for mm_lmm / mm_glmm without a
# hard dependency on the broom stack.
#
# Output shapes follow broom.mixed::tidy.merMod conventions:
#   tidy   -> one row per term; effect in {"fixed","ran_pars","ran_vals"};
#             ran_pars terms are sd__<term> / cor__<a>.<b> / sd__Observation.
#   glance -> nobs, sigma, logLik, AIC, BIC, deviance, df.residual.
#   augment-> the model frame plus .fitted and .resid.

#' Tidy, glance, and augment methods for mixeff fits
#'
#' These implement the \pkg{broom} / \pkg{broom.mixed} generics for
#' [`mm_lmm`][lmm] and [`mm_glmm`][glmm] fits, so `tidy()`, `glance()`, and
#' `augment()` work on mixeff models the same way they do on `lme4` fits.
#'
#' `tidy()` returns one row per model term. `effects = "fixed"` yields the
#' fixed-effect coefficients (`estimate`, `std.error`, `statistic`, and, for
#' GLMMs, a Wald `p.value`); `effects = "ran_pars"` yields the variance-
#' component standard deviations (`sd__<term>`), correlations
#' (`cor__<a>.<b>`), and the residual SD (`sd__Observation`); `effects =
#' "ran_vals"` yields the conditional modes. `glance()` returns a one-row
#' model-summary frame; `augment()` returns the model frame with `.fitted`
#' and `.resid` columns.
#'
#' These methods are registered with \pkg{generics} when the package is loaded;
#' call them via `broom::tidy()` / `broom.mixed::tidy()` etc.
#'
#' @param x,data A fitted `mm_lmm` or `mm_glmm` (and, for `augment()`, optional
#'   data to augment; defaults to the model frame).
#' @param effects Which terms to return: any of `"fixed"`, `"ran_pars"`,
#'   `"ran_vals"`.
#' @param conf.int Logical; add Wald `conf.low`/`conf.high` for fixed effects.
#' @param conf.level Confidence level for `conf.int`.
#' @param ... Unused; for generic compatibility.
#'
#' @return A data frame.
#'
#' @name mm_broom
#' @keywords internal
NULL

# ---- tidy ----------------------------------------------------------------

tidy.mm_lmm <- function(x, effects = c("ran_pars", "fixed"),
                        conf.int = FALSE, conf.level = 0.95, ...) {
  mm_tidy_impl(x, effects = effects, conf.int = conf.int,
               conf.level = conf.level, glmm = FALSE)
}

tidy.mm_glmm <- function(x, effects = c("ran_pars", "fixed"),
                         conf.int = FALSE, conf.level = 0.95, ...) {
  mm_tidy_impl(x, effects = effects, conf.int = conf.int,
               conf.level = conf.level, glmm = TRUE)
}

mm_tidy_impl <- function(x, effects, conf.int, conf.level, glmm) {
  effects <- match.arg(effects, c("fixed", "ran_pars", "ran_vals"),
                       several.ok = TRUE)
  parts <- list()
  if ("fixed" %in% effects) {
    parts$fixed <- mm_tidy_fixed(x, conf.int, conf.level, glmm)
  }
  if ("ran_pars" %in% effects) {
    parts$ran_pars <- mm_tidy_ran_pars(x)
  }
  if ("ran_vals" %in% effects) {
    parts$ran_vals <- mm_tidy_ran_vals(x)
  }
  # Preserve fixed -> ran_pars -> ran_vals ordering.
  ordered <- parts[intersect(c("fixed", "ran_pars", "ran_vals"), names(parts))]
  mm_rbind_fill(ordered)
}

mm_tidy_fixed <- function(x, conf.int, conf.level, glmm) {
  beta <- fixef(x)
  se <- x$std_errors
  if (is.null(se) || length(se) != length(beta)) {
    se <- rep(NA_real_, length(beta))
  }
  statistic <- beta / se
  out <- data.frame(
    effect = "fixed",
    group = NA_character_,
    term = names(beta),
    estimate = as.numeric(beta),
    std.error = as.numeric(se),
    statistic = as.numeric(statistic),
    stringsAsFactors = FALSE
  )
  if (glmm) {
    # GLMM fixed effects use an asymptotic Wald z-test (matching glmer's tidy).
    out$p.value <- 2 * stats::pnorm(-abs(statistic))
  }
  if (isTRUE(conf.int)) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- out$estimate - z * out$std.error
    out$conf.high <- out$estimate + z * out$std.error
  }
  out
}

mm_tidy_ran_pars <- function(x) {
  vc <- as.data.frame(VarCorr(x))
  if (!nrow(vc)) {
    return(data.frame(
      effect = character(0), group = character(0), term = character(0),
      estimate = numeric(0), stringsAsFactors = FALSE
    ))
  }
  term <- character(nrow(vc))
  for (i in seq_len(nrow(vc))) {
    if (identical(vc$grp[i], "Residual")) {
      term[i] <- "sd__Observation"
    } else if (is.na(vc$var2[i])) {
      term[i] <- paste0("sd__", vc$var1[i])
    } else {
      term[i] <- paste0("cor__", vc$var1[i], ".", vc$var2[i])
    }
  }
  data.frame(
    effect = "ran_pars",
    group = vc$grp,
    term = term,
    estimate = vc$sdcor,
    stringsAsFactors = FALSE
  )
}

mm_tidy_ran_vals <- function(x) {
  rv <- as.data.frame(ranef(x, condVar = TRUE))
  if (!nrow(rv)) {
    return(data.frame(
      effect = character(0), group = character(0), level = character(0),
      term = character(0), estimate = numeric(0), std.error = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    effect = "ran_vals",
    group = rv$grpvar,
    level = rv$grp,
    term = rv$term,
    estimate = rv$condval,
    std.error = rv$condsd,
    stringsAsFactors = FALSE
  )
}

# rbind a list of data frames with differing columns, filling missing columns
# with type-appropriate NA. Character columns get NA_character_; everything
# else NA_real_.
mm_rbind_fill <- function(dfs) {
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  dfs <- dfs[vapply(dfs, function(d) nrow(d) > 0L, logical(1))]
  if (!length(dfs)) {
    return(data.frame(
      effect = character(0), group = character(0), term = character(0),
      estimate = numeric(0), stringsAsFactors = FALSE
    ))
  }
  char_cols <- c("effect", "group", "term", "level")
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  filled <- lapply(dfs, function(d) {
    for (m in setdiff(cols, names(d))) {
      d[[m]] <- if (m %in% char_cols) NA_character_ else NA_real_
    }
    d[cols]
  })
  out <- do.call(rbind, filled)
  rownames(out) <- NULL
  out
}

# ---- glance --------------------------------------------------------------

glance.mm_lmm <- function(x, ...) {
  mm_glance_impl(x)
}

glance.mm_glmm <- function(x, ...) {
  mm_glance_impl(x)
}

mm_glance_impl <- function(x) {
  data.frame(
    nobs = nobs(x),
    sigma = tryCatch(as.numeric(sigma(x)), error = function(e) NA_real_),
    logLik = as.numeric(x$logLik),
    AIC = as.numeric(x$AIC),
    BIC = as.numeric(x$BIC),
    deviance = as.numeric(x$deviance),
    df.residual = df.residual(x),
    stringsAsFactors = FALSE
  )
}

# ---- augment -------------------------------------------------------------

augment.mm_lmm <- function(x, data = stats::model.frame(x), ...) {
  mm_augment_impl(x, data)
}

augment.mm_glmm <- function(x, data = stats::model.frame(x), ...) {
  mm_augment_impl(x, data)
}

mm_augment_impl <- function(x, data) {
  out <- as.data.frame(data, stringsAsFactors = FALSE)
  fitted_vals <- tryCatch(as.numeric(fitted(x)), error = function(e) NULL)
  resid_vals <- tryCatch(as.numeric(stats::residuals(x)),
                         error = function(e) NULL)
  if (!is.null(fitted_vals) && length(fitted_vals) == nrow(out)) {
    out$.fitted <- fitted_vals
  }
  if (!is.null(resid_vals) && length(resid_vals) == nrow(out)) {
    out$.resid <- resid_vals
  }
  out
}
