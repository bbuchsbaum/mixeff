#' Predict from a fitted mixeff LMM
#'
#' Predictions follow the lme4 generic shape. In-sample predictions reuse
#' the cached fitted/fixed values; new-data predictions are dispatched
#' through the Rust `predict_new` contract.
#'
#' @param object A fitted `mm_lmm` object.
#' @param newdata Optional new data. Must be a `data.frame` containing every
#'   variable the requested prediction needs (the fixed-effect variables, plus
#'   the grouping and slope variables of the random terms used); other
#'   columns are ignored. Grouping columns may be factor, character, or
#'   integer (coerced the way they were at fit time); fixed-effect factor and
#'   character columns are matched to the training levels, and a level not
#'   seen when fitting is refused. Stateful terms (`poly()`, `scale()`,
#'   `ns()`, ...) are evaluated with the training basis, as in [stats::lm()].
#' @param re.form Random-effects conditioning, following lme4. `NULL`
#'   returns conditional predictions; `NA` (or `~0`) returns population-level
#'   (fixed-effect) predictions; a one-sided formula naming some of the
#'   model's random terms (e.g. `~ (1 | subject)`) conditions on those terms'
#'   conditional modes only. Standard errors and intervals are not available
#'   for such partial conditioning.
#' @param allow.new.levels When `FALSE` (default), unseen grouping levels in
#'   `newdata` raise `mm_inference_unavailable` through the Rust
#'   `NewReLevels::Error` policy. When `TRUE`, unseen levels are replaced by
#'   the population mean (zero random effect), matching
#'   `lme4::predict(allow.new.levels = TRUE)`.
#' @param type Prediction scale. Gaussian LMMs use the same values for
#'   `"response"` and `"link"`. For `residuals()`, the residual type, as in
#'   lme4: LMMs default to `"response"` (`y - mu`, also `"working"`), with
#'   `"pearson"`/`"deviance"` the weighted residuals `sqrt(w) * (y - mu)`;
#'   GLMMs default to `"deviance"` and also offer `"pearson"`
#'   (`(y - mu) * sqrt(w / V(mu))`), `"working"` and `"response"`.
#' @param se.fit Logical; when `TRUE`, returns a list with `fit` and `se.fit`.
#'   For population predictions (`re.form = NA`) the standard error is the
#'   Wald SE of the fixed-effect linear predictor, `sqrt(diag(X V X'))`. For
#'   conditional predictions (`re.form = NULL`) the SE comes from the engine
#'   prediction-variance payload, which adds the random-effect contribution
#'   (BLUP variance and the fixed/random covariance). Rows the engine cannot
#'   certify — e.g. unseen grouping levels under `allow.new.levels = TRUE` —
#'   return `NA` with the engine's reason in the `mm_reason` attribute.
#'   (`lme4::predict.merMod` offers no conditional SE at all.)
#' @param interval Interval type: `"confidence"` for the fitted mean or
#'   `"prediction"` for a new observation (adds the residual variance).
#'   Population (`re.form = NA`) intervals are `fit +/- z*se` computed R-side;
#'   conditional (`re.form = NULL`) bounds come from the engine
#'   prediction-variance payload. Returns a matrix with `fit`/`lwr`/`upr`.
#' @param level Confidence level for `interval` / `se.fit` intervals.
#' @param na.action Missing-value handling for `newdata`, as in
#'   `lme4::predict.merMod()`: the default [stats::na.pass] returns `NA` for
#'   rows with a missing value in a variable the prediction needs;
#'   [stats::na.omit] drops them and [stats::na.exclude] pads them back as
#'   `NA`. In-sample predictions of a fit made with `na.action = na.exclude`
#'   are padded to the original rows, like [fitted()] and [residuals()].
#' @param offset Offset for `newdata` rows when the model was fitted with the
#'   `offset =` argument: a numeric vector (one value per `newdata` row) or an
#'   expression evaluated in `newdata`. Offsets written in the formula
#'   (`offset(log(t))`) are evaluated from `newdata` automatically. (lme4
#'   silently drops an `offset =` argument offset for new data; mixeff refuses
#'   instead of guessing.)
#' @param scaled Logical; when `TRUE`, residuals are divided by `sigma(object)`
#'   (once), as in lme4.
#' @param ... Reserved for generic compatibility.
#'
#' @return A numeric vector, or a list with `fit` and `se.fit` when
#'   `se.fit = TRUE`.
#'
#' @export
predict.mm_lmm <- function(object,
                           newdata = NULL,
                           re.form = NULL,
                           allow.new.levels = FALSE,
                           type = c("response", "link"),
                           se.fit = FALSE,
                           interval = c("none", "confidence", "prediction"),
                           level = 0.95,
                           na.action = stats::na.pass,
                           offset = NULL,
                           ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  mm_reject_unsupported_dots(
    list(...), "predict",
    c(random.only = "use `re.form = NA` for population-level (fixed-only) predictions or `re.form = NULL` for conditional predictions.")
  )
  mm_validate_allow_new_levels(allow.new.levels)
  rf <- mm_resolve_re_form(object, re.form)
  want_se <- isTRUE(se.fit) || !identical(interval, "none")
  mm_refuse_partial_se(rf, want_se)
  offset_new <- if (!is.null(newdata) && is.data.frame(newdata)) {
    eval(substitute(offset), newdata, parent.frame())
  } else {
    offset
  }

  if (is.null(newdata)) {
    if (!is.null(offset_new)) mm_refuse_insample_offset()
    res <- mm_lmm_predict_rows(object, NULL, rf, allow.new.levels, want_se,
                               se.fit, interval, level)
    return(mm_prediction_map(res, function(x) mm_napredict(object, x)))
  }
  nd <- mm_prepare_newdata(object, newdata, rf, na.action, offset_new)
  res <- mm_lmm_predict_rows(object, nd$data, rf, allow.new.levels, want_se,
                             se.fit, interval, level)
  mm_prediction_map(res, function(x) mm_newdata_fill(x, nd))
}

# Rows of `data` (NULL = the training frame) are complete and normalised.
mm_lmm_predict_rows <- function(object, data, rf, allow_new_levels, want_se,
                                se.fit, interval, level) {
  target <- rf$target
  insample <- is.null(data)
  frame <- if (insample) object$model_frame else data
  nm <- rownames(frame)
  if (!nrow(frame)) {
    pred <- stats::setNames(numeric(), character())
  } else if (insample) {
    pred <- switch(
      target,
      conditional = object$fitted,
      population  = object$fixed_fitted,
      partial     = object$fixed_fitted +
        mm_re_eta(object, frame, rf$formula, allow_new_levels)
    )
  } else {
    off <- mm_newdata_offset(object, data)
    pred <- switch(
      target,
      conditional = mm_predict_conditional_newdata(object, data,
                                                   allow_new_levels) + off,
      population  = mm_predict_fixed_only(object, data) + off,
      partial     = mm_predict_fixed_only(object, data) + off +
        mm_re_eta(object, data, rf$formula, allow_new_levels)
    )
  }
  pred <- as.numeric(pred)
  names(pred) <- nm
  if (!want_se) {
    return(pred)
  }
  if (!nrow(frame)) {
    mm_abort(
      message = "`newdata` has no complete rows to compute standard errors for.",
      class = "mm_data_error"
    )
  }
  off <- if (insample) object$offset %||% 0 else mm_newdata_offset(object, data)

  # Population (fixed-effect-only) SEs/intervals are computed R-side below from
  # sqrt(diag(X V X')). Conditional (`re.form = NULL`) SEs and intervals come
  # from the engine prediction-variance payload, which adds the random-effect
  # contribution (per-row se_fit and confidence/prediction bounds; new grouping
  # levels return NA with a reason).
  if (!identical(target, "population")) {
    pv <- mm_lmm_prediction_variance(object, frame, allow_new_levels, level)
    se <- pv$se_fit
    names(se) <- names(pred)
    if (anyNA(se)) attr(se, "mm_reason") <- pv$reason
    if (!identical(interval, "none")) {
      lwr <- if (identical(interval, "prediction")) pv$prediction_lower else pv$confidence_lower
      upr <- if (identical(interval, "prediction")) pv$prediction_upper else pv$confidence_upper
      # The engine works on the offset-adjusted response scale.
      out <- cbind(fit = pred, lwr = lwr + off, upr = upr + off)
      rownames(out) <- names(pred)
      attr(out, "interval") <- interval
      attr(out, "level") <- level
      if (anyNA(out)) attr(out, "mm_reason") <- pv$reason
      if (isTRUE(se.fit)) {
        return(list(fit = out, se.fit = se))
      }
      return(out)
    }
    return(list(fit = pred, se.fit = se))
  }

  se <- mm_fixed_prediction_se(object, frame)
  if (!identical(interval, "none")) {
    crit <- stats::qnorm((1 + level) / 2)
    extra <- if (identical(interval, "prediction")) object$sigma^2 else 0
    half <- crit * sqrt(se^2 + extra)
    out <- cbind(fit = pred, lwr = pred - half, upr = pred + half)
    rownames(out) <- names(pred)
    attr(out, "interval") <- interval
    attr(out, "level") <- level
    if (isTRUE(se.fit)) {
      return(list(fit = out, se.fit = se))
    }
    return(out)
  }
  names(se) <- names(pred)
  list(fit = pred, se.fit = se)
}

mm_validate_allow_new_levels <- function(allow.new.levels) {
  if (!is.logical(allow.new.levels) || length(allow.new.levels) != 1L ||
      is.na(allow.new.levels)) {
    mm_abort(
      message = "`allow.new.levels` must be TRUE or FALSE.",
      class = "mm_arg_error",
      input = allow.new.levels
    )
  }
  invisible(TRUE)
}

mm_refuse_partial_se <- function(rf, want_se) {
  if (identical(rf$target, "partial") && isTRUE(want_se)) {
    mm_abort(
      message = paste(
        "Standard errors and intervals are available for `re.form = NULL`",
        "(conditional) or `re.form = NA` (population) predictions only, not",
        "for conditioning on a subset of the random terms."
      ),
      class = "mm_inference_unavailable",
      reason_code = "partial_re_form_se_unavailable"
    )
  }
  invisible(TRUE)
}

mm_refuse_insample_offset <- function() {
  mm_abort(
    message = paste(
      "`offset` applies to `newdata` rows only; in-sample predictions use the",
      "offset the model was fitted with."
    ),
    class = "mm_arg_error"
  )
}

# Apply `f` to every per-row component of a prediction result (vector,
# interval matrix, or se.fit list), keeping attributes such as mm_reason.
mm_prediction_map <- function(res, f) {
  apply_one <- function(x) {
    if (is.null(x)) return(x)
    keep <- attributes(x)[intersect(names(attributes(x)),
                                    c("mm_reason", "interval", "level"))]
    out <- f(x)
    for (a in names(keep)) {
      if (identical(a, "mm_reason")) {
        r <- keep[[a]]
        # Reasons are per row: map them alongside the values.
        attr(out, a) <- if (length(r) == NROW(x)) f(r) else r
      } else {
        attr(out, a) <- keep[[a]]
      }
    }
    out
  }
  if (is.list(res) && !is.data.frame(res)) {
    lapply(res, apply_one)
  } else {
    apply_one(res)
  }
}

# Spread values computed for the complete rows of a prepared newdata back over
# the output rows (NA where a needed variable was missing), then pad for
# na.exclude.
mm_newdata_fill <- function(x, nd) {
  n <- length(nd$out_names)
  if (is.matrix(x)) {
    out <- matrix(NA_real_, nrow = n, ncol = ncol(x),
                  dimnames = list(nd$out_names, colnames(x)))
    out[nd$ok, ] <- x
  } else if (is.character(x)) {
    out <- rep(NA_character_, n)
    out[nd$ok] <- x
    names(out) <- nd$out_names
  } else {
    out <- rep(NA_real_, n)
    out[nd$ok] <- as.numeric(x)
    names(out) <- nd$out_names
  }
  if (!is.null(nd$na_action)) out <- stats::napredict(nd$na_action, out)
  out
}

# Resolve lme4's re.form into a prediction target: "conditional" (all random
# terms), "population" (none), or "partial" (the bars of a one-sided formula
# that name a strict subset of the model's random terms).
mm_resolve_re_form <- function(fit, re.form) {
  target <- mm_prediction_target(re.form)
  if (!identical(target, "unsupported")) {
    return(list(target = target, formula = NULL))
  }
  if (!inherits(re.form, "formula")) {
    mm_abort(
      message = "`re.form` must be NULL, NA, ~0, or a one-sided formula of random-effect terms.",
      class = "mm_arg_error",
      input = re.form
    )
  }
  bars <- mm_find_bars(re.form[[length(re.form)]])
  if (!length(bars)) {
    mm_abort(
      message = paste(
        "`re.form` must name random-effect terms such as `~ (1 | g)`; use",
        "`re.form = NA` for population-level predictions."
      ),
      class = "mm_arg_error",
      input = re.form
    )
  }
  key <- function(specs) {
    vapply(specs, function(s) paste(deparse1(s$lhs), s$label, sep = "|"),
           character(1))
  }
  wanted <- mm_lme4_group_specs(re.form)
  model <- mm_lme4_group_specs(fit$formula)
  if (setequal(key(wanted), key(model))) {
    return(list(target = "conditional", formula = NULL))
  }
  list(target = "partial", formula = re.form)
}

# Formula holding only the random terms of `formula` (lme4::reOnly()).
mm_reonly_formula <- function(formula) {
  bars <- mm_find_bars(formula[[length(formula)]])
  stats::as.formula(
    paste("~", paste0("(", vapply(bars, deparse1, character(1)), ")",
                      collapse = " + ")),
    env = environment(formula)
  )
}

# Validate and normalise `newdata` for prediction (lme4 semantics): keep only
# the variables the requested prediction needs, coerce grouping and factor
# columns the way the training frame was coerced, and apply `na.action`.
# Returns the complete rows to predict plus the bookkeeping that maps them
# back to the output rows.
mm_prepare_newdata <- function(fit, newdata, rf, na.action, offset_new = NULL) {
  if (!is.data.frame(newdata)) {
    mm_abort(
      message = "`newdata` must be a data.frame.",
      class = "mm_data_error",
      input = newdata
    )
  }
  fixed_vars <- all.vars(mm_user_fixed_terms(fit))
  re_vars <- switch(
    rf$target,
    population = character(),
    conditional = all.vars(mm_reonly_formula(fit$formula)),
    partial = all.vars(rf$formula)
  )
  needed <- unique(c(fixed_vars, re_vars))
  # Variables from the formula environment (e.g. `k` in poly(x, k)) are not
  # data columns.
  needed <- intersect(needed, c(names(fit$model_frame), re_vars))
  missing_vars <- setdiff(needed, names(newdata))
  if (length(missing_vars)) {
    mm_abort(
      message = sprintf(
        "`newdata` is missing variable(s) required by the model formula: %s.",
        paste(missing_vars, collapse = ", ")
      ),
      class = "mm_data_error",
      input = newdata,
      missing = missing_vars
    )
  }
  nd <- newdata[, needed, drop = FALSE]

  grouping <- mm_formula_grouping_vars(fit$formula)
  for (v in needed) {
    train <- fit$model_frame[[v]]
    if (!is.factor(train)) next
    chr <- as.character(nd[[v]])
    lv <- levels(train)
    if (v %in% grouping && !(v %in% fixed_vars)) {
      new_lv <- setdiff(unique(chr[!is.na(chr)]), lv)
      nd[[v]] <- factor(chr, levels = c(lv, new_lv))
    } else {
      bad <- !is.na(chr) & !(chr %in% lv)
      if (any(bad)) {
        mm_abort(
          message = sprintf(
            "Factor `%s` has level(s) in `newdata` not seen when fitting: %s.",
            v, paste(unique(chr[bad]), collapse = ", ")
          ),
          class = "mm_data_error",
          input = unique(chr[bad])
        )
      }
      nd[[v]] <- factor(chr, levels = lv, ordered = is.ordered(train))
    }
  }

  offset_arg_fit <- mm_fit_offset_arg(fit)
  if (!is.null(offset_new)) {
    if (is.null(offset_arg_fit)) {
      mm_abort(
        message = paste(
          "`offset` was supplied, but the model was not fitted with an",
          "`offset =` argument; formula offset() terms are taken from",
          "`newdata` automatically."
        ),
        class = "mm_arg_error"
      )
    }
    if (!is.numeric(offset_new) || length(offset_new) != nrow(newdata)) {
      mm_abort(
        message = sprintf(
          "`offset` must be numeric with one value per row of `newdata` (%d).",
          nrow(newdata)
        ),
        class = "mm_arg_error",
        input = offset_new
      )
    }
    nd[[".mm_offset_arg"]] <- as.numeric(offset_new)
  } else if (!is.null(offset_arg_fit)) {
    mm_abort(
      message = paste(
        "This model was fitted with an `offset =` argument, which cannot be",
        "reconstructed for `newdata`. Supply it with predict(..., offset = ),",
        "or write the offset in the formula (`+ offset(log(t))`) so it is",
        "evaluated from `newdata`. (lme4 silently predicts new data with a",
        "zero offset in this case.)"
      ),
      class = "mm_inference_unavailable",
      reason_code = "newdata_offset_unavailable",
      input = "offset"
    )
  }

  if (!is.function(na.action)) na.action <- match.fun(na.action)
  out_idx <- seq_len(nrow(nd))
  na_out <- attr(na.action(nd), "na.action")
  if (!is.null(na_out)) out_idx <- out_idx[-as.integer(na_out)]
  kept <- nd[out_idx, , drop = FALSE]
  ok <- stats::complete.cases(kept)
  list(
    data = kept[ok, , drop = FALSE],
    out_names = rownames(nd)[out_idx],
    ok = ok,
    na_action = na_out
  )
}

# The `offset =` part of a fit's offset (NULL when none). Fits saved before
# formula offsets existed stored only the argument offset in `$offset`.
mm_fit_offset_arg <- function(fit) {
  if ("offset_arg" %in% names(fit)) fit$offset_arg else fit$offset
}

# Total offset for prepared newdata rows: formula offset() terms evaluated in
# the data plus the supplied `offset =` values.
mm_newdata_offset <- function(fit, data) {
  off <- mm_formula_offset(fit, data)
  arg <- data[[".mm_offset_arg"]]
  if (!is.null(arg)) off <- if (is.null(off)) arg else off + arg
  if (is.null(off)) 0 else as.numeric(off)
}

# Newdata restricted to (and augmented with) the columns the engine formula
# uses: expanded synthetic columns are rebuilt with the training basis.
mm_engine_newdata <- function(fit, data) {
  data <- mm_expand_newdata(fit, data)
  vars <- setdiff(intersect(all.vars(mm_engine_formula(fit)), names(data)),
                  all.vars(mm_engine_formula(fit)[[2L]]))
  data[, vars, drop = FALSE]
}

# Standard error of the population (fixed-effect-only) linear predictor for the
# rows of `data`: sqrt(diag(X V X')) with X the engine-basis fixed-effect design
# and V the stored fixed-effect covariance. Used by predict() for se.fit and
# population (re.form = NA) confidence/prediction intervals.
mm_fixed_prediction_se <- function(fit, data) {
  X <- mm_engine_fixed_matrix(fit, data)
  V <- as.matrix(unclass(fit$fixed_effect_vcov))
  dimnames(V) <- list(names(fit$beta), names(fit$beta))
  V <- V[colnames(X), colnames(X), drop = FALSE]
  sqrt(pmax(0, rowSums((X %*% V) * X)))
}

#' @rdname predict.mm_lmm
#' @export
fitted.mm_lmm <- function(object, ...) {
  out <- object$fitted
  names(out) <- rownames(object$model_frame)
  mm_napredict(object, out)
}

#' @rdname predict.mm_lmm
#' @export
residuals.mm_lmm <- function(object,
                             type = c("response", "pearson", "deviance",
                                      "working"),
                             scaled = FALSE, ...) {
  type <- match.arg(type)
  out <- object$residuals
  # lme4 (residuals.lmResp): "working"/"response" are y - mu; "pearson" and
  # "deviance" are the weighted residuals sqrt(w) * (y - mu). Only
  # `scaled = TRUE` divides by sigma, exactly once.
  if (type %in% c("pearson", "deviance")) {
    out <- out * sqrt(mm_prior_weights(object))
  }
  if (isTRUE(scaled)) {
    out <- out / sigma(object)
  }
  names(out) <- rownames(object$model_frame)
  mm_naresid(object, out)
}

#' @rdname predict.mm_lmm
#' @export
fitted.mm_glmm <- fitted.mm_lmm

#' @rdname predict.mm_lmm
#' @export
residuals.mm_glmm <- function(object,
                              type = c("deviance", "pearson", "working",
                                       "response"),
                              scaled = FALSE, ...) {
  # lme4 (residuals.glmResp): deviance residuals by default, computed from
  # the family at the conditional mean mu, the response y (a proportion for
  # binomial fits) and the prior weights (trial counts).
  type <- match.arg(type)
  out <- mm_glmm_residuals(object, type)
  if (isTRUE(scaled)) {
    out <- out / sigma(object)
  }
  names(out) <- rownames(object$model_frame)
  mm_naresid(object, out)
}

#' Predict from a fitted mixeff GLMM
#'
#' GLMM predictions are computed on the R side from the stored fixed effects
#' (population, `re.form = NA`) or fixed effects plus conditional modes
#' (`re.form = NULL`, or a subset of the random terms via a one-sided
#' `re.form` formula), then mapped through the family link. This mirrors
#' `lme4::predict.merMod` for generalized models: `type = "link"` returns the
#' linear predictor and `type = "response"` the mean. Offsets written in the
#' formula (`offset(log(t))`) are evaluated from `newdata`; an `offset =`
#' argument offset must be supplied again through `predict(offset = )`.
#'
#' In-sample response predictions reuse the engine's certified fitted values.
#'
#' Standard errors and confidence intervals: population (`re.form = NA`)
#' SEs are the fixed-effect Wald SE mapped through the link by the delta
#' method; conditional (`re.form = NULL`) SEs and confidence bounds come from
#' the engine prediction-variance payload. The engine certifies these rows
#' for `method = "joint_laplace"` fits and for `pirls_profiled` fits
#' whose post-fit profiled-optimum certificate is issued (per-row status
#' `"available"`). Uncertified fits (e.g. singular fits, or fits whose
#' certificate fails) keep status `"degraded"`, and their conditional SEs and
#' bounds are withheld as `NA` with the engine's reason in the `mm_reason`
#' attribute — consistent with the package's "no fake certainty" contract.
#' For a fit with an offset, conditional SEs/intervals on `newdata` are
#' available on the link scale only.
#'
#' Prediction (future-observation) intervals (`interval = "prediction"`) are
#' available for conditional, response-scale predictions: the engine returns
#' quantiles of the plug-in predictive distribution (the family conditional
#' distribution mixed over link-scale fitted-mean uncertainty), so bounds are
#' integers for count families and support points for Bernoulli. They are
#' refused with a typed condition on the link scale (future observations are
#' response-scale objects), for population-level requests, and for grouped
#' binomial fits (the future trial count is not representable in `newdata`).
#'
#' @inheritParams predict.mm_lmm
#' @param object A fitted `mm_glmm` object.
#' @return A numeric vector, or a list with `fit` and `se.fit` when
#'   `se.fit = TRUE`.
#' @method predict mm_glmm
#' @export
predict.mm_glmm <- function(object,
                            newdata = NULL,
                            re.form = NULL,
                            allow.new.levels = FALSE,
                            type = c("response", "link"),
                            se.fit = FALSE,
                            interval = c("none", "confidence", "prediction"),
                            level = 0.95,
                            na.action = stats::na.pass,
                            offset = NULL,
                            ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  mm_reject_unsupported_dots(
    list(...), "predict",
    c(random.only = "use `re.form = NA` for population-level (fixed-only) predictions or `re.form = NULL` for conditional predictions.")
  )
  mm_validate_allow_new_levels(allow.new.levels)
  rf <- mm_resolve_re_form(object, re.form)
  target <- rf$target
  want_se <- isTRUE(se.fit) || !identical(interval, "none")
  mm_refuse_partial_se(rf, want_se)
  if (identical(interval, "prediction")) {
    if (identical(type, "link")) {
      mm_abort(
        message = paste(
          "GLMM prediction (future-observation) intervals are response-scale",
          "objects; the engine refuses them on the link scale. Use",
          "`type = \"response\"`, or `interval = \"confidence\"` for",
          "fitted-mean intervals on the link scale."
        ),
        class = "mm_inference_unavailable",
        input = interval
      )
    }
    if (!identical(target, "conditional")) {
      mm_abort(
        message = paste(
          "GLMM prediction (future-observation) intervals are only available",
          "for conditional predictions (`re.form = NULL`); the population",
          "Wald path carries no family variance term."
        ),
        class = "mm_inference_unavailable",
        input = re.form
      )
    }
    if (identical(mm_glmm_engine_family(object), "binomial")) {
      mm_abort(
        message = paste(
          "GLMM prediction intervals are unavailable for grouped binomial",
          "fits: the future trial count is not representable in `newdata`,",
          "so the engine refuses the rows."
        ),
        class = "mm_inference_unavailable",
        input = interval
      )
    }
  }
  offset_new <- if (!is.null(newdata) && is.data.frame(newdata)) {
    eval(substitute(offset), newdata, parent.frame())
  } else {
    offset
  }
  if (is.null(newdata)) {
    if (!is.null(offset_new)) mm_refuse_insample_offset()
    res <- mm_glmm_predict_rows(object, NULL, rf, allow.new.levels, type,
                                want_se, se.fit, interval, level)
    return(mm_prediction_map(res, function(x) mm_napredict(object, x)))
  }
  nd <- mm_prepare_newdata(object, newdata, rf, na.action, offset_new)
  if (want_se && !identical(target, "population") &&
      !is.null(object$offset) && identical(type, "response")) {
    mm_abort(
      message = paste(
        "Conditional standard errors and intervals on `newdata` for a GLMM",
        "with an offset are available on the link scale only",
        "(`type = \"link\"`): the engine's response-scale variance is computed",
        "without the new rows' offset."
      ),
      class = "mm_inference_unavailable",
      reason_code = "newdata_offset_response_se_unavailable"
    )
  }
  res <- mm_glmm_predict_rows(object, nd$data, rf, allow.new.levels, type,
                              want_se, se.fit, interval, level)
  mm_prediction_map(res, function(x) mm_newdata_fill(x, nd))
}

mm_glmm_predict_rows <- function(object, data, rf, allow_new_levels, type,
                                 want_se, se.fit, interval, level) {
  target <- rf$target
  family <- mm_glmm_family_from_info(object$family)
  insample <- is.null(data)
  frame <- if (insample) object$model_frame else data
  nm <- rownames(frame)

  if (!nrow(frame)) {
    eta <- mu <- numeric()
  } else if (insample && identical(target, "conditional")) {
    # Engine-certified conditional mean; eta via the link.
    mu <- object$fitted
    eta <- family$linkfun(mu)
  } else {
    off <- if (insample) object$offset %||% 0 else mm_newdata_offset(object, data)
    eta <- as.numeric(mm_predict_fixed_only(object, frame)) + off
    if (identical(target, "conditional")) {
      eta <- eta + mm_re_eta(object, frame, mm_reonly_formula(object$formula),
                             allow_new_levels,
                             missing_class = "mm_inference_unavailable")
    } else if (identical(target, "partial")) {
      eta <- eta + mm_re_eta(object, frame, rf$formula, allow_new_levels)
    }
    mu <- family$linkinv(eta)
  }

  pred <- if (identical(type, "response")) as.numeric(mu) else as.numeric(eta)
  names(pred) <- nm
  if (!want_se) {
    return(pred)
  }
  if (!nrow(frame)) {
    mm_abort(
      message = "`newdata` has no complete rows to compute standard errors for.",
      class = "mm_data_error"
    )
  }

  if (identical(target, "population")) {
    # Population (fixed-only) Wald SE of the linear predictor; for
    # type = "response" map through the link with the delta method. Confidence
    # intervals are estimate +/- z * SE on the requested scale.
    se_link <- mm_fixed_prediction_se(object, frame)
    se <- if (identical(type, "response")) se_link * abs(family$mu.eta(eta)) else se_link
    names(se) <- nm
    if (!identical(interval, "none")) {
      crit <- stats::qnorm((1 + level) / 2)
      out <- cbind(fit = pred, lwr = pred - crit * se, upr = pred + crit * se)
      rownames(out) <- nm
      attr(out, "interval") <- interval
      attr(out, "level") <- level
      if (isTRUE(se.fit)) return(list(fit = out, se.fit = se))
      return(out)
    }
    return(list(fit = pred, se.fit = se))
  }

  # Conditional (re.form = NULL): engine prediction-variance payload on the
  # requested scale (the engine propagates variance through the link). Certified
  # (joint_laplace) fits return finite se_fit / confidence bounds; uncertified
  # fast-PIRLS or new-grouping-level rows return NA with a reason. On newdata
  # the engine sees no offset, so link-scale bounds are shifted by it here.
  pv <- mm_glmm_prediction_variance(object, frame, type, allow_new_levels, level)
  shift <- if (!insample && identical(type, "link")) mm_newdata_offset(object, frame) else 0
  se <- pv$se_fit
  names(se) <- nm
  if (anyNA(se)) attr(se, "mm_reason") <- pv$reason
  if (!identical(interval, "none")) {
    lwr <- if (identical(interval, "prediction")) pv$prediction_lower else pv$confidence_lower
    upr <- if (identical(interval, "prediction")) pv$prediction_upper else pv$confidence_upper
    out <- cbind(fit = pred, lwr = lwr + shift, upr = upr + shift)
    rownames(out) <- nm
    attr(out, "interval") <- interval
    attr(out, "level") <- level
    if (anyNA(out)) attr(out, "mm_reason") <- pv$reason
    if (isTRUE(se.fit)) return(list(fit = out, se.fit = se))
    return(out)
  }
  list(fit = pred, se.fit = se)
}

# Random-effect contribution to the GLMM linear predictor for arbitrary `data`
# from the stored conditional modes (all of the model's random terms).
mm_glmm_re_eta <- function(fit, data, allow_new_levels) {
  mm_re_eta(fit, data, mm_reonly_formula(fit$formula), allow_new_levels,
            missing_class = "mm_inference_unavailable")
}

mm_prediction_target <- function(re.form) {
  if (is.null(re.form)) {
    return("conditional")
  }
  if (length(re.form) == 1L && is.na(re.form)) {
    return("population")
  }
  # ~0 / ~ 0 evaluates to a one-sided formula with character "~0"; treat as
  # explicit population request.
  if (inherits(re.form, "formula") &&
      identical(length(re.form), 2L) &&
      identical(trimws(deparse1(re.form[[2L]])), "0")) {
    return("population")
  }
  "unsupported"
}


# Conditional newdata prediction via the Rust `predict_new` FFI. The Rust
# contract refits the model from formula + training data, so the call is
# self-contained.
mm_predict_conditional_newdata <- function(fit, newdata, allow_new_levels) {
  policy <- if (isTRUE(allow_new_levels)) "population" else "error"

  spec_data <- mm_translate_data(mm_engine_frame(fit))
  formula_string <- mm_coerce_formula_string(mm_engine_formula(fit))
  control_json <- mm_refit_control_json(fit)
  new_data <- mm_translate_data(mm_engine_newdata(fit, newdata))

  json <- tryCatch(
    .Call(
      wrap__mm_lmm_predict_new_json,
      formula_string,
      isTRUE(fit$REML),
      spec_data$column_order,
      spec_data$numeric_columns,
      spec_data$categorical_values,
      spec_data$categorical_levels,
      spec_data$categorical_ordered,
      mm_bridge_weights(fit$weights),
      as.character(control_json),
      new_data$column_order,
      new_data$numeric_columns,
      new_data$categorical_values,
      new_data$categorical_levels,
      # Ordered-factor coding for newdata is intentionally empty: the engine
      # re-derives newdata contrasts from the fitted training snapshot
      # (align_newdata_to_training), so flagging ordered columns here is a
      # no-op that could only spuriously error on an ordered newdata column
      # observed at <2 levels. Training flags (above) drive the coding.
      character(0),
      policy
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, newdata = newdata)
  }

  payload <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse predict_new JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        parent = cnd
      )
    }
  )
  schema <- payload$schema %||% list()
  if (!identical(as.character(schema$schema_name), "mixeff.lmm_predict_new") ||
      !identical(as.character(schema$schema_version), "1")) {
    mm_abort(
      message = "predict_new payload has an unknown schema header.",
      class = "mm_schema_error",
      input = payload
    )
  }

  preds_raw <- payload$predictions %||% list()
  out <- vapply(preds_raw, function(v) {
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  if (length(out) != nrow(newdata)) {
    mm_abort(
      message = sprintf(
        "predict_new returned %d prediction(s) for %d input row(s).",
        length(out), nrow(newdata)
      ),
      class = "mm_schema_error",
      input = payload
    )
  }
  names(out) <- rownames(newdata)
  out
}

# Parse an engine `PredictionVariancePayload` (schema
# "mixedmodels.prediction_variance") into a per-row data.frame. Unavailable /
# degraded rows carry `se_fit = NA` (encoded as JSON null by the engine) plus a
# stable `reason`, preserving the no-fake-certainty contract.
mm_parse_prediction_variance <- function(json, n) {
  payload <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse prediction-variance JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error", parent = cnd
      )
    }
  )
  if (!identical(as.character(payload$schema_name), "mixedmodels.prediction_variance")) {
    mm_abort(
      message = "prediction-variance payload has an unknown schema header.",
      class = "mm_schema_error", input = payload
    )
  }
  rows <- payload$rows %||% list()
  if (length(rows) != n) {
    mm_abort(
      message = sprintf(
        "prediction-variance returned %d row(s) for %d input row(s).",
        length(rows), n
      ),
      class = "mm_schema_error", input = payload
    )
  }
  num <- function(field) {
    vapply(rows, function(r) {
      v <- r[[field]]
      if (is.null(v)) NA_real_ else as.numeric(v)
    }, numeric(1))
  }
  chr <- function(field) {
    vapply(rows, function(r) {
      v <- r[[field]]
      if (is.null(v)) NA_character_ else as.character(v)
    }, character(1))
  }
  out <- data.frame(
    se_fit = num("se_fit"),
    fixed_variance = num("fixed_variance"),
    confidence_lower = num("confidence_lower"),
    confidence_upper = num("confidence_upper"),
    prediction_lower = num("prediction_lower"),
    prediction_upper = num("prediction_upper"),
    status = chr("status"),
    reason = chr("reason"),
    stringsAsFactors = FALSE
  )
  # Only engine-certified rows cross the boundary with numbers. `degraded`
  # rows (e.g. the uncertified fast-PIRLS working-delta variance) carry finite
  # values the engine does not certify, so they are withheld here with the
  # engine's reason rather than reported — no fake certainty.
  masked <- !is.na(out$status) & out$status != "available"
  if (any(masked)) {
    numeric_fields <- c("se_fit", "fixed_variance", "confidence_lower",
                        "confidence_upper", "prediction_lower", "prediction_upper")
    out[masked, numeric_fields] <- NA_real_
    out$reason[masked & is.na(out$reason)] <-
      "the engine did not certify the prediction variance for this row"
  }
  out
}

# Conditional prediction variance for an LMM via the engine
# `predict_new_variance_with_level` FFI. The engine refits from the stored
# training frame and computes per-row se_fit / confidence / prediction bounds
# including the random-effect contribution. `se_data` is the frame to predict on
# (the training model frame for in-sample, or user `newdata`).
mm_lmm_prediction_variance <- function(fit, se_data, allow_new_levels, level) {
  policy <- if (isTRUE(allow_new_levels)) "population" else "error"
  spec_data <- mm_translate_data(mm_engine_frame(fit))
  new_data <- mm_translate_data(mm_engine_newdata(fit, se_data))
  control_json <- mm_refit_control_json(fit)
  json <- tryCatch(
    .Call(
      wrap__mm_lmm_predict_new_variance_json,
      mm_coerce_formula_string(mm_engine_formula(fit)),
      isTRUE(fit$REML),
      spec_data$column_order,
      spec_data$numeric_columns,
      spec_data$categorical_values,
      spec_data$categorical_levels,
      spec_data$categorical_ordered,
      mm_bridge_weights(fit$weights),
      as.character(control_json),
      new_data$column_order,
      new_data$numeric_columns,
      new_data$categorical_values,
      new_data$categorical_levels,
      # Ordered-factor coding for newdata is intentionally empty: the engine
      # re-derives newdata contrasts from the fitted training snapshot
      # (align_newdata_to_training), so flagging ordered columns here is a
      # no-op that could only spuriously error on an ordered newdata column
      # observed at <2 levels. Training flags (above) drive the coding.
      character(0),
      policy,
      as.numeric(level)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) mm_abort_from_bridge(json, newdata = se_data)
  mm_parse_prediction_variance(json, nrow(se_data))
}

# Engine family string as used at fit time (binomial resolves to "bernoulli"
# for 0/1 responses and "binomial" when trial weights are present), recomputed
# from the stored family + prepped weights so the variance refit reproduces the
# original fit exactly.
mm_glmm_engine_family <- function(fit) {
  fam <- fit$family$family
  if (identical(fam, "binomial")) {
    if (is.null(fit$weights)) "bernoulli" else "binomial"
  } else if (identical(fam, "negative_binomial")) {
    mm_glmm_nb_engine_spec(fit$family)
  } else {
    fam
  }
}

# Conditional prediction variance for a GLMM via the engine
# `predict_new_variance_with_level` FFI. `scale` is "link" or "response"; the
# engine propagates the variance through the link, so no R-side delta method is
# needed. Certified (joint_laplace) fits return available rows; fast-PIRLS /
# new-level rows come back with se_fit = NA and a reason.
mm_glmm_prediction_variance <- function(fit, se_data, scale, allow_new_levels, level) {
  policy <- if (isTRUE(allow_new_levels)) "population" else "error"
  spec_data <- mm_translate_data(mm_engine_frame(fit))
  new_data <- mm_translate_data(mm_engine_newdata(fit, se_data))
  control_json <- mm_refit_control_json(fit)
  json <- tryCatch(
    .Call(
      wrap__mm_glmm_predict_new_variance_json,
      mm_coerce_formula_string(mm_engine_formula(fit)),
      mm_glmm_engine_family(fit),
      fit$family$link,
      fit$method,
      as.integer(fit$nAGQ),
      spec_data$column_order,
      spec_data$numeric_columns,
      spec_data$categorical_values,
      spec_data$categorical_levels,
      spec_data$categorical_ordered,
      mm_bridge_weights(fit$weights),
      mm_bridge_weights(fit$offset),
      as.character(control_json),
      new_data$column_order,
      new_data$numeric_columns,
      new_data$categorical_values,
      new_data$categorical_levels,
      # Ordered-factor coding for newdata is intentionally empty: the engine
      # re-derives newdata contrasts from the fitted training snapshot
      # (align_newdata_to_training), so flagging ordered columns here is a
      # no-op that could only spuriously error on an ordered newdata column
      # observed at <2 levels. Training flags (above) drive the coding.
      character(0),
      scale,
      policy,
      as.numeric(level)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) mm_abort_from_bridge(json, newdata = se_data)
  mm_parse_prediction_variance(json, nrow(se_data))
}

# Fixed-effect-only prediction (re.form = NA path). Build the FE design
# matrix from the fit's stored fixed formula, align column names to beta,
# and multiply. Treating columns absent from `newdata` as zero matches the
# upstream `predict_new` behavior for partial X matrices.
mm_predict_fixed_only <- function(fit, newdata) {
  mm_new <- tryCatch(
    mm_engine_fixed_matrix(fit, newdata),
    error = function(cnd) cnd
  )
  if (inherits(mm_new, "condition")) {
    # model.frame failures (e.g. a missing predictor column) surface as a
    # data error; the engine-basis alignment failure already carries its own
    # typed class and is re-raised as-is.
    if (inherits(mm_new, "mm_inference_unavailable")) stop(mm_new)
    mm_abort(
      message = sprintf("Failed to build model frame for newdata: %s",
                        conditionMessage(mm_new)),
      class = "mm_data_error",
      parent = mm_new
    )
  }
  # Columns are aligned to names(beta) by mm_engine_fixed_matrix().
  pred <- as.numeric(mm_new %*% fit$beta)
  names(pred) <- rownames(newdata)
  pred
}

mm_training_xlevels <- function(fit) {
  rhs <- stats::delete.response(stats::terms(mm_fixed_formula(fit)))
  fe_vars <- all.vars(rhs)
  factor_cols <- vapply(fit$model_frame, is.factor, logical(1))
  factor_names <- intersect(names(factor_cols)[factor_cols], fe_vars)
  if (!length(factor_names)) return(list())
  lapply(fit$model_frame[factor_names], levels)
}

mm_training_contrasts <- function(fit) {
  # Restrict to factors that appear in the fixed-effect part; passing
  # contrasts for random-effect grouping factors makes model.matrix warn
  # ("variable ... is absent, its contrast will be ignored").
  fe_vars <- all.vars(stats::delete.response(stats::terms(mm_fixed_formula(fit))))
  is_fac <- vapply(fit$model_frame, is.factor, logical(1))
  factor_vars <- intersect(names(is_fac)[is_fac], fe_vars)
  if (!length(factor_vars)) return(NULL)
  # The engine codes unordered factors with treatment (dummy) contrasts and
  # ordered factors with contr.poly (orthonormal polynomial trends), matching
  # lme4/R defaults. We must force the SAME per-factor coding so the
  # reconstructed design matches `beta`'s basis; see mm_engine_fixed_matrix()
  # for the full rationale.
  codings <- vapply(
    factor_vars,
    function(v) if (is.ordered(fit$model_frame[[v]])) "contr.poly" else "contr.treatment",
    character(1)
  )
  stats::setNames(as.list(codings), factor_vars)
}

# Build the fixed-effect design matrix for arbitrary `data` (newdata or a
# reference grid) in the SAME basis the fit exposes, with columns named and
# ordered to match `fit$beta`.
#
# fit$beta carries lme4/R-style names in model.matrix() column order (see
# mm_apply_lme4_coef_naming()), so reconstructing with the fit-time xlevels
# and per-factor contrasts (treatment for unordered, contr.poly for ordered —
# see mm_training_contrasts()) yields columns that align by NAME. Positional
# alignment or blind renames silently corrupt X %*% beta whenever an ordered
# factor or an interaction is present, so any unmatched coefficient aborts.
mm_engine_fixed_matrix <- function(fit, data) {
  rhs <- stats::delete.response(stats::terms(mm_fixed_formula(fit)))
  # Stateful formula terms: rebuild the synthetic columns from the user's
  # transforms with the training basis (predvars), then name the columns as
  # model.matrix() names the user's formula.
  data <- mm_expand_newdata(fit, data)
  mf <- stats::model.frame(rhs, data = data, na.action = stats::na.pass,
                           xlev = mm_training_xlevels(fit))
  X <- stats::model.matrix(rhs, data = mf,
                           contrasts.arg = mm_training_contrasts(fit))
  if (!is.null(fit$expansion)) {
    colnames(X) <- mm_expansion_user_names(colnames(X), fit$expansion)
  }
  beta_names <- names(fit$beta)

  missing <- setdiff(beta_names, colnames(X))
  if (length(missing)) {
    mm_abort(
      message = paste0(
        "Could not reconstruct the fixed-effect design in the engine's ",
        "coefficient basis; missing column(s): ",
        paste(missing, collapse = ", "), "."
      ),
      class = "mm_inference_unavailable",
      expected = beta_names,
      observed = colnames(X)
    )
  }
  X[, beta_names, drop = FALSE]
}
