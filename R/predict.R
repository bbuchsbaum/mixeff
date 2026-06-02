#' Predict from a fitted mixeff LMM
#'
#' Predictions follow the lme4 generic shape. In-sample predictions reuse
#' the cached fitted/fixed values; new-data predictions are dispatched
#' through the Rust `predict_new` contract.
#'
#' @param object A fitted `mm_lmm` object.
#' @param newdata Optional new data. Must be a `data.frame` containing every
#'   variable referenced by the model's formula. Categorical levels must
#'   either match the training factor levels or trigger the `allow.new.levels`
#'   policy.
#' @param re.form Random-effects conditioning, following lme4's basic
#'   convention. `NULL` returns conditional predictions; `NA` (or `~0`)
#'   returns population-level (fixed-effect) predictions. Conditioning on a
#'   subset of grouping factors via a one-sided formula is not supported by
#'   the current Rust contract and raises `mm_inference_unavailable`.
#' @param allow.new.levels When `FALSE` (default), unseen grouping levels in
#'   `newdata` raise `mm_inference_unavailable` through the Rust
#'   `NewReLevels::Error` policy. When `TRUE`, unseen levels are replaced by
#'   the population mean (zero random effect), matching
#'   `lme4::predict(allow.new.levels = TRUE)`.
#' @param type Prediction scale. Gaussian LMMs use the same values for
#'   `"response"` and `"link"`.
#' @param se.fit Logical; when `TRUE`, returns a list with `fit` and `se.fit`.
#'   For population predictions (`re.form = NA`) the standard error is the
#'   Wald SE of the fixed-effect linear predictor, `sqrt(diag(X V X'))`. For
#'   conditional predictions (`re.form = NULL`) the SE is `NA` with a reason
#'   (random-effect posterior variances are not certified; lme4 likewise
#'   returns no conditional SE).
#' @param interval Interval type for population (`re.form = NA`) predictions:
#'   `"confidence"` for the mean (`fit +/- z*se`) or `"prediction"` for a new
#'   observation (adds the residual variance). Conditional intervals raise
#'   `mm_inference_unavailable`. Returns a matrix with `fit`/`lwr`/`upr`.
#' @param level Confidence level for `interval` / `se.fit` intervals.
#' @param scaled Logical; when `TRUE`, residuals are divided by the residual
#'   scale.
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
                           ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  mm_reject_unsupported_dots(
    list(...), "predict",
    c(random.only = "use `re.form = NA` for population-level (fixed-only) predictions or `re.form = NULL` for conditional predictions.")
  )
  if (!is.logical(allow.new.levels) || length(allow.new.levels) != 1L ||
      is.na(allow.new.levels)) {
    mm_abort(
      message = "`allow.new.levels` must be TRUE or FALSE.",
      class = "mm_arg_error",
      input = allow.new.levels
    )
  }
  target <- mm_prediction_target(re.form)
  if (identical(target, "unsupported")) {
    mm_abort(
      message = paste(
        "`re.form` values other than NULL, NA, or ~0 are not supported by the",
        "current Rust prediction contract."
      ),
      class = "mm_inference_unavailable",
      input = re.form
    )
  }

  if (is.null(newdata)) {
    pred <- switch(
      target,
      conditional = object$fitted,
      population  = object$fixed_fitted
    )
    names(pred) <- rownames(object$model_frame)
    se_data <- object$model_frame
  } else {
    pred <- mm_predict_newdata(object, newdata, target, allow.new.levels)
    se_data <- newdata
  }

  want_se <- isTRUE(se.fit) || !identical(interval, "none")
  if (!want_se) {
    return(pred)
  }

  # Standard errors / intervals are defined here only for the population
  # (fixed-effect-only) target: se = sqrt(diag(X V X')). Conditional SEs need
  # random-effect posterior variances, which the contract does not certify
  # (lme4's predict() likewise returns no conditional SE) -> NA / refuse.
  if (!identical(target, "population")) {
    if (!identical(interval, "none")) {
      mm_abort(
        message = paste(
          "Conditional (`re.form = NULL`) prediction intervals require",
          "random-effect posterior variances, which are not certified. Use",
          "`re.form = NA` for population-level intervals."
        ),
        class = "mm_inference_unavailable",
        input = interval
      )
    }
    se <- rep(NA_real_, length(pred))
    attr(se, "mm_unavailable_reason") <- "conditional_prediction_se_unavailable"
    out <- list(fit = pred, se.fit = se)
    attr(out, "mm_unavailable_reason") <- "conditional_prediction_se_unavailable"
    return(out)
  }

  se <- mm_fixed_prediction_se(object, se_data)
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
  out
}

#' @rdname predict.mm_lmm
#' @export
residuals.mm_lmm <- function(object,
                             type = c("response", "pearson", "deviance",
                                      "working"),
                             scaled = FALSE, ...) {
  type <- match.arg(type)
  out <- object$residuals
  # For a Gaussian LMM with iid residuals, the working and deviance residuals
  # equal the response residuals; the Pearson residual divides by the residual
  # scale (unit working weights).
  if (identical(type, "pearson")) {
    out <- out / object$sigma
  }
  if (isTRUE(scaled)) {
    out <- out / object$sigma
  }
  names(out) <- rownames(object$model_frame)
  out
}

#' @rdname predict.mm_lmm
#' @export
fitted.mm_glmm <- fitted.mm_lmm

#' @rdname predict.mm_lmm
#' @export
residuals.mm_glmm <- function(object, type = c("response"), ...) {
  # The engine returns response-scale residuals for GLMMs; Pearson/deviance
  # residuals are not certified by the current contract, so refuse rather than
  # silently return response residuals under a different label.
  type <- match.arg(type)
  out <- object$residuals
  names(out) <- rownames(object$model_frame)
  out
}

#' Predict from a fitted mixeff GLMM
#'
#' GLMM predictions are computed on the R side from the stored fixed effects
#' (population, `re.form = NA`) or fixed effects plus conditional modes
#' (`re.form = NULL`), then mapped through the family link. This mirrors
#' `lme4::predict.merMod` for generalized models: `type = "link"` returns the
#' linear predictor and `type = "response"` the mean.
#'
#' In-sample response predictions reuse the engine's certified fitted values;
#' all other paths are an explicit R-side plug-in (not engine-certified), so
#' `se.fit` returns `NA` with a reason and `interval` is refused. An
#' engine-certified `GeneralizedLinearMixedModel::predict_new` (symmetry with
#' the LMM `predict_new` path) is tracked upstream in mixeff-rs as
#' bd-01KT3XJD6SE54G7GT3Z2T6HMEX; when it lands the R-side reconstruction below
#' can be replaced by an FFI call.
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
                            ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  mm_reject_unsupported_dots(
    list(...), "predict",
    c(random.only = "use `re.form = NA` for population-level (fixed-only) predictions or `re.form = NULL` for conditional predictions.")
  )
  if (!is.logical(allow.new.levels) || length(allow.new.levels) != 1L ||
      is.na(allow.new.levels)) {
    mm_abort(
      message = "`allow.new.levels` must be TRUE or FALSE.",
      class = "mm_arg_error",
      input = allow.new.levels
    )
  }
  if (!identical(interval, "none")) {
    mm_abort(
      message = paste(
        "`interval` prediction requires prediction standard errors, which are",
        "not certified for GLMMs by the current Rust inference contract."
      ),
      class = "mm_inference_unavailable",
      input = interval
    )
  }
  target <- mm_prediction_target(re.form)
  if (identical(target, "unsupported")) {
    mm_abort(
      message = paste(
        "`re.form` values other than NULL, NA, or ~0 are not supported by the",
        "current GLMM prediction path."
      ),
      class = "mm_inference_unavailable",
      input = re.form
    )
  }
  if (!is.null(object$offset) && !is.null(newdata)) {
    mm_abort(
      message = paste(
        "Prediction on `newdata` is not supported for a GLMM fitted with an",
        "`offset`; the offset cannot be reconstructed for new rows."
      ),
      class = "mm_inference_unavailable",
      input = "offset"
    )
  }

  family <- mm_glmm_family_from_info(object$family)

  if (is.null(newdata)) {
    if (identical(target, "conditional")) {
      # Engine-certified conditional mean; eta via the link.
      mu <- object$fitted
      eta <- family$linkfun(mu)
      nm <- rownames(object$model_frame)
    } else {
      eta <- as.numeric(mm_predict_fixed_only(object, object$model_frame))
      if (!is.null(object$offset)) eta <- eta + as.numeric(object$offset)
      mu <- family$linkinv(eta)
      nm <- rownames(object$model_frame)
    }
  } else {
    needed <- all.vars(stats::delete.response(stats::terms(object$formula)))
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
    eta <- as.numeric(mm_predict_fixed_only(object, newdata))
    if (identical(target, "conditional")) {
      eta <- eta + mm_glmm_re_eta(object, newdata, allow.new.levels)
    }
    mu <- family$linkinv(eta)
    nm <- rownames(newdata)
  }

  pred <- if (identical(type, "response")) as.numeric(mu) else as.numeric(eta)
  names(pred) <- nm

  if (isTRUE(se.fit)) {
    # Population (fixed-only) Wald SE of the link-scale linear predictor; for
    # type="response" map through the link with the delta method. Conditional
    # SEs need RE posterior variances (not certified) -> NA with a reason.
    if (identical(target, "population")) {
      se_data <- if (is.null(newdata)) object$model_frame else newdata
      se_link <- mm_fixed_prediction_se(object, se_data)
      se <- if (identical(type, "response")) {
        se_link * abs(family$mu.eta(eta))
      } else {
        se_link
      }
    } else {
      se <- rep(NA_real_, length(pred))
      attr(se, "mm_unavailable_reason") <- "conditional_prediction_se_unavailable"
    }
    names(se) <- nm
    out <- list(fit = pred, se.fit = se)
    return(out)
  }
  pred
}

# Random-effect contribution to the GLMM linear predictor for arbitrary `data`,
# reconstructed from the stored conditional modes (BLUPs). For each random
# term, look up each row's grouping level in the BLUP table and accumulate
# basis_value * mode over the term's basis columns. Unseen levels contribute
# zero when allow_new_levels = TRUE, otherwise raise a typed condition (matching
# lme4's allow.new.levels semantics).
mm_glmm_re_eta <- function(fit, data, allow_new_levels) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  eta <- numeric(nrow(data))
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    group_label <- mm_random_term_group_label(fit, term, i)
    group <- as.character(mm_group_factor(data, group_label))
    blup <- fit$random_effects[[group_label]]
    if (is.null(blup)) {
      mm_abort(
        message = sprintf("No conditional modes stored for grouping factor `%s`.",
                          group_label),
        class = "mm_inference_unavailable",
        input = group_label
      )
    }
    idx <- match(group, rownames(blup))
    unseen <- is.na(idx)
    if (any(unseen) && !isTRUE(allow_new_levels)) {
      mm_abort(
        message = sprintf(
          paste0("`newdata` has %d row(s) with grouping level(s) of `%s` not ",
                 "seen during fitting. Set allow.new.levels = TRUE to predict ",
                 "them at the population mean (zero random effect)."),
          sum(unseen), group_label
        ),
        class = "mm_inference_unavailable",
        group = group_label
      )
    }
    basis <- term$basis %||% list()
    basis_labels <- if (length(basis)) {
      vapply(basis, mm_basis_label, character(1))
    } else {
      "(Intercept)"
    }
    basis_values <- if (length(basis)) {
      lapply(basis, mm_basis_values, frame = data)
    } else {
      list(rep(1, nrow(data)))
    }
    for (j in seq_along(basis_values)) {
      lab <- basis_labels[[j]]
      if (!lab %in% names(blup)) next
      modes <- blup[[lab]][idx]
      modes[is.na(modes)] <- 0  # unseen levels -> zero RE (allow.new.levels)
      eta <- eta + basis_values[[j]] * modes
    }
  }
  eta
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

# Dispatch to the upstream predict_new contract (conditional path) or
# compute fixed-only predictions via the cached fixed-effect design matrix
# (population path). lme4 mirrors this split: re.form = NA bypasses the
# RE bridge entirely, so we do the same on the R side without paying the
# refit cost in the Rust FFI.
mm_predict_newdata <- function(fit, newdata, target, allow_new_levels) {
  if (!is.data.frame(newdata)) {
    mm_abort(
      message = "`newdata` must be a data.frame.",
      class = "mm_data_error",
      input = newdata
    )
  }
  needed <- all.vars(stats::delete.response(stats::terms(fit$formula)))
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

  if (identical(target, "population")) {
    return(mm_predict_fixed_only(fit, newdata))
  }

  mm_predict_conditional_newdata(fit, newdata, allow_new_levels)
}

# Conditional newdata prediction via the Rust `predict_new` FFI. The Rust
# contract refits the model from formula + training data, so the call is
# self-contained.
mm_predict_conditional_newdata <- function(fit, newdata, allow_new_levels) {
  policy <- if (isTRUE(allow_new_levels)) "population" else "error"

  spec_data <- mm_translate_data(fit$model_frame)
  formula_string <- mm_coerce_formula_string(fit$formula)
  control_json <- jsonlite::toJSON(unclass(fit$control %||% mm_control()),
                                   auto_unbox = TRUE, null = "null")
  new_data <- mm_translate_data(newdata)

  json <- tryCatch(
    .Call(
      wrap__mm_lmm_predict_new_json,
      formula_string,
      isTRUE(fit$REML),
      spec_data$column_order,
      spec_data$numeric_columns,
      spec_data$categorical_values,
      spec_data$categorical_levels,
      mm_bridge_weights(fit$weights),
      as.character(control_json),
      new_data$column_order,
      new_data$numeric_columns,
      new_data$categorical_values,
      new_data$categorical_levels,
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
  # The Rust engine codes EVERY factor with treatment (dummy) contrasts,
  # including ordered factors (which R would otherwise give contr.poly). We
  # must force contr.treatment so the reconstructed design matches `beta`'s
  # basis; see mm_engine_fixed_matrix() for the full rationale.
  stats::setNames(rep(list("contr.treatment"), length(factor_vars)), factor_vars)
}

# Build the fixed-effect design matrix for arbitrary `data` (newdata or a
# reference grid) in the SAME basis the engine used at fit time, with columns
# named and ordered to match `fit$beta`.
#
# Why this exists: the engine codes all factors with treatment contrasts and
# labels coefficients in a mixeff-specific encoding (e.g. "recipe: B",
# "recipe: B:temperature: 215"). R's model.matrix() instead (a) uses
# contr.poly for ordered factors and (b) orders interaction columns with the
# first factor varying fastest, whereas the engine varies the last factor
# fastest. Either mismatch silently corrupts X %*% beta. The previous code
# relied on positional alignment + a blind colnames<- rename, which produced
# wrong predictions/marginal means whenever an ordered factor or an
# interaction was present (e.g. lme4::cake). We instead force treatment
# contrasts, translate R-style column names into the engine encoding, and
# align by name.
mm_engine_fixed_matrix <- function(fit, data) {
  rhs <- stats::delete.response(stats::terms(mm_fixed_formula(fit)))
  mf <- stats::model.frame(rhs, data = data, na.action = stats::na.pass,
                           xlev = mm_training_xlevels(fit))
  X <- stats::model.matrix(rhs, data = mf,
                           contrasts.arg = mm_training_contrasts(fit))
  beta_names <- names(fit$beta)

  # Factor variables present in the fixed part, longest name first so that a
  # factor whose name is a prefix of another (e.g. "rec" vs "recipe") is
  # matched greedily.
  fe_vars <- all.vars(rhs)
  is_fac <- vapply(fit$model_frame, is.factor, logical(1))
  factor_vars <- intersect(names(is_fac)[is_fac], fe_vars)
  factor_vars <- factor_vars[order(nchar(factor_vars), decreasing = TRUE)]

  translate_component <- function(comp) {
    for (v in factor_vars) {
      if (startsWith(comp, v)) {
        lev <- substring(comp, nchar(v) + 1L)
        if (nzchar(lev)) return(paste0(v, ": ", lev))
      }
    }
    comp
  }
  translate_name <- function(nm) {
    if (nm == "(Intercept)") return(nm)
    parts <- strsplit(nm, ":", fixed = TRUE)[[1L]]
    paste(vapply(parts, translate_component, character(1)), collapse = ":")
  }
  colnames(X) <- vapply(colnames(X), translate_name, character(1))

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
