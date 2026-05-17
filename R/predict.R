#' Predict from a fitted mixeff LMM
#'
#' Phase 2 supports in-sample predictions for [lmm()] fits. `re.form = NULL`
#' returns conditional fitted values, while `re.form = NA` returns the
#' population-level fixed-effect part. New data and prediction standard errors
#' require the later lazy-handle prediction layer.
#'
#' @param object A fitted `mm_lmm` object.
#' @param newdata Optional new data. Not supported: predictions are
#'   returned for the fitted data only.
#' @param re.form Random-effects conditioning, following lme4's basic
#'   convention: `NULL` for conditional predictions and `NA` for population
#'   predictions. Other values are not supported.
#' @param allow.new.levels Reserved for future new-data prediction support.
#' @param type Prediction scale. Gaussian LMMs use the same values for
#'   `"response"` and `"link"`.
#' @param se.fit Logical; when `TRUE`, returns `NA` standard errors with an
#'   unavailable-reason attribute.
#' @param interval Prediction interval type. Intervals require prediction
#'   standard errors and therefore raise `mm_inference_unavailable` until Rust
#'   certifies them.
#' @param level Confidence level for future interval support.
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
  if (!is.null(newdata)) {
    mm_abort(
      message = paste(
        "`newdata` prediction is not supported by the current Rust",
        "inference contract; predictions are returned for the fitted data."
      ),
      class = "mm_inference_unavailable",
      input = newdata
    )
  }
  if (!is.logical(allow.new.levels) || length(allow.new.levels) != 1L ||
      is.na(allow.new.levels)) {
    mm_abort(
      message = "`allow.new.levels` must be TRUE or FALSE.",
      class = "mm_arg_error",
      input = allow.new.levels
    )
  }

  target <- mm_prediction_target(re.form)
  pred <- switch(
    target,
    conditional = object$fitted,
    population = object$fixed_fitted,
    mm_abort(
      message = "`re.form` values other than NULL or NA are not supported.",
      class = "mm_inference_unavailable",
      input = re.form
    )
  )
  names(pred) <- rownames(object$model_frame)

  if (!identical(interval, "none")) {
    mm_abort(
      message = paste(
        "`interval` prediction requires prediction standard errors, which are",
        "not certified by the current Rust inference contract."
      ),
      class = "mm_inference_unavailable",
      input = interval
    )
  }

  if (isTRUE(se.fit)) {
    se <- rep(NA_real_, length(pred))
    attr(se, "mm_unavailable_reason") <- "prediction_se_unavailable_phase_2"
    out <- list(fit = pred, se.fit = se)
    attr(out, "mm_unavailable_reason") <- "prediction_se_unavailable_phase_2"
    return(out)
  }
  pred
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
residuals.mm_lmm <- function(object, type = c("response"), ...) {
  type <- match.arg(type)
  out <- object$residuals
  names(out) <- rownames(object$model_frame)
  out
}

mm_prediction_target <- function(re.form) {
  if (is.null(re.form)) {
    return("conditional")
  }
  if (length(re.form) == 1L && is.na(re.form)) {
    return("population")
  }
  "unsupported"
}
