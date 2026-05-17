#' Extract components from a fitted mixeff LMM
#'
#' These methods provide the common lme4-style extractor surface for
#' [lmm()] fits. The required values are stored directly on the R object or
#' rebuilt lazily from the serialized artifact, so these methods do not require
#' a live Rust handle after `saveRDS()` / `readRDS()`.
#'
#' @param object,x,formula A fitted `mm_lmm` object.
#' @param condVar Logical; when `TRUE`, Phase 2 returns the random-effects
#'   tables with an `NA` `postVar` array and an `mm_unavailable_reason`
#'   attribute rather than fabricating conditional variances.
#' @param k Penalty per parameter for [AIC()].
#' @param ... Reserved for generic compatibility.
#'
#' @name mm_lmm-methods
NULL

#' @rdname mm_lmm-methods
#' @export
fixef <- function(object, ...) {
  UseMethod("fixef")
}

#' @rdname mm_lmm-methods
#' @export
fixef.mm_lmm <- function(object, ...) {
  object$beta
}

#' @rdname mm_lmm-methods
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' @rdname mm_lmm-methods
#' @export
ranef.mm_lmm <- function(object, condVar = FALSE, ...) {
  if (isTRUE(condVar)) {
    out <- lapply(object$random_effects, mm_attach_ranef_postvar)
    class(out) <- class(object$random_effects)
    attr(out, "mm_unavailable_reason") <- "random_effect_conditional_variance_unavailable"
    return(out)
  }
  object$random_effects
}

mm_attach_ranef_postvar <- function(df) {
  p <- ncol(df)
  n <- nrow(df)
  postvar <- array(
    NA_real_,
    dim = c(p, p, n),
    dimnames = list(names(df), names(df), rownames(df))
  )
  attr(df, "postVar") <- postvar
  attr(df, "mm_unavailable_reason") <- "random_effect_conditional_variance_unavailable"
  df
}

#' @rdname mm_lmm-methods
#' @export
coef.mm_lmm <- function(object, ...) {
  re <- ranef(object)
  fixed <- fixef(object)
  out <- lapply(re, function(df) {
    df_out <- df
    for (nm in names(df_out)) {
      if (nm %in% names(fixed)) {
        df_out[[nm]] <- df_out[[nm]] + fixed[[nm]]
      }
    }
    df_out
  })
  class(out) <- c("mm_coef", "list")
  out
}

#' @rdname mm_lmm-methods
#' @export
VarCorr <- function(x, ...) {
  UseMethod("VarCorr")
}

#' @rdname mm_lmm-methods
#' @export
VarCorr.mm_lmm <- function(x, ...) {
  x$varcorr
}

#' @rdname mm_lmm-methods
#' @export
sigma.mm_lmm <- function(object, ...) {
  object$sigma
}

#' @rdname mm_lmm-methods
#' @export
logLik.mm_lmm <- function(object, ...) {
  structure(
    object$logLik,
    df = object$dof,
    nobs = object$nobs,
    class = "logLik"
  )
}

#' @rdname mm_lmm-methods
#' @export
deviance.mm_lmm <- function(object, ...) {
  object$deviance
}

#' @rdname mm_lmm-methods
#' @export
AIC.mm_lmm <- function(object, ..., k = 2) {
  dots <- list(...)
  if (length(dots)) {
    mm_abort(
      message = paste(
        "AIC comparison across multiple objects is not supported;",
        "call `AIC()` per model and compare with `compare()`."
      ),
      class = "mm_inference_unavailable",
      input = dots
    )
  }
  -2 * object$logLik + k * object$dof
}

#' @rdname mm_lmm-methods
#' @export
BIC.mm_lmm <- function(object, ...) {
  dots <- list(...)
  if (length(dots)) {
    mm_abort(
      message = paste(
        "BIC comparison across multiple objects is not supported;",
        "call `BIC()` per model and compare with `compare()`."
      ),
      class = "mm_inference_unavailable",
      input = dots
    )
  }
  object$BIC
}

#' @rdname mm_lmm-methods
#' @export
nobs.mm_lmm <- function(object, ...) {
  object$nobs
}

#' @rdname mm_lmm-methods
#' @export
df.residual.mm_lmm <- function(object, ...) {
  object$df_residual
}

#' @rdname mm_lmm-methods
#' @export
formula.mm_lmm <- function(x, ...) {
  x$formula
}

#' @rdname mm_lmm-methods
#' @export
model.frame.mm_lmm <- function(formula, ...) {
  formula$model_frame
}

mm_ranef_from_terms <- function(terms) {
  out <- list()
  for (term in terms %||% list()) {
    group <- as.character(term$group)
    names <- as.character(unlist(term$names, use.names = FALSE))
    levels <- as.character(unlist(term$levels, use.names = FALSE))
    values <- term$values %||% list()
    mat <- if (length(values)) {
      do.call(rbind, lapply(values, function(row) {
        as.numeric(unlist(row, use.names = FALSE))
      }))
    } else {
      matrix(numeric(), nrow = 0L, ncol = length(names))
    }
    if (is.null(dim(mat))) {
      mat <- matrix(mat, ncol = length(names), byrow = TRUE)
    }
    df <- as.data.frame(mat, check.names = FALSE)
    names(df) <- names
    rownames(df) <- levels

    if (group %in% names(out)) {
      existing <- out[[group]]
      if (identical(sort(rownames(existing)), sort(rownames(df)))) {
        df <- df[rownames(existing), , drop = FALSE]
      }
      out[[group]] <- cbind(existing, df)
    } else {
      out[[group]] <- df
    }
  }
  class(out) <- c("mm_ranef", "list")
  out
}

mm_varcorr_from_result <- function(varcorr) {
  components <- lapply(varcorr$components %||% list(), function(component) {
    names <- as.character(unlist(component$names, use.names = FALSE))
    std_dev <- as.numeric(unlist(component$std_dev, use.names = FALSE))
    correlations <- as.numeric(unlist(component$correlations, use.names = FALSE))
    rows <- vector("list", length(names))
    for (i in seq_along(names)) {
      rows[[i]] <- data.frame(
        group = as.character(component$group),
        name = names[[i]],
        variance = std_dev[[i]]^2,
        std_dev = std_dev[[i]],
        correlation = mm_varcorr_correlation_text(correlations, i),
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  })
  table <- if (length(components)) {
    out <- do.call(rbind, components)
    rownames(out) <- NULL
    out
  } else {
    data.frame(
      group = character(),
      name = character(),
      variance = numeric(),
      std_dev = numeric(),
      correlation = character(),
      stringsAsFactors = FALSE
    )
  }
  residual_sd <- as.numeric(varcorr$residual_sd %||% NA_real_)
  table$boundary <- mm_varcorr_boundary_flag(table$std_dev, residual_sd)

  out <- list(
    table = table,
    residual_sd = residual_sd
  )
  class(out) <- c("mm_varcorr", "list")
  out
}

# PRD §9.4: boundary fits are reported model state, not a warning.
# Flag a variance component as on-boundary when its std_dev is
# negligibly small relative to the residual scale. Threshold:
# std_dev <= max(1e-8, 1e-6 * residual_sd). The constant floor protects
# the no-residual case (REML on a saturated design); the relative
# component protects models with very small or very large residuals.
mm_varcorr_boundary_flag <- function(std_dev, residual_sd) {
  if (!length(std_dev)) return(logical())
  scale <- if (is.finite(residual_sd) && residual_sd > 0) residual_sd else 1
  threshold <- max(1e-8, 1e-6 * scale)
  is.finite(std_dev) & std_dev <= threshold
}

mm_varcorr_correlation_text <- function(correlations, row_index) {
  if (row_index <= 1L || !length(correlations)) {
    return("")
  }
  offset <- (row_index - 1L) * (row_index - 2L) / 2L
  vals <- correlations[seq.int(offset + 1L, offset + row_index - 1L)]
  paste(sprintf("%+.2f", vals), collapse = " ")
}
