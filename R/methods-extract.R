#' Extract components from a fitted mixeff LMM
#'
#' These methods provide the common lme4-style extractor surface for
#' [lmm()] fits. The required values are stored directly on the R object or
#' rebuilt lazily from the serialized artifact, so these methods do not require
#' a live Rust handle after `saveRDS()` / `readRDS()`.
#'
#' @param object,x,formula,fit A fitted `mm_lmm` or `mm_glmm` object.
#' @param REML Ignored; included for S3 compatibility with likelihood and
#'   deviance generics.
#' @param scale Ignored; included for S3 compatibility with [extractAIC()].
#' @param correlation Logical; accepted for S3 compatibility with [vcov()].
#' @param condVar Logical; when `TRUE`, Phase 2 returns the random-effects
#'   tables with an `NA` `postVar` array and an `mm_unavailable_reason`
#'   attribute rather than fabricating conditional variances.
#' @param k Penalty per parameter for [AIC()].
#' @param ... Reserved for generic compatibility.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   y = rnorm(60), x = rnorm(60),
#'   g = factor(rep(seq_len(10), each = 6))
#' )
#' fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
#' fixef(fit)
#' VarCorr(fit)
#' head(ranef(fit)$g)
#' sigma(fit)
#' logLik(fit)
#' nobs(fit)
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
fixef.mm_glmm <- fixef.mm_lmm

#' @rdname mm_lmm-methods
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' @rdname mm_lmm-methods
#' @export
ranef.mm_lmm <- function(object, condVar = FALSE, ...) {
  if (isTRUE(condVar)) {
    postvars <- tryCatch(
      mm_cond_var_postvars(object),
      error = function(cnd) cnd
    )
    if (inherits(postvars, "condition")) {
      out <- lapply(object$random_effects, mm_attach_ranef_postvar_unavailable,
                    reason = conditionMessage(postvars))
      class(out) <- class(object$random_effects)
      attr(out, "mm_unavailable_reason") <-
        "random_effect_conditional_variance_unavailable"
      return(out)
    }
    out <- mm_attach_ranef_postvars(object$random_effects, postvars)
    class(out) <- class(object$random_effects)
    return(out)
  }
  object$random_effects
}

#' @rdname mm_lmm-methods
#' @export
ranef.mm_glmm <- function(object, condVar = FALSE, ...) {
  if (isTRUE(condVar)) {
    out <- lapply(
      object$random_effects,
      mm_attach_ranef_postvar_unavailable,
      reason = "random_effect_conditional_variance_unavailable_for_glmm"
    )
    class(out) <- class(object$random_effects)
    attr(out, "mm_unavailable_reason") <-
      "random_effect_conditional_variance_unavailable_for_glmm"
    return(out)
  }
  object$random_effects
}

# Fallback postvar attachment used when the Rust cond_var() bridge fails
# (typed-refusal path). We keep the array-shaped attribute so downstream
# callers that test attr(df, "postVar") never see NULL, but every entry is
# NA and an explicit `mm_unavailable_reason` attribute records why.
mm_attach_ranef_postvar_unavailable <- function(df, reason) {
  p <- ncol(df)
  n <- nrow(df)
  postvar <- array(
    NA_real_,
    dim = c(p, p, n),
    dimnames = list(names(df), names(df), rownames(df))
  )
  attr(df, "postVar") <- postvar
  attr(df, "mm_unavailable_reason") <-
    "random_effect_conditional_variance_unavailable"
  attr(df, "mm_cond_var_error") <- reason
  df
}

# Back-compat alias; some callers (e.g. revive() helpers) still reference
# the old name. New code should use mm_attach_ranef_postvar_unavailable.
mm_attach_ranef_postvar <- function(df) {
  mm_attach_ranef_postvar_unavailable(
    df,
    reason = "random_effect_conditional_variance_unavailable"
  )
}

mm_attach_ranef_postvars <- function(ranef_list, postvars) {
  for (group in names(ranef_list)) {
    df <- ranef_list[[group]]
    pv <- postvars[[group]]
    if (is.null(pv)) {
      ranef_list[[group]] <- mm_attach_ranef_postvar_unavailable(
        df,
        reason = sprintf("cond_var bridge returned no payload for group `%s`",
                         group)
      )
      next
    }
    df_cols <- names(df)
    df_levels <- rownames(df)
    if (!setequal(df_cols, dimnames(pv)[[1L]]) ||
        !setequal(df_levels, dimnames(pv)[[3L]])) {
      ranef_list[[group]] <- mm_attach_ranef_postvar_unavailable(
        df,
        reason = sprintf(
          "cond_var bridge returned mismatched names/levels for group `%s`",
          group
        )
      )
      next
    }
    pv_aligned <- pv[df_cols, df_cols, df_levels, drop = FALSE]
    attr(df, "postVar") <- pv_aligned
    ranef_list[[group]] <- df
  }
  ranef_list
}

# Round-trip a cond_var FFI call for `fit` and return a named list keyed by
# grouping factor; each element is a `p × p × n` PSD array (dimnames =
# slope-names × slope-names × levels). Cached on fit$lazy_cache so repeated
# ranef(condVar=TRUE) calls do not pay the ~refit cost again.
mm_cond_var_postvars <- function(fit) {
  .mm_lazy(fit, "cond_var", mm_compute_cond_var_postvars)
}

mm_compute_cond_var_postvars <- function(fit) {
  spec_data <- mm_translate_data(fit$model_frame)
  formula_string <- mm_coerce_formula_string(fit$formula)
  control_json <- jsonlite::toJSON(unclass(fit$control %||% mm_control()),
                                   auto_unbox = TRUE, null = "null")

  json <- tryCatch(
    .Call(
      wrap__mm_lmm_cond_var_json,
      formula_string,
      isTRUE(fit$REML),
      spec_data$column_order,
      spec_data$numeric_columns,
      spec_data$categorical_values,
      spec_data$categorical_levels,
      spec_data$categorical_ordered,
      mm_bridge_weights(fit$weights),
      as.character(control_json)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json)
  }
  payload <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse cond_var JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        parent = cnd
      )
    }
  )
  schema <- payload$schema %||% list()
  if (!identical(as.character(schema$schema_name), "mixeff.lmm_cond_var") ||
      !identical(as.character(schema$schema_version), "1")) {
    mm_abort(
      message = "cond_var payload has an unknown schema header.",
      class = "mm_schema_error",
      input = payload
    )
  }

  out <- list()
  for (term in payload$terms %||% list()) {
    group <- as.character(term$group)
    slope_names <- as.character(unlist(term$names, use.names = FALSE))
    level_names <- as.character(unlist(term$levels, use.names = FALSE))
    flat <- as.numeric(unlist(term$postvar, use.names = FALSE))
    dims <- as.integer(unlist(term$dim, use.names = FALSE))
    arr <- array(flat, dim = dims,
                 dimnames = list(slope_names, slope_names, level_names))
    if (is.null(out[[group]])) {
      out[[group]] <- arr
    } else {
      # `(1|g) + (0+t|g)` collapses to one ranef(fit)$g frame whose columns
      # are the union of the per-term slopes. Mirror that here as a block-
      # diagonal postVar (the two RE blocks are independent in this
      # parameterisation, so off-diagonal blocks are zero).
      out[[group]] <- mm_merge_block_diag_postvar(out[[group]], arr, group)
    }
  }
  out
}

mm_merge_block_diag_postvar <- function(existing, incoming, group) {
  existing_levels <- dimnames(existing)[[3L]]
  incoming_levels <- dimnames(incoming)[[3L]]
  if (!setequal(existing_levels, incoming_levels)) {
    mm_abort(
      message = sprintf(
        "cond_var: cannot combine RE terms for group `%s` because their level sets differ.",
        group
      ),
      class = "mm_schema_error"
    )
  }
  incoming_aligned <- incoming[, , existing_levels, drop = FALSE]
  existing_names <- dimnames(existing)[[1L]]
  incoming_names <- dimnames(incoming_aligned)[[1L]]
  new_names <- c(existing_names, incoming_names)
  new_p <- length(new_names)
  n <- length(existing_levels)
  combined <- array(0, dim = c(new_p, new_p, n),
                    dimnames = list(new_names, new_names, existing_levels))
  old_p <- length(existing_names)
  combined[seq_len(old_p), seq_len(old_p), ] <- existing
  combined[(old_p + 1L):new_p, (old_p + 1L):new_p, ] <- incoming_aligned
  combined
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
coef.mm_glmm <- coef.mm_lmm

#' @rdname mm_lmm-methods
#' @export
VarCorr <- function(x, ...) {
  UseMethod("VarCorr")
}

#' @rdname mm_lmm-methods
#' @export
VarCorr.mm_lmm <- function(x, ...) {
  mm_varcorr_with_design_notes(x)
}

#' @rdname mm_lmm-methods
#' @export
VarCorr.mm_glmm <- VarCorr.mm_lmm

mm_varcorr_with_design_notes <- function(fit) {
  out <- fit$varcorr
  groups <- mm_design_weak_identifiability_groups(fit)
  attr(out, "mm_design_weak_identifiability_groups") <- groups
  if (!length(groups) || !is.data.frame(out$table) || !nrow(out$table)) {
    return(out)
  }
  if (!"note" %in% names(out$table)) {
    out$table$note <- ""
  }
  hit <- out$table$group %in% groups & out$table$name == "(Intercept)"
  out$table$note[hit] <- mm_note_append(
    out$table$note[hit],
    "[design_weak_identifiability]"
  )
  out
}

mm_note_append <- function(existing, note) {
  existing <- as.character(existing)
  if (!length(existing)) return(character())
  existing[is.na(existing)] <- ""
  note <- as.character(note)
  if (length(note) == 1L) {
    note <- rep(note, length(existing))
  }
  if (length(note) != length(existing)) {
    note <- rep(note, length.out = length(existing))
  }
  note[is.na(note)] <- ""
  for (i in seq_along(existing)) {
    if (!nzchar(note[[i]])) next
    if (!nzchar(existing[[i]])) {
      existing[[i]] <- note[[i]]
    } else if (!grepl(note[[i]], existing[[i]], fixed = TRUE)) {
      existing[[i]] <- paste(existing[[i]], note[[i]])
    }
  }
  existing
}

#' @rdname mm_lmm-methods
#' @export
sigma.mm_lmm <- function(object, ...) {
  object$sigma
}

#' @rdname mm_lmm-methods
#' @export
sigma.mm_glmm <- sigma.mm_lmm

# Refuse recognized-but-unsupported lme4 arguments that `...` would otherwise
# swallow silently, returning a plausible-but-wrong value (no silent surgery;
# PRD §8.1). `reject` maps an argument name to an actionable explanation.
mm_reject_unsupported_dots <- function(dots, method, reject) {
  if (!length(dots) || !length(reject)) return(invisible(NULL))
  nms <- names(dots)
  hit <- intersect(nms[nzchar(nms)], names(reject))
  if (length(hit)) {
    arg <- hit[[1L]]
    mm_abort(
      message = sprintf("`%s()` does not support the `%s` argument: %s",
                        method, arg, reject[[arg]]),
      class = "mm_arg_error",
      argument = arg,
      method = method
    )
  }
  invisible(NULL)
}

# The fit carries a single (RE)ML criterion. logLik()/deviance() cannot switch
# it after the fact, so honor a matching `REML=` and refuse a mismatching one
# (rather than silently returning the wrong-criterion value).
mm_check_reml_request <- function(object, REML, method) {
  if (is.null(REML)) return(invisible(NULL))
  if (!identical(isTRUE(REML), isTRUE(object$REML))) {
    mm_abort(
      message = sprintf(
        paste0("This model was fitted with REML = %s; `%s(REML = %s)` cannot ",
               "change the estimation criterion after the fact. Refit with ",
               "lmm(..., REML = %s) to obtain the %s value."),
        isTRUE(object$REML), method, isTRUE(REML), isTRUE(REML),
        if (isTRUE(REML)) "REML" else "ML"
      ),
      class = "mm_inference_unavailable",
      requested_REML = isTRUE(REML),
      fit_REML = isTRUE(object$REML)
    )
  }
  invisible(NULL)
}

#' @rdname mm_lmm-methods
#' @export
logLik.mm_lmm <- function(object, REML = NULL, ...) {
  mm_check_reml_request(object, REML, "logLik")
  structure(
    object$logLik,
    df = object$dof,
    nobs = object$nobs,
    class = "logLik"
  )
}

#' @rdname mm_lmm-methods
#' @export
logLik.mm_glmm <- logLik.mm_lmm

#' @rdname mm_lmm-methods
#' @export
deviance.mm_lmm <- function(object, REML = NULL, ...) {
  mm_check_reml_request(object, REML, "deviance")
  object$deviance
}

#' @rdname mm_lmm-methods
#' @export
deviance.mm_glmm <- deviance.mm_lmm

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
AIC.mm_glmm <- AIC.mm_lmm

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
BIC.mm_glmm <- BIC.mm_lmm

#' @rdname mm_lmm-methods
#' @export
nobs.mm_lmm <- function(object, ...) {
  object$nobs
}

#' @rdname mm_lmm-methods
#' @export
nobs.mm_glmm <- nobs.mm_lmm

#' @rdname mm_lmm-methods
#' @export
df.residual.mm_lmm <- function(object, ...) {
  object$df_residual
}

#' @rdname mm_lmm-methods
#' @export
df.residual.mm_glmm <- df.residual.mm_lmm

#' @rdname mm_lmm-methods
#' @export
formula.mm_lmm <- function(x, ...) {
  x$formula
}

#' @rdname mm_lmm-methods
#' @export
formula.mm_glmm <- formula.mm_lmm

#' @rdname mm_lmm-methods
#' @export
model.frame.mm_lmm <- function(formula, ...) {
  formula$model_frame
}

#' @rdname mm_lmm-methods
#' @export
model.frame.mm_glmm <- model.frame.mm_lmm

#' Number of groups per random-effect grouping factor
#'
#' `ngrps()` returns a named integer vector giving the number of levels of each
#' random-effect grouping factor, mirroring `lme4::ngrps()`.
#'
#' @return A named integer vector of group counts.
#' @rdname mm_lmm-methods
#' @export
ngrps <- function(object, ...) {
  UseMethod("ngrps")
}

#' @rdname mm_lmm-methods
#' @export
ngrps.default <- function(object, ...) {
  mm_abort(
    message = "`ngrps()` has no method for this object.",
    class = "mm_arg_error",
    input = object
  )
}

#' @rdname mm_lmm-methods
#' @export
ngrps.mm_lmm <- function(object, ...) {
  re <- object$random_effects %||% list()
  vapply(re, nrow, integer(1L))
}

#' @rdname mm_lmm-methods
#' @export
ngrps.mm_glmm <- ngrps.mm_lmm

#' @rdname mm_lmm-methods
#' @importFrom stats weights
#' @export
weights.mm_lmm <- function(object, ...) {
  object$weights
}

#' @rdname mm_lmm-methods
#' @export
weights.mm_glmm <- weights.mm_lmm

#' @rdname mm_lmm-methods
#' @importFrom stats extractAIC
#' @export
extractAIC.mm_lmm <- function(fit, scale, k = 2, ...) {
  edf <- fit$dof
  c(edf, -2 * fit$logLik + k * edf)
}

#' @rdname mm_lmm-methods
#' @export
extractAIC.mm_glmm <- extractAIC.mm_lmm

#' @rdname mm_lmm-methods
#' @importFrom stats terms
#' @export
terms.mm_lmm <- function(x, ...) {
  stats::terms(mm_fixed_formula(x))
}

#' @rdname mm_lmm-methods
#' @export
terms.mm_glmm <- terms.mm_lmm

#' Coerce variance components to an lme4-style data frame
#'
#' Produces the long form returned by `as.data.frame(lme4::VarCorr(.))`:
#' one row per variance (`var2 = NA`) and one row per covariance
#' (`var1`, `var2` both set), with a final `Residual` row for LMMs. `vcov`
#' holds the (co)variance and `sdcor` the standard deviation (diagonal) or
#' correlation (off-diagonal). This is the shape `broom.mixed::tidy()` expects.
#'
#' @param row.names,optional Ignored; present for S3 consistency.
#' @rdname mm_lmm-methods
#' @export
as.data.frame.mm_varcorr <- function(x, row.names = NULL, optional = FALSE,
                                     ...) {
  grp <- character(0)
  var1 <- character(0)
  var2 <- character(0)
  vcov <- numeric(0)
  sdcor <- numeric(0)
  for (comp in x$components_raw %||% list()) {
    nm <- gsub(": *", "", comp$names)  # lme4 concatenates: "sex: female" -> "sexfemale"
    sd <- comp$std_dev
    corr <- comp$correlations
    p <- length(nm)
    # Variance (diagonal) rows first, matching lme4's ordering.
    for (i in seq_len(p)) {
      grp <- c(grp, comp$group)
      var1 <- c(var1, nm[i])
      var2 <- c(var2, NA_character_)
      vcov <- c(vcov, sd[i]^2)
      sdcor <- c(sdcor, sd[i])
    }
    # Covariance (off-diagonal) rows; correlations are stored row-major in the
    # strict lower triangle (see mm_varcorr_correlation_values()).
    for (i in seq_len(p)) {
      offset <- (i - 1L) * (i - 2L) / 2L
      for (j in seq_len(i - 1L)) {
        r <- corr[offset + j]
        grp <- c(grp, comp$group)
        var1 <- c(var1, nm[j])
        var2 <- c(var2, nm[i])
        vcov <- c(vcov, r * sd[i] * sd[j])
        sdcor <- c(sdcor, r)
      }
    }
  }
  if (!is.null(x$residual_sd) && length(x$residual_sd) == 1L &&
      is.finite(x$residual_sd)) {
    grp <- c(grp, "Residual")
    var1 <- c(var1, NA_character_)
    var2 <- c(var2, NA_character_)
    vcov <- c(vcov, x$residual_sd^2)
    sdcor <- c(sdcor, x$residual_sd)
  }
  data.frame(grp = grp, var1 = var1, var2 = var2, vcov = vcov, sdcor = sdcor,
             stringsAsFactors = FALSE)
}

#' Coerce conditional modes to an lme4-style data frame
#'
#' Produces the long form returned by `as.data.frame(lme4::ranef(.))`: columns
#' `grpvar`, `term`, `grp`, `condval`, and `condsd`. `condsd` is the
#' conditional standard deviation, taken from the `postVar` attribute when the
#' modes were extracted with `condVar = TRUE`, and `NA` otherwise.
#'
#' @rdname mm_lmm-methods
#' @export
as.data.frame.mm_ranef <- function(x, row.names = NULL, optional = FALSE, ...) {
  grpvar <- character(0)
  term <- character(0)
  grp <- character(0)
  condval <- numeric(0)
  condsd <- numeric(0)
  for (g in names(x)) {
    df <- x[[g]]
    pv <- attr(df, "postVar")
    terms <- colnames(df)
    levs <- rownames(df)
    for (ti in seq_along(terms)) {
      for (li in seq_along(levs)) {
        sdv <- if (!is.null(pv)) sqrt(pv[ti, ti, li]) else NA_real_
        grpvar <- c(grpvar, g)
        term <- c(term, terms[ti])
        grp <- c(grp, levs[li])
        condval <- c(condval, df[li, ti])
        condsd <- c(condsd, sdv)
      }
    }
  }
  data.frame(grpvar = grpvar, term = term, grp = grp, condval = condval,
             condsd = condsd, stringsAsFactors = FALSE)
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

mm_varcorr_from_result <- function(varcorr, artifact = NULL) {
  components_in <- varcorr$components %||% list()
  group_labels <- mm_varcorr_group_labels(components_in, artifact)
  # Keep the raw (full-precision) per-group covariance pieces so
  # as.data.frame.mm_varcorr() can reconstruct the lme4 long form
  # (grp/var1/var2/vcov/sdcor) without reshaping `table`.
  components_raw <- lapply(seq_along(components_in), function(i) {
    component <- components_in[[i]]
    list(
      group = group_labels[[i]],
      names = as.character(unlist(component$names, use.names = FALSE)),
      std_dev = as.numeric(unlist(component$std_dev, use.names = FALSE)),
      correlations = as.numeric(unlist(component$correlations, use.names = FALSE))
    )
  })
  # Correlations are stored as full-precision numerics, one column per
  # preceding term within the group: `correlation` holds each term's
  # correlation with the group's first term, `correlation2` with the second,
  # and so on (NA where no such pair exists). Rounding to 2 d.p. is a
  # presentation concern handled by mm_varcorr_correlation_display().
  n_corr_cols <- max(
    1L,
    vapply(
      components_in,
      function(component) {
        length(unlist(component$names, use.names = FALSE)) - 1L
      },
      integer(1),
      USE.NAMES = FALSE
    ),
    na.rm = TRUE
  )
  corr_col_names <- mm_varcorr_correlation_col_names(n_corr_cols)
  components <- lapply(seq_along(components_in), function(component_index) {
    component <- components_in[[component_index]]
    names <- as.character(unlist(component$names, use.names = FALSE))
    std_dev <- as.numeric(unlist(component$std_dev, use.names = FALSE))
    correlations <- as.numeric(unlist(component$correlations, use.names = FALSE))
    rows <- vector("list", length(names))
    for (i in seq_along(names)) {
      row <- data.frame(
        group = group_labels[[component_index]],
        name = names[[i]],
        variance = std_dev[[i]]^2,
        std_dev = std_dev[[i]],
        stringsAsFactors = FALSE
      )
      vals <- mm_varcorr_correlation_values(correlations, i)
      row[corr_col_names] <- as.list(c(vals, rep(NA_real_, n_corr_cols - length(vals))))
      rows[[i]] <- row
    }
    do.call(rbind, rows)
  })
  table <- if (length(components)) {
    out <- do.call(rbind, components)
    rownames(out) <- NULL
    out
  } else {
    out <- data.frame(
      group = character(),
      name = character(),
      variance = numeric(),
      std_dev = numeric(),
      stringsAsFactors = FALSE
    )
    out$correlation <- numeric()
    out
  }
  residual_sd <- as.numeric(varcorr$residual_sd %||% NA_real_)
  table$boundary <- mm_varcorr_boundary_flag(table$std_dev, residual_sd)

  out <- list(
    table = table,
    residual_sd = residual_sd,
    components_raw = components_raw
  )
  class(out) <- c("mm_varcorr", "list")
  out
}

# Map each varcorr component's OWN group to its lme4-style label. The engine
# spells interaction/cell groups "a & b" (GroupingFactor::Interaction joins
# with " & "); lme4 writes "a:b". Labels must come from the component itself:
# the engine emits varcorr components in optimizer order, which need not
# match `semantic_model$random_terms` (formula) order, so any positional
# pairing of components with semantic terms swaps labels between crossed
# groups (caught by the test-bw-lme-tutorial.R crossed intercept cases).
mm_varcorr_group_labels <- function(components, artifact = NULL) {
  vapply(
    components,
    function(component) {
      mm_group_lme4_label(component$group,
                          fallback = as.character(component$group %||% ""))
    },
    character(1)
  )
}

mm_group_lme4_label <- function(group, fallback = "") {
  if (is.character(group) && length(group) == 1L && nzchar(group)) {
    return(gsub(" & ", ":", group, fixed = TRUE))
  }
  if (is.list(group)) {
    if (!is.null(group$single$name)) {
      return(as.character(group$single$name))
    }
    if (!is.null(group$cell$names)) {
      return(paste(as.character(unlist(group$cell$names, use.names = FALSE)),
                   collapse = ":"))
    }
    if (!is.null(group$interaction$names)) {
      return(paste(as.character(unlist(group$interaction$names, use.names = FALSE)),
                   collapse = ":"))
    }
  }
  as.character(fallback %||% "")
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

# Row `row_index` of a k-term group correlates with terms 1..row_index-1;
# the engine stores those values row-major in the strict lower triangle.
mm_varcorr_correlation_values <- function(correlations, row_index) {
  if (row_index <= 1L || !length(correlations)) {
    return(numeric())
  }
  offset <- (row_index - 1L) * (row_index - 2L) / 2L
  correlations[seq.int(offset + 1L, offset + row_index - 1L)]
}

mm_varcorr_correlation_col_names <- function(n) {
  if (n < 1L) {
    return(character())
  }
  c("correlation", if (n > 1L) paste0("correlation", seq.int(2L, n)))
}

# Presentation form of the numeric correlation columns: each row's non-NA
# correlations rendered "%+.2f" and space-joined, "" when the row has none.
# Shared by print.mm_varcorr() and the reporting layer so the stored table
# can stay full-precision numeric.
mm_varcorr_correlation_display <- function(table) {
  cols <- grep("^correlation[0-9]*$", names(table), value = TRUE)
  if (!length(cols) || !nrow(table)) {
    return(character(nrow(table)))
  }
  vapply(
    seq_len(nrow(table)),
    function(i) {
      vals <- as.numeric(table[i, cols])
      vals <- vals[!is.na(vals)]
      if (!length(vals)) "" else paste(sprintf("%+.2f", vals), collapse = " ")
    },
    character(1)
  )
}
