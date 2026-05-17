#' Compare fitted mixeff models
#'
#' `compare()` is the namespace-qualified model-comparison front door. For LMMs
#' it reports likelihood, information criteria, and asymptotic likelihood-ratio
#' comparisons. REML fits are refit by ML when `refit_for_comparison = "auto"`
#' or `"ml"`; `"error"` refuses that comparison.
#'
#' @param object A fitted `mm_lmm`.
#' @param ... Additional fitted `mm_lmm` objects.
#' @param target Comparison target label.
#' @param method `"auto"` / `"lrt"` for asymptotic likelihood-ratio rows,
#'   `"aic"` for information criteria only, or `"bootstrap"` for a small
#'   parametric-bootstrap LRT when `nsim > 0`.
#' @param refit_for_comparison How to handle REML fits.
#' @param nsim Number of bootstrap simulations for `method = "bootstrap"`.
#' @param seed Optional bootstrap seed.
#'
#' @return An `mm_model_comparison` object with a data-frame `table`.
#'
#' @importFrom stats anova drop1
#' @export
compare <- function(object, ...) {
  UseMethod("compare")
}

#' @rdname compare
#' @export
compare.mm_lmm <- function(object,
                           ...,
                           target = c("fixed_effects", "random_effects", "prediction"),
                           method = c("auto", "lrt", "bootstrap", "aic"),
                           refit_for_comparison = c("auto", "error", "ml"),
                           nsim = 0L,
                           seed = NULL) {
  target <- match.arg(target)
  method <- match.arg(method)
  refit_for_comparison <- match.arg(refit_for_comparison)
  fits <- c(list(object), list(...))
  if (!all(vapply(fits, inherits, logical(1), what = "mm_lmm"))) {
    mm_abort(
      message = "`compare()` requires fitted `mm_lmm` objects.",
      class = "mm_inference_unavailable",
      input = fits
    )
  }
  mm_assert_comparable_lmm(fits)
  prepared <- mm_prepare_comparison_fits(fits, refit_for_comparison)
  fits <- prepared$fits
  table <- mm_compare_table(fits, method, prepared$refit)
  bootstrap <- NULL
  if (identical(method, "bootstrap") && nsim > 0L && length(fits) == 2L) {
    bootstrap <- parametric_bootstrap(
      fits[[1L]],
      fits[[2L]],
      nsim = nsim,
      seed = seed
    )
    table$p_value[nrow(table)] <- bootstrap$p_value
    table$status[nrow(table)] <- "parametric_bootstrap"
    table$reason[nrow(table)] <- sprintf("parametric bootstrap with nsim=%d", nsim)
  } else if (identical(method, "bootstrap")) {
    table$status <- "bootstrap_not_run"
    table$reason <- "set nsim > 0 and compare exactly two models to run bootstrap"
  }
  obj <- list(
    table = table,
    fits = fits,
    target = target,
    method = method,
    refit_for_comparison = refit_for_comparison,
    bootstrap = bootstrap
  )
  class(obj) <- "mm_model_comparison"
  obj
}

#' @method print mm_model_comparison
#' @export
print.mm_model_comparison <- function(x, ...) {
  cat("Model comparison:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Parametric bootstrap likelihood-ratio comparison
#'
#' The bootstrap simulates from the smaller model, refits both models to each
#' simulated response, and compares simulated likelihood-ratio statistics with
#' the observed statistic.
#'
#' @param null,alternative Fitted `mm_lmm` objects.
#' @param nsim Number of simulations.
#' @param seed Optional random seed.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_parametric_bootstrap` object.
#'
#' @export
parametric_bootstrap <- function(null, alternative, nsim = 100L, seed = NULL, ...) {
  if (!inherits(null, "mm_lmm") || !inherits(alternative, "mm_lmm")) {
    mm_abort(
      message = "`parametric_bootstrap()` requires two fitted `mm_lmm` objects.",
      class = "mm_inference_unavailable",
      input = list(null = null, alternative = alternative)
    )
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    mm_abort(
      message = "`nsim` must be a positive integer.",
      class = "mm_inference_unavailable",
      input = nsim
    )
  }
  nsim <- as.integer(nsim)
  observed <- mm_lrt_stat(null, alternative)
  stats <- mm_with_seed(seed, {
    vapply(seq_len(nsim), function(i) {
      y <- stats::simulate(null, nsim = 1L)[[1L]]
      ref_null <- refit(null, y)
      ref_alt <- refit(alternative, y)
      mm_lrt_stat(ref_null, ref_alt)
    }, numeric(1))
  })
  p_value <- mean(stats >= observed)
  out <- list(
    observed = observed,
    simulated = stats,
    p_value = p_value,
    nsim = nsim,
    seed = seed
  )
  class(out) <- "mm_parametric_bootstrap"
  out
}

#' @method print mm_parametric_bootstrap
#' @export
print.mm_parametric_bootstrap <- function(x, ...) {
  cat("Parametric bootstrap LRT:\n")
  cat(sprintf("  observed: %.6g\n", x$observed))
  cat(sprintf("  nsim:     %d\n", x$nsim))
  cat(sprintf("  p.value:  %.6g\n", x$p_value))
  invisible(x)
}

#' @method anova mm_lmm
#' @export
anova.mm_lmm <- function(object, ..., type = c("III", "II", "I"),
                         method = c("auto", "satterthwaite", "kenward_roger",
                                    "bootstrap", "asymptotic", "none"),
                         refit_for_comparison = c("auto", "error", "ml")) {
  dots <- list(...)
  if (length(dots)) {
    cmp_method <- if (identical(match.arg(method), "bootstrap")) "bootstrap" else "auto"
    return(compare(
      object,
      ...,
      method = cmp_method,
      refit_for_comparison = match.arg(refit_for_comparison)
    ))
  }
  type <- match.arg(type)
  method <- match.arg(method)
  refit_for_comparison <- match.arg(refit_for_comparison)
  terms <- setdiff(mm_fixed_effect_terms(object), "1")
  if (identical(method, "none")) {
    table <- mm_unavailable_effect_table(terms, method)
  } else {
    parsed <- mm_rust_term_table(object, method)
    table <- parsed$table[parsed$table$term %in% terms, , drop = FALSE]
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reason", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
  }
  table$type <- type
  table <- table[, c("term", "type", setdiff(names(table), c("term", "type"))),
                 drop = FALSE]
  obj <- list(
    table = table,
    type = type,
    requested_method = method,
    refit_for_comparison = refit_for_comparison
  )
  class(obj) <- "mm_anova"
  obj
}

#' Drop one fixed-effect term at a time
#'
#' `drop1.mm_lmm()` refits reduced fixed-effect models and compares them to the
#' original fit. It is conservative: random-effect terms are preserved exactly,
#' and the reduced formulas are reported in the result table.
#'
#' @param object A fitted `mm_lmm`.
#' @param scope Optional character vector of fixed-effect terms to drop.
#' @param test Comparison test label. `"Chisq"` reports asymptotic LRT rows;
#'   `"none"` reports information criteria only.
#' @param refit_for_comparison How to handle REML fits.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_drop1` object.
#'
#' @method drop1 mm_lmm
#' @export
drop1.mm_lmm <- function(object,
                         scope = NULL,
                         test = c("none", "Chisq"),
                         refit_for_comparison = c("auto", "error", "ml"),
                         ...) {
  test <- match.arg(test)
  refit_for_comparison <- match.arg(refit_for_comparison)
  terms <- setdiff(mm_fixed_effect_terms(object), "1")
  if (!is.null(scope)) {
    terms <- intersect(terms, as.character(scope))
  }
  full <- mm_prepare_comparison_fits(list(object), refit_for_comparison)$fits[[1L]]
  rows <- lapply(terms, function(term) {
    reduced_formula <- mm_drop_fixed_term_formula(full, term)
    reduced <- lmm(reduced_formula, full$model_frame, REML = isTRUE(full$REML),
                   weights = full$weights,
                   control = mm_control(verbose = -1))
    stat <- mm_lrt_stat(reduced, full)
    df <- full$dof - reduced$dof
    data.frame(
      dropped = term,
      formula = deparse1(reduced_formula),
      df = df,
      logLik = as.numeric(logLik(reduced)),
      AIC = AIC(reduced),
      BIC = BIC(reduced),
      LRT = if (identical(test, "Chisq")) stat else NA_real_,
      p_value = if (identical(test, "Chisq") && df > 0) {
        stats::pchisq(stat, df = df, lower.tail = FALSE)
      } else {
        NA_real_
      },
      method = if (identical(test, "Chisq")) "asymptotic_lrt" else "none",
      stringsAsFactors = FALSE
    )
  })
  table <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      dropped = character(),
      formula = character(),
      df = numeric(),
      logLik = numeric(),
      AIC = numeric(),
      BIC = numeric(),
      LRT = numeric(),
      p_value = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    )
  }
  rownames(table) <- NULL
  obj <- list(table = table, full = full)
  class(obj) <- "mm_drop1"
  obj
}

#' @method print mm_drop1
#' @export
print.mm_drop1 <- function(x, ...) {
  cat("Single-term deletion table:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

mm_assert_comparable_lmm <- function(fits) {
  n <- vapply(fits, nobs, integer(1))
  if (length(unique(n)) != 1L) {
    mm_abort(
      message = "Compared models must have the same number of observations.",
      class = "mm_inference_unavailable",
      input = n
    )
  }
  responses <- vapply(fits, mm_response_name, character(1))
  if (length(unique(responses)) != 1L) {
    mm_abort(
      message = "Compared models must use the same response variable.",
      class = "mm_inference_unavailable",
      input = responses
    )
  }
}

mm_prepare_comparison_fits <- function(fits, refit_for_comparison) {
  has_reml <- vapply(fits, function(x) isTRUE(x$REML), logical(1))
  refit <- rep(FALSE, length(fits))
  if (any(has_reml)) {
    if (identical(refit_for_comparison, "error")) {
      mm_abort(
        message = "REML fits require `refit_for_comparison = \"auto\"` or `\"ml\"` for likelihood comparison.",
        class = "mm_inference_unavailable",
        input = has_reml
      )
    }
    fits <- lapply(fits, function(fit) {
      if (!isTRUE(fit$REML)) return(fit)
      lmm(fit$formula, fit$model_frame, REML = FALSE,
          weights = fit$weights,
          control = mm_control(verbose = -1))
    })
    refit <- has_reml
  }
  list(fits = fits, refit = refit)
}

mm_compare_table <- function(fits, method, refit) {
  ord <- order(vapply(fits, function(x) x$dof, numeric(1)))
  fits <- fits[ord]
  refit <- refit[ord]
  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    data.frame(
      model = paste0("m", i),
      formula = deparse1(fit$formula),
      nobs = nobs(fit),
      df = fit$dof,
      logLik = as.numeric(logLik(fit)),
      deviance = deviance(fit),
      AIC = AIC(fit),
      BIC = BIC(fit),
      REML = isTRUE(fit$REML),
      refit = isTRUE(refit[[i]]),
      delta_df = NA_real_,
      LRT = NA_real_,
      p_value = NA_real_,
      method = if (identical(method, "aic")) "none" else "asymptotic_lrt",
      status = if (identical(method, "aic")) "information_criteria" else "asymptotic_uncertified",
      reason = if (identical(method, "aic")) "" else "ordinary LRT is not a finite-sample mixed-model certificate",
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  if (!identical(method, "aic") && nrow(table) > 1L) {
    for (i in 2:nrow(table)) {
      table$delta_df[[i]] <- table$df[[i]] - table$df[[i - 1L]]
      table$LRT[[i]] <- pmax(0, table$deviance[[i - 1L]] - table$deviance[[i]])
      if (table$delta_df[[i]] > 0) {
        table$p_value[[i]] <- stats::pchisq(table$LRT[[i]],
                                            df = table$delta_df[[i]],
                                            lower.tail = FALSE)
      }
    }
  }
  rownames(table) <- NULL
  table
}

mm_lrt_stat <- function(null, alternative) {
  pmax(0, deviance(null) - deviance(alternative))
}

mm_drop_fixed_term_formula <- function(fit, term) {
  response <- mm_response_name(fit)
  fixed <- setdiff(mm_fixed_effect_terms(fit), c("1", term))
  fixed_rhs <- if (length(fixed)) paste(fixed, collapse = " + ") else "1"
  random <- vapply(
    fit$artifact$semantic_model$random_terms %||% list(),
    function(x) x$source_syntax$text %||% "",
    character(1)
  )
  random <- random[nzchar(random)]
  rhs <- paste(c(fixed_rhs, random), collapse = " + ")
  stats::as.formula(paste(response, "~", rhs), env = environment(fit$formula))
}
