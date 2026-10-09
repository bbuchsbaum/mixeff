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
#'   `"aic"` for information criteria only, `"bootstrap"` for a small
#'   parametric-bootstrap LRT when `nsim > 0`, or `"kenward_roger"` /
#'   `"satterthwaite"` for pbkrtest's `KRmodcomp()` / `SATmodcomp()`: an F
#'   test of the larger model's fixed effects restricted to the smaller
#'   model's (two models with identical random effects, nested in their fixed
#'   effects). The restriction matrix is built from the two design matrices
#'   as pbkrtest does; the test runs on the REML fit of the larger model
#'   (refitted by REML if it was fitted by ML, as pbkrtest does). The
#'   engine's Kenward-Roger F is the unscaled statistic (pbkrtest's `FtestU`
#'   row: same F and denominator df; the KR scaling factor is not applied
#'   and is reported as `f_scaling = 1`). The F row
#'   replaces the LRT row (`statistic_name = "F"`, `df` = numerator df,
#'   `den_df` = denominator df); the full result is in `$fixed_f`.
#' @param refit_for_comparison How to handle REML fits.
#' @param nsim Number of bootstrap simulations for `method = "bootstrap"`.
#' @param seed Optional bootstrap seed.
#' @param threads Worker threads for the bootstrap refits (default 1; see
#'   [bootstrap_control()]). The result is identical for every value.
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
                           method = c("auto", "lrt", "bootstrap", "aic",
                                      "kenward_roger", "satterthwaite"),
                           refit_for_comparison = c("auto", "error", "ml"),
                           nsim = 0L,
                           seed = NULL,
                           threads = 1L) {
  target <- match.arg(target)
  method <- match.arg(method)
  refit_for_comparison <- match.arg(refit_for_comparison)
  if (identical(method, "bootstrap")) {
    if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 0) {
      mm_abort(
        message = "`nsim` must be a non-negative integer for `compare(method = \"bootstrap\")`.",
        class = "mm_arg_error",
        input = nsim
      )
    }
    nsim <- as.integer(nsim)
  }
  fits <- c(list(object), list(...))
  if (!all(vapply(fits, inherits, logical(1), what = "mm_lmm"))) {
    mm_abort(
      message = "`compare()` requires fitted `mm_lmm` objects.",
      class = "mm_arg_error",
      input = fits
    )
  }
  mm_assert_comparable_lmm(fits)
  fixed_f <- NULL
  if (method %in% c("kenward_roger", "satterthwaite")) {
    if (length(fits) != 2L) {
      mm_abort(
        message = sprintf("`method = \"%s\"` compares exactly two nested models.", method),
        class = "mm_arg_error",
        input = length(fits)
      )
    }
    ord <- order(vapply(fits, function(x) length(x$beta), numeric(1)))
    fixed_f <- mm_compare_fixed_f(fits[[ord[[1L]]]], fits[[ord[[2L]]]], method)
  }
  prepared <- mm_prepare_comparison_fits(fits, refit_for_comparison)
  fits <- prepared$fits
  table <- mm_compare_table(
    fits, if (is.null(fixed_f)) method else "auto", prepared$refit
  )
  if (!is.null(fixed_f)) {
    table <- mm_apply_fixed_f_to_table(table, fixed_f, method)
  }
  bootstrap <- NULL
  if (identical(method, "bootstrap") && nsim > 0L && length(fits) == 2L) {
    bootstrap <- parametric_bootstrap(
      fits[[1L]],
      fits[[2L]],
      nsim = nsim,
      seed = seed,
      threads = threads
    )
    last <- nrow(table)
    table$p_value[last] <- bootstrap$p_value
    table$method[last] <- "parametric_bootstrap_lrt"
    table$status[last] <- bootstrap$status
    succ <- bootstrap$successful_replicates
    table$reason[last] <- if (identical(bootstrap$status, "available")) {
      sprintf("parametric bootstrap LRT (%s/%d replicates%s)",
              if (is.na(succ)) "?" else as.character(succ),
              nsim,
              if (is.na(bootstrap$mcse)) "" else sprintf(", MCSE=%.4g", bootstrap$mcse))
    } else {
      bootstrap$reason %||% "parametric bootstrap LRT did not certify a p-value"
    }
  } else if (identical(method, "bootstrap")) {
    # Bootstrap was requested but cannot run: do not let the asymptotic LRT
    # p-value masquerade as the requested bootstrap result.
    table$p_value <- NA_real_
    table$method <- "bootstrap_not_run"
    table$status <- "bootstrap_not_run"
    table$reason <- "set nsim > 0 and compare exactly two models to run bootstrap"
  }
  ledger <- mm_comparison_ledger(
    table,
    target = target,
    requested_method = method,
    refit_for_comparison = refit_for_comparison,
    source = "mixeff.compare"
  )
  obj <- list(
    table = table,
    ledger = ledger,
    fits = fits,
    target = target,
    method = method,
    refit_for_comparison = refit_for_comparison,
    bootstrap = bootstrap,
    fixed_f = fixed_f
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
#' Runs the engine-certified parametric-bootstrap likelihood-ratio test
#' between two nested ML-fitted LMMs through the Rust
#' `mm_bootstrap_lrt_json` entry point. The smaller model (fewer estimated
#' parameters) is the reduced model; the larger is the alternative. The
#' returned object carries the engine's replicate accounting (successful and
#' completed replicates, boundary count, Monte-Carlo standard error, seed)
#' rather than a bare `mean()` p-value, so every reported number traces back
#' to a versioned Rust payload.
#'
#' The engine refuses REML fits: refit with `lmm(..., REML = FALSE)` before
#' calling. (`compare(method = "bootstrap")` refits REML to ML automatically.)
#'
#' @param null,alternative Fitted `mm_lmm` objects. Order is irrelevant; the
#'   model with fewer parameters is treated as the reduced model.
#' @param nsim Number of bootstrap replicates.
#' @param seed Optional bootstrap seed.
#' @param threads Worker threads for the replicate refits (default 1). The
#'   responses are simulated serially and each replicate refits fresh copies
#'   of both models, so the result is identical for every value; see
#'   [bootstrap_control()] for the threading contract.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_parametric_bootstrap` object.
#'
#' @export
parametric_bootstrap <- function(null, alternative, nsim = 100L, seed = NULL,
                                 threads = 1L, ...) {
  if (!inherits(null, "mm_lmm") || !inherits(alternative, "mm_lmm")) {
    mm_abort(
      message = "`parametric_bootstrap()` requires two fitted `mm_lmm` objects.",
      class = "mm_arg_error",
      input = list(null = null, alternative = alternative)
    )
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    mm_abort(
      message = "`nsim` must be a positive integer.",
      class = "mm_arg_error",
      input = nsim
    )
  }
  nsim <- as.integer(nsim)
  # Reduced = fewer estimated parameters; the engine LRT is direction-aware
  # so we order explicitly rather than trust the call order.
  if (alternative$dof < null$dof) {
    tmp <- null
    null <- alternative
    alternative <- tmp
  }
  if (isTRUE(null$REML) || isTRUE(alternative$REML)) {
    mm_abort(
      message = paste(
        "parametric bootstrap likelihood-ratio test requires ML-fitted",
        "models; refit with `lmm(..., REML = FALSE)` and retry"
      ),
      class = "mm_inference_unavailable",
      reason_code = "bootstrap_lrt_requires_ml",
      input = list(null_reml = isTRUE(null$REML),
                   alternative_reml = isTRUE(alternative$REML))
    )
  }
  mm_assert_bootstrap_lrt_pair(null, alternative)
  bootstrap <- bootstrap_control(nsim = nsim, seed = seed, threads = threads)
  bridge <- mm_rust_fit_bridge_payload(alternative)
  bootstrap_json <- jsonlite::toJSON(
    mm_bootstrap_wire(bootstrap),
    auto_unbox = TRUE,
    null = "null"
  )
  json <- tryCatch(
    mm_bootstrap_lrt_json(
      deparse1(null$formula),
      bridge$formula_string,
      bridge$spec_data$column_order,
      bridge$spec_data$numeric_columns,
      bridge$spec_data$categorical_values,
      bridge$spec_data$categorical_levels,
      bridge$spec_data$categorical_ordered,
      bridge$weights,
      bridge$control_json,
      as.character(bootstrap_json)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort(
      message = conditionMessage(json),
      class = "mm_inference_unavailable",
      reason_code = "bootstrap_lrt_engine_refused",
      parent = json
    )
  }
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  payload <- parsed$payload %||% list()
  meta <- payload$metadata %||% list()
  simulated <- as.numeric(unlist(payload$replicate_statistics %||% list(),
                                 use.names = FALSE))
  certified <- !is.null(parsed$p_value)
  out <- list(
    observed = parsed$observed_statistic %||% mm_lrt_stat(null, alternative),
    simulated = simulated,
    p_value = parsed$p_value %||% NA_real_,
    nsim = nsim,
    threads = bootstrap$threads,
    successful_replicates = meta$successful_replicates %||% NA_integer_,
    completed_replicates = meta$completed_replicates %||% NA_integer_,
    boundary_count = meta$boundary_count %||% NA_integer_,
    mcse = parsed$mcse %||% NA_real_,
    seed = meta$seed_record$seed %||% seed,
    status = if (certified) "available" else "not_assessed",
    reason = if (certified) {
      NA_character_
    } else {
      "parametric bootstrap likelihood-ratio test did not certify a p-value"
    },
    notes = as.character(unlist(parsed$notes %||% list(), use.names = FALSE)),
    reduced_formula = deparse1(null$formula),
    alternative_formula = bridge$formula_string
  )
  class(out) <- "mm_parametric_bootstrap"
  out
}

#' @method print mm_parametric_bootstrap
#' @export
print.mm_parametric_bootstrap <- function(x, ...) {
  fmt_int <- function(v) if (is.null(v) || is.na(v)) "NA" else format(as.integer(v))
  fmt_num <- function(v) if (is.null(v) || !is.finite(v)) "NA" else sprintf("%.6g", v)
  cat("Parametric bootstrap likelihood-ratio test:\n")
  cat(sprintf("  status:    %s\n", x$status %||% "unknown"))
  cat(sprintf("  observed:  %s\n", fmt_num(x$observed)))
  cat(sprintf("  requested replicates: %s\n", fmt_int(x$nsim)))
  cat(sprintf("  successful / completed: %s / %s\n",
              fmt_int(x$successful_replicates),
              fmt_int(x$completed_replicates)))
  if (!is.null(x$boundary_count) && !is.na(x$boundary_count) &&
      x$boundary_count > 0L) {
    cat(sprintf("  boundary replicates: %s\n", fmt_int(x$boundary_count)))
  }
  cat(sprintf("  MCSE:      %s\n", fmt_num(x$mcse)))
  cat(sprintf("  seed:      %s\n", fmt_int(x$seed)))
  if (identical(x$status, "available")) {
    cat(sprintf("  p.value:   %s\n", fmt_num(x$p_value)))
  } else {
    cat(sprintf("  p.value:   not certified -- %s\n",
                x$reason %||% "no reason recorded"))
  }
  notes <- x$notes %||% character()
  if (length(notes)) {
    cat("  notes:\n")
    for (n in notes) cat(sprintf("    - %s\n", n))
  }
  invisible(x)
}

#' Analysis of variance for a mixeff LMM
#'
#' With no extra models, `anova()` returns a term-level F-test table for
#' `object` shaped like lmerTest's: a data frame of class `"anova"` with
#' columns `Sum Sq`, `Mean Sq`, `NumDF`, `DenDF`, `F value` and `Pr(>F)`
#' (rows named by term; `Mean Sq = F * sigma^2` and `Sum Sq = NumDF * Mean
#' Sq`, as in lmerTest). With further fitted models in `...`, it returns
#' lme4's likelihood-ratio table (`npar`, `AIC`, `BIC`, `logLik`,
#' `-2*log(L)`, `Chisq`, `Df`, `Pr(>Chisq)`, rows named after the arguments;
#' lme4 before 2.0 called the `-2*log(L)` column `deviance`),
#' computed by [compare()] (REML fits are refit by ML, as in lme4).
#'
#' mixeff's provenance stays on the result: `x$table` is the full
#' term-level (or [compare()]) table with method, status, reliability and
#' reason columns, and `x$type` / `x$requested_method` record the request.
#' Rows the engine could not certify keep `NA` statistics and their reason
#' is printed under the table.
#'
#' @param object A fitted `mm_lmm`.
#' @param ... Optional additional fitted models; triggers [compare()].
#' @param type Term-hypothesis type. `"III"` (default) tests marginal
#'   Type III hypotheses: a term's own contrast columns plus the
#'   equally-weighted average over the levels of every term that contains
#'   it, which makes the test invariant to which factor level is the
#'   reference and matches SAS / `car` / `lmerTest` Type III. `"II"`
#'   respects marginality (each term adjusted for all terms that do not
#'   contain it). `"I"` is sequential in `terms()` order (main effects,
#'   then two-way interactions, and so on). `"block"` tests the raw
#'   coefficient block for each term — under treatment coding that is the
#'   simple effect at the other factors' reference levels, which is a
#'   legitimate quantity but is *not* Type III on unbalanced designs; it
#'   is the hypothesis mixeff computed for `"III"` before engine
#'   `1f3f689`. lmerTest's spellings `3`, `2`, `1` (numeric or character)
#'   are accepted.
#' @param method Degrees-of-freedom / test method.
#' @param refit_for_comparison Passed to [compare()] when `...` is used.
#' @param ddf lmerTest's spelling of the method: `"Satterthwaite"`
#'   (`method = "satterthwaite"`), `"Kenward-Roger"`
#'   (`method = "kenward_roger"`), or `"lme4"`, which returns lme4's own
#'   single-model table (sequential `npar`, `Sum Sq`, `Mean Sq`, `F value`,
#'   no denominator df or p-values). Supply `ddf` or `method`, not both.
#'
#' @return A data frame of class `c("mm_anova", "anova", "data.frame")` for
#'   one model, or `c("mm_anova_comparison", "mm_model_comparison",
#'   "anova", "data.frame")` for several.
#'
#' @method anova mm_lmm
#' @export
anova.mm_lmm <- function(object, ..., type = c("III", "II", "I", "block"),
                         method = c("auto", "satterthwaite", "kenward_roger",
                                    "bootstrap", "asymptotic", "none"),
                         refit_for_comparison = c("auto", "error", "ml"),
                         ddf = NULL) {
  dots <- list(...)
  lme4_table <- FALSE
  if (!is.null(ddf)) {
    if (!missing(method)) {
      mm_abort(
        message = "Supply either `method` or lmerTest's `ddf`, not both.",
        class = "mm_arg_error",
        input = ddf
      )
    }
    mapped <- mm_ddf_to_method(ddf)
    method <- mapped$method
    lme4_table <- mapped$lme4
  }
  if (length(dots)) {
    cmp_method <- if (identical(match.arg(method), "bootstrap")) "bootstrap" else "auto"
    labels <- vapply(as.list(match.call(expand.dots = FALSE)$...),
                     deparse1, character(1))
    labels <- c(deparse1(substitute(object)), labels)
    cmp <- compare(
      object,
      ...,
      method = cmp_method,
      refit_for_comparison = match.arg(refit_for_comparison)
    )
    return(mm_anova_comparison_frame(cmp, c(list(object), dots), labels))
  }
  type <- mm_anova_type(type)
  method <- match.arg(method)
  refit_for_comparison <- match.arg(refit_for_comparison)
  if (lme4_table) {
    return(mm_lme4_sequential_anova(object))
  }
  terms <- setdiff(mm_fixed_effect_terms(object), "1")
  if (identical(method, "none")) {
    table <- mm_unavailable_effect_table(mm_user_term_label(object, terms), method)
  } else {
    parsed <- mm_rust_term_table(object, method, type = type)
    table <- parsed$table[parsed$table$term %in% terms, , drop = FALSE]
    # Expanded formula terms: report the user's term names (a multi-column
    # term such as poly(x, 2) is tested per column by the engine).
    table$term <- mm_user_term_label(object, table$term)
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reason", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
    # Single-df terms come back as t-form rows; lme4's anova presents the
    # mathematically identical F form (F(1, nu) = t(nu)^2). Convert so the
    # ANOVA table reads like lme4's and num_df is never a bare NA for an
    # available row.
    t_rows <- !is.na(table$statistic_name) & table$statistic_name == "t"
    if (any(t_rows)) {
      table$statistic[t_rows] <- table$statistic[t_rows]^2
      table$statistic_name[t_rows] <- "f"
      table$num_df[t_rows] <- 1
    }
    # Asymptotic single-df rows are Wald z: F(1, Inf) = z^2.
    z_rows <- !is.na(table$statistic_name) & table$statistic_name == "z"
    if (any(z_rows)) {
      table$statistic[z_rows] <- table$statistic[z_rows]^2
      table$statistic_name[z_rows] <- "f"
      table$num_df[z_rows] <- 1
      table$den_df[z_rows] <- Inf
    }
  }
  if (!identical(method, "none")) {
    table <- mm_anova_group_expanded_terms(object, table, type)
  }
  table$type <- type
  table <- table[, c("term", "type", setdiff(names(table), c("term", "type"))),
                 drop = FALSE]
  rownames(table) <- NULL
  mm_anova_frame(object, table, type, method, refit_for_comparison)
}

# A user term that expands into several engine columns (poly(x, 2),
# ns(x, 3), ...) is tested per column by the engine. lmerTest tests it as one
# multi-df term. For a term no other term contains (so Type I/II/III all test
# exactly its coefficients), replace the per-column rows by one joint F test
# of all the term's coefficients, with the same df method as the rows.
mm_anova_group_expanded_terms <- function(fit, table, type) {
  ex <- fit$expansion
  if (is.null(ex) || !length(ex$term_map)) return(table)
  multi <- names(ex$term_map)[lengths(ex$term_map) > 1L]
  if (!length(multi)) return(table)
  droppable <- mm_droppable_terms(fit)
  X <- stats::model.matrix(fit, type = "fixed")
  for (lab in intersect(multi, droppable)) {
    pieces <- which(startsWith(table$term, lab))
    if (length(pieces) < 2L) next
    method <- unique(table$method[pieces])
    if (length(method) != 1L ||
        !method %in% c("satterthwaite", "kenward_roger")) next
    # Each per-column row is labelled with its coefficient's name.
    cols <- match(table$term[pieces], colnames(X))
    if (anyNA(cols)) next
    L <- matrix(0, nrow = length(cols), ncol = ncol(X),
                dimnames = list(colnames(X)[cols], colnames(X)))
    L[cbind(seq_along(cols), cols)] <- 1
    joint <- tryCatch(mm_anova_joint_f(fit, L, method, lab),
                      error = function(cnd) NULL)
    if (is.null(joint)) next
    row <- table[pieces[[1L]], , drop = FALSE]
    row$term <- lab
    row$num_df <- joint$num_df
    row$den_df <- joint$den_df
    row$statistic <- joint$statistic
    row$statistic_name <- "f"
    row$p_value <- joint$p_value
    row$status <- joint$status
    table <- rbind(table[seq_len(pieces[[1L]] - 1L), , drop = FALSE], row,
                   table[-seq_len(max(pieces)), , drop = FALSE])
  }
  table
}

mm_anova_joint_f <- function(fit, L, method, label) {
  bridge <- mm_rust_fit_bridge_payload(fit)
  json <- mm_fixed_effect_joint_test_json(
    bridge$formula_string,
    isTRUE(fit$REML),
    bridge$spec_data$column_order,
    bridge$spec_data$numeric_columns,
    bridge$spec_data$categorical_values,
    bridge$spec_data$categorical_levels,
    bridge$spec_data$categorical_ordered,
    bridge$weights,
    bridge$control_json,
    as.numeric(t(mm_coef_l_to_engine(L, fit))),
    as.integer(nrow(L)),
    as.integer(ncol(L)),
    label,
    rep(0, nrow(L)),
    method
  )
  row <- mm_json_parse_fixed_effect_inference_table(
    jsonlite::fromJSON(json, simplifyVector = FALSE)
  )$table[1L, , drop = FALSE]
  list(
    statistic = as.numeric(row$statistic),
    num_df = as.numeric(row$numerator_df),
    den_df = as.numeric(row$denominator_df),
    p_value = as.numeric(row$p_value),
    status = as.character(row$status)
  )
}

# lmerTest's ddf= spellings -> mixeff method names.
mm_ddf_to_method <- function(ddf) {
  choices <- c("Satterthwaite", "Kenward-Roger", "lme4")
  if (!is.character(ddf) || length(ddf) != 1L || is.na(ddf)) {
    mm_abort(
      message = "`ddf` must be one of \"Satterthwaite\", \"Kenward-Roger\", or \"lme4\".",
      class = "mm_arg_error",
      input = ddf
    )
  }
  hit <- choices[pmatch(tolower(ddf), tolower(choices))]
  if (is.na(hit)) {
    mm_abort(
      message = sprintf(
        "Unknown `ddf = \"%s\"`; use \"Satterthwaite\", \"Kenward-Roger\", or \"lme4\".",
        ddf
      ),
      class = "mm_arg_error",
      input = ddf
    )
  }
  switch(
    hit,
    Satterthwaite = list(method = "satterthwaite", lme4 = FALSE),
    `Kenward-Roger` = list(method = "kenward_roger", lme4 = FALSE),
    lme4 = list(method = "asymptotic", lme4 = TRUE)
  )
}

mm_anova_type <- function(type) {
  if (length(type) > 1L) return("III")
  type <- as.character(type)
  map <- c("3" = "III", "2" = "II", "1" = "I", III = "III", II = "II",
           I = "I", block = "block", marginal = "III")
  if (!type %in% names(map)) {
    mm_abort(
      message = sprintf(
        "Unknown anova `type = \"%s\"`; use \"III\", \"II\", \"I\" (or 3, 2, 1) or \"block\".",
        type
      ),
      class = "mm_arg_error",
      input = type
    )
  }
  unname(map[[type]])
}

mm_anova_method_label <- function(table, method) {
  used <- unique(stats::na.omit(table$method[table$status %in% "available"]))
  m <- if (length(used) == 1L) used else method
  switch(
    m,
    satterthwaite = "Satterthwaite's method",
    kenward_roger = "Kenward-Roger's method",
    asymptotic_wald_z = "asymptotic Wald (chi-square/df) method",
    asymptotic = "asymptotic Wald (chi-square/df) method",
    sprintf("method `%s`", m)
  )
}

# lmerTest-shaped single-model table; the full provenance table rides along
# as an attribute (reachable as x$table).
mm_anova_frame <- function(object, table, type, method, refit_for_comparison) {
  F <- as.numeric(table$statistic)
  F[!(table$statistic_name %in% "f")] <- NA_real_
  num <- as.numeric(table$num_df)
  s2 <- as.numeric(sigma(object))^2
  ms <- F * s2
  out <- data.frame(
    `Sum Sq` = ms * num,
    `Mean Sq` = ms,
    NumDF = num,
    DenDF = as.numeric(table$den_df),
    `F value` = F,
    `Pr(>F)` = as.numeric(table$p_value),
    check.names = FALSE
  )
  row.names(out) <- table$term
  heading <- sprintf("Type %s Analysis of Variance Table with %s",
                     type, mm_anova_method_label(table, method))
  structure(
    out,
    heading = heading,
    mm_table = table,
    mm_type = type,
    mm_requested_method = method,
    mm_refit_for_comparison = refit_for_comparison,
    class = c("mm_anova", "anova", "data.frame")
  )
}

# lme4::anova(<lmerMod>) with one model (lmerTest's ddf = "lme4"): sequential
# sums of squares from the RX factor, no denominator df.
mm_lme4_sequential_anova <- function(object) {
  X <- stats::model.matrix(object, type = "fixed")
  asgn <- attr(X, "assign")
  beta <- object$beta
  idx <- match(names(beta), colnames(X))
  if (is.null(asgn) || anyNA(idx)) {
    mm_abort(
      message = "lme4's sequential ANOVA table cannot be rebuilt for this fit.",
      class = "mm_inference_unavailable",
      reason_code = "lme4_anova_unavailable"
    )
  }
  asgn <- asgn[idx]
  RX <- tryCatch(getME(object, "RX"), error = function(cnd) cnd)
  if (inherits(RX, "condition")) {
    mm_abort(
      message = paste("lme4's sequential ANOVA table needs a full-rank",
                      "fixed-effect design:", conditionMessage(RX)),
      class = "mm_inference_unavailable",
      reason_code = "lme4_anova_unavailable"
    )
  }
  ss <- as.vector(RX %*% beta)^2
  labels <- attr(stats::terms(mm_fixed_formula(object)), "term.labels")
  nm <- c(if (any(asgn == 0L)) "(Intercept)", labels[unique(asgn[asgn > 0L])])
  ss <- unlist(lapply(split(ss, asgn), sum))
  df <- lengths(split(asgn, asgn))
  ms <- ss / df
  table <- data.frame(npar = as.integer(df), `Sum Sq` = ss, `Mean Sq` = ms,
                      `F value` = ms / as.numeric(sigma(object))^2,
                      check.names = FALSE)
  row.names(table) <- nm
  if ("(Intercept)" %in% nm) {
    table <- table[-match("(Intercept)", nm), , drop = FALSE]
  }
  structure(table, heading = "Analysis of Variance Table",
            class = c("anova", "data.frame"))
}

# lme4-shaped multi-model LRT table built from a compare() result.
mm_anova_comparison_frame <- function(cmp, fits, labels) {
  tbl <- cmp$table
  formulas <- vapply(fits, function(f) deparse1(f$formula), character(1))
  # compare() orders models by complexity; map rows back to argument labels.
  row_labels <- labels
  if (length(labels) == nrow(tbl) && "formula" %in% names(tbl)) {
    fit_forms <- vapply(cmp$fits, function(f) mm_compare_formula_label(f),
                        character(1))
    orig_forms <- vapply(fits, function(f) mm_compare_formula_label(f),
                         character(1))
    hit <- match(fit_forms, orig_forms)
    if (!anyNA(hit) && !anyDuplicated(hit)) row_labels <- labels[hit]
    formulas <- formulas[if (!anyNA(hit) && !anyDuplicated(hit)) hit else
      seq_along(formulas)]
  }
  if (anyDuplicated(row_labels)) row_labels <- paste0("MODEL", seq_along(row_labels))
  loglik <- as.numeric(tbl$logLik)
  out <- data.frame(
    npar = as.numeric(tbl$df),
    AIC = as.numeric(tbl$AIC),
    BIC = as.numeric(tbl$BIC),
    logLik = loglik,
    `-2*log(L)` = -2 * loglik,
    Chisq = as.numeric(tbl$LRT),
    Df = as.numeric(tbl$delta_df),
    `Pr(>Chisq)` = as.numeric(tbl$p_value),
    check.names = FALSE
  )
  row.names(out) <- row_labels
  data_expr <- fits[[1L]]$call$data
  heading <- c(
    mm_anova_data_line(data_expr),
    "Models:",
    paste(row_labels, formulas, sep = ": ")
  )
  structure(
    out,
    heading = heading,
    mm_comparison = cmp,
    # mm_model_comparison stays in the class so reporting_table() and code
    # reading `$table` / `$ledger` keep working.
    class = c("mm_anova_comparison", "mm_model_comparison", "anova",
              "data.frame")
  )
}

# lme4's "Data: <expr>" heading line. A call holding literal data (as
# update() leaves after refitting the stored model frame) is not shown.
mm_anova_data_line <- function(data_expr) {
  if (is.null(data_expr)) return(NULL)
  if (!is.name(data_expr) && !is.call(data_expr)) return(NULL)
  text <- deparse1(data_expr)
  if (nchar(text) > 80L) return(NULL)
  paste("Data:", text)
}

mm_compare_formula_label <- function(fit) {
  paste(deparse1(fit$formula), length(fit$beta), fit$dof)
}

#' @export
`$.mm_anova` <- function(x, name) {
  if (name %in% names(x)) return(.subset2(x, name))
  attr(x, paste0("mm_", name), exact = TRUE)
}

#' @export
`$.mm_anova_comparison` <- function(x, name) {
  if (name %in% names(x)) return(.subset2(x, name))
  cmp <- attr(x, "mm_comparison", exact = TRUE)
  if (is.null(cmp)) return(NULL)
  cmp[[name]]
}

#' @method print mm_anova
#' @export
print.mm_anova <- function(x, ...) {
  NextMethod()
  table <- x$table
  if (!is.null(table) && "status" %in% names(table)) {
    bad <- !is.na(table$status) & table$status != "available"
    if (any(bad)) {
      cat("\nNot certified (statistics withheld):\n")
      for (i in which(bad)) {
        cat(sprintf("  %s: %s%s\n", table$term[[i]], table$status[[i]],
                    if (!is.na(table$reason[[i]]) && nzchar(table$reason[[i]]))
                      paste0(" -- ", table$reason[[i]]) else ""))
      }
    }
  }
  invisible(x)
}

#' @method print mm_anova_comparison
#' @export
print.mm_anova_comparison <- function(x, ...) {
  # Print as stats' anova table (skip the mm_*_comparison print methods,
  # which are kept in the class for `$table` / reporting_table()).
  shown <- x
  attr(shown, "mm_comparison") <- NULL
  class(shown) <- c("anova", "data.frame")
  print(shown, ...)
  table <- x$table
  if (!is.null(table) && "status" %in% names(table)) {
    bad <- !is.na(table$status) &
      !table$status %in% c("available", "reference_model")
    if (any(bad)) {
      cat("\nNot certified:\n")
      for (i in which(bad)) {
        cat(sprintf("  row %d: %s%s\n", i, table$status[[i]],
                    if (!is.na(table$reason[[i]]) && nzchar(table$reason[[i]]))
                      paste0(" -- ", table$reason[[i]]) else ""))
      }
    }
  }
  invisible(x)
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
  terms <- mm_drop1_terms(object)
  if (!is.null(scope)) {
    terms <- intersect(terms, as.character(scope))
  } else {
    # Match stats::drop1 marginality semantics: by default only terms not
    # contained in a higher-order term are droppable (you don't test a main
    # effect while its interaction is in the model). An explicit `scope` may
    # still request a non-marginal drop -- the engine fits those correctly
    # (matching lme4's full-dummy expansion), so they produce ordinary rows.
    terms <- intersect(terms, mm_droppable_terms(object))
  }
  prepared_full <- mm_prepare_comparison_fits(list(object), refit_for_comparison)
  full <- prepared_full$fits[[1L]]
  full_refit <- isTRUE(prepared_full$refit[[1L]])
  rows <- lapply(terms, function(term) {
    reduced_formula <- mm_drop_fixed_term_formula(full, term)
    reduced <- tryCatch(
      mm_internal_lmm(reduced_formula, full$model_frame,
                      REML = isTRUE(full$REML), weights = full$weights,
                      control = mm_internal_control(full, keep_start = FALSE)),
      error = function(cnd) cnd
    )
    if (inherits(reduced, "condition")) {
      # Explicit-scope non-marginal drops (or any refit refusal) surface as an
      # unavailable row instead of aborting the whole table.
      return(data.frame(
        dropped = term,
        formula = mm_drop_formula_label(full, term, reduced_formula),
        df = NA_real_,
        logLik = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        LRT = NA_real_,
        p_value = NA_real_,
        method = "unavailable",
        status = "unavailable",
        reason = conditionMessage(reduced),
        stringsAsFactors = FALSE
      ))
    }
    stat <- mm_lrt_stat(reduced, full)
    df <- full$dof - reduced$dof
    data.frame(
      dropped = term,
      formula = mm_drop_formula_label(full, term, reduced_formula),
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
      status = "available",
      reason = NA_character_,
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
  ledger <- mm_drop1_comparison_ledger(
    full = full,
    table = table,
    test = test,
    refit_for_comparison = refit_for_comparison,
    full_refit = full_refit
  )
  obj <- list(table = table, ledger = ledger, full = full)
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
    # lme4's anova() refuses the same way ("models were not all fitted to
    # the same size of dataset"); typically different missing-value patterns
    # dropped different rows under `na.action`.
    mm_abort(
      message = sprintf(
        paste0(
          "Models were not all fitted to the same size of dataset ",
          "(nobs: %s). Likelihood comparisons need the same rows; this usually ",
          "means `na.action` dropped different incomplete rows for different ",
          "models. Refit all models to the same complete-case data (e.g. ",
          "na.omit() over the union of their variables)."
        ),
        paste(n, collapse = ", ")
      ),
      class = "mm_arg_error",
      reason_code = "different_nobs",
      input = n
    )
  }
  responses <- vapply(fits, mm_response_name, character(1))
  if (length(unique(responses)) != 1L) {
    mm_abort(
      message = "Compared models must use the same response variable.",
      class = "mm_arg_error",
      input = responses
    )
  }
}

mm_assert_bootstrap_lrt_pair <- function(null, alternative) {
  if (!isTRUE(alternative$dof > null$dof)) {
    mm_abort(
      message = paste(
        "Parametric bootstrap LRT requires nested models with the",
        "alternative estimating more parameters than the reduced model."
      ),
      class = "mm_arg_error",
      reason_code = "bootstrap_lrt_requires_nested_models",
      input = list(null_df = null$dof, alternative_df = alternative$dof)
    )
  }

  null_vars <- all.vars(null$formula)
  alt_names <- names(alternative$model_frame %||% data.frame())
  missing <- setdiff(null_vars, alt_names)
  if (length(missing)) {
    mm_abort(
      message = sprintf(
        "Parametric bootstrap LRT requires the reduced model variables to be present in the alternative model frame; missing: %s.",
        paste(missing, collapse = ", ")
      ),
      class = "mm_arg_error",
      reason_code = "bootstrap_lrt_requires_nested_model_frames",
      input = missing
    )
  }

  shared <- intersect(names(null$model_frame %||% data.frame()), alt_names)
  mismatched <- shared[!vapply(shared, function(nm) {
    identical(null$model_frame[[nm]], alternative$model_frame[[nm]])
  }, logical(1))]
  if (length(mismatched)) {
    mm_abort(
      message = sprintf(
        "Parametric bootstrap LRT requires compared fits to share identical model-frame values; mismatched column(s): %s.",
        paste(mismatched, collapse = ", ")
      ),
      class = "mm_arg_error",
      reason_code = "bootstrap_lrt_requires_same_observations",
      input = mismatched
    )
  }

  if (!identical(null$weights, alternative$weights)) {
    mm_abort(
      message = "Parametric bootstrap LRT requires compared fits to use identical case weights.",
      class = "mm_arg_error",
      reason_code = "bootstrap_lrt_requires_same_weights",
      input = list(null = null$weights, alternative = alternative$weights)
    )
  }
  invisible(TRUE)
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
      # Same model, ML criterion: keep the user's mm_control() (optimizer,
      # tolerances, start), silenced.
      mm_internal_lmm(fit$formula, fit$model_frame, REML = FALSE,
                      weights = fit$weights, offset = mm_fit_offset_arg(fit),
                      control = mm_internal_control(fit))
    })
    refit <- has_reml
  }
  list(fits = fits, refit = refit)
}

mm_compare_table <- function(fits, method, refit) {
  ord <- order(vapply(fits, function(x) x$dof, numeric(1)))
  fits <- fits[ord]
  refit <- refit[ord]
  payloads <- lapply(fits, function(fit) {
    payload <- mm_rust_fit_bridge_payload(fit)
    payload$REML <- isTRUE(fit$REML)
    payload
  })
  json <- tryCatch(
    mm_compare_models_json(payloads, method, "never"),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = method)
  }
  parsed <- mm_json_parse_model_comparison_table(json)
  table <- mm_compare_table_from_rust_payload(parsed$payload, fits, refit, method)
  rownames(table) <- NULL
  table
}

mm_json_parse_model_comparison_table <- function(json) {
  if (!is.character(json) || length(json) != 1L || is.na(json) || !nzchar(json)) {
    mm_abort(
      message = "`json` must be a single non-empty character string.",
      class = "mm_schema_error",
      input = json
    )
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse model-comparison JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  schema <- parsed$schema
  if (!is.list(schema) ||
      !identical(as.character(schema$schema_name), "mixedmodels.model_comparison_table") ||
      !identical(as.character(schema$schema_version), "1.0.0")) {
    mm_abort(
      message = "Model-comparison JSON has an unknown schema header.",
      class = "mm_schema_error",
      input = parsed
    )
  }
  if (!is.list(parsed$payload) || !is.list(parsed$payload$rows)) {
    mm_abort(
      message = "Model-comparison JSON is missing its row payload.",
      class = "mm_schema_error",
      input = parsed
    )
  }
  parsed
}

mm_compare_table_from_rust_payload <- function(payload, fits, refit, method) {
  rows <- payload$rows
  n <- length(rows)
  if (n != length(fits)) {
    mm_abort(
      message = "Model-comparison row count does not match compared fits.",
      class = "mm_schema_error",
      input = list(rows = n, fits = length(fits))
    )
  }

  scalar <- function(row, field, default) {
    value <- row[[field]]
    if (is.null(value)) default else value
  }
  chr <- function(field, default = NA_character_) {
    vapply(rows, function(row) as.character(scalar(row, field, default)),
           character(1))
  }
  num <- function(field, default = NA_real_) {
    vapply(rows, function(row) as.numeric(scalar(row, field, default)),
           numeric(1))
  }
  int <- function(field, default = NA_integer_) {
    as.integer(vapply(rows, function(row) as.integer(scalar(row, field, default)),
                      integer(1)))
  }
  bool <- function(field, default = FALSE) {
    vapply(rows, function(row) isTRUE(scalar(row, field, default)), logical(1))
  }

  lrt_available <- bool("lrt_available")
  requires_ml_refit <- bool("requires_ml_refit")
  status <- rep("not_available", n)
  status[seq_len(n) == 1L] <- "reference_model"
  status[lrt_available] <- "available"
  status[requires_ml_refit] <- "ml_refit_required"
  if (identical(method, "aic")) {
    status[] <- "information_criteria"
  }

  method_col <- rep("not_available", n)
  method_col[seq_len(n) == 1L | lrt_available] <- "asymptotic_lrt"
  if (identical(method, "aic")) {
    method_col[] <- "none"
  }

  reason <- chr("reason", "")
  reason[is.na(reason)] <- ""

  data.frame(
    model = paste0("m", seq_len(n)),
    formula = chr("label", ""),
    nobs = int("nobs"),
    df = int("dof"),
    logLik = num("loglik"),
    deviance = num("deviance"),
    AIC = num("aic"),
    BIC = num("bic"),
    delta_aic = num("delta_aic"),
    delta_bic = num("delta_bic"),
    REML = vapply(fits, function(fit) isTRUE(fit$REML), logical(1)),
    refit = vapply(refit, isTRUE, logical(1)),
    fit_status = vapply(fits, mm_fit_status_label, character(1)),
    delta_df = num("chisq_dof"),
    LRT = num("chisq"),
    p_value = num("pvalue"),
    method = method_col,
    status = status,
    reason = reason,
    reason_code = chr("reason_code"),
    comparison_class = chr("comparison_class"),
    lrt_available = lrt_available,
    information_criteria_available = bool("information_criteria_available", TRUE),
    requires_ml_refit = requires_ml_refit,
    loglik_within_optimizer_tol = bool("loglik_within_optimizer_tol", NA),
    rust_method = as.character(payload$method %||% NA_character_),
    rust_refit_policy = as.character(payload$refit_policy %||% NA_character_),
    stringsAsFactors = FALSE
  )
}

mm_comparison_ledger <- function(table, target, requested_method,
                                 refit_for_comparison,
                                 source = "mixeff.compare") {
  n <- nrow(table)
  if (!n) {
    return(mm_comparison_ledger_empty())
  }
  status <- mm_table_col(table, "status", "not_available")
  reason <- mm_comparison_reason(
    status,
    mm_table_col(table, "reason", NA_character_),
    mm_table_col(table, "reason_code", NA_character_)
  )
  reml <- mm_logical_col(table, "REML", FALSE)
  refit <- mm_logical_col(table, "refit", FALSE)
  data.frame(
    comparison_id = rep(mm_comparison_id(table$formula, target, requested_method), n),
    model_id = mm_table_col(table, "model", paste0("m", seq_len(n))),
    model_index = seq_len(n),
    model_role = ifelse(seq_len(n) == 1L, "reference", "candidate"),
    formula = as.character(table$formula),
    fit_method = ifelse(reml, "REML", "ML"),
    original_fit_method = ifelse(refit, "REML", ifelse(reml, "REML", "ML")),
    refit = refit,
    refit_policy = rep(refit_for_comparison, n),
    comparison_target = rep(target, n),
    requested_method = rep(requested_method, n),
    comparison_method = mm_table_col(table, "method", requested_method),
    statistic = as.numeric(mm_table_col(table, "LRT", NA_real_)),
    statistic_name = ifelse(is.na(mm_table_col(table, "LRT", NA_real_)),
                            NA_character_, "LRT"),
    df = as.numeric(mm_table_col(table, "delta_df", NA_real_)),
    p_value = as.numeric(mm_table_col(table, "p_value", NA_real_)),
    nobs = as.integer(mm_table_col(table, "nobs", NA_integer_)),
    logLik = as.numeric(mm_table_col(table, "logLik", NA_real_)),
    AIC = as.numeric(mm_table_col(table, "AIC", NA_real_)),
    BIC = as.numeric(mm_table_col(table, "BIC", NA_real_)),
    dof = as.numeric(mm_table_col(table, "df", NA_real_)),
    delta_aic = as.numeric(mm_table_col(table, "delta_aic", NA_real_)),
    delta_bic = as.numeric(mm_table_col(table, "delta_bic", NA_real_)),
    fit_status = mm_table_col(table, "fit_status", "not_assessed"),
    validity_status = status,
    status = status,
    reason = reason,
    reason_code = mm_table_col(table, "reason_code", NA_character_),
    comparison_class = mm_table_col(table, "comparison_class", NA_character_),
    lrt_available = mm_logical_col(table, "lrt_available", FALSE),
    information_criteria_available = mm_logical_col(
      table, "information_criteria_available", TRUE
    ),
    requires_ml_refit = mm_logical_col(table, "requires_ml_refit", FALSE),
    loglik_within_optimizer_tol = mm_table_col(table, "loglik_within_optimizer_tol", NA),
    report_row = seq_len(n),
    source = rep(source, n),
    stringsAsFactors = FALSE
  )
}

mm_drop1_comparison_ledger <- function(full, table, test, refit_for_comparison,
                                       full_refit) {
  n <- nrow(table)
  if (!n) {
    return(mm_comparison_ledger_empty())
  }
  status <- if (identical(test, "Chisq")) {
    ifelse(is.finite(table$LRT) & !is.na(table$p_value), "available",
           "not_available")
  } else {
    rep("information_criteria", n)
  }
  reason <- ifelse(
    identical(test, "Chisq") | status != "information_criteria",
    NA_character_,
    "single-term deletion reported without likelihood-ratio test"
  )
  data.frame(
    comparison_id = rep(mm_comparison_id(c(deparse1(full$formula), table$formula),
                                         "fixed_effects", test), n),
    model_id = paste0("drop_", seq_len(n)),
    model_index = seq_len(n),
    model_role = "reduced_candidate",
    dropped = as.character(table$dropped),
    formula = as.character(table$formula),
    reference_formula = rep(deparse1(full$formula), n),
    fit_method = ifelse(isTRUE(full$REML), "REML", "ML"),
    original_fit_method = ifelse(isTRUE(full_refit), "REML",
                                 ifelse(isTRUE(full$REML), "REML", "ML")),
    refit = rep(isTRUE(full_refit), n),
    refit_policy = rep(refit_for_comparison, n),
    comparison_target = "fixed_effects",
    requested_method = rep(test, n),
    comparison_method = as.character(table$method),
    statistic = as.numeric(table$LRT),
    statistic_name = ifelse(identical(test, "Chisq"), "LRT", NA_character_),
    df = as.numeric(table$df),
    p_value = as.numeric(table$p_value),
    nobs = rep(nobs(full), n),
    logLik = as.numeric(table$logLik),
    AIC = as.numeric(table$AIC),
    BIC = as.numeric(table$BIC),
    dof = NA_real_,
    delta_aic = NA_real_,
    delta_bic = NA_real_,
    fit_status = rep(mm_fit_status_label(full), n),
    validity_status = status,
    status = status,
    reason = mm_comparison_reason(status, reason, NA_character_),
    reason_code = NA_character_,
    comparison_class = "drop1_fixed_effect",
    lrt_available = identical(test, "Chisq") & is.finite(table$LRT),
    information_criteria_available = TRUE,
    requires_ml_refit = rep(isTRUE(full_refit), n),
    loglik_within_optimizer_tol = NA,
    report_row = seq_len(n),
    source = "mixeff.drop1",
    stringsAsFactors = FALSE
  )
}

mm_comparison_ledger_empty <- function() {
  data.frame(
    comparison_id = character(),
    model_id = character(),
    model_index = integer(),
    model_role = character(),
    formula = character(),
    fit_method = character(),
    original_fit_method = character(),
    refit = logical(),
    refit_policy = character(),
    comparison_target = character(),
    requested_method = character(),
    comparison_method = character(),
    statistic = numeric(),
    statistic_name = character(),
    df = numeric(),
    p_value = numeric(),
    nobs = integer(),
    logLik = numeric(),
    AIC = numeric(),
    BIC = numeric(),
    dof = numeric(),
    delta_aic = numeric(),
    delta_bic = numeric(),
    fit_status = character(),
    validity_status = character(),
    status = character(),
    reason = character(),
    reason_code = character(),
    comparison_class = character(),
    lrt_available = logical(),
    information_criteria_available = logical(),
    requires_ml_refit = logical(),
    loglik_within_optimizer_tol = logical(),
    report_row = integer(),
    source = character(),
    stringsAsFactors = FALSE
  )
}

mm_comparison_reason <- function(status, reason, reason_code) {
  status <- as.character(status)
  reason <- as.character(reason)
  reason_code <- rep(as.character(reason_code), length.out = length(status))

  missing <- is.na(reason) | !nzchar(reason)
  reason[missing & status == "reference_model"] <- "baseline model for comparison"
  reason[missing & status == "information_criteria"] <- "information criteria row"

  missing <- is.na(reason) | !nzchar(reason)
  use_code <- missing & !is.na(reason_code) & nzchar(reason_code)
  reason[use_code] <- reason_code[use_code]

  missing <- is.na(reason) | !nzchar(reason)
  needs_reason <- !status %in% c("available", "reference_model",
                                 "information_criteria")
  reason[missing & needs_reason] <- "comparison status not available"
  reason[reason == ""] <- NA_character_
  reason
}

mm_table_col <- function(table, col, default) {
  if (col %in% names(table)) {
    table[[col]]
  } else {
    rep(default, nrow(table))
  }
}

mm_logical_col <- function(table, col, default) {
  vapply(mm_table_col(table, col, default), isTRUE, logical(1))
}

mm_comparison_id <- function(formulas, target, method) {
  text <- paste(c(target, method, formulas), collapse = "\r")
  codes <- utf8ToInt(enc2utf8(text))
  if (!length(codes)) {
    return("cmp_00000000")
  }
  checksum <- sum((seq_along(codes) %% 997L) * codes) %% 1000000007
  sprintf("cmp_%08x", as.integer(checksum))
}

mm_fit_status_label <- function(fit) {
  as.character(
    fit$fit_status %||%
      fit$artifact$optimizer_certificate$status %||%
      "not_assessed"
  )
}

mm_lrt_stat <- function(null, alternative) {
  # -2 log-likelihood difference (deviance() of a GLMM is the residual
  # deviance, not the Laplace objective).
  pmax(0, 2 * (as.numeric(alternative$logLik) - as.numeric(null$logLik)))
}

# Fixed-effect right-hand side that preserves the fit's intercept choice: a
# no-intercept model stays no-intercept (`0 + ...`) after a term is dropped.
mm_fixed_terms_have_intercept <- function(terms) {
  any(c("1", "(Intercept)") %in% terms)
}

mm_fixed_rhs_text <- function(terms, has_intercept) {
  terms <- setdiff(terms, c("1", "0", "(Intercept)"))
  if (has_intercept) {
    if (length(terms)) paste(terms, collapse = " + ") else "1"
  } else {
    paste(c("0", terms), collapse = " + ")
  }
}

mm_drop_fixed_term_formula <- function(fit, term) {
  response <- mm_response_name(fit)
  drop <- term
  if (!is.null(fit$expansion) && term %in% names(fit$expansion$term_map)) {
    # A user term (e.g. poly(x, 2)) spans several engine terms.
    engine <- setdiff(mm_fixed_effect_terms(fit), "1")
    drop <- engine[mm_term_key(engine) %in%
                     mm_term_key(fit$expansion$term_map[[term]])]
  }
  all_fixed <- mm_fixed_effect_terms(fit)
  fixed_rhs <- mm_fixed_rhs_text(setdiff(all_fixed, drop),
                                 mm_fixed_terms_have_intercept(all_fixed))
  random <- vapply(
    fit$artifact$semantic_model$random_terms %||% list(),
    function(x) x$source_syntax$text %||% "",
    character(1)
  )
  random <- random[nzchar(random)]
  rhs <- paste(c(fixed_rhs, random), collapse = " + ")
  stats::as.formula(paste(response, "~", rhs), env = environment(fit$formula))
}

# Terms droppable under stats::drop1 marginality rules: a term is droppable
# iff no OTHER term contains all of its variables (e.g. `recipe` is not
# droppable from `recipe * temperature`).
mm_droppable_terms <- function(fit) {
  tt <- mm_user_fixed_terms(fit)
  labels <- attr(tt, "term.labels")
  fac <- attr(tt, "factors")
  if (!length(labels) || is.null(dim(fac))) return(labels)
  vars_of <- lapply(seq_along(labels), function(i) rownames(fac)[fac[, i] > 0])
  droppable <- vapply(seq_along(labels), function(i) {
    !any(vapply(seq_along(labels)[-i], function(j) {
      all(vars_of[[i]] %in% vars_of[[j]])
    }, logical(1)))
  }, logical(1))
  labels[droppable]
}
