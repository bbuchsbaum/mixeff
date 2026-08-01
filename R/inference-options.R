#' Inspect inference methods available for this fit
#'
#' `inference_options()` is the audit verb for fixed-effect inference. It does
#' not run any test; it predicts, from the fit's metadata, which inference
#' methods will succeed on this fit and at what approximate cost. The goal is
#' to remove trial-and-error: a user reading the table can see which routes
#' are immediately available, which will refuse and why, and which require a
#' bootstrap.
#'
#' Like [random_options()], this function does not rank or recommend. There is
#' no "best method" row.
#'
#' @param fit A fitted `mm_lmm`.
#' @param term Optional fixed-effect term name. Reserved for future
#'   per-term refinement; currently unused (the table is fit-level).
#' @param nsim Bootstrap replicate count to use when estimating cost. Used
#'   only to format the `approx_cost` column.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_inference_options` object with a `table` data frame of one
#'   row per candidate method.
#'
#' @export
inference_options <- function(fit, term = NULL, nsim = 1000L, ...) {
  UseMethod("inference_options")
}

#' @rdname inference_options
#' @export
inference_options.mm_lmm <- function(fit, term = NULL, nsim = 1000L, ...) {
  if (!is.null(term)) {
    fixed_terms <- mm_fixed_effect_terms(fit)
    if (!term %in% fixed_terms) {
      mm_abort(
        message = sprintf(
          "Unknown fixed-effect term: %s. Known terms: %s",
          term, paste(fixed_terms, collapse = ", ")),
        class = "mm_arg_error",
        input = term
      )
    }
  }

  fs <- fit$fit_status %||% "unknown"
  is_boundary <- fs %in% c("converged_boundary", "converged_reduced_rank")
  is_reml <- isTRUE(fit$REML)

  # Number of grouping levels in the largest random factor (proxy for cluster
  # bootstrap viability)
  n_groups_max <- mm_inference_options_n_groups_max(fit)

  # What auto landed on for this fit
  current_inf <- inference_table(fit)$table
  current_method <- if (nrow(current_inf)) current_inf$method[[1L]] else NA_character_

  rows <- list(
    mm_inference_options_row_wald(fit, current_inf, is_boundary),
    mm_inference_options_row_satterthwaite(is_boundary),
    mm_inference_options_row_kenward_roger(is_boundary),
    mm_inference_options_row_bootstrap(fit, nsim),
    mm_inference_options_row_bootstrap_lrt(fit, is_reml, nsim),
    mm_inference_options_row_cluster_bootstrap(n_groups_max, nsim),
    mm_inference_options_row_profile_ci(fit, is_boundary, is_reml)
  )
  tab <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  tab$current <- tab$method == current_method
  tab <- mm_inference_options_add_display(tab)
  rownames(tab) <- NULL

  obj <- list(
    table = tab,
    fit_status = fs,
    is_reml = is_reml,
    n_groups_max = n_groups_max,
    term = term
  )
  class(obj) <- "mm_inference_options"
  obj
}

#' @rdname inference_options
#' @export
inference_options.mm_glmm <- function(fit, term = NULL, nsim = 1000L, ...) {
  # The GLMM route menu (Audit1 DEC-1 tier 2): this table is where the
  # working-Hessian approximation is DISCOVERED -- refusal messages name only
  # the certified paths plus a pointer here, so the caveat and the invocation
  # syntax travel together with the option.
  cap <- mm_glmm_inference_capability(fit)
  certified <- isTRUE(cap$wald) &&
    identical(cap$source, "certified_inference_table")
  opted_in <- identical(as.character(fit$inference_request %||% "auto"),
                        "working_hessian")

  is_nb <- identical(as.character(fit$family$family %||% ""),
                     "negative_binomial")
  rows <- list(
    list(
      method = "asymptotic_wald_z",
      expected_status = if (certified) "available" else "not_assessed",
      expected_reliability_reason = if (certified) {
        "glmm_certified_wald"
      } else if (is_nb) {
        "nb_joint_laplace_unavailable"
      } else {
        "glmm_wald_uncertified_for_profiled_estimator"
      },
      r_verb = if (certified) {
        "summary(fit)"
      } else if (is_nb) {
        # The joint route refuses NB, so advising it here would be a dead
        # end; the bootstrap row below is NB's interval path.
        "confint(fit, method = \"bootstrap\")"
      } else {
        "glmm(..., method = \"joint_laplace\")"
      },
      approx_cost = if (certified) "immediate" else "one refit (slower estimator)",
      notes = if (certified) {
        "certified Wald z from the engine inference table"
      } else if (is_nb) {
        "certified Wald is unavailable for negative binomial (no joint route); use the bootstrap"
      } else {
        "the default profiled estimator withholds Wald inference; the joint route certifies it"
      }
    ),
    list(
      method = "glmm_parametric_bootstrap",
      expected_status = "available",
      expected_reliability_reason = "bootstrap_monte_carlo_replicates",
      r_verb = "confint(fit, method = \"bootstrap\", nsim = 999)",
      approx_cost = "nsim model refits",
      notes = paste0(
        "percentile intervals from simulate-and-refit replicates of the ",
        "effective estimator; every supported family incl. negative ",
        "binomial; fixed-theta NB conditions replicates on that theta, ",
        "accounting and MCSE attached to the result"
      )
    ),
    list(
      method = "wald_z_working_hessian",
      expected_status = if (opted_in) "available" else "opt_in",
      expected_reliability_reason = "working_hessian_uncertified_approximation",
      r_verb = "glmm(..., inference = \"working_hessian\")",
      approx_cost = "immediate",
      notes = paste0(
        "UNCERTIFIED approximation, reliability moderate; SEs ran ~11% ",
        "smaller than glmer on the package's reference data ",
        "(anti-conservative). Exploration and screening, not reporting."
      )
    )
  )
  tab <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  tab$current <- tab$method == cap$method
  tab <- mm_inference_options_add_display(tab)
  rownames(tab) <- NULL

  obj <- list(
    table = tab,
    fit_status = fit$fit_status %||% "unknown",
    is_reml = FALSE,
    n_groups_max = mm_inference_options_n_groups_max(fit),
    term = term
  )
  class(obj) <- "mm_inference_options"
  obj
}

#' @method print mm_inference_options
#' @export
print.mm_inference_options <- function(x, ...) {
  cat(sprintf("Inference options (fit_status: %s, REML: %s):\n",
              x$fit_status, x$is_reml))
  print(x$table[, c("method", "display_status", "display_reason",
                    "what_to_do_next", "approx_cost", "current")],
        row.names = FALSE)
  cat("\nUse `<obj>$table` for raw enum columns (`expected_status`, `expected_reliability_reason`) and notes.\n")
  invisible(x)
}

# ---- helpers ----

mm_inference_options_n_groups_max <- function(fit) {
  re_levels <- vapply(ranef(fit) %||% list(), function(x) {
    if (is.data.frame(x)) nrow(x) else NA_integer_
  }, integer(1L))
  if (!length(re_levels) || all(is.na(re_levels))) return(NA_integer_)
  max(re_levels, na.rm = TRUE)
}

mm_inference_options_row_wald <- function(fit, current_inf, is_boundary) {
  predicted <- if (is_boundary)
    "asymptotic_wald_z_at_boundary"
  else
    "interior_converged"
  observed <- if (nrow(current_inf) &&
                  identical(current_inf$method[[1L]], "asymptotic_wald_z") &&
                  !is.na(current_inf$reliability_reason[[1L]])) {
    current_inf$reliability_reason[[1L]]
  } else {
    NA_character_
  }
  list(
    method = "asymptotic_wald_z",
    expected_status = "available",
    expected_reliability_reason = if (!is.na(observed)) observed else predicted,
    r_verb = "summary(fit)",
    approx_cost = "immediate",
    notes = if (is_boundary)
      "boundary fit; SEs may understate uncertainty"
    else
      "default fast route"
  )
}

mm_inference_options_row_satterthwaite <- function(is_boundary) {
  list(
    method = "satterthwaite",
    expected_status = if (is_boundary) "not_assessed" else "available",
    expected_reliability_reason = if (is_boundary)
      "satterthwaite_unavailable_at_boundary"
    else
      "satterthwaite_finite_difference_approximation",
    r_verb = "test_effect(fit, term, method = 'satterthwaite')",
    approx_cost = "immediate",
    notes = if (is_boundary)
      "boundary varpar; df derivative undefined"
    else
      "df-corrected t test"
  )
}

mm_inference_options_row_kenward_roger <- function(is_boundary) {
  list(
    method = "kenward_roger",
    expected_status = if (is_boundary) "not_assessed" else "available",
    expected_reliability_reason = if (is_boundary)
      "kenward_roger_unavailable_at_boundary"
    else
      "kenward_roger_approximation",
    r_verb = "test_effect(fit, term, method = 'kenward_roger')",
    approx_cost = "immediate",
    notes = if (is_boundary)
      "boundary varpar; df derivative undefined"
    else
      "df + small-sample variance adjustment"
  )
}

mm_inference_options_row_bootstrap <- function(fit, nsim) {
  list(
    method = "bootstrap",
    expected_status = "available",
    expected_reliability_reason = "bootstrap_monte_carlo_replicates",
    r_verb = sprintf(
      "test_effect(fit, term, method = 'bootstrap', bootstrap = bootstrap_control(nsim = %d))",
      nsim),
    approx_cost = mm_inference_options_format_cost(fit, nsim, factor = 1),
    notes = paste(
      "model-conditioned; respects fitted covariance.",
      "Single-df terms produce a t row; multi-df terms produce a joint Wald-F",
      "row with `num_df = effective restriction rank`."
    )
  )
}

mm_inference_options_row_bootstrap_lrt <- function(fit, is_reml, nsim) {
  list(
    method = "bootstrap_lrt",
    expected_status = if (is_reml) "not_assessed" else "available",
    expected_reliability_reason = if (is_reml)
      "bootstrap_lrt_requires_ml"
    else
      mm_inference_options_bootstrap_lrt_reliability_reason(nsim),
    r_verb = sprintf(
      "test_effect(fit, term, method = 'bootstrap_lrt', bootstrap = bootstrap_control(nsim = %d))",
      nsim),
    approx_cost = mm_inference_options_format_cost(fit, nsim, factor = 2),
    notes = if (is_reml)
      "ML required; refit with `lmm(..., REML = FALSE)` to enable bootstrap LRT"
    else
      "model-comparison LRT; refits reduced and alternative per replicate"
  )
}

mm_inference_options_bootstrap_lrt_reliability_reason <- function(nsim) {
  if (is.numeric(nsim) && length(nsim) == 1L && !is.na(nsim) &&
      nsim >= mm_bootstrap_moderate_min()) {
    "bootstrap_monte_carlo_replicates"
  } else {
    "bootstrap_insufficient_replicates"
  }
}

mm_inference_options_row_cluster_bootstrap <- function(n_groups_max, nsim) {
  too_few <- !is.na(n_groups_max) && n_groups_max < 10L
  list(
    method = "cluster_bootstrap",
    expected_status = "not_assessed",
    expected_reliability_reason = "bootstrap_cluster_resample_p_value_unavailable",
    r_verb = "test_effect(fit, term, method = 'cluster_bootstrap', group = '<group>')",
    approx_cost = "-",
    notes = if (too_few)
      sprintf(paste0("heuristic (rule of thumb, not a validated threshold): ",
                     "resampling only %d clusters gives the cluster ",
                     "bootstrap little granularity"),
              n_groups_max)
    else
      "cluster_resample is an estimator-distribution target; fixed-effect p-values are not certified"
  )
}

mm_inference_options_row_profile_ci <- function(fit, is_boundary, is_reml) {
  refused <- is_boundary || is_reml
  reason <- if (is_boundary) {
    "profile_ci_unavailable_at_boundary"
  } else if (is_reml) {
    "profile_beta_unavailable_under_reml"
  } else {
    "profile_likelihood_ci"
  }
  list(
    method = "profile_ci",
    expected_status = if (refused) "not_assessed" else "available",
    expected_reliability_reason = reason,
    r_verb = "confint(fit, method = 'profile')",
    approx_cost = "slow profile refits",
    notes = if (is_boundary) {
      "profile intervals are not certified for boundary or reduced-rank fits"
    } else if (is_reml) {
      "REML profile payloads omit fixed-effect beta intervals; refit ML for beta profile CIs"
    } else {
      "profile-likelihood confidence intervals; slower than Wald"
    }
  )
}

mm_inference_options_add_display <- function(tab) {
  tab$display_status <- vapply(
    tab$expected_status,
    mm_inference_options_display_status,
    character(1)
  )
  tab$display_reason <- mapply(
    mm_inference_options_display_reason,
    tab$expected_reliability_reason,
    tab$expected_status,
    tab$method,
    USE.NAMES = FALSE
  )
  tab$what_to_do_next <- mapply(
    mm_inference_options_next_step,
    tab$method,
    tab$expected_status,
    tab$expected_reliability_reason,
    tab$r_verb,
    USE.NAMES = FALSE
  )
  tab
}

mm_inference_options_display_status <- function(status) {
  switch(
    status,
    available = "runs now",
    not_assessed = "refused on this fit",
    not_yet_wired = "available upstream; R bridge pending",
    opt_in = "opt-in (uncertified)",
    "check route status"
  )
}

mm_inference_options_display_reason <- function(reason, status, method) {
  if (is.na(reason) || !nzchar(reason) || identical(reason, "not_available")) {
    return("-")
  }
  switch(
    reason,
    asymptotic_wald_z_at_boundary =
      "Wald route runs, but boundary fits can understate uncertainty",
    interior_converged =
      "interior fit; default fast route",
    satterthwaite_unavailable_at_boundary =
      "variance-parameter derivative undefined at boundary",
    kenward_roger_unavailable_at_boundary =
      "variance-parameter derivative undefined at boundary",
    satterthwaite_finite_difference_approximation =
      "finite-sample t test via finite-difference degrees of freedom",
    kenward_roger_approximation =
      "finite-sample degrees of freedom with covariance adjustment",
    bootstrap_monte_carlo_replicates =
      "calibrated by nsim and Monte Carlo error",
    bootstrap_insufficient_replicates =
      "increase nsim for a moderate reliability grade",
    bootstrap_lrt_requires_ml =
      "requires an ML fit; refit with REML = FALSE",
    bootstrap_cluster_resample_p_value_unavailable =
      "cluster resampling reports estimator distributions, not fixed-effect p-values",
    profile_likelihood_ci =
      "profile-likelihood confidence interval",
    profile_beta_unavailable_under_reml =
      "fixed-effect profile intervals require an ML fit",
    profile_ci_unavailable_at_boundary =
      "profile intervals are not certified at the boundary",
    gsub("_", " ", reason, fixed = TRUE)
  )
}

mm_inference_options_next_step <- function(method, status, reason, r_verb) {
  if (identical(status, "available")) {
    return(r_verb)
  }
  if (identical(reason, "bootstrap_lrt_requires_ml")) {
    return(sprintf("Refit with lmm(..., REML = FALSE), then run %s", r_verb))
  }
  if (reason %in% c("satterthwaite_unavailable_at_boundary",
                    "kenward_roger_unavailable_at_boundary",
                    "profile_ci_unavailable_at_boundary")) {
    return("Use asymptotic_wald_z or bootstrap; simplify the random-effects structure if the boundary is unintended")
  }
  if (identical(reason, "bootstrap_cluster_resample_p_value_unavailable")) {
    return("Use bootstrap or bootstrap_lrt for fixed-effect p-values")
  }
  if (identical(reason, "profile_beta_unavailable_under_reml")) {
    return("Refit with lmm(..., REML = FALSE), then run confint(fit, method = 'profile')")
  }
  if (identical(status, "not_yet_wired")) {
    return("Track the R bridge wiring before using this route")
  }
  "Inspect the raw status columns before using this route"
}

mm_inference_options_format_cost <- function(fit, nsim, factor = 1) {
  # Report WORKLOAD, not invented seconds: per-replicate cost depends on the
  # model, the machine, and the build, and a two-branch constant dressed as
  # a runtime reads far more precise than it is (Audit1 WI-5.3). One
  # bootstrap replicate costs one model refit; the bootstrap LRT refits the
  # reduced and alternative models per replicate.
  refits <- as.integer(round(nsim * factor))
  sprintf("%s model refits (nsim = %d)",
          format(refits, big.mark = ","), as.integer(nsim))
}
