#' Contrast fixed effects
#'
#' Note: this is not `emmeans::contrast`. `contrast()` is mixeff's fixed-effect
#' contrast front door. R validates the contrast matrix shape, then asks Rust
#' to evaluate estimability, method prerequisites, standard errors, degrees of
#' freedom, statistics, p-values, reliability, and unavailable reasons.
#'
#' @param fit A fitted `mm_lmm`.
#' @param L A numeric contrast vector or matrix with one column per fixed
#'   effect.
#' @param rhs Numeric right-hand side, recycled to the number of contrasts.
#' @param method Requested inference method.
#' @param bootstrap Optional [bootstrap_control()] object for
#'   `method = "bootstrap"`.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_contrast` object with a data-frame `table`. The `estimate`
#'   column is the tested difference, `L beta_hat - rhs`.
#'
#' @export
contrast <- function(fit, L, rhs = 0, method = c("auto", "satterthwaite",
                                                 "kenward_roger", "bootstrap",
                                                 "asymptotic", "boundary_lrt",
                                                 "none"),
                     bootstrap = NULL, ...) {
  UseMethod("contrast")
}

#' @rdname contrast
#' @export
contrast.mm_lmm <- function(fit, L, rhs = 0, method = c("auto", "satterthwaite",
                                                        "kenward_roger", "bootstrap",
                                                        "asymptotic", "boundary_lrt",
                                                        "none"),
                            bootstrap = NULL, ...) {
  method <- match.arg(method)
  L <- mm_contrast_matrix(L, fit)
  rhs <- rep(as.numeric(rhs), length.out = nrow(L))
  if (identical(method, "boundary_lrt")) {
    table <- mm_boundary_lrt_fixed_effect_contrast_table(L, rhs)
    obj <- list(
      table = table,
      L = L,
      rhs = rhs,
      requested_method = method,
      raw = NULL
    )
    class(obj) <- "mm_contrast"
    return(obj)
  }

  if (identical(method, "none")) {
    estimate <- as.numeric(L %*% fit$beta) - rhs
    table <- data.frame(
      contrast = rownames(L),
      estimate = estimate,
      rhs = rhs,
      std_error = NA_real_,
      df = NA_real_,
      statistic = NA_real_,
      statistic_name = NA_character_,
      p_value = NA_real_,
      method = "not_computed",
      requested_method = method,
      status = "not_assessed",
      reliability = "not_available",
      estimability = "not_assessed",
      reason = mm_inference_unavailable_reason(method),
      notes = I(rep(list(character()), nrow(L))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    obj <- list(table = table, L = L, rhs = rhs, requested_method = method, raw = NULL)
    class(obj) <- "mm_contrast"
    return(obj)
  }

  if (mm_boundary_df_method_unavailable(fit, method)) {
    table <- mm_boundary_df_unavailable_contrast_table(fit, L, rhs, method)
    obj <- list(
      table = table,
      L = L,
      rhs = rhs,
      requested_method = method,
      raw = NULL
    )
    class(obj) <- "mm_contrast"
    return(obj)
  }

  parsed <- mm_rust_contrast_table(fit, L, rhs, method, bootstrap = bootstrap)
  table <- parsed$table
  table$contrast <- table$label
  table$rhs <- rhs
  table$requested_method <- method
  table <- table[, c("contrast", "estimate", "rhs", "std_error", "df",
                     "statistic", "statistic_name", "p_value", "method",
                     "requested_method", "status", "reliability",
                     "reliability_reason", "estimability", "reason",
                     "reason_code", "reason_detail", "details", "notes"),
                 drop = FALSE]
  obj <- list(
    table = table,
    L = L,
    rhs = rhs,
    requested_method = method,
    raw = parsed$raw
  )
  class(obj) <- "mm_contrast"
  obj
}

#' Fixed-effect bootstrap control
#'
#' @param nsim Requested bootstrap replicate count.
#' @param seed Optional integer seed. `NULL` (the default) draws a seed from
#'   R's random number generator when the bootstrap runs, so `set.seed()`
#'   makes the bootstrap reproducible; the seed used is recorded in the
#'   result.
#' @param failed_refit_policy How failed refits are accounted for. Stable Rust
#'   wire labels are `"exclude"`, `"count_extreme"`, and `"abort"`.
#' @param threads Number of worker threads for the replicate refits
#'   (default `1`, serial). See the "Threads" section.
#'
#' @section Threads:
#' The engine simulates every replicate's response serially, in the serial
#' random-number order, and refits the replicates on `threads` worker
#' threads. The result is identical for every `threads` value; only the
#' wall-clock time changes. Workers never call into R: the R thread alone
#' polls for interrupts, so Ctrl-C / Esc still stops the run. The default is
#' `1` because CRAN policy limits packages to at most two threads unless the
#' user explicitly asks for more; passing `threads` is that request (for
#' example `threads = parallel::detectCores()`).
#'
#' @return A list used by `contrast(..., method = "bootstrap")`,
#'   `test_effect()`, and `confint(..., method = "bootstrap")`.
#'
#' @examples
#' bootstrap_control(nsim = 199, seed = 1, threads = 2)
#'
#' @export
bootstrap_control <- function(nsim = 999L,
                              seed = NULL,
                              failed_refit_policy = c("exclude", "count_extreme", "abort"),
                              threads = 1L) {
  failed_refit_policy <- match.arg(failed_refit_policy)
  threads <- mm_check_threads(threads)
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    mm_abort(
      message = "`nsim` must be a positive integer.",
      class = "mm_arg_error",
      input = nsim
    )
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0)) {
    mm_abort(
      message = "`seed` must be `NULL` or a non-negative integer.",
      class = "mm_arg_error",
      input = seed
    )
  }
  structure(
    list(
      requested_replicates = as.integer(nsim),
      seed = if (is.null(seed)) NULL else as.integer(seed),
      failed_refit_policy = failed_refit_policy,
      threads = threads
    ),
    class = "mm_bootstrap_control"
  )
}

# Validate a worker-thread count (bootstrap_control(), confint()/profile()
# `threads`). Any positive whole number is honoured as asked: the default is
# 1 and a larger value is the user's explicit opt-in (CRAN policy: <= 2
# unless the user asks for more).
mm_check_threads <- function(threads) {
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
      threads < 1 || threads != round(threads) || !is.finite(threads)) {
    mm_abort(
      message = "`threads` must be a single positive whole number.",
      class = "mm_arg_error",
      input = threads
    )
  }
  as.integer(threads)
}

# Wire form of a bootstrap control. A NULL seed is resolved from R's RNG at
# call time (as the GLMM bootstrap does) so `set.seed()` governs the engine's
# otherwise entropy-seeded StdRng.
mm_bootstrap_wire <- function(bootstrap) {
  out <- unclass(bootstrap)
  out$threads <- out$threads %||% 1L
  if (is.null(out$seed)) {
    out$seed <- sample.int(.Machine$integer.max, 1L)
  }
  out
}

#' @method print mm_contrast
#' @export
print.mm_contrast <- function(x, ...) {
  cat("Fixed-effect contrasts:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Test a fixed-effect term
#'
#' `test_effect()` asks Rust to construct fixed-effect term hypotheses and
#' returns the corresponding fixed-effect inference rows.
#'
#' @param fit A fitted `mm_lmm`.
#' @param term A fixed-effect term label.
#' @param method Requested inference method.
#' @param bootstrap Optional [bootstrap_control()] object for bootstrap-backed
#'   methods.
#' @param group Optional grouping factor for `method = "cluster_bootstrap"`.
#'   Required for crossed or multi-grouping-factor models. In schema 1.0.0,
#'   cluster resampling is an estimator-distribution target and term-level
#'   p-values return `not_assessed` with a stable reason code.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_effect_test` object.
#'
#' @export
test_effect <- function(fit, term, method = c("auto", "satterthwaite",
                                             "kenward_roger", "bootstrap",
                                             "bootstrap_lrt",
                                             "cluster_bootstrap",
                                             "asymptotic", "boundary_lrt",
                                             "none"),
                        bootstrap = NULL, group = NULL, ...) {
  UseMethod("test_effect")
}

#' @rdname test_effect
#' @export
test_effect.mm_lmm <- function(fit, term, method = c("auto", "satterthwaite",
                                                     "kenward_roger", "bootstrap",
                                                     "bootstrap_lrt",
                                                     "cluster_bootstrap",
                                                     "asymptotic", "boundary_lrt",
                                                     "none"),
                               bootstrap = NULL, group = NULL, ...) {
  method <- match.arg(method)
  if (!is.character(term) || !length(term)) {
    mm_abort(
      message = "`term` must be a non-empty character vector.",
      class = "mm_arg_error",
      input = term
    )
  }
  fixed_terms <- mm_fixed_effect_terms(fit)
  unknown <- setdiff(term, fixed_terms)
  if (length(unknown)) {
    mm_abort(
      message = sprintf("Unknown fixed-effect term(s): %s.",
                        paste(unknown, collapse = ", ")),
      class = "mm_arg_error",
      input = unknown
    )
  }
  if (identical(method, "boundary_lrt")) {
    table <- mm_boundary_lrt_fixed_effect_table(term)
  } else if (identical(method, "none")) {
    table <- mm_unavailable_effect_table(term, method)
  } else if (identical(method, "bootstrap")) {
    rows <- lapply(term, function(t) {
      mm_rust_term_bootstrap_row(fit, t, bootstrap)
    })
    table <- do.call(rbind, rows)
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reliability_reason", "reason",
                       "reason_code", "reason_detail", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
  } else if (identical(method, "bootstrap_lrt")) {
    rows <- lapply(term, function(t) {
      mm_rust_term_bootstrap_lrt_row(fit, t, bootstrap)
    })
    table <- do.call(rbind, rows)
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reliability_reason", "reason",
                       "reason_code", "reason_detail", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
  } else if (identical(method, "cluster_bootstrap")) {
    table <- mm_cluster_bootstrap_unavailable_effect_table(fit, term, group)
  } else if (mm_boundary_df_method_unavailable(fit, method)) {
    table <- mm_boundary_df_unavailable_effect_table(term, method)
  } else {
    parsed <- mm_rust_term_table(fit, method)
    table <- parsed$table[parsed$table$term %in% term, , drop = FALSE]
    table$requested_method <- method
    table <- table[, c("term", "numerator_df", "denominator_df", "statistic",
                       "statistic_name", "p_value", "method", "requested_method",
                       "status", "reliability", "reliability_reason", "reason",
                       "reason_code", "reason_detail", "details", "notes"),
                   drop = FALSE]
    names(table)[names(table) == "numerator_df"] <- "num_df"
    names(table)[names(table) == "denominator_df"] <- "den_df"
  }
  obj <- list(table = table, requested_method = method)
  class(obj) <- "mm_effect_test"
  obj
}

#' Test a random-effect variance component
#'
#' `test_random_effect()` exposes the boundary-aware likelihood-ratio route for
#' random-effect variance components. The v1 certified route is a nested ML
#' comparison that adds exactly one variance parameter and reports the
#' Self-Liang 50:50 mixture reference distribution; a comparison that adds
#' only correlation parameters (e.g. `||` versus `|`) is referred to the
#' ordinary chi-square, since a zero correlation is interior to the
#' parameter space (`method = "boundary_lrt_ordinary_chisq"`). It is intentionally
#' separate from [test_effect()], which tests fixed effects.
#'
#' @param fit A fitted `mm_lmm`.
#' @param term Random-effect term to test. This can be the term id (`"r0"`),
#'   the original random-effect fragment such as `"(1 | subject)"`, or a
#'   unique grouping factor name such as `"subject"`.
#' @param method Currently `"boundary_lrt"`.
#' @param refit_for_comparison How to handle REML fits. `"auto"` and `"ml"`
#'   refit to ML; `"error"` refuses.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_random_effect_test` object with a one-row `table`.
#'
#' @export
test_random_effect <- function(fit, term,
                               method = c("boundary_lrt"),
                               refit_for_comparison = c("auto", "error", "ml"),
                               ...) {
  UseMethod("test_random_effect")
}

#' @rdname test_random_effect
#' @export
test_random_effect.mm_lmm <- function(fit, term,
                                      method = c("boundary_lrt"),
                                      refit_for_comparison = c("auto", "error", "ml"),
                                      ...) {
  method <- match.arg(method)
  refit_for_comparison <- match.arg(refit_for_comparison)
  if (!is.character(term) || length(term) != 1L || is.na(term) || !nzchar(term)) {
    mm_abort(
      message = "`term` must be one random-effect term id, fragment, or group.",
      class = "mm_arg_error",
      input = term
    )
  }

  full <- mm_boundary_lrt_ml_fit(fit, refit_for_comparison)
  terms <- mm_random_effect_term_table(full)
  selected <- mm_match_random_effect_term(terms, term)
  reduced_formula <- mm_drop_random_term_formula(full, selected$index)
  remaining <- nrow(terms) - 1L
  reduced <- NULL
  reduced_payload <- NULL
  if (remaining > 0L) {
    reduced <- mm_internal_lmm(reduced_formula, full$model_frame, REML = FALSE,
                               weights = full$weights,
                               control = mm_internal_control(full, keep_start = FALSE))
    reduced_payload <- mm_boundary_lrt_bridge_payload(reduced)
  }

  json <- tryCatch(
    mm_boundary_lrt_json(
      reduced_payload,
      mm_boundary_lrt_bridge_payload(full),
      deparse1(reduced_formula)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = method, term = term)
  }
  payload <- mm_json_parse_boundary_lrt(json)
  table <- mm_boundary_lrt_table(
    payload = payload,
    term = selected,
    reduced_formula = reduced_formula,
    full = full,
    reduced = reduced,
    refit = !identical(full, fit)
  )
  obj <- list(
    table = table,
    requested_method = method,
    term = selected,
    reduced = reduced,
    full = full,
    raw = payload
  )
  class(obj) <- "mm_random_effect_test"
  obj
}

#' @method print mm_random_effect_test
#' @export
print.mm_random_effect_test <- function(x, ...) {
  cat("Random-effect variance-component test:\n")
  cols <- intersect(
    c("term", "group", "statistic", "statistic_name", "p_value",
      "reference_distribution", "status", "reason_code"),
    names(x$table)
  )
  print(x$table[, cols, drop = FALSE], row.names = FALSE)
  invisible(x)
}

mm_boundary_lrt_ml_fit <- function(fit, refit_for_comparison) {
  if (!inherits(fit, "mm_lmm")) {
    mm_abort(
      message = "`test_random_effect()` requires a fitted `mm_lmm` object.",
      class = "mm_arg_error",
      input = fit
    )
  }
  if (!isTRUE(fit$REML)) {
    return(fit)
  }
  if (identical(refit_for_comparison, "error")) {
    mm_abort(
      message = "Boundary LRT requires ML-fitted models; use `refit_for_comparison = \"auto\"` or refit with `REML = FALSE`.",
      class = "mm_inference_unavailable",
      reason_code = "boundary_lrt_requires_ml",
      input = list(REML = TRUE)
    )
  }
  mm_internal_lmm(fit$formula, fit$model_frame, REML = FALSE,
                  weights = fit$weights, offset = mm_fit_offset_arg(fit),
                  control = mm_internal_control(fit))
}

mm_boundary_lrt_bridge_payload <- function(fit) {
  payload <- mm_rust_fit_bridge_payload(fit)
  payload$REML <- isTRUE(fit$REML)
  payload
}

mm_random_effect_term_table <- function(fit) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  if (!length(terms)) {
    mm_abort(
      message = "No random-effect terms are available to test.",
      class = "mm_inference_unavailable",
      reason_code = "boundary_lrt_requires_variance_component_comparison"
    )
  }
  rows <- lapply(seq_along(terms), function(i) {
    term <- terms[[i]]
    basis <- term$basis %||% list()
    p <- length(basis)
    if (!p) p <- 1L
    covariance <- mm_scalar_text(term$covariance, "full")
    theta_parameters <- switch(
      covariance,
      scalar = 1L,
      diagonal = p,
      diag = p,
      p * (p + 1L) / 2L
    )
    data.frame(
      index = i,
      term_id = term$id %||% sprintf("r%d", i - 1L),
      term = term$source_syntax$text %||% sprintf("random term %d", i),
      group = mm_random_term_group_label(fit, term, i),
      basis = paste(vapply(basis, mm_basis_label, character(1)), collapse = ", "),
      covariance = covariance,
      theta_parameters = as.integer(theta_parameters),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

mm_match_random_effect_term <- function(terms, term) {
  hits <- terms$term_id == term | terms$term == term | terms$group == term
  if (!any(hits)) {
    mm_abort(
      message = sprintf(
        "Unknown random-effect term `%s`. Available terms: %s.",
        term,
        paste(unique(c(terms$term_id, terms$term, terms$group)), collapse = ", ")
      ),
      class = "mm_arg_error",
      input = term
    )
  }
  if (sum(hits) > 1L) {
    mm_abort(
      message = sprintf(
        "Random-effect term `%s` is ambiguous; use the term id or exact fragment.",
        term
      ),
      class = "mm_arg_error",
      input = term
    )
  }
  terms[which(hits)[[1L]], , drop = FALSE]
}

mm_drop_random_term_formula <- function(fit, drop_index) {
  response <- mm_response_name(fit)
  fixed <- mm_fixed_effect_terms(fit)
  fixed_rhs <- mm_fixed_rhs_text(fixed, mm_fixed_terms_have_intercept(fixed))
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  random <- vapply(terms, function(x) x$source_syntax$text %||% "", character(1))
  random <- random[nzchar(random)]
  if (drop_index <= length(random)) {
    random <- random[-drop_index]
  }
  rhs <- paste(c(fixed_rhs, random), collapse = " + ")
  stats::as.formula(paste(response, "~", rhs), env = environment(fit$formula))
}

mm_json_parse_boundary_lrt <- function(json) {
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
        message = sprintf("Failed to parse boundary-LRT JSON: %s", conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  if (!identical(as.character(parsed$schema_name), "mixedmodels.boundary_lrt") ||
      !identical(as.character(parsed$schema_version), "1.0.0")) {
    mm_abort(
      message = "Boundary-LRT JSON has an unknown schema header.",
      class = "mm_schema_error",
      input = parsed
    )
  }
  parsed
}

mm_boundary_lrt_table <- function(payload, term, reduced_formula, full, reduced, refit) {
  status <- as.character(payload$status %||% "not_assessed")
  available <- identical(status, "available")
  reference_distribution <- mm_boundary_lrt_reference_label(payload$mixture)
  # The engine uses the 50:50 boundary mixture only when exactly one
  # variance is added; a comparison adding only correlations (`||` vs `|`)
  # is referred to the ordinary chi-square (a single weight-1 component).
  mixture <- payload$mixture %||% list()
  is_mixture <- any(vapply(mixture, function(component) {
    !is.null(component$point_mass_at) ||
      !isTRUE(all.equal(as.numeric(component$weight %||% NA_real_), 1))
  }, logical(1)))
  reason <- payload$reason %||% if (available) NA_character_ else
    "boundary LRT did not certify a p-value"
  reason_code <- payload$reason_code %||% if (available) NA_character_ else
    "boundary_lrt_not_assessed"
  details <- list(
    schema_name = payload$schema_name,
    schema_version = payload$schema_version,
    mixture = payload$mixture %||% list(),
    references = payload$references %||% character(),
    notes = payload$notes %||% character(),
    comparison_class = payload$comparison_class %||% NA_character_,
    ordinary_chisq_dof = payload$ordinary_chisq_dof %||% NA_integer_,
    loglik_within_optimizer_tol = payload$loglik_within_optimizer_tol %||% NA,
    reduced_formula = deparse1(reduced_formula),
    full_formula = deparse1(full$formula),
    reduced_class = if (is.null(reduced)) "lm_fixed_effect_submodel" else "mm_lmm",
    shape_restrictions = c(
      "nested ML fits",
      "same fixed-effect column space",
      "exactly one tested variance/covariance parameter",
      "no extra nuisance parameter on the boundary"
    )
  )
  data.frame(
    term = term$term,
    term_id = term$term_id,
    group = term$group,
    theta_parameters = term$theta_parameters,
    statistic = payload$statistic %||% NA_real_,
    statistic_name = if (is_mixture || !available) "chi_bar_square" else "chi_square",
    ordinary_chisq_dof = payload$ordinary_chisq_dof %||% NA_integer_,
    p_value = payload$pvalue %||% NA_real_,
    method = if (!available) {
      "boundary_lrt"
    } else if (is_mixture) {
      "boundary_lrt_self_liang_mixture"
    } else {
      "boundary_lrt_ordinary_chisq"
    },
    requested_method = "boundary_lrt",
    status = status,
    reason = reason,
    reason_code = reason_code,
    reference_distribution = reference_distribution,
    refit = isTRUE(refit),
    details = I(list(details)),
    notes = I(list(payload$notes %||% character())),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mm_boundary_lrt_reference_label <- function(mixture) {
  if (is.null(mixture) || !length(mixture)) {
    return(NA_character_)
  }
  parts <- vapply(mixture, function(component) {
    weight <- as.numeric(component$weight %||% NA_real_)
    df <- component$chisq_df
    if (is.null(df) && identical(as.numeric(component$point_mass_at %||% NA_real_), 0)) {
      df <- 0L
    }
    if (is.null(df)) {
      sprintf("%.3g * point-mass(%s)", weight,
              as.character(component$point_mass_at %||% NA_real_))
    } else {
      sprintf("%.3g * chi-square(%s)", weight, as.integer(df))
    }
  }, character(1))
  paste(parts, collapse = " + ")
}

mm_boundary_lrt_fixed_effect_reason <- function() {
  "boundary_lrt is a variance-component route, not a fixed-effect p-value method"
}

mm_boundary_lrt_fixed_effect_contrast_table <- function(L, rhs) {
  reason <- mm_boundary_lrt_fixed_effect_reason()
  reason_code <- "boundary_lrt_not_applicable_to_fixed_effects"
  data.frame(
    contrast = rownames(L),
    estimate = NA_real_,
    rhs = rhs,
    std_error = NA_real_,
    df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = "not_applicable",
    requested_method = "boundary_lrt",
    status = "unsupported",
    reliability = "not_available",
    reliability_reason = "not_applicable_to_fixed_effects",
    estimability = "not_assessed",
    reason = reason,
    reason_code = reason_code,
    reason_detail = reason,
    details = I(rep(list(list(boundary_lrt = list(
      route = "variance_component",
      fixed_effect_applicable = FALSE
    ))), nrow(L))),
    notes = I(rep(list(character()), nrow(L))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mm_boundary_lrt_fixed_effect_table <- function(term) {
  reason <- mm_boundary_lrt_fixed_effect_reason()
  reason_code <- "boundary_lrt_not_applicable_to_fixed_effects"
  data.frame(
    term = term,
    num_df = NA_real_,
    den_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = "not_applicable",
    requested_method = "boundary_lrt",
    status = "unsupported",
    reliability = "not_available",
    reliability_reason = "not_applicable_to_fixed_effects",
    reason = reason,
    reason_code = reason_code,
    reason_detail = reason,
    details = I(rep(list(list(boundary_lrt = list(
      route = "variance_component",
      fixed_effect_applicable = FALSE
    ))), length(term))),
    notes = I(rep(list(character()), length(term))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' @method print mm_effect_test
#' @export
print.mm_effect_test <- function(x, ...) {
  cat("Effect tests:\n")
  cols <- intersect(
    c("term", "num_df", "den_df", "statistic", "statistic_name", "p_value",
      "method", "status", "reason_code"),
    names(x$table)
  )
  show <- x$table[, cols, drop = FALSE]
  show <- show[, vapply(show, function(col) {
    if (is.numeric(col)) return(any(!is.na(col)))
    values <- as.character(col)
    any(!is.na(values) & nzchar(values))
  }, logical(1)), drop = FALSE]
  print(show, row.names = FALSE)
  hidden <- setdiff(names(x$table), names(show))
  if (length(hidden)) {
    cat(sprintf("Full audit columns available in `x$table` (%d hidden).\n",
                length(hidden)))
  }
  invisible(x)
}

#' Assess contrast estimability
#'
#' Routes each requested contrast row through the Rust fixed-effect inference
#' bridge and reports the upstream estimability assessment verbatim.
#' Returned rows carry `status` (the closed enum from upstream:
#' `estimable`, `not_estimable`, `aliased`, ...), a boolean `estimable`
#' convenience flag, the contrast `rank` and `requested_rank`, and a
#' stable `reason` populated only when the engine refuses the contrast.
#'
#' @param fit A fitted `mm_lmm`.
#' @param L Optional contrast vector or matrix. Defaults to the fixed-effect
#'   coefficient basis.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_estimability` object.
#'
#' @export
estimability <- function(fit, L = NULL, ...) {
  UseMethod("estimability")
}

#' @rdname estimability
#' @export
estimability.mm_lmm <- function(fit, L = NULL, ...) {
  if (is.null(L)) {
    L <- diag(length(fit$beta))
    colnames(L) <- names(fit$beta)
    rownames(L) <- names(fit$beta)
  } else {
    L <- mm_contrast_matrix(L, fit)
  }

  parsed <- tryCatch(
    mm_rust_contrast_table(fit, L, rep(0, nrow(L)), "auto", bootstrap = NULL),
    error = function(cnd) cnd
  )

  if (inherits(parsed, "condition")) {
    table <- data.frame(
      contrast = rownames(L),
      estimable = NA,
      status = "not_assessed",
      rank = NA_integer_,
      requested_rank = NA_integer_,
      reason = conditionMessage(parsed),
      stringsAsFactors = FALSE
    )
  } else {
    table <- mm_estimability_table_from_inference(parsed$table, rownames(L))
  }

  obj <- list(table = table, L = L)
  class(obj) <- "mm_estimability"
  obj
}

mm_estimability_table_from_inference <- function(inference_table, contrast_names) {
  rows <- lapply(seq_len(nrow(inference_table)), function(i) {
    payload <- inference_table$estimability[[i]]
    assessment <- payload$assessment %||% list()
    status <- as.character(assessment$status %||% "not_assessed")
    diagnostics <- assessment$diagnostics %||% list()
    reason <- if (identical(status, "estimable")) {
      NA_character_
    } else if (length(diagnostics)) {
      paste(vapply(diagnostics, function(d) as.character(d$message %||% d), character(1)),
            collapse = "; ")
    } else if (identical(status, "not_assessed")) {
      "rust_estimability_certificate_unavailable"
    } else {
      sprintf("contrast %s under the fitted fixed-effect design", status)
    }
    data.frame(
      contrast = contrast_names[[i]],
      estimable = if (identical(status, "estimable")) TRUE else
                  if (identical(status, "not_assessed")) NA else FALSE,
      status = status,
      rank = as.integer(assessment$rank %||% NA_integer_),
      requested_rank = as.integer(assessment$requested_rank %||% NA_integer_),
      reason = reason,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @method print mm_estimability
#' @export
print.mm_estimability <- function(x, ...) {
  cat("Estimability:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Degrees of freedom for a contrast
#'
#' @param fit A fitted `mm_lmm`.
#' @param L A contrast vector or matrix.
#' @param method Requested degrees-of-freedom method.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_df_for_contrast` object with `$table` (one row per
#'   contrast: `contrast`, `df`, `method`, `requested_method`, `reason`),
#'   `$df` (the named numeric vector), and `$method`. When the method is
#'   `"none"` or the engine refuses, `df` is `NA` and `reason` records why.
#'
#' @export
df_for_contrast <- function(fit, L, method = c("auto", "satterthwaite",
                                              "kenward_roger", "bootstrap",
                                              "asymptotic", "none"), ...) {
  UseMethod("df_for_contrast")
}

#' @rdname df_for_contrast
#' @export
df_for_contrast.mm_lmm <- function(fit, L, method = c("auto", "satterthwaite",
                                                      "kenward_roger", "bootstrap",
                                                      "asymptotic", "none"), ...) {
  method <- match.arg(method)
  L <- mm_contrast_matrix(L, fit)
  n <- nrow(L)
  out <- rep(NA_real_, n)
  names(out) <- rownames(L)

  if (identical(method, "none")) {
    return(mm_new_df_for_contrast(out, method_label = "not_requested",
                                  requested = method,
                                  reason = mm_inference_unavailable_reason(method)))
  }

  parsed <- tryCatch(
    mm_rust_contrast_table(fit, L, rep(0, n), method, bootstrap = NULL),
    error = function(cnd) cnd
  )
  if (inherits(parsed, "condition")) {
    return(mm_new_df_for_contrast(out, method_label = "unavailable",
                                  requested = method,
                                  reason = conditionMessage(parsed)))
  }

  table <- parsed$table
  df_values <- if ("df" %in% names(table)) {
    as.numeric(table$df)
  } else if ("denominator_df" %in% names(table)) {
    as.numeric(table$denominator_df)
  } else {
    rep(NA_real_, nrow(table))
  }
  out[] <- df_values[seq_along(out)]

  resolved_methods <- unique(stats::na.omit(as.character(table$method)))
  method_label <- if (length(resolved_methods) == 1L) resolved_methods else method
  reasons <- stats::na.omit(as.character(table$reason))
  reason_label <- if (all(is.na(out))) {
    if (length(reasons)) reasons[[1L]] else mm_inference_unavailable_reason(method)
  } else {
    NA_character_
  }

  mm_new_df_for_contrast(out, method_label = method_label,
                         requested = method, reason = reason_label)
}

# Sibling-consistent constructor: an mm_* object carrying $table plus the
# named df vector for programmatic use.
mm_new_df_for_contrast <- function(df, method_label, requested, reason) {
  if (is.null(names(df))) names(df) <- paste0("c", seq_along(df))
  table <- data.frame(
    contrast = names(df),
    df = as.numeric(df),
    method = method_label,
    requested_method = requested,
    reason = if (length(reason) == 1L && !is.na(reason)) reason else NA_character_,
    stringsAsFactors = FALSE
  )
  obj <- list(
    table = table,
    df = df,
    method = method_label,
    requested_method = requested
  )
  class(obj) <- "mm_df_for_contrast"
  obj
}

#' @method print mm_df_for_contrast
#' @export
print.mm_df_for_contrast <- function(x, ...) {
  cat("Degrees of freedom for contrasts:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

# lme4's confint.merMod spellings ("Wald", "boot") map onto mixeff's names.
# Only exact lme4 spellings are translated; the default vector passes
# through untouched so match.arg() picks the first entry.
mm_confint_method_spelling <- function(method) {
  if (length(method) != 1L) return(method)
  switch(method, Wald = "wald", boot = "bootstrap", method)
}

#' Confidence intervals for a mixeff LMM
#'
#' `method = "asymptotic"` (the default; synonyms `"wald"` and lme4's
#' `"Wald"`) gives Wald intervals for the fixed effects from the stored
#' standard errors. `method = "profile"` gives profile-likelihood intervals
#' with lme4's rows and names: `.sig01`, `.sig02`, ... (random-effect
#' standard deviations and correlations, numbered as lme4 >= 2.0 does: term
#' by term, each term's standard deviations before its correlations),
#' `.sigma`, then the fixed effects. As in lme4, every fit is profiled on the
#' ML deviance: a REML fit's intervals are those of its ML refit. The
#' engine's relative-Cholesky `theta1`, ... intervals (on the fit's own
#' criterion) are kept in `attr(ci, "mm_profile")$table` and can be
#' requested by name through `parm`. When one parameter's profile is
#' irregular (non-monotone, not bracketing the estimate, or a failed refit),
#' its undetermined bounds are `NA` -- lme4 would warn and interpolate
#' linearly, often reporting `[-1, 1]` for a correlation -- and the reason is
#' in `attr(ci, "mm_profile")$table` (`status`, `reason_code`, `reason`);
#' printing the result adds a note. See [profile.mm_lmm()].
#' `method = "bootstrap"` (lme4's
#' `"boot"`) gives parametric-bootstrap intervals for the fixed effects,
#' controlled by [bootstrap_control()].
#'
#' lme4's `confint.merMod()` defaults to `method = "profile"`; mixeff keeps
#' the Wald default because the profile refits the model many times
#' (about the cost of lme4's profile; several seconds for a modest model)
#' and can refuse on near-singular fits, while the Wald route is free.
#' Request `method = "profile"` to reproduce lme4's default output.
#'
#' @param object A fitted `mm_lmm`.
#' @param parm Parameters to report: names (fixed-effect names, and for
#'   `"profile"` also `".sig01"`, ..., `".sigma"`, `"theta1"`, ...), lme4's
#'   `"theta_"` / `"beta_"` shortcuts (`"profile"` only), or indices.
#'   Defaults to all fixed effects (all lme4 rows for `"profile"`).
#' @param level Confidence level.
#' @param method `"asymptotic"`, `"wald"`, `"bootstrap"`, or `"profile"`.
#' @param bootstrap Optional [bootstrap_control()] for `method = "bootstrap"`
#'   (its `threads` field sets the bootstrap's worker threads).
#' @param interval Bootstrap interval type, `"percentile"` or `"basic"`.
#' @param threads Worker threads for `method = "profile"` (default 1). The
#'   per-parameter profiles run concurrently; the result is identical for
#'   every value. The default is 1 because CRAN policy limits packages to at
#'   most two threads unless the user explicitly asks for more; passing
#'   `threads` is that request.
#' @param ... Unused.
#'
#' @return An `mm_confint` matrix of lower/upper bounds.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(60))
#' d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(60)
#' fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
#' confint(fit)
#' confint(fit, method = "profile")
#'
#' @method confint mm_lmm
#' @importFrom stats confint
#' @export
confint.mm_lmm <- function(object, parm, level = 0.95,
                           method = c("asymptotic", "wald", "bootstrap",
                                      "profile"),
                           bootstrap = NULL,
                           interval = c("percentile", "basic"),
                           threads = 1L, ...) {
  method <- match.arg(mm_confint_method_spelling(method),
                      c("asymptotic", "wald", "bootstrap", "profile"))
  # `"asymptotic"` is the package-wide canonical name for the closed-form
  # Wald interval (shared with contrast()/test_effect()/inference_table());
  # `"wald"` stays accepted as a synonym. Internals compute under "wald".
  if (identical(method, "asymptotic")) method <- "wald"
  interval <- match.arg(interval)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    mm_abort(
      message = "`level` must be a single number between 0 and 1.",
      class = "mm_arg_error",
      input = level
    )
  }
  if (identical(method, "profile")) {
    return(mm_profile_confint(object, if (missing(parm)) NULL else parm,
                              level, threads = threads))
  }
  terms <- names(object$beta)
  if (missing(parm)) {
    parm <- terms
  }
  if (is.numeric(parm)) {
    parm <- terms[parm]
  }
  unknown <- setdiff(parm, terms)
  if (length(unknown)) {
    mm_abort(
      message = sprintf("Unknown fixed-effect parameter(s): %s.",
                        paste(unknown, collapse = ", ")),
      class = "mm_arg_error",
      input = unknown
    )
  }
  if (identical(method, "bootstrap")) {
    return(mm_bootstrap_confint(object, parm, level, bootstrap, interval))
  }
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)
  est <- object$beta[parm]
  se <- object$std_errors[parm]
  out <- cbind(est - crit * se, est + crit * se)
  colnames(out) <- c(sprintf("%.1f %%", 100 * alpha / 2),
                     sprintf("%.1f %%", 100 * (1 - alpha / 2)))
  attr(out, "method") <- "wald_asymptotic_from_stored_standard_errors"
  attr(out, "status") <- "not_certified_by_rust_inference_contract"
  mm_new_confint(out)
}

#' Confidence intervals for fixed effects of a mixeff GLMM
#'
#' Two routes are available. `method = "asymptotic"` builds Wald intervals
#' (`estimate +/- z * SE`) from the Rust fixed-effect inference table; it
#' requires a certified fit (`method = "joint_laplace"`) or the explicit
#' working-Hessian opt-in, and refuses the profiled (`pirls_profiled`) estimator with
#' a typed reason. `method = "bootstrap"` runs a parametric bootstrap —
#' each replicate simulates a response from the fitted model under fresh
#' random-effect draws and refits the same (effective) estimator — and
#' returns percentile intervals; it works on profiled-estimator fits for
#' every supported family, including negative binomial (fixed-theta fits
#' condition replicates on that theta; estimated-theta fits re-estimate
#' theta per replicate, propagating its uncertainty). On `joint_laplace`
#' fits the bootstrap is refused at this engine pin — replicate refits
#' cannot re-run the joint estimator — and those fits have certified Wald
#' intervals instead.
#' The bootstrap result carries replicate accounting as attributes:
#' `mm_bootstrap` (requested, successful, failed, seed, per-parameter
#' standard errors, Monte Carlo SE, and a reliability grade) and
#' `mm_estimator`. Failed replicates are counted, never silently dropped.
#' Profile intervals are not certified for GLMMs and are refused.
#'
#' @param object A fitted `mm_glmm`.
#' @param parm Optional fixed-effect names or indices; defaults to all.
#' @param level Confidence level.
#' @param method `"asymptotic"` (the default; the package-wide name for the
#'   closed-form Wald interval), its synonyms `"wald"` and lme4's `"Wald"`, or
#'   `"bootstrap"` (lme4's spelling `"boot"` is accepted). lme4 defaults to
#'   `"profile"`; mixeff keeps the Wald default because profile intervals are
#'   refused for GLMMs. The same spellings are accepted by the LMM method
#'   ([confint.mm_lmm()]), whose `"profile"` intervals use lme4's `.sig01`,
#'   ..., `.sigma` rows.
#' @param ... For `method = "bootstrap"`: `nsim` (number of replicates,
#'   default 999), `seed` (non-negative integer; when omitted one is
#'   drawn from R's RNG so `set.seed()` governs reproducibility), and
#'   `threads` (worker threads for the replicate refits, default 1; the
#'   result is identical for every value -- see [bootstrap_control()] for
#'   the threading contract and CRAN's two-thread default policy).
#'
#' @return An `mm_confint` matrix of lower/upper bounds (`"asymptotic"`), or
#'   a plain matrix of percentile bounds with bootstrap-accounting
#'   attributes (`"bootstrap"`).
#'
#' @method confint mm_glmm
#' @export
confint.mm_glmm <- function(object, parm, level = 0.95,
                            method = c("asymptotic", "wald", "profile",
                                       "bootstrap"), ...) {
  method <- match.arg(mm_confint_method_spelling(method),
                      c("asymptotic", "wald", "profile", "bootstrap"))
  if (identical(method, "profile")) {
    # Advice must be actionable on THIS fit: asymptotic Wald intervals are
    # certified only for joint_laplace fits.
    advice <- if (identical(object$method, "joint_laplace")) {
      # The bootstrap cannot refit a joint template at this engine pin, so
      # it must not be advised here.
      "Use method = \"asymptotic\"."
    } else {
      paste0("Use method = \"bootstrap\", or re-fit with ",
             "glmm(..., method = \"joint_laplace\") and use ",
             "method = \"asymptotic\".")
    }
    mm_abort(
      message = sprintf(
        paste0("`confint(method = \"profile\")` is not available for GLMM ",
               "fits; profile likelihood is not certified for generalized ",
               "models. %s"),
        advice
      ),
      class = "mm_inference_unavailable",
      reason_code = "glmm_confint_method_unavailable",
      input = method
    )
  }
  if (identical(method, "asymptotic")) method <- "wald"
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    mm_abort(
      message = "`level` must be a single number between 0 and 1.",
      class = "mm_arg_error",
      input = level
    )
  }
  terms <- names(object$beta)
  if (missing(parm)) {
    parm <- terms
  }
  if (is.numeric(parm)) {
    parm <- terms[parm]
  }
  unknown <- setdiff(parm, terms)
  if (length(unknown)) {
    mm_abort(
      message = sprintf("Unknown fixed-effect parameter(s): %s.",
                        paste(unknown, collapse = ", ")),
      class = "mm_arg_error",
      input = unknown
    )
  }
  if (identical(method, "bootstrap")) {
    return(mm_glmm_parametric_bootstrap_confint(object, parm = parm,
                                                level = level, ...))
  }
  inference <- mm_glmm_wald_z_inference(object)
  rows <- inference$table
  rows <- rows[match(parm, rows$label), , drop = FALSE]
  # `available` = certified rows; `available_noninferential` = the labelled
  # working-Hessian rows produced only under the explicit fit-level opt-in
  # (the builder never emits them otherwise; WI-2.9).
  bad <- is.na(rows$label) |
    is.na(rows$status) |
    !(rows$status %in% c("available", "available_noninferential")) |
    !is.finite(rows$std_error) |
    rows$std_error <= 0
  if (any(bad)) {
    mm_abort(
      message = mm_glmm_wald_refusal_message("confint() Wald intervals"),
      class = "mm_inference_unavailable",
      reason_code = "glmm_wald_confint_unavailable",
      input = parm[bad]
    )
  }
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)
  est <- stats::setNames(rows$estimate, rows$label)
  se <- stats::setNames(rows$std_error, rows$label)
  out <- cbind(est - crit * se, est + crit * se)
  colnames(out) <- c(sprintf("%.1f %%", 100 * alpha / 2),
                     sprintf("%.1f %%", 100 * (1 - alpha / 2)))
  # The attrs must carry the rows' own certificate, not a blanket one: for
  # an inference = "working_hessian" fit the intervals come from the
  # labelled uncertified approximation, and stamping them
  # available/from-the-inference-table would be a false certificate
  # (Audit1 WI-2.9).
  if (any(rows$method == "wald_z_working_hessian", na.rm = TRUE)) {
    attr(out, "method") <- "wald_z_working_hessian"
    attr(out, "status") <- "available_noninferential"
  } else {
    attr(out, "method") <- "wald_asymptotic_from_rust_inference_table"
    attr(out, "status") <- "available"
  }
  mm_new_confint(out)
}

#' @method print mm_confint
#' @export
print.mm_confint <- function(x, ...) {
  cat("Confidence intervals:\n")
  show <- matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )
  print(show, ...)

  method <- attr(x, "method") %||% "not_recorded"
  status <- attr(x, "status") %||% "not_recorded"
  cat(sprintf("method: %s\n", method))
  interval <- attr(x, "interval")
  if (!is.null(interval)) {
    cat(sprintf("interval: %s\n", interval))
  }
  # Display translation for the internal status label (the stored attribute
  # is unchanged): the raw contract label reads as an alarm without telling
  # the user anything actionable.
  status_shown <- if (identical(status, "not_certified_by_rust_inference_contract")) {
    paste0("Wald (asymptotic) intervals from stored standard errors ",
           "(engine-certified profile intervals: method = \"profile\")")
  } else {
    status
  }
  cat(sprintf("status: %s\n", status_shown))
  mm_print_profile_irregular_note(x)

  bootstrap <- attr(x, "bootstrap")
  if (length(bootstrap)) {
    summary <- mm_bootstrap_confint_summary(x, bootstrap)
    cat("\nBootstrap run:\n")
    print(summary, row.names = FALSE)
    notes <- unique(unlist(lapply(bootstrap, function(payload) {
      payload$metadata$notes %||% character()
    }), use.names = FALSE))
    if (length(notes)) {
      cat("notes:\n")
      for (note in notes) {
        wrapped <- strwrap(note, width = 78, exdent = 4)
        cat(sprintf("  - %s\n", paste(wrapped, collapse = "\n    ")))
      }
    }
    cat("Full bootstrap payload available in `attr(x, \"bootstrap\")`.\n")
  }
  invisible(x)
}

mm_contrast_matrix <- function(L, fit) {
  if (is.null(L)) {
    mm_abort(
      message = "`L` must be a numeric contrast vector or matrix.",
      class = "mm_arg_error",
      input = L
    )
  }
  if (is.vector(L)) {
    L <- matrix(as.numeric(L), nrow = 1L)
  }
  if (!is.matrix(L) || !is.numeric(L)) {
    mm_abort(
      message = "`L` must be a numeric contrast vector or matrix.",
      class = "mm_arg_error",
      input = L
    )
  }
  if (ncol(L) != length(fit$beta)) {
    mm_abort(
      message = sprintf("`L` must have %d column(s), one for each fixed effect.",
                        length(fit$beta)),
      class = "mm_arg_error",
      input = L
    )
  }
  colnames(L) <- colnames(L) %||% names(fit$beta)
  if (is.null(rownames(L))) {
    rownames(L) <- paste0("c", seq_len(nrow(L)))
  }
  L
}


# One-line scale notice before a long silent bootstrap FFI call: the engine
# runs every replicate refit inside a single native call with no progress
# output, so say up front roughly how much work was requested. Suppressed by
# mm_control(verbose = -1).
mm_inform_bootstrap_scale <- function(fit, bootstrap) {
  verbose <- fit$control$verbose %||% 0L
  nsim <- bootstrap$nsim %||% NA_integer_
  if (verbose >= 0L && is.finite(nsim) && nsim >= 200) {
    mm_inform(
      sprintf(paste0(
        "Running %d bootstrap refits in a single native call with no ",
        "progress output; expect roughly %d x the cost of one fit. Silence ",
        "this message with mm_control(verbose = -1)."
      ), as.integer(nsim), as.integer(nsim)),
      class = "mm_runtime_notice"
    )
  }
  invisible(NULL)
}

mm_rust_contrast_table <- function(fit, L, rhs, method, bootstrap = NULL) {
  bridge <- mm_rust_fit_bridge_payload(fit)
  if (identical(method, "bootstrap") && !is.null(bootstrap)) {
    if (!inherits(bootstrap, "mm_bootstrap_control")) {
      bootstrap <- do.call(bootstrap_control, as.list(bootstrap))
    }
    mm_inform_bootstrap_scale(fit, bootstrap)
    bootstrap_json <- jsonlite::toJSON(
      mm_bootstrap_wire(bootstrap),
      auto_unbox = TRUE,
      null = "null"
    )
    json <- tryCatch(
      mm_fixed_effect_bootstrap_contrast_json(
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
        as.character(rownames(L)),
        as.numeric(rhs),
        as.character(bootstrap_json)
      ),
      error = function(cnd) cnd
    )
  } else {
    json <- tryCatch(
      mm_fixed_effect_contrast_json(
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
        as.character(rownames(L)),
        as.numeric(rhs),
        method
      ),
      error = function(cnd) cnd
    )
  }
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = method)
  }
  table <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse fixed-effect contrast JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  mm_json_parse_fixed_effect_inference_table(table)
}

mm_rust_term_table <- function(fit, method, type = "III") {
  bridge <- mm_rust_fit_bridge_payload(fit)
  json <- tryCatch(
    mm_fixed_effect_term_json(
      bridge$formula_string,
      isTRUE(fit$REML),
      bridge$spec_data$column_order,
      bridge$spec_data$numeric_columns,
      bridge$spec_data$categorical_values,
      bridge$spec_data$categorical_levels,
      bridge$spec_data$categorical_ordered,
      bridge$weights,
      bridge$control_json,
      method,
      mm_fixed_effect_term_type_label(type)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = method)
  }
  table <- tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse fixed-effect term JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  mm_json_parse_fixed_effect_inference_table(table)
}

mm_fixed_effect_term_type_label <- function(type) {
  type <- as.character(type)
  if (length(type) != 1L || is.na(type) || !nzchar(type)) {
    mm_abort(
      message = "`type` must be one of `I`, `II`, `III`, or `block`.",
      class = "mm_arg_error",
      input = type
    )
  }
  switch(
    toupper(type),
    I = "type_i",
    II = "type_ii",
    III = "type_iii",
    # Simple effects at reference levels (the pre-1f3f689 "Type III"),
    # kept available under its honest name (engine coefficient_block).
    BLOCK = "coefficient_block",
    mm_abort(
      message = "`type` must be one of `I`, `II`, `III`, or `block`.",
      class = "mm_arg_error",
      input = type
    )
  )
}

mm_rust_fit_bridge_payload <- function(fit, warm_start = FALSE) {
  spec_data <- mm_translate_data(mm_engine_frame(fit))
  formula_string <- mm_coerce_formula_string(mm_engine_formula(fit))
  list(
    spec_data = spec_data,
    formula_string = formula_string,
    weights = mm_bridge_weights(fit$weights),
    control_json = mm_refit_control_json(fit, warm_start = warm_start)
  )
}

# Control JSON for an engine-side refit of `fit` itself. The user's
# mm_control() is kept. Refits that recompute quantities of the stored fit
# (contrasts, predictions, profiles) deliberately do NOT warm start: a cold
# start with the same control replays the original optimizer path and lands
# on the identical optimum, whereas a warm start from the fitted theta stops
# at a point that differs at the optimizer tolerance (~1e-6 relative), so
# summary()/predict(newdata) would no longer reproduce fixef()/fitted(); it
# also saved < 10% of the refit time on the benchmarks. `warm_start = TRUE`
# (start from the fitted theta unless the user set `start`) is for resampling
# refits of the same model.
mm_refit_control_json <- function(fit, warm_start = FALSE) {
  control <- unclass(fit$control %||% mm_control())
  if (isTRUE(warm_start) && inherits(fit, "mm_lmm") && is.null(control$start)) {
    theta <- as.numeric(unlist(fit$theta %||% numeric(), use.names = FALSE))
    if (length(theta) && all(is.finite(theta))) {
      control$start <- theta
    }
  }
  as.character(jsonlite::toJSON(control, auto_unbox = TRUE, null = "null",
                                digits = NA))
}

mm_unavailable_effect_table <- function(term, method) {
  data.frame(
    term = term,
    num_df = NA_real_,
    den_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = "not_computed",
    requested_method = method,
    status = "not_assessed",
    reliability = "not_available",
    reason = mm_inference_unavailable_reason(method),
    details = I(rep(list(NULL), length(term))),
    notes = I(rep(list(character()), length(term))),
    stringsAsFactors = FALSE
  )
}

mm_boundary_df_method_unavailable <- function(fit, method) {
  method %in% c("satterthwaite", "kenward_roger") && isTRUE(is_singular(fit))
}

mm_boundary_df_reason_code <- function(method) {
  switch(
    method,
    satterthwaite = "satterthwaite_unavailable_at_boundary",
    kenward_roger = "kenward_roger_unavailable_at_boundary",
    "finite_sample_df_unavailable_at_boundary"
  )
}

mm_boundary_df_reason <- function(method) {
  sprintf(
    "%s degrees-of-freedom inference is not certified for boundary or reduced-rank variance-parameter fits.",
    method
  )
}

mm_boundary_df_unavailable_contrast_table <- function(fit, L, rhs, method) {
  reason_code <- mm_boundary_df_reason_code(method)
  reason <- mm_boundary_df_reason(method)
  n <- nrow(L)
  data.frame(
    contrast = rownames(L),
    estimate = as.numeric(L %*% fit$beta) - rhs,
    rhs = rhs,
    std_error = NA_real_,
    df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = method,
    requested_method = method,
    status = "not_assessed",
    reliability = "not_available",
    reliability_reason = reason_code,
    estimability = I(rep(list(NULL), n)),
    reason = reason,
    reason_code = reason_code,
    reason_detail = fit_status(fit),
    details = I(rep(list(list(
      boundary_policy = reason_code,
      fit_status = fit_status(fit)
    )), n)),
    notes = I(rep(list(character()), n)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mm_boundary_df_unavailable_effect_table <- function(term, method) {
  reason_code <- mm_boundary_df_reason_code(method)
  reason <- mm_boundary_df_reason(method)
  n <- length(term)
  data.frame(
    term = term,
    num_df = NA_real_,
    den_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = method,
    requested_method = method,
    status = "not_assessed",
    reliability = "not_available",
    reliability_reason = reason_code,
    reason = reason,
    reason_code = reason_code,
    reason_detail = NA_character_,
    details = I(rep(list(list(boundary_policy = reason_code)), n)),
    notes = I(rep(list(character()), n)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mm_cluster_bootstrap_unavailable_effect_table <- function(fit, term, group = NULL) {
  groups <- names(ranef(fit))
  if (is.null(group)) {
    if (length(groups) == 1L) {
      group <- groups[[1L]]
      reason_code <- "bootstrap_cluster_resample_p_value_unavailable"
      reason <- paste(
        "cluster_resample is an estimator-distribution target in schema 1.0.0;",
        "it does not certify fixed-effect hypothesis-test p-values"
      )
    } else {
      group <- NA_character_
      reason_code <- "cluster_bootstrap_multifactor_ambiguous"
      reason <- paste(
        "cluster bootstrap requires an explicit `group` when a model has",
        "multiple random-effect grouping factors"
      )
    }
  } else {
    if (!is.character(group) || length(group) != 1L || is.na(group) || !nzchar(group)) {
      mm_abort(
        message = "`group` must be a single grouping-factor name.",
        class = "mm_arg_error",
        input = group
      )
    }
    if (!group %in% groups) {
      mm_abort(
        message = sprintf("Unknown random-effect grouping factor `%s`. Known groups: %s.",
                          group, paste(groups, collapse = ", ")),
        class = "mm_arg_error",
        input = group
      )
    }
    reason_code <- "bootstrap_cluster_resample_p_value_unavailable"
    reason <- paste(
      "cluster_resample is an estimator-distribution target in schema 1.0.0;",
      "it does not certify fixed-effect hypothesis-test p-values"
    )
  }

  data.frame(
    term = term,
    num_df = NA_real_,
    den_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = "cluster_bootstrap",
    requested_method = "cluster_bootstrap",
    status = "not_assessed",
    reliability = "not_available",
    reliability_reason = "not_available",
    reason = reason,
    reason_code = reason_code,
    reason_detail = reason,
    details = I(rep(list(list(
      bootstrap = list(
        target_kind = "cluster_resample",
        group = group,
        p_value_certified = FALSE
      )
    )), length(term))),
    notes = I(rep(list(character()), length(term))),
    stringsAsFactors = FALSE
  )
}

mm_fixed_effect_terms <- function(fit) {
  terms <- as.character(unlist(
    fit$artifact$semantic_model$fixed_terms %||% list(),
    use.names = FALSE
  ))
  if (!length(terms)) {
    terms <- names(fit$beta)
  }
  terms
}

mm_inference_unavailable_reason <- function(method) {
  if (identical(method, "none")) {
    return("inference_not_requested")
  }
  "rust_inference_certificate_unavailable"
}

mm_bootstrap_confint <- function(fit, parm, level, bootstrap, interval) {
  if (is.null(bootstrap)) bootstrap <- bootstrap_control()
  if (!inherits(bootstrap, "mm_bootstrap_control")) {
    bootstrap <- do.call(bootstrap_control, as.list(bootstrap))
  }
  rows <- lapply(parm, function(parameter) {
    payload <- mm_full_model_bootstrap_payload(fit, parameter, level, bootstrap)
    interval_row <- mm_select_bootstrap_interval(payload, level, interval)
    list(
      parameter = parameter,
      lower = interval_row$lower,
      upper = interval_row$upper,
      payload = payload
    )
  })

  out <- cbind(
    vapply(rows, `[[`, numeric(1), "lower"),
    vapply(rows, `[[`, numeric(1), "upper")
  )
  rownames(out) <- vapply(rows, `[[`, character(1), "parameter")
  alpha <- 1 - level
  colnames(out) <- c(sprintf("%.1f %%", 100 * alpha / 2),
                     sprintf("%.1f %%", 100 * (1 - alpha / 2)))
  attr(out, "method") <- "bootstrap_full_model_distribution"
  attr(out, "interval") <- interval
  attr(out, "status") <- "available"
  attr(out, "bootstrap") <- lapply(rows, `[[`, "payload")
  mm_new_confint(out)
}

# Rows of a profile interval matrix whose per-parameter profile was
# irregular (engine status non_monotone / not_bracketing / failed). Their
# undetermined bounds are NA; the reason is in attr(x, "mm_profile")$table.
mm_profile_irregular_rows <- function(x) {
  tab <- attr(x, "mm_profile")$table
  if (is.null(tab) || !"status" %in% names(tab)) return(NULL)
  tab <- tab[tab$parameter %in% rownames(x), , drop = FALSE]
  tab[!is.na(tab$status) & tab$status != "ok", , drop = FALSE]
}

mm_print_profile_irregular_note <- function(x) {
  bad <- mm_profile_irregular_rows(x)
  if (is.null(bad) || !nrow(bad)) return(invisible(NULL))
  na_rows <- is.na(bad$lower) | is.na(bad$upper)
  what <- sprintf("%s (%s)", bad$parameter, gsub("_", " ", bad$reason_code,
                                                 fixed = TRUE))
  cat(strwrap(
    paste0(
      "Note: irregular profile for ", paste(what, collapse = ", "), ". ",
      if (any(na_rows)) {
        paste0("Bounds the profile could not determine are reported as NA ",
               "(lme4 instead falls back to linear interpolation, often ",
               "giving [-1, 1] for a correlation). ")
      } else {
        "The reported bounds come from the regular part of the profile. "
      },
      "Reasons: attr(x, \"mm_profile\")$table$reason."
    ),
    width = 78, exdent = 2
  ), sep = "\n")
  invisible(NULL)
}

mm_new_confint <- function(x) {
  class(x) <- c("mm_confint", "matrix")
  x
}

mm_bootstrap_confint_summary <- function(x, bootstrap) {
  parameters <- rownames(x)
  if (is.null(parameters)) {
    parameters <- rep(NA_character_, length(bootstrap))
  }
  rows <- lapply(seq_along(bootstrap), function(i) {
    metadata <- bootstrap[[i]]$metadata %||% list()
    seed <- metadata$seed_record$seed %||% NA
    data.frame(
      parameter = parameters[[i]],
      requested = as.integer(metadata$requested_replicates %||% NA_integer_),
      successful = as.integer(metadata$successful_replicates %||% NA_integer_),
      failed_refits = as.integer(metadata$failed_refits %||% NA_integer_),
      boundary_rate = round(as.numeric(metadata$boundary_rate %||% NA_real_), 3),
      seed = as.integer(seed),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mm_full_model_bootstrap_payload <- function(fit, parameter, level, bootstrap) {
  L <- matrix(0, nrow = 1L, ncol = length(fit$beta))
  colnames(L) <- names(fit$beta)
  rownames(L) <- parameter
  L[1L, parameter] <- 1
  bridge <- mm_rust_fit_bridge_payload(fit)
  bootstrap_json <- jsonlite::toJSON(
    mm_bootstrap_wire(bootstrap),
    auto_unbox = TRUE,
    null = "null"
  )
  json <- tryCatch(
    mm_full_model_bootstrap_contrast_json(
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
      as.character(rownames(L)),
      0,
      as.character(bootstrap_json),
      as.numeric(level)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = "bootstrap")
  }
  tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf("Failed to parse full-model bootstrap JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
}

mm_select_bootstrap_interval <- function(payload, level, interval) {
  intervals <- payload$intervals %||% list()
  hits <- vapply(intervals, function(x) {
    identical(mm_scalar_text(x$method), interval) &&
      isTRUE(all.equal(as.numeric(x$level), level, tolerance = 1e-12))
  }, logical(1))
  if (!any(hits)) {
    mm_abort(
      message = sprintf("Bootstrap payload did not contain a %s interval at level %.3f.",
                        interval, level),
      class = "mm_schema_error",
      input = payload
    )
  }
  row <- intervals[[which(hits)[[1L]]]]
  list(
    lower = as.numeric(row$lower),
    upper = as.numeric(row$upper)
  )
}

# Map a term name to a row-permutation L matrix that picks the coefficients
# belonging to the term. Selection goes through the coefficient map's
# `assign` vector (model.matrix's column -> term index), not name parsing:
# with lme4-style names ("recipeB") the level is fused to the variable name,
# so no prefix regex can distinguish term "SNR" from a coefficient of a term
# named "SNRHard". The intercept term `1` selects assign == 0.
mm_term_to_l_matrix <- function(fit, term) {
  coef_names <- names(fit$beta)
  k <- length(coef_names)
  map <- mm_fit_coef_map(fit)
  hits <- if (identical(term, "1")) {
    map$assign == 0L
  } else {
    map$assign == match(term, map$term_labels)
  }
  hits[is.na(hits)] <- FALSE
  if (!any(hits)) {
    mm_abort(
      message = sprintf(
        "Term '%s' did not match any fixed-effect coefficient. Coefficients: %s",
        term, paste(coef_names, collapse = ", ")),
      class = "mm_inference_unavailable",
      input = term
    )
  }
  idx <- which(hits)
  L <- matrix(0, nrow = length(idx), ncol = k)
  for (i in seq_along(idx)) L[i, idx[[i]]] <- 1
  rownames(L) <- coef_names[idx]
  colnames(L) <- coef_names
  L
}

regex_escape <- function(s) {
  gsub("([.\\\\+*?\\[\\^\\]$(){}=!<>|:#-])", "\\\\\\1", s, perl = TRUE)
}

# Dispatch a term-level bootstrap test through the new Rust bridge entrypoint.
# Single-df terms produce a t-form row; multi-df terms produce an F-form row
# with `numerator_df` set to the effective restriction rank.
mm_rust_term_bootstrap_row <- function(fit, term, bootstrap) {
  if (is.null(bootstrap)) bootstrap <- bootstrap_control()
  if (!inherits(bootstrap, "mm_bootstrap_control")) {
    bootstrap <- do.call(bootstrap_control, as.list(bootstrap))
  }
  L <- mm_term_to_l_matrix(fit, term)
  rhs <- rep(0, nrow(L))
  mm_inform_bootstrap_scale(fit, bootstrap)
  bridge <- mm_rust_fit_bridge_payload(fit)
  bootstrap_json <- jsonlite::toJSON(
    mm_bootstrap_wire(bootstrap),
    auto_unbox = TRUE,
    null = "null"
  )
  json <- tryCatch(
    mm_fixed_effect_bootstrap_term_json(
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
      as.character(term),
      as.numeric(rhs),
      as.character(bootstrap_json)
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json, method = "bootstrap")
  }
  table <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  parsed <- mm_json_parse_fixed_effect_inference_table(table)
  row <- parsed$table[1L, , drop = FALSE]
  row$term <- term
  row
}

# Dispatch a bootstrap likelihood-ratio test for `term`. Builds the reduced
# formula (drops the term from the fixed effects), refuses REML with a stable
# reason, and routes through the new `mm_bootstrap_lrt_json` bridge entrypoint.
# Returns a single inference-table-shaped row.
mm_rust_term_bootstrap_lrt_row <- function(fit, term, bootstrap) {
  if (isTRUE(fit$REML)) {
    return(mm_inference_row_unavailable(
      term = term,
      method = "bootstrap_lrt",
      reason = paste(
        "bootstrap likelihood-ratio test requires ML-fitted models;",
        "this fit is REML -- refit with `lmm(..., REML = FALSE)` and retry"
      ),
      reason_code = "bootstrap_lrt_requires_ml"
    ))
  }
  if (is.null(bootstrap)) bootstrap <- bootstrap_control()
  if (!inherits(bootstrap, "mm_bootstrap_control")) {
    bootstrap <- do.call(bootstrap_control, as.list(bootstrap))
  }
  reduced_formula <- tryCatch(
    update(fit$formula, stats::as.formula(sprintf(". ~ . - %s", term))),
    error = function(e) e
  )
  if (inherits(reduced_formula, "condition")) {
    return(mm_inference_row_unavailable(
      term = term,
      method = "bootstrap_lrt",
      reason = sprintf("could not drop term '%s' from formula: %s",
                       term, conditionMessage(reduced_formula)),
      reason_code = "bootstrap_lrt_reduced_formula_failed"
    ))
  }
  bridge <- mm_rust_fit_bridge_payload(fit)
  bootstrap_json <- jsonlite::toJSON(
    mm_bootstrap_wire(bootstrap),
    auto_unbox = TRUE,
    null = "null"
  )
  json <- tryCatch(
    mm_bootstrap_lrt_json(
      deparse1(reduced_formula),
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
    return(mm_inference_row_unavailable(
      term = term,
      method = "bootstrap_lrt",
      reason = conditionMessage(json),
      reason_code = "bootstrap_lrt_engine_refused"
    ))
  }
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  payload <- parsed$payload %||% list()
  grade <- mm_bootstrap_reliability(
    certified = !is.null(parsed$p_value),
    successful = payload$metadata$successful_replicates %||% NA_integer_,
    mcse = parsed$mcse %||% NA_real_
  )
  data.frame(
    term = term,
    label = term,
    kind = "coefficient",
    estimate = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    numerator_df = NA_real_,
    denominator_df = NA_real_,
    statistic = parsed$observed_statistic %||% NA_real_,
    statistic_name = "chi_square",
    p_value = parsed$p_value %||% NA_real_,
    method = "bootstrap_lrt",
    status = if (is.null(parsed$p_value)) "not_assessed" else "available",
    reliability = grade$reliability,
    reliability_reason = grade$reason,
    reason = NA_character_,
    reason_code = NA_character_,
    reason_detail = NA_character_,
    estimability = I(list(NULL)),
    details = I(list(list(
      bootstrap = list(
        target_kind = payload$metadata$target$kind %||%
          payload$target_kind %||% "likelihood_ratio",
        observed_statistic = parsed$observed_statistic %||% NA_real_,
        mcse = parsed$mcse %||% NA_real_,
        successful_replicates = payload$metadata$successful_replicates %||% NA_integer_,
        completed_replicates = payload$metadata$completed_replicates %||% NA_integer_,
        replicate_statistics = payload$replicate_statistics %||% list(),
        boundary_count = payload$metadata$boundary_count %||% NA_integer_,
        seed = payload$metadata$seed_record$seed %||% NA_integer_
      ),
      reduced_formula = deparse1(reduced_formula),
      alternative_formula = bridge$formula_string
    ))),
    notes = I(list(as.character(unlist(parsed$notes %||% list(),
                                       use.names = FALSE)))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# Replicate floor for a `moderate` bootstrap grade -- ONE definition,
# consumed by mm_bootstrap_reliability() and the inference_options route
# map (previously the literal 999 was duplicated in two files).
#
# Rationale for 999 (documented in the glossary vignette's reliability
# rubric): (i) with B = 999 the achievable two-sided bootstrap p-values are
# exact multiples of 1/(B+1) = 0.001, so a reported p of 0.05 is resolved
# to +/- one grid step; (ii) the Monte Carlo SE of an estimated rejection/
# coverage rate near 0.05 is sqrt(.05*.95/999) ~ 0.0069 -- small relative
# to the 0.05 target; below ~1000 replicates that MC noise becomes a
# material share of the quantity being reported. The step at the threshold
# is a grading convention, not a claim that B = 998 differs materially
# from B = 999; the accounting fields (successful replicates, MCSE) ship
# alongside the grade so readers can judge the evidence directly.
mm_bootstrap_moderate_min <- function() 999L

# Grade bootstrap reliability from the replicate accounting. `moderate`
# requires a certified payload, a finite Monte-Carlo standard error, and at
# least mm_bootstrap_moderate_min() successful replicates; fewer replicates
# or a non-finite MCSE are `low`; an uncertified payload is `not_available`.
mm_bootstrap_reliability <- function(certified, successful, mcse,
                                     min_moderate = mm_bootstrap_moderate_min()) {
  if (!isTRUE(certified)) {
    return(list(reliability = "not_available",
                reason = "bootstrap_not_certified"))
  }
  succ_ok <- is.numeric(successful) && length(successful) == 1L &&
    !is.na(successful) && successful >= min_moderate
  mcse_ok <- is.numeric(mcse) && length(mcse) == 1L && is.finite(mcse)
  if (succ_ok && mcse_ok) {
    list(reliability = "moderate",
         reason = "bootstrap_monte_carlo_replicates")
  } else {
    list(reliability = "low",
         reason = "bootstrap_insufficient_replicates")
  }
}

# Construct a one-row "unavailable" inference-table-shaped data frame.
mm_inference_row_unavailable <- function(term, method, reason, reason_code = NA_character_) {
  data.frame(
    term = term,
    label = term,
    kind = "coefficient",
    estimate = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    numerator_df = NA_real_,
    denominator_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = method,
    status = "not_assessed",
    reliability = "not_available",
    reliability_reason = "inference_unavailable_by_policy",
    reason = reason,
    reason_code = reason_code,
    reason_detail = NA_character_,
    estimability = I(list(NULL)),
    details = I(list(NULL)),
    notes = I(list(character())),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Profile a fitted linear mixed model
#'
#' Computes profile-likelihood intervals for the model's parameters via the
#' engine's certified profile payload and returns them as an `mm_profile`
#' object: `$table` has one row per profiled parameter (`parameter`,
#' `parameter_kind`, `estimate`, `lower`, `upper`, `regularity`, `status`,
#' `reason_code`, `reason`, `profiled_criterion`). Rows use lme4's names: `.sig01`,
#' ... (kind `"sd"` or `"cor"`, numbered in lme4's term order), `.sigma`,
#' and the fixed effects, all profiled on the ML deviance as lme4 does (a
#' REML fit is refitted by ML); the engine's relative-Cholesky `theta1`, ...
#' rows follow, profiled on the fit's own criterion. Use [confint()] with
#' `method = "profile"` for the matrix form.
#'
#' A parameter whose profile is irregular keeps its row: `status` is
#' `"non_monotone"` (the signed-root deviance is not monotone),
#' `"not_bracketing"` (the profile does not reach the cutoff on one side),
#' or `"failed"` (a constrained refit failed); `reason_code` holds the
#' matching `regularity` code (`"non_monotone_profile"`,
#' `"profile_not_bracketing"`, `"profile_failed"`) and `reason` the engine's
#' explanation. Bounds the profile cannot determine are `NA`; a bound inside
#' the regular part of the profile is still reported. lme4 instead warns
#' and falls back to linear interpolation (often giving `[-1, 1]` for a
#' correlation); mixeff reports `NA` plus the reason rather than an
#' interpolated bound. The other parameters' intervals are unaffected.
#'
#' @param fitted A fitted `mm_lmm`.
#' @param which Optional character vector of parameter names to keep
#'   (coefficient names, `".sig01"`, `".sigma"`, `"theta1"`, ...).
#' @param level Confidence level for the reported interval endpoints.
#' @param threads Worker threads (default 1); see [confint.mm_lmm()].
#' @param ... Unused; for generic consistency.
#' @return An `mm_profile` object with `$table`, `$level`, `$fit_criterion`,
#'   and `$notes`.
#' @seealso [confint()] with `method = "profile"` for the matrix form.
#' @export
profile.mm_lmm <- function(fitted, which = NULL, level = 0.95,
                           threads = 1L, ...) {
  ci <- mm_profile_confint(fitted, parm = which, level = level,
                           threads = threads)
  prof <- attr(ci, "mm_profile")
  out <- list(
    table = prof$table,
    level = prof$level,
    fit_criterion = prof$fit_criterion,
    notes = prof$notes
  )
  class(out) <- "mm_profile"
  out
}

#' @method print mm_profile
#' @export
print.mm_profile <- function(x, ...) {
  cat(sprintf("Profile-likelihood intervals (%s, level %.2f)\n",
              x$fit_criterion, x$level))
  tbl <- x$table
  show <- tbl[, c("parameter", "estimate", "lower", "upper"), drop = FALSE]
  print(format(show, digits = 4), row.names = FALSE)
  unavailable <- tbl[!is.na(tbl$reason_code), , drop = FALSE]
  if (nrow(unavailable)) {
    cat(sprintf("Irregular profile: %s.\n",
                paste(sprintf("%s (%s)", unavailable$parameter,
                              unavailable$reason_code),
                      collapse = ", ")))
  }
  for (note in x$notes) cat("Note:", note, "\n")
  invisible(x)
}

#' @export
confint.mm_profile <- function(object, parm, level = NULL, ...) {
  if (!is.null(level) && !isTRUE(all.equal(level, object$level))) {
    mm_abort(
      message = sprintf(
        "This profile was computed at level %.3f; re-run profile() with level = %.3f.",
        object$level, level
      ),
      class = "mm_arg_error",
      input = level
    )
  }
  tbl <- object$table
  if (!missing(parm) && !is.null(parm)) {
    tbl <- tbl[tbl$parameter %in% as.character(parm), , drop = FALSE]
  }
  mat <- as.matrix(tbl[, c("lower", "upper"), drop = FALSE])
  rownames(mat) <- tbl$parameter
  alpha <- 1 - object$level
  colnames(mat) <- c(sprintf("%.1f %%", 100 * alpha / 2),
                     sprintf("%.1f %%", 100 * (1 - alpha / 2)))
  mat
}

# Stage D.3 (bd-01KRFGFSK4A0MGPFQVCNY5SYFK): wrapper for the upstream
# `profile_confint_payload` FFI. Refits the model from the bridge payload
# (no persistent Rust handle), parses the schema-versioned JSON, and
# returns a list with the rendered confint matrix plus the raw payload.

mm_profile_confint <- function(fit, parm = NULL, level = 0.95, threads = 1L) {
  threads <- mm_check_threads(threads)
  reml <- isTRUE(fit$REML)
  # lme4 profiles every fit on the ML deviance (a REML fit's confint()
  # equals its ML refit's), so the lme4-scale rows (.sigNN, .sigma, beta)
  # come from the ML profile. The relative-Cholesky theta rows stay on the
  # fit's own criterion.
  payload_ml <- mm_profile_confint_payload(fit, level, threads, reml = FALSE)
  payload <- if (reml) {
    mm_profile_confint_payload(fit, level, threads, reml = TRUE)
  } else {
    payload_ml
  }
  fit_criterion <- as.character(payload$fit_criterion %||% "")

  translate <- function(p) {
    rows <- lapply(p$intervals %||% list(), mm_translate_profile_row, fit = fit)
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) do.call(rbind, rows) else mm_empty_profile_table()
  }
  tab_ml <- translate(payload_ml)
  tab_fit <- if (reml) translate(payload) else tab_ml
  lme4_kinds <- c("sd", "cor", "sigma", "beta")
  sig <- tab_ml[tab_ml$parameter_kind %in% c("sd", "cor"), , drop = FALSE]
  sig <- sig[order(sig$parameter), , drop = FALSE]
  sigma_row <- tab_ml[tab_ml$parameter_kind == "sigma", , drop = FALSE]
  beta <- tab_ml[tab_ml$parameter_kind == "beta", , drop = FALSE]
  beta <- beta[order(match(beta$parameter, names(fit$beta))), , drop = FALSE]
  other <- tab_fit[!tab_fit$parameter_kind %in% lme4_kinds, , drop = FALSE]
  table <- rbind(sig, sigma_row, beta, other)
  table$profiled_criterion <- c(rep("ML", nrow(sig) + nrow(sigma_row) + nrow(beta)),
                                rep(fit_criterion, nrow(other)))

  # parm: names (lme4's ".sig01", ".sigma", coefficient names, or the theta
  # rows "theta1", ...), lme4's "theta_" (all variance parameters) and
  # "beta_" (all fixed effects), or indices into the lme4 rows.
  is_lme4 <- table$parameter_kind %in% lme4_kinds
  if (!is.null(parm)) {
    if (is.numeric(parm)) {
      keep <- table$parameter[is_lme4][parm]
    } else {
      parm <- as.character(parm)
      keep <- c(
        setdiff(parm, c("theta_", "beta_")),
        if ("theta_" %in% parm) {
          table$parameter[table$parameter_kind %in% c("sd", "cor", "sigma")]
        },
        if ("beta_" %in% parm) table$parameter[table$parameter_kind == "beta"]
      )
    }
    table <- table[table$parameter %in% keep, , drop = FALSE]
  }

  rownames(table) <- NULL
  shown <- if (is.null(parm)) {
    table[table$parameter_kind %in% lme4_kinds, , drop = FALSE]
  } else {
    table
  }
  mat <- as.matrix(shown[, c("lower", "upper"), drop = FALSE])
  rownames(mat) <- shown$parameter
  alpha <- 1 - level
  colnames(mat) <- c(sprintf("%.1f %%", 100 * alpha / 2),
                     sprintf("%.1f %%", 100 * (1 - alpha / 2)))
  attr(mat, "method") <- "profile_likelihood"
  attr(mat, "status") <- "available"
  attr(mat, "fit_criterion") <- fit_criterion
  notes <- unique(c(
    as.character(unlist(payload$notes %||% list(), use.names = FALSE)),
    as.character(unlist(payload_ml$notes %||% list(), use.names = FALSE))
  ))
  # The engine's REML-payload note on beta no longer applies: the beta rows
  # come from the ML profile, as in lme4.
  notes <- notes[!grepl("omit fixed-effect beta profiles", notes, fixed = TRUE)]
  if (reml) {
    notes <- c(notes, paste0(
      "REML fit: the .sigNN, .sigma and fixed-effect intervals profile the ",
      "ML deviance (the model refitted by ML), as lme4's confint() does; ",
      "theta rows profile the REML criterion."
    ))
  }
  attr(mat, "mm_profile") <- list(
    schema = list(
      schema_name = as.character(payload$schema_name %||% NA_character_),
      schema_version = as.character(payload$schema_version %||% NA_character_)
    ),
    fit_criterion = fit_criterion,
    level = as.numeric(payload$level %||% level),
    threads = threads,
    notes = notes,
    table = table
  )
  mm_new_confint(mat)
}

mm_profile_confint_payload <- function(fit, level, threads = 1L,
                                       reml = isTRUE(fit$REML)) {
  bridge <- mm_rust_fit_bridge_payload(fit)
  json <- tryCatch(
    .Call(
      wrap__mm_lmm_profile_confint_json,
      bridge$formula_string,
      isTRUE(reml),
      bridge$spec_data$column_order,
      bridge$spec_data$numeric_columns,
      bridge$spec_data$categorical_values,
      bridge$spec_data$categorical_levels,
      bridge$spec_data$categorical_ordered,
      bridge$weights,
      bridge$control_json,
      as.numeric(level),
      as.integer(threads)
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
        message = sprintf("Failed to parse profile CI JSON: %s",
                          conditionMessage(cnd)),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
  schema_name <- as.character(payload$schema_name %||% "")
  if (!nzchar(schema_name)) {
    mm_abort(
      message = "Profile CI payload is missing its schema header.",
      class = "mm_schema_error",
      input = payload
    )
  }
  mm_json_negotiate(list(
    schema_name = schema_name,
    schema_version = as.character(payload$schema_version %||% NA_character_)
  ))
  payload
}

# Translate one upstream profile row from its parameter-name encoding
# ("β1", "β2", "σ", "θ1", ...) into mixeff-side coordinates so the matrix
# row names align with fit$beta and the rest of the wrapper.
mm_translate_profile_row <- function(row, fit) {
  upstream <- as.character(row$parameter %||% "")
  mapped <- mm_map_profile_parameter(upstream, fit)
  if (is.null(mapped)) return(NULL)
  # Irregular per-parameter profiles (engine 26f9690+) arrive as rows with
  # `status` != "ok", a `reason`, and JSON-null (-> NA) bounds that could not
  # be determined; they are reported, never dropped.
  status <- as.character(row$status %||% "ok")
  regularity <- as.character(
    row$regularity %||%
      mm_profile_status_regularity[status] %||%
      "regular_profile_likelihood"
  )
  if (is.na(regularity)) regularity <- "regular_profile_likelihood"
  data.frame(
    parameter = mapped$name,
    parameter_kind = mapped$kind,
    estimate = as.numeric(row$estimate %||% NA_real_),
    lower = as.numeric(row$lower %||% NA_real_),
    upper = as.numeric(row$upper %||% NA_real_),
    method = as.character(row$method %||% "profile_likelihood"),
    regularity = regularity,
    boundary_clamped_lower = isTRUE(row$boundary_clamped_lower),
    status = status,
    reason_code = if (identical(status, "ok")) NA_character_ else regularity,
    reason = mm_profile_reason_text(row$reason, upstream, mapped$name),
    stringsAsFactors = FALSE
  )
}

# The engine prefixes a failed profile's reason with its own parameter name
# ("profile_.sig02: ..."; engine theta order). Relabel that prefix with the
# lme4 name the row carries so the reason matches the row.
mm_profile_reason_text <- function(reason, upstream, mapped) {
  reason <- mm_optional_text(reason)
  if (is.na(reason) || identical(upstream, mapped)) return(reason)
  prefix <- paste0("profile_", upstream, ":")
  if (startsWith(reason, prefix)) {
    reason <- paste0("profile_", mapped, ":",
                     substr(reason, nchar(prefix) + 1L, nchar(reason)))
  }
  reason
}

# Engine (mixeff-rs 26f9690+) per-parameter profile status -> the closed
# `regularity` code it travels with. A payload from an older engine has no
# `status`; its rows are all regular.
mm_profile_status_regularity <- c(
  non_monotone = "non_monotone_profile",
  not_bracketing = "profile_not_bracketing",
  failed = "profile_failed"
)

# lme4 (>= 2.0) numbers profile parameters term by term: a term's standard
# deviations first, then its correlations (lower triangle, column-major).
# `pos` is the lme4 theta index of the slot.
mm_profile_sig_slot <- function(st, pos) {
  offset <- 0L
  for (tm in st$terms) {
    cells <- which(tm$Tidx > 0L, arr.ind = TRUE)
    if (!nrow(cells)) next
    is_diag <- cells[, 1L] == cells[, 2L]
    ord <- order(!is_diag, cells[, 2L], cells[, 1L])
    hit <- match(pos, tm$Tidx[cells[ord, , drop = FALSE]])
    if (!is.na(hit)) {
      return(list(index = offset + hit, diag = is_diag[ord][[hit]]))
    }
    offset <- offset + nrow(cells)
  }
  NULL
}

mm_map_profile_parameter <- function(upstream, fit) {
  if (upstream %in% c("\u03C3", "<U+03C3>", ".sigma")) {
    return(list(name = ".sigma", kind = "sigma"))
  }
  sig_idx <- if (grepl("^\\.sig[0-9]+$", upstream)) {
    as.integer(sub("^\\.sig", "", upstream))
  } else {
    NA_integer_
  }
  if (!is.na(sig_idx)) {
    # Engine `.sigNN` = theta slot NN in ENGINE theta order (SD for a
    # diagonal slot, correlation otherwise); lme4 numbers the same
    # quantities in its own term order (mm_re_structure()$perm).
    st <- tryCatch(mm_re_structure(fit), error = function(cnd) NULL)
    if (is.null(st)) return(list(name = upstream, kind = "unknown"))
    pos <- match(sig_idx, st$perm)
    if (is.na(pos)) return(list(name = upstream, kind = "unknown"))
    sig <- mm_profile_sig_slot(st, pos)
    if (is.null(sig)) return(list(name = upstream, kind = "unknown"))
    return(list(name = sprintf(".sig%02d", sig$index),
                kind = if (sig$diag) "sd" else "cor"))
  }
  beta_idx <- mm_profile_parameter_index(upstream, c("\u03B2", "<U+03B2>"))
  if (!is.na(beta_idx)) {
    beta_names <- names(fit$beta)
    if (beta_idx >= 1L && beta_idx <= length(beta_names)) {
      # The engine indexes beta in ITS column order; fit$beta is stored in
      # lme4 order, so route the positional index through the map (NA when
      # no map is stored, e.g. legacy fits: keep the raw index).
      lme4_pos <- match(beta_idx, fit$coef_map$perm)
      if (is.na(lme4_pos)) lme4_pos <- beta_idx
      return(list(name = beta_names[[lme4_pos]], kind = "beta"))
    }
    return(NULL)
  }
  theta_idx <- mm_profile_parameter_index(upstream, c("\u03B8", "<U+03B8>"))
  if (!is.na(theta_idx)) {
    return(list(name = sprintf("theta%d", theta_idx), kind = "theta"))
  }
  # Unrecognized parameter name; do not silently drop — surface a stable
  # synthetic label so callers see the row in attr(ci, "mm_profile")$table.
  list(name = upstream, kind = "unknown")
}

mm_profile_parameter_index <- function(upstream, prefix) {
  for (one_prefix in prefix) {
    if (!startsWith(upstream, one_prefix)) next
    suffix <- substr(upstream, nchar(one_prefix) + 1L, nchar(upstream))
    if (!grepl("^[0-9]+$", suffix)) return(NA_integer_)
    return(suppressWarnings(as.integer(suffix)))
  }
  NA_integer_
}

mm_profile_beta_subset <- function(fit, parm) {
  beta_names <- names(fit$beta)
  if (is.null(parm)) return(beta_names)
  parm <- as.character(parm)
  intersect(beta_names, parm)
}

mm_empty_profile_table <- function() {
  data.frame(
    parameter = character(),
    parameter_kind = character(),
    estimate = numeric(),
    lower = numeric(),
    upper = numeric(),
    method = character(),
    regularity = character(),
    boundary_clamped_lower = logical(),
    status = character(),
    reason_code = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

#' Wald inference on a linear combination of fixed effects
#'
#' Convenience helper for the common case of testing
#' \eqn{H_0:\, c^\top \beta = 0}{H0: c' beta = 0} where `c` is a sparse,
#' named weight vector. The estimate is \eqn{c^\top \hat\beta}{c' beta_hat},
#' the standard error is \eqn{\sqrt{c^\top V c}}{sqrt(c' V c)} where `V`
#' is the model's fixed-effect covariance, the statistic is the Wald
#' ratio, and the interval is the symmetric Wald CI at `level`.
#'
#' For `mm_glmm`, the statistic is the asymptotic Wald *z* (no df). For
#' `mm_lmm`, the default is Satterthwaite-approximated *t* via
#' [df_for_contrast()]; pass `method = "asymptotic"` to force Wald *z*.
#'
#' Weight names must be a subset of `names(fixef(fit))`. Coefficients
#' not named in `weights` contribute zero. Pass the long-form
#' [contrast()] front door if you need multiple contrasts or non-default
#' rhs.
#'
#' @param fit A fitted `mm_lmm` or `mm_glmm`.
#' @param weights A named numeric vector (or named list / single-row
#'   data.frame coercible to one). Names must match
#'   `names(fixef(fit))` exactly.
#' @param level Confidence level for the Wald interval. Default 0.95.
#' @param method For `mm_lmm`, the degrees-of-freedom method passed to
#'   [df_for_contrast()]. Defaults to `"auto"` (Satterthwaite when
#'   available). For `mm_glmm`, only `"asymptotic"` is accepted.
#' @param ... Reserved for future methods.
#'
#' @return A single-row data.frame with columns `estimate`, `std_error`,
#'   `statistic`, `statistic_name` (`"t"` or `"z"`), `df`, `p_value`,
#'   `lower`, `upper`, and `method`. The result carries an `"mm_status"`
#'   attribute reflecting the underlying vcov reliability (`status`,
#'   `method`, `reliability`, `reason`).
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   g   = factor(rep(1:12, each = 8)),
#'   trt = factor(rep(c("ctrl", "drug"), 48)),
#'   x   = rnorm(96)
#' )
#' d$y <- 1 + 0.5 * (d$trt == "drug") + 0.3 * d$x +
#'   0.4 * (d$trt == "drug") * d$x + rnorm(12)[d$g] + rnorm(96)
#' fit <- lmm(y ~ trt * x + (1 | g), d, control = mm_control(verbose = -1))
#' # Treatment effect at x = 1 (weights use lme4 coefficient names):
#' mm_lincomb(fit, c(trtdrug = 1, "trtdrug:x" = 1))
#'
#' @seealso [contrast()] for the long-form, Rust-routed contrast surface
#'   with full estimability / reliability reporting.
#'
#' @export
mm_lincomb <- function(fit, weights, level = 0.95, method = NULL, ...) {
  UseMethod("mm_lincomb")
}

#' @rdname mm_lincomb
#' @export
mm_lincomb.default <- function(fit, weights, level = 0.95, method = NULL, ...) {
  mm_abort(
    message = "`mm_lincomb()` is only defined for `mm_lmm` and `mm_glmm` fits.",
    class = "mm_arg_error",
    input = fit
  )
}

#' @rdname mm_lincomb
#' @export
mm_lincomb.mm_glmm <- function(fit, weights, level = 0.95, method = NULL, ...) {
  if (!is.null(method) && !identical(method, "asymptotic") && !identical(method, "auto")) {
    mm_abort(
      message = "`mm_lincomb()` on `mm_glmm` only supports `method = \"asymptotic\"` (Wald z).",
      class = "mm_arg_error",
      input = method
    )
  }
  # Wald capability comes from the single arbiter (WI-2.1) -- never from the
  # covariance payload's own status, which reads available_noninferential for
  # the uncertified working Hessian (see mm_glmm_inference_capability).
  cap <- mm_glmm_inference_capability(fit)
  if (!isTRUE(cap$wald)) {
    mm_abort(
      message = mm_glmm_wald_refusal_message("mm_lincomb() standard errors and Wald intervals"),
      class = "mm_inference_unavailable",
      reason_code = cap$reason_code,
      input = fit
    )
  }
  beta_named <- mm_lincomb_fixef(fit)
  beta  <- as.numeric(beta_named)
  bnms  <- names(beta_named)
  Vfull <- stats::vcov(fit)
  V     <- as.matrix(unclass(Vfull))
  dimnames(V) <- list(bnms, bnms)
  w     <- mm_lincomb_weights_vector(weights, bnms)
  est   <- sum(w * beta)
  se    <- sqrt(drop(t(w) %*% V %*% w))
  z     <- est / se
  q     <- stats::qnorm((1 + level) / 2)
  p     <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  out <- data.frame(
    estimate       = est,
    std_error      = se,
    statistic      = z,
    statistic_name = "z",
    df             = NA_real_,
    p_value        = p,
    lower          = est - q * se,
    upper          = est + q * se,
    method         = cap$method,
    stringsAsFactors = FALSE,
    check.names    = FALSE
  )
  attr(out, "mm_status") <- mm_glmm_capability_status(cap)
  out
}

#' @rdname mm_lincomb
#' @export
mm_lincomb.mm_lmm <- function(fit, weights,
                              level = 0.95,
                              method = c("auto", "satterthwaite",
                                         "kenward_roger", "asymptotic"),
                              ...) {
  method <- match.arg(method)
  beta_named <- mm_lincomb_fixef(fit)
  beta  <- as.numeric(beta_named)
  bnms  <- names(beta_named)
  Vfull <- stats::vcov(fit)
  V     <- as.matrix(unclass(Vfull))
  dimnames(V) <- list(bnms, bnms)
  w     <- mm_lincomb_weights_vector(weights, bnms)
  est   <- sum(w * beta)
  se    <- sqrt(drop(t(w) %*% V %*% w))

  if (identical(method, "asymptotic")) {
    df <- NA_real_
    stat_name <- "z"
    statistic <- est / se
    p <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
    q <- stats::qnorm((1 + level) / 2)
    method_used <- "asymptotic"
  } else {
    L <- matrix(w, nrow = 1L, dimnames = list(NULL, bnms))
    df_tbl <- df_for_contrast(fit, L, method = method)
    df <- as.numeric(df_tbl$df[[1L]])
    method_used <- df_tbl$method %||% method
    statistic <- est / se
    stat_name <- "t"
    if (is.finite(df) && df > 0) {
      p <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
      q <- stats::qt((1 + level) / 2, df = df)
    } else {
      p <- NA_real_
      q <- NA_real_
    }
  }

  out <- data.frame(
    estimate       = est,
    std_error      = se,
    statistic      = statistic,
    statistic_name = stat_name,
    df             = df,
    p_value        = p,
    lower          = if (is.finite(q)) est - q * se else NA_real_,
    upper          = if (is.finite(q)) est + q * se else NA_real_,
    method         = method_used,
    stringsAsFactors = FALSE,
    check.names    = FALSE
  )
  attr(out, "mm_status") <- mm_lincomb_status_from_vcov(Vfull)
  out
}

mm_lincomb_weights_vector <- function(weights, fixef_names) {
  if (is.null(weights)) {
    mm_abort(
      message = "`weights` must be a named numeric vector.",
      class = "mm_arg_error",
      input = weights
    )
  }
  if (is.data.frame(weights)) {
    if (nrow(weights) != 1L) {
      mm_abort(
        message = "`weights` data.frame must have exactly one row.",
        class = "mm_arg_error",
        input = weights
      )
    }
    nms <- names(weights)
    weights <- stats::setNames(as.numeric(unlist(weights[1L, ], use.names = FALSE)), nms)
  } else if (is.list(weights)) {
    nms <- names(weights)
    weights <- stats::setNames(as.numeric(unlist(weights, use.names = FALSE)), nms)
  }
  if (!is.numeric(weights) || is.null(names(weights)) || any(!nzchar(names(weights)))) {
    mm_abort(
      message = "`weights` must be a named numeric vector with non-empty names.",
      class = "mm_arg_error",
      input = weights
    )
  }
  if (anyNA(weights)) {
    mm_abort(
      message = "`weights` must not contain NA.",
      class = "mm_arg_error",
      input = weights
    )
  }
  if (anyDuplicated(names(weights))) {
    mm_abort(
      message = "`weights` names must be unique.",
      class = "mm_arg_error",
      input = weights
    )
  }
  unknown <- setdiff(names(weights), fixef_names)
  if (length(unknown)) {
    mm_abort(
      message = sprintf(
        "Unknown coefficient name(s) in `weights`: %s. Valid names: %s.",
        paste(sprintf("`%s`", unknown), collapse = ", "),
        paste(sprintf("`%s`", fixef_names), collapse = ", ")
      ),
      class = "mm_arg_error",
      input = weights
    )
  }
  out <- stats::setNames(rep(0, length(fixef_names)), fixef_names)
  out[names(weights)] <- as.numeric(weights)
  out
}

mm_lincomb_fixef <- function(fit) {
  ## fixef() returns a named numeric vector for both mm_lmm and mm_glmm;
  ## stats::coef() on mm_glmm returns the random-effects list, so we
  ## resolve fixed effects explicitly here.
  if (utils::isS3method("fixef", class = class(fit)[1L], envir = asNamespace("mixeff")) ||
      !is.null(utils::getS3method("fixef", class(fit)[1L], optional = TRUE))) {
    out <- fixef(fit)
  } else {
    out <- fit$beta
  }
  if (is.null(names(out))) {
    mm_abort(
      message = "Fitted model has unnamed fixed effects; mm_lincomb() needs names.",
      class = "mm_arg_error",
      input = fit
    )
  }
  out
}

mm_lincomb_status_from_vcov <- function(V) {
  list(
    status      = attr(V, "mm_status")      %||% "available",
    method      = attr(V, "mm_method")      %||% NA_character_,
    reliability = attr(V, "mm_reliability") %||% NA_character_,
    reason      = attr(V, "mm_reason")      %||% attr(V, "mm_unavailable_reason") %||% NA_character_
  )
}
