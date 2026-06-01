aphantasia_fixture_dir <- function() {
  candidates <- c(
    system.file("extdata", "aphantasia", package = "mixeff"),
    file.path("inst", "extdata", "aphantasia"),
    testthat::test_path("..", "..", "inst", "extdata", "aphantasia")
  )
  hit <- candidates[dir.exists(candidates)][1L]
  if (is.na(hit) || !nzchar(hit)) {
    testthat::skip("aphantasia fixture is unavailable")
  }
  hit
}

aphantasia_reference <- function() {
  jsonlite::fromJSON(
    file.path(aphantasia_fixture_dir(), "reference.json"),
    simplifyVector = FALSE
  )
}

aphantasia_trials <- function() {
  readRDS(file.path(aphantasia_fixture_dir(), "trials.rds"))
}

aphantasia_metadata <- function() {
  readRDS(file.path(aphantasia_fixture_dir(), "metadata.rds"))
}

aphantasia_run_full <- function() {
  identical(tolower(Sys.getenv("MIXEFF_RUN_APHANTASIA")), "true")
}

aphantasia_run_stress <- function() {
  identical(tolower(Sys.getenv("MIXEFF_RUN_APHANTASIA_STRESS")), "true")
}

aphantasia_run_joint_proof <- function() {
  identical(tolower(Sys.getenv("MIXEFF_APHANTASIA_JOINT_PROOF")), "true")
}

aphantasia_joint_budget <- function() {
  raw <- Sys.getenv("MIXEFF_APHANTASIA_JOINT_BUDGET", unset = "40")
  value <- suppressWarnings(as.integer(raw))
  if (is.na(value) || value < 1L) {
    stop("MIXEFF_APHANTASIA_JOINT_BUDGET must be a positive integer")
  }
  value
}

aphantasia_prepare_model_data <- function(trials, stimtype = FALSE) {
  out <- transform(
    trials,
    participant = factor(participant),
    item = factor(trial_image),
    group = factor(ifelse(aphantasia == "yes", "aphant", "control"),
                   levels = c("control", "aphant")),
    mask = factor(ifelse(back_masked == "yes", "masked", "unmasked"),
                  levels = c("unmasked", "masked")),
    block = factor(block_num),
    soa_log = log(SOA)
  )
  out$soa_s <- as.numeric(scale(out$soa_log))
  if (isTRUE(stimtype)) {
    out$stimtype <- factor(
      ifelse(out$bubbled == "yes", "occluded", "intact"),
      levels = c("intact", "occluded")
    )
  }
  out
}

aphantasia_data_sets <- function(ref, trials) {
  excluded <- unlist(ref$excluded_participants, use.names = FALSE)
  primary_raw <- subset(
    trials,
    bubbled == "yes" & !is.na(correct) & !participant %in% excluded
  )
  sensitivity_raw <- subset(trials, bubbled == "yes" & !is.na(correct))
  sensitivity_raw$aphantasia[
    sensitivity_raw$participant %in% excluded
  ] <- "no"
  intact_raw <- subset(
    trials,
    bubbled == "no" & !is.na(correct) & !participant %in% excluded
  )
  combined_raw <- subset(trials, !is.na(correct) & !participant %in% excluded)

  primary <- aphantasia_prepare_model_data(primary_raw)
  age <- subset(primary, !is.na(age))
  age$age_z <- as.numeric(scale(age$age))

  matched_ids <- unique(c(
    as.character(primary$participant[primary$aphantasia == "yes"]),
    as.character(primary$participant[
      primary$aphantasia != "yes" &
        primary$source_folder == "prolific_control_age_match"
    ])
  ))
  matched <- droplevels(primary[primary$participant %in% matched_ids, ])
  matched_age <- matched
  matched_age$age_z <- as.numeric(scale(matched_age$age))

  ## The manuscript's RT model (fit_rt_lmm) re-scales soa_s *within* the
  ## correct-trial RT subset, not on the full occluded set. Match that or
  ## the soa_s covariate differs and the LMM coefficients drift ~5e-3.
  rt <- subset(primary, correct == 1 & is.finite(rt) & rt > 0)
  rt$soa_s <- as.numeric(scale(rt$soa_log))
  rt$log_rt <- log(rt$rt)

  list(
    primary = primary,
    sensitivity = aphantasia_prepare_model_data(sensitivity_raw),
    intact = aphantasia_prepare_model_data(intact_raw),
    combined = aphantasia_prepare_model_data(combined_raw, stimtype = TRUE),
    rt = rt,
    age = age,
    matched = matched,
    matched_age = matched_age
  )
}

aphantasia_fit_cases <- function(ref, data_sets) {
  models <- ref$models
  list(
    primary = list(data = data_sets$primary, family = stats::binomial()),
    sensitivity = list(data = data_sets$sensitivity, family = stats::binomial()),
    intact = list(data = data_sets$intact, family = stats::binomial()),
    combined = list(data = data_sets$combined, family = stats::binomial()),
    rt = list(data = data_sets$rt),
    S1_intercept_only = list(data = data_sets$primary, family = stats::binomial()),
    S1_current_uncorrelated_slopes = list(data = data_sets$primary,
                                          family = stats::binomial()),
    S1_correlated_slopes = list(data = data_sets$primary,
                                family = stats::binomial()),
    S1_item_mask_slope = list(data = data_sets$primary,
                              family = stats::binomial()),
    S1_maximal = list(data = data_sets$primary, family = stats::binomial()),
    S7_age_covariate = list(data = data_sets$age, family = stats::binomial()),
    S9_age_matched_subset = list(data = data_sets$matched,
                                 family = stats::binomial()),
    S9_age_matched_subset_age_covariate = list(
      data = data_sets$matched_age,
      family = stats::binomial()
    )
  )[names(models)]
}

aphantasia_core_case_ids <- function(ref) {
  setdiff(names(ref$models), aphantasia_s1_case_ids(ref))
}

aphantasia_s1_case_ids <- function(ref) {
  grep("^S1_", names(ref$models), value = TRUE)
}

aphantasia_use_joint_glmm <- function() {
  identical(tolower(Sys.getenv("MIXEFF_APHANTASIA_JOINT")), "true")
}

aphantasia_glmm_method <- function(id) {
  if (id %in% c("intact", "combined") && aphantasia_use_joint_glmm()) {
    return("joint_laplace")
  }
  "pirls_profiled"
}

aphantasia_glmm_control <- function(id) {
  if (id %in% c("intact", "combined") && aphantasia_use_joint_glmm()) {
    return(mixeff::mm_control(verbose = -1, max_feval = aphantasia_joint_budget()))
  }
  mixeff::mm_control(verbose = -1)
}

aphantasia_lme4_key <- function(x) {
  gsub(": ", "", as.character(x), fixed = TRUE)
}

aphantasia_reference_rows <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x) && length(x) &&
      all(vapply(x, is.list, logical(1)))) {
    rows <- lapply(x, function(row) {
      as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }
  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

aphantasia_expect_fit_matches_reference <- function(fit, ref, id) {
  tol <- if (identical(ref$model_type, "lmm")) {
    unlist(aphantasia_reference()$tolerances$lmm, use.names = TRUE)
  } else {
    unlist(aphantasia_reference()$tolerances$glmm, use.names = TRUE)
  }

  observed <- mixeff::fixef(fit)
  names(observed) <- aphantasia_lme4_key(names(observed))
  expected <- unlist(ref$fixef, use.names = TRUE)
  common <- intersect(names(expected), names(observed))
  expect_equal(length(common), length(expected),
               info = sprintf("not all fixed effects aligned for `%s`", id))
  case_id <- paste0("aphantasia_", id)
  mm_assert_parity(
    observed[common],
    expected[common],
    case_id,
    "fixef",
    as.numeric(tol["fixef_abs"]),
    sprintf("aphantasia `%s` fixed effects", id),
    mode = "absolute"
  )

  loglik_ref <- ref$logLik
  aic_ref <- ref$AIC
  mm_assert_parity(
    as.numeric(stats::logLik(fit)),
    loglik_ref,
    case_id,
    "logLik",
    as.numeric(tol["logLik_rel"]),
    sprintf("aphantasia `%s` logLik", id)
  )
  mm_assert_parity(
    stats::AIC(fit),
    aic_ref,
    case_id,
    "AIC",
    as.numeric(tol["AIC_rel"]),
    sprintf("aphantasia `%s` AIC", id)
  )
}

aphantasia_has_glmm_full_vcov <- function(fit) {
  V <- stats::vcov(fit)
  identical(attr(V, "mm_status"), "available") &&
    is.matrix(V) &&
    all(is.finite(V))
}

## Map lme4-style coefficient names (e.g. "groupaphant:maskmasked") to
## mixeff's "group: aphant:mask: masked" naming so user code written
## against lme4 conventions can drive mm_lincomb() on an mm_glmm fit.
aphantasia_to_mixeff_key <- function(fit, lme4_keys) {
  mix_names <- names(mixeff::fixef(fit))
  lme4_to_mix <- setNames(mix_names, aphantasia_lme4_key(mix_names))
  unname(lme4_to_mix[lme4_keys])
}

aphantasia_lincomb <- function(fit, weights) {
  mix_keys <- aphantasia_to_mixeff_key(fit, names(weights))
  stopifnot(!anyNA(mix_keys))
  named <- setNames(as.numeric(weights), mix_keys)
  out <- mixeff::mm_lincomb(fit, named)
  ## Adapt to the historical column shape this helper used to return
  ## (estimate, SE, z, p) so existing assertions don't need rewriting.
  data.frame(
    estimate = out$estimate,
    SE       = out$std_error,
    z        = out$statistic,
    p        = out$p_value
  )
}

aphantasia_diagnostic_payloads <- function(fit) {
  cert <- fit$artifact$optimizer_certificate %||% list()
  diags <- c(fit$artifact$diagnostics %||% list(), cert$diagnostics %||% list())
  lapply(diags, function(diagnostic) diagnostic$payload %||% list())
}

test_that("aphantasia fixture has anonymized data and frozen references", {
  ref <- aphantasia_reference()
  trials <- aphantasia_trials()
  metadata <- aphantasia_metadata()

  expect_identical(ref$schema$name, "mixeff.aphantasia_fixture_reference")
  expect_equal(nrow(trials), ref$counts$trials)
  expect_equal(nrow(metadata), ref$counts$metadata_rows)
  expect_true(all(grepl("^p_[0-9a-f]{16}$", trials$participant)))
  expect_true(all(grepl("^p_[0-9a-f]{16}$", metadata$participant)))
  expect_false(any(grepl("^[0-9a-f]{24}$", trials$participant)))
  expect_true(all(c("primary", "sensitivity", "intact", "combined", "rt") %in%
                    names(ref$models)))
  expect_equal(length(aphantasia_s1_case_ids(ref)), 5L)
  expect_equal(ref$counts$primary_trials, 17280L)
  expect_equal(ref$counts$combined_trials, 23040L)
})

test_that("aphantasia core fit-side reproduction matches cached lme4 references when enabled", {
  testthat::skip_on_cran()
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    aphantasia_run_full(),
    "Set MIXEFF_RUN_APHANTASIA=true to run the core aphantasia reproduction."
  )

  ref <- aphantasia_reference()
  data_sets <- aphantasia_data_sets(ref, aphantasia_trials())
  cases <- aphantasia_fit_cases(ref, data_sets)

  for (id in aphantasia_core_case_ids(ref)) {
    model_ref <- ref$models[[id]]
    form <- stats::as.formula(model_ref$formula)
    fit <- if (identical(model_ref$model_type, "lmm")) {
      mixeff::lmm(form, cases[[id]]$data, REML = FALSE,
                  control = mixeff::mm_control(verbose = -1))
    } else {
      mixeff::glmm(form, cases[[id]]$data, family = cases[[id]]$family,
                   method = aphantasia_glmm_method(id), nAGQ = 1L,
                   control = aphantasia_glmm_control(id))
    }
    aphantasia_expect_fit_matches_reference(fit, model_ref, id)
  }
})

test_that("aphantasia intact budgeted joint Laplace proof improves the profiled gap when enabled", {
  testthat::skip_on_cran()
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    aphantasia_run_joint_proof(),
    paste(
      "Set MIXEFF_APHANTASIA_JOINT_PROOF=true to run the budgeted intact",
      "joint-Laplace timing/parity proof."
    )
  )

  budget <- aphantasia_joint_budget()
  ref <- aphantasia_reference()
  data_sets <- aphantasia_data_sets(ref, aphantasia_trials())
  cases <- aphantasia_fit_cases(ref, data_sets)
  model_ref <- ref$models$intact
  fit <- mixeff::glmm(
    stats::as.formula(model_ref$formula),
    cases$intact$data,
    family = cases$intact$family,
    method = "joint_laplace",
    nAGQ = 1L,
    control = mixeff::mm_control(verbose = -1, max_feval = budget)
  )

  observed <- mixeff::fixef(fit)
  names(observed) <- aphantasia_lme4_key(names(observed))
  expected <- unlist(model_ref$fixef, use.names = TRUE)
  common <- intersect(names(expected), names(observed))
  max_fixef_drift <- max(abs(observed[common] - expected[common]))
  loglik_gap <- abs(as.numeric(stats::logLik(fit)) - model_ref$logLik)
  payloads <- aphantasia_diagnostic_payloads(fit)
  fit_modes <- unlist(lapply(payloads, function(payload) {
    payload$fit_mode %||% NA_character_
  }), use.names = FALSE)
  scorecard_classes <- unlist(lapply(payloads, function(payload) {
    payload$scorecard_class %||% NA_character_
  }), use.names = FALSE)

  expect_lte(fit$fit$optimizer$function_evaluations, budget)
  expect_match(fit$fit$optimizer$return_value, "^JOINT_LAPLACE:")
  expect_lt(loglik_gap, 0.5)
  expect_lt(max_fixef_drift, 0.1)
  expect_true("uncertified_joint_candidate" %in% fit_modes)
  expect_true("budget_limited_joint_candidate" %in% scorecard_classes)
})

test_that("aphantasia S1 random-effects stability fits run in the stress tier", {
  testthat::skip_on_cran()
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    aphantasia_run_stress(),
    "Set MIXEFF_RUN_APHANTASIA_STRESS=true to run S1 random-effects stability fits."
  )

  ref <- aphantasia_reference()
  data_sets <- aphantasia_data_sets(ref, aphantasia_trials())
  cases <- aphantasia_fit_cases(ref, data_sets)

  for (id in aphantasia_s1_case_ids(ref)) {
    model_ref <- ref$models[[id]]
    fit <- mixeff::glmm(
      stats::as.formula(model_ref$formula),
      cases[[id]]$data,
      family = cases[[id]]$family,
      control = mixeff::mm_control(verbose = -1)
    )
    aphantasia_expect_fit_matches_reference(fit, model_ref, id)
  }
})

test_that("aphantasia GLMM inference checks are gated on full vcov support", {
  testthat::skip_on_cran()
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    aphantasia_run_full(),
    "Set MIXEFF_RUN_APHANTASIA=true to run the core aphantasia reproduction."
  )

  ref <- aphantasia_reference()
  data_sets <- aphantasia_data_sets(ref, aphantasia_trials())
  primary_ref <- ref$models$primary
  fit <- mixeff::glmm(
    stats::as.formula(primary_ref$formula),
    data_sets$primary,
    family = stats::binomial(),
    control = mixeff::mm_control(verbose = -1)
  )
  testthat::skip_if_not(
    aphantasia_has_glmm_full_vcov(fit),
    "GLMM fixed-effect covariance payload is not yet available."
  )

  s25 <- (log(0.025) - mean(data_sets$primary$soa_log)) /
    stats::sd(data_sets$primary$soa_log)
  observed <- rbind(
    cbind(where = "centered_soa",
          aphantasia_lincomb(fit, c("groupaphant:maskmasked" = 1))),
    cbind(where = "25_ms",
          aphantasia_lincomb(
            fit,
            c("groupaphant:maskmasked" = 1,
              "groupaphant:maskmasked:soa_s" = s25)
          ))
  )
  expected <- aphantasia_reference_rows(ref$inference$primary_dd)
  ## Absolute tolerance: the DiD log-odds estimate drifts ~3-4% relative
  ## (~0.007-0.011 absolute) vs lme4 from PIRLS/optimizer convergence,
  ## so a relative 2.5e-2 over-fails. The qualitative result (sign,
  ## significance, CI excluding zero) is unchanged. See
  ## bd-01KS61JMB85G8DC3BXVSQW48MA.
  max_dd_abs <- max(abs(observed$estimate - unlist(expected$estimate)))
  expect_true(
    max_dd_abs < 2e-2,
    info = sprintf("DiD estimate drift vs lme4 too large: max_abs=%s",
                   signif(max_dd_abs, 4))
  )
  expect_equal(observed$where, unlist(expected$where))
})
