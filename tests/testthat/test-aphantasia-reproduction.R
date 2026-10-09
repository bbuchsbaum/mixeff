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

# Estimator routing for the core GLMM cases. The headline cases (primary,
# intact, combined) run glmm()'s default estimator, joint Laplace (glmer
# nAGQ = 1); the remaining sensitivity/subset variants (sensitivity, S7,
# S9 x2) run the profiled fast path, which also clears the strict case
# tolerances there and keeps real-data coverage of pirls_profiled at ~10x
# lower cost. Native `||` expands like lme4 (the factor term `mask` keeps
# its own correlated block: 6 theta, df identical to glmer), so both routes
# fit lme4's exact model family. Measured 2026-10-09 (engine b0037a4,
# release .so, single core) vs inst/extdata/aphantasia/reference.json:
#   case         route    max|dfixef|  logLik gap  AIC gap  df  time
#   primary      joint    5e-5         +0.0008     -0.002   15  ~119s
#   intact       joint    3.9e-3       -0.0022     +0.004   15  ~104s
#   combined     joint    1.27e-2      +0.137      -0.274   23  ~449s
#   sensitivity  profiled 1.27e-2      -0.145      +0.291   15  ~13s
#   S7           profiled 1.28e-2      -0.172      +0.344   16  ~10s
#   S9           profiled 1.32e-2      -0.185      +0.370   15  ~22s
#   S9 + age     profiled 1.57e-2      -0.194      +0.388   16  ~7s
#   rt (lmm)     --       <1e-5        <1e-10      <1e-10   12  <1s
# (Joint on the four variants: max|dfixef| <= 1.1e-4, logLik within 2e-3,
# but 82-203s per fit.) Profiled misses strict parity on intact
# (max|dfixef| 0.47) and combined (0.24), so those two need the joint route;
# combined's joint optimum sits 0.14 logLik ABOVE the glmer reference.
# MIXEFF_APHANTASIA_JOINT=false is the debug-build escape hatch (joint fits
# are many times slower under an unoptimized build): it routes every case
# through pirls_profiled, so expect strict-tolerance failures for intact
# and combined when it is engaged.
aphantasia_use_joint_glmm <- function() {
  !identical(tolower(Sys.getenv("MIXEFF_APHANTASIA_JOINT")), "false")
}

aphantasia_default_route_cases <- c("primary", "intact", "combined")

aphantasia_use_joint_for_case <- function(id) {
  id %in% aphantasia_default_route_cases && aphantasia_use_joint_glmm()
}

aphantasia_glmm_method <- function(id) {
  if (aphantasia_use_joint_for_case(id)) {
    return("joint_laplace")
  }
  "pirls_profiled"
}

aphantasia_glmm_control <- function(id) {
  mixeff::mm_control(verbose = -1)
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
  # Parameter count parity: lme4's df recovered from the frozen AIC/logLik.
  expect_identical(
    as.integer(attr(stats::logLik(fit), "df")),
    as.integer(round((aic_ref + 2 * loglik_ref) / 2)),
    info = sprintf("aphantasia `%s` df (parameter count) differs from lme4", id)
  )
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
  # Profiled fits carry `available_noninferential` since covariance schema
  # 1.1.0 (engine f82c646): geometry decodes, Wald certification stays with
  # the inference table. Joint-Laplace fits still say `available`.
  attr(V, "mm_status") %in% c("available", "available_noninferential") &&
    is.matrix(V) &&
    all(is.finite(V))
}

## mixeff coefficient names are lme4-identical (the coef-map contract), so
## lme4-style weight names drive mm_lincomb() directly — no key translation.
aphantasia_lincomb <- function(fit, weights) {
  stopifnot(all(names(weights) %in% names(mixeff::fixef(fit))))
  out <- mixeff::mm_lincomb(fit, weights)
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
    if (identical(id, "primary") && aphantasia_use_joint_for_case(id)) {
      # The manuscript's DiD point estimates on the default (joint Laplace)
      # path: a fixed linear combination of fixef, so it inherits the 5e-5
      # fixef parity (measured |dEst| < 1e-4). Wald SEs are exercised in the
      # working-Hessian block below.
      s25 <- (log(0.025) - mean(cases$primary$data$soa_log)) /
        stats::sd(cases$primary$data$soa_log)
      b <- mixeff::fixef(fit)
      did <- c(b[["groupaphant:maskmasked"]],
               b[["groupaphant:maskmasked"]] +
                 s25 * b[["groupaphant:maskmasked:soa_s"]])
      expected_dd <- aphantasia_reference_rows(ref$inference$primary_dd)
      mm_assert_parity(
        did, unlist(expected_dd$estimate),
        case_id = "aphantasia_primary",
        field = "inference.primary_dd.estimate_default_route",
        tolerance = 0.02,
        label = "aphantasia primary DiD estimate (default route)",
        mode = "absolute"
      )
      # Engine afc7c36: the active Hessian is not positive definite here
      # (participant mask block, corr ~0.98), so the Wald rows fall back to
      # RX conditional on theta -- glmer's vcov(use.hessian = FALSE) -- and
      # are labelled as such (reliability low, notes naming the Hessian).
      rows <- mixeff:::mm_glmm_coefficient_inference_rows(fit)
      expect_true(all(rows$status == "available"))
      expect_true(all(rows$reliability == "low"))
      expect_true(all(rows$reliability_reason ==
                        "glmm_laplace_rx_conditional_on_theta_wald"))
      expect_warning(stats::vcov(fit), class = "mm_vcov_rx_fallback")
      out <- paste(capture.output(print(suppressWarnings(summary(fit)))),
                   collapse = "\n")
      expect_match(out, "fixed-effect block RX", fixed = TRUE)
      # The joint stop is FTOL with an eager objective gap of ~4e-4 deviance
      # units (> the 1e-6 certification policy; logLik 0.001 above glmer):
      # presented as "convergence not certified", a warning, not a failure.
      if (identical(fit$fit_status, "not_optimized")) {
        conv <- mixeff:::mm_glmm_convergence_assessment(fit)
        expect_true(conv$not_certified_small_gap)
        expect_match(out, "convergence not certified", fixed = TRUE)
      }
    }
  }
})

test_that("native || fits the same model as lme4's explicit expansion on aphantasia data", {
  # lme4 expands (1 + mask + soa_s || participant) into
  # (1 | participant) + (0 + mask | participant) + (0 + soa_s | participant):
  # the factor term keeps its own correlated 2x2 block. The native `||` now
  # follows the same rule, so writing the expansion out by hand must fit the
  # identical model (same 6 theta, same df, same optimum). Profiled fits keep
  # this cheap (~10s each, measured 2026-10-09); the estimator does not
  # matter for an equivalence check of the random-effects structure.
  testthat::skip_on_cran()
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    aphantasia_run_full(),
    "Set MIXEFF_RUN_APHANTASIA=true to run the core aphantasia reproduction."
  )

  ref <- aphantasia_reference()
  data_sets <- aphantasia_data_sets(ref, aphantasia_trials())
  model_ref <- ref$models$primary
  native <- stats::as.formula(model_ref$formula)
  expanded <- correct ~ group * mask * soa_s + block +
    (1 | participant) + (0 + mask | participant) +
    (0 + soa_s | participant) + (1 | item)

  fit_native <- mixeff::glmm(native, data_sets$primary,
                             family = stats::binomial(),
                             method = "pirls_profiled",
                             control = mixeff::mm_control(verbose = -1))
  fit_expanded <- mixeff::glmm(expanded, data_sets$primary,
                               family = stats::binomial(),
                               method = "pirls_profiled",
                               control = mixeff::mm_control(verbose = -1))

  # lme4's parameter count: 9 fixef + 6 theta
  expect_identical(as.integer(attr(stats::logLik(fit_native), "df")), 15L)
  expect_identical(as.integer(attr(stats::logLik(fit_expanded), "df")), 15L)
  expect_equal(length(fit_native$theta), length(unlist(model_ref$theta)))
  expect_equal(as.numeric(stats::logLik(fit_native)),
               as.numeric(stats::logLik(fit_expanded)), tolerance = 1e-8)
  # Same model, slightly different optimizer path (term order differs):
  # measured fixef agreement ~3e-5 relative, logLik to 1e-8.
  expect_equal(mixeff::fixef(fit_native), mixeff::fixef(fit_expanded),
               tolerance = 1e-4)
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

  # Profiled fits clear the strict tolerances for every S1 variant and keep
  # this tier at ~100s total. Measured 2026-10-09 (max|dfixef|, logLik gap,
  # df = lme4's, time): intercept_only 1.2e-2 / -0.16 / 11 / 2s;
  # current_uncorrelated 1.3e-2 / -0.15 / 15 / 10s; correlated 1.3e-2 /
  # -0.13 / 16 / 18s; item_mask 1.6e-2 / -0.20 / 17 / 21s; maximal
  # 1.2e-2 / -0.14 / 22 / 48s. (The default joint route reaches
  # max|dfixef| <= 8e-5 and logLik within 1.3e-3 on all five, at 39-303s
  # per fit, ~10 min in total.)
  for (id in aphantasia_s1_case_ids(ref)) {
    model_ref <- ref$models[[id]]
    fit <- mixeff::glmm(
      stats::as.formula(model_ref$formula),
      cases[[id]]$data,
      family = cases[[id]]$family,
      method = "pirls_profiled",
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
  # inference = "working_hessian": this block DELIBERATELY exercises the
  # labelled working-Hessian approximation on the profiled estimator -- its
  # point is to pin how far the profiled DiD estimates and SEs sit from
  # glmer on real data (the evidence behind the opt-in's documented caveat).
  # Without the opt-in mm_lincomb() refuses SEs here: the profiled
  # covariance is uncertified. (On this dataset the default joint fit's
  # Hessian is not positive definite on the active parameter space -- the
  # fitted mask block sits near a boundary, corr ~0.98 -- so since engine
  # afc7c36 its Wald rows come from RX conditional on theta, graded low.)
  # ~10s per fit (measured 2026-10-09).
  fit <- mixeff::glmm(
    stats::as.formula(primary_ref$formula),
    data_sets$primary,
    family = stats::binomial(),
    method = "pirls_profiled",
    inference = "working_hessian",
    control = mixeff::mm_control(verbose = -1)
  )
  # Assert, don't skip: the payload is available on the current pin, and a
  # pin bump that regressed it would otherwise silently delete the DiD
  # estimate + SE parity assertions below from the opt-in run.
  testthat::expect_true(
    aphantasia_has_glmm_full_vcov(fit),
    info = "GLMM fixed-effect covariance payload regressed to unavailable"
  )
  if (!aphantasia_has_glmm_full_vcov(fit)) {
    return(invisible(NULL))
  }

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

  # With lme4-style `||` the profiled fit is in glmer's model family, and the
  # DiD rows sit inside the strict case tolerances (no ledger entry):
  # measured 2026-10-09 |dEst| = 0.0039 / 0.0083 (centered / 25 ms) against
  # tol 0.02, |dSE| = 8e-5 / 1.4e-4 (SE within 0.1% of glmer) against tol
  # 0.002 (~2.5% of the SE). Routing through mm_assert_parity records a
  # parity-scoreboard row per field.
  mm_assert_parity(
    observed$estimate, unlist(expected$estimate),
    case_id = "aphantasia_primary",
    field = "inference.primary_dd.estimate",
    tolerance = 0.02,
    label = "aphantasia primary DiD estimate",
    mode = "absolute"
  )
  mm_assert_parity(
    observed$SE, unlist(expected$SE),
    case_id = "aphantasia_primary",
    field = "inference.primary_dd.SE",
    tolerance = 0.002,
    label = "aphantasia primary DiD SE",
    mode = "absolute"
  )
  expect_equal(observed$where, unlist(expected$where))
})
