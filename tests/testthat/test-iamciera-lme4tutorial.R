# Integration cases extracted from https://github.com/iamciera/lme4tutorial.
# Plotting, QQ plots, modelcheck(), languageR, lmerTest, car, and confidence
# interval examples are intentionally excluded; this file covers the lme4 LMM
# fitting and model-comparison examples.

iamciera_stomata_path <- function() {
  candidates <- c(
    file.path("tests", "fixtures", "iamciera_modeling_example.txt"),
    testthat::test_path("..", "fixtures", "iamciera_modeling_example.txt")
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip("iamciera lme4tutorial stomata fixture is unavailable")
  }
  hit
}

iamciera_stomata_data <- function() {
  dat <- utils::read.delim(iamciera_stomata_path(), stringsAsFactors = TRUE)
  dat$trans_abs_stom <- sqrt(dat$abs_stom)
  dat$trans_epi_count <- sqrt(dat$epi_count)
  dat
}

iamciera_cases <- function() {
  list(
    max_model = list(
      formula = trans_abs_stom ~ il + (1 | tray) + (1 | row) + (1 | col),
      slow = TRUE,
      tolerance = list(fixef = 1e-4, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3)
    ),
    drop_col = list(
      formula = trans_abs_stom ~ il + (1 | tray) + (1 | row),
      slow = TRUE,
      tolerance = list(fixef = 1e-4, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3)
    ),
    drop_row = list(
      formula = trans_abs_stom ~ il + (1 | tray) + (1 | col),
      slow = TRUE,
      tolerance = list(fixef = 1e-4, scalar = 1e-4, fitted = 2e-4,
                       varcorr = 3e-3)
    ),
    drop_tray = list(
      formula = trans_abs_stom ~ il + (1 | row) + (1 | col),
      slow = TRUE,
      tolerance = list(fixef = 1e-4, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3)
    ),
    drop_il = list(
      formula = trans_abs_stom ~ 1 + (1 | tray) + (1 | row) + (1 | col),
      tolerance = list(fixef = 1e-6, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3)
    )
  )
}

# Gate status (2026-07-09 baseline at engine pin 3b6ec69): all fits are fast
# (0.02-0.07s); the sole opt-in failure is the `drop_tray` LRT (63.75 vs
# lme4 63.02). Root cause: the reduced model's ML cold start stops
# non-stationary (certificate free_gradient_norm = 17.7) yet reports
# converged_interior; the SAME model under REML matches lme4 to 3e-9, and
# warm-start/cobyla/pattern_search recover lme4's ML optimum exactly. Filed
# upstream as mixeff-rs bd-01KX33ZEQHHE8CWV5Z1KA7EG5G; reproducer in
# planning/probes/MINIMAL_case2_iamciera_drop_tray.R. Ungate when it lands.
iamciera_run_slow_parity <- function() {
  identical(tolower(Sys.getenv("MIXEFF_RUN_SLOW_PARITY")), "true")
}

iamciera_label_key <- function(x) {
  x <- as.character(x)
  x <- gsub(": ", "", x, fixed = TRUE)
  x
}

iamciera_case_id <- function(label) {
  paste0("iamciera_", label)
}

iamciera_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- iamciera_stomata_data()
  list(
    mixeff = mixeff::lmm(case$formula, dat, REML = TRUE,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = dat, REML = TRUE)
    ))
  )
}

iamciera_expect_fixef_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::fixef(fit)
  expected <- lme4::fixef(ref)
  names(observed) <- iamciera_label_key(names(observed))
  names(expected) <- iamciera_label_key(names(expected))
  common <- intersect(names(observed), names(expected))
  mm_assert_parity(observed[common], expected[common],
                   case_id = iamciera_case_id(label),
                   field = "fixef",
                   tolerance = tolerance,
                   label = "fixed effects")
  expect_equal(length(common), length(expected),
               info = sprintf("not all fixed-effect labels aligned for `%s`", label))
}

iamciera_expect_varcorr_close <- function(fit, ref, tolerance, label) {
  case_id <- iamciera_case_id(label)
  observed <- mixeff::VarCorr(fit)
  expected_full <- as.data.frame(lme4::VarCorr(ref))
  expected <- expected_full[is.na(expected_full$var2) &
                              expected_full$grp != "Residual", ,
                            drop = FALSE]

  for (i in seq_len(nrow(expected))) {
    grp <- expected$grp[[i]]
    name <- expected$var1[[i]]
    hit <- mm_lme4_group_key(observed$table$group) == mm_lme4_group_key(grp) &
      observed$table$name == name
    expect_true(any(hit), info = sprintf(
      "variance component `%s/%s` missing for `%s`", grp, name, label
    ))
    if (any(hit)) {
      mm_assert_parity(
        observed$table$std_dev[hit][1L],
        expected$sdcor[[i]],
        case_id = case_id,
        field = sprintf("varcorr.std_dev.%s.%s", grp, name),
        tolerance = tolerance,
        label = sprintf("VarCorr std_dev[%s/%s]", grp, name)
      )
    }
  }

  residual_sd <- expected_full$sdcor[expected_full$grp == "Residual"][1L]
  mm_assert_parity(observed$residual_sd, residual_sd,
                   case_id = case_id,
                   field = "varcorr.residual_sd",
                   tolerance = tolerance,
                   label = "VarCorr residual SD")
}

test_that("iamciera stomata fixture preserves tutorial data shape", {
  dat <- iamciera_stomata_data()

  expect_equal(nrow(dat), 727L)
  expect_true(all(c("plant", "abs_stom", "epi_count", "il", "row", "tray", "col") %in%
                    names(dat)))
  expect_equal(length(unique(dat$il)), 75L)
  expect_equal(length(unique(dat$row)), 10L)
  expect_equal(length(unique(dat$tray)), 16L)
  expect_equal(length(unique(dat$col)), 5L)
  expect_equal(dat$trans_abs_stom, sqrt(dat$abs_stom))
})

test_that("iamciera lme4tutorial models match core lme4 outputs", {
  mm_skip_if_no_lme4()

  for (label in names(iamciera_cases())) {
    case <- iamciera_cases()[[label]]
    if (isTRUE(case$slow)) {
      next
    }
    pair <- iamciera_fit_pair(case)
    tol <- case$tolerance
    case_id <- iamciera_case_id(label)

    iamciera_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
    mm_assert_parity(sigma(pair$mixeff), sigma(pair$lme4),
                     case_id = case_id, field = "sigma",
                     tolerance = tol$scalar, label = "sigma")
    mm_assert_parity(as.numeric(logLik(pair$mixeff)),
                     as.numeric(stats::logLik(pair$lme4)),
                     case_id = case_id, field = "logLik",
                     tolerance = tol$scalar, label = "logLik")
    mm_assert_parity(AIC(pair$mixeff), AIC(pair$lme4),
                     case_id = case_id, field = "AIC",
                     tolerance = tol$scalar, label = "AIC")
    mm_assert_parity(BIC(pair$mixeff), BIC(pair$lme4),
                     case_id = case_id, field = "BIC",
                     tolerance = tol$scalar, label = "BIC")
    mm_assert_parity(fitted(pair$mixeff), fitted(pair$lme4),
                     case_id = case_id, field = "fitted",
                     tolerance = tol$fitted, label = "fitted")
    mm_assert_parity(residuals(pair$mixeff), residuals(pair$lme4),
                     case_id = case_id, field = "residuals",
                     tolerance = tol$fitted, label = "residuals")
    iamciera_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
  }
})

test_that("slow iamciera lme4tutorial models match lme4 when enabled", {
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    iamciera_run_slow_parity(),
    "Set MIXEFF_RUN_SLOW_PARITY=true to run slow iamciera lme4tutorial parity cases."
  )

  for (label in names(iamciera_cases())) {
    case <- iamciera_cases()[[label]]
    if (!isTRUE(case$slow)) {
      next
    }
    pair <- iamciera_fit_pair(case)
    tol <- case$tolerance
    case_id <- iamciera_case_id(label)

    iamciera_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
    mm_assert_parity(sigma(pair$mixeff), sigma(pair$lme4),
                     case_id = case_id, field = "sigma",
                     tolerance = tol$scalar, label = "sigma")
    mm_assert_parity(as.numeric(logLik(pair$mixeff)),
                     as.numeric(stats::logLik(pair$lme4)),
                     case_id = case_id, field = "logLik",
                     tolerance = tol$scalar, label = "logLik")
    mm_assert_parity(AIC(pair$mixeff), AIC(pair$lme4),
                     case_id = case_id, field = "AIC",
                     tolerance = tol$scalar, label = "AIC")
    mm_assert_parity(BIC(pair$mixeff), BIC(pair$lme4),
                     case_id = case_id, field = "BIC",
                     tolerance = tol$scalar, label = "BIC")
    mm_assert_parity(fitted(pair$mixeff), fitted(pair$lme4),
                     case_id = case_id, field = "fitted",
                     tolerance = tol$fitted, label = "fitted")
    mm_assert_parity(residuals(pair$mixeff), residuals(pair$lme4),
                     case_id = case_id, field = "residuals",
                     tolerance = tol$fitted, label = "residuals")
    iamciera_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
  }
})

test_that("iamciera backwards-selection model comparisons agree with lme4", {
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    iamciera_run_slow_parity(),
    "Set MIXEFF_RUN_SLOW_PARITY=true to run slow iamciera model-comparison parity cases."
  )
  dat <- iamciera_stomata_data()
  cases <- iamciera_cases()
  comparisons <- list(
    drop_col = c("drop_col", "max_model"),
    drop_row = c("drop_row", "max_model"),
    drop_tray = c("drop_tray", "max_model"),
    drop_il = c("drop_il", "max_model")
  )

  mm_fits <- lapply(cases, function(case) {
    mixeff::lmm(case$formula, dat, REML = TRUE,
                control = mixeff::mm_control(verbose = -1))
  })
  lme4_fits <- lapply(cases, function(case) {
    suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = dat, REML = TRUE)
    ))
  })

  for (label in names(comparisons)) {
    pair <- comparisons[[label]]
    mm_cmp <- mixeff::compare(mm_fits[[pair[[1L]]]], mm_fits[[pair[[2L]]]])$table
    lme4_cmp <- stats::anova(lme4_fits[[pair[[1L]]]], lme4_fits[[pair[[2L]]]])

    expect_true(all(mm_cmp$refit), info = sprintf("comparison `%s` was not ML-refit", label))
    expect_equal(mm_cmp$LRT[[2L]], lme4_cmp$Chisq[[2L]], tolerance = 1e-4,
                 info = sprintf("LRT statistic mismatch for `%s`", label))
    expect_equal(mm_cmp$p_value[[2L]], lme4_cmp$`Pr(>Chisq)`[[2L]],
                 tolerance = 1e-4,
                 info = sprintf("LRT p-value mismatch for `%s`", label))
  }
})
