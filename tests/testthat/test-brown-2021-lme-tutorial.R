# Integration cases extracted from:
# Brown (2021), "An Introduction to Linear Mixed-Effects Modeling in R",
# planning/Brown 2021 - An Introduction to Linear Mixed-Effects Modeling in R.pdf
#
# Associated data and code are from OSF project https://osf.io/v6qag/:
# - /Data/rt_dummy_data.csv
# - /Data/rt_dummy_data_interaction.csv
# - /Data/figure_data.csv
# - /Code/intro_to_lmer.Rmd
#
# The paper also includes binomial glmer() examples using acc_dummy_data.csv.
# Those are intentionally not included here because mixeff::glmm() is currently
# a typed Phase 4 boundary, not a fitted GLMM bridge.

brown_fixture_path <- function(file) {
  candidates <- c(
    file.path("tests", "fixtures", file),
    testthat::test_path("..", "fixtures", file)
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip(sprintf("Brown 2021 fixture `%s` is unavailable", file))
  }
  hit
}

brown_rt_data <- function() {
  dat <- utils::read.csv(brown_fixture_path("brown_rt_dummy_data.csv"),
                         stringsAsFactors = FALSE)
  dat$PID <- factor(dat$PID)
  dat$stim <- factor(dat$stim)
  dat$modality <- ifelse(dat$modality == "Audio-only", 0, 1)
  dat
}

brown_rt_interaction_data <- function() {
  dat <- utils::read.csv(
    brown_fixture_path("brown_rt_dummy_data_interaction.csv"),
    stringsAsFactors = FALSE
  )
  dat$PID <- factor(dat$PID)
  dat$stim <- factor(dat$stim)
  dat$modality <- ifelse(dat$modality == "Audio-only", 0, 1)
  dat$SNR <- ifelse(dat$SNR == "Easy", 0, 1)
  dat
}

brown_figure_data <- function() {
  dat <- utils::read.csv(brown_fixture_path("brown_figure_data.csv"),
                         stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM")
  dat$PID <- factor(dat$PID)
  dat
}

brown_cases <- function() {
  list(
    rt_modality_full = list(
      data = brown_rt_data,
      formula = RT ~ 1 + modality + (1 + modality | PID) +
        (1 + modality | stim),
      REML = TRUE,
      slow = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa"),
      tolerance = list(fixef = 2e-4, scalar = 2e-3, fitted = 2e-3,
                       varcorr = 3e-1)
    ),
    rt_modality_reduced = list(
      data = brown_rt_data,
      formula = RT ~ 1 + (1 + modality | stim) + (1 + modality | PID),
      REML = TRUE,
      slow = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa"),
      tolerance = list(fixef = 1, scalar = 2e-3, fitted = 1,
                       varcorr = 3e-1)
    ),
    rt_interaction = list(
      data = brown_rt_interaction_data,
      formula = RT ~ 1 + modality + SNR + modality:SNR +
        (0 + modality | stim) + (1 | stim) +
        (1 + modality + SNR | PID),
      REML = TRUE,
      slow = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa"),
      tolerance = list(fixef = 5e-4, scalar = 5e-3, fitted = 5e-3,
                       varcorr = 2e-2)
    ),
    figure_random_intercepts = list(
      data = brown_figure_data,
      formula = yvar ~ 1 + xvar + (1 | PID),
      REML = TRUE,
      tolerance = list(fixef = 1e-3, scalar = 1e-2, fitted = 1e-2,
                       varcorr = 5e-2)
    ),
    figure_random_slopes = list(
      data = brown_figure_data,
      formula = yvar ~ 1 + xvar + (1 + xvar | PID),
      REML = TRUE,
      tolerance = list(fixef = 1e-3, scalar = 1e-2, fitted = 1e-2,
                       varcorr = 5e-2)
    )
  )
}

brown_run_slow_rt_parity <- function() {
  identical(tolower(Sys.getenv("MIXEFF_RUN_SLOW_PARITY")), "true")
}

brown_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- case$data()
  ref_args <- list(
    formula = case$formula,
    data = dat,
    REML = case$REML
  )
  if (!is.null(case$control)) {
    ref_args$control <- case$control
  }
  list(
    mixeff = mixeff::lmm(case$formula, dat, REML = case$REML,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      do.call(lme4::lmer, ref_args)
    ))
  )
}

brown_fixef_key <- function(x) {
  gsub(": ", "", as.character(x), fixed = TRUE)
}

brown_expect_fixef_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::fixef(fit)
  expected <- lme4::fixef(ref)
  names(observed) <- brown_fixef_key(names(observed))
  names(expected) <- brown_fixef_key(names(expected))
  common <- intersect(names(observed), names(expected))
  expect_equal(observed[common], expected[common], tolerance = tolerance,
               info = sprintf("fixed-effect parity failed for `%s`", label))
  expect_equal(length(common), length(expected),
               info = sprintf("not all fixed-effect labels aligned for `%s`",
                              label))
}

brown_group_key <- function(x) {
  sub("\\.[0-9]+$", "", mm_lme4_group_key(x))
}

brown_expect_varcorr_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::VarCorr(fit)
  expected <- as.data.frame(lme4::VarCorr(ref))
  expected <- expected[is.na(expected$var2) & expected$grp != "Residual", ,
                       drop = FALSE]

  for (i in seq_len(nrow(expected))) {
    hit <- brown_group_key(observed$table$group) ==
      brown_group_key(expected$grp[[i]]) &
      brown_fixef_key(observed$table$name) == brown_fixef_key(expected$var1[[i]])
    expect_true(any(hit), info = sprintf(
      "variance component `%s/%s` missing for `%s`",
      expected$grp[[i]], expected$var1[[i]], label
    ))
    if (any(hit)) {
      expect_equal(observed$table$std_dev[hit][1L], expected$sdcor[[i]],
                   tolerance = tolerance,
                   info = sprintf("variance component parity failed for `%s`",
                                  label))
    }
  }

  residual_sd <- as.data.frame(lme4::VarCorr(ref))$sdcor[
    as.data.frame(lme4::VarCorr(ref))$grp == "Residual"
  ][1L]
  expect_equal(observed$residual_sd, residual_sd, tolerance = tolerance,
               info = sprintf("residual SD parity failed for `%s`", label))
}

test_that("Brown 2021 OSF fixtures preserve tutorial data shapes", {
  rt <- brown_rt_data()
  rti <- brown_rt_interaction_data()
  fig <- brown_figure_data()

  expect_equal(nrow(rt), 21679L)
  expect_equal(length(unique(rt$PID)), 53L)
  expect_equal(length(unique(rt$stim)), 543L)
  expect_equal(sort(unique(rt$modality)), c(0, 1))

  expect_equal(nrow(rti), 21679L)
  expect_equal(length(unique(rti$PID)), 53L)
  expect_equal(length(unique(rti$stim)), 543L)
  expect_equal(sort(unique(rti$modality)), c(0, 1))
  expect_equal(sort(unique(rti$SNR)), c(0, 1))

  expect_equal(nrow(fig), 16L)
  expect_equal(length(unique(fig$PID)), 4L)
  expect_equal(length(unique(fig$xvar)), 16L)
})

test_that("Brown 2021 LMM tutorial examples match core lme4 outputs", {
  mm_skip_if_no_lme4()

  for (label in names(brown_cases())) {
    case <- brown_cases()[[label]]
    if (isTRUE(case$slow)) {
      next
    }
    pair <- brown_fit_pair(case)
    tol <- case$tolerance

    brown_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
    expect_equal(sigma(pair$mixeff), sigma(pair$lme4), tolerance = tol$scalar,
                 info = sprintf("sigma parity failed for `%s`", label))
    expect_equal(as.numeric(logLik(pair$mixeff)),
                 as.numeric(stats::logLik(pair$lme4)),
                 tolerance = tol$scalar,
                 info = sprintf("logLik parity failed for `%s`", label))
    expect_equal(AIC(pair$mixeff), AIC(pair$lme4), tolerance = tol$scalar,
                 info = sprintf("AIC parity failed for `%s`", label))
    expect_equal(BIC(pair$mixeff), BIC(pair$lme4), tolerance = tol$scalar,
                 info = sprintf("BIC parity failed for `%s`", label))
    expect_equal(fitted(pair$mixeff), fitted(pair$lme4),
                 tolerance = tol$fitted,
                 info = sprintf("fitted-value parity failed for `%s`", label))
    expect_equal(residuals(pair$mixeff), residuals(pair$lme4),
                 tolerance = tol$fitted,
                 info = sprintf("residual parity failed for `%s`", label))
    brown_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
  }
})

test_that("Brown 2021 large RT tutorial examples match lme4 when slow parity is enabled", {
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    brown_run_slow_rt_parity(),
    "Set MIXEFF_RUN_SLOW_PARITY=true to run large Brown 2021 RT mixed-model parity cases."
  )

  for (label in names(brown_cases())) {
    case <- brown_cases()[[label]]
    if (!isTRUE(case$slow)) {
      next
    }
    pair <- brown_fit_pair(case)
    tol <- case$tolerance

    brown_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
    expect_equal(sigma(pair$mixeff), sigma(pair$lme4), tolerance = tol$scalar,
                 info = sprintf("sigma parity failed for `%s`", label))
    expect_equal(as.numeric(logLik(pair$mixeff)),
                 as.numeric(stats::logLik(pair$lme4)),
                 tolerance = tol$scalar,
                 info = sprintf("logLik parity failed for `%s`", label))
    expect_equal(AIC(pair$mixeff), AIC(pair$lme4), tolerance = tol$scalar,
                 info = sprintf("AIC parity failed for `%s`", label))
    expect_equal(BIC(pair$mixeff), BIC(pair$lme4), tolerance = tol$scalar,
                 info = sprintf("BIC parity failed for `%s`", label))
    expect_equal(fitted(pair$mixeff), fitted(pair$lme4),
                 tolerance = tol$fitted,
                 info = sprintf("fitted-value parity failed for `%s`", label))
    expect_equal(residuals(pair$mixeff), residuals(pair$lme4),
                 tolerance = tol$fitted,
                 info = sprintf("residual parity failed for `%s`", label))
    brown_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
  }
})

test_that("Brown 2021 RT likelihood-ratio example agrees with lme4", {
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    brown_run_slow_rt_parity(),
    "Set MIXEFF_RUN_SLOW_PARITY=true to run the large Brown 2021 RT likelihood-ratio case."
  )
  dat <- brown_rt_data()
  full <- RT ~ 1 + modality + (1 + modality | PID) + (1 + modality | stim)
  reduced <- RT ~ 1 + (1 + modality | stim) + (1 + modality | PID)

  mm_full <- mixeff::lmm(full, dat, REML = FALSE,
                         control = mixeff::mm_control(verbose = -1))
  mm_reduced <- mixeff::lmm(reduced, dat, REML = FALSE,
                            control = mixeff::mm_control(verbose = -1))
  lme4_full <- suppressMessages(suppressWarnings(
    lme4::lmer(full, data = dat, REML = FALSE,
               control = lme4::lmerControl(optimizer = "bobyqa"))
  ))
  lme4_reduced <- suppressMessages(suppressWarnings(
    lme4::lmer(reduced, data = dat, REML = FALSE,
               control = lme4::lmerControl(optimizer = "bobyqa"))
  ))

  mm_cmp <- mixeff::compare(mm_reduced, mm_full)$table
  lme4_cmp <- stats::anova(lme4_reduced, lme4_full)
  expect_equal(mm_cmp$LRT[[2L]], lme4_cmp$Chisq[[2L]], tolerance = 5e-3)
  expect_equal(mm_cmp$p_value[[2L]], lme4_cmp$`Pr(>Chisq)`[[2L]],
               tolerance = 1e-4)
})
