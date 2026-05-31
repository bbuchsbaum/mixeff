# Integration cases extracted from planning/bw_LME_tutorial.pdf.
# The source CSV is the tutorial's politeness_data.csv, normalized to LF line
# endings in tests/fixtures/bw_politeness_data.csv.

bw_politeness_path <- function() {
  candidates <- c(
    file.path("tests", "fixtures", "bw_politeness_data.csv"),
    testthat::test_path("..", "fixtures", "bw_politeness_data.csv")
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip("Bodo Winter politeness fixture is unavailable")
  }
  hit
}

bw_politeness_data <- function(complete = TRUE) {
  dat <- utils::read.csv(bw_politeness_path(), stringsAsFactors = TRUE)
  dat$scenario <- factor(dat$scenario)
  if (complete) {
    dat <- stats::na.omit(dat)
  }
  dat
}

bw_tutorial_cases <- function() {
  list(
    random_intercepts_attitude = list(
      formula = frequency ~ attitude + (1 | subject) + (1 | scenario),
      REML = TRUE,
      tolerance = list(fixef = 2e-5, scalar = 1e-3, fitted = 1e-3,
                       varcorr = 1e-2)
    ),
    random_intercepts_attitude_gender = list(
      formula = frequency ~ attitude + gender + (1 | subject) + (1 | scenario),
      REML = TRUE,
      tolerance = list(fixef = 1e-5, scalar = 1e-3, fitted = 1e-3,
                       varcorr = 1e-2)
    ),
    ml_null_gender = list(
      formula = frequency ~ gender + (1 | subject) + (1 | scenario),
      REML = FALSE,
      tolerance = list(fixef = 1e-5, scalar = 1e-3, fitted = 1e-3,
                       varcorr = 1e-2)
    ),
    ml_full_attitude_gender = list(
      formula = frequency ~ attitude + gender + (1 | subject) + (1 | scenario),
      REML = FALSE,
      tolerance = list(fixef = 2e-5, scalar = 1e-3, fitted = 2e-3,
                       varcorr = 1e-2)
    ),
    ml_interaction = list(
      formula = frequency ~ attitude * gender + (1 | subject) + (1 | scenario),
      REML = FALSE,
      tolerance = list(fixef = 5e-5, scalar = 1e-3, fitted = 5e-3,
                       varcorr = 1e-2)
    ),
    random_slopes_attitude = list(
      formula = frequency ~ attitude + gender +
        (1 + attitude | subject) + (1 + attitude | scenario),
      REML = FALSE,
      tolerance = list(fixef = 2e-3, scalar = 1e-3, fitted = 6e-3,
                       varcorr = 1e-2)
    ),
    random_slopes_null_gender = list(
      formula = frequency ~ gender +
        (1 + attitude | subject) + (1 + attitude | scenario),
      REML = FALSE,
      tolerance = list(fixef = 2e-2, scalar = 1e-3, fitted = 1e-2,
                       varcorr = 2e-2)
    )
  )
}

bw_label_key <- function(x) {
  gsub(": ", "", as.character(x), fixed = TRUE)
}

bw_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- bw_politeness_data()
  list(
    mixeff = mixeff::lmm(case$formula, dat, REML = case$REML,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = dat, REML = case$REML)
    ))
  )
}

bw_expect_fixef_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::fixef(fit)
  expected <- lme4::fixef(ref)
  names(observed) <- bw_label_key(names(observed))
  names(expected) <- bw_label_key(names(expected))
  common <- intersect(names(observed), names(expected))
  expect_equal(observed[common], expected[common], tolerance = tolerance,
               info = sprintf("fixed-effect parity failed for `%s`", label))
}

bw_expect_varcorr_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::VarCorr(fit)
  expected <- as.data.frame(lme4::VarCorr(ref))
  expected <- expected[is.na(expected$var2) & expected$grp != "Residual", ,
                       drop = FALSE]

  for (i in seq_len(nrow(expected))) {
    hit <- mm_lme4_group_key(observed$table$group) ==
      mm_lme4_group_key(expected$grp[[i]]) &
      bw_label_key(observed$table$name) == bw_label_key(expected$var1[[i]])
    expect_true(any(hit), info = sprintf(
      "variance component `%s/%s` missing for `%s`",
      expected$grp[[i]], expected$var1[[i]], label
    ))
    if (any(hit)) {
      field <- sprintf(
        "varcorr.std_dev.%s.%s",
        mm_lme4_group_key(expected$grp[[i]]),
        bw_label_key(expected$var1[[i]])
      )
      mm_assert_parity(
        observed$table$std_dev[hit][1L],
        expected$sdcor[[i]],
        label,
        field,
        tolerance,
        "variance component"
      )
    }
  }

  residual_sd <- as.data.frame(lme4::VarCorr(ref))$sdcor[
    as.data.frame(lme4::VarCorr(ref))$grp == "Residual"
  ][1L]
  mm_assert_parity(observed$residual_sd, residual_sd, label,
                   "varcorr.residual_sd", tolerance, "residual SD")
}

test_that("Bodo Winter politeness fixture preserves tutorial data shape", {
  raw <- bw_politeness_data(complete = FALSE)
  dat <- bw_politeness_data()

  expect_equal(nrow(raw), 84L)
  expect_equal(sum(is.na(raw$frequency)), 1L)
  expect_equal(nrow(dat), 83L)
  expect_equal(length(unique(dat$subject)), 6L)
  expect_equal(length(unique(dat$scenario)), 7L)
  expect_equal(table(dat$attitude, dat$gender)[["inf", "F"]], 21L)
})

test_that("Bodo Winter LME tutorial models match lme4 core outputs", {
  mm_skip_if_no_lme4()

  for (label in names(bw_tutorial_cases())) {
    case <- bw_tutorial_cases()[[label]]
    pair <- bw_fit_pair(case)
    tol <- case$tolerance

    bw_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
    mm_assert_parity(sigma(pair$mixeff), sigma(pair$lme4), label,
                     "sigma", tol$scalar, "sigma")
    mm_assert_parity(as.numeric(logLik(pair$mixeff)),
                     as.numeric(stats::logLik(pair$lme4)), label,
                     "logLik", tol$scalar, "logLik")
    mm_assert_parity(AIC(pair$mixeff), AIC(pair$lme4), label,
                     "AIC", tol$scalar, "AIC")
    mm_assert_parity(BIC(pair$mixeff), BIC(pair$lme4), label,
                     "BIC", tol$scalar, "BIC")
    mm_assert_parity(fitted(pair$mixeff), fitted(pair$lme4), label,
                     "fitted", tol$fitted, "fitted-value")
    mm_assert_parity(residuals(pair$mixeff), residuals(pair$lme4), label,
                     "residuals", tol$fitted, "residual")
    bw_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
  }
})

test_that("Bodo Winter likelihood-ratio examples agree with lme4", {
  mm_skip_if_no_lme4()
  dat <- bw_politeness_data()

  comparisons <- list(
    attitude = list(
      null = frequency ~ gender + (1 | subject) + (1 | scenario),
      full = frequency ~ attitude + gender + (1 | subject) + (1 | scenario)
    ),
    interaction = list(
      null = frequency ~ attitude + gender + (1 | subject) + (1 | scenario),
      full = frequency ~ attitude * gender + (1 | subject) + (1 | scenario)
    ),
    random_slope_attitude = list(
      null = frequency ~ gender +
        (1 + attitude | subject) + (1 + attitude | scenario),
      full = frequency ~ attitude + gender +
        (1 + attitude | subject) + (1 + attitude | scenario)
    )
  )

  for (label in names(comparisons)) {
    spec <- comparisons[[label]]
    mm_null <- mixeff::lmm(spec$null, dat, REML = FALSE,
                           control = mixeff::mm_control(verbose = -1))
    mm_full <- mixeff::lmm(spec$full, dat, REML = FALSE,
                           control = mixeff::mm_control(verbose = -1))
    lme4_null <- suppressMessages(suppressWarnings(
      lme4::lmer(spec$null, data = dat, REML = FALSE)
    ))
    lme4_full <- suppressMessages(suppressWarnings(
      lme4::lmer(spec$full, data = dat, REML = FALSE)
    ))

    mm_cmp <- mixeff::compare(mm_null, mm_full)$table
    lme4_cmp <- stats::anova(lme4_null, lme4_full)
    mm_assert_parity(mm_cmp$LRT[[2L]], lme4_cmp$Chisq[[2L]], label,
                     "lrt.statistic", 1e-3, "LRT statistic")
    mm_assert_parity(mm_cmp$p_value[[2L]], lme4_cmp$`Pr(>Chisq)`[[2L]],
                     label, "lrt.p_value", 1e-4, "LRT p-value")
  }
})
