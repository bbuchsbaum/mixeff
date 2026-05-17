# Integration cases extracted from:
# https://ourcodingclub.github.io/tutorials/mixed-models/
# Plotting, QQ plots, residual diagnostics, ggeffects, sjPlot, and stargazer
# examples are intentionally excluded; this file covers the lme4 model examples.

cc_dragons_path <- function() {
  candidates <- c(
    file.path("tests", "fixtures", "codingclub_dragons.RData"),
    testthat::test_path("..", "fixtures", "codingclub_dragons.RData")
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip("Coding Club dragons fixture is unavailable")
  }
  hit
}

cc_dragons_data <- function() {
  env <- new.env(parent = emptyenv())
  load(cc_dragons_path(), envir = env)
  dat <- env$dragons
  dat$bodyLength2 <- as.numeric(scale(dat$bodyLength, center = TRUE, scale = TRUE))
  dat$sample <- factor(dat$mountainRange:dat$site)
  dat
}

cc_tutorial_cases <- function() {
  list(
    first_random_intercept = list(
      formula = testScore ~ bodyLength2 + (1 | mountainRange),
      REML = TRUE,
      tolerance = list(fixef = 1e-5, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3)
    ),
    crossed_site_warning_example = list(
      formula = testScore ~ bodyLength2 + (1 | mountainRange) + (1 | site),
      REML = TRUE,
      tolerance = list(fixef = 5e-4, scalar = 1e-4, fitted = 1e-3,
                       varcorr = 1e-2)
    ),
    nested_sample_explicit = list(
      formula = testScore ~ bodyLength2 + (1 | mountainRange) + (1 | sample),
      REML = TRUE,
      tolerance = list(fixef = 2e-3, scalar = 1e-3, fitted = 5e-3,
                       varcorr = 3e-2)
    ),
    nested_slash_equivalent = list(
      formula = testScore ~ bodyLength2 + (1 | mountainRange / site),
      REML = TRUE,
      tolerance = list(fixef = 5e-4, scalar = 1e-4, fitted = 1e-3,
                       varcorr = 1e-2)
    ),
    nested_random_slope = list(
      formula = testScore ~ bodyLength2 +
        (1 + bodyLength2 | mountainRange / site),
      REML = TRUE,
      tolerance = list(fixef = 2e-3, scalar = 1e-4, fitted = 3e-3,
                       varcorr = 2e-2)
    ),
    ml_full_body_length = list(
      formula = testScore ~ bodyLength2 + (1 | mountainRange) + (1 | sample),
      REML = FALSE,
      tolerance = list(fixef = 1e-3, scalar = 1e-4, fitted = 2e-3,
                       varcorr = 1e-2)
    ),
    ml_reduced_intercept = list(
      formula = testScore ~ 1 + (1 | mountainRange) + (1 | sample),
      REML = FALSE,
      tolerance = list(fixef = 1e-8, scalar = 1e-4, fitted = 1e-3,
                       varcorr = 1e-2)
    )
  )
}

cc_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- cc_dragons_data()
  list(
    mixeff = mixeff::lmm(case$formula, dat, REML = case$REML,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = dat, REML = case$REML)
    ))
  )
}

cc_expect_fixef_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::fixef(fit)
  expected <- lme4::fixef(ref)
  expect_equal(observed[names(expected)], expected, tolerance = tolerance,
               info = sprintf("fixed-effect parity failed for `%s`", label))
}

cc_group_key <- function(x) {
  vapply(strsplit(mm_lme4_group_key(x), ":", fixed = TRUE), function(parts) {
    paste(sort(parts), collapse = ":")
  }, character(1))
}

cc_name_key <- function(x) {
  gsub(": ", "", as.character(x), fixed = TRUE)
}

cc_expect_varcorr_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::VarCorr(fit)
  expected <- as.data.frame(lme4::VarCorr(ref))
  expected <- expected[is.na(expected$var2) & expected$grp != "Residual", ,
                       drop = FALSE]

  for (i in seq_len(nrow(expected))) {
    hit <- cc_group_key(observed$table$group) == cc_group_key(expected$grp[[i]]) &
      cc_name_key(observed$table$name) == cc_name_key(expected$var1[[i]])
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

test_that("Coding Club dragons fixture preserves tutorial data shape", {
  dat <- cc_dragons_data()

  expect_equal(nrow(dat), 480L)
  expect_equal(length(unique(dat$mountainRange)), 8L)
  expect_equal(length(unique(dat$site)), 3L)
  expect_equal(length(unique(dat$sample)), 24L)
  expect_true(all(table(dat$mountainRange, dat$site) == 20L))
  expect_equal(mean(dat$bodyLength2), 0, tolerance = 1e-12)
})

test_that("Coding Club lme4 tutorial models match core lme4 outputs", {
  mm_skip_if_no_lme4()

  for (label in names(cc_tutorial_cases())) {
    case <- cc_tutorial_cases()[[label]]
    pair <- cc_fit_pair(case)
    tol <- case$tolerance

    cc_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
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
    cc_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
  }
})

test_that("Coding Club ML likelihood-ratio example agrees with lme4", {
  mm_skip_if_no_lme4()
  dat <- cc_dragons_data()

  full <- testScore ~ bodyLength2 + (1 | mountainRange) + (1 | sample)
  reduced <- testScore ~ 1 + (1 | mountainRange) + (1 | sample)

  mm_full <- mixeff::lmm(full, dat, REML = FALSE,
                         control = mixeff::mm_control(verbose = -1))
  mm_reduced <- mixeff::lmm(reduced, dat, REML = FALSE,
                            control = mixeff::mm_control(verbose = -1))
  lme4_full <- suppressMessages(suppressWarnings(
    lme4::lmer(full, data = dat, REML = FALSE)
  ))
  lme4_reduced <- suppressMessages(suppressWarnings(
    lme4::lmer(reduced, data = dat, REML = FALSE)
  ))

  mm_cmp <- mixeff::compare(mm_reduced, mm_full)$table
  lme4_cmp <- stats::anova(lme4_reduced, lme4_full)
  expect_equal(mm_cmp$LRT[[2L]], lme4_cmp$Chisq[[2L]], tolerance = 1e-4)
  expect_equal(mm_cmp$p_value[[2L]], lme4_cmp$`Pr(>Chisq)`[[2L]],
               tolerance = 1e-4)
})
