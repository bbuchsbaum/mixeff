# Integration cases extracted from the UBC Biology 501 mixed-effects workshop.
# Source page: https://www.zoology.ubc.ca/~bio501/R/workshops/lme.html
# Source data:
# - https://www.zoology.ubc.ca/~bio501/R/data/flycatcher.csv
# - https://www.zoology.ubc.ca/~bio501/R/data/goldfish.csv
# - https://www.zoology.ubc.ca/~bio501/R/data/kluane.csv

ubc_bio501_fixture_path <- function(name) {
  candidates <- c(
    file.path("tests", "fixtures", name),
    testthat::test_path("..", "fixtures", name)
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip(sprintf("UBC Biology 501 fixture `%s` is unavailable", name))
  }
  hit
}

ubc_bio501_data <- function(name) {
  dat <- utils::read.csv(ubc_bio501_fixture_path(name), stringsAsFactors = TRUE)
  if ("bird" %in% names(dat)) {
    dat$bird <- factor(dat$bird)
    dat$year <- factor(dat$year)
  }
  if ("plot" %in% names(dat)) {
    dat$plot <- factor(dat$plot)
    dat$treatment <- factor(dat$treatment,
                            levels = c("control", "exclosure", "fertilizer", "both"))
  }
  dat
}

ubc_bio501_cases <- function() {
  list(
    flycatcher_repeatability = list(
      file = "ubc_bio501_flycatcher.csv",
      formula = patch ~ 1 + (1 | bird),
      REML = TRUE,
      expected_status = "match",
      tolerance = list(fixef = 1e-5, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3, ranef = 1e-4)
    ),
    goldfish_singular_wavelength = list(
      file = "ubc_bio501_goldfish.csv",
      formula = sensitivity ~ wavelength + (1 | fish),
      REML = TRUE,
      expected_status = "known_boundary",
      tolerance = list(fixef = 1e-5, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3, ranef = 1e-4),
      boundary_group = "fish"
    ),
    kluane_yarrow_no_interaction = list(
      file = "ubc_bio501_kluane.csv",
      formula = log(phen.ach) ~ treatment + duration + (1 | plot),
      REML = TRUE,
      expected_status = "known_boundary",
      tolerance = list(fixef = 1e-5, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3, ranef = 1e-4),
      boundary_group = "plot"
    ),
    kluane_yarrow_interaction = list(
      file = "ubc_bio501_kluane.csv",
      formula = log(phen.ach) ~ treatment * duration + (1 | plot),
      REML = TRUE,
      expected_status = "match",
      tolerance = list(fixef = 1e-5, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3, ranef = 1e-4)
    )
  )
}

ubc_bio501_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- ubc_bio501_data(case$file)
  list(
    data = dat,
    mixeff = mixeff::lmm(case$formula, dat, REML = case$REML,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = dat, REML = case$REML)
    ))
  )
}

ubc_bio501_name_key <- function(x) {
  gsub(": ", "", as.character(x), fixed = TRUE)
}

ubc_bio501_expect_fixef_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::fixef(fit)
  expected <- lme4::fixef(ref)
  names(observed) <- ubc_bio501_name_key(names(observed))
  names(expected) <- ubc_bio501_name_key(names(expected))
  common <- intersect(names(observed), names(expected))
  expect_equal(observed[common], expected[common], tolerance = tolerance,
               info = sprintf("fixed-effect parity failed for `%s`", label))
}

ubc_bio501_expect_varcorr_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::VarCorr(fit)
  expected <- as.data.frame(lme4::VarCorr(ref))
  expected <- expected[is.na(expected$var2) & expected$grp != "Residual", ,
                       drop = FALSE]

  for (i in seq_len(nrow(expected))) {
    hit <- mm_lme4_group_key(observed$table$group) ==
      mm_lme4_group_key(expected$grp[[i]]) &
      ubc_bio501_name_key(observed$table$name) ==
      ubc_bio501_name_key(expected$var1[[i]])
    expect_true(any(hit), info = sprintf(
      "variance component `%s/%s` missing for `%s`",
      expected$grp[[i]], expected$var1[[i]], label
    ))
    if (any(hit)) {
      mm_assert_parity(
        observed$table$std_dev[hit][1L],
        expected$sdcor[[i]],
        label,
        sprintf("varcorr.std_dev.%s.%s",
                mm_lme4_group_key(expected$grp[[i]]),
                ubc_bio501_name_key(expected$var1[[i]])),
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

ubc_bio501_expect_ranef_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::ranef(fit)
  expected <- lme4::ranef(ref)
  common_groups <- intersect(mm_lme4_group_key(names(observed)),
                             mm_lme4_group_key(names(expected)))
  expect_true(length(common_groups) > 0L,
              info = sprintf("no common random-effect groups for `%s`", label))

  for (group in common_groups) {
    obs_group <- names(observed)[mm_lme4_group_key(names(observed)) == group][[1L]]
    exp_group <- names(expected)[mm_lme4_group_key(names(expected)) == group][[1L]]
    obs <- observed[[obs_group]]
    exp <- expected[[exp_group]]
    obs_rows <- mm_lme4_level_key(rownames(obs))
    exp_rows <- mm_lme4_level_key(rownames(exp))
    rows <- intersect(obs_rows, exp_rows)
    obs <- obs[match(rows, obs_rows), , drop = FALSE]
    exp <- exp[match(rows, exp_rows), , drop = FALSE]
    mm_assert_parity(as.matrix(obs), as.matrix(exp), label,
                     sprintf("ranef.%s", exp_group), tolerance,
                     "random effects")
  }
}

ubc_bio501_expect_boundary_state <- function(fit, ref, case, label) {
  expected_boundary <- identical(case$expected_status, "known_boundary")

  expect_identical(lme4::isSingular(ref), expected_boundary,
                   info = sprintf("lme4 singular status mismatch for `%s`", label))
  expect_identical(mixeff::is_singular(fit), expected_boundary,
                   info = sprintf("mixeff singular status mismatch for `%s`", label))

  vc <- mixeff::VarCorr(fit)
  expect_true("boundary" %in% names(vc$table),
              info = sprintf("VarCorr boundary column missing for `%s`", label))
  expect_identical(any(vc$table$boundary), expected_boundary,
                   info = sprintf("VarCorr boundary flag mismatch for `%s`", label))

  if (expected_boundary) {
    expected_group <- mm_lme4_group_key(case$boundary_group)
    boundary_groups <- mm_lme4_group_key(vc$table$group[vc$table$boundary])
    expect_true(expected_group %in% boundary_groups,
                info = sprintf("expected boundary group `%s` not flagged for `%s`",
                               case$boundary_group, label))
  }
}

test_that("UBC Biology 501 fixtures preserve workshop data shapes", {
  flycatcher <- ubc_bio501_data("ubc_bio501_flycatcher.csv")
  goldfish <- ubc_bio501_data("ubc_bio501_goldfish.csv")
  kluane <- ubc_bio501_data("ubc_bio501_kluane.csv")

  expect_equal(nrow(flycatcher), 60L)
  expect_equal(length(unique(flycatcher$bird)), 30L)
  expect_equal(as.numeric(table(flycatcher$year)), c(30, 30))

  expect_equal(nrow(goldfish), 45L)
  expect_equal(length(unique(goldfish$fish)), 5L)
  expect_equal(length(unique(goldfish$wavelength)), 9L)

  expect_equal(nrow(kluane), 32L)
  expect_equal(length(unique(kluane$plot)), 16L)
  expect_equal(levels(kluane$treatment),
               c("control", "exclosure", "fertilizer", "both"))
  expect_equal(as.numeric(table(kluane$duration)), c(16, 16))
})

test_that("UBC Biology 501 fixture cases declare boundary intent", {
  for (label in names(ubc_bio501_cases())) {
    case <- ubc_bio501_cases()[[label]]
    expect_true(case$expected_status %in% c("match", "known_boundary"),
                info = sprintf("unexpected expected_status for `%s`", label))
    if (identical(case$expected_status, "known_boundary")) {
      expect_false(is.null(case$boundary_group),
                   info = sprintf("known-boundary case `%s` needs boundary_group metadata",
                                  label))
    }
  }
})

test_that("UBC Biology 501 LMM workshop models match lme4 core outputs", {
  mm_skip_if_no_lme4()

  for (label in names(ubc_bio501_cases())) {
    case <- ubc_bio501_cases()[[label]]
    pair <- ubc_bio501_fit_pair(case)
    tol <- case$tolerance

    ubc_bio501_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
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
                     "fitted", tol$fitted, "fitted values")
    mm_assert_parity(residuals(pair$mixeff), residuals(pair$lme4), label,
                     "residuals", tol$fitted, "residuals")
    ubc_bio501_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
    ubc_bio501_expect_ranef_close(pair$mixeff, pair$lme4, tol$ranef, label)
    ubc_bio501_expect_boundary_state(pair$mixeff, pair$lme4, case, label)
  }
})

test_that("UBC flycatcher repeatability matches workshop calculation", {
  pair <- ubc_bio501_fit_pair(ubc_bio501_cases()$flycatcher_repeatability)
  vc <- as.data.frame(lme4::VarCorr(pair$lme4))$vcov
  repeatability <- vc[[1L]] / (vc[[1L]] + vc[[2L]])

  expect_equal(repeatability, 0.7764331, tolerance = 1e-6)
  expect_equal(
    pair$mixeff$varcorr$table$variance[[1L]] /
      (pair$mixeff$varcorr$table$variance[[1L]] + pair$mixeff$sigma^2),
    repeatability,
    tolerance = 1e-4
  )
})
