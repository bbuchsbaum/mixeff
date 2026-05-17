# Parity test against the pw2 / Patrick Ward "Crossed vs Nested Random Effects"
# tutorial (https://github.com/pw2/crossed_vs_nested_effects). The dataset is
# a 60-row simulated sports-team example with 6 teams and 37 unique players,
# of whom 17 play on multiple teams during the season (max 4 teams per
# player). This produces the partially-crossed-and-nested topology that the
# tutorial uses to contrast `(1 | team) + (1 | player_id)` (crossed) with
# `(1 | team) + (1 | team:player_id)` (nested).
#
# The CSV at tests/fixtures/pw2_crossed_nested.csv is the deterministic
# output of the tutorial's seed = 225544 simulation. See the Rmd source in
# the repo for the data-generation code reproduced verbatim.

pw2_path <- function() {
  candidates <- c(
    file.path("tests", "fixtures", "pw2_crossed_nested.csv"),
    testthat::test_path("..", "fixtures", "pw2_crossed_nested.csv")
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip("pw2 crossed/nested fixture is unavailable")
  }
  hit
}

pw2_data <- function() {
  dat <- utils::read.csv(pw2_path(), stringsAsFactors = FALSE)
  dat$team <- factor(dat$team)
  dat$player_id <- factor(dat$player_id)
  dat
}

pw2_cases <- function() {
  list(
    crossed = list(
      formula = player_value ~ 1 + (1 | team) + (1 | player_id),
      REML    = TRUE,
      tolerance = list(fixef = 1e-4, scalar = 1e-3, fitted = 5e-3,
                       varcorr = 5e-2)
    ),
    nested = list(
      formula = player_value ~ 1 + (1 | team) + (1 | team:player_id),
      REML    = TRUE,
      tolerance = list(fixef = 1e-4, scalar = 1e-3, fitted = 5e-3,
                       varcorr = 5e-2)
    )
  )
}

pw2_label_key <- function(x) {
  gsub(": ", "", as.character(x), fixed = TRUE)
}

pw2_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- pw2_data()
  list(
    case = case,
    data = dat,
    mixeff = mixeff::lmm(case$formula, dat, REML = case$REML,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = dat, REML = case$REML)
    ))
  )
}

pw2_expect_fixef_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::fixef(fit)
  expected <- lme4::fixef(ref)
  names(observed) <- pw2_label_key(names(observed))
  names(expected) <- pw2_label_key(names(expected))
  common <- intersect(names(observed), names(expected))
  expect_equal(observed[common], expected[common], tolerance = tolerance,
               info = sprintf("fixed-effect parity failed for `%s`", label))
}

pw2_expect_varcorr_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::VarCorr(fit)
  expected <- as.data.frame(lme4::VarCorr(ref))
  expected <- expected[is.na(expected$var2) & expected$grp != "Residual", ,
                       drop = FALSE]
  for (i in seq_len(nrow(expected))) {
    hit <- mm_lme4_group_key(observed$table$group) ==
      mm_lme4_group_key(expected$grp[[i]]) &
      pw2_label_key(observed$table$name) == pw2_label_key(expected$var1[[i]])
    expect_true(any(hit), info = sprintf(
      "variance component `%s/%s` missing for `%s`",
      expected$grp[[i]], expected$var1[[i]], label
    ))
    if (any(hit)) {
      expect_equal(observed$table$std_dev[hit][1L], expected$sdcor[[i]],
                   tolerance = tolerance,
                   info = sprintf("variance-component parity failed for `%s/%s` (%s)",
                                  expected$grp[[i]], expected$var1[[i]], label))
    }
  }
  residual_sd <- as.data.frame(lme4::VarCorr(ref))$sdcor[
    as.data.frame(lme4::VarCorr(ref))$grp == "Residual"
  ][1L]
  expect_equal(observed$residual_sd, residual_sd, tolerance = tolerance,
               info = sprintf("residual SD parity failed for `%s`", label))
}

test_that("pw2 crossed/nested fixture has the expected shape", {
  dat <- pw2_data()
  expect_equal(nrow(dat), 60L)
  expect_equal(length(levels(dat$team)), 6L)
  expect_equal(length(levels(dat$player_id)), 37L)
  # At least one player plays on more than one team -- the crossed signal.
  team_counts <- vapply(split(dat$team, dat$player_id),
                       function(t) length(unique(t)), integer(1))
  expect_true(any(team_counts > 1L),
              info = "fixture should contain players on multiple teams")
})

test_that("pw2 crossed and nested random-intercept models match lme4", {
  mm_skip_if_no_lme4()
  for (label in names(pw2_cases())) {
    case <- pw2_cases()[[label]]
    pair <- pw2_fit_pair(case)
    tol <- case$tolerance

    pw2_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
    expect_equal(sigma(pair$mixeff), sigma(pair$lme4), tolerance = tol$scalar,
                 info = sprintf("sigma parity failed for `%s`", label))
    expect_equal(as.numeric(logLik(pair$mixeff)),
                 as.numeric(stats::logLik(pair$lme4)),
                 tolerance = tol$scalar,
                 info = sprintf("logLik parity failed for `%s`", label))
    expect_equal(AIC(pair$mixeff), AIC(pair$lme4), tolerance = tol$scalar,
                 info = sprintf("AIC parity failed for `%s`", label))
    expect_equal(fitted(pair$mixeff), fitted(pair$lme4),
                 tolerance = tol$fitted,
                 info = sprintf("fitted-value parity failed for `%s`", label))
    pw2_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
  }
})

test_that("pw2 nested form expands to the (1 | team:player_id) grouping in the artifact", {
  case <- pw2_cases()$nested
  fit <- mixeff::lmm(case$formula, pw2_data(), REML = case$REML,
                     control = mixeff::mm_control(verbose = -1))
  ranef_groups <- names(mixeff::ranef(fit))
  # The two grouping factors should be `team` and the interaction `team:player_id`.
  expect_true("team" %in% ranef_groups,
              info = "nested model should expose a `team` random effect")
  expect_true(any(grepl("team", ranef_groups) & grepl("player_id", ranef_groups)),
              info = "nested model should expose a `team:player_id` random effect")
})
