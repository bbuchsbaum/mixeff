# Integration cases extracted from:
# https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html
# Plotting, QQ plots, afex::mixed(), buildmer(), confidence intervals, and the
# reduced correlated speeddate example are intentionally excluded; this file
# covers lme4::lmer examples with stable mixeff::lmm() parity.

sdamr_fixture_path <- function(file) {
  candidates <- c(
    file.path("tests", "fixtures", file),
    testthat::test_path("..", "fixtures", file)
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip(sprintf("sdamr fixture `%s` is unavailable", file))
  }
  hit
}

sdamr_anchoring_data <- function() {
  dat <- utils::read.csv(sdamr_fixture_path("sdamr_anchoring.csv"),
                         stringsAsFactors = TRUE)
  dat$anchor <- factor(dat$anchor)
  contrasts(dat$anchor) <- c(1 / 2, -1 / 2)
  dat$referrer <- factor(dat$referrer)
  dat
}

sdamr_speeddate_data <- function(complete = TRUE) {
  dat <- utils::read.csv(sdamr_fixture_path("sdamr_speeddate_lmm.csv"),
                         stringsAsFactors = TRUE)
  dat$iid <- factor(dat$iid)
  dat$pid <- factor(dat$pid)
  if (complete) {
    dat <- stats::na.omit(dat)
  }
  dat
}

sdamr_cases <- function() {
  list(
    anchoring_random_intercept = list(
      data = sdamr_anchoring_data,
      formula = everest_feet ~ anchor + (1 | referrer),
      REML = TRUE,
      fixed_basis = "anchor_treatment",
      tolerance = list(fixef = 1e-3, scalar = 1e-2, fitted = 0.1,
                       varcorr = 0.1)
    ),
    anchoring_random_slope = list(
      data = sdamr_anchoring_data,
      formula = everest_feet ~ anchor + (1 + anchor | referrer),
      REML = TRUE,
      fixed_basis = "anchor_treatment",
      varcorr_basis = "anchor_treatment",
      tolerance = list(fixef = 1e-3, scalar = 2e-2, fitted = 2,
                       varcorr = 2e-2),
      varcorr_mode = "relative"
    ),
    anchoring_uncorrelated_numeric_slope = list(
      data = sdamr_anchoring_data,
      formula = everest_feet ~ anchor_contrast +
        (1 | referrer) + (0 + anchor_contrast | referrer),
      REML = TRUE,
      tolerance = list(fixef = 0.5, scalar = 0.2, fitted = 10,
                       varcorr = 1e-2),
      varcorr_mode = "relative"
    ),
    speeddate_maximal_crossed = list(
      data = sdamr_speeddate_data,
      formula = other_like ~ other_attr_c * other_intel_c +
        (1 + other_attr_c * other_intel_c | iid) +
        (1 + other_attr_c * other_intel_c | pid),
      REML = TRUE,
      slow = TRUE,
      tolerance = list(fixef = 1e-4, scalar = 1e-4, fitted = 2e-3,
                       varcorr = 1e-3)
    ),
    speeddate_uncorrelated_crossed = list(
      data = sdamr_speeddate_data,
      formula = other_like ~ other_attr_c * other_intel_c +
        (1 + other_attr_c + other_intel_c || iid) +
        (1 + other_attr_c + other_intel_c || pid),
      REML = TRUE,
      slow = TRUE,
      tolerance = list(fixef = 1e-4, scalar = 1e-4, fitted = 1e-3,
                       varcorr = 1e-3)
    )
  )
}

sdamr_run_slow_parity <- function() {
  identical(tolower(Sys.getenv("MIXEFF_RUN_SLOW_PARITY")), "true")
}

sdamr_fixef_key <- function(x) {
  x <- as.character(x)
  x <- gsub(": ", "", x, fixed = TRUE)
  x
}

sdamr_case_id <- function(label) {
  paste0("sdamr_", label)
}

sdamr_expected_fixef <- function(ref, basis = NULL) {
  expected <- lme4::fixef(ref)
  names(expected) <- sdamr_fixef_key(names(expected))

  if (identical(basis, "anchor_treatment") &&
      all(c("(Intercept)", "anchor1") %in% names(expected))) {
    expected <- c(
      "(Intercept)" = unname(expected[["(Intercept)"]] +
        0.5 * expected[["anchor1"]]),
      anchorlow = unname(-expected[["anchor1"]])
    )
  }

  expected
}

sdamr_group_key <- function(x) {
  sub("\\.[0-9]+$", "", mm_lme4_group_key(x))
}

sdamr_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- case$data()
  list(
    mixeff = mixeff::lmm(case$formula, dat, REML = case$REML,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = dat, REML = case$REML)
    ))
  )
}

sdamr_expect_fixef_close <- function(fit, ref, tolerance, label,
                                     basis = NULL) {
  observed <- mixeff::fixef(fit)
  expected <- sdamr_expected_fixef(ref, basis)
  names(observed) <- sdamr_fixef_key(names(observed))
  common <- intersect(names(observed), names(expected))
  mm_assert_parity(observed[common], expected[common],
                   case_id = sdamr_case_id(label),
                   field = "fixef",
                   tolerance = tolerance,
                   label = "fixed effects")
  expect_equal(length(common), length(expected),
               info = sprintf("not all fixed-effect labels aligned for `%s`", label))
}

sdamr_expected_varcorr <- function(ref, basis = NULL) {
  expected <- as.data.frame(lme4::VarCorr(ref))

  if (identical(basis, "anchor_treatment")) {
    referrer <- expected$grp == "referrer"
    intercept <- referrer & expected$var1 == "(Intercept)" & is.na(expected$var2)
    slope <- referrer & expected$var1 == "anchor1" & is.na(expected$var2)
    covariance <- referrer & expected$var1 == "(Intercept)" &
      !is.na(expected$var2) & expected$var2 == "anchor1"

    if (any(intercept) && any(slope) && any(covariance)) {
      high_variance <- expected$vcov[intercept][1L] +
        0.25 * expected$vcov[slope][1L] + expected$vcov[covariance][1L]
      expected$vcov[intercept] <- high_variance
      expected$sdcor[intercept] <- sqrt(high_variance)
      expected$var1[slope] <- "anchorlow"
      expected$sdcor[slope] <- sqrt(expected$vcov[slope][1L])
    }
  }

  expected
}

sdamr_expect_varcorr_close <- function(fit, ref, tolerance, label,
                                       basis = NULL,
                                       mode = c("absolute", "relative")) {
  mode <- match.arg(mode)
  case_id <- sdamr_case_id(label)
  observed <- mixeff::VarCorr(fit)
  expected_full <- sdamr_expected_varcorr(ref, basis)
  expected <- expected_full[is.na(expected_full$var2) &
                              expected_full$grp != "Residual", ,
                            drop = FALSE]

  for (i in seq_len(nrow(expected))) {
    grp <- expected$grp[[i]]
    name <- expected$var1[[i]]
    hit <- sdamr_group_key(observed$table$group) == sdamr_group_key(grp) &
      sdamr_fixef_key(observed$table$name) == sdamr_fixef_key(name)
    expect_true(any(hit), info = sprintf(
      "variance component `%s/%s` missing for `%s`", grp, name, label
    ))
    if (any(hit)) {
      mm_assert_parity(
        as.numeric(observed$table$std_dev[hit][1L]),
        as.numeric(expected$sdcor[[i]]),
        case_id = case_id,
        field = sprintf("varcorr.std_dev.%s.%s",
                        sdamr_group_key(grp), sdamr_fixef_key(name)),
        tolerance = tolerance,
        label = sprintf("VarCorr std_dev[%s/%s]", grp, name),
        mode = mode
      )
    }
  }

  residual_expected <- sdamr_expected_varcorr(ref)
  residual_sd <- residual_expected$sdcor[residual_expected$grp == "Residual"][1L]
  mm_assert_parity(
    as.numeric(observed$residual_sd),
    as.numeric(residual_sd),
    case_id = case_id,
    field = "varcorr.residual_sd",
    tolerance = tolerance,
    label = "VarCorr residual SD",
    mode = mode
  )
}

test_that("sdamr fixtures preserve tutorial data shapes", {
  anchoring <- sdamr_anchoring_data()
  speeddate_raw <- sdamr_speeddate_data(complete = FALSE)
  speeddate <- sdamr_speeddate_data()

  expect_equal(nrow(anchoring), 4632L)
  expect_equal(length(unique(anchoring$referrer)), 31L)
  expect_equal(sort(levels(anchoring$anchor)), c("high", "low"))
  expect_true(all(c("anchor_contrast", "everest_feet") %in% names(anchoring)))

  expect_equal(nrow(speeddate_raw), 1562L)
  expect_equal(nrow(speeddate), 1509L)
  expect_equal(length(unique(speeddate$iid)), 102L)
  expect_equal(length(unique(speeddate$pid)), 102L)
})

test_that("sdamr companion lme4 examples match core lme4 outputs", {
  mm_skip_if_no_lme4()

  for (label in names(sdamr_cases())) {
    case <- sdamr_cases()[[label]]
    if (isTRUE(case$slow)) {
      next
    }
    pair <- sdamr_fit_pair(case)
    tol <- case$tolerance
    case_id <- sdamr_case_id(label)

    sdamr_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label,
                             case$fixed_basis)
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
    varcorr_mode <- case$varcorr_mode
    if (is.null(varcorr_mode)) {
      varcorr_mode <- "absolute"
    }
    sdamr_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label,
                               case$varcorr_basis, mode = varcorr_mode)
  }
})

test_that("slow sdamr speeddate crossed models match lme4 when enabled", {
  mm_skip_if_no_lme4()
  testthat::skip_if_not(
    sdamr_run_slow_parity(),
    "Set MIXEFF_RUN_SLOW_PARITY=true to run slow SDAMR speeddate crossed-model parity cases."
  )

  for (label in names(sdamr_cases())) {
    case <- sdamr_cases()[[label]]
    if (!isTRUE(case$slow)) {
      next
    }
    pair <- sdamr_fit_pair(case)
    tol <- case$tolerance
    case_id <- sdamr_case_id(label)

    sdamr_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label,
                             case$fixed_basis)
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
    varcorr_mode <- case$varcorr_mode
    if (is.null(varcorr_mode)) {
      varcorr_mode <- "absolute"
    }
    sdamr_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label,
                               case$varcorr_basis, mode = varcorr_mode)
  }
})

test_that("sdamr companion random-slope LRT agrees with lme4", {
  mm_skip_if_no_lme4()
  dat <- sdamr_anchoring_data()
  correlated <- everest_feet ~ anchor + (1 + anchor | referrer)
  uncorrelated <- everest_feet ~ anchor_contrast +
    (1 | referrer) + (0 + anchor_contrast | referrer)

  mm_correlated <- mixeff::lmm(correlated, dat, REML = FALSE,
                               control = mixeff::mm_control(verbose = -1))
  mm_uncorrelated <- mixeff::lmm(uncorrelated, dat, REML = FALSE,
                                 control = mixeff::mm_control(verbose = -1))
  lme4_correlated <- suppressMessages(suppressWarnings(
    lme4::lmer(correlated, data = dat, REML = FALSE)
  ))
  lme4_uncorrelated <- suppressMessages(suppressWarnings(
    lme4::lmer(uncorrelated, data = dat, REML = FALSE)
  ))

  mm_lrt <- 2 * (as.numeric(logLik(mm_correlated)) -
    as.numeric(logLik(mm_uncorrelated)))
  lme4_cmp <- stats::anova(lme4_uncorrelated, lme4_correlated)
  expect_equal(mm_lrt, lme4_cmp$Chisq[[2L]], tolerance = 1e-1)
  expect_equal(stats::pchisq(mm_lrt, df = 1L, lower.tail = FALSE),
               lme4_cmp$`Pr(>Chisq)`[[2L]], tolerance = 1e-3)
})
