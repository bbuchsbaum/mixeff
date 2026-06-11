# Integration cases from Phillip M. Alday's RPubs example:
# https://rpubs.com/palday/mixed-interactions
#
# The example uses lme4::cake directly, so there is no external data file to
# vendor. The useful fixture is the ordered-factor interaction analysis path:
# angle ~ recipe * temperature + (1 | recipe:replicate).

palday_cake_data <- function() {
  mm_skip_if_no_lme4()
  env <- new.env(parent = emptyenv())
  utils::data("cake", package = "lme4", envir = env)
  if (!exists("cake", envir = env, inherits = FALSE)) {
    testthat::skip("lme4::cake is unavailable")
  }
  get("cake", envir = env, inherits = FALSE)
}

palday_formula <- function(interaction = TRUE) {
  formula <- if (interaction) {
    "angle ~ recipe * temperature + (1 | recipe:replicate)"
  } else {
    "angle ~ recipe + temperature + (1 | recipe:replicate)"
  }
  stats::as.formula(formula, env = parent.frame())
}

palday_fit_pair <- function(reml = FALSE, interaction = TRUE) {
  dat <- palday_cake_data()
  formula <- palday_formula(interaction)
  list(
    data = dat,
    mixeff = mixeff::lmm(formula, dat, REML = reml,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::lmer(formula, data = dat, REML = reml)
    ))
  )
}

palday_case_metadata <- function() {
  list(
    source = "https://rpubs.com/palday/mixed-interactions",
    dataset = "lme4::cake",
    formula = "angle ~ recipe * temperature + (1 | recipe:replicate)",
    ordered_factor_contrast_policy = "contrast_basis_divergence",
    contrast_note = paste(
      "lme4 uses contr.poly for ordered factors; mixeff currently uses",
      "treatment coding. Comparable tests avoid coefficient-name/value parity",
      "until bd-01KTVF5NFWB4ERGCXVNEP2FWCV is implemented."
    )
  )
}

palday_lmerTest_fit <- function(reml) {
  mm_skip_if_no_lmerTest()
  dat <- palday_cake_data()
  suppressMessages(suppressWarnings(
    lmerTest::lmer(palday_formula(TRUE), data = dat, REML = reml)
  ))
}

palday_lmerTest_anova_row <- function(ref, term, ddf) {
  row <- as.data.frame(stats::anova(ref, type = 3, ddf = ddf))
  row$term <- rownames(row)
  row <- row[row$term == term, , drop = FALSE]
  expect_equal(nrow(row), 1L,
               info = sprintf("lmerTest did not return term `%s`", term))
  row
}

palday_mixeff_term_row <- function(fit, term, method) {
  row <- stats::anova(fit, type = "III", method = method)$table
  row <- row[row$term == term, , drop = FALSE]
  expect_equal(nrow(row), 1L,
               info = sprintf("mixeff did not return term `%s`", term))
  row
}

palday_expect_term_parity <- function(fit, ref, term, method, ddf) {
  observed <- palday_mixeff_term_row(fit, term, method)
  expected <- palday_lmerTest_anova_row(ref, term, ddf)

  expect_identical(observed$status, "available")
  expect_equal(observed$num_df, expected$NumDF, tolerance = 1e-12)
  expect_equal(observed$den_df, expected$DenDF, tolerance = 1e-2)
  expect_equal(observed$statistic, expected$`F value`, tolerance = 1e-3)
  if (is.finite(expected$`Pr(>F)`) && expected$`Pr(>F)` < 1e-12) {
    expect_true(observed$p_value <= 1e-12)
  } else {
    expect_equal(observed$p_value, expected$`Pr(>F)`, tolerance = 1e-3)
  }
}

palday_drop_row <- function(drop_obj, term) {
  table <- drop_obj$table %||% as.data.frame(drop_obj)
  row <- table[table$dropped == term, , drop = FALSE]
  expect_equal(nrow(row), 1L,
               info = sprintf("drop1() did not return dropped term `%s`", term))
  row
}

test_that("Palday cake fixture uses the lme4 bundled split-plot data", {
  dat <- palday_cake_data()

  expect_equal(nrow(dat), 270L)
  expect_equal(ncol(dat), 5L)
  expect_equal(names(dat),
               c("replicate", "recipe", "temperature", "angle", "temp"))
  expect_true(is.ordered(dat$temperature))
  expect_equal(levels(dat$temperature), c("175", "185", "195", "205", "215", "225"))
  expect_equal(length(unique(dat$recipe)), 3L)
  expect_equal(nlevels(interaction(dat$recipe, dat$replicate)), 45L)
})

test_that("Palday fixture metadata records ordered-factor contrast policy", {
  metadata <- palday_case_metadata()

  expect_identical(metadata$ordered_factor_contrast_policy,
                   "contrast_basis_divergence")
  expect_match(metadata$contrast_note, "contr.poly", fixed = TRUE)
  expect_match(metadata$contrast_note, "treatment coding", fixed = TRUE)
})

test_that("Palday ML interaction model matches lme4 on comparable fit quantities", {
  pair <- palday_fit_pair(reml = FALSE, interaction = TRUE)
  tol <- list(sigma = 1e-4, fitted = 1e-4, ranef = 1e-4, varcorr = 1e-3)

  mm_assert_parity(sigma(pair$mixeff), sigma(pair$lme4), "palday_cake_ml",
                   "sigma", tol$sigma, "sigma")
  mm_assert_parity(stats::predict(pair$mixeff), stats::predict(pair$lme4),
                   "palday_cake_ml", "fitted", tol$fitted, "fitted values")
  mm_assert_parity(residuals(pair$mixeff), residuals(pair$lme4),
                   "palday_cake_ml", "residuals", tol$fitted, "residuals")
  mm_expect_varcorr_lme4_parity(pair$mixeff, pair$lme4, tol$varcorr,
                                list(id = "palday_cake_ml"))
  expect_true("recipe:replicate" %in% pair$mixeff$varcorr$table$group)
  expect_false("recipe & replicate" %in% pair$mixeff$varcorr$table$group)
  expect_true("recipe:replicate" %in% as.data.frame(pair$mixeff$varcorr)$grp)
  expect_true("recipe:replicate" %in%
                reporting_table(pair$mixeff, "random_effects")$group)
  mm_expect_ranef_lme4_parity(list(
    id = "cake_ordered_temperature_interaction_ml",
    dataset = "cake",
    package = "lme4",
    formula = "angle ~ recipe * temperature + (1 | recipe:replicate)",
    reml = FALSE,
    default_tolerances = list(ranef = tol$ranef)
  ))

  expect_false(lme4::isSingular(pair$lme4))
  expect_false(mixeff::is_singular(pair$mixeff))
})

test_that("Palday drop1 interaction LRT matches lme4", {
  dat <- palday_cake_data()
  fit <- mixeff::lmm(palday_formula(TRUE), dat, REML = FALSE,
                     control = mixeff::mm_control(verbose = -1))
  ref <- suppressMessages(suppressWarnings(
    lme4::lmer(palday_formula(TRUE), data = dat, REML = FALSE)
  ))
  mm_drop <- stats::drop1(fit, test = "Chisq")
  lme4_drop <- stats::drop1(
    ref,
    scope = c("recipe", "temperature", "recipe:temperature"),
    test = "Chisq"
  )
  observed <- palday_drop_row(mm_drop, "recipe:temperature")
  expected <- lme4_drop["recipe:temperature", , drop = FALSE]

  expect_s3_class(mm_drop, "mm_drop1")
  expect_equal(observed$df, expected$npar, tolerance = 1e-12)
  expect_equal(observed$AIC, expected$AIC, tolerance = 1e-4)
  expect_equal(observed$LRT, expected$LRT, tolerance = 1e-4)
  expect_equal(observed$p_value, expected$`Pr(Chi)`, tolerance = 1e-4)
})

test_that("Palday nested-model interaction comparison matches lme4", {
  full <- palday_fit_pair(reml = FALSE, interaction = TRUE)
  reduced <- palday_fit_pair(reml = FALSE, interaction = FALSE)

  observed <- mixeff::compare(reduced$mixeff, full$mixeff)$table
  expected <- stats::anova(reduced$lme4, full$lme4)

  expect_identical(observed$status[[2L]], "available")
  expect_identical(observed$comparison_class[[2L]], "nested_fixed_effects")
  expect_false(observed$requires_ml_refit[[2L]])
  expect_equal(observed$df[[2L]], expected$npar[[2L]], tolerance = 1e-12)
  expect_equal(observed$logLik[[2L]], expected$logLik[[2L]], tolerance = 1e-4)
  expect_equal(observed$AIC[[2L]], expected$AIC[[2L]], tolerance = 1e-4)
  expect_equal(observed$BIC[[2L]], expected$BIC[[2L]], tolerance = 1e-4)
  expect_equal(observed$LRT[[2L]], expected$Chisq[[2L]], tolerance = 1e-4)
  expect_equal(observed$p_value[[2L]], expected$`Pr(>Chisq)`[[2L]],
               tolerance = 1e-4)
})

test_that("Palday interaction finite-sample term tests match lmerTest", {
  mm_skip_if_no_lmerTest()
  mm_skip_if_no_pbkrtest()

  ml <- palday_fit_pair(reml = FALSE, interaction = TRUE)$mixeff
  reml <- palday_fit_pair(reml = TRUE, interaction = TRUE)$mixeff
  lmerTest_ml <- palday_lmerTest_fit(reml = FALSE)
  lmerTest_reml <- palday_lmerTest_fit(reml = TRUE)

  palday_expect_term_parity(
    ml, lmerTest_ml, "recipe:temperature", "satterthwaite", "Satterthwaite"
  )
  palday_expect_term_parity(
    reml, lmerTest_reml, "recipe:temperature", "kenward_roger", "Kenward-Roger"
  )
})

test_that("Palday ordered-factor coefficient caveat is explicit", {
  pair <- palday_fit_pair(reml = FALSE, interaction = TRUE)
  mm_names <- names(mixeff::fixef(pair$mixeff))
  lme4_names <- names(lme4::fixef(pair$lme4))

  expect_true(is.ordered(pair$data$temperature))
  expect_true(any(grepl("temperature: 185", mm_names, fixed = TRUE)))
  expect_true(any(grepl("temperature.L", lme4_names, fixed = TRUE)))
  expect_false(identical(mm_names, lme4_names))
})
