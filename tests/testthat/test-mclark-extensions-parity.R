# Parity tests mined from Michael Clark's "Mixed Models with R" extensions
# chapter (https://m-clark.github.io/mixed-models-with-R/extensions.html).
# Datasets are downloaded from the companion GitHub repo
# (github.com/m-clark/mixed-models-with-R/tree/master/data) and stored as
# tests/fixtures/mclark_*.RData.
#
# Structural patterns covered:
#   nurses  — nested hierarchical RE: (1|hospital) + (1|hospital:ward)
#             and slash-notation alias: (1|hospital/ward)
#   pupils  — truly crossed RE: (1|primary_school_id) + (1|secondary_school_id)
#             with continuous + ordered-factor fixed covariates
#   gpa     — longitudinal random intercept: (1|student) with numeric occasion
#   patents — Poisson GLMM: (1|year)  [profiled-PIRLS vs lme4's joint-Laplace;
#             small algorithmic offset expected, tolerances set accordingly]

# ---------------------------------------------------------------------------
# Dataset loaders
# ---------------------------------------------------------------------------

mclark_path <- function(name) {
  fname <- paste0("mclark_", name, ".RData")
  candidates <- c(
    file.path("tests", "fixtures", fname),
    testthat::test_path("..", "fixtures", fname)
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip(sprintf("m-clark fixture `%s` unavailable", fname))
  }
  hit
}

mclark_load <- function(name) {
  e <- new.env(parent = emptyenv())
  load(mclark_path(name), envir = e)
  get(ls(e)[1L], envir = e)
}

mclark_nurses <- function() {
  dat <- mclark_load("nurses")
  dat$hospital <- factor(dat$hospital)
  dat$ward     <- factor(dat$ward)
  dat$wardid   <- factor(dat$wardid)
  dat
}

mclark_pupils <- function() {
  dat <- mclark_load("pupils")
  dat$primary_school_id   <- factor(dat$primary_school_id)
  dat$secondary_school_id <- factor(dat$secondary_school_id)
  dat
}

mclark_gpa <- function() {
  mclark_load("gpa")   # student already factor; occasion numeric 0–5
}

mclark_patents <- function() {
  dat <- mclark_load("patents")
  dat$year <- factor(dat$year)
  dat
}

# ---------------------------------------------------------------------------
# LMM case list
# ---------------------------------------------------------------------------

mclark_lmm_cases <- function() {
  list(
    nurses_nested = list(
      data_fn  = mclark_nurses,
      formula  = stress ~ age + sex + experience + treatment + wardtype +
                   hospsize + (1 | hospital) + (1 | hospital:ward),
      REML     = TRUE,
      # Observed at authoring: logLik diff 3.7e-5, fixef max diff 6.8e-6
      tolerance = list(fixef = 1e-4, scalar = 1e-3, fitted = 1e-3,
                       varcorr = 5e-3)
    ),
    pupils_crossed = list(
      data_fn  = mclark_pupils,
      formula  = achievement ~ sex + ses +
                   (1 | primary_school_id) + (1 | secondary_school_id),
      REML     = TRUE,
      # Observed at authoring: logLik diff 3.2e-7, fixef max diff 1.7e-6
      tolerance = list(fixef = 1e-4, scalar = 1e-3, fitted = 1e-3,
                       varcorr = 5e-3)
    ),
    gpa_random_intercept = list(
      data_fn  = mclark_gpa,
      formula  = gpa ~ occasion + (1 | student),
      REML     = TRUE,
      # Observed at authoring: logLik diff 5.3e-8, fixef max diff 1.7e-13
      tolerance = list(fixef = 1e-5, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3)
    ),
    gpa_with_sex = list(
      data_fn  = mclark_gpa,
      formula  = gpa ~ occasion + sex + (1 | student),
      REML     = TRUE,
      # Cluster-level covariate; observed: fixef ~1e-13, logLik ~2e-9
      tolerance = list(fixef = 1e-5, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3)
    ),
    gpa_random_slopes = list(
      data_fn  = mclark_gpa,
      formula  = gpa ~ occasion + (1 + occasion | student),
      REML     = TRUE,
      # Correlated random intercept+slope; observed: fixef ~1e-13, logLik ~2e-9
      tolerance = list(fixef = 1e-5, scalar = 1e-4, fitted = 1e-4,
                       varcorr = 1e-3)
    ),
    pupils_re_slopes_sex = list(
      data_fn  = mclark_pupils,
      formula  = achievement ~ sex + ses + (1 + sex | primary_school_id),
      REML     = TRUE,
      # Binary categorical random slope; observed: fixef ~1e-6, fitted ~4e-5
      tolerance = list(fixef = 1e-4, scalar = 1e-3, fitted = 1e-3,
                       varcorr = 5e-3)
    ),
    pupils_interaction_re = list(
      data_fn  = mclark_pupils,
      formula  = achievement ~ sex + ses +
                   (1 | primary_school_id) + (1 | primary_school_id:sex),
      REML     = TRUE,
      # Approach 3 from Clark categorical section; observed: fixef ~3e-6
      tolerance = list(fixef = 1e-4, scalar = 1e-3, fitted = 1e-3,
                       varcorr = 5e-3)
    )
  )
}

# ---------------------------------------------------------------------------
# Fit pair helpers
# ---------------------------------------------------------------------------

mclark_lmm_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- case$data_fn()
  list(
    case    = case,
    data    = dat,
    mixeff  = mixeff::lmm(case$formula, dat, REML = case$REML,
                          control = mixeff::mm_control(verbose = -1)),
    lme4    = suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = dat, REML = case$REML)
    ))
  )
}

mclark_glmm_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- case$data_fn()
  list(
    case   = case,
    data   = dat,
    mixeff = mixeff::glmm(case$formula, dat, family = case$family,
                          control = mixeff::mm_control(verbose = -1)),
    lme4   = suppressMessages(suppressWarnings(
      lme4::glmer(case$formula, data = dat, family = case$family)
    ))
  )
}

# ---------------------------------------------------------------------------
# VarCorr parity helper (LMM only; no residual for Poisson GLMM)
# ---------------------------------------------------------------------------

# Normalize RE term names: lme4 concatenates factor variable+level ("sexfemale")
# while mixeff formats them as "sex: female".  Stripping ": " aligns them.
mm_re_term_key <- function(x) gsub(": *", "", x)

mclark_expect_varcorr_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::VarCorr(fit)
  expected <- as.data.frame(lme4::VarCorr(ref))
  diag_rows <- expected[is.na(expected$var2) & expected$grp != "Residual", ,
                        drop = FALSE]
  for (i in seq_len(nrow(diag_rows))) {
    grp  <- mm_lme4_group_key(diag_rows$grp[[i]])
    var1 <- diag_rows$var1[[i]]
    hit  <- mm_lme4_group_key(observed$table$group) == grp &
              mm_re_term_key(observed$table$name) == mm_re_term_key(var1)
    testthat::expect_true(
      any(hit),
      info = sprintf("VarCorr component `%s/%s` missing for `%s`",
                     diag_rows$grp[[i]], var1, label)
    )
    if (any(hit)) {
      testthat::expect_equal(
        observed$table$std_dev[hit][1L],
        diag_rows$sdcor[[i]],
        tolerance = tolerance,
        info = sprintf("VarCorr std_dev[%s/%s] parity failed for `%s`",
                       diag_rows$grp[[i]], var1, label)
      )
    }
  }
  residual_row <- expected[expected$grp == "Residual", , drop = FALSE]
  if (nrow(residual_row)) {
    testthat::expect_equal(
      observed$residual_sd,
      residual_row$sdcor[[1L]],
      tolerance = tolerance,
      info = sprintf("VarCorr residual_sd parity failed for `%s`", label)
    )
  }
}

# ---------------------------------------------------------------------------
# Shape/integrity tests
# ---------------------------------------------------------------------------

test_that("mclark nurses fixture has expected shape", {
  dat <- mclark_nurses()
  expect_equal(nrow(dat), 1000L)
  expect_equal(nlevels(dat$hospital), 25L)
  # ward is locally coded 1-4; wardid is globally unique
  expect_true(nlevels(dat$wardid) > nlevels(dat$ward))
  expect_true(all(c("stress", "age", "sex", "experience",
                    "treatment", "wardtype", "hospsize") %in% names(dat)))
})

test_that("mclark pupils fixture has expected shape", {
  dat <- mclark_pupils()
  expect_equal(nrow(dat), 1000L)
  expect_true(nlevels(dat$primary_school_id) > 1L)
  expect_true(nlevels(dat$secondary_school_id) > 1L)
  # confirm crossing: not all students attend the same secondary school
  tab <- table(dat$primary_school_id, dat$secondary_school_id)
  expect_true(sum(tab > 0L) > max(dim(tab)))
})

test_that("mclark gpa fixture has expected shape", {
  dat <- mclark_gpa()
  expect_equal(nrow(dat), 1200L)
  expect_equal(nlevels(dat$student), 200L)   # 200 students × 6 occasions
  expect_true(all(c("gpa", "occasion", "student") %in% names(dat)))
  expect_equal(sort(unique(dat$occasion)), 0:5)
})

test_that("mclark patents fixture has expected shape", {
  dat <- mclark_patents()
  expect_equal(nrow(dat), 4809L)
  expect_true(all(c("ncit", "opposition", "biopharm", "year") %in% names(dat)))
  # year is now a factor; sensible number of levels (late 1970s – mid 1990s)
  expect_true(nlevels(dat$year) >= 15L)
})

# ---------------------------------------------------------------------------
# LMM parity tests
# ---------------------------------------------------------------------------

test_that("mclark nurses (nested) matches lme4 core outputs", {
  mm_skip_if_no_lme4()
  case  <- mclark_lmm_cases()$nurses_nested
  pair  <- mclark_lmm_fit_pair(case)
  tol   <- case$tolerance
  label <- "nurses_nested"

  fit <- pair$mixeff
  ref <- pair$lme4

  fe_obs <- unname(mixeff::fixef(fit))
  fe_exp <- unname(lme4::fixef(ref))
  expect_equal(length(fe_obs), length(fe_exp),
               info = "nurses_nested: fixed-effect vector length mismatch")
  expect_equal(fe_obs, fe_exp, tolerance = tol$fixef,
               info = "nurses_nested: fixed-effect parity failed")

  mm_assert_parity(sigma(fit), sigma(ref), label, "sigma",
                   tol$scalar, "sigma")
  mm_assert_parity(as.numeric(logLik(fit)),
                   as.numeric(stats::logLik(ref)),
                   label, "logLik", tol$scalar, "logLik")
  mm_assert_parity(AIC(fit), AIC(ref), label, "AIC", tol$scalar, "AIC")
  mm_assert_parity(BIC(fit), BIC(ref), label, "BIC", tol$scalar, "BIC")
  mm_assert_parity(fitted(fit), fitted(ref), label, "fitted",
                   tol$fitted, "fitted-value")
  mm_assert_parity(residuals(fit), residuals(ref), label, "residuals",
                   tol$fitted, "residual")
  mclark_expect_varcorr_close(fit, ref, tol$varcorr, label)
})

test_that("mclark nurses slash notation (1|hospital/ward) equals explicit nesting", {
  mm_skip_if_no_lme4()
  dat <- mclark_nurses()

  form_explicit <- stress ~ age + sex + experience + treatment + wardtype +
                     hospsize + (1 | hospital) + (1 | hospital:ward)
  form_slash    <- stress ~ age + sex + experience + treatment + wardtype +
                     hospsize + (1 | hospital/ward)

  fit_exp   <- mixeff::lmm(form_explicit, dat, REML = TRUE,
                            control = mixeff::mm_control(verbose = -1))
  fit_slash <- mixeff::lmm(form_slash,    dat, REML = TRUE,
                            control = mixeff::mm_control(verbose = -1))

  expect_equal(unname(mixeff::fixef(fit_exp)),
               unname(mixeff::fixef(fit_slash)),
               tolerance = 1e-8,
               info = "slash and explicit nesting should give identical fixef")
  expect_equal(as.numeric(logLik(fit_exp)),
               as.numeric(logLik(fit_slash)),
               tolerance = 1e-8,
               info = "slash and explicit nesting should give identical logLik")
  expect_equal(sigma(fit_exp), sigma(fit_slash), tolerance = 1e-8,
               info = "slash and explicit nesting should give identical sigma")
})

test_that("mclark pupils (crossed RE) matches lme4 core outputs", {
  mm_skip_if_no_lme4()
  case  <- mclark_lmm_cases()$pupils_crossed
  pair  <- mclark_lmm_fit_pair(case)
  tol   <- case$tolerance
  label <- "pupils_crossed"

  fit <- pair$mixeff
  ref <- pair$lme4

  fe_obs <- unname(mixeff::fixef(fit))
  fe_exp <- unname(lme4::fixef(ref))
  expect_equal(length(fe_obs), length(fe_exp),
               info = "pupils_crossed: fixed-effect vector length mismatch")
  expect_equal(fe_obs, fe_exp, tolerance = tol$fixef,
               info = "pupils_crossed: fixed-effect parity failed")

  mm_assert_parity(sigma(fit), sigma(ref), label, "sigma",
                   tol$scalar, "sigma")
  mm_assert_parity(as.numeric(logLik(fit)),
                   as.numeric(stats::logLik(ref)),
                   label, "logLik", tol$scalar, "logLik")
  mm_assert_parity(fitted(fit), fitted(ref), label, "fitted",
                   tol$fitted, "fitted-value")
  mclark_expect_varcorr_close(fit, ref, tol$varcorr, label)
})

test_that("mclark GPA (longitudinal random intercept) matches lme4", {
  mm_skip_if_no_lme4()
  case  <- mclark_lmm_cases()$gpa_random_intercept
  pair  <- mclark_lmm_fit_pair(case)
  tol   <- case$tolerance
  label <- "gpa_random_intercept"

  fit <- pair$mixeff
  ref <- pair$lme4

  expect_equal(unname(mixeff::fixef(fit)), unname(lme4::fixef(ref)),
               tolerance = tol$fixef,
               info = "gpa_ri: fixed-effect parity failed")
  mm_assert_parity(sigma(fit), sigma(ref), label, "sigma",
                   tol$scalar, "sigma")
  mm_assert_parity(as.numeric(logLik(fit)),
                   as.numeric(stats::logLik(ref)),
                   label, "logLik", tol$scalar, "logLik")
  mm_assert_parity(fitted(fit), fitted(ref), label, "fitted",
                   tol$fitted, "fitted-value")
  mclark_expect_varcorr_close(fit, ref, tol$varcorr, label)
})

# ---------------------------------------------------------------------------
# GLMM parity — patents Poisson
# Observed at authoring (profiled-PIRLS vs lme4 joint-Laplace):
#   logLik diff 4.5e-4, fixef max diff 2.6e-3
# Tolerances sit above those margins; deviance is not compared because
# lme4::deviance.glmer reports the GLM deviance-residual scale, not -2*logLik.
# ---------------------------------------------------------------------------

test_that("mclark patents (Poisson GLMM, 1|year) matches lme4 within PIRLS tolerance", {
  mm_skip_if_no_lme4()

  dat <- mclark_patents()
  form <- ncit ~ opposition + biopharm + (1 | year)

  fit <- mixeff::glmm(form, dat, family = poisson(),
                      control = mixeff::mm_control(verbose = -1))
  ref <- suppressMessages(suppressWarnings(
    lme4::glmer(form, data = dat, family = poisson)
  ))

  fe_obs <- unname(mixeff::fixef(fit))
  fe_exp <- unname(lme4::fixef(ref))
  expect_equal(length(fe_obs), length(fe_exp),
               info = "patents: fixed-effect vector length mismatch")
  # 5e-3 sits comfortably above the observed 2.6e-3 profiled-vs-Laplace offset
  expect_lt(max(abs(fe_obs - fe_exp)), 5e-3,
            label = "patents: max|fixef diff| vs lme4")

  loglik_diff <- abs(as.numeric(logLik(fit)) - as.numeric(stats::logLik(ref)))
  expect_lt(loglik_diff, 5e-3,
            label = "patents: |logLik diff| vs lme4")

  # year random-effect SD: check lme4 VarCorr has a year component and compare
  vc_ref <- as.data.frame(lme4::VarCorr(ref))
  year_sd_ref <- vc_ref$sdcor[vc_ref$grp == "year" & is.na(vc_ref$var2)][1L]
  vc_fit <- mixeff::VarCorr(fit)
  year_hit <- mm_lme4_group_key(vc_fit$table$group) == "year"
  expect_true(any(year_hit), info = "patents: year RE component missing from VarCorr")
  if (any(year_hit)) {
    year_sd_fit <- vc_fit$table$std_dev[year_hit][1L]
    expect_lt(abs(year_sd_fit - year_sd_ref), 5e-2,
              label = "patents: year RE std_dev vs lme4")
  }
})

# ---------------------------------------------------------------------------
# LMM parity — random slopes chapter cases
# ---------------------------------------------------------------------------

test_that("mclark GPA with sex cluster-level covariate matches lme4", {
  mm_skip_if_no_lme4()
  case  <- mclark_lmm_cases()$gpa_with_sex
  pair  <- mclark_lmm_fit_pair(case)
  tol   <- case$tolerance
  label <- "gpa_with_sex"

  fit <- pair$mixeff
  ref <- pair$lme4

  expect_equal(unname(mixeff::fixef(fit)), unname(lme4::fixef(ref)),
               tolerance = tol$fixef,
               info = "gpa_with_sex: fixed-effect parity failed")
  mm_assert_parity(sigma(fit), sigma(ref), label, "sigma",
                   tol$scalar, "sigma")
  mm_assert_parity(as.numeric(logLik(fit)),
                   as.numeric(stats::logLik(ref)),
                   label, "logLik", tol$scalar, "logLik")
  mm_assert_parity(fitted(fit), fitted(ref), label, "fitted",
                   tol$fitted, "fitted-value")
  mclark_expect_varcorr_close(fit, ref, tol$varcorr, label)
})

test_that("mclark GPA random slopes (1 + occasion | student) matches lme4", {
  mm_skip_if_no_lme4()
  case  <- mclark_lmm_cases()$gpa_random_slopes
  pair  <- mclark_lmm_fit_pair(case)
  tol   <- case$tolerance
  label <- "gpa_random_slopes"

  fit <- pair$mixeff
  ref <- pair$lme4

  expect_equal(unname(mixeff::fixef(fit)), unname(lme4::fixef(ref)),
               tolerance = tol$fixef,
               info = "gpa_random_slopes: fixed-effect parity failed")
  mm_assert_parity(sigma(fit), sigma(ref), label, "sigma",
                   tol$scalar, "sigma")
  mm_assert_parity(as.numeric(logLik(fit)),
                   as.numeric(stats::logLik(ref)),
                   label, "logLik", tol$scalar, "logLik")
  mm_assert_parity(fitted(fit), fitted(ref), label, "fitted",
                   tol$fitted, "fitted-value")
  mclark_expect_varcorr_close(fit, ref, tol$varcorr, label)

  # Check the intercept-slope correlation. $table$correlation is stored as a
  # full-precision numeric (bd-01KTQHCZ0K), so parity holds at the varcorr
  # tolerance with no rounding slack.
  vc_fit <- mixeff::VarCorr(fit)
  vc_ref <- as.data.frame(lme4::VarCorr(ref))
  corr_ref <- vc_ref$sdcor[!is.na(vc_ref$var2) & vc_ref$grp == "student"][1L]
  corr_fit <- vc_fit$table$correlation[
    vc_fit$table$group == "student" & vc_fit$table$name == "occasion"
  ][1L]
  expect_true(is.numeric(corr_fit),
              info = "gpa_random_slopes: $table$correlation must be numeric")
  if (!is.na(corr_fit) && !is.na(corr_ref)) {
    expect_equal(corr_fit, corr_ref, tolerance = tol$varcorr,
                 info = "gpa_random_slopes: intercept-slope correlation parity failed")
  }
})

test_that("mclark pupils binary random slope (1 + sex | school) matches lme4", {
  mm_skip_if_no_lme4()
  case  <- mclark_lmm_cases()$pupils_re_slopes_sex
  pair  <- mclark_lmm_fit_pair(case)
  tol   <- case$tolerance
  label <- "pupils_re_slopes_sex"

  fit <- pair$mixeff
  ref <- pair$lme4

  fe_obs <- unname(mixeff::fixef(fit))
  fe_exp <- unname(lme4::fixef(ref))
  expect_equal(length(fe_obs), length(fe_exp),
               info = "pupils_re_slopes_sex: fixed-effect vector length mismatch")
  expect_equal(fe_obs, fe_exp, tolerance = tol$fixef,
               info = "pupils_re_slopes_sex: fixed-effect parity failed")
  mm_assert_parity(sigma(fit), sigma(ref), label, "sigma",
                   tol$scalar, "sigma")
  mm_assert_parity(as.numeric(logLik(fit)),
                   as.numeric(stats::logLik(ref)),
                   label, "logLik", tol$scalar, "logLik")
  mm_assert_parity(fitted(fit), fitted(ref), label, "fitted",
                   tol$fitted, "fitted-value")
  mclark_expect_varcorr_close(fit, ref, tol$varcorr, label)
})

test_that("mclark pupils interaction RE (1|school) + (1|school:sex) matches lme4", {
  mm_skip_if_no_lme4()
  case  <- mclark_lmm_cases()$pupils_interaction_re
  pair  <- mclark_lmm_fit_pair(case)
  tol   <- case$tolerance
  label <- "pupils_interaction_re"

  fit <- pair$mixeff
  ref <- pair$lme4

  fe_obs <- unname(mixeff::fixef(fit))
  fe_exp <- unname(lme4::fixef(ref))
  expect_equal(length(fe_obs), length(fe_exp),
               info = "pupils_interaction_re: fixed-effect vector length mismatch")
  expect_equal(fe_obs, fe_exp, tolerance = tol$fixef,
               info = "pupils_interaction_re: fixed-effect parity failed")
  mm_assert_parity(sigma(fit), sigma(ref), label, "sigma",
                   tol$scalar, "sigma")
  mm_assert_parity(as.numeric(logLik(fit)),
                   as.numeric(stats::logLik(ref)),
                   label, "logLik", tol$scalar, "logLik")
  mm_assert_parity(fitted(fit), fitted(ref), label, "fitted",
                   tol$fitted, "fitted-value")
  mclark_expect_varcorr_close(fit, ref, tol$varcorr, label)
})
