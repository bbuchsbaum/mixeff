# Integration case extracted from OSF project https://osf.io/h6rfv/:
# - Study 2 dataset.csv, https://osf.io/download/b2wyt/
# - Study 2 glmm codes.R, https://osf.io/download/tujz4/
# - Syntax (Study2) anova.sps, https://osf.io/download/psy3k/
#
# The supplied R script names the objects "glmm", but the fitted models are
# Gaussian lmer() fits. The raw CSV is wide; the test derives the long
# repeated-measures table implied by the SPSS 3 x 3 x 3 within-subject design.

osf_study2_fixture_path <- function() {
  candidates <- c(
    file.path("tests", "fixtures", "osf_study2_dataset.csv"),
    testthat::test_path("..", "fixtures", "osf_study2_dataset.csv")
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip("OSF Study 2 fixture is unavailable")
  }
  hit
}

osf_study2_wide_data <- function() {
  utils::read.csv(
    osf_study2_fixture_path(),
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM",
    check.names = TRUE
  )
}

osf_study2_prefix <- function(time, stim) {
  if (identical(time, "X")) {
    return(sprintf("X.%d", stim))
  }
  sprintf("%s_%d", time, stim)
}

osf_study2_long_data <- function() {
  wide <- osf_study2_wide_data()
  times <- c("pre", "post", "X")
  rows <- vector("list", nrow(wide) * 9L * length(times))
  row_i <- 0L

  for (sub_i in seq_len(nrow(wide))) {
    for (stim in seq_len(9L)) {
      attractiveness <- ((stim - 1L) %/% 3L) + 1L
      condition <- ((stim - 1L) %% 3L) + 1L
      for (time in times) {
        prefix <- osf_study2_prefix(time, stim)
        row_i <- row_i + 1L
        rows[[row_i]] <- data.frame(
          sub = wide$id[[sub_i]],
          stim = stim,
          con = condition,
          cute_level = attractiveness,
          time = time,
          cute = wide[[paste0(prefix, "_cute")]][[sub_i]],
          warmth = wide[[paste0(prefix, "_warmth")]][[sub_i]],
          competence = wide[[paste0(prefix, "_competence")]][[sub_i]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  out$sub <- factor(out$sub)
  out$stim <- factor(out$stim)
  out$con_f <- factor(out$con)
  out$time_f <- factor(out$time, levels = times)
  out$level_f <- factor(out$cute_level)
  out
}

osf_study2_cases <- function() {
  base <- list(
    formula_rhs = "~ con_f * time_f * level_f + (1 | stim) + (1 | sub) +
      (1 | sub:con_f) + (1 | sub:time_f)",
    REML = TRUE
  )
  list(
    cute_full_factorial = utils::modifyList(base, list(
      outcome = "cute",
      tolerance = list(fixef = 1e-7, scalar = 1e-3, fitted = 3e-3,
                       varcorr = 2e-3)
    )),
    warmth_full_factorial = utils::modifyList(base, list(
      outcome = "warmth",
      tolerance = list(fixef = 1e-7, scalar = 1e-4, fitted = 5e-4,
                       varcorr = 5e-4)
    )),
    competence_full_factorial = utils::modifyList(base, list(
      outcome = "competence",
      tolerance = list(fixef = 1e-7, scalar = 2.5e-1, fitted = 6e-2,
                       varcorr = 3e-2)
    ))
  )
}

osf_study2_case_id <- function(label) {
  paste0("osf_study2_", label)
}

osf_study2_formula <- function(case) {
  stats::as.formula(paste(case$outcome, case$formula_rhs))
}

osf_study2_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  dat <- osf_study2_long_data()
  formula <- osf_study2_formula(case)
  list(
    mixeff = mixeff::lmm(formula, dat, REML = case$REML,
                         control = mixeff::mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::lmer(formula, data = dat, REML = case$REML)
    ))
  )
}

osf_study2_fixef_key <- function(x) {
  gsub(": ", "", as.character(x), fixed = TRUE)
}

osf_study2_expect_fixef_close <- function(fit, ref, tolerance, label) {
  observed <- mixeff::fixef(fit)
  expected <- lme4::fixef(ref)
  names(observed) <- osf_study2_fixef_key(names(observed))
  names(expected) <- osf_study2_fixef_key(names(expected))
  common <- intersect(names(observed), names(expected))
  case_id <- osf_study2_case_id(label)

  mm_assert_parity(observed[common], expected[common],
                   case_id = case_id, field = "fixef",
                   tolerance = tolerance, label = "fixed effects",
                   mode = "absolute")
  expect_equal(length(common), length(expected),
               info = sprintf("not all fixed-effect labels aligned for `%s`",
                              label))
}

osf_study2_expect_varcorr_close <- function(fit, ref, tolerance, label) {
  case_id <- osf_study2_case_id(label)
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
        label = sprintf("VarCorr std_dev[%s/%s]", grp, name),
        mode = "absolute"
      )
    }
  }

  residual_sd <- expected_full$sdcor[expected_full$grp == "Residual"][1L]
  mm_assert_parity(observed$residual_sd, residual_sd,
                   case_id = case_id,
                   field = "varcorr.residual_sd",
                   tolerance = tolerance,
                   label = "VarCorr residual SD",
                   mode = "absolute")
}

osf_study2_expect_stim_design_note <- function(fit) {
  diag <- mixeff::diagnostics(fit)
  expect_true("design_weak_identifiability" %in% diag$table$code)
  expect_match(
    diag$table$message[diag$table$code == "design_weak_identifiability"][[1L]],
    "random intercept `stim` is aliased",
    fixed = TRUE
  )

  warnings <- mixeff::diagnostics(fit, severity = "warning")
  expect_true("design_weak_identifiability" %in% warnings$table$code)

  audit_text <- mixeff::audit(fit)$text
  expect_match(audit_text, "Wrapper design notes:", fixed = TRUE)
  expect_match(audit_text, "design_weak_identifiability", fixed = TRUE)

  explain_text <- mixeff::explain_model(fit)$text
  expect_match(explain_text, "Design notes:", fixed = TRUE)
  expect_match(explain_text, "design_weak_identifiability", fixed = TRUE)

  vc <- mixeff::VarCorr(fit)
  expect_true("note" %in% names(vc$table))
  stim_hit <- vc$table$group == "stim" & vc$table$name == "(Intercept)"
  expect_true(any(stim_hit))
  expect_match(
    vc$table$note[stim_hit][[1L]],
    "[design_weak_identifiability]",
    fixed = TRUE
  )

  summary_text <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(summary_text, "design_weak_identifiability", fixed = TRUE)
}


test_that("OSF Study 2 fixture preserves source data shape", {
  wide <- osf_study2_wide_data()
  long <- osf_study2_long_data()

  expect_equal(dim(wide), c(86L, 84L))
  expect_true(all(c("id", "study", "age", "pre_1_cute", "post_9_warmth",
                    "X.9_competence") %in% names(wide)))
  expect_equal(as.integer(table(wide$study)), c(44L, 42L))

  expect_equal(nrow(long), 86L * 9L * 3L)
  expect_equal(length(unique(long$sub)), 86L)
  expect_equal(length(unique(long$stim)), 9L)
  expect_equal(sort(levels(long$con_f)), c("1", "2", "3"))
  expect_equal(levels(long$time_f), c("pre", "post", "X"))
  expect_equal(sort(levels(long$level_f)), c("1", "2", "3"))
  expect_true(all(table(long$con_f, long$time_f, long$level_f) == 86L))
})

test_that("OSF Study 2 lmer examples match core lme4 outputs", {
  mm_skip_if_no_lme4()

  for (label in names(osf_study2_cases())) {
    case <- osf_study2_cases()[[label]]
    pair <- osf_study2_fit_pair(case)
    tol <- case$tolerance
    case_id <- osf_study2_case_id(label)

    osf_study2_expect_fixef_close(pair$mixeff, pair$lme4, tol$fixef, label)
    mm_assert_parity(sigma(pair$mixeff), sigma(pair$lme4),
                     case_id = case_id, field = "sigma",
                     tolerance = tol$scalar, label = "sigma",
                     mode = "absolute")
    mm_assert_parity(as.numeric(logLik(pair$mixeff)),
                     as.numeric(stats::logLik(pair$lme4)),
                     case_id = case_id, field = "logLik",
                     tolerance = tol$scalar, label = "logLik",
                     mode = "absolute")
    mm_assert_parity(AIC(pair$mixeff), AIC(pair$lme4),
                     case_id = case_id, field = "AIC",
                     tolerance = tol$scalar, label = "AIC",
                     mode = "absolute")
    mm_assert_parity(BIC(pair$mixeff), BIC(pair$lme4),
                     case_id = case_id, field = "BIC",
                     tolerance = tol$scalar, label = "BIC",
                     mode = "absolute")
    mm_assert_parity(fitted(pair$mixeff), fitted(pair$lme4),
                     case_id = case_id, field = "fitted",
                     tolerance = tol$fitted, label = "fitted",
                     mode = "absolute")
    mm_assert_parity(residuals(pair$mixeff), residuals(pair$lme4),
                     case_id = case_id, field = "residuals",
                     tolerance = tol$fitted, label = "residuals",
                     mode = "absolute")
    osf_study2_expect_varcorr_close(pair$mixeff, pair$lme4, tol$varcorr, label)
    if (identical(label, "cute_full_factorial")) {
      osf_study2_expect_stim_design_note(pair$mixeff)
    }
  }
})
