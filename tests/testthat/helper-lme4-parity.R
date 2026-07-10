mm_lme4_parity_manifest_path <- function() {
  candidates <- c(
    file.path("tests", "fixtures", "lme4_baseline_cases.json"),
    testthat::test_path("..", "fixtures", "lme4_baseline_cases.json")
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    testthat::skip("lme4 parity fixture manifest is unavailable")
  }
  hit
}

mm_lme4_parity_manifest <- function() {
  path <- mm_lme4_parity_manifest_path()
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

mm_lme4_parity_cases <- function(ids = NULL, model = "lmm") {
  manifest <- mm_lme4_parity_manifest()
  cases <- manifest$cases
  if (!is.null(ids)) {
    cases <- Filter(function(x) x$id %in% ids, cases)
  } else if (!is.null(model)) {
    cases <- Filter(function(x) identical(x$model %||% "lmm", model), cases)
  }
  lapply(cases, function(case) {
    case$default_tolerances <- manifest$default_tolerances
    case
  })
}

mm_skip_if_no_lme4 <- function() {
  testthat::skip_if_not_installed("lme4")
}

mm_skip_if_no_lmerTest <- function() {
  testthat::skip_if_not_installed("lmerTest")
}

mm_skip_if_no_pbkrtest <- function() {
  testthat::skip_if_not_installed("pbkrtest")
}

mm_reference_versions <- function() {
  refs <- c("lme4", "lmerTest", "pbkrtest")
  versions <- vapply(refs, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(utils::packageVersion(pkg))
    } else {
      NA_character_
    }
  }, character(1))
  versions
}

mm_lme4_case_data <- function(case) {
  env <- new.env(parent = emptyenv())
  utils::data(list = case$dataset, package = case$package %||% "lme4", envir = env)
  if (!exists(case$dataset, envir = env, inherits = FALSE)) {
    testthat::skip(sprintf("Dataset `%s` is unavailable", case$dataset))
  }
  get(case$dataset, envir = env, inherits = FALSE)
}

mm_lme4_case_formula <- function(case) {
  stats::as.formula(case$formula, env = parent.frame())
}

mm_lme4_case_tolerances <- function(case) {
  defaults <- case$default_tolerances %||% list()
  overrides <- case$tolerances %||% list()
  utils::modifyList(defaults, overrides)
}

mm_lme4_group_key <- function(x) {
  gsub("\\s+", "", gsub("&", ":", as.character(x)))
}

mm_lme4_level_key <- function(x) {
  gsub("_", ":", as.character(x), fixed = TRUE)
}

mm_lme4_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  data <- mm_lme4_case_data(case)
  formula <- mm_lme4_case_formula(case)
  reml <- isTRUE(case$reml)
  list(
    case = case,
    data = data,
    formula = formula,
    mixeff = lmm(formula, data, REML = reml, control = mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(lme4::lmer(formula, data = data, REML = reml)))
  )
}

mm_lmerTest_fit_pair <- function(case) {
  mm_skip_if_no_lmerTest()
  data <- mm_lme4_case_data(case)
  formula <- mm_lme4_case_formula(case)
  reml <- isTRUE(case$reml)
  list(
    case = case,
    data = data,
    formula = formula,
    mixeff = lmm(formula, data, REML = reml, control = mm_control(verbose = -1)),
    lmerTest = suppressMessages(suppressWarnings(
      lmerTest::lmer(formula, data = data, REML = reml)
    ))
  )
}

mm_numeric_payload <- function(x) {
  if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  } else if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  dims <- dim(x)
  out <- as.numeric(x)
  if (!is.null(dims)) {
    dim(out) <- dims
  }
  out
}

mm_expect_numeric_close <- function(observed, expected, tolerance, label, case) {
  testthat::expect_equal(
    mm_numeric_payload(observed),
    mm_numeric_payload(expected),
    tolerance = tolerance,
    info = sprintf(
      "%s parity failed for case `%s` using %s; tolerance=%s; reference versions: %s",
      label,
      case$id,
      case$formula,
      format(tolerance, scientific = TRUE),
      paste(names(mm_reference_versions()), mm_reference_versions(), sep = "=", collapse = ", ")
    )
  )
}

mm_lmerTest_p_value <- function(row) {
  p_col <- grep("^Pr\\(", names(row), value = TRUE)
  if (!length(p_col)) {
    p_col <- grep("^p\\.value$", names(row), value = TRUE)
  }
  row[[p_col[[1L]]]][[1L]]
}

mm_scalar_contrast_index <- function(ref) {
  if (length(lme4::fixef(ref)) > 1L) {
    2L
  } else {
    1L
  }
}

mm_scalar_contrast_vector <- function(ref, index = mm_scalar_contrast_index(ref)) {
  L <- rep(0, length(lme4::fixef(ref)))
  L[[index]] <- 1
  names(L) <- names(lme4::fixef(ref))
  L
}

mm_lmerTest_reference_contrast <- function(ref, L, rhs, method) {
  if (identical(method, "satterthwaite")) {
    row <- lmerTest::contest1D(ref, L, rhs = rhs, ddf = "Satterthwaite")
    p_value <- mm_lmerTest_p_value(row)
    return(data.frame(
      estimate = row[["Estimate"]][[1L]] - rhs,
      std_error = row[["Std. Error"]][[1L]],
      df = row[["df"]][[1L]],
      statistic = row[["t value"]][[1L]],
      p_value = p_value,
      statistic_name = "t",
      stringsAsFactors = FALSE
    ))
  }
  if (identical(method, "kenward_roger")) {
    mm_skip_if_no_pbkrtest()
    row <- lmerTest::contest1D(ref, L, rhs = rhs, ddf = "Kenward-Roger")
    p_value <- mm_lmerTest_p_value(row)
    return(data.frame(
      estimate = row[["Estimate"]][[1L]] - rhs,
      std_error = row[["Std. Error"]][[1L]],
      df = row[["df"]][[1L]],
      statistic = row[["t value"]][[1L]],
      p_value = p_value,
      statistic_name = "t",
      stringsAsFactors = FALSE
    ))
  }
  if (identical(method, "asymptotic")) {
    beta <- lme4::fixef(ref)
    estimate <- as.numeric(crossprod(L, beta)) - rhs
    std_error <- sqrt(as.numeric(t(L) %*% stats::vcov(ref) %*% L))
    statistic <- estimate / std_error
    return(data.frame(
      estimate = estimate,
      std_error = std_error,
      df = NA_real_,
      statistic = statistic,
      p_value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
      statistic_name = "z",
      stringsAsFactors = FALSE
    ))
  }
  stop(sprintf("Unsupported scalar contrast reference method `%s`.", method),
       call. = FALSE)
}

mm_lmerTest_scalar_rhs <- function(ref, L) {
  row <- lmerTest::contest1D(ref, L, rhs = 0, ddf = "Satterthwaite")
  row[["Estimate"]][[1L]] - row[["Std. Error"]][[1L]]
}

mm_expect_scalar_lmerTest_parity <- function(case, method) {
  pair <- mm_lmerTest_fit_pair(case)
  L <- mm_scalar_contrast_vector(pair$lmerTest)
  rhs <- mm_lmerTest_scalar_rhs(pair$lmerTest, L)
  observed <- contrast(pair$mixeff, L, rhs = rhs, method = method)$table
  expected <- mm_lmerTest_reference_contrast(pair$lmerTest, L, rhs, method)
  label <- sprintf("%s scalar contrast", method)
  tol <- list(
    estimate = 1e-6,
    std_error = 1e-3,
    df = 1e-2,
    statistic = 1e-3,
    p_value = 1e-3
  )

  mm_expect_numeric_close(observed$estimate, expected$estimate, tol$estimate,
                          sprintf("%s estimate", label), case)
  mm_expect_numeric_close(observed$std_error, expected$std_error, tol$std_error,
                          sprintf("%s std_error", label), case)
  if (is.na(expected$df)) {
    testthat::expect_true(is.na(observed$df),
                          info = sprintf("%s df should be NA for case `%s`",
                                         label, case$id))
  } else {
    mm_expect_numeric_close(observed$df, expected$df, tol$df,
                            sprintf("%s df", label), case)
  }
  mm_expect_numeric_close(observed$statistic, expected$statistic, tol$statistic,
                          sprintf("%s statistic", label), case)
  mm_expect_numeric_close(observed$p_value, expected$p_value, tol$p_value,
                          sprintf("%s p_value", label), case)
  testthat::expect_identical(observed$statistic_name, expected$statistic_name)
  testthat::expect_identical(observed$status, "available")
  testthat::expect_identical(observed$requested_method, method)
  if (identical(method, "asymptotic")) {
    testthat::expect_identical(observed$method, "asymptotic_wald_z")
  } else {
    testthat::expect_identical(observed$method, method)
  }
  testthat::expect_false(is.null(observed$details[[1L]]$contrast_family))
  invisible(pair)
}

mm_lmerTest_anova_table <- function(ref, method) {
  ddf <- switch(
    method,
    satterthwaite = "Satterthwaite",
    kenward_roger = {
      mm_skip_if_no_pbkrtest()
      "Kenward-Roger"
    },
    stop(sprintf("Unsupported lmerTest ANOVA method `%s`.", method), call. = FALSE)
  )
  out <- as.data.frame(stats::anova(ref, type = 3, ddf = ddf))
  out$term <- rownames(out)
  rownames(out) <- NULL
  out
}

mm_lmerTest_anova_row <- function(ref, term, method) {
  table <- mm_lmerTest_anova_table(ref, method)
  row <- table[table$term == term, , drop = FALSE]
  if (!nrow(row)) {
    testthat::skip(sprintf("lmerTest ANOVA did not return term `%s`", term))
  }
  row
}

mm_mixeff_term_row <- function(fit, term, method, source = c("test_effect", "anova")) {
  source <- match.arg(source)
  table <- if (identical(source, "test_effect")) {
    test_effect(fit, term, method = method)$table
  } else {
    stats::anova(fit, type = "III", method = method)$table
  }
  row <- table[table$term == term, , drop = FALSE]
  testthat::expect_equal(nrow(row), 1L,
                         info = sprintf("%s returned %d row(s) for term `%s`",
                                        source, nrow(row), term))
  row
}

mm_term_statistic_as_f <- function(row) {
  if (identical(row$statistic_name[[1L]], "t")) {
    row$statistic[[1L]]^2
  } else {
    row$statistic[[1L]]
  }
}

mm_expect_p_value_close <- function(observed, expected, tolerance, label, case) {
  if (is.finite(expected) && expected < 1e-12) {
    testthat::expect_true(
      is.finite(observed) && observed <= 1e-12,
      info = sprintf(
        "%s parity failed for case `%s` using %s; expected tiny p=%s, observed=%s",
        label,
        case$id,
        case$formula,
        format(expected, scientific = TRUE),
        format(observed, scientific = TRUE)
      )
    )
  } else {
    mm_expect_numeric_close(observed, expected, tolerance, label, case)
  }
}

mm_expect_term_lmerTest_parity <- function(case, term, method, source) {
  pair <- mm_lmerTest_fit_pair(case)
  observed <- mm_mixeff_term_row(pair$mixeff, term, method, source = source)
  expected <- mm_lmerTest_anova_row(pair$lmerTest, term, method)
  label <- sprintf("%s %s term `%s`", source, method, term)
  expected_num_df <- expected[["NumDF"]][[1L]]
  expected_den_df <- expected[["DenDF"]][[1L]]
  expected_f <- expected[["F value"]][[1L]]
  expected_p <- expected[["Pr(>F)"]][[1L]]
  observed_f <- mm_term_statistic_as_f(observed)
  details <- observed$details[[1L]]
  contrast_family <- details$contrast_family

  testthat::expect_identical(observed$status, "available")
  testthat::expect_identical(observed$requested_method, method)
  testthat::expect_identical(observed$method, method)
  testthat::expect_false(is.null(contrast_family))
  testthat::expect_equal(contrast_family$effective_rank, expected_num_df,
                         info = sprintf("%s contrast-family rank mismatch", label))

  if (identical(observed$statistic_name[[1L]], "f")) {
    mm_expect_numeric_close(observed$num_df, expected_num_df, 1e-12,
                            sprintf("%s numerator df", label), case)
  } else {
    testthat::expect_true(is.na(observed$num_df),
                          info = sprintf("%s scalar t row should not report top-level num_df",
                                         label))
  }
  mm_expect_numeric_close(observed$den_df, expected_den_df, 1e-2,
                          sprintf("%s denominator df", label), case)
  mm_expect_numeric_close(observed_f, expected_f, 1e-3,
                          sprintf("%s F-equivalent statistic", label), case)
  mm_expect_p_value_close(observed$p_value, expected_p, 1e-3,
                          sprintf("%s p_value", label), case)

  if (identical(method, "kenward_roger")) {
    testthat::expect_false(is.null(details$kenward_roger),
                           info = sprintf("%s missing KR detail metadata", label))
  }
  invisible(pair)
}

mm_expect_core_lme4_parity <- function(case) {
  pair <- mm_lme4_fit_pair(case)
  tol <- mm_lme4_case_tolerances(case)
  fit <- pair$mixeff
  ref <- pair$lme4

  # Names AND order must be lme4-identical (the coef-map contract): every
  # parity case doubles as a regression test for the renaming layer.
  testthat::expect_identical(
    names(fixef(fit)), names(lme4::fixef(ref)),
    info = sprintf("fixef name/order parity failed for case `%s`", case$id)
  )
  mm_assert_parity(fixef(fit), lme4::fixef(ref), case$id, "fixef", tol$fixef,
                  "fixef")
  mm_assert_parity(sigma(fit), sigma(ref), case$id, "sigma", tol$sigma,
                  "sigma")
  mm_assert_parity(as.numeric(logLik(fit)), as.numeric(stats::logLik(ref)),
                  case$id, "logLik", tol$logLik,
                  "logLik")
  mm_assert_parity(deviance(fit), suppressWarnings(deviance(ref)),
                  case$id, "deviance", tol$deviance,
                  "deviance")
  mm_assert_parity(AIC(fit), AIC(ref), case$id, "AIC", tol$AIC,
                  "AIC")
  mm_assert_parity(BIC(fit), BIC(ref), case$id, "BIC", tol$BIC,
                  "BIC")
  mm_assert_parity(fitted(fit), fitted(ref), case$id, "fitted", tol$fitted,
                  "fitted")
  mm_assert_parity(residuals(fit), residuals(ref), case$id, "residuals", tol$residuals,
                  "residuals")

  testthat::expect_identical(nobs(fit), stats::nobs(ref),
                             info = sprintf("nobs parity failed for case `%s`", case$id))
  testthat::expect_identical(df.residual(fit), df.residual(ref),
                             info = sprintf("df.residual parity failed for case `%s`", case$id))

  mm_assert_parity(stats::model.matrix(fit), lme4::getME(ref, "X"),
                  case$id, "model_matrix", tol$model_matrix,
                  "fixed model matrix", mode = "absolute")
  mm_expect_varcorr_lme4_parity(fit, ref, tol$varcorr, case)

  theta_fit <- tryCatch(getME(fit, "theta"), error = function(cnd) NULL)
  theta_ref <- lme4::getME(ref, "theta")
  if (!is.null(theta_fit) && length(theta_fit) == length(theta_ref)) {
    mm_assert_parity(theta_fit, theta_ref, case$id, "theta", tol$theta,
                    "theta")
  }

  invisible(pair)
}

mm_expect_varcorr_lme4_parity <- function(fit, ref, tolerance, case) {
  observed <- VarCorr(fit)
  expected <- as.data.frame(lme4::VarCorr(ref))
  expected_diagonal <- expected[is.na(expected$var2) & expected$grp != "Residual", ,
                                drop = FALSE]

  for (i in seq_len(nrow(expected_diagonal))) {
    hit <- mm_lme4_group_key(observed$table$group) ==
      mm_lme4_group_key(expected_diagonal$grp[[i]]) &
      observed$table$name == expected_diagonal$var1[[i]]
    testthat::expect_true(
      any(hit),
      info = sprintf(
        "VarCorr component `%s/%s` missing for case `%s`",
        expected_diagonal$grp[[i]],
        expected_diagonal$var1[[i]],
        case$id
      )
    )
    if (any(hit)) {
      mm_assert_parity(
        observed$table$std_dev[hit][1L],
        expected_diagonal$sdcor[[i]],
        case$id,
        sprintf("varcorr.std_dev.%s.%s",
                mm_lme4_group_key(expected_diagonal$grp[[i]]),
                expected_diagonal$var1[[i]]),
        tolerance,
        sprintf("VarCorr std_dev[%s/%s]",
                expected_diagonal$grp[[i]],
                expected_diagonal$var1[[i]])
      )
    }
  }

  residual_sd <- expected$sdcor[expected$grp == "Residual"][1L]
  mm_assert_parity(
    observed$residual_sd,
    residual_sd,
    case$id,
    "varcorr.residual_sd",
    tolerance,
    "VarCorr residual_sd"
  )
}

mm_expect_ranef_lme4_parity <- function(case) {
  pair <- mm_lme4_fit_pair(case)
  tol <- mm_lme4_case_tolerances(case)
  fit_re <- ranef(pair$mixeff)
  ref_re <- lme4::ranef(pair$lme4)
  fit_keys <- stats::setNames(names(fit_re), mm_lme4_group_key(names(fit_re)))
  ref_keys <- stats::setNames(names(ref_re), mm_lme4_group_key(names(ref_re)))
  common_keys <- intersect(names(fit_keys), names(ref_keys))
  if (!length(common_keys)) {
    testthat::skip(sprintf(
      "No common ranef labels for case `%s`; extractor parity is covered by scalar/fitted fields",
      case$id
    ))
  }
  for (key in common_keys) {
    fit_group <- fit_keys[[key]]
    ref_group <- ref_keys[[key]]
    fit_df <- fit_re[[fit_group]]
    ref_df <- ref_re[[ref_group]]
    fit_row_keys <- mm_lme4_level_key(rownames(fit_df))
    ref_row_keys <- mm_lme4_level_key(rownames(ref_df))
    common_rows <- intersect(fit_row_keys, ref_row_keys)
    testthat::expect_true(length(common_rows) > 0L,
                          info = sprintf("No common ranef levels for case `%s` group `%s`",
                                         case$id, ref_group))
    fit_df <- fit_df[match(common_rows, fit_row_keys), , drop = FALSE]
    ref_df <- ref_df[match(common_rows, ref_row_keys), , drop = FALSE]
    mm_assert_parity(as.matrix(fit_df), as.matrix(ref_df),
                    case$id, sprintf("ranef.%s", ref_group),
                    tol$ranef, sprintf("ranef[%s]", ref_group))
  }
  invisible(pair)
}

mm_expect_prediction_lme4_parity <- function(case) {
  pair <- mm_lme4_fit_pair(case)
  tol <- mm_lme4_case_tolerances(case)
  mm_assert_parity(stats::predict(pair$mixeff),
                   stats::predict(pair$lme4),
                   case$id, "predict_conditional", tol$fitted,
                   "predict conditional fitted-data")
  mm_assert_parity(stats::predict(pair$mixeff, re.form = NA),
                   stats::predict(pair$lme4, re.form = NA),
                   case$id, "predict_fixed_only", tol$fitted,
                   "predict fixed-only fitted-data")
  testthat::expect_equal(stats::predict(pair$mixeff), fitted(pair$mixeff),
                         ignore_attr = TRUE,
                         info = sprintf("predict()/fitted() mismatch for case `%s`",
                                        case$id))
  # lme4::predict.merMod offers no conditional se.fit, so there is no lme4
  # reference value here. Assert the engine prediction-variance contract
  # instead: finite positive SEs whose fixed component reproduces the
  # population (re.form = NA) Wald SE.
  se <- stats::predict(pair$mixeff, se.fit = TRUE)
  testthat::expect_equal(se$fit, stats::predict(pair$mixeff),
                         ignore_attr = TRUE,
                         info = sprintf("predict(se.fit=TRUE)$fit mismatch for case `%s`",
                                        case$id))
  testthat::expect_true(all(is.finite(se$se.fit)) && all(se$se.fit > 0),
                        info = sprintf("conditional se.fit should be finite and positive for case `%s`",
                                       case$id))
  pop_se <- stats::predict(pair$mixeff, re.form = NA, se.fit = TRUE)$se.fit
  pv <- mm_lmm_prediction_variance(pair$mixeff, pair$mixeff$model_frame,
                                   FALSE, 0.95)
  testthat::expect_equal(sqrt(pv$fixed_variance), unname(pop_se),
                         tolerance = 1e-6,
                         info = sprintf("engine fixed_variance should match population se for case `%s`",
                                        case$id))

  # newdata prediction is wired through mm_lmm_predict_new_json (Stage C.1,
  # bd-01KRCKCZJ5B5AQS5BV77VMM8ZF). Re-running on the training rows must
  # reproduce the in-sample fitted values within the case tolerance.
  newdata_conditional <- stats::predict(pair$mixeff, newdata = pair$data,
                                        re.form = NULL)
  mm_assert_parity(newdata_conditional, stats::predict(pair$lme4),
                   case$id, "predict_newdata_conditional", tol$fitted,
                   "predict conditional with newdata=training")
  newdata_population <- stats::predict(pair$mixeff, newdata = pair$data,
                                       re.form = NA)
  mm_assert_parity(newdata_population,
                   stats::predict(pair$lme4, re.form = NA),
                   case$id, "predict_newdata_population", tol$fitted,
                   "predict population with newdata=training")
  # `~0` is the explicit population form; behavior must mirror `re.form = NA`.
  newdata_zero <- stats::predict(pair$mixeff,
                                 re.form = stats::as.formula("~0"))
  mm_assert_parity(newdata_zero,
                   stats::predict(pair$lme4, re.form = NA),
                   case$id, "predict_re_form_zero", tol$fitted,
                   "predict re.form=~0 on in-sample data")

  ci <- stats::predict(pair$mixeff, interval = "confidence")
  testthat::expect_true(is.matrix(ci) && all(is.finite(ci)),
                        info = sprintf("conditional confidence interval should be finite for case `%s`",
                                       case$id))
  testthat::expect_true(all(ci[, "lwr"] < ci[, "fit"]) && all(ci[, "fit"] < ci[, "upr"]),
                        info = sprintf("conditional confidence interval should bracket fit for case `%s`",
                                       case$id))
  testthat::expect_equal(unname(ci[, "upr"] - ci[, "fit"]),
                         unname(stats::qnorm(0.975) * as.numeric(se$se.fit)),
                         tolerance = 1e-6,
                         info = sprintf("confidence half-width should equal z*se for case `%s`",
                                        case$id))
  invisible(pair)
}
