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
# Those remain outside this LMM parity file; mixeff::glmm() now has its own
# profiled-PIRLS bridge tests and still labels unavailable joint backends
# explicitly.

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
  dat <- utils::read.csv(brown_fixture_path("brown_rt_dummy_data.csv.gz"),
                         stringsAsFactors = FALSE)
  dat$PID <- factor(dat$PID)
  dat$stim <- factor(dat$stim)
  dat$modality <- ifelse(dat$modality == "Audio-only", 0, 1)
  dat
}

brown_rt_interaction_data <- function() {
  dat <- utils::read.csv(
    brown_fixture_path("brown_rt_dummy_data_interaction.csv.gz"),
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

brown_mkl_endpoint <- function() {
  jsonlite::read_json(
    brown_fixture_path("brown-mkl-reference-endpoint.json"),
    simplifyVector = FALSE
  )
}

brown_cases <- function() {
  list(
    rt_modality_full = list(
      data = brown_rt_data,
      formula = RT ~ 1 + modality + (1 + modality | PID) +
        (1 + modality | stim),
      REML = TRUE,
      slow = TRUE,
      # Default BOBYQA can stop at an inferior nearly singular endpoint on
      # MKL. Tight controls escape the captured endpoint (regression below).
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(rhobeg = 0.01, rhoend = 1e-9, maxfun = 20000)
      ),
      # NOTE (pin 4a2abb3, 2026-07-10): the reduced-rank/optimizer bugs that
      # used to fail rt_interaction are fixed; this maximal crossed model now
      # converges (both engines converged_interior, REML deviance 0.05 apart
      # ~= 3.6e-7 relative). Only per-observation residual parity misses at
      # fitted=2e-3 (max |resid| diff ~2.3 on a ~254 scale). Tight-tolerance
      # refit (ftol_rel=1e-12) narrows the REML gap to 0.014 but does not
      # close it: the two engines reach a REML-equivalent point whose
      # random-effect BLUPs differ slightly on the flat ridge of a maximal
      # crossed fit (PRD Decision D: statistical equivalence, not bit-exact).
      # RESOLVED (user decision, 2026-08-01): residual parity is judged on
      # the response scale (see the loop below), same tolerance as fitted;
      # the sdamr_speeddate precedent (commit 96076e4) carried the same
      # logic. Fitted tolerance stays at 2e-3.
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

# Former engine-defect gates were removed after unchanged parity checks
# passed at engine pin 977ef11 (2026-10-06).

brown_diag_dir <- function() {
  out <- Sys.getenv("MIXEFF_DIAG_OUT", unset = "")
  if (!nzchar(out)) out <- testthat::test_path("_problems")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

brown_save_diagnostics <- function(label, payload) {
  path <- tempfile(
    sprintf("brown-parity-%s-", gsub("[^[:alnum:]_-]", "-", label)),
    tmpdir = brown_diag_dir(), fileext = ".rds"
  )
  payload$session <- utils::sessionInfo()
  payload$written_at <- as.character(Sys.time())
  saveRDS(payload, path)
  path
}

brown_reference_fit <- function(ref_args, label, mixeff_fit = NULL) {
  ref_control <- ref_args$control %||% lme4::lmerControl()
  ref_control$calc.derivs <- TRUE
  # lme4 otherwise skips derivative checking above these size limits even
  # when calc.derivs is TRUE. Keep older lme4 releases compatible by only
  # overriding limits that their control object exposes.
  if (!is.null(ref_control$checkConv$check.conv.nobsmax)) {
    ref_control$checkConv$check.conv.nobsmax <- Inf
  }
  if (!is.null(ref_control$checkConv$check.conv.nparmax)) {
    ref_control$checkConv$check.conv.nparmax <- Inf
  }
  ref_args$control <- ref_control
  warnings <- character()
  messages <- character()
  ref <- tryCatch(
    withCallingHandlers(
      do.call(lme4::lmer, ref_args),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) e
  )
  if (inherits(ref, "error")) {
    receipt <- brown_save_diagnostics(label, list(
      stage = "reference_fit_error", reference_error = conditionMessage(ref),
      reference_warnings = warnings, reference_messages = messages,
      reference_args = ref_args, mixeff = mixeff_fit
    ))
    stop(sprintf("lme4 reference failed for `%s`; diagnostics saved to %s: %s",
                 label, receipt, conditionMessage(ref)), call. = FALSE)
  }
  optinfo <- ref@optinfo
  lme4_messages <- optinfo$conv$lme4$messages %||% character()
  all_messages <- unique(c(messages, lme4_messages))
  opt_code <- optinfo$conv$opt
  derivs <- optinfo$derivs
  expected_boundary <- grepl("^boundary \\(singular\\) fit", all_messages)
  unexpected_messages <- all_messages[!expected_boundary]
  known_warning <- grepl("^boundary \\(singular\\) fit", warnings)
  unexpected_warnings <- warnings[!known_warning]
  derivatives_ok <- !is.null(derivs$gradient) && !is.null(derivs$Hessian) &&
    all(is.finite(derivs$gradient)) && all(is.finite(derivs$Hessian))
  # optinfo$val is the optimizer parameter vector, not its objective. The
  # fitted likelihood supplies the finite REML/ML deviance criterion.
  criterion <- -2 * as.numeric(stats::logLik(ref))
  criterion_ok <- length(criterion) == 1L && is.finite(criterion)
  # A boundary/singular fit is a valid lme4 optimum. Optimizer convergence
  # messages and nonzero optimizer codes are not an oracle for parity.
  optimizer_failure <- is.null(opt_code) || length(opt_code) != 1L ||
    !is.finite(opt_code) || opt_code != 0L || !derivatives_ok || !criterion_ok ||
    length(unexpected_messages) > 0L || length(unexpected_warnings) > 0L
  diagnostics <- list(
    reference_warnings = warnings, emitted_messages = messages,
    reference_messages = lme4_messages, all_messages = all_messages,
    optinfo = optinfo, reference_control = ref_control,
    singular = lme4::isSingular(ref),
    expected_boundary_messages = all_messages[expected_boundary],
    unexpected_messages = unexpected_messages,
    unexpected_warnings = unexpected_warnings,
    derivatives_ok = derivatives_ok, criterion = criterion,
    criterion_ok = criterion_ok,
    optimizer_failure = optimizer_failure
  )
  if (optimizer_failure) {
    receipt <- brown_save_diagnostics(label, c(list(
      stage = "reference_optimizer_failure", reference = ref,
      reference_args = ref_args, mixeff = mixeff_fit
    ), diagnostics))
    stop(sprintf("lme4 reference for `%s` did not converge; diagnostics saved to %s",
                 label, receipt), call. = FALSE)
  }
  attr(ref, "brown_reference_diagnostics") <- diagnostics
  ref
}

brown_with_diagnostics <- function(label, pair, expr) {
  tryCatch(
    force(expr),
    error = function(e) {
      receipt <- brown_save_diagnostics(label, list(
        stage = "parity_failure", error = conditionMessage(e),
        mixeff = pair$mixeff, lme4 = pair$lme4,
        reference_diagnostics = attr(pair$lme4, "brown_reference_diagnostics")
      ))
      stop(sprintf("%s (Brown diagnostics: %s)", conditionMessage(e), receipt),
           call. = FALSE)
    }
  )
}

brown_fit_pair <- function(case, label = "unlabelled") {
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
  # Derivatives are needed for lme4's convergence assessment. Preserve a
  # caller's other controls while making that verification explicit.
  mixeff_fit <- mixeff::lmm(case$formula, dat, REML = case$REML,
                            control = mixeff::mm_control(verbose = -1))
  list(mixeff = mixeff_fit,
       lme4 = brown_reference_fit(ref_args, label, mixeff_fit))
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

brown_expect_pair <- function(pair, tolerance, label) {
  brown_with_diagnostics(label, pair, {
    brown_expect_fixef_close(pair$mixeff, pair$lme4, tolerance$fixef, label)
    expect_equal(sigma(pair$mixeff), sigma(pair$lme4), tolerance = tolerance$scalar,
                 info = sprintf("sigma parity failed for `%s`", label))
    expect_equal(as.numeric(logLik(pair$mixeff)),
                 as.numeric(stats::logLik(pair$lme4)), tolerance = tolerance$scalar,
                 info = sprintf("logLik parity failed for `%s`", label))
    expect_equal(AIC(pair$mixeff), AIC(pair$lme4), tolerance = tolerance$scalar,
                 info = sprintf("AIC parity failed for `%s`", label))
    expect_equal(BIC(pair$mixeff), BIC(pair$lme4), tolerance = tolerance$scalar,
                 info = sprintf("BIC parity failed for `%s`", label))
    expect_equal(fitted(pair$mixeff), fitted(pair$lme4), tolerance = tolerance$fitted,
                 info = sprintf("fitted-value parity failed for `%s`", label))
    # Residual and fitted gaps agree on the shared response scale.
    resid_gap <- max(abs(residuals(pair$mixeff) - residuals(pair$lme4)))
    resid_scale <- max(abs(fitted(pair$lme4)))
    expect_lt(resid_gap / resid_scale, tolerance$fitted,
              label = sprintf("response-scaled residual gap for `%s`", label))
    brown_expect_varcorr_close(pair$mixeff, pair$lme4, tolerance$varcorr, label)
  })
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
    pair <- brown_fit_pair(case, label)
    brown_expect_pair(pair, case$tolerance, label)
  }
})

test_that("Brown 2021 large RT tutorial examples match lme4", {
  mm_skip_if_no_lme4()

  for (label in names(brown_cases())) {
    case <- brown_cases()[[label]]
    if (!isTRUE(case$slow)) {
      next
    }
    pair <- brown_fit_pair(case, label)
    brown_expect_pair(pair, case$tolerance, label)
  }
})

test_that("Brown 2021 RT likelihood-ratio example agrees with lme4", {
  mm_skip_if_no_lme4()
  dat <- brown_rt_data()
  full <- RT ~ 1 + modality + (1 + modality | PID) + (1 + modality | stim)
  reduced <- RT ~ 1 + (1 + modality | stim) + (1 + modality | PID)

  mm_full <- mixeff::lmm(full, dat, REML = FALSE,
                         control = mixeff::mm_control(verbose = -1))
  mm_reduced <- mixeff::lmm(reduced, dat, REML = FALSE,
                            control = mixeff::mm_control(verbose = -1))
  lme4_full <- brown_reference_fit(
    list(formula = full, data = dat, REML = FALSE,
         control = lme4::lmerControl(optimizer = "bobyqa", calc.derivs = TRUE)),
    "lrt_full", mm_full
  )
  lme4_reduced <- brown_reference_fit(
    list(formula = reduced, data = dat, REML = FALSE,
         control = lme4::lmerControl(optimizer = "bobyqa", calc.derivs = TRUE)),
    "lrt_reduced", mm_reduced
  )

  mm_cmp <- mixeff::compare(mm_reduced, mm_full)$table
  lme4_cmp <- stats::anova(lme4_reduced, lme4_full)
  brown_with_diagnostics("lrt", list(mixeff = list(full = mm_full, reduced = mm_reduced),
                                      lme4 = list(full = lme4_full, reduced = lme4_reduced)), {
    expect_equal(mm_cmp$LRT[[2L]], lme4_cmp$Chisq[[2L]], tolerance = 5e-3)
    expect_equal(mm_cmp$p_value[[2L]], lme4_cmp$`Pr(>Chisq)`[[2L]], tolerance = 1e-4)
  })
})

test_that("Brown captured lme4 endpoint is improved by a tight validated refit", {
  mm_skip_if_no_lme4()
  endpoint <- brown_mkl_endpoint()
  dat <- brown_rt_data()
  expect_identical(
    unname(tools::md5sum(brown_fixture_path("brown_rt_dummy_data.csv.gz"))),
    endpoint$fixture_md5
  )
  theta <- unlist(endpoint$theta, use.names = TRUE)
  captured_control <- lme4::lmerControl(
    optimizer = endpoint$captured_control$optimizer,
    calc.derivs = TRUE
  )
  tight_control <- lme4::lmerControl(
    optimizer = endpoint$tight_refit$control$optimizer,
    optCtrl = list(
      rhobeg = endpoint$tight_refit$control$rhobeg,
      rhoend = endpoint$tight_refit$control$rhoend,
      maxfun = endpoint$tight_refit$control$maxfun
    ),
    calc.derivs = TRUE
  )
  formula <- brown_cases()$rt_modality_full$formula
  expected_cnms <- lapply(endpoint$cnms, unlist, use.names = FALSE)
  # lme4's term order defines the theta slots. Check it before evaluating the
  # captured numeric vector or using it as an optimizer start.
  expect_identical(lme4::lFormula(formula, data = dat)$reTrms$cnms, expected_cnms)
  # This evaluates the recorded endpoint directly; it does not rely on a
  # current default optimizer run reproducing the historical failure.
  devfun <- lme4::lmer(formula, data = dat, REML = TRUE,
                        control = captured_control, devFunOnly = TRUE)
  captured_objective <- as.numeric(devfun(theta))
  expect_lt(abs(captured_objective - endpoint$common_objective), 1e-6)
  mixeff_fit <- mixeff::lmm(formula, dat, REML = TRUE,
                            control = mixeff::mm_control(verbose = -1))
  refit <- brown_reference_fit(
    list(formula = formula, data = dat, REML = TRUE, control = tight_control,
         start = list(theta = theta)),
    "mkl-captured-endpoint", mixeff_fit
  )
  expect_identical(names(lme4::getME(refit, "theta")), names(theta))
  refit_objective <- -2 * as.numeric(stats::logLik(refit))
  expect_gt(captured_objective - refit_objective, 0.1)
  expect_lt(abs(refit_objective - endpoint$tight_refit$objective), 1e-5)
  brown_expect_pair(
    list(mixeff = mixeff_fit, lme4 = refit),
    brown_cases()$rt_modality_full$tolerance, "rt_modality_full_mkl_refit"
  )
})

test_that("Brown reference harness records a real optimizer nonconvergence", {
  mm_skip_if_no_lme4()
  old <- Sys.getenv("MIXEFF_DIAG_OUT", unset = NA_character_)
  out <- tempfile("brown-diag-")
  Sys.setenv(MIXEFF_DIAG_OUT = out)
  on.exit(if (is.na(old)) Sys.unsetenv("MIXEFF_DIAG_OUT") else Sys.setenv(MIXEFF_DIAG_OUT = old), add = TRUE)
  dat <- brown_figure_data()
  mixeff_fit <- mixeff::lmm(yvar ~ 1 + xvar + (1 | PID), dat, REML = TRUE,
                            control = mixeff::mm_control(verbose = -1))
  expect_error(
    brown_reference_fit(
      list(formula = yvar ~ 1 + xvar + (1 | PID), data = dat, REML = TRUE,
           control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1),
                                       calc.derivs = TRUE)),
      "low-maxfun", mixeff_fit
    ),
    "did not converge"
  )
  receipt <- list.files(out, pattern = "low-maxfun.*rds$", full.names = TRUE)
  expect_length(receipt, 1L)
  saved <- readRDS(receipt[[1L]])
  expect_identical(saved$stage, "reference_optimizer_failure")
  expect_s3_class(saved$mixeff, "mm_lmm")
  expect_s4_class(saved$reference, "lmerMod")
})

test_that("Brown reference harness accepts a converged singular boundary fit", {
  mm_skip_if_no_lme4()
  d <- data.frame(g = gl(8, 5), x = rep(seq(-1, 1, length.out = 5), 8))
  # The repeated within-group residual pattern gives every group exactly the
  # same mean, producing a genuine zero random-intercept boundary solution.
  d$y <- 1 + 0.7 * d$x + rep(c(-0.1, 0.1, -0.05, 0.05, 0), 8)
  ref <- brown_reference_fit(
    list(formula = y ~ x + (1 | g), data = d, REML = TRUE,
         control = lme4::lmerControl(optimizer = "bobyqa", calc.derivs = TRUE)),
    "singular-boundary"
  )
  diag <- attr(ref, "brown_reference_diagnostics")
  expect_true(lme4::isSingular(ref))
  expect_false(diag$optimizer_failure)
  expect_true(any(grepl("boundary \\(singular\\) fit", diag$expected_boundary_messages)))
})

test_that("Brown parity failures save both fits and diagnostics", {
  old <- Sys.getenv("MIXEFF_DIAG_OUT", unset = NA_character_)
  out <- tempfile("brown-diag-")
  Sys.setenv(MIXEFF_DIAG_OUT = out)
  on.exit(if (is.na(old)) Sys.unsetenv("MIXEFF_DIAG_OUT") else Sys.setenv(MIXEFF_DIAG_OUT = old), add = TRUE)
  case <- brown_cases()$figure_random_intercepts
  pair <- brown_fit_pair(case, "deliberate-source")
  expect_error(brown_with_diagnostics(
    "deliberate", pair,
    expect_equal(mixeff::fixef(pair$mixeff), rep(NA_real_, length(mixeff::fixef(pair$mixeff))))
  ),
               "Brown diagnostics")
  receipt <- list.files(out, pattern = "deliberate.*rds$", full.names = TRUE)
  expect_length(receipt, 1L)
  saved <- readRDS(receipt[[1L]])
  expect_s3_class(saved$mixeff, "mm_lmm")
  expect_s4_class(saved$lme4, "lmerMod")
  expect_equal(mixeff::fixef(saved$mixeff), mixeff::fixef(pair$mixeff))
  expect_equal(lme4::fixef(saved$lme4), lme4::fixef(pair$lme4))
})
