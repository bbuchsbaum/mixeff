# Convergence-honesty invariants (Audit1 remediation WI-1.4): a converged_*
# status must never coexist with non-finite claimed quantities, on any model
# class or estimator. See helper-convergence-invariants.R for the checker.

mm_ci_lmm_data <- function() {
  set.seed(42)
  g <- factor(rep(1:12, each = 8))
  x <- rnorm(96)
  y <- 2 + 0.5 * x + rep(rnorm(12, sd = 0.7), each = 8) + rnorm(96, sd = 0.4)
  data.frame(y = y, x = x, g = g)
}

mm_ci_glmm_data <- function() {
  set.seed(43)
  g <- factor(rep(1:15, each = 12))
  x <- rnorm(180)
  eta <- -0.3 + 0.8 * x + rep(rnorm(15, sd = 0.6), each = 12)
  data.frame(y = rbinom(180, 1, plogis(eta)), x = x, g = g)
}

test_that("converged LMM fits claim only finite quantities", {
  fit <- lmm(y ~ x + (1 | g), mm_ci_lmm_data(),
             control = mm_control(verbose = -1))
  expect_true(startsWith(as.character(fit$fit_status), "converged"))
  mm_expect_convergence_honesty(fit, "interior LMM")
})

test_that("boundary LMM fits claim only finite quantities", {
  # Near-zero group variance drives theta to the boundary; the certificate
  # may say converged_boundary, and every claimed number must still be finite.
  set.seed(44)
  d <- mm_ci_lmm_data()
  d$y <- 2 + 0.5 * d$x + rnorm(nrow(d), sd = 0.4) # no true group effect
  fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
  mm_expect_convergence_honesty(fit, "boundary LMM")
})

test_that("profiled and joint GLMM fits claim only finite quantities", {
  d <- mm_ci_glmm_data()
  fp <- glmm(y ~ x + (1 | g), d, family = binomial(),
             control = mm_control(verbose = -1))
  mm_expect_convergence_honesty(fp, "profiled GLMM")

  fj <- glmm(y ~ x + (1 | g), d, family = binomial(),
             method = "joint_laplace", control = mm_control(verbose = -1))
  mm_expect_convergence_honesty(fj, "joint GLMM")
})

test_that("the honesty checker fires on a corrupted converged payload", {
  # Prove the guard can fail: a fake fit claiming converged_interior with a
  # non-finite objective/logLik must be flagged, not passed through.
  fake <- list(
    fit_status = "converged_interior",
    logLik = -Inf,
    beta = c(1, NaN),
    theta = 0.5,
    fitted = c(0.2, 0.8),
    artifact = list(optimizer_certificate = list(objective_value = Inf))
  )
  problems <- mm_check_convergence_honesty(fake)
  expect_setequal(problems, c("logLik", "beta", "certificate_objective"))

  # And a non-converged status is out of scope for this invariant (the
  # refusal/labelling contracts cover it), regardless of field contents.
  fake$fit_status <- "not_optimized"
  expect_identical(mm_check_convergence_honesty(fake), character(0))

  # A corrupted NA status must fail cleanly (out of scope), not crash.
  fake$fit_status <- NA_character_
  expect_identical(mm_check_convergence_honesty(fake), character(0))
})

test_that("estimator-substitution surfaces name both estimators and the evidence", {
  # The suite's only real substituted fit lives in the minute-class OSF file;
  # cover the notice/note text cheaply against a faked typed record with the
  # exact fields the engine emits (EstimatorSubstitution, engine f82c646).
  sub <- list(
    requested_method = "joint_laplace",
    effective_method = "fast_pirls_profiled",
    requested_fit_status = "not_optimized",
    requested_return_code = "JOINT_LAPLACE:XTOL_REACHED",
    requested_free_gradient_norm = NA_real_,
    reason = "joint GLMM did not certify; returning labelled fast-PIRLS fallback"
  )
  msg <- mixeff:::mm_glmm_substitution_notice(sub)
  expect_match(msg, "joint_laplace", fixed = TRUE)
  expect_match(msg, "fast_pirls_profiled", fixed = TRUE)
  expect_match(msg, "verify_convergence", fixed = TRUE)

  note <- mixeff:::mm_estimator_substitution_note(sub)
  expect_match(note, "estimator substitution", fixed = TRUE)
  expect_match(note, "fallback", fixed = TRUE)
})

test_that("mm_estimator_substitution is NULL on a fit with no substitution", {
  fit <- glmm(y ~ x + (1 | g), mm_ci_glmm_data(), family = binomial(),
              control = mm_control(verbose = -1))
  expect_null(mixeff:::mm_estimator_substitution(fit))
})
