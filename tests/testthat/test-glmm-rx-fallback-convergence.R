# Engine follow-ups (mixeff-rs afc7c36): a joint-Laplace GLMM whose active
# Hessian is not usable reports Wald rows from RX conditional on theta
# (covariance_method laplace_rx_conditional_on_theta, reliability low,
# reliability_reason glmm_laplace_rx_conditional_on_theta_wald). The only
# known trigger is the aphantasia primary model (~2 min; covered by the
# MIXEFF_RUN_APHANTASIA block in test-aphantasia-reproduction.R), so these
# tests rewrite a fast cbpp joint fit's payloads into the RX shape the
# engine emits and check the R presentation. The same fixture covers the
# "convergence not certified" presentation of a not_optimized joint fit
# with a small Newton-decrement objective gap.

rxf_base_fit <- local({
  cache <- NULL
  function() {
    skip_if_not_installed("lme4")
    if (is.null(cache)) {
      env <- new.env()
      utils::data("cbpp", package = "lme4", envir = env)
      cache <<- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
                     env$cbpp, family = stats::binomial(),
                     control = mm_control(verbose = -1))
    }
    cache
  }
})

rxf_notes <- c(
  paste0("joint-laplace active Hessian unusable for Wald inference: ",
         "joint-laplace GLMM Hessian is not positive definite on the active ",
         "parameter space: min eigenvalue -8.8e-1 <= tolerance 4.9e-4"),
  paste0("fixed-effect covariance falls back to the fixed-effect block (RX) ",
         "of the Laplace PLS factorization at the joint optimum, conditional ",
         "on theta (as lme4 vcov(fit, use.hessian = FALSE))")
)

# Rewrite the fit into the engine's RX-fallback payload shape.
rxf_as_rx_fallback <- function(fit) {
  rows <- fit$artifact$fixed_effect_inference_table$rows
  fit$artifact$fixed_effect_inference_table$rows <- lapply(rows, function(r) {
    r$covariance_method <- "laplace_rx_conditional_on_theta"
    r$reliability <- "low"
    r$reliability_reason <- "glmm_laplace_rx_conditional_on_theta_wald"
    r$notes <- as.list(rxf_notes)
    r
  })
  cov <- fit$artifact$fixed_effect_covariance_matrix
  cov$method <- "laplace_rx_conditional_on_theta"
  cov$reliability <- "low"
  cov$notes <- as.list(rxf_notes)
  fit$artifact$fixed_effect_covariance_matrix <- cov
  if (!is.null(fit$fixed_effect_vcov)) {
    attr(fit$fixed_effect_vcov, "mm_method") <- "laplace_rx_conditional_on_theta"
    attr(fit$fixed_effect_vcov, "mm_reliability") <- "low"
    attr(fit$fixed_effect_vcov, "mm_notes") <- rxf_notes
  }
  fit
}

# Rewrite the certificate into a not_optimized joint stop with the given
# eager objective gap (the aphantasia primary shape).
rxf_as_not_optimized <- function(fit, gap, budget_exhausted = FALSE) {
  cert <- fit$artifact$optimizer_certificate
  cert$status <- "not_optimized"
  cert$stationarity_decrement$eager <- list(
    variant = "beta_block_theta_diagonal",
    verdict = "exceeds_tolerance",
    objective_gap = gap
  )
  cert$stationarity_decrement$gap_tolerance <- 1e-6
  cert$evidence$optimizer_stop$acceptable_stop <- TRUE
  cert$evidence$optimizer_stop$budget_exhausted <- budget_exhausted
  fit$artifact$optimizer_certificate <- cert
  fit$fit_status <- "not_optimized"
  fit
}

test_that("RX-fallback Wald rows are labelled in summary(), vcov() and inference_options()", {
  base <- rxf_base_fit()
  expect_identical(attr(base$fixed_effect_vcov, "mm_method"),
                   "joint_laplace_active_hessian")
  expect_no_warning(vcov(base))
  fit <- rxf_as_rx_fallback(base)

  # summary(): the rows stay available (glmer reports RX SEs too), graded
  # low, with the Hessian failure in the notes.
  s <- summary(fit)
  expect_true(all(s$inference$table$status == "available"))
  expect_true(all(s$inference$table$reliability_reason ==
                    "glmm_laplace_rx_conditional_on_theta_wald"))
  expect_equal(unname(s$coefficients[, "Std. Error"]),
               unname(summary(base)$coefficients[, "Std. Error"]))
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "Wald-z reliability: low", fixed = TRUE)
  expect_match(out, "fixed-effect block RX of the Laplace factorization",
               fixed = TRUE)
  expect_match(out, "use.hessian = TRUE", fixed = TRUE)
  expect_match(out, "not positive definite on the active parameter space",
               fixed = TRUE)
  expect_false(grepl("RX of the Laplace", paste(capture.output(print(
    summary(base))), collapse = "\n"), fixed = TRUE))

  # vcov(): glmer-style fallback warning, typed, with the engine notes kept.
  expect_warning(V <- vcov(fit), class = "mm_vcov_rx_fallback")
  expect_identical(attr(V, "mm_method"), "laplace_rx_conditional_on_theta")
  expect_identical(attr(V, "mm_reliability"), "low")
  expect_identical(attr(V, "mm_notes"), rxf_notes)
  w <- tryCatch(vcov(fit), warning = function(w) conditionMessage(w))
  expect_match(w, "use.hessian = FALSE", fixed = TRUE)
  expect_match(w, "not positive definite on the active", fixed = TRUE)

  # inference_options(): the current Wald route names the RX fallback, not a
  # plain certified grade, and carries the engine notes.
  io <- inference_options(fit)
  wald <- io$table[io$table$method == "asymptotic_wald_z", ]
  expect_identical(wald$expected_status, "available")
  expect_identical(wald$expected_reliability_reason,
                   "glmm_laplace_rx_conditional_on_theta_wald")
  expect_match(wald$display_reason, "RX", fixed = TRUE)
  expect_match(wald$notes, "not positive definite", fixed = TRUE)
  # The joint route's bootstrap refusal is reported, not "runs now".
  boot <- io$table[io$table$method == "glmm_parametric_bootstrap", ]
  expect_identical(boot$expected_status, "not_assessed")
  expect_identical(boot$expected_reliability_reason,
                   "glmm_bootstrap_joint_laplace_unavailable")
  expect_error(confint(fit, method = "bootstrap", nsim = 10),
               class = "mm_inference_unavailable")
})

test_that("a not_optimized joint fit with a small objective gap reads as 'convergence not certified'", {
  fit <- rxf_as_not_optimized(rxf_as_rx_fallback(rxf_base_fit()), gap = 4.4e-5)
  conv <- mixeff:::mm_glmm_convergence_assessment(fit)
  expect_true(conv$not_certified_small_gap)
  expect_equal(conv$objective_gap, 4.4e-5)
  expect_equal(conv$gap_tolerance, 1e-6)
  expect_identical(fit$fit_status, "not_optimized")  # raw label untouched

  out <- paste(capture.output(print(suppressWarnings(summary(fit)))),
               collapse = "\n")
  expect_match(out, "Fit status: not_optimized (convergence not certified; warning)",
               fixed = TRUE)
  expect_match(out, "Warning: convergence not certified", fixed = TRUE)
  expect_match(out, "objective gap 4.4e-05 deviance units", fixed = TRUE)
  expect_match(out, "verify_convergence(fit)", fixed = TRUE)
  expect_false(grepl("completed its optimization budget", out, fixed = TRUE))
  expect_match(paste(capture.output(print(fit)), collapse = "\n"),
               "convergence not certified", fixed = TRUE)

  io <- inference_options(fit)
  expect_match(io$fit_status, "convergence not certified", fixed = TRUE)
  wald <- io$table[io$table$method == "asymptotic_wald_z", ]
  expect_identical(wald$expected_status, "available")
  expect_false(identical(wald$expected_reliability_reason,
                         "glmm_wald_uncertified_for_profiled_estimator"))
  expect_match(wald$notes, "a warning, not a refusal", fixed = TRUE)
})

test_that("large gaps, exhausted budgets and converged fits keep their labels", {
  base <- rxf_base_fit()
  expect_null(mixeff:::mm_glmm_convergence_assessment(base))
  big <- rxf_as_not_optimized(base, gap = 0.5)
  expect_false(mixeff:::mm_glmm_convergence_assessment(big)$not_certified_small_gap)
  out <- paste(capture.output(print(summary(big))), collapse = "\n")
  expect_match(out, "Fit status: not_optimized\n", fixed = TRUE)
  expect_false(grepl("convergence not certified", out, fixed = TRUE))
  budget <- rxf_as_not_optimized(base, gap = 1e-5, budget_exhausted = TRUE)
  expect_false(
    mixeff:::mm_glmm_convergence_assessment(budget)$not_certified_small_gap
  )
  # Exactly at the presentation limit still counts as small.
  edge <- rxf_as_not_optimized(base, gap = mixeff:::mm_joint_gap_warning_limit)
  expect_true(mixeff:::mm_glmm_convergence_assessment(edge)$not_certified_small_gap)
})

test_that("a joint fit whose Wald rows the engine withheld is not labelled as the profiled estimator", {
  fit <- rxf_as_not_optimized(rxf_base_fit(), gap = 0.5)
  rows <- fit$artifact$fixed_effect_inference_table$rows
  fit$artifact$fixed_effect_inference_table$rows <- lapply(rows, function(r) {
    r[c("std_error", "statistic", "p_value", "statistic_name")] <- list(NULL)
    r$status <- "not_assessed"
    r$reliability <- "not_available"
    r$reliability_reason <- "glmm_joint_laplace_not_optimized"
    r$reason <- "joint-Laplace fit is not optimized; Wald inference withheld"
    r
  })
  io <- inference_options(fit)
  wald <- io$table[io$table$method == "asymptotic_wald_z", ]
  expect_identical(wald$expected_status, "not_assessed")
  expect_identical(wald$expected_reliability_reason,
                   "glmm_joint_laplace_not_optimized")
  expect_identical(wald$r_verb, "verify_convergence(fit)")
  expect_match(wald$notes, "withheld Wald inference for this joint-Laplace fit",
               fixed = TRUE)
  expect_false(grepl("profiled", wald$notes, fixed = TRUE))
})
