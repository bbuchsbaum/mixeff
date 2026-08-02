# Closed-world GLMM inference contract (Audit1 WI-2.6).
#
# THE release property, from the audit:
#   "For a pirls_profiled fit, no public path can return a finite standard
#    error, inferential covariance, test statistic, p-value, or confidence
#    interval."
# ...unless the user explicitly opted in to the labelled working-Hessian
# approximation at fit time, in which case every emitted row must carry the
# `wald_z_working_hessian` label (both directions are asserted here).
#
# vcov() note: vcov.mm_glmm returns the covariance GEOMETRY with its wire
# status attached (`available_noninferential` for profiled fits). That is
# deliberate -- prediction and diagnostics need the matrix -- and it is not
# an "inferential covariance" in the audit's sense because its status
# attribute machine-readably says so and no public route computes Wald
# quantities from it without consulting the capability arbiter. The arbiter
# <-> wire consistency test below pins that agreement.

mm_contract_data <- function(seed = 99L) {
  set.seed(seed)
  g <- factor(rep(1:14, each = 10))
  x <- rnorm(140)
  trt <- factor(rep(c("a", "b"), length.out = 140))
  eta <- -0.2 + 0.6 * x + 0.5 * (trt == "b") + rep(rnorm(14, sd = 0.5), each = 10)
  data.frame(y = rbinom(140, 1, plogis(eta)), x = x, trt = trt, g = g)
}

mm_contract_fit <- local({
  cache <- list()
  function(method = "pirls_profiled", inference = "auto") {
    key <- paste(method, inference, sep = "/")
    if (is.null(cache[[key]])) {
      cache[[key]] <<- glmm(y ~ x + trt + (1 | g), mm_contract_data(),
                            family = binomial(), method = method,
                            inference = inference,
                            control = mm_control(verbose = -1))
    }
    cache[[key]]
  }
})

test_that("no public route yields finite Wald inference from a default profiled GLMM", {
  fit <- mm_contract_fit()

  # summary(): estimates present, inference withheld.
  ct <- summary(fit, tests = "coefficients")$coefficients
  expect_true(all(is.finite(ct[, "Estimate"])))
  expect_false(any(is.finite(ct[, "Std. Error"])))
  expect_false(any(is.finite(ct[, "statistic"])))
  expect_false(any(is.finite(ct[, "p.value"])))

  # confint(): typed refusal.
  expect_error(confint(fit), class = "mm_inference_unavailable")

  # contrast(): estimate reported, SE/z/p withheld.
  w <- stats::setNames(numeric(length(fixef(fit))), names(fixef(fit)))
  w["trtb"] <- 1
  ctr <- contrast(fit, rbind(w))
  expect_true(all(is.finite(ctr$table$estimate)))
  expect_false(any(is.finite(ctr$table$std_error)))
  expect_false(any(is.finite(ctr$table$p_value)))
  expect_true(all(ctr$table$method == "not_computed"))

  # mm_lincomb(): typed refusal.
  expect_error(mm_lincomb(fit, c(x = 1)), class = "mm_inference_unavailable")

  # inference_table(): rows carry no finite SE/statistic/p.
  it <- inference_table(fit)$table
  expect_false(any(is.finite(it$std_error)))
  expect_false(any(is.finite(it$p_value)))

  # reporting tables inherit the gate through the inference table.
  # Reporting tables inherit the gate: the compact fixed-effects section on
  # a withheld fit carries estimates with method not_computed and NO
  # standard-error / p-value columns at all (columns without content are
  # dropped, which is stronger than printing NA).
  fe <- reporting_table(fit)$sections$fixed_effects
  expect_true(is.data.frame(fe) && nrow(fe) > 0L)
  expect_true(all(fe$method == "not_computed"))
  se_col <- intersect(c("std_error", "Std. Error", "std.error"), names(fe))
  if (length(se_col)) {
    expect_false(any(is.finite(fe[[se_col[1L]]])))
  }
  p_col <- intersect(c("p_value", "p.value"), names(fe))
  if (length(p_col)) {
    expect_false(any(is.finite(fe[[p_col[1L]]])))
  }

  # broom: NA inference.
  skip_if_not_installed("generics")
  td <- generics::tidy(fit, effects = "fixed")
  expect_false(any(is.finite(td$std.error)))
  expect_false(any(is.finite(td$p.value)))
})

test_that("GLMM verbs without a validated route refuse with typed errors", {
  fit <- mm_contract_fit()
  for (call in list(
    quote(test_effect(fit, "trt")),
    quote(mm_means(fit, ~trt)),
    quote(mm_comparisons(fit, ~trt)),
    quote(mm_predictions(fit)),
    quote(mm_grid(fit, ~trt))
  )) {
    err <- tryCatch(eval(call), error = function(e) e)
    expect_s3_class(err, "mm_inference_unavailable")
    expect_false(inherits(err, "simpleError"))
  }
})

test_that("emmeans surface refuses a default profiled GLMM", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("estimability")
  fit <- mm_contract_fit()
  # emmeans catches recover_data errors with try() -- our refusal text is
  # printed to the console and emmeans then throws its own generic error --
  # so the pipeline is asserted to fail with our message visible, and the
  # typed class is asserted on the direct entry point.
  printed <- capture.output(
    expect_error(emmeans::emmeans(fit, ~trt)),
    type = "message"
  )
  expect_true(any(grepl("withheld for this GLMM estimator", printed,
                        fixed = TRUE)))
  expect_error(recover_data.mm_glmm(fit), class = "mm_inference_unavailable")
})

test_that("the working-Hessian opt-in labels every route it unlocks", {
  fit <- mm_contract_fit(inference = "working_hessian")

  ct <- summary(fit, tests = "coefficients")$coefficients
  expect_true(all(is.finite(ct[, "Std. Error"])))
  s <- summary(fit, tests = "coefficients")
  expect_true(all(s$inference$table$method == "wald_z_working_hessian"))
  expect_true(all(s$inference$table$status == "available_noninferential"))
  out <- capture.output(print(s))
  expect_true(any(grepl("UNCERTIFIED working-Hessian", out, fixed = TRUE)))

  ci <- confint(fit)
  expect_true(all(is.finite(as.matrix(ci))))
  # The interval object itself must carry the approximation's certificate,
  # not the certified-table one (review finding: false certificate).
  expect_identical(attr(ci, "method"), "wald_z_working_hessian")
  expect_identical(attr(ci, "status"), "available_noninferential")

  w <- stats::setNames(numeric(length(fixef(fit))), names(fixef(fit)))
  w["trtb"] <- 1
  ctr <- contrast(fit, rbind(w))
  expect_true(all(is.finite(ctr$table$std_error)))
  expect_true(all(ctr$table$method == "wald_z_working_hessian"))

  lc <- mm_lincomb(fit, c(x = 1))
  expect_true(is.finite(lc$std_error))
  expect_identical(lc$method, "wald_z_working_hessian")

  skip_if_not_installed("generics")
  td <- generics::tidy(fit, effects = "fixed")
  expect_true(all(is.finite(td$std.error)))
})

test_that("a certified joint fit reports certified Wald, not the approximation", {
  fit <- mm_contract_fit(method = "joint_laplace")
  cap <- mixeff:::mm_glmm_inference_capability(fit)
  expect_true(cap$wald)
  expect_identical(cap$source, "certified_inference_table")
  expect_identical(cap$method, "asymptotic_wald_z")
  ct <- summary(fit, tests = "coefficients")$coefficients
  expect_true(all(is.finite(ct[, "Std. Error"])))
})

test_that("arbiter verdicts agree with the covariance wire status (WI-2.1)", {
  profiled <- mm_contract_fit()
  cap_p <- mixeff:::mm_glmm_inference_capability(profiled)
  expect_false(cap_p$wald)
  expect_true(attr(vcov(profiled), "mm_status") %in%
                c("available_noninferential", "unavailable"))

  joint <- mm_contract_fit(method = "joint_laplace")
  cap_j <- mixeff:::mm_glmm_inference_capability(joint)
  expect_true(cap_j$wald)
  expect_identical(attr(vcov(joint), "mm_status"), "available")

  # Opt-in changes capability, never the wire status: the covariance payload
  # stays noninferential on the wire; the unlock is an R-side, user-visible
  # decision (WI-2.9).
  opted <- mm_contract_fit(inference = "working_hessian")
  cap_o <- mixeff:::mm_glmm_inference_capability(opted)
  expect_true(cap_o$wald)
  expect_identical(cap_o$source, "working_hessian_opt_in")
  expect_identical(attr(vcov(opted), "mm_status"), "available_noninferential")
})

test_that("inference_options(mm_glmm) carries the route menu with the opt-in caveat", {
  fit <- mm_contract_fit()
  opts <- inference_options(fit)
  tab <- opts$table
  expect_true("asymptotic_wald_z" %in% tab$method)
  expect_true("wald_z_working_hessian" %in% tab$method)
  wh <- tab[tab$method == "wald_z_working_hessian", ]
  expect_identical(wh$expected_status, "opt_in")
  expect_match(wh$notes, "UNCERTIFIED", fixed = TRUE)
  expect_match(wh$r_verb, "working_hessian", fixed = TRUE)
  # Refusal messages point here; the pointer must never dangle.
  expect_s3_class(opts, "mm_inference_options")
  expect_output(print(opts), "wald_z_working_hessian")
})