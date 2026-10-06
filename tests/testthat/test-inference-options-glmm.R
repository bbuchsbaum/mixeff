route_glmm_data <- function() {
  set.seed(7)
  g <- factor(rep(1:12, each = 12))
  x <- rnorm(144)
  eta <- -0.3 + 0.8 * x + rep(rnorm(12, sd = 0.5), each = 12)
  data.frame(y = rbinom(144, 1, plogis(eta)), x = x, g = g)
}

test_that("joint-fit bootstrap execution agrees with the discovery map", {
  fit <- glmm(y ~ x + (1 | g), route_glmm_data(), family = binomial(),
               method = "joint_laplace", control = mm_control(verbose = -1))
  expect_identical(mixeff:::mm_glmm_effective_method(fit), "joint_laplace")
  route <- subset(inference_options(fit)$table, method == "glmm_parametric_bootstrap")
  expect_identical(route$expected_status, "available")
  ci <- confint(fit, method = "bootstrap", nsim = 12L, seed = 42L)
  expect_true(all(is.finite(ci)))
  expect_identical(attr(ci, "mm_estimator"), "joint_laplace")
  account <- attr(ci, "mm_bootstrap")
  expect_identical(account$requested, account$successful + account$failed)
  expect_identical(account$replicate_estimator_policy, "preserve_effective_or_fail")
  expect_equal(account$template$n_agq, fit$nAGQ)
  expect_equal(account$template_objective, fit$fit$optimizer$objective)
})

test_that("profiled and opted-in GLMM bootstrap routes execute as advertised", {
  for (inference in c("auto", "working_hessian")) {
    fit <- glmm(y ~ x + (1 | g), route_glmm_data(), family = binomial(),
                 inference = inference, control = mm_control(verbose = -1))
    route <- subset(inference_options(fit)$table, method == "glmm_parametric_bootstrap")
    expect_identical(route$expected_status, "available")
    ci <- confint(fit, method = "bootstrap", nsim = 12L, seed = 42L)
    expect_true(all(is.finite(ci)))
    expect_identical(attr(ci, "mm_estimator"), "pirls_profiled")
    acct <- attr(ci, "mm_bootstrap")
    expect_identical(acct$requested, acct$successful + acct$failed)
  }
})

test_that("a substituted fit's route follows its recorded effective estimator", {
  fit <- glmm(y ~ x + (1 | g), route_glmm_data(), family = binomial(),
               control = mm_control(verbose = -1))
  # Contract fixture: a profiled result with a recorded failed joint request.
  fit$artifact$optimizer_certificate$estimator_substitution <- list(
    requested_method = "joint_laplace", effective_method = "fast_pirls_profiled",
    requested_fit_status = "not_optimized", requested_return_code = "test_fixture"
  )
  fit$method <- "joint_laplace"
  route <- subset(inference_options(fit)$table, method == "glmm_parametric_bootstrap")
  expect_identical(route$expected_status, "available")
  ci <- confint(fit, method = "bootstrap", nsim = 12L, seed = 42L)
  expect_true(all(is.finite(ci)))
  expect_identical(attr(ci, "mm_estimator"), "pirls_profiled")
})

test_that("negative-binomial and missing-table profiled fits retain bootstrap routes", {
  d <- route_glmm_data()
  set.seed(102)
  d$y <- rnbinom(nrow(d), mu = exp(.2 + .3*d$x), size = 3)
  nb <- glmm(y ~ x + (1 | g), d, family = mm_negative_binomial(theta=3),
             control = mm_control(verbose=-1))
  legacy <- glmm(y ~ x + (1 | g), route_glmm_data(), family=binomial(),
                 control=mm_control(verbose=-1))
  legacy$artifact$fixed_effect_inference_table <- NULL
  for (fit in list(nb, legacy)) {
    route <- subset(inference_options(fit)$table, method == "glmm_parametric_bootstrap")
    expect_identical(route$expected_status, "available")
    ci <- confint(fit, method="bootstrap", nsim=12L, seed=42L)
    expect_true(all(is.finite(ci)))
    acct <- attr(ci, "mm_bootstrap")
    expect_identical(acct$requested, acct$successful + acct$failed)
    expect_error(confint(fit, method="asymptotic"), class="mm_inference_unavailable")
  }
})
