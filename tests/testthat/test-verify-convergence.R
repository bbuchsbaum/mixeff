# verify_convergence(): bounded engine-side verification of a fitted optimum
# (restart from the optimum, jittered restarts, opt-in optimizer consensus).
# The verdict and per-run deltas are engine-owned; these tests lock the R
# binding's shape, argument validation, and print rendering.

mk_verify_fit <- function() {
  set.seed(42)
  df <- data.frame(
    t = rep(seq_len(6), 12),
    s = factor(rep(seq_len(12), each = 6))
  )
  df$y <- 2 + 0.5 * df$t + rnorm(12)[as.integer(df$s)] * 2 + rnorm(72)
  lmm(y ~ t + (1 | s), df, control = mm_control(verbose = -1))
}

test_that("verify_convergence() returns the engine verification payload", {
  fit <- mk_verify_fit()
  v <- verify_convergence(fit)

  expect_s3_class(v, "mm_convergence_verification")
  expect_true(v$status %in% c("not_run", "restart_agrees",
                              "optimizer_consensus", "fragile", "unstable"))
  expect_false(identical(v$status, "not_run"))
  expect_true(nzchar(v$message))

  expect_s3_class(v$table, "data.frame")
  expect_named(v$table,
               c("label", "optimizer", "return_code", "objective_value",
                 "objective_delta", "theta_delta", "beta_delta", "agrees",
                 "diagnostics"))
  # default options: restart from optimum + one jittered start
  expect_gte(nrow(v$table), 2L)
  expect_true(all(nzchar(v$table$label)))
  expect_true(is.logical(v$table$agrees))

  expect_true(is.numeric(v$reference$theta))
  expect_length(v$reference$theta, length(fit$theta))
  expect_equal(v$tolerances$objective, 1e-5)
})

test_that("verify_convergence() run count follows restart/jitter options", {
  fit <- mk_verify_fit()
  v <- verify_convergence(fit, jitter_starts = 0L)
  expect_identical(nrow(v$table), 1L)
  expect_true(all(!grepl("consensus", v$table$label, ignore.case = TRUE)))
})

test_that("a well-behaved fit verifies as agreeing on restart", {
  fit <- mk_verify_fit()
  v <- verify_convergence(fit, jitter_starts = 0L)
  expect_identical(v$status, "restart_agrees")
  expect_true(all(v$table$agrees))
})

test_that("verify_convergence() validates its arguments", {
  fit <- mk_verify_fit()
  expect_error(verify_convergence(fit, jitter_starts = -1),
               class = "mm_arg_error")
  expect_error(verify_convergence(fit, jitter_scale = 0),
               class = "mm_arg_error")
  expect_error(verify_convergence(fit, max_feval = 0),
               class = "mm_arg_error")
  expect_error(verify_convergence(fit, objective_tolerance = -1),
               class = "mm_arg_error")
  expect_error(verify_convergence(fit, restart = NA),
               class = "mm_arg_error")
})

test_that("verify_convergence() refuses non-LMM input with a typed error", {
  expect_error(verify_convergence(list()), class = "mm_schema_error")
  expect_error(verify_convergence(lm(mpg ~ wt, mtcars)),
               class = "mm_schema_error")
})

test_that("print.mm_convergence_verification renders status, runs, tolerances", {
  fit <- mk_verify_fit()
  v <- verify_convergence(fit)
  printed <- paste(capture.output(print(v)), collapse = "\n")
  expect_match(printed, "Convergence verification (status:", fixed = TRUE)
  expect_match(printed, "Runs:", fixed = TRUE)
  expect_match(printed, "Tolerances:", fixed = TRUE)
})
