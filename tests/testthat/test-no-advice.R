collect_printed_output <- function(...) {
  paste(capture.output(...), collapse = "\n")
}

expect_no_forbidden_advice <- function(output) {
  output <- tolower(output)
  forbidden <- c(
    "suggested starting model",
    "we recommend",
    "you should",
    "try .* instead",
    "drop the random slope"
  )
  for (pattern in forbidden) {
    expect_false(grepl(pattern, output, perl = TRUE),
                 info = sprintf("forbidden phrase matched: %s", pattern))
  }
}

test_that("Phase 1 verb surface keeps R9 forbidden advice phrases out", {
  set.seed(12)
  df <- data.frame(
    y = rnorm(36),
    t = rep(seq_len(6), 6),
    s = factor(rep(seq_len(6), each = 6))
  )
  spec <- compile_model(y ~ t + (1 | s), df)
  fit <- lmm(y ~ t + (1 | s), df, control = mm_control(verbose = -1))

  output <- paste(
    collect_printed_output(print(explain_model(spec))),
    collect_printed_output(print(random_options(spec, "s", "t"))),
    collect_printed_output(print(compare_covariance(spec))),
    collect_printed_output(print(audit(fit))),
    collect_printed_output(print(changes(fit))),
    collect_printed_output(print(diagnostics(fit))),
    collect_printed_output(print(parameterization(fit))),
    collect_printed_output(print(roles(fit))),
    collect_printed_output(print(summary(fit))),
    collect_printed_output(print(VarCorr(fit))),
    sep = "\n"
  )

  expect_no_forbidden_advice(output)
})
