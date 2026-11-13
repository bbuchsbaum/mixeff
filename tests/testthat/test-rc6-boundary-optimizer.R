# Simulated 816-row regression from mixeff-rs issue #19 (R seed 24).
# The old engine collapsed both variances and returned criterion 2313.00455.
# Independent lme4 2.0.1 reference: 2311.4456015819; tolerance allows native
# optimizer termination differences while rejecting the 1.56-unit defect.
test_that("a zero intercept variance does not collapse a supported weak slope", {
  dat <- read.csv(test_path("..", "fixtures", "boundary_collapse.csv"))
  dat$subject <- factor(dat$subject)
  reduced <- lmm(y ~ x + (0 + x | subject), dat,
                 control = mm_control(verbose = -1))
  for (f in list(y ~ x + (1 | subject) + (0 + x | subject),
                 y ~ x + (1 + x || subject),
                 y ~ x + (0 + x | subject) + (1 | subject))) {
    fit <- lmm(f, dat, control = mm_control(verbose = -1))
    criterion <- -2 * as.numeric(logLik(fit))
    expect_lt(abs(criterion - 2311.4456015819), 1e-6)
    expect_lte(criterion, -2 * as.numeric(logLik(reduced)) + 1e-6)
    expect_true(any(getME(fit, "theta") > 0.1))
    expect_true(fit_status(fit) %in% c("converged_boundary", "converged_reduced_rank"))
    expect_true(all(is.finite(summary(fit)$coefficients[, "Estimate"])))
  }
})
