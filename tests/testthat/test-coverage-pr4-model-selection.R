# PR-4 coverage: compare(method = "kenward_roger") on ML fits (REML refit),
# its data guard, step() defaults/printing, and step.default().

cp4s_ctrl <- function() mm_control(verbose = -1)

cp4s_data <- function(seed = 461L) {
  set.seed(seed)
  g <- gl(12, 8)
  x <- rep(0:7, 12)
  z <- rnorm(96)
  w <- rnorm(96)
  y <- 0.3 * w + 2 + 0.6 * x + 0.4 * z + rnorm(12, sd = 1.5)[g] +
    rnorm(12, sd = 0.5)[g] * x + rnorm(96, sd = 0.6)
  data.frame(y = y, x = x, z = z, w = w, g = g)
}

test_that("Kenward-Roger comparison of ML fits refits the larger model by REML", {
  d <- cp4s_data()
  small <- lmm(y ~ x + (1 | g), d, REML = FALSE, control = cp4s_ctrl())
  big <- lmm(y ~ x + z + w + (1 | g), d, REML = FALSE, control = cp4s_ctrl())
  cmp <- compare(small, big, method = "kenward_roger")
  last <- cmp$table[nrow(cmp$table), ]
  expect_identical(last$statistic_name, "F")
  expect_match(last$reason, "REML refit")
  expect_true(cmp$fixed_f$refit_reml)
  expect_equal(cmp$fixed_f$num_df, 2)
  skip_if_not_installed("pbkrtest")
  skip_if_not_installed("lme4")
  rb <- eval(bquote(lme4::lmer(y ~ x + z + w + (1 | g), data = .(d))))
  rs <- eval(bquote(lme4::lmer(y ~ x + (1 | g), data = .(d))))
  kr <- pbkrtest::KRmodcomp(rb, rs)$test
  expect_equal(cmp$fixed_f$statistic, kr["FtestU", "stat"], tolerance = 1e-4)
  expect_equal(cmp$fixed_f$den_df, kr["FtestU", "ddf"], tolerance = 1e-3)
  expect_equal(last$p_value, kr["FtestU", "p.value"], tolerance = 1e-3)
})

test_that("F comparisons refuse models fitted to different data", {
  d <- cp4s_data()
  small <- lmm(y ~ x + (1 | g), d[-1, ], control = cp4s_ctrl())
  big <- lmm(y ~ x + z + (1 | g), d, control = cp4s_ctrl())
  expect_error(mixeff:::mm_compare_fixed_f(small, big, "kenward_roger"),
               class = "mm_inference_unavailable")
})

test_that("step() keeps significant terms and prints empty tables", {
  d <- cp4s_data()
  fit <- lmm(y ~ x + z + (1 | g) + (0 + x | g), d, control = cp4s_ctrl())
  s <- step(fit, keep = c("x", "z"))
  # Both random terms are clearly needed: nothing is eliminated.
  expect_true(all(s$random$eliminated == 0L))
  expect_identical(nrow(s$random), 2L)
  expect_null(s$fixed)
  expect_identical(deparse1(s$model$formula), deparse1(fit$formula))
  out <- capture.output(print(s))
  expect_true(any(grepl("Backward reduced random-effect table", out)))
  expect_true(any(grepl("(none tested)", out, fixed = TRUE)))
  expect_true(any(grepl("Model found", out)))

  s2 <- step(fit, reduce.random = FALSE, reduce.fixed = FALSE)
  expect_null(s2$random)
  expect_null(s2$fixed)
  out2 <- capture.output(print(s2))
  expect_identical(sum(grepl("(none tested)", out2, fixed = TRUE)), 2L)
  s3 <- step(fit, reduce.random = FALSE, keep = "x")
  expect_identical(s3$fixed$term[[1]], "z")
  expect_output(print(s3), "Backward reduced fixed-effect table \\(ddf: Satterthwaite\\)")
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    expect_identical(lmerTest::get_model(s2), s2$model)
  }
})

test_that("step.default() is stats::step()", {
  d <- cp4s_data()
  m <- lm(y ~ x + z, d)
  ours <- step(m, trace = 0)
  ref <- stats::step(m, trace = 0)
  expect_equal(coef(ours), coef(ref))
})
