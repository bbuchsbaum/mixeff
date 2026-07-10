mk_boundary_lrt_data <- function(seed = 617L, slope = FALSE) {
  set.seed(seed)
  n_subjects <- 12L
  n_per <- 6L
  subject <- factor(rep(seq_len(n_subjects), each = n_per))
  x <- rep(seq(-1, 1, length.out = n_per), n_subjects)
  b0 <- rnorm(n_subjects, sd = 0.7)
  y <- 1 + 0.4 * x + b0[as.integer(subject)] + rnorm(length(x), sd = 0.35)
  if (slope) {
    b1 <- rnorm(n_subjects, sd = 0.25)
    y <- y + b1[as.integer(subject)] * x
  }
  data.frame(y = y, x = x, subject = subject)
}

test_that("test_random_effect() exposes certified one-component boundary LRT", {
  df <- mk_boundary_lrt_data()
  fit <- lmm(y ~ x + (1 | subject), df, REML = FALSE,
             control = mm_control(verbose = -1))

  out <- test_random_effect(fit, "subject", method = "boundary_lrt")
  row <- out$table[1L, ]

  expect_s3_class(out, "mm_random_effect_test")
  expect_identical(row$status, "available")
  expect_identical(row$method, "boundary_lrt_self_liang_mixture")
  expect_identical(row$statistic_name, "chi_bar_square")
  expect_true(is.finite(row$statistic))
  expect_true(is.finite(row$p_value))
  expect_gte(row$p_value, 0)
  expect_lte(row$p_value, 1)
  expect_identical(row$ordinary_chisq_dof, 1L)
  expect_identical(row$reference_distribution,
                   "0.5 * chi-square(0) + 0.5 * chi-square(1)")

  details <- row$details[[1L]]
  expect_length(details$mixture, 2L)
  expect_equal(vapply(details$mixture, `[[`, numeric(1), "weight"),
               c(0.5, 0.5))
  expect_true(any(grepl("Self and Liang", details$references, fixed = TRUE)))
  expect_true(any(grepl("nested ML fits", details$shape_restrictions,
                        fixed = TRUE)))

  printed <- paste(capture.output(print(out)), collapse = "\n")
  expect_match(printed, "chi-square(0)", fixed = TRUE)
  report <- reporting_table(out)
  expect_s3_class(report, "mm_reporting_table")
  expect_identical(report$table$reference_distribution[[1L]],
                   "0.5 * chi-square(0) + 0.5 * chi-square(1)")
})

test_that("boundary_lrt is refused on fixed-effect inference surfaces", {
  df <- mk_boundary_lrt_data()
  fit <- lmm(y ~ x + (1 | subject), df, REML = FALSE,
             control = mm_control(verbose = -1))

  effect <- test_effect(fit, "x", method = "boundary_lrt")
  expect_s3_class(effect, "mm_effect_test")
  expect_identical(effect$table$status[[1L]], "unsupported")
  expect_identical(effect$table$reason_code[[1L]],
                   "boundary_lrt_not_applicable_to_fixed_effects")

  con <- contrast(fit, c(0, 1), method = "boundary_lrt")
  expect_s3_class(con, "mm_contrast")
  expect_identical(con$table$status[[1L]], "unsupported")
  expect_identical(con$table$reason_code[[1L]],
                   "boundary_lrt_not_applicable_to_fixed_effects")
})

test_that("boundary_lrt refuses multi-parameter random-effect geometry", {
  df <- mk_boundary_lrt_data(slope = TRUE)
  fit <- lmm(y ~ x + (1 + x | subject), df, REML = FALSE,
             control = mm_control(verbose = -1))

  out <- test_random_effect(fit, "subject", method = "boundary_lrt")
  row <- out$table[1L, ]

  expect_identical(row$status, "not_assessed")
  expect_identical(row$reason_code,
                   "boundary_lrt_mixture_weights_not_certified")
  expect_true(is.na(row$p_value))
  expect_match(row$reason, "certifies only one added boundary", fixed = TRUE)
  expect_identical(row$theta_parameters, 3L)
})
