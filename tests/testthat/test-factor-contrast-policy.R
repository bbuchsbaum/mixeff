contrast_policy_data <- function() {
  set.seed(392)
  d <- expand.grid(g = factor(1:10), f = factor(c("a", "b", "c")), rep = 1:3)
  d$x <- rnorm(nrow(d))
  d$y <- 2 + c(0, 0.4, 0.9)[d$f] + rnorm(10, sd = 0.6)[d$g] +
    rnorm(nrow(d), sd = 0.3)
  d$binary <- rbinom(nrow(d), 1, plogis(-0.4 + 0.5 * d$x))
  d
}

test_that("unsupported attached contrasts are refused at all fitting entry points", {
  d <- contrast_policy_data()
  contrasts(d$f) <- contr.sum(3)
  calls <- list(
    function() compile_model(y ~ f + (1 | g), d),
    function() lmm(y ~ f + (1 | g), d, control = mm_control(verbose = -1)),
    function() glmm(binary ~ f + (1 | g), d, family = binomial(),
                    control = mm_control(verbose = -1)),
    function() compile_model(y ~ x + (f | g), d)
  )
  for (call in calls) {
    err <- tryCatch(call(), error = identity)
    expect_s3_class(err, "mm_arg_error")
    expect_identical(err$reason_code, "unsupported_factor_contrasts")
    expect_identical(err$column, "f")
  }
})

test_that("named and unnamed contrast options do not silently change a design", {
  old <- options("contrasts")
  on.exit(options(old), add = TRUE)
  d <- contrast_policy_data()
  for (opt in list(c("contr.sum", "contr.poly"),
                   c(unordered = "contr.sum", ordered = "contr.poly"))) {
    options(contrasts = opt)
    expect_error(compile_model(y ~ f + (1 | g), d), class = "mm_arg_error")
    expect_error(lmm(y ~ f + (1 | g), d, control = mm_control(verbose = -1)),
                 class = "mm_arg_error")
    expect_error(glmm(binary ~ f + (1 | g), d, family = binomial(),
                      control = mm_control(verbose = -1)), class = "mm_arg_error")
    # No factor predictor: a grouping label does not use these contrasts.
    expect_s3_class(compile_model(y ~ x + (1 | g), d), "mm_spec")
  }
})

test_that("unused and grouping-only contrast attributes do not constrain a fit", {
  d <- contrast_policy_data()
  contrasts(d$g) <- contr.sum(nlevels(d$g))
  d$unused <- d$f
  contrasts(d$unused) <- contr.sum(3)
  expect_s3_class(compile_model(y ~ f + (1 | g), d), "mm_spec")
  expect_s3_class(lmm(y ~ f + (1 | g), d, control = mm_control(verbose = -1)),
                   "mm_lmm")
  # The same factor used as a predictor must be checked even if it also groups.
  expect_error(compile_model(y ~ g + (1 | g), d), class = "mm_arg_error")
})

test_that("explicit supported coding overrides options and survives persistence", {
  old <- options("contrasts")
  on.exit(options(old), add = TRUE)
  d <- contrast_policy_data()
  expected <- model.matrix(~ f, d, contrasts.arg = list(f = "contr.treatment"))
  contrasts(d$f) <- contr.sum(3)
  options(contrasts = c("contr.sum", "contr.poly"))
  fit <- lmm(y ~ f + (1 | g), d, contrasts = list(f = "contr.treatment"),
             control = mm_control(verbose = -1))
  expect_equal(unname(model.matrix(fit)), unname(expected))
  expect_identical(colnames(model.matrix(fit)), colnames(expected))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_equal(fixef(restored), fixef(fit))
  expect_equal(predict(restored, re.form = NA), predict(fit, re.form = NA))
  # Newdata carries the old sum-coding attribute: R announces that it is
  # discarded, and predictions still use the recorded training basis.
  expect_warning(
    expect_equal(predict(restored, newdata = d, re.form = NA),
                 as.numeric(expected %*% fixef(fit)), ignore_attr = TRUE),
    "contrasts dropped from factor f", fixed = TRUE
  )
})

test_that("canonical treatment attributes work but SAS coding is not treatment coding", {
  d <- contrast_policy_data()
  contrasts(d$f) <- contr.treatment(levels(d$f))
  expect_s3_class(compile_model(y ~ f + (1 | g), d), "mm_spec")
  expect_false(identical(contr.SAS(3), contr.treatment(3)))
  expect_error(lmm(y ~ f + (1 | g), d, contrasts = list(f = "contr.SAS")),
               class = "mm_arg_error")
  expect_error(lmm(y ~ f + (1 | g), d, contrasts = list(f = "contr.sum")),
               class = "mm_arg_error")
  expect_error(lmm(y ~ f + (1 | g), d, contrasts = list(missing = "contr.treatment")),
               class = "mm_arg_error")
})
