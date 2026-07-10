# User-facing error/notice UX regression battery.
#
# The package's UX contract (never more obscure than lme4) rests on a set of
# structured conditions with plain-language, actionable messages. This file
# protects that surface against silent regression during refactors: for each
# critical path it asserts the CONDITION CLASS and a stable MESSAGE FRAGMENT
# (not the exact full string, which is allowed to evolve). Runs in the fast
# default suite.

mk_ux_lmm_data <- function(n = 60L, seed = 1L) {
  set.seed(seed)
  g <- factor(rep(seq_len(n / 4L), each = 4L))
  data.frame(y = rnorm(n), x = rnorm(n), g = g)
}

test_that("unsupported GLMM family names supported families and points to glmer", {
  d <- mk_ux_lmm_data()
  d$yb <- rbinom(nrow(d), 1, 0.4)
  err <- tryCatch(
    glmm(yb ~ x + (1 | g), d, family = inverse.gaussian(),
         control = mm_control(verbose = -1)),
    condition = function(e) e
  )
  expect_s3_class(err, "mm_inference_unavailable")
  expect_match(conditionMessage(err), "not supported", fixed = TRUE)
  expect_match(conditionMessage(err), "glmer", fixed = TRUE)
  # No raw engine-contract jargon in the user message.
  expect_false(grepl("certified upstream contract", conditionMessage(err),
                     fixed = TRUE))
})

test_that("multivariate cbind LMM response is refused in plain language", {
  d <- data.frame(y1 = rnorm(20), y2 = rnorm(20), x = rnorm(20),
                  g = factor(rep(1:5, 4)))
  err <- tryCatch(
    lmm(cbind(y1, y2) ~ x + (1 | g), d, control = mm_control(verbose = -1)),
    condition = function(e) e
  )
  expect_s3_class(err, "mm_inference_unavailable")
  expect_match(conditionMessage(err), "Multivariate responses", fixed = TRUE)
  expect_match(conditionMessage(err), "own model", fixed = TRUE)
})

test_that("no-random-effects formula error points to stats::lm()", {
  d <- mk_ux_lmm_data()
  err <- tryCatch(lmm(y ~ x, d, control = mm_control(verbose = -1)),
                  condition = function(e) e)
  expect_s3_class(err, "mm_fit_error")
  expect_match(conditionMessage(err), "not a mixed model", fixed = TRUE)
  expect_match(conditionMessage(err), "stats::lm()", fixed = TRUE)
  # The duplicated 'Caused by doTryCatch' chain must not leak.
  expect_false(grepl("doTryCatch", conditionMessage(err), fixed = TRUE))
})

test_that("NA in a model variable is a typed data error by default", {
  d <- mk_ux_lmm_data()
  d$x[3] <- NA
  err <- tryCatch(lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1)),
                  condition = function(e) e)
  expect_s3_class(err, "mm_data_error")
  expect_match(conditionMessage(err), "NA", fixed = TRUE)
})

test_that("new grouping level in predict gives R-level remedies, not Rust enums", {
  d <- mk_ux_lmm_data()
  fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
  nd <- data.frame(x = 0, g = factor("999"))
  err <- tryCatch(predict(fit, nd), condition = function(e) e)
  expect_s3_class(err, "condition")
  msg <- conditionMessage(err)
  # R-facing guidance, not the engine's NewReLevels:: enum syntax.
  expect_false(grepl("NewReLevels::", msg, fixed = TRUE))
  expect_true(grepl("re.form", msg, fixed = TRUE) ||
                grepl("allow.new.levels", msg, fixed = TRUE) ||
                grepl("new level", msg, ignore.case = TRUE))
})

test_that("GLMM default-method summary explains withheld inference plainly", {
  d <- mk_ux_lmm_data()
  d$yb <- rbinom(nrow(d), 1, 0.4)
  fit <- glmm(yb ~ x + (1 | g), d, family = binomial(),
              control = mm_control(verbose = -1))
  s <- summary(fit)
  txt <- paste(capture.output(print(s)), collapse = "\n")
  if (any(is.na(s$coefficients[["Std. Error"]]))) {
    # When SE/z/p are withheld the note must name the remedy, not the
    # working-Hessian geometry.
    expect_match(txt, "joint_laplace", fixed = TRUE)
    expect_false(grepl("active-subspace Hessian", txt, fixed = TRUE))
  } else {
    succeed("this fit certified inference; withheld-note path not exercised")
  }
})

test_that("predictor far from unit scale triggers a rescaling notice", {
  set.seed(1)
  g <- factor(rep(1:20, each = 6)); x <- runif(120, 0, 1e5)
  y <- 1 + 2e-4 * x + rnorm(20, 0, 2)[g] + rnorm(120)
  expect_message(
    lmm(y ~ x + (1 | g), data.frame(y, x, g)),
    "scale",
    class = "mm_scaling_notice"
  )
})

test_that("integer grouping variable is coerced with an announced notice", {
  d <- mk_ux_lmm_data()
  d$g <- as.integer(as.character(d$g))
  expect_message(
    lmm(y ~ x + (1 | g), d),
    class = "mm_grouping_coercion_notice"
  )
})

test_that("mm_negative_binomial rejects bad theta with a typed arg error", {
  expect_error(mm_negative_binomial(theta = -1), class = "mm_arg_error")
  expect_error(mm_negative_binomial(theta = c(1, 2)), class = "mm_arg_error")
})
