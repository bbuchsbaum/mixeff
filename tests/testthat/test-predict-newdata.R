# Stage C.1 (bd-01KRCKCZJ5B5AQS5BV77VMM8ZF): predict(newdata=) through the
# Rust `predict_new` contract.
#
# Done condition from the bead: predict(fit, newdata, re.form = NULL / NA)
# matches lme4::predict within ledger tolerance; parity test green.

mm_skip_if_no_lme4_local <- function() {
  testthat::skip_if_not_installed("lme4")
}

mm_sleepstudy_dataset <- function() {
  mm_skip_if_no_lme4_local()
  env <- new.env(parent = emptyenv())
  utils::data("sleepstudy", package = "lme4", envir = env)
  if (!exists("sleepstudy", envir = env, inherits = FALSE)) {
    testthat::skip("sleepstudy dataset unavailable")
  }
  get("sleepstudy", envir = env, inherits = FALSE)
}

mm_sleepstudy_pair <- function(data) {
  fit <- lmm(
    Reaction ~ Days + (1 + Days | Subject),
    data = data,
    REML = TRUE,
    control = mm_control(verbose = -1)
  )
  ref <- suppressMessages(suppressWarnings(
    lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
               data = data, REML = TRUE)
  ))
  list(fit = fit, ref = ref)
}

test_that("predict(newdata=) re.form=NULL agrees with lme4 on held-in rows", {
  data <- mm_sleepstudy_dataset()
  pair <- mm_sleepstudy_pair(data)

  obs <- stats::predict(pair$fit, newdata = data, re.form = NULL)
  ref <- stats::predict(pair$ref, newdata = data, re.form = NULL)
  expect_equal(unname(obs), unname(ref), tolerance = 1e-4)
})

test_that("predict(newdata=) re.form=NA agrees with lme4 on held-in rows", {
  data <- mm_sleepstudy_dataset()
  pair <- mm_sleepstudy_pair(data)

  obs <- stats::predict(pair$fit, newdata = data, re.form = NA)
  ref <- stats::predict(pair$ref, newdata = data, re.form = NA)
  expect_equal(unname(obs), unname(ref), tolerance = 1e-6)
})

test_that("predict(newdata=) allow.new.levels routes to NewReLevels policy", {
  data <- mm_sleepstudy_dataset()
  pair <- mm_sleepstudy_pair(data)

  # Hold out one Subject; the remaining 17 train the model and the held-out
  # rows are the unseen-level test case.
  held_out_subject <- levels(droplevels(data$Subject))[[1L]]
  train <- droplevels(data[data$Subject != held_out_subject, , drop = FALSE])
  test  <- droplevels(data[data$Subject == held_out_subject, , drop = FALSE])

  fit <- lmm(Reaction ~ Days + (1 + Days | Subject), data = train,
             REML = TRUE, control = mm_control(verbose = -1))

  # Default policy (error) must refuse unseen levels.
  expect_error(
    stats::predict(fit, newdata = test, re.form = NULL,
                   allow.new.levels = FALSE),
    class = "mm_inference_unavailable"
  )

  # Population fallback (allow.new.levels = TRUE) returns finite values
  # whose RE contribution is zero; agree with lme4(allow.new.levels=TRUE).
  obs <- stats::predict(fit, newdata = test, re.form = NULL,
                        allow.new.levels = TRUE)
  ref_fit <- suppressMessages(suppressWarnings(
    lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
               data = train, REML = TRUE)
  ))
  ref <- stats::predict(ref_fit, newdata = test, re.form = NULL,
                        allow.new.levels = TRUE)
  expect_true(all(is.finite(obs)))
  expect_equal(unname(obs), unname(ref), tolerance = 1e-4)
})

test_that("predict(newdata=) refuses unsupported re.form values", {
  data <- mm_sleepstudy_dataset()
  pair <- mm_sleepstudy_pair(data)
  expect_error(
    stats::predict(pair$fit, newdata = data,
                   re.form = stats::as.formula("~(1|Subject)")),
    class = "mm_inference_unavailable"
  )
})

test_that("predict(newdata=) refuses missing required variables", {
  data <- mm_sleepstudy_dataset()
  pair <- mm_sleepstudy_pair(data)
  bad <- data[, "Days", drop = FALSE]  # drop Subject
  expect_error(
    stats::predict(pair$fit, newdata = bad, re.form = NULL),
    class = "mm_data_error"
  )
})

test_that("predict(newdata=) accepts ~0 as the population form", {
  data <- mm_sleepstudy_dataset()
  pair <- mm_sleepstudy_pair(data)
  obs_na <- stats::predict(pair$fit, newdata = data, re.form = NA)
  obs_zero <- stats::predict(pair$fit, newdata = data,
                             re.form = stats::as.formula("~0"))
  expect_equal(unname(obs_na), unname(obs_zero), tolerance = 1e-12)
})

test_that("predict(newdata=) with an ordered fixed factor reuses the training contr.poly basis", {
  mm_skip_if_no_lme4_local()
  set.seed(21)
  group_count <- 16L
  reps <- 8L
  dat <- expand.grid(
    g = factor(seq_len(group_count)),
    o = factor(c("lo", "mid", "hi"), levels = c("lo", "mid", "hi"),
               ordered = TRUE),
    rep = seq_len(reps)
  )
  dat$y <- 1 + 0.7 * as.integer(dat$o) +
    rnorm(group_count, sd = 0.3)[as.integer(dat$g)] + rnorm(nrow(dat), sd = 0.4)
  fit <- lmm(y ~ o + (1 | g), dat, REML = TRUE, control = mm_control(verbose = -1))
  ref <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ o + (1 | g), dat, REML = TRUE)
  ))

  # Normal newdata (all ordered levels present): population + conditional match.
  nd <- dat[dat$o == "mid", , drop = FALSE]
  expect_equal(as.numeric(stats::predict(fit, newdata = nd, re.form = NA)),
               as.numeric(stats::predict(ref, newdata = nd, re.form = NA)),
               tolerance = 1e-6)
  expect_equal(as.numeric(stats::predict(fit, newdata = nd, re.form = NULL)),
               as.numeric(stats::predict(ref, newdata = nd, re.form = NULL)),
               tolerance = 1e-6)

  # Regression guard: newdata whose ordered factor declares a SINGLE level must
  # still predict, because the engine re-derives newdata contrasts from the
  # fitted training snapshot. Eagerly poly-coding the newdata column would abort
  # here (contr.poly needs >= 2 levels); the predict wrappers therefore send an
  # empty ordered-flag set for newdata.
  nd1 <- data.frame(o = factor("mid", levels = "mid", ordered = TRUE),
                    g = factor("1", levels = levels(dat$g)))
  expect_identical(nlevels(nd1$o), 1L)
  expect_equal(as.numeric(stats::predict(fit, newdata = nd1, re.form = NA)),
               as.numeric(stats::predict(ref, newdata = nd1, re.form = NA)),
               tolerance = 1e-6)
})
