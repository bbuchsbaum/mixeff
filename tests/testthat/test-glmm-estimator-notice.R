# The default GLMM estimator is joint_laplace (glmer's nAGQ = 1). When the
# user does not choose a method and the request needs the profiled path
# (negative-binomial family, nAGQ > 1), glmm() must announce the switch with
# a typed mm_estimator_notice -- never swap silently.

mm_notice_data <- function() {
  set.seed(808)
  ng <- 12L
  per <- 10L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  re <- rnorm(ng, sd = 0.6)[as.integer(g)]
  y <- rbinom(n, 1, plogis(-0.2 + 0.6 * x + re))
  data.frame(y = y, x = x, g = g)
}

# Run `thunk` and report whether an mm_estimator_notice was signalled, while
# muffling it and suppressing the explain_model stdout print.
mm_fired_notice <- function(thunk) {
  fired <- FALSE
  invisible(utils::capture.output(withCallingHandlers(
    thunk(),
    mm_estimator_notice = function(cnd) {
      fired <<- TRUE
      invokeRestart("muffleMessage")
    }
  )))
  fired
}

mm_nb_notice_data <- function() {
  set.seed(809)
  g <- factor(rep(seq_len(12L), each = 10L))
  x <- rnorm(120)
  y <- rnbinom(120, size = 2, mu = exp(0.5 + 0.4 * x + rnorm(12, sd = 0.4)[g]))
  data.frame(y = y, x = x, g = g)
}

test_that("default glmm() is joint_laplace and emits no estimator notice", {
  df <- mm_notice_data()
  fit <- NULL
  expect_false(mm_fired_notice(function() {
    fit <<- glmm(y ~ x + (1 | g), df, family = binomial(),
                 control = mm_control(verbose = 0))
  }))
  expect_identical(fit$method, "joint_laplace")
})

test_that("default-method NB glmm() uses the profiled path with a notice", {
  df <- mm_nb_notice_data()
  fit <- NULL
  expect_true(mm_fired_notice(function() {
    fit <<- glmm(y ~ x + (1 | g), df, family = mm_negative_binomial(theta = 2),
                 control = mm_control(verbose = 0))
  }))
  expect_identical(fit$method, "pirls_profiled")
  # verbose = -1 silences the notice but not the resolution.
  expect_false(mm_fired_notice(function() {
    glmm(y ~ x + (1 | g), df, family = mm_negative_binomial(theta = 2),
         control = mm_control(verbose = -1))
  }))
  # An explicit joint_laplace request with NB is refused, not swapped.
  expect_error(
    glmm(y ~ x + (1 | g), df, family = mm_negative_binomial(theta = 2),
         method = "joint_laplace", control = mm_control(verbose = -1)),
    class = "mm_inference_unavailable"
  )
})

test_that("default-method nAGQ > 1 uses the profiled path with a notice", {
  df <- mm_notice_data()
  fit <- NULL
  expect_true(mm_fired_notice(function() {
    fit <<- glmm(y ~ x + (1 | g), df, family = binomial(), nAGQ = 5,
                 control = mm_control(verbose = 0))
  }))
  expect_identical(fit$method, "pirls_profiled")
  expect_identical(fit$nAGQ, 5L)
})

test_that("explicit method suppresses the estimator notice", {
  df <- mm_notice_data()
  expect_false(mm_fired_notice(function() {
    glmm(y ~ x + (1 | g), df, family = binomial(),
         method = "pirls_profiled", control = mm_control(verbose = 0))
  }))
  expect_false(mm_fired_notice(function() {
    glmm(y ~ x + (1 | g), df, family = binomial(), nAGQ = 5,
         method = "pirls_profiled", control = mm_control(verbose = 0))
  }))
})

test_that("nAGQ = 0 is lme4's fast estimate: the profiled path", {
  df <- mm_notice_data()
  fit <- glmm(y ~ x + (1 | g), df, family = binomial(), nAGQ = 0,
              control = mm_control(verbose = -1))
  expect_identical(fit$method, "pirls_profiled")
  expect_identical(fit$nAGQ_requested, 0L)
  expect_error(
    glmm(y ~ x + (1 | g), df, family = binomial(), nAGQ = 0,
         method = "joint_laplace", control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  skip_if_not_installed("lme4")
  ref <- lme4::glmer(y ~ x + (1 | g), df, family = binomial(), nAGQ = 0)
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)), tolerance = 1e-4)
  expect_equal(unname(fit$theta), unname(lme4::getME(ref, "theta")),
               tolerance = 1e-3)
})
