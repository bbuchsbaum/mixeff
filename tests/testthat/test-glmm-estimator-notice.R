# The default GLMM estimator (pirls_profiled) is not glmer's; glmm() must
# surface that when the user does not explicitly choose a method.

mm_notice_data <- function() {
  set.seed(808)
  ng <- 12L; per <- 10L
  g <- factor(rep(seq_len(ng), each = per)); n <- ng * per
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

test_that("default-method glmm() emits the estimator notice when verbose >= 0", {
  df <- mm_notice_data()
  expect_true(mm_fired_notice(function() {
    glmm(y ~ x + (1 | g), df, family = binomial(),
         control = mm_control(verbose = 0))
  }))
})

test_that("explicit method suppresses the estimator notice", {
  df <- mm_notice_data()
  expect_false(mm_fired_notice(function() {
    glmm(y ~ x + (1 | g), df, family = binomial(),
         method = "pirls_profiled", control = mm_control(verbose = 0))
  }))
})

test_that("verbose = -1 suppresses the estimator notice", {
  df <- mm_notice_data()
  expect_false(mm_fired_notice(function() {
    glmm(y ~ x + (1 | g), df, family = binomial(),
         control = mm_control(verbose = -1))
  }))
})

test_that("joint_laplace (explicit) does not emit the profiled-default notice", {
  df <- mm_notice_data()
  expect_false(mm_fired_notice(function() {
    glmm(y ~ x + (1 | g), df, family = binomial(),
         method = "joint_laplace",
         control = mm_control(verbose = 0, max_feval = 50000L))
  }))
})
