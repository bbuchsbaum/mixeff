# Grouping-variable coercion (lme4 parity). lmm()/glmm() coerce a
# non-categorical grouping variable (integer / numeric / logical) to a factor
# for the random-effects structure, the way lme4/nlme/glmmTMB do, instead of
# letting the native fit refuse it with "grouping factor not categorical".
# The coercion is announced via a suppressible mm_grouping_coercion_notice and
# never changes the fit. Surfaced by the OSF "Willingness to wait" reproduction
# (mote bd-01KT3ZRCKWRZQFA4W7TXTGWAZ0).

# --- helper: formula grouping-variable extraction (internal) -----------------

test_that("mm_formula_grouping_vars extracts grouping vars across RE forms", {
  expect_equal(mm_formula_grouping_vars(y ~ x + (1 | g)), "g")
  expect_equal(mm_formula_grouping_vars(y ~ x + (1 + x | ID) + (1 + x | Title)),
               c("ID", "Title"))
  expect_equal(mm_formula_grouping_vars(y ~ x + (1 | a / b)), c("a", "b"))
  expect_equal(mm_formula_grouping_vars(y ~ x + (1 + x || g)), "g")
  expect_equal(mm_formula_grouping_vars(y ~ x + (1 | s) + (1 | s:cond)),
               c("s", "cond"))
  expect_equal(mm_formula_grouping_vars(y ~ x + (x | g)), "g")  # slope-only
  expect_equal(mm_formula_grouping_vars(y ~ x), character(0))    # no RE
})

test_that("mm_coerce_grouping_factors coerces only non-categorical grouping cols", {
  d <- data.frame(
    y = 1:6, x = rnorm(6),
    int_g = c(10L, 10L, 2L, 2L, 7L, 7L),     # integer grouping  -> coerce
    num_g = c(1.5, 1.5, 2.5, 2.5, 3.5, 3.5), # double grouping    -> coerce
    chr_g = rep(c("a", "b", "c"), each = 2),  # character grouping -> leave
    fct_g = factor(rep(c("p", "q", "r"), each = 2)) # factor       -> leave
  )
  res <- mm_coerce_grouping_factors(
    y ~ x + (1 | int_g) + (1 | num_g) + (1 | chr_g) + (1 | fct_g), d)
  expect_setequal(res$coerced, c("int_g", "num_g"))
  expect_true(is.factor(res$data$int_g))
  expect_true(is.factor(res$data$num_g))
  expect_true(is.character(res$data$chr_g))           # untouched
  expect_true(is.factor(res$data$fct_g))              # already a factor
  # integer grouping factors get numerically-sorted levels, matching lme4
  expect_identical(levels(res$data$int_g), c("2", "7", "10"))
  # a grouping variable absent from the formula is never touched
  res2 <- mm_coerce_grouping_factors(y ~ x + (1 | chr_g), d)
  expect_length(res2$coerced, 0L)
})

# --- lmm() ------------------------------------------------------------------

test_that("lmm() fits an integer grouping variable (previously refused) + announces it", {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = rep(1:12, each = 5)) # g integer
  expect_true(is.integer(d$g) || is.numeric(d$g))

  # default verbose: coercion announced
  expect_message(
    utils::capture.output(lmm(y ~ x + (1 | g), d)),
    class = "mm_grouping_coercion_notice"
  )
  # verbose = -1: silent, still fits
  expect_no_message(
    fit_int <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1)),
    class = "mm_grouping_coercion_notice"
  )
  expect_s3_class(fit_int, "mm_lmm")

  # coercion does not change the fit: identical to pre-factored grouping
  d2 <- d
  d2$g <- factor(d2$g)
  fit_fac <- lmm(y ~ x + (1 | g), d2, control = mm_control(verbose = -1))
  expect_equal(unname(fixef(fit_int)), unname(fixef(fit_fac)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit_int)), as.numeric(logLik(fit_fac)),
               tolerance = 1e-8)
})

test_that("lmm() leaves a factor grouping variable alone (no notice)", {
  set.seed(2)
  d <- data.frame(y = rnorm(60), x = rnorm(60),
                  g = factor(rep(letters[1:12], each = 5)))
  expect_no_message(
    # default verbose (explain_model message fires); assert no coercion notice
    utils::capture.output(fit <- lmm(y ~ x + (1 | g), d)),
    class = "mm_grouping_coercion_notice"
  )
  expect_s3_class(fit, "mm_lmm")
})

# --- glmm() -----------------------------------------------------------------

test_that("glmm() fits an integer grouping variable + announces it; -1 silences", {
  set.seed(3)
  ng <- 12L
  d <- data.frame(
    y = rbinom(120, 1, 0.5),
    x = rnorm(120),
    g = rep(seq_len(ng), each = 10)   # integer grouping
  )
  expect_message(
    utils::capture.output(
      glmm(y ~ x + (1 | g), d, family = binomial(), method = "joint_laplace")),
    class = "mm_grouping_coercion_notice"
  )
  expect_no_message(
    fit_int <- glmm(y ~ x + (1 | g), d, family = binomial(),
                    method = "joint_laplace", control = mm_control(verbose = -1)),
    class = "mm_grouping_coercion_notice"
  )
  expect_s3_class(fit_int, "mm_glmm")

  d2 <- d
  d2$g <- factor(d2$g)
  fit_fac <- glmm(y ~ x + (1 | g), d2, family = binomial(),
                  method = "joint_laplace", control = mm_control(verbose = -1))
  expect_equal(unname(fixef(fit_int)), unname(fixef(fit_fac)), tolerance = 1e-6)
})

test_that("glmm() leaves a character grouping variable alone (no notice)", {
  set.seed(4)
  d <- data.frame(
    y = rbinom(120, 1, 0.5),
    x = rnorm(120),
    g = rep(paste0("grp", seq_len(12)), each = 10)  # character grouping
  )
  expect_true(is.character(d$g))
  expect_no_message(
    fit <- glmm(y ~ x + (1 | g), d, family = binomial(),
                method = "joint_laplace", control = mm_control(verbose = -1)),
    class = "mm_grouping_coercion_notice"
  )
  expect_s3_class(fit, "mm_glmm")
})
