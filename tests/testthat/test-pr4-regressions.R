# Regressions for bugs found while raising PR-4 coverage.

test_that("synthetic names terminate when a taken name prefixes the base", {
  set.seed(1)
  d <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
  d$y <- d$x + rnorm(12)[d$g] + rnorm(120)
  # `.pmax_x_0` is a prefix of `.pmax_x_0_5`: this used to loop forever.
  expect_identical(
    mixeff:::mm_synthetic_base("pmax(x, 0.5)", ".pmax_x_0"),
    ".v2_pmax_x_0_5"
  )
  fit <- lmm(y ~ pmax(x, 0) + pmax(x, 0.5) + (1 | g), d,
             control = mm_control(verbose = -1))
  expect_identical(names(fixef(fit)),
                   c("(Intercept)", "pmax(x, 0)", "pmax(x, 0.5)"))
  skip_if_not_installed("lme4")
  ref <- lme4::lmer(y ~ pmax(x, 0) + pmax(x, 0.5) + (1 | g), d)
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)), tolerance = 1e-5)
})

test_that("a one-restriction F comparison reports F = t^2 with 1 numerator df", {
  skip_if_not_installed("pbkrtest")
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4", envir = environment())
  set.seed(4)
  d <- sleepstudy
  d$z <- rnorm(nrow(d))
  big <- lmm(Reaction ~ Days + z + (Days | Subject), d,
             control = mm_control(verbose = -1))
  small <- lmm(Reaction ~ Days + (Days | Subject), d,
               control = mm_control(verbose = -1))
  B <- eval(bquote(lme4::lmer(Reaction ~ Days + z + (Days | Subject),
                              data = .(d))))
  A <- eval(bquote(lme4::lmer(Reaction ~ Days + (Days | Subject), data = .(d))))
  kr <- compare(small, big, method = "kenward_roger")
  ref <- pbkrtest::KRmodcomp(B, A)$test
  expect_equal(kr$fixed_f$num_df, 1)
  expect_equal(kr$fixed_f$statistic, ref["Ftest", "stat"], tolerance = 1e-4)
  expect_equal(kr$fixed_f$den_df, ref["Ftest", "ddf"], tolerance = 1e-4)
  expect_equal(kr$fixed_f$p_value, ref["Ftest", "p.value"], tolerance = 1e-3)
  last <- kr$table[nrow(kr$table), ]
  expect_identical(last$statistic_name, "F")
  expect_equal(last$df, 1)
  sat <- compare(small, big, method = "satterthwaite")
  sref <- pbkrtest::SATmodcomp(B, A)$test
  expect_equal(sat$fixed_f$statistic, sref$statistic, tolerance = 1e-4)
  expect_equal(sat$fixed_f$num_df, 1)
})

test_that("hatvalues() and cooks.distance() work for negative-binomial GLMMs", {
  set.seed(5)
  d <- data.frame(g = factor(rep(1:15, each = 8)), x = rnorm(120))
  d$y <- rnbinom(120, mu = exp(1 + 0.4 * d$x + rnorm(15, sd = 0.5)[d$g]),
                 size = 3)
  fit <- glmm(y ~ x + (1 | g), d, family = mm_negative_binomial(),
              control = mm_control(verbose = -1))
  h <- hatvalues(fit)
  expect_length(h, 120)
  expect_true(all(h > 0 & h < 1))
  # trace of the hat matrix lies between p and p + q
  expect_gt(sum(h), 2)
  expect_lt(sum(h), 2 + 15)
  # NB working weights: mu^2 / V(mu) with V(mu) = mu + mu^2 / theta
  parts <- mixeff:::mm_influence_parts(fit)
  mu <- parts$mu
  expect_equal(parts$working, mu^2 / (mu + mu^2 / fit$family$nb_theta))
  cd <- cooks.distance(fit)
  expect_true(all(is.finite(cd) & cd >= 0))
})

test_that("dfbeta() on a do.coef = FALSE influence result is a typed error", {
  set.seed(6)
  d <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
  d$y <- d$x + rnorm(6)[d$g] + rnorm(30)
  fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
  inf <- influence(fit, groups = "g", do.coef = FALSE)
  expect_error(dfbeta(inf), class = "mm_arg_error")
})

test_that("glmm() accepts a computed binomial response such as I(y > c)", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4", envir = environment())
  fit <- glmm(I(Reaction > 300) ~ Days + (1 | Subject), sleepstudy,
              family = binomial(), control = mm_control(verbose = -1))
  ref <- lme4::glmer(I(Reaction > 300) ~ Days + (1 | Subject), sleepstudy,
                     family = binomial())
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  expect_length(predict(fit, newdata = sleepstudy[1:3, ], type = "response"), 3)
})

test_that("anova() heading omits a Data: line holding literal data", {
  expect_identical(mixeff:::mm_anova_data_line(quote(sleepstudy)),
                   "Data: sleepstudy")
  expect_null(mixeff:::mm_anova_data_line(NULL))
  big <- as.call(c(as.name("structure"), list(as.list(seq_len(100)))))
  expect_null(mixeff:::mm_anova_data_line(big))
  expect_null(mixeff:::mm_anova_data_line(data.frame(x = 1)))
})
