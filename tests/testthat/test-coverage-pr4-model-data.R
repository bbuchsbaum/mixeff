# PR-4 coverage: R-side model-data preparation (formula expansion, subset,
# na.action, offsets) and its formula helpers.

cp4m_data <- function(seed = 411L) {
  set.seed(seed)
  g <- gl(10, 6)
  x <- rnorm(60)
  f <- factor(rep(c("p", "q", "r"), 20))
  o <- runif(60, -0.5, 0.5)
  y <- 1 + 0.5 * x + o + rnorm(10, sd = 0.8)[g] + rnorm(60, sd = 0.6)
  data.frame(y = y, x = x, g = g, f = f, o = o)
}

cp4m_ctrl <- function() mm_control(verbose = -1)

test_that("bar removal handles parentheses and unary plus", {
  nb <- mixeff:::mm_nobars_term
  expect_identical(nb(quote((x + (1 | g)))), quote((x)))
  expect_null(nb(quote(((1 | g)))))
  expect_identical(nb(quote(+x)), quote(+x))
  expect_null(nb(quote(+(1 | g))))
  expect_identical(nb(quote((1 | g) + x)), quote(x))
  expect_identical(nb(quote(x + (1 | g))), quote(x))
  expect_identical(mixeff:::mm_fixed_only_formula(y ~ (1 | g)), y ~ 1)
  expect_identical(mixeff:::mm_lme4_group_specs("y ~ x"), list())
  specs <- mixeff:::mm_lme4_group_specs(y ~ x + (1 | a / b) + (x || h))
  expect_identical(vapply(specs, `[[`, "", "label"), c("b:a", "a", "h"))
  expect_identical(vapply(specs, `[[`, NA, "double"), c(FALSE, FALSE, TRUE))
})

test_that("synthetic column names avoid collisions", {
  sb <- mixeff:::mm_synthetic_base
  expect_identical(sb("!!", character()), ".term")
  expect_identical(sb("log(x)", character()), ".log_x")
  # A name in play that extends the base forces a numbered variant.
  expect_identical(sb("log(x)", ".log_x_extra"), ".log_x_2")
  expect_identical(sb("log(x)", c(".log_x_extra", ".log_x_2_b")), ".log_x_3")
  sfx <- mixeff:::mm_numeric_component_suffixes
  expect_identical(sfx(1:3), "")
  expect_identical(sfx(matrix(1:3, ncol = 1)), "")
  expect_identical(sfx(matrix(1:4, ncol = 2)), c("1", "2"))
  expect_identical(sfx(cbind(a = 1:2, b = 3:4)), c("a", "b"))
  expect_identical(mixeff:::mm_bare_term_vars(y ~ x + log(z) + (w | g)),
                   c("x", "w"))
  # A bar side terms() cannot parse contributes no bare variables.
  expect_identical(mixeff:::mm_bare_term_vars(y ~ x + (. | g)), "x")
})

test_that("character-valued and offset-only formulas match lme4", {
  skip_if_not_installed("lme4")
  d <- cp4m_data()
  fit <- lmm(y ~ ifelse(x > 0, "hi", "lo") + (1 | g), d, control = cp4m_ctrl())
  ref <- lme4::lmer(y ~ ifelse(x > 0, "hi", "lo") + (1 | g), d)
  expect_identical(names(fixef(fit)), names(lme4::fixef(ref)))
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)), tolerance = 1e-5)

  off <- lmm(y ~ offset(o) + (1 | g), d, control = cp4m_ctrl())
  roff <- lme4::lmer(y ~ offset(o) + (1 | g), d)
  expect_identical(names(fixef(off)), "(Intercept)")
  expect_equal(unname(fixef(off)), unname(lme4::fixef(roff)), tolerance = 1e-5)
  expect_equal(as.numeric(fitted(off)), as.numeric(fitted(roff)),
               tolerance = 1e-4)
})

test_that("colliding synthetic names still give lme4's coefficients", {
  skip_if_not_installed("lme4")
  d <- cp4m_data()
  d$x <- runif(60)
  fit <- lmm(y ~ pmax(x, 0.5) + pmax(x, 0) + (1 | g), d, control = cp4m_ctrl())
  ref <- lme4::lmer(y ~ pmax(x, 0.5) + pmax(x, 0) + (1 | g), d)
  expect_identical(names(fixef(fit)), names(lme4::fixef(ref)))
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)), tolerance = 1e-5)
  cols <- unlist(lapply(fit$expansion$records, `[[`, "columns"))
  expect_false(anyDuplicated(cols) > 0)
})

test_that("unsupported transform types are refused", {
  d <- cp4m_data()
  expect_error(
    lmm(y ~ as.complex(x) + (1 | g), d, control = cp4m_ctrl()),
    class = "mm_formula_error"
  )
})

test_that("data preparation validates its inputs", {
  d <- cp4m_data()
  prep <- function(formula = y ~ x + (1 | g), data = d, subset = NULL,
                   na.action = NULL, offset = NULL, lmm = TRUE) {
    mixeff:::mm_prepare_model_data(formula, data, subset, na.action,
                                   weights = NULL, offset_arg = offset,
                                   enclos = environment(), verbose = FALSE,
                                   lmm = lmm)
  }
  expect_error(prep(data = as.matrix(d[, 1:2])), class = "mm_data_error")
  expect_error(prep(formula = ~ x + (1 | g)), class = "mm_formula_error")
  expect_error(prep(formula = "y ~ x"), class = "mm_formula_error")
  expect_error(prep(offset = 1:3), class = "mm_arg_error")
  expect_error(prep(offset = c(NA, rep(0, 59))), class = "mm_arg_error")
  expect_error(prep(subset = quote(c(TRUE, FALSE))), class = "mm_arg_error")
  expect_error(prep(subset = quote("a")), class = "mm_arg_error")
  neg <- prep(subset = quote(-(1:6)))
  expect_identical(nrow(neg$data), 54L)
  expect_identical(rownames(neg$data)[[1]], "7")
  idx <- prep(subset = quote(c(1:10, NA)))
  expect_identical(nrow(idx$data), 10L)

  dn <- d
  dn$x[c(2, 5)] <- NA
  om <- prep(data = dn, na.action = "na.omit")
  expect_identical(nrow(om$data), 58L)
  expect_s3_class(om$na_action, "omit")

  # A factor response cannot carry an LMM offset.
  df <- d
  df$y <- factor(df$y > 1)
  expect_error(prep(formula = y ~ x + (1 | g), data = df, offset = d$o),
               class = "mm_data_error")
})

test_that("dropping unused levels keeps a character contrasts attribute", {
  d <- cp4m_data()
  attr(d$f, "contrasts") <- "contr.treatment"
  fit <- lmm(y ~ x + f + (1 | g), d, subset = f != "r",
             control = cp4m_ctrl())
  expect_identical(levels(fit$model_frame$f), c("p", "q"))
  expect_identical(attr(fit$model_frame$f, "contrasts"), "contr.treatment")
  expect_identical(nobs(fit), 40L)
})

test_that("negative subsets and string na.action work through lmm()", {
  skip_if_not_installed("lme4")
  d <- cp4m_data()
  d$x[3] <- NA
  fit <- lmm(y ~ x + (1 | g), d, subset = -(1:6), na.action = "na.exclude",
             control = cp4m_ctrl())
  ref <- lme4::lmer(y ~ x + (1 | g), d, subset = -(1:6),
                    na.action = na.exclude)
  expect_identical(nobs(fit), nobs(ref))
  expect_equal(unname(fixef(fit)), unname(lme4::fixef(ref)), tolerance = 1e-5)
})

test_that("expanded factor terms check newdata levels", {
  d <- cp4m_data()
  fit <- lmm(y ~ ifelse(x > 0, "hi", "lo") + (1 | g), d, control = cp4m_ctrl())
  # Data that already carries the synthetic columns is returned untouched.
  same <- mixeff:::mm_expand_newdata(fit, fit$model_frame)
  expect_identical(same, fit$model_frame)
  nd <- data.frame(x = c(-1, 1), g = factor(c("1", "2"), levels = levels(d$g)))
  p <- predict(fit, newdata = nd)
  expect_length(p, 2L)
  expect_true(p[[2]] != p[[1]])

  rec <- fit$expansion$records[[1]]
  bad <- fit
  # Pretend the training data only saw one value of the transform.
  bad$model_frame[[rec$columns]] <- factor(rep("lo", nrow(fit$model_frame)))
  expect_error(mixeff:::mm_expand_newdata(bad, nd), class = "mm_data_error")
})

test_that("random-effect predictions need the grouping variables", {
  d <- cp4m_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4m_ctrl())
  expect_error(
    mixeff:::mm_re_eta(fit, data.frame(x = 1), ~ (1 | g),
                       allow_new_levels = FALSE),
    class = "mm_data_error"
  )
  expect_identical(mixeff:::mm_rename_ranef_groups(list(a = 1), list()),
                   list(a = 1))
})
