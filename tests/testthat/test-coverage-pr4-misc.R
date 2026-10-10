# PR-4 coverage: assorted small paths -- optCtrl translation, lme4-style
# model frames, ranef() ordering fallbacks, VarCorr printing variants,
# confint() spellings and profile parm selection, GLMM drop1() with weights
# and offsets, and summary printing variants.

cp4x_ctrl <- function() mm_control(verbose = -1)

cp4x_data <- function(seed = 501L) {
  set.seed(seed)
  g <- gl(10, 6)
  x <- rnorm(60)
  w <- rep(c(1, 2, 3), 20)
  y <- 1 + 0.5 * x + rnorm(10, sd = 0.8)[g] + rnorm(60, sd = 0.5)
  data.frame(y = y, x = x, g = g, w = w)
}

test_that("mm_control(optCtrl = ) maps every tolerance and validates names", {
  ctl <- mm_control(optCtrl = list(ftol_rel = 1e-9, ftol_abs = 1e-10,
                                   xtol_rel = 1e-7))
  expect_identical(ctl$ftol_rel, 1e-9)
  expect_identical(ctl$ftol_abs, 1e-10)
  expect_identical(ctl$xtol_rel, 1e-7)
  # Explicit arguments win over optCtrl.
  ctl2 <- mm_control(ftol_rel = 1e-6, optCtrl = list(ftol_rel = 1e-9))
  expect_identical(ctl2$ftol_rel, 1e-6)
  expect_error(mm_control(optCtrl = list(1e-9)), class = "mm_arg_error")
  expect_error(mm_control(optCtrl = "maxfun"), class = "mm_arg_error")
  expect_error(mm_control(optCtrl = stats::setNames(list(1, 2), c("maxfun", ""))),
               class = "mm_arg_error")
})

test_that("lme4-style model frames carry weights/offset and survive bad formulas", {
  skip_if_not_installed("lme4")
  d <- cp4x_data()
  d$o <- runif(60)
  fit <- lmm(y ~ x + (1 | g), d, weights = w, offset = o, control = cp4x_ctrl())
  mf <- model.frame(fit)
  ref <- model.frame(lme4::lmer(y ~ x + (1 | g), d, weights = w, offset = o))
  expect_equal(mf[["(weights)"]], ref[["(weights)"]])
  expect_equal(mf[["(offset)"]], ref[["(offset)"]])
  broken <- fit
  broken$formula <- y ~ not_a_column + (1 | g)
  fr <- mixeff:::mm_lme4_model_frame(broken)
  # model.frame() fails on the stored rows, so the stored frame is used.
  expect_identical(rownames(fr), rownames(fit$model_frame))
  expect_identical(fr$x, fit$model_frame$x)
  expect_equal(fr[["(weights)"]], d$w)
})

test_that("ranef() ordering falls back when the structure is unavailable", {
  d <- cp4x_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4x_ctrl())
  re <- fit$random_effects
  no_map <- fit
  no_map$artifact$theta_maps <- list()
  no_map$lazy_cache <- new.env(parent = emptyenv())
  expect_identical(mixeff:::mm_ranef_lme4_order(no_map, re), re)
  renamed <- re
  names(renamed) <- "other"
  expect_identical(mixeff:::mm_ranef_lme4_order(fit, renamed), renamed)
  expect_identical(mixeff:::mm_ranef_lme4_order(fit, list()), list())
})

test_that("GLMM condVar degrades to an annotated result when unavailable", {
  set.seed(502)
  d <- data.frame(x = rnorm(64), g = gl(8, 8))
  d$y <- rpois(64, exp(0.3 + 0.2 * d$x + rnorm(8, sd = 0.4)[d$g]))
  fit <- glmm(y ~ x + (1 | g), d, family = poisson(), control = cp4x_ctrl())
  broken <- fit
  broken$artifact$theta_maps <- list()
  broken$lazy_cache <- new.env(parent = emptyenv())
  re <- ranef(broken, condVar = TRUE)
  expect_identical(attr(re, "mm_unavailable_reason"),
                   "random_effect_conditional_variance_unavailable")
  expect_equal(re$g[[1]], fit$random_effects$g[[1]])
})

test_that("VarCorr print variants: internal record and table boundary marks", {
  internal <- structure(
    list(table = data.frame(group = "g", sd = 1.5), residual_sd = 0.5),
    class = "mm_varcorr"
  )
  out <- capture.output(print(internal))
  expect_identical(out[[1]], "Variance components:")
  expect_true(any(grepl("Residual std. dev.: 0.5", out, fixed = TRUE)))
  none <- structure(list(table = data.frame(group = character()),
                         residual_sd = NA_real_), class = "mm_varcorr")
  expect_output(print(none), "none")
  expect_false(any(grepl("Residual", capture.output(print(none)))))

  d <- cp4x_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4x_ctrl())
  vc <- VarCorr(fit)
  attr(vc, "mm_boundary") <- NULL
  tbl <- attr(vc, "mm_table")
  tbl$boundary <- c(TRUE, rep(FALSE, nrow(tbl) - 1L))
  attr(vc, "mm_table") <- tbl
  expect_output(print(vc), "\\[boundary\\]")
  tbl$boundary <- FALSE
  attr(vc, "mm_table") <- tbl
  expect_false(any(grepl("boundary", capture.output(print(vc)))))
})

test_that("confint() accepts lme4's method spellings and profile parm forms", {
  sp <- mixeff:::mm_confint_method_spelling
  expect_identical(sp("Wald"), "wald")
  expect_identical(sp("boot"), "bootstrap")
  expect_identical(sp("profile"), "profile")
  expect_identical(sp(c("a", "b")), c("a", "b"))
  d <- cp4x_data()
  fit <- lmm(y ~ x + (1 | g), d, REML = FALSE, control = cp4x_ctrl())
  expect_equal(confint(fit, method = "Wald"), confint(fit, method = "wald"))
  pr <- confint(fit, method = "profile")
  by_index <- confint(fit, method = "profile", parm = 1:2)
  expect_identical(rownames(by_index), rownames(pr)[1:2])
  th <- confint(fit, method = "profile", parm = "theta_")
  expect_identical(rownames(th), c(".sig01", ".sigma"))
  expect_equal(unclass(th), unclass(pr)[c(".sig01", ".sigma"), ],
               ignore_attr = TRUE)
  skip_if_not_installed("lme4")
  ref <- suppressMessages(confint(lme4::lmer(y ~ x + (1 | g), d, REML = FALSE),
                                  parm = "theta_", method = "profile"))
  expect_equal(unname(unclass(th)[, 1:2]), unname(ref), tolerance = 1e-3,
               ignore_attr = TRUE)
})

test_that("profile parameter mapping reports unknown slots", {
  d <- cp4x_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4x_ctrl())
  mp <- mixeff:::mm_map_profile_parameter
  expect_identical(mp(".sig01", fit), list(name = ".sig01", kind = "sd"))
  expect_identical(mp(".sig09", fit)$kind, "unknown")
  no_map <- fit
  no_map$artifact$theta_maps <- list()
  no_map$lazy_cache <- new.env(parent = emptyenv())
  expect_identical(mp(".sig01", no_map)$kind, "unknown")
  st <- mixeff:::mm_re_structure(fit)
  expect_null(mixeff:::mm_profile_sig_slot(st, 5L))
  empty <- st
  empty$terms[[1]]$Tidx[] <- 0L
  expect_null(mixeff:::mm_profile_sig_slot(empty, 1L))
})

test_that("GLMM drop1() refits with the fit's weights and offset (vs lme4)", {
  set.seed(503)
  d <- data.frame(x = rnorm(64), z = rnorm(64), g = gl(8, 8),
                  t = runif(64, 1, 3))
  d$y <- rpois(64, d$t * exp(0.2 + 0.3 * d$x + rnorm(8, sd = 0.4)[d$g]))
  fit <- glmm(y ~ x + z + offset(log(t)) + (1 | g), d, family = poisson(),
              control = cp4x_ctrl())
  dr <- drop1(fit, test = "Chisq")
  expect_identical(dr$table$dropped, c("x", "z"))
  skip_if_not_installed("lme4")
  ref <- drop1(lme4::glmer(y ~ x + z + offset(log(t)) + (1 | g), d,
                           family = poisson()), test = "Chisq")
  expect_equal(dr$table$LRT, ref$LRT[-1], tolerance = 1e-3)
})

test_that("GLMM drop1() keeps binomial trial weights", {
  set.seed(504)
  d <- data.frame(x = rnorm(60), g = gl(10, 6), n = rep(c(4L, 7L), 30))
  d$k <- rbinom(60, d$n, plogis(0.2 + 0.5 * d$x + rnorm(10, sd = 0.4)[d$g]))
  fit <- glmm(cbind(k, n - k) ~ x + (1 | g), d, family = binomial(),
              control = cp4x_ctrl())
  expect_false(is.null(fit$weights))
  dr <- drop1(fit, test = "Chisq")
  skip_if_not_installed("lme4")
  ref <- drop1(lme4::glmer(cbind(k, n - k) ~ x + (1 | g), d,
                           family = binomial()), test = "Chisq")
  expect_equal(dr$table$LRT, ref$LRT[-1], tolerance = 1e-3)
})

test_that("small helpers: character columns, legacy bridge results", {
  df <- data.frame(a = c("b", "a"), n = 1:2, stringsAsFactors = FALSE)
  out <- mixeff:::mm_factor_character_columns(df)
  expect_identical(levels(out$a), c("a", "b"))
  expect_identical(out$n, 1:2)
  expect_identical(mixeff:::mm_bridge_fit_result("{}", function(s) list(raw = s)),
                   list(raw = "{}"))
})

test_that("summary printing: lme4 ddf table and pirls_profiled note", {
  d <- cp4x_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4x_ctrl())
  s <- summary(fit, ddf = "lme4")
  out <- capture.output(print(s))
  expect_true(any(grepl("Fixed effects", out)))
  expect_false(any(grepl("Pr\\(>\\|t\\|\\)", out)))
  set.seed(505)
  d$b <- rbinom(60, 1, plogis(0.3 * d$x + rnorm(10)[d$g]))
  gp <- glmm(b ~ x + (1 | g), d, family = binomial(), method = "pirls_profiled",
             control = cp4x_ctrl())
  expect_output(print(summary(gp)), "pirls_profiled")
})
