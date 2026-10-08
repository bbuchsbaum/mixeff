# Tests for lmm() data preparation: subset, na.action, contrasts.

mm_dh_data <- function() {
  set.seed(909)
  ng <- 12L
  per <- 8L
  g <- factor(rep(seq_len(ng), each = per))
  n <- ng * per
  x <- rnorm(n)
  f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  re <- rnorm(ng, sd = 1.0)[as.integer(g)]
  y <- 1 + 0.5 * x + as.integer(f) * 0.2 + re + rnorm(n, sd = 0.7)
  data.frame(y = y, x = x, f = f, g = g)
}

test_that("subset selects rows and matches a pre-filtered fit", {
  df <- mm_dh_data()
  keep <- df$g %in% as.character(1:8)
  fit <- lmm(y ~ x + (1 | g), df, subset = g %in% as.character(1:8),
             control = mm_control(verbose = -1))
  manual <- lmm(y ~ x + (1 | g), droplevels(df[keep, ]),
                control = mm_control(verbose = -1))
  expect_equal(nobs(fit), sum(keep))
  expect_equal(unname(fixef(fit)), unname(fixef(manual)), tolerance = 1e-6)
})

test_that("subset accepts a numeric index and keeps weights aligned", {
  df <- mm_dh_data()
  w <- runif(nrow(df), 0.5, 2)
  fit <- lmm(y ~ x + (1 | g), df, weights = w, subset = 1:48,
             control = mm_control(verbose = -1))
  expect_equal(nobs(fit), 48L)
  expect_equal(length(weights(fit)), 48L)
})

test_that("na.action = na.omit drops incomplete rows; default refuses NA", {
  df <- mm_dh_data()
  df$x[c(3, 17, 40)] <- NA
  expect_error(lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1)),
               class = "mm_data_error")
  fit <- lmm(y ~ x + (1 | g), df, na.action = na.omit,
             control = mm_control(verbose = -1))
  expect_equal(nobs(fit), nrow(df) - 3L)
  manual <- lmm(y ~ x + (1 | g), na.omit(df[c("y", "x", "g")]),
                control = mm_control(verbose = -1))
  expect_equal(unname(fixef(fit)), unname(fixef(manual)), tolerance = 1e-6)
})

test_that("na.action = na.fail errors on missing data", {
  df <- mm_dh_data()
  df$y[5] <- NA
  expect_error(lmm(y ~ x + (1 | g), df, na.action = na.fail,
                   control = mm_control(verbose = -1)))
})

test_that("contrasts refuses non-treatment coding but accepts treatment", {
  df <- mm_dh_data()
  expect_error(
    lmm(y ~ f + (1 | g), df, contrasts = list(f = "contr.sum"),
        control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  expect_no_error(
    lmm(y ~ f + (1 | g), df, contrasts = list(f = "contr.treatment"),
        control = mm_control(verbose = -1))
  )
})

test_that("character fixed-effect predictors fit like factors with sorted levels (lme4)", {
  df <- mm_dh_data()
  # First-appearance order differs from sorted order, so a first-appearance
  # reference level would give different coefficients than lme4.
  df$fc <- c("c", "a", "b")[as.integer(df$f)]
  fit_chr <- lmm(y ~ x + fc + (1 | g), df, control = mm_control(verbose = -1))
  df$ff <- factor(df$fc)
  fit_fac <- lmm(y ~ x + ff + (1 | g), df, control = mm_control(verbose = -1))
  expect_identical(names(fixef(fit_chr)), c("(Intercept)", "x", "fcb", "fcc"))
  expect_equal(unname(fixef(fit_chr)), unname(fixef(fit_fac)), tolerance = 1e-8)
  # Prediction on the training rows reproduces the fitted values.
  expect_equal(unname(predict(fit_chr, newdata = df)), unname(fitted(fit_chr)),
               tolerance = 1e-8)
  if (requireNamespace("lme4", quietly = TRUE)) {
    ref <- lme4::lmer(y ~ x + fc + (1 | g), df, REML = TRUE)
    expect_equal(fixef(fit_chr), lme4::fixef(ref), tolerance = 1e-5)
  }
})

test_that("a numeric predictor whose name extends a factor's name is not mangled", {
  df <- mm_dh_data()
  set.seed(4242)
  df$group <- factor(sample(c("u", "v"), nrow(df), replace = TRUE))
  df$group_size <- rnorm(nrow(df))
  df$f_x <- rnorm(nrow(df))
  fit <- lmm(y ~ group + group_size + f + f_x + (1 | g), df,
             control = mm_control(verbose = -1))
  expect_identical(
    names(fixef(fit)),
    c("(Intercept)", "groupv", "group_size", "fb", "fc", "f_x")
  )
})

test_that("contr.SAS is refused: the engine always uses the first level as reference", {
  df <- mm_dh_data()
  expect_error(
    lmm(y ~ f + (1 | g), df, contrasts = list(f = "contr.SAS"),
        control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
})

test_that("y ~ 0 + f estimates one mean per level, matching lme4 and y ~ f", {
  df <- mm_dh_data()
  cell <- lmm(y ~ 0 + f + x + (1 | g), df, REML = FALSE,
              control = mm_control(verbose = -1))
  base <- lmm(y ~ f + x + (1 | g), df, REML = FALSE,
              control = mm_control(verbose = -1))
  expect_identical(names(fixef(cell)), c("fa", "fb", "fc", "x"))
  expect_equal(as.numeric(logLik(cell)), as.numeric(logLik(base)),
               tolerance = 1e-6)
  b <- fixef(base)
  expect_equal(unname(fixef(cell)[c("fa", "fb", "fc")]),
               unname(b[["(Intercept)"]] + c(0, b[["fb"]], b[["fc"]])),
               tolerance = 1e-5)
  if (requireNamespace("lme4", quietly = TRUE)) {
    ref <- lme4::lmer(y ~ 0 + f + x + (1 | g), df, REML = FALSE)
    expect_equal(fixef(cell), lme4::fixef(ref), tolerance = 1e-5)
  }
})

test_that("predict(newdata) reproduces fitted values for nested fixed effects", {
  df <- mm_dh_data()
  set.seed(31)
  df$h <- factor(sample(c("u", "v"), nrow(df), replace = TRUE))
  fit <- lmm(y ~ f / h + (1 | g), df, control = mm_control(verbose = -1))
  expect_equal(unname(predict(fit, newdata = df)), unname(fitted(fit)),
               tolerance = 1e-8)
})

test_that("simulate() scales residual noise by prior weights (sigma / sqrt(w))", {
  df <- mm_dh_data()
  w <- rep(c(1, 100), length.out = nrow(df))
  fit <- lmm(y ~ x + (1 | g), df, weights = w,
             control = mm_control(verbose = -1))
  sims <- as.matrix(simulate(fit, nsim = 200, seed = 1, re.form = NA))
  resid_sd <- apply(sims - fit$fixed_fitted, 1, stats::sd)
  heavy <- w == 100
  expect_gt(mean(resid_sd[!heavy]), 5 * mean(resid_sd[heavy]))
})

test_that("simulate() keeps a zero-variance crossed term at zero", {
  set.seed(2)
  d <- expand.grid(subj = factor(1:15), item = factor(1:8), rep = 1:2)
  d$y <- rnorm(15, sd = 2)[as.integer(d$subj)] + rnorm(nrow(d), sd = 0.5)
  fit <- lmm(y ~ 1 + (1 | subj) + (1 | item), d,
             control = mm_control(verbose = -1))
  vc <- VarCorr(fit)$table
  skip_if_not(any(vc$group == "item" & vc$variance == 0),
              "item variance is not exactly zero for this draw")
  # Before the fix the zero SD was treated as "missing" and replaced by the
  # first "(Intercept)" variance in VarCorr -- the subject variance.
  terms <- fit$artifact$semantic_model$random_terms
  labels <- vapply(seq_along(terms), function(i) {
    mixeff:::mm_random_term_group_label(fit, terms[[i]], i)
  }, character(1))
  i <- match("item", labels)
  Sigma <- mixeff:::mm_random_term_covariance(
    fit, if (is.null(terms[[i]]$id)) sprintf("r%d", i - 1L) else terms[[i]]$id,
    "(Intercept)", "item"
  )
  expect_equal(Sigma[1, 1], 0)
})
