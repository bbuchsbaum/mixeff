# Tests for lmm() data preparation: subset, na.action, contrasts.

mm_dh_data <- function() {
  set.seed(909)
  ng <- 12L; per <- 8L
  g <- factor(rep(seq_len(ng), each = per)); n <- ng * per
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
