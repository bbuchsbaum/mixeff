# PR-4 coverage: anova()/compare() argument handling, lmerTest's ddf=
# spellings, lme4's sequential table, and the comparison print paths.

cp4c_data <- function(seed = 421L) {
  set.seed(seed)
  g <- gl(10, 6)
  x <- rnorm(60)
  z <- rnorm(60)
  f <- factor(rep(c("p", "q", "r"), 20))
  y <- 1 + 0.5 * x + 0.3 * z + c(0, 0.4, -0.3)[f] + rnorm(10, sd = 0.7)[g] +
    rnorm(60, sd = 0.5)
  data.frame(y = y, x = x, z = z, f = f, g = g)
}

cp4c_ctrl <- function() mm_control(verbose = -1)

test_that("compare() with a df method needs exactly two models", {
  d <- cp4c_data()
  a <- lmm(y ~ x + (1 | g), d, REML = FALSE, control = cp4c_ctrl())
  b <- lmm(y ~ x + z + (1 | g), d, REML = FALSE, control = cp4c_ctrl())
  c3 <- lmm(y ~ x + z + f + (1 | g), d, REML = FALSE, control = cp4c_ctrl())
  expect_error(compare(a, b, c3, method = "kenward_roger"),
               class = "mm_arg_error")
  expect_error(compare(a, b, c3, method = "satterthwaite"),
               class = "mm_arg_error")
})

test_that("anova(ddf = ) follows lmerTest's spellings", {
  d <- cp4c_data()
  fit <- lmm(y ~ x + z + f + (1 | g), d, control = cp4c_ctrl())
  expect_error(anova(fit, method = "satterthwaite", ddf = "lme4"),
               class = "mm_arg_error")
  expect_error(anova(fit, ddf = 1), class = "mm_arg_error")
  expect_error(anova(fit, ddf = NA_character_), class = "mm_arg_error")
  expect_error(anova(fit, ddf = "bogus"), class = "mm_arg_error")
  sat <- anova(fit, ddf = "Satterthwaite")
  expect_s3_class(sat, "mm_anova")
  expect_match(attr(sat, "heading"), "Satterthwaite")
  expect_identical(row.names(sat), c("x", "z", "f"))
  expect_identical(sat$NumDF, c(1, 1, 2))
  expect_s3_class(sat$table, "data.frame")
  expect_identical(sat$type, "III")
  kr <- anova(fit, ddf = "Kenward")
  expect_match(attr(kr, "heading"), "Kenward-Roger")
  skip_if_not_installed("lmerTest")
  ref <- anova(lmerTest::lmer(y ~ x + z + f + (1 | g), d), ddf = "Satterthwaite")
  expect_equal(sat[["F value"]], ref[["F value"]], tolerance = 1e-4)
  expect_equal(sat$DenDF, ref$DenDF, tolerance = 1e-3)
})

test_that("anova(ddf = \"lme4\") rebuilds lme4's sequential table", {
  d <- cp4c_data()
  fit <- lmm(y ~ x + z + f + (1 | g), d, control = cp4c_ctrl())
  seq_tab <- anova(fit, ddf = "lme4")
  expect_identical(row.names(seq_tab), c("x", "z", "f"))
  expect_identical(seq_tab$npar, c(1L, 1L, 2L))
  skip_if_not_installed("lme4")
  ref <- anova(lme4::lmer(y ~ x + z + f + (1 | g), d))
  expect_equal(seq_tab[["F value"]], ref[["F value"]], tolerance = 1e-4)
  expect_equal(seq_tab[["Sum Sq"]], ref[["Sum Sq"]], tolerance = 1e-4)
})

test_that("the sequential table refuses fits it cannot rebuild", {
  d <- cp4c_data()
  fit <- lmm(y ~ x + z + (1 | g), d, control = cp4c_ctrl())
  renamed <- fit
  names(renamed$beta)[[2]] <- "not_a_column"
  expect_error(anova(renamed, ddf = "lme4"),
               class = "mm_inference_unavailable")
  no_map <- fit
  no_map$artifact$theta_maps <- list()
  no_map$lazy_cache <- new.env(parent = emptyenv())
  expect_error(anova(no_map, ddf = "lme4"),
               class = "mm_inference_unavailable")
})

test_that("method labels and expanded single-column terms", {
  lab <- mixeff:::mm_anova_method_label
  none <- data.frame(method = NA_character_, status = "unavailable")
  expect_identical(lab(none, "asymptotic"),
                   "asymptotic Wald (chi-square/df) method")
  expect_identical(lab(none, "bootstrap"), "method `bootstrap`")
  d <- cp4c_data()
  # An expanded term with a single column needs no joint test.
  fit <- lmm(y ~ ifelse(x > 0, "hi", "lo") + z + (1 | g), d,
             control = cp4c_ctrl())
  a <- anova(fit)
  expect_identical(row.names(a), c("ifelse(x > 0, \"hi\", \"lo\")", "z"))
  expect_true(all(a$NumDF == 1))
  # Asymptotic rows are not grouped into a joint F.
  pf <- lmm(y ~ poly(x, 2) + (1 | g), d, control = cp4c_ctrl())
  aa <- anova(pf, method = "asymptotic")
  expect_identical(nrow(aa), 2L)
  expect_true(all(is.infinite(aa$DenDF)))
  sa <- anova(pf, method = "satterthwaite")
  expect_identical(row.names(sa), "poly(x, 2)")
  expect_identical(sa$NumDF, 2)
})

test_that("anova comparison objects expose and print their provenance", {
  d <- cp4c_data()
  a <- lmm(y ~ x + (1 | g), d, REML = FALSE, control = cp4c_ctrl())
  b <- lmm(y ~ x + z + (1 | g), d, REML = FALSE, control = cp4c_ctrl())
  tab <- anova(a, b)
  expect_s3_class(tab, "mm_anova_comparison")
  expect_identical(tab$npar, c(4, 5))
  expect_s3_class(tab$table, "data.frame")
  bare <- tab
  attr(bare, "mm_comparison") <- NULL
  expect_null(bare$ledger)
  out <- capture.output(print(tab))
  expect_false(any(grepl("Not certified", out)))
  cmp <- attr(tab, "mm_comparison")
  cmp$table$status[[2]] <- "unavailable"
  cmp$table$reason[[2]] <- "pretend refusal"
  attr(tab, "mm_comparison") <- cmp
  expect_output(print(tab), "row 2: unavailable -- pretend refusal")
  cmp$table$reason[[2]] <- NA_character_
  attr(tab, "mm_comparison") <- cmp
  expect_output(print(tab), "row 2: unavailable$")
})
