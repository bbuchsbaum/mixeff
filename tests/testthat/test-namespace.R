test_that("mixeff does not mask lme4::lmer or lme4::glmer on attach", {
  # PRD §11: attaching mixeff must not displace lme4's lmer/glmer for users
  # who have both packages loaded. We don't *export* a colliding symbol from
  # mixeff, so this is a one-line guard against accidental name reuse.
  exported <- getNamespaceExports("mixeff")
  expect_false("lmer"  %in% exported)
  expect_false("glmer" %in% exported)
  expect_true("lmm" %in% exported)
  expect_true("random_options" %in% exported)
  expect_true("compare_covariance" %in% exported)
  expect_true("fixef" %in% exported)
  expect_true("ranef" %in% exported)
  expect_true("VarCorr" %in% exported)
})

test_that("mm_lmm methods dispatch through lme4 generics when lme4 masks mixeff", {
  testthat::skip_if_not_installed("lme4")
  suppressPackageStartupMessages(suppressWarnings(library(lme4)))
  df <- data.frame(
    y = c(1.0, 1.2, 2.0, 2.3, 3.1, 3.3, 4.0, 4.1),
    x = rep(0:1, 4),
    subject = factor(rep(1:4, each = 2))
  )
  fit <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))
  expect_named(fixef(fit), c("(Intercept)", "x"))
  expect_s3_class(ranef(fit), "mm_ranef")
  expect_s3_class(VarCorr(fit), "mm_varcorr")
})

test_that("typed condition base class is documented", {
  # `mm-conditions.Rd` is generated from the roxygen block in R/conditions.R.
  # Verify the help topic is reachable in whichever mode the test runs in:
  #
  #   - devtools::test()     -> source tree; man/ has the .Rd file
  #   - R CMD check          -> installed package; help db has the topic
  #
  # The installed package does not retain man/ as a directory, so we try the
  # help db first; if that's absent (load_all-style sourcing without an
  # install), fall back to the source-tree man/ file.
  rd_db <- tryCatch(tools::Rd_db("mixeff"), error = function(e) NULL)
  if (!is.null(rd_db)) {
    expect_true(any(grepl("mm-conditions\\.Rd$", names(rd_db))))
  } else {
    rd_path <- file.path(testthat::test_path("..", ".."), "man",
                         "mm-conditions.Rd")
    expect_true(file.exists(rd_path))
  }
})
