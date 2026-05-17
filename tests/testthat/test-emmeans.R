mk_emmeans_fit <- function(seed = 7L) {
  set.seed(seed)
  subject <- factor(rep(seq_len(10L), each = 6L))
  trt <- factor(rep(c("a", "b", "c"), 20L))
  x <- rep(c(-0.5, 0, 0.5, 0.25, -0.25, 0.75), 10L)
  b0 <- rnorm(nlevels(subject), sd = 0.3)
  y <- 1 +
    c(a = 0, b = 0.4, c = 1)[trt] +
    0.2 * x +
    b0[as.integer(subject)] +
    rnorm(length(subject), sd = 0.2)
  lmm(
    y ~ trt + x + (1 | subject),
    data.frame(y = y, trt = trt, x = x, subject = subject),
    control = mm_control(verbose = -1)
  )
}

test_that("optional emmeans methods reproduce native marginal mean estimates", {
  testthat::skip_if_not_installed("emmeans")
  testthat::skip_if_not_installed("estimability")
  fit <- mk_emmeans_fit()

  eg <- emmeans::emmeans(fit, ~ trt, at = list(x = 0), method = "asymptotic")
  em_df <- as.data.frame(eg)
  native <- mm_means(fit, ~ trt, at = list(x = 0), method = "asymptotic")
  trt_levels <- levels(fit$model_frame$trt)
  basis <- emm_basis.mm_lmm(
    fit,
    stats::delete.response(stats::terms(mm_fixed_formula(fit))),
    xlev = list(trt = trt_levels),
    grid = data.frame(
      trt = factor(trt_levels, levels = trt_levels),
      x = 0
    ),
    method = "asymptotic"
  )

  expect_equal(em_df$emmean, native$table$estimate, tolerance = 1e-8)
  expect_true(all(is.finite(em_df$SE)))
  expect_true(all(is.infinite(em_df$df)))
  expect_identical(attr(basis$V, "mm_method"), "model_based")
  expect_identical(attr(basis$V, "mm_status"), "available")
  expect_identical(attr(basis$V, "mm_schema_name"),
                   "mixedmodels.fixed_effect_covariance_matrix")
  expect_true(any(grepl("fixed-effect covariance from mixedmodels.fixed_effect_covariance_matrix",
                        attr(em_df, "mesg"), fixed = TRUE)))
})

test_that("optional emmeans pairwise estimates agree with fixed-effect differences", {
  testthat::skip_if_not_installed("emmeans")
  testthat::skip_if_not_installed("estimability")
  fit <- mk_emmeans_fit()

  eg <- emmeans::emmeans(fit, ~ trt, at = list(x = 0), method = "asymptotic")
  pairs <- as.data.frame(emmeans::contrast(eg, method = "pairwise"))
  beta <- fixef(fit)

  expect_equal(pairs$estimate[pairs$contrast == "a - b"],
               -unname(beta[["trt: b"]]), tolerance = 1e-8)
  expect_equal(pairs$estimate[pairs$contrast == "a - c"],
               -unname(beta[["trt: c"]]), tolerance = 1e-8)
  expect_true(all(is.finite(pairs$SE)))
})

test_that("emmeans support methods are exported for conditional registration", {
  testthat::skip_if_not_installed("emmeans")
  exported <- getNamespaceExports("mixeff")
  expect_true("recover_data.mm_lmm" %in% exported)
  expect_true("emm_basis.mm_lmm" %in% exported)
  expect_false(is.null(getS3method("recover_data", "mm_lmm",
                                   optional = TRUE,
                                   envir = asNamespace("emmeans"))))
  expect_false(is.null(getS3method("emm_basis", "mm_lmm",
                                   optional = TRUE,
                                   envir = asNamespace("emmeans"))))
})
