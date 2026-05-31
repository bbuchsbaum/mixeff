mk_marginal_fit <- function(seed = 42L, weighted = FALSE) {
  set.seed(seed)
  subject <- factor(rep(seq_len(12L), each = 6L))
  trt <- factor(rep(c("a", "b", "c"), times = 24L))
  block <- factor(rep(rep(c("low", "high"), each = 3L), times = 12L))
  x <- rep(c(-0.5, 0, 0.5, 0.25, -0.25, 0.75), times = 12L)
  b0 <- rnorm(nlevels(subject), sd = 0.4)
  y <- 2 +
    c(a = 0, b = 0.5, c = 1.25)[trt] +
    c(low = -0.2, high = 0.3)[block] +
    0.25 * x +
    b0[as.integer(subject)] +
    rnorm(length(subject), sd = 0.25)
  df <- data.frame(y = y, trt = trt, block = block, x = x, subject = subject)
  if (isTRUE(weighted)) {
    df$w <- rep(c(1, 2, 0.75, 1.5, 3, 0.5), length.out = nrow(df))
    return(lmm(
      y ~ trt + block + x + (1 | subject),
      df,
      weights = w,
      control = mm_control(verbose = -1)
    ))
  }
  lmm(y ~ trt + block + x + (1 | subject), df,
      control = mm_control(verbose = -1))
}

test_that("mm_grid() builds a fixed-effect reference grid", {
  fit <- mk_marginal_fit()
  grid <- mm_grid(fit, ~ trt, at = list(x = 0))

  expect_s3_class(grid, "mm_grid")
  expect_named(grid$grid, c("trt", "x", "block"))
  expect_equal(nrow(grid$grid), 6L)
  expect_equal(ncol(grid$X), length(fixef(fit)))
  expect_identical(colnames(grid$X), names(fixef(fit)))
  expect_setequal(as.character(grid$grid$trt), levels(fit$model_frame$trt))
  expect_setequal(as.character(grid$grid$block), levels(fit$model_frame$block))
  expect_true(all(grid$grid$x == 0))
})

test_that("mm_predictions() returns contract-shaped population prediction rows", {
  fit <- mk_marginal_fit()
  grid <- mm_grid(fit, ~ trt, at = list(x = 0))
  pred <- mm_predictions(fit, grid = grid, method = "asymptotic")

  expect_s3_class(pred, "mm_marginal_quantity")
  expect_equal(nrow(pred$table), nrow(grid$grid))
  expect_true(all(c(
    "quantity", "label", "estimate", "std_error", "df", "statistic",
    "rhs", "statistic_name", "p_value", "conf_low", "conf_high", "method",
    "requested_method", "status", "reliability", "estimability",
    "reason", "target", "scale", "weights", "comparison", "by",
    "specs", "grid_id", "details", "notes"
  ) %in% names(pred$table)))
  expect_identical(unique(pred$table$quantity), "prediction")
  expect_identical(unique(pred$table$target), "population")
  expect_identical(unique(pred$table$scale), "response")
  expect_identical(unique(pred$table$comparison), "identity")
  expect_identical(pred$table$method, rep("asymptotic_wald_z", nrow(grid$grid)))
  expect_equal(length(unique(pred$table$label)), nrow(grid$grid))
  expect_true(all(is.finite(pred$table$estimate)))
  expect_true(all(is.finite(pred$table$conf_low)))
  expect_true(all(vapply(pred$table$specs, identical, logical(1), "trt")))
})

test_that("mm_means() averages nuisance fixed-factor cells through contrast()", {
  fit <- mk_marginal_fit()
  means <- mm_means(fit, ~ trt, at = list(x = 0), method = "asymptotic")
  grid <- means$grid

  expect_equal(nrow(means$table), nlevels(fit$model_frame$trt))
  expect_identical(unique(means$table$quantity), "mean")
  expect_identical(unique(means$table$weights), "equal")
  expect_true(all(means$table$status == "available"))

  for (i in seq_len(nrow(means$table))) {
    keep <- as.character(grid$grid$trt) == levels(fit$model_frame$trt)[[i]]
    expected <- colMeans(grid$X[keep, , drop = FALSE])
    expect_equal(unname(means$L[i, ]), unname(expected))
    expect_equal(means$table$estimate[[i]],
                 as.numeric(expected %*% fixef(fit)),
                 tolerance = 1e-8)
  }
})

test_that("mm_comparisons() returns pairwise differences within by groups", {
  fit <- mk_marginal_fit()
  cmp <- mm_comparisons(
    fit,
    ~ trt | block,
    at = list(x = 0),
    method = "asymptotic"
  )

  expect_s3_class(cmp, "mm_marginal_quantity")
  expect_equal(nrow(cmp$table), 6L)
  expect_identical(unique(cmp$table$quantity), "comparison")
  expect_identical(unique(cmp$table$comparison), "difference")
  expect_true(all(cmp$table$status == "available"))
  expect_true(all(grepl(" - ", cmp$table$label, fixed = TRUE)))
  expect_true(all(vapply(cmp$table$by, identical, logical(1), "block")))
  expect_true(all(vapply(cmp$table$specs, identical, logical(1), "trt")))
})

test_that("weighted LMM covariance feeds marginal standard errors", {
  fit <- mk_marginal_fit(weighted = TRUE)
  unweighted <- mk_marginal_fit()
  V <- stats::vcov(fit, type = "fixed")
  grid <- mm_grid(fit, ~ trt, at = list(x = 0))

  pred <- mm_predictions(fit, grid = grid, method = "asymptotic")
  means <- mm_means(fit, ~ trt, at = list(x = 0), method = "asymptotic")
  expected_pred_se <- sqrt(rowSums((grid$X %*% V) * grid$X))
  expected_mean_se <- sqrt(rowSums((means$L %*% V) * means$L))
  expected_weights <- rep(c(1, 2, 0.75, 1.5, 3, 0.5),
                          length.out = nrow(fit$model_frame))

  expect_equal(fit$weights, expected_weights)
  expect_identical(attr(V, "mm_status"), "available")
  expect_identical(attr(V, "mm_schema_name"),
                   "mixedmodels.fixed_effect_covariance_matrix")
  expect_gt(max(abs(unname(V) - unname(stats::vcov(unweighted, type = "fixed")))),
            1e-6)
  expect_equal(pred$table$std_error, unname(expected_pred_se), tolerance = 1e-8)
  expect_equal(means$table$std_error, unname(expected_mean_se), tolerance = 1e-8)
  expect_true(all(pred$table$status == "available"))
  expect_true(all(means$table$status == "available"))
})

test_that("unsupported marginal quantity requests are typed", {
  fit <- mk_marginal_fit()
  expect_error(mm_grid(fit, ~ missing), class = "mm_arg_error")
  expect_error(
    mm_comparisons(fit, ~ trt, comparison = "ratio"),
    class = "mm_inference_unavailable"
  )
})
