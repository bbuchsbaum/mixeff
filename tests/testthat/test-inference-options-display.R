mm_inference_options_display_fit <- function(REML = TRUE) {
  set.seed(31)
  n_subj <- 10L
  days <- as.numeric(0:4)
  b0 <- rnorm(n_subj, sd = 0.5)
  d <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
    data.frame(
      subj = factor(i),
      days = days,
      y = 1 + 0.5 * days + b0[i] + rnorm(length(days), sd = 0.3)
    )
  }))
  lmm(y ~ days + (1 | subj), d, REML = REML,
      control = mm_control(verbose = -1))
}

mm_raw_enum_like <- function(x) {
  grepl("^[a-z][a-z0-9]*(_[a-z0-9]+)+$", x)
}

test_that("inference_options display columns are populated and readable", {
  fit <- mm_inference_options_display_fit(REML = TRUE)
  opt <- inference_options(fit, nsim = 50)
  table <- opt$table

  expect_true(all(c("display_status", "display_reason", "what_to_do_next") %in%
                    names(table)))
  expect_true(all(c("asymptotic_wald_z", "satterthwaite", "kenward_roger",
                    "bootstrap", "bootstrap_lrt", "cluster_bootstrap",
                    "profile_ci") %in% table$method))

  display_cols <- c("display_status", "display_reason", "what_to_do_next")
  for (col in display_cols) {
    values <- as.character(table[[col]])
    expect_true(all(!is.na(values) & nzchar(values)),
                info = sprintf("display column `%s` must be populated", col))
    expect_false(any(values %in% table$expected_status),
                 info = sprintf("display column `%s` leaked status enums", col))
    expect_false(any(values %in% table$expected_reliability_reason),
                 info = sprintf("display column `%s` leaked reason enums", col))
  }
  expect_false(any(mm_raw_enum_like(table$display_status)))
  expect_false(any(mm_raw_enum_like(table$display_reason)))
  expect_false(any(mm_raw_enum_like(table$what_to_do_next)))
})

test_that("inference_options print uses display columns by default", {
  fit <- mm_inference_options_display_fit(REML = TRUE)
  opt <- inference_options(fit, nsim = 50)
  printed <- paste(capture.output(print(opt)), collapse = "\n")

  expect_match(printed, "runs now", fixed = TRUE)
  expect_match(printed, "refused on this fit", fixed = TRUE)
  expect_match(printed, "what_to_do_next", fixed = TRUE)
  expect_match(printed, "raw enum columns", fixed = TRUE)
})

test_that("route-table refusal reasons match the verbs they advertise", {
  fit_reml <- mm_inference_options_display_fit(REML = TRUE)
  opt_reml <- inference_options(fit_reml, nsim = 20)$table

  cluster_row <- opt_reml[opt_reml$method == "cluster_bootstrap", ,
                          drop = FALSE]
  cluster_te <- test_effect(fit_reml, "days", method = "cluster_bootstrap",
                            group = "subj")
  expect_identical(cluster_row$expected_status, "not_assessed")
  expect_identical(cluster_row$expected_reliability_reason,
                   cluster_te$table$reason_code)

  lrt_row <- opt_reml[opt_reml$method == "bootstrap_lrt", , drop = FALSE]
  lrt_te <- test_effect(fit_reml, "days", method = "bootstrap_lrt",
                        bootstrap = bootstrap_control(nsim = 20, seed = 1))
  expect_identical(lrt_row$expected_status, "not_assessed")
  expect_identical(lrt_row$expected_reliability_reason,
                   lrt_te$table$reason_code)
})

test_that("profile_ci route follows ML and REML profile contracts", {
  fit_reml <- mm_inference_options_display_fit(REML = TRUE)
  fit_ml <- mm_inference_options_display_fit(REML = FALSE)

  reml_row <- inference_options(fit_reml)$table
  reml_row <- reml_row[reml_row$method == "profile_ci", , drop = FALSE]
  expect_identical(reml_row$expected_status, "not_assessed")
  expect_identical(reml_row$expected_reliability_reason,
                   "profile_beta_unavailable_under_reml")

  ml_row <- inference_options(fit_ml)$table
  ml_row <- ml_row[ml_row$method == "profile_ci", , drop = FALSE]
  expect_identical(ml_row$expected_status, "available")
  expect_identical(ml_row$expected_reliability_reason,
                   "profile_likelihood_ci")
})
