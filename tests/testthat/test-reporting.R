mk_reporting_fit <- function(seed = 41L) {
  set.seed(seed)
  n_subjects <- 10L
  n_per <- 5L
  subject <- factor(rep(seq_len(n_subjects), each = n_per))
  x <- rep(seq_len(n_per) - 1, n_subjects)
  b0 <- rnorm(n_subjects, sd = 0.5)
  y <- 1.5 + 0.4 * x + b0[as.integer(subject)] + rnorm(length(x), sd = 0.25)
  lmm(
    y ~ x + (1 | subject),
    data.frame(y = y, x = x, subject = subject),
    control = mm_control(verbose = -1)
  )
}

test_that("model_report() assembles required reporting sections", {
  fit <- mk_reporting_fit()
  report <- model_report(fit)

  expect_s3_class(report, "mm_model_report")
  expect_named(report$sections, c(
    "overview",
    "model_specification",
    "data_design",
    "random_terms",
    "random_effects",
    "fixed_effects",
    "fit_statistics",
    "optimizer",
    "comparison_ledger",
    "reproducibility",
    "unavailable"
  ))
  expect_true(all(c("metadata", "sections", "unavailable", "provenance") %in%
                    names(report)))
  expect_true(nrow(report$metadata) > 0L)
  expect_true(nrow(report$provenance) == length(report$sections))
})

test_that("model_report() leaves no claim-level provenance holes", {
  fit <- mk_reporting_fit()
  report <- model_report(fit)

  for (section in names(report$sections)) {
    table <- report$sections[[section]]
    expect_true(
      all(c("source", "status", "reason") %in% names(table)),
      info = sprintf("missing provenance columns in section %s", section)
    )
    if (nrow(table)) {
      expect_false(
        any(is.na(table$source) | !nzchar(table$source)),
        info = sprintf("missing source in section %s", section)
      )
      non_available <- !table$status %in% c("available", "available_from_varcorr")
      expect_false(
        any(non_available & (is.na(table$reason) | !nzchar(table$reason))),
        info = sprintf("missing reason for caveated rows in section %s", section)
      )
    }
  }
})

test_that("model_report() unavailable ledger is explicit", {
  fit <- mk_reporting_fit()
  unavailable <- model_report(fit)$unavailable

  expect_true(all(c("section", "field", "status", "reason", "source",
                    "action_taken") %in% names(unavailable)))
  expect_false(any(is.na(unavailable$reason) | !nzchar(unavailable$reason)))
  expect_false(any(is.na(unavailable$source) | !nzchar(unavailable$source)))
  expect_true(all(unavailable$action_taken == "reported"))
})

test_that("overview and reporting_table() expose compact report fields", {
  fit <- mk_reporting_fit()
  overview <- reporting_table(fit, "overview")

  expect_identical(names(overview), c("field", "value"))
  expect_true(all(c("formula", "fit_method", "nobs", "fit_status",
                    "inference") %in% overview$field))
  expect_identical(
    reporting_table(model_report(fit), "overview"),
    overview
  )
})

test_that("reporting_table() supports compact and audit views", {
  fit <- mk_reporting_fit()
  compact <- reporting_table(fit, "fixed_effects")
  audit <- reporting_table(fit, "fixed_effects", view = "audit")

  expect_false("source" %in% names(compact))
  expect_false("details" %in% names(compact))
  expect_true(all(c("source", "reason", "details", "notes") %in% names(audit)))
  expect_true(all(names(compact) %in% names(audit)))
  expect_identical(
    reporting_table(model_report(fit), "fixed_effects", view = "audit"),
    audit
  )
})

test_that("fixed-effect report rows preserve Rust inference status", {
  fit <- mk_reporting_fit()
  fixed <- reporting_table(fit, "fixed_effects", view = "audit")
  inf <- inference_table(fit)$table

  expect_identical(fixed$method, inf$method)
  expect_identical(fixed$status, inf$status)
  expect_identical(fixed$reason, inf$reason)
  expect_true("source" %in% names(fixed))
  expect_true(any(fixed$status == "available"))
})

test_that("random-term report rows preserve Rust-authored cards", {
  fit <- mk_reporting_fit()
  random_terms <- reporting_table(fit, "random_terms", view = "audit")

  expect_true(all(c("term_id", "original_fragment", "canonical_fragment",
                    "group", "basis", "covariance", "english",
                    "design_status", "source") %in% names(random_terms)))
  expect_true(any(random_terms$group == "subject"))
  expect_true(any(nzchar(random_terms$english)))
  expect_true(all(random_terms$source %in%
                    c("random_term_cards", "cross_card_constraints")))
})

test_that("data-design report includes grouping unit counts", {
  fit <- mk_reporting_fit()
  design <- reporting_table(fit, "data_design")

  expect_true("subject" %in% design$group)
  row <- design[design$group == "subject", , drop = FALSE]
  expect_equal(row$group_levels[[1L]], 10L)
  expect_equal(row$min_rows_per_group[[1L]], 5L)
  expect_equal(row$max_rows_per_group[[1L]], 5L)
  expect_true("status" %in% names(row))
})

test_that("random-effect report includes variance rows and explicit schema gap", {
  fit <- mk_reporting_fit()
  random_effects <- reporting_table(fit, "random_effects")
  unavailable <- model_report(fit)$unavailable

  expect_true(any(random_effects$kind == "variance"))
  expect_true(any(random_effects$group == "Residual"))
  expect_true(any(unavailable$field ==
                    "stable_random_effect_variance_covariance_payload"))
  expect_true(any(unavailable$reason ==
                    "using_fit_varcorr_until_rust_report_payload_is_available"))
})

test_that("report output avoids recommendation language", {
  fit <- mk_reporting_fit()
  txt <- paste(capture.output(print(model_report(fit))), collapse = "\n")

  forbidden <- c(
    "we recommend",
    "you should",
    "try .* instead",
    "drop the random slope",
    "suggested starting model"
  )
  for (pattern in forbidden) {
    expect_false(grepl(pattern, txt, ignore.case = TRUE))
  }
})

test_that("saved fits preserve report sections backed by stored artifacts", {
  fit <- mk_reporting_fit()
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  revived <- readRDS(path)

  expect_identical(
    reporting_table(revived, "fixed_effects", view = "audit"),
    reporting_table(fit, "fixed_effects", view = "audit")
  )
  expect_identical(
    reporting_table(revived, "random_effects", view = "audit"),
    reporting_table(fit, "random_effects", view = "audit")
  )
})
