# reporting_table() returns an mm_reporting_table object (API stabilization
# 2026-07-09); most assertions below inspect the section table itself.
rt_tbl <- function(...) reporting_table(...)$table
rt_sections <- function(...) reporting_table(...)$sections

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
  overview <- rt_tbl(fit, "overview")

  expect_identical(names(overview), c("field", "value"))
  expect_true(all(c("formula", "fit_method", "nobs", "fit_status",
                    "inference") %in% overview$field))
  expect_identical(
    rt_tbl(model_report(fit), "overview"),
    overview
  )
})

test_that("reporting_table() supports compact and audit views", {
  fit <- mk_reporting_fit()
  compact <- rt_tbl(fit, "fixed_effects")
  audit <- rt_tbl(fit, "fixed_effects", view = "audit")

  expect_false("source" %in% names(compact))
  expect_false("details" %in% names(compact))
  expect_true(all(c("source", "reason", "details", "notes") %in% names(audit)))
  expect_true(all(names(compact) %in% names(audit)))
  expect_identical(
    rt_tbl(model_report(fit), "fixed_effects", view = "audit"),
    audit
  )
})

test_that("fixed-effect report rows preserve Rust inference status", {
  fit <- mk_reporting_fit()
  fixed <- rt_tbl(fit, "fixed_effects", view = "audit")
  inf <- inference_table(fit)$table

  expect_identical(fixed$method, inf$method)
  expect_identical(fixed$status, inf$status)
  expect_identical(fixed$reason, inf$reason)
  expect_true("source" %in% names(fixed))
  expect_true(any(fixed$status == "available"))
})

test_that("random-term report rows preserve Rust-authored cards", {
  fit <- mk_reporting_fit()
  random_terms <- rt_tbl(fit, "random_terms", view = "audit")

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
  design <- rt_tbl(fit, "data_design")

  expect_true("subject" %in% design$group)
  row <- design[design$group == "subject", , drop = FALSE]
  expect_equal(row$group_levels[[1L]], 10L)
  expect_equal(row$min_rows_per_group[[1L]], 5L)
  expect_equal(row$max_rows_per_group[[1L]], 5L)
  expect_true("status" %in% names(row))
})

test_that("random-effect report consumes Rust fit-summary VarCorr payload", {
  fit <- mk_reporting_fit()
  random_effects <- rt_tbl(fit, "random_effects")
  unavailable <- model_report(fit)$unavailable

  expect_identical(fit$fit_summary$schema_name, "mixedmodels.fit_summary")
  expect_identical(fit$fit_summary$schema_version, "1.0.0")
  expect_true(any(random_effects$kind == "variance"))
  expect_true(any(random_effects$group == "Residual"))
  expect_true(all(random_effects$status == "available"))
  expect_true(all(random_effects$source == "mixedmodels.fit_summary.varcorr"))
  expect_false(any(unavailable$field ==
                     "stable_random_effect_variance_covariance_payload"))
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

test_that("reporting_table() exposes durable comparison ledgers", {
  full <- mk_reporting_fit()
  reduced <- lmm(y ~ 1 + (1 | subject), full$model_frame,
                 control = mm_control(verbose = -1))

  cmp <- compare(reduced, full)
  compact <- rt_tbl(cmp, "comparison_ledger")
  audit <- rt_tbl(cmp, "comparison_ledger", view = "audit")
  all_sections <- rt_sections(cmp, "all")

  expect_equal(nrow(audit), nrow(cmp$ledger))
  expect_true(all(c("comparison_id", "formula", "comparison_method",
                    "statistic", "p_value", "status", "reason") %in%
                    names(compact)))
  expect_false("source" %in% names(compact))
  expect_identical(all_sections$comparison_ledger, compact)
  expect_identical(audit$source, cmp$ledger$source)

  dropped <- stats::drop1(full, test = "Chisq")
  dropped_ledger <- rt_tbl(dropped, "comparison_ledger", view = "audit")

  expect_equal(nrow(dropped_ledger), nrow(dropped$ledger))
  expect_identical(dropped_ledger$reference_formula, dropped$ledger$reference_formula)
  expect_error(
    rt_tbl(cmp, "fixed_effects"),
    class = "mm_schema_error"
  )
})

test_that("saved fits preserve report sections backed by stored artifacts", {
  fit <- mk_reporting_fit()
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  revived <- readRDS(path)

  expect_identical(
    rt_tbl(revived, "fixed_effects", view = "audit"),
    rt_tbl(fit, "fixed_effects", view = "audit")
  )
  expect_identical(
    rt_tbl(revived, "random_effects", view = "audit"),
    rt_tbl(fit, "random_effects", view = "audit")
  )
})
