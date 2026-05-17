# JSON schema round-trip tests. Every artifact the bridge serializes must
# validate against a committed schema in inst/schemas/, so drift in the Rust
# contract surfaces here rather than silently misparsing downstream.
#
# Phase 1 / Stage A.4 coverage: random_term_card. Subsequent stages extend
# this file with compiled_model_artifact, model_audit_report, lmm_fit_result,
# and inference_table schemas as additional `validate_*` test_that blocks.

schema_path <- function(name) {
  candidates <- c(
    system.file("schemas", name, package = "mixeff"),
    file.path("inst", "schemas", name),
    testthat::test_path("..", "..", "inst", "schemas", name)
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    testthat::skip(sprintf("Schema `%s` is not installed", name))
  }
  candidates[[1L]]
}

mk_schema_design <- function() {
  set.seed(7L)
  df <- expand.grid(
    s = factor(seq_len(8)),
    i = factor(seq_len(4))
  )
  df$t <- rep(seq_len(4), length.out = nrow(df))
  df$y <- rnorm(nrow(df))
  df
}

cards_from <- function(formula) {
  audit <- audit_design(compile_model(formula, mk_schema_design()))
  audit$random_term_cards %||% list()
}

validate_against_random_term_card_schema <- function(card_json) {
  schema_file <- schema_path("mixedmodels.random_term_card.schema.json")
  validator <- jsonvalidate::json_validator(schema_file, engine = "ajv")
  result <- validator(card_json, verbose = TRUE, greedy = TRUE)
  if (!isTRUE(result)) {
    errors <- attr(result, "errors")
    if (is.data.frame(errors)) {
      cat(sprintf(
        "\nSchema validation errors:\n%s\n",
        paste(capture.output(print(errors)), collapse = "\n")
      ))
    }
  }
  expect_true(isTRUE(result),
              info = sprintf("random_term_card validates: %s", schema_file))
}

test_that("random_term_card JSON validates against the v1 schema for the §9.5.7 patterns", {
  skip_if_not_installed("jsonvalidate")

  patterns <- list(
    punt        = y ~ t + (1 | s),
    slope       = y ~ t + (0 + t | s),
    full        = y ~ t + (1 + t | s),
    double_bar  = y ~ t + (1 + t || s),
    split       = y ~ t + (1 | s) + (0 + t | s),
    interaction = y ~ t + (1 | s:i),
    crossed     = y ~ t + (1 | s) + (1 | i)
  )

  for (key in names(patterns)) {
    cards <- cards_from(patterns[[key]])
    expect_true(length(cards) > 0L,
                info = sprintf("pattern `%s` produced at least one card", key))
    for (card in cards) {
      card_json <- jsonlite::toJSON(card, auto_unbox = TRUE, null = "null")
      validate_against_random_term_card_schema(card_json)
    }
  }
})

test_that("random_term_card schema header is the v1 contract", {
  skip_if_not_installed("jsonlite")
  schema_file <- schema_path("mixedmodels.random_term_card.schema.json")
  schema <- jsonlite::fromJSON(schema_file, simplifyVector = FALSE)
  expect_equal(schema$properties$schema_name$const,
               "mixedmodels.random_term_card")
  expect_equal(schema$properties$schema_version$minimum, 1L)
})
