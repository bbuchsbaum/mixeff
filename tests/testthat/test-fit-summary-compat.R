test_that("fit summaries accept known additive versions and preserve provenance", {
  for (version in c("1.0.0", "1.1.0")) {
    payload <- list(
      schema_name = "mixedmodels.fit_summary",
      schema_version = version,
      varcorr = list(),
      coefficients = list(covariance_methods = list("model_based"))
    )
    expect_identical(mm_json_parse_fit_summary(payload), payload)
  }
})

test_that("fit summaries reject unknown or malformed versions", {
  for (version in list("1.2.0", "2.0.0", "", NULL, NA_character_,
                       c("1.0.0", "1.1.0"))) {
    payload <- list(
      schema_name = "mixedmodels.fit_summary",
      schema_version = version,
      varcorr = list()
    )
    expect_error(mm_json_parse_fit_summary(payload), class = "mm_schema_error")
  }
})

test_that("fit summaries retain schema-name and required-payload validation", {
  payload <- list(
    schema_name = "mixedmodels.fit_summary",
    schema_version = "1.1.0",
    varcorr = list()
  )
  malformed <- payload
  malformed$schema_name <- "unknown.summary"
  expect_error(mm_json_parse_fit_summary(malformed), class = "mm_schema_error")
  malformed <- payload
  malformed$varcorr <- NULL
  expect_error(mm_json_parse_fit_summary(malformed), class = "mm_schema_error")
  expect_error(mm_json_parse_fit_summary(NULL), class = "mm_schema_error")
})
