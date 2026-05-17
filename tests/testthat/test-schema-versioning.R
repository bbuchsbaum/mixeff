test_that("mm_json_negotiate accepts the formula v0 schema", {
  expect_true(
    mm_json_negotiate(list(schema_name = "formula", schema_version = "v0"))
  )
})

test_that("mm_json_negotiate ignores extra header fields", {
  # Real artifact headers carry crate_version + package_version too; the
  # negotiator must accept them rather than refusing.
  expect_true(
    mm_json_negotiate(list(
      schema_name     = "formula",
      schema_version  = "v0",
      crate_version   = "0.1.0",
      package_version = "0.0.0.9000"
    ))
  )
})

test_that("mm_json_negotiate raises mm_schema_error on version mismatch", {
  err <- expect_error(
    mm_json_negotiate(list(schema_name = "formula", schema_version = "v99")),
    class = "mm_schema_error"
  )
  expect_match(conditionMessage(err), "version mismatch")
})

test_that("mm_json_negotiate raises mm_schema_error on unknown schema", {
  err <- expect_error(
    mm_json_negotiate(list(schema_name = "not_a_schema",
                            schema_version = "v0")),
    class = "mm_schema_error"
  )
  expect_match(conditionMessage(err), "unknown schema")
})

test_that("mm_json_negotiate raises mm_schema_error on malformed input", {
  expect_error(mm_json_negotiate(NULL),                         class = "mm_schema_error")
  expect_error(mm_json_negotiate("not a list"),                 class = "mm_schema_error")
  expect_error(mm_json_negotiate(list(schema_name = "formula")),
               class = "mm_schema_error")
  expect_error(mm_json_negotiate(list(schema_name = c("a", "b"),
                                       schema_version = "v0")),
               class = "mm_schema_error")
  expect_error(mm_json_negotiate(list(schema_name = NA_character_,
                                       schema_version = "v0")),
               class = "mm_schema_error")
})

test_that("mm_schema_error condition object carries the offending header", {
  bad <- list(schema_name = "formula", schema_version = "v99")
  caught <- tryCatch(
    mm_json_negotiate(bad),
    mm_schema_error = function(cnd) cnd
  )
  expect_s3_class(caught, "mm_schema_error")
  expect_s3_class(caught, "mm_condition")
  # Stored under `input`, not `header`, because rlang reserves the
  # `header` slot on conditions for cnd_header() formatting.
  expect_identical(caught$input, bad)
})

test_that("mm_json_known_schemas returns a data frame with at least 'formula'", {
  ks <- mm_json_known_schemas()
  expect_s3_class(ks, "data.frame")
  expect_named(ks, c("name", "version"), ignore.order = TRUE)
  expect_true("formula" %in% ks$name)
  expect_identical(ks$version[ks$name == "formula"], "v0")
  expect_identical(
    ks$version[ks$name == "mixedmodels.model_audit_report"], "2"
  )
  expect_identical(
    ks$version[ks$name == "mixedmodels.random_term_card"], "1"
  )
  expect_identical(
    ks$version[ks$name == "mixedmodels.fixed_effect_inference_table"], "1.0.0"
  )
  expect_identical(
    ks$version[ks$name == "mixedmodels.marginal_quantity_table"], "1.0.0"
  )
})

test_that("marginal quantity schema document exists and pins required row fields", {
  path <- system.file(
    "schemas/mixedmodels.marginal_quantity_table.schema.json",
    package = "mixeff",
    mustWork = TRUE
  )
  schema <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expect_identical(schema$properties$schema_name$const,
                   "mixedmodels.marginal_quantity_table")
  expect_identical(schema$properties$schema_version$const, "1.0.0")

  row_required <- schema$`$defs`$marginal_quantity_row$required
  expect_true(all(c(
    "quantity", "label", "estimate", "std_error", "df", "statistic",
    "rhs", "statistic_name", "p_value", "conf_low", "conf_high", "method",
    "requested_method", "status", "reliability", "estimability",
    "reason", "target", "scale", "weights", "comparison", "by",
    "specs", "grid_id", "details", "notes"
  ) %in% row_required))
})

test_that("mm_json_known_schemas content matches the manifest's schema_versions", {
  # The two surfaces must agree — if a schema appears in the manifest, the
  # negotiator must know about it.
  ks <- mm_json_known_schemas()
  manifest_schemas <- mm_formula_manifest()$schema_versions
  for (nm in names(manifest_schemas)) {
    expect_true(
      nm %in% ks$name,
      info = sprintf("manifest schema '%s' missing from negotiator's known set",
                     nm)
    )
    expect_identical(
      ks$version[ks$name == nm],
      manifest_schemas[[nm]],
      info = sprintf("schema '%s' version disagreement between manifest and negotiator",
                     nm)
    )
  }
})
