test_that("mm_formula_manifest returns the documented top-level shape", {
  m <- mm_formula_manifest()
  expect_type(m, "list")
  expect_named(
    m,
    c("mixeff_rust_version", "crate_version", "schema_versions",
      "formula_features", "capabilities"),
    ignore.order = TRUE
  )
})

test_that("manifest version fields are non-empty strings", {
  m <- mm_formula_manifest()
  for (field in c("mixeff_rust_version", "crate_version")) {
    v <- m[[field]]
    expect_type(v, "character")
    expect_length(v, 1L)
    expect_true(nzchar(v))
  }
})

test_that("manifest exposes the formula schema version", {
  m <- mm_formula_manifest()
  expect_type(m$schema_versions, "list")
  expect_true("formula" %in% names(m$schema_versions))
  expect_identical(m$schema_versions$formula, "v0")
})

test_that("manifest formula_features includes the operators we parse", {
  feat <- mm_formula_manifest()$formula_features
  expect_named(
    feat,
    c("operators", "intercept_forms", "random_term_forms",
      "transformations"),
    ignore.order = TRUE
  )
  # The lme4-style operator surface — must include the bars and
  # double-bar that distinguish correlated vs uncorrelated random terms.
  for (op in c("+", "*", ":", "/", "&", "|", "||")) {
    expect_true(op %in% feat$operators,
                info = sprintf("operator %s missing from manifest", op))
  }
  expect_true("(1 | g)" %in% feat$random_term_forms)
  expect_true("(1 + x || g)" %in% feat$random_term_forms)
})

test_that("manifest capabilities accurately advertise the current phase", {
  cap <- mm_formula_manifest()$capabilities
  expect_named(
    cap,
    c("parse_formula", "compile_model", "audit_design", "explain_model",
      "random_options", "compare_covariance", "fit_lmm", "audit",
      "changes", "diagnostics", "fit_status", "parameterization",
      "roles", "as_json", "fit_glmm", "simulate", "inference",
      "fixed_effect_inference_table", "satterthwaite",
      "kenward_roger_explicit", "bootstrap_fixed_effect_payload",
      "model_comparison_table", "fit_summary_payload", "marginal_quantity_table",
      "marginal_quantities"),
    ignore.order = TRUE
  )
  # Shipped: Phase 0 (parse_formula) and Phase 1.A-F (compile/audit,
  # random-effects explanation helpers, LMM fitting, and audit verbs).
  expect_true(cap$parse_formula)
  expect_true(cap$compile_model)
  expect_true(cap$audit_design)
  expect_true(cap$explain_model)
  expect_true(cap$random_options)
  expect_true(cap$compare_covariance)
  expect_true(cap$fit_lmm)
  expect_true(cap$audit)
  expect_true(cap$changes)
  expect_true(cap$diagnostics)
  expect_true(cap$fit_status)
  expect_true(cap$parameterization)
  expect_true(cap$roles)
  expect_true(cap$as_json)
  expect_true(cap$fit_glmm)
  expect_true(cap$simulate)
  expect_true(cap$inference)
  expect_true(cap$fixed_effect_inference_table)
  expect_true(cap$satterthwaite)
  expect_true(cap$kenward_roger_explicit)
  expect_true(cap$bootstrap_fixed_effect_payload)
  expect_true(cap$model_comparison_table)
  expect_true(cap$fit_summary_payload)
  expect_true(cap$marginal_quantity_table)
  expect_true(cap$marginal_quantities)
})

test_that("manifest exposes the compiled-model-artifact schema version", {
  m <- mm_formula_manifest()
  # The negotiator and the manifest must agree — see the agreement test
  # in test-schema-versioning.R. The artifact's actual schema_name is
  # the upstream fully-qualified form `mixedmodels.compiled_model_artifact`.
  expect_true("mixedmodels.compiled_model_artifact" %in% names(m$schema_versions))
  expect_identical(m$schema_versions[["mixedmodels.compiled_model_artifact"]], "1")
})

test_that("manifest exposes the audit report and random term card schemas", {
  schemas <- mm_formula_manifest()$schema_versions
  expect_identical(schemas[["mixedmodels.model_audit_report"]], "2")
  expect_identical(schemas[["mixedmodels.random_term_card"]], "1")
  expect_identical(schemas[["mixedmodels.fixed_effect_inference_table"]], "1.0.0")
  expect_identical(schemas[["mixedmodels.model_comparison_table"]], "1.0.0")
  expect_identical(schemas[["mixedmodels.fit_summary"]], "1.0.0")
  expect_identical(schemas[["mixedmodels.marginal_quantity_table"]], "1.0.0")
})
