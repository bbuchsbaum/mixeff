# Audit1 WI-9.6: deliberate FFI containment tests. CRAN forbids compiled code
# terminating the R session; the audit requires that this be demonstrated, not
# assumed. Every test below drives an INTERNAL Rust-bridge entry point (the
# extendr wrappers in R/extendr-wrappers.R) with hostile-but-type-correct
# input and asserts the outcome is a catchable R condition -- never a crash.
#
# Entry-point inventory (R/extendr-wrappers.R) and the surfaces exercised:
#   formula/schema : mm_parse_formula (raw wrap__ symbol), mm_json_negotiate_one
#   compile/audit  : mm_compile_model_json, mm_audit_report_text
#   fit            : mm_fit_lmm_json, mm_fit_glmm_json
#   verify         : mm_verify_convergence_json, mm_verify_convergence_glmm_json
#   bootstrap      : mm_glmm_parametric_bootstrap_json,
#                    mm_fixed_effect_bootstrap_contrast_json
#   inference table: mm_fixed_effect_contrast_json
#   predict        : mm_lmm_predict_new_json
# (Not covered here: mm_formula_manifest / mm_json_known_schemas take no
# input; the remaining *_json entry points share the same bridge payload
# decoding paths as the ones tested.)
#
# Authoring evidence: each hostile call in this file was first executed in an
# isolated Rscript subprocess (18 probes, 2026-08); every probe exited rc=0
# with a catchable condition -- no R session terminated. No crash was found,
# so there is no excluded reproducer.
#
# Two deliberate scope limits:
#   * No huge-allocation inputs (no memory bombs on CRAN); hostile dimensions
#     here are absurd-but-SMALL (0 rows, negative nAGQ, nsim = 0).
#   * Wrong-TYPE arguments may be rejected R-side by extendr's coercion layer
#     before Rust runs; that is still containment and is asserted as such.
#
# Error-shape note: some refusals cross the boundary as typed mixeff
# conditions (e.g. mm_formula_error), while raw bridge errors surface as
# plain simpleError strings carrying an "mm_*_error:" routing prefix that the
# R callers translate via mm_abort_from_bridge(). Where a typed class is part
# of the contract we assert it; otherwise we assert the routing prefix (or
# just catchability, for extendr coercion errors).

## Small clean fixtures reused across payload-shaped probes -------------------

mm_ffi_lmm_fit <- function() {
  set.seed(1)
  d <- data.frame(y = rnorm(60), t = rnorm(60), s = factor(rep(1:10, each = 6)))
  lmm(y ~ t + (1 | s), d, control = mm_control(verbose = -1))
}

mm_ffi_glmm_fit <- function() {
  set.seed(2)
  d <- data.frame(y = rbinom(80, 1, 0.5), x = rnorm(80),
                  g = factor(rep(1:10, each = 8)))
  glmm(y ~ x + (1 | g), d, family = binomial(),
       control = mm_control(verbose = -1))
}

# A well-formed 20-row spec_data payload; individual tests mutate it into
# hostile shapes.
mm_ffi_spec_data <- function() {
  set.seed(3)
  list(
    column_order = c("y", "x", "g"),
    numeric_columns = list(y = rnorm(20), x = rnorm(20)),
    categorical_values = list(g = as.character(rep(1:4, each = 5))),
    categorical_levels = list(g = as.character(1:4)),
    categorical_ordered = character(0)
  )
}

## formula / schema surfaces ---------------------------------------------------

test_that("malformed formula strings are contained as typed mm_formula_error", {
  expect_error(
    mixeff:::mm_parse_formula("y ~~~ (1 |"),
    class = "mm_formula_error"
  )
})

test_that("wrong-type Robj at the raw parse-formula symbol is contained", {
  # Bypass the R-side shim entirely: hand extendr a list where its signature
  # says string. The coercion layer refuses ("Expected Strings got List");
  # that refusal IS the containment under test.
  expect_error(
    .Call(mixeff:::wrap__mm_parse_formula, list(a = 1))
  )
})

test_that("unknown schema negotiation is contained with the schema routing prefix", {
  expect_error(
    mixeff:::mm_json_negotiate_one("no.such.schema", "9.9.9"),
    regexp = "^mm_schema_error:"
  )
  # Wrong-type / zero-length Robj arguments are refused by extendr coercion.
  expect_error(mixeff:::mm_json_negotiate_one(list(), 1L))
})

## compile / audit surfaces ----------------------------------------------------

test_that("a spec_data payload with a missing column is contained at compile", {
  s <- mm_ffi_spec_data()
  s$numeric_columns$x <- NULL # column_order still names "x"
  expect_error(
    mixeff:::mm_compile_model_json(
      "y ~ x + (1 | g)", s$column_order, s$numeric_columns,
      s$categorical_values, s$categorical_levels, s$categorical_ordered
    ),
    regexp = "^mm_data_error:"
  )
})

test_that("malformed JSON handed to the audit renderer is contained", {
  expect_error(
    mixeff:::mm_audit_report_text("{ this is not json"),
    regexp = "^mm_schema_error:"
  )
})

## fit surfaces ----------------------------------------------------------------

test_that("malformed control JSON in the LMM fit primitive is contained", {
  s <- mm_ffi_spec_data()
  expect_error(
    mixeff:::mm_fit_lmm_json(
      "y ~ x + (1 | g)", TRUE, s$column_order, s$numeric_columns,
      s$categorical_values, s$categorical_levels, s$categorical_ordered,
      numeric(), "{ not json"
    ),
    regexp = "^mm_fit_error:"
  )
})

test_that("a 0-row fit payload is contained", {
  s <- list(
    column_order = c("y", "x", "g"),
    numeric_columns = list(y = numeric(0), x = numeric(0)),
    categorical_values = list(g = character(0)),
    categorical_levels = list(g = c("1", "2")),
    categorical_ordered = character(0)
  )
  # Surfaces as a plain error without the mm_* routing prefix ("Matrix index
  # out of bounds.") -- i.e. an engine-level refusal caught at the boundary
  # rather than a pre-validated one. Still a catchable condition, which is
  # the containment contract; the message shape is deliberately not pinned.
  expect_error(
    mixeff:::mm_fit_lmm_json(
      "y ~ x + (1 | g)", TRUE, s$column_order, s$numeric_columns,
      s$categorical_values, s$categorical_levels, s$categorical_ordered,
      numeric(), "{}"
    )
  )
})

test_that("non-positive nAGQ in the GLMM fit primitive is contained", {
  s <- mm_ffi_spec_data()
  s$numeric_columns$y <- as.numeric(rbinom(20, 1, 0.5))
  for (bad_nagq in c(0, -3)) {
    expect_error(
      mixeff:::mm_fit_glmm_json(
        "y ~ x + (1 | g)", "bernoulli", "logit", "joint_laplace", bad_nagq,
        s$column_order, s$numeric_columns, s$categorical_values,
        s$categorical_levels, s$categorical_ordered,
        numeric(), numeric(), "{}"
      ),
      regexp = "^mm_arg_error:"
    )
  }
})

test_that("an unknown engine family string is contained", {
  s <- mm_ffi_spec_data()
  s$numeric_columns$y <- as.numeric(rbinom(20, 1, 0.5))
  expect_error(
    mixeff:::mm_fit_glmm_json(
      "y ~ x + (1 | g)", "cauchy", "logit", "pirls_profiled", 1,
      s$column_order, s$numeric_columns, s$categorical_values,
      s$categorical_levels, s$categorical_ordered,
      numeric(), numeric(), "{}"
    ),
    regexp = "^mm_fit_error:"
  )
})

## verify-convergence surfaces -------------------------------------------------

test_that("an empty verify payload list is contained", {
  expect_error(
    mixeff:::mm_verify_convergence_json(list(), "{}"),
    regexp = "^mm_schema_error:"
  )
})

test_that("malformed verification options JSON is contained", {
  fit <- mm_ffi_lmm_fit()
  payload <- mixeff:::mm_rust_fit_bridge_payload(fit)
  payload$REML <- isTRUE(fit$REML)
  expect_error(
    mixeff:::mm_verify_convergence_json(payload, "{ nope"),
    regexp = "^mm_arg_error:"
  )
})

test_that("a negative n_agq smuggled into the GLMM verify payload is contained", {
  fit <- mm_ffi_glmm_fit()
  payload <- mixeff:::mm_rust_glmm_refit_payload(fit, "joint_laplace")
  payload$n_agq <- -1
  expect_error(
    mixeff:::mm_verify_convergence_glmm_json(payload, "{}"),
    regexp = "^mm_schema_error:"
  )
})

## bootstrap surfaces ----------------------------------------------------------

test_that("nsim = 0 in the GLMM parametric-bootstrap options is contained", {
  fit <- mm_ffi_glmm_fit()
  payload <- mixeff:::mm_rust_glmm_refit_payload(fit, "pirls_profiled")
  expect_error(
    mixeff:::mm_glmm_parametric_bootstrap_json(payload, '{"nsim":0,"seed":42}'),
    regexp = "^mm_arg_error:"
  )
})

test_that("malformed bootstrap options JSON in the contrast bootstrap is contained", {
  fit <- mm_ffi_lmm_fit()
  b <- mixeff:::mm_rust_fit_bridge_payload(fit)
  expect_error(
    mixeff:::mm_fixed_effect_bootstrap_contrast_json(
      b$formula_string, FALSE, b$spec_data$column_order,
      b$spec_data$numeric_columns, b$spec_data$categorical_values,
      b$spec_data$categorical_levels, b$spec_data$categorical_ordered,
      b$weights, b$control_json,
      c(0, 1), 1L, 2L, "r1", 0, "{ definitely not json"
    ),
    regexp = "^mm_inference_unavailable:"
  )
})

## inference-table surfaces ----------------------------------------------------

test_that("an inconsistent contrast-matrix payload is contained", {
  fit <- mm_ffi_lmm_fit()
  b <- mixeff:::mm_rust_fit_bridge_payload(fit)
  # 3 values declared as a 2x2 matrix.
  expect_error(
    mixeff:::mm_fixed_effect_contrast_json(
      b$formula_string, TRUE, b$spec_data$column_order,
      b$spec_data$numeric_columns, b$spec_data$categorical_values,
      b$spec_data$categorical_levels, b$spec_data$categorical_ordered,
      b$weights, b$control_json,
      c(0, 1, 1), 2L, 2L, c("r1", "r2"), c(0, 0), "satterthwaite"
    ),
    regexp = "^mm_inference_unavailable:"
  )
})

## predict surfaces ------------------------------------------------------------

test_that("wrong-type Robj arguments to predict_new are contained", {
  # A list where the formula string is expected: extendr coercion refuses
  # ("Expected Scalar, got List") before Rust logic runs. That R-side refusal
  # is containment.
  expect_error(
    mixeff:::mm_lmm_predict_new_json(
      list(1, 2), TRUE, "y", list(y = 1), list(), list(), character(0),
      numeric(), "{}", "y", list(y = 1), list(), list(), character(0),
      "error"
    )
  )
})
