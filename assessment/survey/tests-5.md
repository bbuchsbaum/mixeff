# Test Survey: tests-5 (predict-newdata, ranef-condvar, phase2-revive, schema-versioning)

Survey date: 2026-05-31
Files:
- tests/testthat/test-predict-newdata.R
- tests/testthat/test-ranef-condvar.R
- tests/testthat/test-phase2-revive.R
- tests/testthat/test-schema-versioning.R

---

## test-predict-newdata.R

### What is tested

Six `test_that` blocks all exercised on `sleepstudy` + random-slope model
`Reaction ~ Days + (1 + Days | Subject)`.

| # | Test | Assertion / tolerance |
|---|------|-----------------------|
| 1 | `predict(newdata=, re.form=NULL)` (held-in rows) agrees with lme4 | `expect_equal(..., tolerance=1e-4)` |
| 2 | `predict(newdata=, re.form=NA)` (fixed-only, held-in) agrees with lme4 | `expect_equal(..., tolerance=1e-6)` |
| 3 | `allow.new.levels=FALSE` refuses unseen grouping levels with `mm_inference_unavailable`; `allow.new.levels=TRUE` returns finite values matching lme4 on held-out subject | `expect_error(..., class="mm_inference_unavailable")`, `expect_equal(..., tolerance=1e-4)` |
| 4 | Unsupported `re.form` formula (e.g. `~(1|Subject)`) raises `mm_inference_unavailable` | `expect_error(..., class="mm_inference_unavailable")` |
| 5 | Missing required variable in `newdata` raises `mm_data_error` | `expect_error(..., class="mm_data_error")` |
| 6 | `re.form = ~0` is equivalent to `re.form = NA` | `expect_equal(..., tolerance=1e-12)` |

### Skips

- `mm_skip_if_no_lme4_local()` (= `skip_if_not_installed("lme4")`) guards all tests.
- One additional runtime skip if the `sleepstudy` data object is unavailable.

### What is NOT tested (gaps)

- **In-sample predict with no `newdata`**: `predict(fit)` (uses cached fitted values) is never exercised from this file. Whether it returns the same values as lme4 for in-sample predictions is untested here.
- **`se.fit = TRUE` stub**: `predict(..., se.fit=TRUE)` should return a list with `fit` and an all-`NA` `se.fit` vector plus an `mm_unavailable_reason` attribute. The attribute contract is untested.
- **`interval` guard**: `predict(..., interval="confidence")` should raise `mm_inference_unavailable`. Not tested.
- **Column-count mismatch in fixed-only path**: `mm_predict_fixed_only` has an explicit guard for when the design matrix width from `newdata` doesn't match the stored coefficient count. No test exercises this code path.
- **Multi-group models** (`(1|g1) + (1|g2)`): `newdata` prediction with two random grouping factors is untested; the Rust contract expects group-level vectors for all groups.
- **Single-level (intercept-only RE) models**: All tests use a random-slope model. An intercept-only model `y ~ x + (1|subject)` is never tested in this file.
- **`re.form=NULL` with `allow.new.levels=TRUE` parity**: The existing test checks that predictions are finite and match lme4, but it does not separately verify that the RE contribution is exactly zero (population fallback semantics).
- **Factor levels in `newdata` that are valid but reordered**: Whether `predict` tolerates a different level ordering in a factor column in `newdata` vs. the training data is not checked.
- **GLMM predict stub**: `predict.mm_glmm` exists and unconditionally raises `mm_inference_unavailable`. No test asserts this stub fires with the expected condition class.

---

## test-ranef-condvar.R

### What is tested

Four `test_that` blocks, all on `sleepstudy` with model
`Reaction ~ Days + (1 + Days | Subject)`.

| # | Test | Assertion / tolerance |
|---|------|-----------------------|
| 1 | `ranef(condVar=TRUE)` returns a finite, symmetric, PSD 3-D `postVar` array with correct dimensions; no `mm_unavailable_reason` on success | structural + `tolerance=1e-9` symmetry |
| 2 | `postVar` values agree with lme4 | `expect_equal(..., tolerance=1e-3)` |
| 3 | Repeated calls return identical array (caching); `cond_var` key exists in `lazy_cache` | `expect_identical` |
| 4 | `ranef(condVar=FALSE)` (default) returns no `postVar` and no `mm_unavailable_reason` | `expect_null` on both attributes |

### Skips

- Same `skip_if_not_installed("lme4")` + `sleepstudy` runtime guard.

### What is NOT tested (gaps)

- **Intercept-only model** (`y ~ x + (1|subject)`): `condVar=TRUE` on a scalar random effect. The `postVar` array should be 1×1×n. This is the most common real-world model structure and is completely untested.
- **Multiple grouping factors** (`(1|g1) + (1|g2)`): `condVar=TRUE` when two separate groups exist. The R bridge has explicit code for merging same-group RE terms and handling group-level mismatch; neither path is tested.
- **`(0 + x | g)` (slope-only, no intercept) model**: `postVar` should be 1×1×n but with slope variance rather than intercept variance. Untested.
- **GLMM `ranef(condVar=TRUE)` stub**: `ranef.mm_glmm(condVar=TRUE)` attaches an all-`NA` `postVar` with `mm_cond_var_error` attribute. Neither the NA structure nor the `mm_cond_var_error` attribute is verified by any test in this file.
- **`mm_cond_var_error` attribute on Rust bridge failure**: If `wrap__mm_lmm_cond_var_json` returns an error, the fallback path in `mm_attach_ranef_postvar_unavailable` should attach `mm_cond_var_error`. This error path is not simulated in the test suite.
- **Cache invalidation / fresh object**: The caching test verifies identity on the same object, but there is no test that a newly revived (deserialized) fit does not carry a stale `cond_var` cache entry.
- **Boundary/singular model condvar**: A singular model (zero RE variance, e.g. variance collapses to boundary) is not tested; the `postVar` diagonal could be exactly zero in such cases.

---

## test-phase2-revive.R

### What is tested

Eight `test_that` blocks using a small synthetic dataset (`n=8` subjects, `n_per=5`), model `y ~ x + (1|subject)`.

| # | Test | Coverage |
|---|------|---------|
| 1 | Fit carries `lazy_cache` env, correct `schema` fields, and `schema_name` value | structural identity |
| 2 | `model.matrix(type="fixed"/"random")`, `getME("Lambda"/"theta"/"beta"/"flist"/"cnms"/"y")` return correct types and dimensions; `X` and `Z` cached | dim checks + `expect_s4_class` + cache existence |
| 3 | `random_blocks`, `optimizer_certificate`, `inference_table`, `reproducibility`, `is_singular` all return correct S3 classes and field names; `inf$table$method` restricted to known values | structural + class |
| 4 | `vcov()` returns symmetric matrix with correct diagonal (`std_errors^2`), correct `mm_schema_*` attributes, no `mm_unavailable_reason` | `tolerance=1e-12` + attribute identity |
| 5 | `vcov()` fallback (missing `fixed_effect_vcov` and artifact entry) returns diagonal-only matrix with `mm_status="unavailable"` | attribute identity |
| 6 | Rank-deficient fit (`x` aliased by `x2`) returns all-`NA` vcov with correct `mm_unavailable_reason` matching `"rank_deficient"`, `details$rank==2`, `details$aliased=="x2"` | regex match + numeric equality |
| 7 | `mm_fixed_effect_vcov_from_payload()` contract drift detection: rejects non-finite matrix, asymmetric matrix, mismatched coef names, unavailable-with-matrix, unavailable-without-reason | `expect_error(..., class="mm_schema_error")` |
| 8 | `saveRDS` + `readRDS` + `revive()` round-trip preserves `fixef`, `predict`, `audit`, `changes`, `random_blocks`; fresh process revival passes `stopifnot` checks | `ignore_attr=TRUE`; cross-process via `system2` |

### Skips

- Test 8 (cross-session): `testthat::skip_on_cran()`.
- No `skip_if_not_installed` for lme4 (uses synthetic data; lme4 not needed).

### What is NOT tested (gaps)

- **GLMM revival**: `revive.mm_fit` is dispatched on class `mm_fit`; there is no test that a `mm_glmm` fit (which also inherits `mm_fit`) round-trips correctly through `saveRDS`/`revive`.
- **`revive.default` error path**: Passing a non-`mm_fit` object to `revive()` should raise `mm_arg_error`. Not tested.
- **Revival with missing `artifact` field**: The `revive.mm_fit` code checks `is.null(fit$artifact)` and raises `mm_arg_error`. That guard path is untested.
- **`fit_handle_alive()` semantics post-revival**: After revival the handle is always dead; the file tests `fit_handle_alive(fit) == FALSE` once on a fresh fit (Phase 2 never holds a live Rust handle), but `fit_handle_alive` on a live handle (Phase 1 after `lmm()` in an active session) is not exercised—though that is technically tested in `test-lmm.R`.
- **`getME` requesting an unrecognised name**: Should raise an error with useful messaging; not covered.
- **`model.matrix` with unknown `type` argument**: The extractor has `match.arg`; the error path is untested.
- **Inference table method enumeration completeness**: The test checks `method %in% c("satterthwaite", "asymptotic_wald_z", "not_computed")`, but if a new Rust method string is added, this will silently pass. No exhaustiveness assertion.
- **`reproducibility` threshold count floor**: `nrow(repro$thresholds) >= 1L` — only checks a floor of 1. Actual threshold rows (fixef tolerance, theta tolerance, logLik tolerance, sigma tolerance) are specified in the PRD but not individually verified.
- **`vcov` PSD guard**: After constructing the vcov from the payload, symmetry is checked (`isSymmetric`) but positive semi-definiteness (all eigenvalues ≥ 0) is not verified.
- **vcov dimnames**: `dimnames(V)` should equal `list(names(fixef(fit)), names(fixef(fit)))`; this is not asserted anywhere in the file.

---

## test-schema-versioning.R

### What is tested

Seven `test_that` blocks.

| # | Test | Coverage |
|---|------|---------|
| 1 | `mm_json_negotiate` accepts the `formula v0` schema | positive path |
| 2 | `mm_json_negotiate` accepts extra header fields (`crate_version`, `package_version`) beyond the required two | extra-fields tolerance |
| 3 | Version mismatch (`formula / v99`) raises `mm_schema_error` with message matching `"version mismatch"` | `expect_error` + `expect_match` |
| 4 | Unknown schema name raises `mm_schema_error` matching `"unknown schema"` | `expect_error` + `expect_match` |
| 5 | Malformed inputs (NULL, string, missing `schema_version`, non-scalar `schema_name`, NA `schema_name`) all raise `mm_schema_error` | five `expect_error` calls |
| 6 | The `mm_schema_error` condition carries the offending header under `$input` (not `$header`) | `expect_identical(caught$input, bad)` |
| 7 | `mm_json_known_schemas()` returns a data frame with `name`/`version` columns; checks exact version strings for `formula`, `mixedmodels.model_audit_report`, `mixedmodels.random_term_card`, `mixedmodels.fixed_effect_inference_table`, `mixedmodels.marginal_quantity_table` | `expect_identical` per schema |
| 8 | `marginal_quantity_table` schema JSON file exists, pins `schema_name` and `schema_version` constants, and contains all required row-level fields | `expect_true(all(...%in% row_required))` |
| 9 | `mm_json_known_schemas()` content is a superset of the manifest's `schema_versions`; versions agree | loop over manifest schemas |

### Skips

- None (`skip_on_cran` / `skip_if_not_installed` absent).

### What is NOT tested (gaps)

- **Round-trip negotiate for every Rust-facing schema**: `mixeff.lmm_predict_new`, `mixeff.lmm_cond_var`, `mixedmodels.compiled_model_artifact`, `mixedmodels.fixed_effect_covariance_matrix` are all negotiated at runtime in the R source but never directly exercised in `test-schema-versioning.R`.
- **`schema_name` as integer/logical (non-character scalar)**: The malformed-input tests cover `NULL`, string, missing field, length-2 character, and NA. A length-1 integer or logical value might slip through validation depending on the Rust side implementation.
- **`schema_version` with empty-string value**: `schema_version = ""` is not in the malformed battery; it is a plausible edge case from downstream JSON parsers.
- **Schema document existence for all registered schemas**: Only `mixedmodels.marginal_quantity_table` has its schema JSON verified via `system.file(...)`. The other registered schemas (`mixedmodels.fixed_effect_inference_table`, `mixedmodels.random_term_card`, etc.) have no corresponding file-existence test.
- **`mm_json_known_schemas` column types**: The test checks column names and specific values but not that `name` and `version` are character vectors (not factors, not lists).
- **Manifest `schema_versions` superset direction**: The existing test checks that every manifest schema appears in the negotiator's known set. The reverse — that every negotiator entry also appears in the manifest — is NOT checked. A schema could be added to the Rust side without being surfaced in the manifest.
- **Future-version tolerance**: There is no test that `mm_json_negotiate` on a schema with a newer version (e.g., `mixedmodels.model_audit_report / 3`) raises an error with the correct message rather than silently accepting. This is important for catching silent forward-compatibility drift.
- **Concurrent / re-entrant negotiate calls**: Not applicable as a unit test, but the registry is a closed Rust table; no smoke test verifies two distinct schemas can be negotiated in the same R session without interference.

---

## Summary table

| File | Tests | Skips | Gap severity |
|------|-------|-------|-------------|
| test-predict-newdata.R | 6 | `skip_if_not_installed("lme4")` on all | Several major gaps: se.fit contract, multi-group, GLMM stub, column-count mismatch |
| test-ranef-condvar.R | 4 | `skip_if_not_installed("lme4")` on all | Major gaps: intercept-only, multi-group, GLMM stub, Rust failure fallback |
| test-phase2-revive.R | 8 | `skip_on_cran` for cross-process test only | Several minor/major gaps: GLMM revival, revive.default error, getME bad name, vcov PSD, dimnames |
| test-schema-versioning.R | 9 | none | Several minor gaps: per-schema negotiate round-trip, JSON file existence coverage, reverse superset check |
