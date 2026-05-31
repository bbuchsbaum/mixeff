# mixeff-src-11: FFI Surface Inventory

**Files surveyed:** `R/extendr-wrappers.R`, `src/rust/src/lib.rs`
**Supporting files consulted:** `R/fit-lmm.R`, `R/glmm.R`, `R/inference.R`, `R/methods-extract.R`, `R/predict.R`, `R/compare.R`, `R/reporting.R`, `NAMESPACE`

---

## 1. Overview

The FFI layer is a pure JSON-over-strings bridge. Every exported Rust function
returns either a JSON string or an error string with a typed prefix
(`mm_formula_error:`, `mm_fit_error:`, `mm_schema_error:`, `mm_data_error:`,
`mm_inference_unavailable:`, `mm_arg_error:`). R never receives a raw Rust
handle. The R wrapper owns the S3 surface; Rust owns all numerical computation
and audit wording.

The `extendr_module!` block registers exactly **20 functions**. All are
`@noRd` (internal primitives); the public API surface lives entirely in R.

---

## 2. Registered FFI Functions

### 2.1 `mm_parse_formula(formula: &str) -> String`

- **R wrapper:** `mm_parse_formula(formula)` — also exported to users as a
  diagnostic.
- **What it does:** Calls `mixeff_rs::formula::parse_formula`, returns the
  canonical `Display` rendering of the parsed formula.
- **On failure:** Returns `mm_formula_error: <message>` as an R error.
- **Notes:** The only direct user-visible FFI primitive (exported in NAMESPACE).

### 2.2 `mm_formula_manifest() -> List`

- **R wrapper:** `mm_formula_manifest()` — exported.
- **What it returns:** Named R list with:
  - `mixeff_rust_version` (CARGO_PKG_VERSION)
  - `crate_version` ("0.1.0")
  - `schema_versions` (named list of all schema name → version pairs)
  - `formula_features` (operators, intercept_forms, random_term_forms,
    transformations)
  - `capabilities` (named logical flags: `parse_formula`, `compile_model`,
    `audit_design`, `explain_model`, `random_options`, `compare_covariance`,
    `fit_lmm`, `audit`, `changes`, `diagnostics`, `fit_status`,
    `parameterization`, `roles`, `as_json`, `fit_glmm`, `simulate`,
    `inference`, `fixed_effect_inference_table`, `satterthwaite`,
    `kenward_roger_explicit`, `bootstrap_fixed_effect_payload`,
    `model_comparison_table`, `fit_summary_payload`, `marginal_quantity_table`,
    `marginal_quantities`)

### 2.3 `mm_json_negotiate_one(name: &str, version: &str) -> bool`

- **R wrapper:** `mm_json_negotiate_one(name, version)` — used internally by
  `mm_json_negotiate()` (exported).
- **What it does:** Checks `(name, version)` against the closed `KNOWN_SCHEMAS`
  table. Returns `TRUE` on exact match; throws `mm_schema_error:` on version
  mismatch or unknown schema.

### 2.4 `mm_json_known_schemas() -> List`

- **R wrapper:** `mm_json_known_schemas()` — exported.
- **What it returns:** List with `$name` and `$version` character vectors
  covering all 9 known schemas.

**Known schemas (as of this survey):**
| Schema name | Version |
|---|---|
| `formula` | `v0` |
| `mixedmodels.compiled_model_artifact` | `1` |
| `mixedmodels.model_audit_report` | `2` |
| `mixedmodels.random_term_card` | `1` |
| `mixedmodels.fixed_effect_inference_table` | (upstream const) |
| `mixedmodels.marginal_quantity_table` | `1.0.0` |
| `mixedmodels.model_comparison_table` | `1.0.0` |
| `mixedmodels.boundary_lrt` | (upstream const) |
| `mixedmodels.fit_summary` | (upstream const) |
| `mixedmodels.profile_likelihood_ci` | (upstream const) |

### 2.5 `mm_compile_model_json(formula, column_order, numeric_columns, categorical_values, categorical_levels) -> String`

- **R wrapper:** Called by `compile_model()`.
- **Pipeline:** `parse_formula` → `compile_formula_ir` → `CompiledModelArtifact::new`
  → `attach_design_audit`.
- **Returns:** Serialized `CompiledModelArtifact` JSON (schema
  `mixedmodels.compiled_model_artifact` v1). No fit.
- **Errors:** `mm_formula_error:`, `mm_data_error:`, `mm_schema_error:`.

### 2.6 `mm_fit_lmm_json(formula, reml, column_order, numeric_columns, categorical_values, categorical_levels, weights, control_json) -> String`

- **R wrapper:** Called directly by `lmm()`.
- **Pipeline:** `parse_formula` → `build_dataframe` → `LinearMixedModel::new` →
  `.fit(reml)`.
- **Returns:** JSON payload with schema `mixeff.lmm_fit_result` v1, containing:
  - `artifact_json` (post-fit `CompiledModelArtifact`)
  - `formula`, `reml`, `beta`, `beta_names`, `theta`, `sigma`
  - `log_likelihood`, `deviance`, `aic`, `bic`
  - `nobs`, `dof`, `df_residual`, `fit_status`
  - `std_errors`, `fixed_fitted`, `fitted`, `residuals`
  - `ranef` (per-group: group, levels, names, values)
  - `varcorr` (components with group/names/std_dev/correlations + residual_sd)
  - `fit_summary` (schema `mixedmodels.fit_summary`)
  - `optimizer` (backend, algorithm, return_value, function_evaluations,
    objective, reml)
- **Weights:** Optional; validated for finiteness and positivity. Empty vector
  means no weights.
- **control_json:** Parsed but fields not yet acted upon in Rust (validated for
  JSON syntax only).
- **Errors:** `mm_fit_error:`, `mm_formula_error:`, `mm_schema_error:`.

### 2.7 `mm_fit_glmm_json(formula, family, link, method, n_agq, column_order, numeric_columns, categorical_values, categorical_levels, control_json) -> String`

- **R wrapper:** Called directly by `glmm()`.
- **Families accepted:** `"bernoulli"` / `"binomial"`, `"poisson"`, `"gamma"`.
- **Links accepted:** `"identity"`, `"log"`, `"logit"`, `"probit"`,
  `"cloglog"`, `"inverse"`, `"sqrt"`.
- **Method:**
  - `"pirls_profiled"` → `fast = true` (available).
  - `"joint_laplace"` → **hard refused**: `"estimation_method_unavailable:
    method='joint_laplace' requires the upstream nlopt backend, which is
    disabled in this vendored build"`. This is the CRAN-oriented build
    without the `nlopt` feature.
- **Returns:** JSON with schema `mixeff.glmm_fit_result` v1 containing same
  fields as the LMM result minus `fixed_fitted` (stored as `NULL` on R side)
  plus `dispersion`, `family`, `link`, `method`, `n_agq`.
- **GLMM prediction:** `predict.mm_glmm` hard-refuses ("GLMM prediction is not
  certified by the current Rust contract").
- **GLMM condVar:** `ranef(condVar=TRUE)` on `mm_glmm` returns unavailable
  arrays with reason `"random_effect_conditional_variance_unavailable_for_glmm"`.

### 2.8 `mm_fixed_effect_contrast_json(..., l_values, nrow, ncol, labels, rhs, method) -> String`

- **R wrapper:** Called by `contrast.mm_lmm` and `df_for_contrast.mm_lmm`
  via `mm_rust_contrast_table()`.
- **What it does:** Refits the LMM from bridge data, then calls
  `LinearMixedModel::fixed_effect_contrast_inference_table(hypotheses, method)`.
- **Methods:** `"auto"`, `"asymptotic"`/`"asymptotic_wald_z"`,
  `"satterthwaite"`, `"kenward_roger"`, `"bootstrap"` (parametric-bootstrap
  null variant).
- **Returns:** Serialized `FixedEffectInferenceTable` JSON.

### 2.9 `mm_fixed_effect_bootstrap_contrast_json(..., bootstrap_options_json) -> String`

- **R wrapper:** Called when `method = "bootstrap"` and `!is.null(bootstrap)`
  in `contrast.mm_lmm`.
- **What it does:** Null-bootstrap contrast inference (parametric bootstrap
  under the null hypothesis that `L*beta = rhs`).
- **Returns:** Serialized `FixedEffectInferenceTable` JSON.

### 2.10 `mm_full_model_bootstrap_contrast_json(..., bootstrap_options_json, levels) -> String`

- **R wrapper:** Called by `mm_full_model_bootstrap_payload()` in
  `confint.mm_lmm(method="bootstrap")`.
- **What it does:** Full-model parametric bootstrap for CI computation. Uses
  `parametricbootstrap()`, computes percentile and basic intervals in-Rust.
- **Restriction:** Only scalar contrasts (single row L matrix) are certified.
  Multi-row requests return `mm_inference_unavailable:
  full-model bootstrap intervals are currently certified only for scalar
  contrasts`.
- **Returns:** JSON with `intervals` (method/level/lower/upper/n),
  `metadata`, `replicate_statistics`, `observed_estimate`, `contrast_label`.

### 2.11 `mm_fixed_effect_bootstrap_term_json(..., label, bootstrap_options_json) -> String`

- **R wrapper:** Called by `mm_rust_term_bootstrap_row()` in
  `test_effect.mm_lmm(method="bootstrap")`.
- **What it does:** Null-bootstrap term inference wrapping
  `LinearMixedModel::fixed_effect_null_bootstrap_inference_row` with
  `kind = Term`. Multi-df terms produce F-form rows.
- **Returns:** Serialized `FixedEffectInferenceTable` JSON (single row).

### 2.12 `mm_fixed_effect_term_json(..., method) -> String`

- **R wrapper:** Called by `mm_rust_term_table()` in `test_effect.mm_lmm`,
  `anova.mm_lmm`.
- **What it does:** Builds Rust-owned term hypotheses and returns the full term
  inference table.
- **Methods:** Same as contrast: auto/asymptotic/satterthwaite/kenward_roger/
  bootstrap.
- **Returns:** Serialized `FixedEffectInferenceTable` JSON.

### 2.13 `mm_bootstrap_lrt_json(reduced_formula, alternative_formula, ..., bootstrap_options_json) -> String`

- **R wrapper:** Called by `parametric_bootstrap()` (and via `compare(method=
  "bootstrap")`).
- **What it does:** Fits both models (both forced to ML — REML is refused by
  the R layer before this call). Simulates from the reduced model, refits both,
  records LRT statistics. Computes a Monte-Carlo p-value.
- **Returns:** JSON with `observed_statistic`, `p_value`, `mcse`, `notes`,
  `payload` (metadata + replicate_statistics).
- **Refusal condition:** `observed_lrt` not finite → `mm_inference_unavailable`.

### 2.14 `mm_compare_models_json(model_payloads, method, refit_policy) -> String`

- **R wrapper:** Called by `mm_compare_table()` in `compare.mm_lmm`.
- **What it does:** Accepts a list of bridge payloads (each with
  `formula_string`, `REML`, `spec_data`, `weights`, `control_json`), fits all
  models, delegates to `ModelComparisonTable::compare_with_options()`.
- **Methods:** `"auto"`, `"lrt"`/`"likelihood_ratio"`, `"bootstrap"` (mapped to
  Auto), `"aic"`/`"information_criteria"`.
- **Refit policies:** `"error"`, `"auto"`/`"ml"`, `"never"`.
- **Returns:** JSON with schema `mixedmodels.model_comparison_table` v1.0.0.
- **Minimum:** Refuses fewer than 2 models.

### 2.15 `mm_boundary_lrt_json(reduced_payload, full_payload, reduced_formula) -> String`

- **R wrapper:** Called by `test_random_effect.mm_lmm`.
- **What it does:** Boundary-aware variance-component LRT. When `reduced_payload`
  is NULL, the reduced model is an ordinary fixed-effect LM (`LinearModelFit::fit`).
  Otherwise fits the reduced LMM from the bridge payload.
  Calls `BoundaryLikelihoodRatioTest::variance_component()`.
- **Returns:** Serialized `BoundaryLikelihoodRatioTest` JSON (schema
  `mixedmodels.boundary_lrt`).

### 2.16 `mm_audit_report_text(artifact_json: &str) -> String`

- **R wrapper:** Called by `print.mm_audit()`.
- **What it does:** Deserializes `CompiledModelArtifact`, returns
  `ModelAuditReport::Display` rendering as a plain-text string.
- **Design note:** All audit English wording is authored in Rust; this enforces
  the "no advice creep" contract (PRD R9).

### 2.17 `mm_audit_report_json(artifact_json: &str) -> String`

- **R wrapper:** Called by `explain_model()`.
- **What it does:** Deserializes `CompiledModelArtifact`, serializes
  `ModelAuditReport` as JSON for structured consumption by R formatters.

### 2.18 `mm_lmm_cond_var_json(...) -> String`

- **R wrapper:** Called by `mm_compute_cond_var_postvars()` in
  `ranef(condVar=TRUE)`.
- **What it does:** Refits the model, calls `model.cond_var()`, serializes
  per-term `p × p × n` PSD arrays column-major (R `array()` compatible).
- **Returns:** JSON with schema `mixeff.lmm_cond_var` v1, fields: `terms`
  (group, names, levels, postvar flat vector, dim `[p,p,n]`).
- **GLMM:** Hard-refused by the R wrapper (`ranef.mm_glmm` returns unavailable
  arrays before calling FFI).

### 2.19 `mm_lmm_predict_new_json(..., allow_new_levels_policy) -> String`

- **R wrapper:** Called by `mm_predict_conditional_newdata()` in
  `predict.mm_lmm(newdata=...)`.
- **What it does:** Refits the model, calls
  `LinearMixedModel::predict_new(&newdf, policy)`.
- **Policies:** `"error"` → `NewReLevels::Error`; `"population"` →
  `NewReLevels::Population`; `"missing"` → `NewReLevels::Missing`.
  (R currently only sends `"error"` or `"population"`; the `"missing"` policy
  is wired but the R predict.mm_lmm does not expose it yet.)
- **Returns:** JSON with schema `mixeff.lmm_predict_new` v1: `predictions`
  (array of numbers or `null`), `policy`, `n_new`. JSON `null` → `NA_real_`.
- **Prediction SE:** Hard-refused by R wrapper with `"prediction_se_unavailable_phase_2"`.
- **Prediction intervals:** Hard-refused before FFI call.

### 2.20 `mm_lmm_profile_confint_json(..., level) -> String`

- **R wrapper:** Called by `mm_profile_confint()` in `confint.mm_lmm(method=
  "profile")`.
- **What it does:** Refits the model, calls `profile_confint_payload(&mut model,
  level)`.
- **Under REML:** The upstream contract omits beta from profiled parameters; R
  translates this absence to a `profile_beta_unavailable_under_reml` reason
  rather than fabricating beta CIs.
- **Returns:** JSON per upstream `PROFILE_LIKELIHOOD_CI_SCHEMA` (contains
  computed intervals + raw profile rows).

### 2.21 `mm_interrupt_demo(iters: i32) -> i32`

- **R wrapper:** `mm_interrupt_demo(iters)` — internal smoke test.
- **What it does:** Loops `iters` times calling `R_CheckUserInterrupt()` then
  returns `iters`. Phase 0 smoke test only; real fit/inference loops do not yet
  call `check_user_interrupt()` from inside Rust.

---

## 3. Internal Rust Helpers (not exported, but load-bearing)

| Helper | Purpose |
|---|---|
| `fit_lmm_from_bridge_data()` | Shared refit primitive used by all inference FFI functions (refits model from scratch each call). |
| `fit_lmm_from_bridge_payload_robj()` | Extracts bridge fields from an R `List` and calls `fit_lmm_from_bridge_data`. Used by `mm_compare_models_json` and `mm_boundary_lrt_json`. |
| `random_effects_json()` / `random_effects_json_glmm()` | Serializes per-term ranef matrices. |
| `varcorr_json()` / `varcorr_json_glmm()` / `varcorr_value()` | Serializes VarCorr structure (components with group/names/std_dev/correlations + residual_sd). |
| `glmm_family()` / `glmm_link()` / `glmm_method()` | Family/link/method string → enum mapping with early refusal for unsupported values. |
| `fixed_effect_hypotheses()` | Constructs `Vec<FixedEffectHypothesis>` from flat L-matrix arguments. |
| `fixed_effect_test_method()` | Maps method string to `FixedEffectTestMethod` enum. |
| `fixed_effect_bootstrap_options()` | Parses bootstrap options JSON (nsim, seed, failed_refit_policy). |
| `optional_case_weights()` | Validates and extracts case weights; empty vector → None. |
| `model_comparison_method()` / `model_comparison_refit_policy()` | Maps comparison method/policy strings. |
| `make_bootstrap_rng()` | Constructs seeded or entropy `StdRng` + `BootstrapSeedRecord`. |
| `quantile_sorted()` | Linear-interpolation quantile on a pre-sorted slice (mirrors upstream helper that is not `pub`). |

---

## 4. R-Side S3 Surface (from NAMESPACE)

### Fit functions (exported)
- `lmm(formula, data, REML, weights, control)` → `mm_lmm` / `mm_fit` / `mm_compiled`
- `glmm(formula, data, family, random, weights, subset, na.action, contrasts, method, nAGQ, inference, control)` → `mm_glmm` / `mm_fit` / `mm_compiled`

### Extractors (lme4-compatible)
`fixef`, `ranef` (condVar supported for LMM; refused for GLMM), `coef`,
`VarCorr`, `sigma`, `logLik`, `deviance`, `AIC`, `BIC`, `nobs`,
`df.residual`, `formula`, `model.frame`, `model.matrix`, `vcov`,
`fitted`, `residuals`

### Inference
- `contrast.mm_lmm` — arbitrary L-matrix fixed-effect tests; methods: auto /
  satterthwaite / kenward_roger / bootstrap / asymptotic / boundary_lrt / none
- `test_effect.mm_lmm` — term-level fixed-effect tests; additionally supports
  bootstrap_lrt, cluster_bootstrap (cluster_bootstrap returns
  `not_assessed` p-values by design in schema 1.0.0)
- `test_random_effect.mm_lmm` — variance-component boundary LRT only
- `estimability.mm_lmm` — contrast estimability assessment
- `df_for_contrast.mm_lmm`
- `confint.mm_lmm` — methods: wald (default), profile, bootstrap; se.fit and
  interval types refused until further Rust certification
- `inference_options.mm_lmm`, `inference_table.mm_lmm`

### Comparison
- `compare.mm_lmm` — LRT / AIC / bootstrap; REML auto-refit
- `parametric_bootstrap(null, alternative, nsim, seed)` — engine-certified
  bootstrap LRT
- `anova.mm_lmm` — dispatches to compare() when multiple models, otherwise
  term table
- `drop1.mm_lmm` — single-term deletion table

### Audit / design
- `compile_model(formula, data)` → `mm_spec` / `mm_compiled`
- `audit_design()`, `audit.mm_fit`
- `explain_model()`
- `diagnostics.mm_compiled`
- `parameterization.mm_compiled`, `roles()`
- `changes.mm_compiled`, `as_json.mm_compiled`
- `optimizer_certificate.mm_compiled`, `fit_status.mm_fit`
- `is_singular.mm_lmm`, `reproducibility.mm_compiled`

### Prediction
- `predict.mm_lmm` — conditional (newdata via Rust) or population (fixed-only
  in R); re.form=NULL → conditional, re.form=NA/~0 → population; partial
  re.form formulas refused; se.fit refused; intervals refused
- `predict.mm_glmm` — hard-refused entirely

### Reporting / marginals
- `model_report.mm_fit`, `reporting_table.mm_fit` / `.mm_model_comparison` /
  `.mm_drop1` / `.mm_random_effect_test` / `.mm_model_report`
- `mm_means.mm_lmm`, `mm_predictions.mm_lmm`, `mm_comparisons.mm_lmm`,
  `mm_lincomb` (LMM + GLMM), `mm_grid.mm_lmm`
- `emm_basis.mm_lmm`, `emm_basis.mm_glmm`, `recover_data.mm_lmm`,
  `recover_data.mm_glmm` (emmeans integration)

### Revive / save
- `revive.mm_fit` — revival from saved JSON artifact without live handle
- `refit.mm_lmm`
- `simulate.mm_lmm`

### Control / options
- `mm_control()`, `bootstrap_control()`, `random_options()`, `compare_covariance()`
- `mm_parse_formula()`, `mm_formula_manifest()`, `mm_json_known_schemas()`,
  `mm_json_negotiate()`

---

## 5. Stubs, Refusals, and Deferred Paths

| Item | Classification | Notes |
|---|---|---|
| `glmm(method="joint_laplace")` | out-of-scope-by-design | Rust hard-refuses: `nlopt` backend not compiled in CRAN build. Error: `estimation_method_unavailable`. |
| `predict.mm_glmm` | out-of-scope-by-design | Hard-refused: `"GLMM prediction is not certified by the current Rust contract."` |
| `ranef(condVar=TRUE)` on `mm_glmm` | out-of-scope-by-design | Returns NA postVar arrays with reason `random_effect_conditional_variance_unavailable_for_glmm`. |
| `se.fit = TRUE` in `predict.mm_lmm` | partial | Returns NA SE vector with attribute `prediction_se_unavailable_phase_2`; no Rust certification yet. |
| `interval != "none"` in `predict.mm_lmm` | out-of-scope-by-design | Hard-refused before FFI call. |
| `re.form = <formula>` in `predict.mm_lmm` | out-of-scope-by-design | Partial re.form formulas (not NULL/NA/~0) raise `mm_inference_unavailable`. |
| `NewReLevels::Missing` policy for predict_new | partial | Wired in Rust FFI, but R currently only sends `"error"` or `"population"`. |
| `confint(method="wald")` Wald asymptotic CIs | partial | Computed in R from stored SE without Rust certification; status = `"not_certified_by_rust_inference_contract"`. |
| `mm_interrupt_demo` | partial | Phase 0 smoke test only. Rust fit/inference loops do not yet call `check_user_interrupt()` internally. |
| `cluster_bootstrap` method in `test_effect` | partial | Accepted as method arg; returns `not_assessed` p-values with reason code `bootstrap_cluster_resample_p_value_unavailable`. Design intent: estimator distribution only. |
| `full-model bootstrap` with multi-row L matrix | in-scope-missing | Rust refuses: "currently certified only for scalar contrasts". Single-row only. |
| `profile confint` beta under REML | partial | Upstream omits beta from profile; R converts absence to `profile_beta_unavailable_under_reml` reason. Not a bug, documented contract. |
| Satterthwaite/Kenward-Roger at boundary | partial | R intercepts and returns `not_assessed` rows before calling Rust when `is_singular(fit)` is TRUE. |
| `control_json` contents used in Rust | test-gap | `control_json` is validated for JSON syntax in Rust but fields (verbose, optimizer settings) are not yet acted upon inside the Rust engine. |
| GLMM families: gaussian / inverse_gaussian | out-of-scope-by-design | Not listed in `mm_glmm_supported_family_links()`; R refuses before FFI call. Rust's `glmm_family()` maps `"bernoulli"/"binomial"`, `"poisson"`, `"gamma"` only. |
| AIC/BIC multi-model comparison via `AIC(m1, m2, ...)` | out-of-scope-by-design | `AIC.mm_lmm` refuses multi-model calls; directs to `compare()`. |

---

## 6. Schema Version Constants (lib.rs)

```
SCHEMA_VERSION_FORMULA                 = "v0"
SCHEMA_VERSION_COMPILED_ARTIFACT       = "1"
SCHEMA_VERSION_MODEL_AUDIT_REPORT      = "2"
SCHEMA_VERSION_RANDOM_TERM_CARD        = "1"
SCHEMA_VERSION_FIXED_EFFECT_INFERENCE_TABLE = FIXED_EFFECT_INFERENCE_TABLE_SCHEMA_VERSION (upstream)
SCHEMA_NAME_MARGINAL_QUANTITY_TABLE    = "mixedmodels.marginal_quantity_table" @ "1.0.0"
SCHEMA_NAME_MODEL_COMPARISON_TABLE     = "mixedmodels.model_comparison_table" @ "1.0.0"
MIXEDMODELS_CRATE_VERSION              = "0.1.0"
```

---

## 7. Wire-format Data Encoding (data.rs)

Data crossing the FFI uses:
- `column_order`: character vector of column names
- `numeric_columns`: named list of numeric vectors
- `categorical_values`: named list of integer-coded values
- `categorical_levels`: named list of character vectors

The R helper `mm_translate_data(model_frame)` constructs these from a
`data.frame`. This encoding is shared across all fit, inference, comparison,
and prediction FFI calls.

---

## 8. Key Architectural Observations

1. **No live handle caching.** Every inference FFI call (contrast, term,
   cond_var, predict_new, profile_confint) calls `fit_lmm_from_bridge_data()`
   which re-parses the formula and re-fits the model. The `rust_handle = NULL`
   field on R objects is a reserved slot; `lazy_cache` caches cond_var results
   in R.

2. **No in-process interrupt hooks in fit loops.** `R_CheckUserInterrupt()` is
   declared and called in `mm_interrupt_demo` but not wired into
   `LinearMixedModel::fit()` or inference loops yet.

3. **All English wording in Rust.** Audit text, random-term card wording, and
   reason strings are authored in Rust and returned verbatim. R formatters
   receive pre-authored strings.

4. **JSON is the source of truth.** Model objects stored in R carry all
   extractable values as plain R vectors/lists. `revive.mm_fit` can
   reconstruct the full object from the artifact JSON without a live handle.

5. **GLMM inference is limited.** Only `summary(mm_glmm, tests="coefficients")`
   (Wald-z table) and `emmeans` integration are available. No fixed-effect term
   tests, no confint, no ranef condVar, no predict for GLMM.
