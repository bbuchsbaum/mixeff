# mixeff Source Survey — R/inference.R

**File:** `R/inference.R`  
**Size:** 2082 lines  
**Survey date:** 2026-05-31  
**Surveyor:** Claude Sonnet 4.6 (subagent)

---

## Overview

`inference.R` is the entire fixed-effect and random-effect inference surface of mixeff. It contains every exported function for hypothesis testing, confidence intervals, estimability assessment, and the convenience linear-combination helper. All work is stateless: the R layer reconstructs a Rust bridge payload from `fit$model_frame` + `fit$formula` + `fit$control` and re-runs the Rust engine for each inference call (no persistent handle is kept after `lmm()` returns).

---

## Exported Functions

### `contrast()` / `contrast.mm_lmm()`

**Signature:**
```r
contrast(fit, L, rhs = 0,
         method = c("auto","satterthwaite","kenward_roger","bootstrap",
                    "asymptotic","boundary_lrt","none"),
         bootstrap = NULL, ...)
```

**Contract:** Fixed-effect contrast front door. Accepts a contrast vector or matrix `L` (one column per fixed effect), validates shape in R, then routes to Rust. Returns `mm_contrast` with a `$table` data-frame.

**Return columns:** `contrast`, `estimate` (= L β̂ − rhs), `rhs`, `std_error`, `df`, `statistic`, `statistic_name`, `p_value`, `method`, `requested_method`, `status`, `reliability`, `reliability_reason`, `estimability`, `reason`, `reason_code`, `reason_detail`, `details`, `notes`.

**Method routing:**
- `"boundary_lrt"` → immediately returns a structured refusal table (not applicable to fixed effects); no Rust call.
- `"none"` → returns `estimate` only; all inference columns are `NA`; no Rust call.
- singular fit + `"satterthwaite"` or `"kenward_roger"` → `mm_boundary_df_unavailable_contrast_table()` refusal; no Rust call.
- All other methods → `mm_rust_contrast_table()` which calls `mm_fixed_effect_contrast_json()` (standard) or `mm_fixed_effect_bootstrap_contrast_json()` (bootstrap).

**Rust FFI calls:**
- `mm_fixed_effect_contrast_json` — standard / asymptotic / Satterthwaite / KR contrast
- `mm_fixed_effect_bootstrap_contrast_json` — parametric bootstrap contrast

**Refusals / NAs:**
- `boundary_lrt` on fixed effects: all inference columns `NA`, `status = "unsupported"`, `reason_code = "boundary_lrt_not_applicable_to_fixed_effects"`.
- `none`: estimates computed, all inference `NA`.
- Singular fit + df method: estimates computed, inference `NA`, `reason_code = "satterthwaite_unavailable_at_boundary"` or `"kenward_roger_unavailable_at_boundary"`.

---

### `bootstrap_control()`

**Signature:**
```r
bootstrap_control(nsim = 999L, seed = NULL,
                  failed_refit_policy = c("exclude","count_extreme","abort"))
```

**Contract:** Constructs a `mm_bootstrap_control` list used by `contrast(..., method="bootstrap")` and `confint(..., method="bootstrap")`. Validates `nsim` (positive integer) and `seed` (NULL or non-negative integer).

**Return:** Named list with class `"mm_bootstrap_control"`.

**Rust FFI:** None — pure R configuration object serialized to JSON before being passed to Rust.

---

### `print.mm_contrast()`

**Contract:** Prints the `$table` data-frame. No hidden columns omitted.

---

### `test_effect()` / `test_effect.mm_lmm()`

**Signature:**
```r
test_effect(fit, term,
            method = c("auto","satterthwaite","kenward_roger","bootstrap",
                       "bootstrap_lrt","cluster_bootstrap","asymptotic",
                       "boundary_lrt","none"),
            bootstrap = NULL, group = NULL, ...)
```

**Contract:** Tests one or more named fixed-effect terms. Validates `term` against `mm_fixed_effect_terms(fit)`. Returns `mm_effect_test` with a `$table` data-frame.

**Return columns:** `term`, `num_df`, `den_df`, `statistic`, `statistic_name`, `p_value`, `method`, `requested_method`, `status`, `reliability`, `reliability_reason`, `reason`, `reason_code`, `reason_detail`, `details`, `notes`.

**Method routing:**
- `"boundary_lrt"` → structured refusal table (not applicable to fixed effects).
- `"none"` → `mm_unavailable_effect_table()`, all inference `NA`.
- `"bootstrap"` → per-term calls to `mm_rust_term_bootstrap_row()` → Rust FFI `mm_fixed_effect_bootstrap_term_json`.
- `"bootstrap_lrt"` → per-term calls to `mm_rust_term_bootstrap_lrt_row()` → Rust FFI `mm_bootstrap_lrt_json`. Refuses REML fits with `reason_code = "bootstrap_lrt_requires_ml"`.
- `"cluster_bootstrap"` → `mm_cluster_bootstrap_unavailable_effect_table()` — structured refusal. In schema 1.0.0 cluster resampling is an estimator-distribution target only; p-values are never certified. Requires `group` for multi-factor models; if `group` is NULL and model is multi-factor, `reason_code = "cluster_bootstrap_multifactor_ambiguous"`.
- singular fit + `"satterthwaite"` or `"kenward_roger"` → `mm_boundary_df_unavailable_effect_table()`.
- All other methods → `mm_rust_term_table()` → Rust FFI `mm_fixed_effect_term_json`.

**Rust FFI calls:**
- `mm_fixed_effect_term_json` — term-level Wald / Satterthwaite / KR F-tests
- `mm_fixed_effect_bootstrap_term_json` — parametric bootstrap term test
- `mm_bootstrap_lrt_json` — bootstrap likelihood-ratio test

**Refusals / NAs:** All refusal paths produce fully structured rows with stable `reason_code` values, never silent `NA` without explanation.

---

### `test_random_effect()` / `test_random_effect.mm_lmm()`

**Signature:**
```r
test_random_effect(fit, term,
                   method = c("boundary_lrt"),
                   refit_for_comparison = c("auto","error","ml"), ...)
```

**Contract:** Tests a single random-effect variance component using a boundary-aware nested-ML likelihood-ratio test (Self-Liang 50:50 mixture reference distribution). The `term` argument accepts term id (`"r0"`), original formula fragment (`"(1 | subject)"`), or grouping factor name (`"subject"`). Returns `mm_random_effect_test` with a one-row `$table`.

**Return columns:** `term`, `term_id`, `group`, `theta_parameters`, `statistic`, `statistic_name` (`"chi_bar_square"`), `ordinary_chisq_dof`, `p_value`, `method` (`"boundary_lrt_self_liang_mixture"` when certified), `requested_method`, `status`, `reason`, `reason_code`, `reference_distribution`, `refit`, `details`, `notes`.

**Implementation steps:**
1. Refits to ML if `fit$REML` is TRUE (unless `refit_for_comparison = "error"`).
2. Builds `mm_random_effect_term_table()` from `fit$artifact$semantic_model$random_terms`.
3. Constructs reduced formula by dropping the selected term.
4. If any terms remain after drop, refits reduced model via `lmm(..., REML=FALSE)`.
5. Calls `mm_boundary_lrt_json(reduced_payload, full_payload, reduced_formula_string)`.
6. Parses schema-versioned JSON (`schema_name = "mixedmodels.boundary_lrt"`, `schema_version = "1.0.0"`).

**Rust FFI calls:**
- `mm_boundary_lrt_json` — boundary LRT engine (stage D.3, bead `bd-01KRFGFSK4A0MGPFQVCNY5SYFK`)

**Refusals:**
- REML fit + `refit_for_comparison = "error"` → `mm_inference_unavailable`, `reason_code = "boundary_lrt_requires_ml"`.
- No random-effect terms → `mm_inference_unavailable`, `reason_code = "boundary_lrt_requires_variance_component_comparison"`.
- Ambiguous / unknown term → `mm_arg_error`.
- Only method currently accepted is `"boundary_lrt"`.

---

### `print.mm_random_effect_test()`

**Contract:** Prints selected columns (`term`, `group`, `statistic`, `statistic_name`, `p_value`, `reference_distribution`, `status`, `reason_code`). Full `$table` always accessible on the object.

---

### `estimability()` / `estimability.mm_lmm()`

**Signature:**
```r
estimability(fit, L = NULL, ...)
```

**Contract:** Assesses estimability of each row of `L` (defaults to the full fixed-effect coefficient basis, i.e., identity matrix). Routes through `mm_rust_contrast_table()` using `method = "auto"` and extracts the estimability sub-payload from each inference row.

**Return:** `mm_estimability` with `$table` columns: `contrast`, `estimable` (logical), `status`, `rank`, `requested_rank`, `reason`.

**Rust FFI:** `mm_fixed_effect_contrast_json` (via `mm_rust_contrast_table`).

**Refusals:** If the Rust call fails, returns `estimable = NA`, `status = "not_assessed"`, `reason = conditionMessage(error)`.

---

### `print.mm_estimability()`

**Contract:** Prints the `$table` data-frame.

---

### `df_for_contrast()` / `df_for_contrast.mm_lmm()`

**Signature:**
```r
df_for_contrast(fit, L,
                method = c("auto","satterthwaite","kenward_roger",
                           "bootstrap","asymptotic","none"), ...)
```

**Contract:** Returns a named numeric vector of denominator degrees of freedom for each contrast row. On success, carries `method` and `requested_method` attributes. On failure, returns all-`NA` with `mm_unavailable_reason` attribute.

**Return:** Class `c("mm_df_for_contrast","numeric")` with attributes `method`, `requested_method`, and (on failure) `mm_unavailable_reason`.

**Rust FFI:** `mm_fixed_effect_contrast_json` (via `mm_rust_contrast_table`).

**Refusals:** `method = "none"` returns `NA` immediately. Rust errors return all-`NA` with reason.

---

### `print.mm_df_for_contrast()`

**Contract:** Prints numeric values plus `method` and `reason` attributes.

---

### `confint.mm_lmm()`

**Signature:**
```r
confint(object, parm, level = 0.95,
        method = c("wald","asymptotic","bootstrap","profile"),
        bootstrap = NULL,
        interval = c("percentile","basic"), ...)
```

**Contract:** Computes confidence intervals for fixed-effect parameters. `"asymptotic"` is accepted as a synonym for `"wald"`.

**Method routing:**
- `"wald"` / `"asymptotic"` → pure R: `β̂ ± z_{α/2} × SE` from stored `fit$std_errors`. Does **not** call Rust. Carries `status = "not_certified_by_rust_inference_contract"`.
- `"bootstrap"` → `mm_bootstrap_confint()` → `mm_full_model_bootstrap_contrast_json` (Rust). Computes per-parameter full-model bootstrap distribution; supports `"percentile"` and `"basic"` interval types.
- `"profile"` → `mm_profile_confint()` → `.Call(wrap__mm_lmm_profile_confint_json, ...)` (Rust). Schema-versioned JSON; `level` passed through. Under REML, beta parameters are explicitly refused with `reason_code = "profile_beta_unavailable_under_reml"` — these appear as typed-refusal rows rather than silent NAs.

**Return:** `mm_confint` (class `c("mm_confint","matrix")`). Attributes: `method`, `status`, and (bootstrap) `interval`, `bootstrap` (list of per-parameter payloads); (profile) `fit_criterion`, `mm_profile` (list with schema, table, notes).

**Rust FFI calls:**
- `mm_full_model_bootstrap_contrast_json` — full-model bootstrap CI
- `wrap__mm_lmm_profile_confint_json` — profile likelihood CI (called via `.Call` directly, not the extendr wrapper alias)

**Refusals / NAs:**
- Wald: certified by stored SE only; status explicitly marked `"not_certified_by_rust_inference_contract"`.
- Profile under REML: beta rows present with `NA` bounds and `reason_code = "profile_beta_unavailable_under_reml"`.

---

### `print.mm_confint()`

**Contract:** Prints the numeric matrix, method, interval type (if present), status. For bootstrap objects also prints a per-parameter summary table (requested/successful/failed replicates, boundary rate, seed) and any notes from the Rust payload.

---

### `mm_lincomb()` / `mm_lincomb.mm_lmm()` / `mm_lincomb.mm_glmm()`

**Signature:**
```r
mm_lincomb(fit, weights, level = 0.95, method = NULL, ...)
# mm_lmm: method = c("auto","satterthwaite","kenward_roger","asymptotic")
# mm_glmm: method = "asymptotic" only
```

**Contract:** Convenience Wald linear combination helper. `weights` is a named numeric vector (or named list / single-row data.frame) mapping coefficient names to weights; unlisted coefficients contribute zero.

**For `mm_lmm`:**
- Computes estimate, SE, and (unless `"asymptotic"`) delegates df computation to `df_for_contrast()` → Rust.
- Returns t-statistic with Satterthwaite df when df is finite and positive; falls back to `NA` p-value when df unavailable.

**For `mm_glmm`:**
- Only `"asymptotic"` accepted. Uses `vcov(fit)` and z-statistic. No df.

**Return:** Single-row data.frame with columns `estimate`, `std_error`, `statistic`, `statistic_name` (`"t"` or `"z"`), `df`, `p_value`, `lower`, `upper`, `method`. Carries `"mm_status"` attribute reflecting vcov reliability.

**Rust FFI:** `mm_fixed_effect_contrast_json` (indirectly via `df_for_contrast` when method ≠ `"asymptotic"`).

---

## Internal Helper Functions (non-exported, load-bearing)

| Function | Purpose |
|---|---|
| `mm_rust_fit_bridge_payload(fit)` | Assembles the stateless bridge payload: `spec_data` (translated data), `formula_string`, `weights`, `control_json`. Called before every Rust inference dispatch. |
| `mm_rust_contrast_table(fit, L, rhs, method, bootstrap)` | Core Rust dispatch for contrast inference; parses JSON via `mm_json_parse_fixed_effect_inference_table`. |
| `mm_rust_term_table(fit, method)` | Rust dispatch for term-level F-test inference. |
| `mm_rust_term_bootstrap_row(fit, term, bootstrap)` | Rust dispatch for a single parametric bootstrap term test. |
| `mm_rust_term_bootstrap_lrt_row(fit, term, bootstrap)` | Rust dispatch for bootstrap LRT; refuses REML with stable reason. |
| `mm_boundary_lrt_bridge_payload(fit)` | Adds `$REML` flag to bridge payload for boundary LRT. |
| `mm_boundary_lrt_ml_fit(fit, refit_for_comparison)` | Refits REML → ML for boundary LRT; refuses if policy = `"error"`. |
| `mm_random_effect_term_table(fit)` | Builds data-frame of random-effect terms from `fit$artifact$semantic_model$random_terms`. Computes `theta_parameters` from covariance structure. |
| `mm_match_random_effect_term(terms, term)` | Matches term by id, fragment text, or group name; errors on ambiguity. |
| `mm_drop_random_term_formula(fit, drop_index)` | Reconstructs formula with specified random term dropped. |
| `mm_json_parse_boundary_lrt(json)` | Validates `schema_name = "mixedmodels.boundary_lrt"`, `schema_version = "1.0.0"` and parses JSON. |
| `mm_boundary_lrt_table(...)` | Builds the one-row result data-frame from parsed boundary-LRT payload. |
| `mm_boundary_lrt_reference_label(mixture)` | Renders the mixture distribution label (e.g., `"0.5 * chi-square(0) + 0.5 * chi-square(1)"`). |
| `mm_profile_confint(fit, parm, level)` | Profile CI wrapper: calls Rust, translates upstream parameter names (β₁→beta name, σ→sigma, θ₁→theta1), appends REML-refusal rows for betas. |
| `mm_profile_confint_payload(fit, level)` | Raw `.Call(wrap__mm_lmm_profile_confint_json, ...)` dispatch + schema negotiation. |
| `mm_map_profile_parameter(upstream, fit)` | Maps Unicode upstream names (β, σ, θ with numeric suffix) to R-side names. Unknown names surface as `kind = "unknown"` rather than being silently dropped. |
| `mm_contrast_matrix(L, fit)` | Validates and normalises L: vector → 1-row matrix; checks column count; assigns rownames. |
| `mm_term_to_l_matrix(fit, term)` | Builds an L matrix for a named term by matching coefficient names with regex. |
| `mm_boundary_df_method_unavailable(fit, method)` | Returns TRUE when method ∈ {satterthwaite, kenward_roger} AND `is_singular(fit)` is TRUE. |
| `mm_bootstrap_confint(fit, parm, level, bootstrap, interval)` | Bootstrap CI: one `mm_full_model_bootstrap_contrast_json` call per parameter. |
| `mm_bootstrap_reliability(certified, successful, mcse, min_moderate)` | Grades reliability: `"moderate"` (≥999 replicates, finite MCSE), `"low"` (fewer/non-finite), `"not_available"` (uncertified). |
| `mm_inference_row_unavailable(term, method, reason, reason_code)` | Canonical one-row refusal frame used by bootstrap_lrt path. |
| `mm_cluster_bootstrap_unavailable_effect_table(fit, term, group)` | Cluster-bootstrap refusal table; always `status = "not_assessed"` in schema 1.0.0. |
| `mm_fixed_effect_terms(fit)` | Returns fixed-effect term labels from `fit$artifact$semantic_model$fixed_terms`, falls back to `names(fit$beta)`. |
| `mm_lincomb_weights_vector(weights, fixef_names)` | Validates and expands sparse weight vector to full coefficient length. |
| `mm_lincomb_fixef(fit)` | Extracts named fixed effects via `fixef()` or `fit$beta`. |
| `mm_lincomb_status_from_vcov(V)` | Extracts reliability metadata from `vcov` attributes. |
| `regex_escape(s)` | Escapes regex metacharacters for term name pattern matching. |

---

## Rust FFI Call Inventory

All calls are stateless (no persistent handle). The bridge payload is rebuilt from `fit$model_frame`, `fit$formula`, `fit$control`, and `fit$weights` on every call.

| Rust FFI function | Called from | Purpose |
|---|---|---|
| `mm_fixed_effect_contrast_json` | `mm_rust_contrast_table` | Contrast inference (Wald/asymptotic/Satterthwaite/KR) |
| `mm_fixed_effect_bootstrap_contrast_json` | `mm_rust_contrast_table` | Parametric bootstrap contrast |
| `mm_fixed_effect_term_json` | `mm_rust_term_table` | Term-level F-test |
| `mm_fixed_effect_bootstrap_term_json` | `mm_rust_term_bootstrap_row` | Parametric bootstrap term test |
| `mm_bootstrap_lrt_json` | `mm_rust_term_bootstrap_lrt_row` | Bootstrap likelihood-ratio test |
| `mm_full_model_bootstrap_contrast_json` | `mm_bootstrap_confint`, `mm_full_model_bootstrap_payload` | Full-model bootstrap CI |
| `mm_boundary_lrt_json` | `test_random_effect.mm_lmm` | Boundary-aware variance component LRT |
| `wrap__mm_lmm_profile_confint_json` | `mm_profile_confint_payload` | Profile likelihood CI |

---

## Structured Refusal Policy

The file consistently returns NA inference columns with diagnostic metadata rather than erroring silently. Every refusal path provides:
- `status` (closed enum: `"not_assessed"`, `"unsupported"`, `"available"`)
- `reliability` (closed enum: `"not_available"`, `"low"`, `"moderate"`, `"high"`)
- `reason` (human-readable string)
- `reason_code` (stable machine-readable label)

Notable reason codes:
- `boundary_lrt_not_applicable_to_fixed_effects`
- `satterthwaite_unavailable_at_boundary` / `kenward_roger_unavailable_at_boundary`
- `boundary_lrt_requires_ml`
- `boundary_lrt_requires_variance_component_comparison`
- `bootstrap_lrt_requires_ml`
- `bootstrap_lrt_reduced_formula_failed`
- `bootstrap_lrt_engine_refused`
- `bootstrap_cluster_resample_p_value_unavailable`
- `cluster_bootstrap_multifactor_ambiguous`
- `profile_beta_unavailable_under_reml`
- `inference_not_requested` (method = "none")

---

## Notable Design Notes

1. **No persistent Rust handle.** Every inference call rebuilds the bridge payload. This is correct by design (JSON artifacts are the source of truth) but means repeated inference calls re-run the optimizer.

2. **Wald CI is not Rust-certified.** `confint(..., method="wald")` uses `fit$std_errors` directly and carries `status = "not_certified_by_rust_inference_contract"`. Profile and bootstrap CIs are Rust-certified.

3. **`cluster_bootstrap` in test_effect is always a refusal** in schema 1.0.0. It is accepted as a valid `method` argument, produces no error, but always returns `status = "not_assessed"` and `p_value = NA` with a stable reason. The group argument is validated but the distribution is estimator-distribution only.

4. **Profile CI under REML explicitly refuses betas.** Rather than silently omitting rows, `mm_profile_confint` appends typed-refusal rows for each requested beta with `reason_code = "profile_beta_unavailable_under_reml"`.

5. **Stage D.3 reference.** The profile CI path is tagged with bead `bd-01KRFGFSK4A0MGPFQVCNY5SYFK`.

6. **No TODOs, FIXMEs, stubs, or deferred markers** anywhere in the file. All code paths are implemented.

---

## Gaps / Issues Found

### G1: `cluster_bootstrap` in `test_effect` always returns refusal
`test_effect(..., method="cluster_bootstrap")` is accepted and dispatched but always returns `status = "not_assessed"` and `p_value = NA` with `reason_code = "bootstrap_cluster_resample_p_value_unavailable"`. This is by design per schema 1.0.0 (cluster resampling is an estimator-distribution target only, not a p-value certification path). **Classification: out-of-scope-by-design.**

### G2: Wald CI not Rust-certified
`confint(..., method="wald")` computes intervals from stored `fit$std_errors` without a Rust call. The `status` attribute is explicitly set to `"not_certified_by_rust_inference_contract"`. **Classification: partial** — the feature works but lacks the audit-trail certification that other methods provide.

### G3: Profile CI betas unavailable under REML
Under REML fits, `confint(..., method="profile")` cannot profile beta parameters (upstream contract). The wrapper surfaces this as typed-refusal rows rather than silently returning `NA`. Users must refit with `REML=FALSE` to get profile CIs for fixed effects. **Classification: upstream-blocked** (Rust upstream contract).

### G4: `test_effect.mm_lmm` only defined for `mm_lmm`
The generic `test_effect` has no `mm_glmm` method. There is no S3 registration for `test_effect.mm_glmm` in NAMESPACE. **Classification: in-scope-missing** (GLMM inference beyond `summary()` Wald-z table is not yet implemented).

### G5: `contrast()` only registered for `mm_lmm`
`contrast.mm_glmm` does not exist. GLMM users must use `mm_lincomb()` (asymptotic Wald z only) or `summary(..., tests="coefficients")`. **Classification: in-scope-missing** for a complete GLMM inference surface; acceptable for current v1 scope.

### G6: `confint.mm_lmm` only; no `confint.mm_glmm`
`confint` is dispatched only for `mm_lmm`. **Classification: in-scope-missing** for GLMM (PRD §3 defers profile-LL CIs for GLMM to v2; asymptotic Wald would be reasonable to add).

### G7: `estimability.mm_lmm` only; no `mm_glmm` method
**Classification: in-scope-missing.**

### G8: No `test_random_effect` for crossed / multi-term models beyond one-at-a-time
`test_random_effect` tests exactly one variance component at a time. Joint tests of multiple random-effect terms are not supported. **Classification: out-of-scope-by-design** for v1 (one parameter per boundary-LRT call is a shape restriction).

### G9: `bootstrap_lrt` refuses REML without an auto-refit option
Unlike `test_random_effect` which has `refit_for_comparison = "auto"`, `bootstrap_lrt` in `test_effect` unconditionally refuses REML with `reason_code = "bootstrap_lrt_requires_ml"` and no auto-refit. The user must manually refit. **Classification: partial** — works correctly but lacks the convenience `refit_for_comparison` ergonomic that `test_random_effect` provides.
