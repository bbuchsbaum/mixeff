# mixeff source survey — file group 7

**Files:** `R/compare.R`, `R/compare-covariance.R`, `R/changes.R`
**Date:** 2026-05-31
**Family:** mixeff-src-7

---

## 1. `R/compare.R`

### 1.1 Exported functions and S3 methods

#### `compare(object, ...)` — generic
- **Exported:** yes (`export(compare)`)
- **Contract:** Namespace-qualified front door for model comparison. Dispatches on class.
- **Returns:** Depends on method dispatch.

#### `compare.mm_lmm(object, ..., target, method, refit_for_comparison, nsim, seed)`
- **Exported:** yes (`S3method(compare, mm_lmm)`)
- **Arguments:**
  - `object`, `...`: one or more fitted `mm_lmm` objects.
  - `target`: `"fixed_effects"` | `"random_effects"` | `"prediction"` (default `"fixed_effects"`).
  - `method`: `"auto"` | `"lrt"` | `"bootstrap"` | `"aic"` (default `"auto"`).
  - `refit_for_comparison`: `"auto"` | `"error"` | `"ml"` (default `"auto"`).
  - `nsim`: bootstrap replicates (integer; only used when `method = "bootstrap"`).
  - `seed`: optional bootstrap seed.
- **Returns:** `mm_model_comparison` list with fields `table` (data frame), `ledger` (provenance data frame), `fits`, `target`, `method`, `refit_for_comparison`, `bootstrap` (NULL or `mm_parametric_bootstrap`).
- **Rust FFI:** `mm_compare_models_json` (via `mm_compare_table` → `mm_rust_fit_bridge_payload`). Optionally `mm_bootstrap_lrt_json` when `method = "bootstrap"` and `nsim > 0` and exactly two models.
- **REML handling:** REML fits are silently refitted to ML when `refit_for_comparison = "auto"` or `"ml"`; `"error"` aborts.
- **Refuses/NA's:**
  - Bootstrap requested with `nsim = 0` or `>2` models: `p_value` set to `NA_real_`, `method` becomes `"bootstrap_not_run"`.
  - Non-finite observed LRT at Rust level: engine returns error string.
  - Non-nested models: Rust `ModelComparisonTable` sets `lrt_available = false` and records a reason code; R propagates to `status = "not_available"`.
- **No TODO/FIXME markers.**

#### `print.mm_model_comparison(x, ...)`
- **Exported:** yes (`S3method(print, mm_model_comparison)`)
- **Behaviour:** Prints header line `"Model comparison:\n"` then the `table` data frame without row names.

#### `parametric_bootstrap(null, alternative, nsim, seed, ...)`
- **Exported:** yes (`export(parametric_bootstrap)`)
- **Arguments:**
  - `null`, `alternative`: fitted `mm_lmm`. Order is corrected internally (fewer parameters = reduced).
  - `nsim`: positive integer (default 999).
  - `seed`: optional non-negative integer.
- **Returns:** `mm_parametric_bootstrap` list with fields: `observed`, `simulated` (replicate LRT statistics), `p_value`, `nsim`, `successful_replicates`, `completed_replicates`, `boundary_count`, `mcse`, `seed`, `status` (`"available"` or `"not_assessed"`), `reason`, `notes`, `reduced_formula`, `alternative_formula`.
- **Rust FFI:** `mm_bootstrap_lrt_json`. The bridge payload is built from the *alternative* model (so column order, spec data, weights come from the larger model). The null formula string comes from `deparse1(null$formula)`. Both models are refitted inside Rust from bridge data; the bootstrap loop is entirely inside the Rust engine.
- **Refuses:**
  - Either model is REML: throws `mm_inference_unavailable` / `"bootstrap_lrt_requires_ml"`.
  - `alternative$dof <= null$dof` (after auto-swap): throws `"bootstrap_lrt_requires_nested_models"`.
  - `null` formula variables missing from `alternative$model_frame`: throws `"bootstrap_lrt_requires_nested_model_frames"`.
  - Mismatched model-frame values for shared columns: throws `"bootstrap_lrt_requires_same_observations"`.
  - Mismatched weights: throws `"bootstrap_lrt_requires_same_weights"`.
- **p_value certification:** `certified = !is.null(parsed$p_value)`. If the engine does not emit a top-level `p_value` field, `p_value = NA_real_` and `status = "not_assessed"`.
- **No TODO/FIXME markers.**

#### `print.mm_parametric_bootstrap(x, ...)`
- **Exported:** yes (`S3method(print, mm_parametric_bootstrap)`)
- **Behaviour:** Prints replicate accounting (requested, successful/completed, boundary, MCSE, seed). Shows `p.value` only when `status == "available"`; otherwise prints the reason.

#### `anova.mm_lmm(object, ..., type, method, refit_for_comparison)`
- **Exported:** yes (`S3method(anova, mm_lmm)`)
- **Arguments:**
  - `type`: `"III"` | `"II"` | `"I"` (default `"III"`); used only in the single-model branch.
  - `method`: `"auto"` | `"satterthwaite"` | `"kenward_roger"` | `"bootstrap"` | `"asymptotic"` | `"none"`.
  - `refit_for_comparison`: `"auto"` | `"error"` | `"ml"`.
- **Dispatch logic:**
  - If `...` contains additional models → delegates to `compare()` (multi-model branch). Returns `mm_model_comparison`.
  - Otherwise → single-model Type I/II/III F-table. Returns `mm_anova`.
- **Single-model path:** Calls `mm_rust_term_table(object, method)` (Rust FFI via `mm_inference_term_table_json`). Columns renamed: `numerator_df → num_df`, `denominator_df → den_df`. A `type` column is prepended.
- **`method = "none"` path:** Returns `mm_unavailable_effect_table` (all NA inference, no Rust call).
- **Returns:** `mm_anova` object with fields `table`, `type`, `requested_method`, `refit_for_comparison`.
- **No TODO/FIXME markers.**

#### `print.mm_anova(x, ...)`
- **Exported:** yes (`S3method(print, mm_anova)`)
- **Behaviour:** Prints header with type and method, then table without row names.

#### `drop1.mm_lmm(object, scope, test, refit_for_comparison, ...)`
- **Exported:** yes (`S3method(drop1, mm_lmm)`)
- **Arguments:**
  - `scope`: optional character vector of fixed-effect terms to consider dropping.
  - `test`: `"none"` | `"Chisq"` (default `"none"`).
  - `refit_for_comparison`: `"auto"` | `"error"` | `"ml"`.
- **Mechanism:** For each droppable fixed-effect term, rebuilds the reduced formula using `mm_drop_fixed_term_formula()` (preserves all random terms verbatim from `artifact$semantic_model$random_terms`), refits via `lmm()` with `verbose = -1`, computes asymptotic LRT from `deviance()` difference. No Rust bootstrap call.
- **Returns:** `mm_drop1` with fields `table` (columns: `dropped`, `formula`, `df`, `logLik`, `AIC`, `BIC`, `LRT`, `p_value`, `method`), `ledger`, `full`.
- **Refuses/NA's:** `LRT` and `p_value` are `NA_real_` when `test = "none"` or `df <= 0`.
- **Limitation (minor gap):** The asymptotic LRT p-value uses `pchisq(stat, df, lower.tail = FALSE)` — standard chi-square, not boundary-aware. This is appropriate for fixed-effect comparisons but the contrast with `test_random_effect`'s boundary correction is not documented.
- **No TODO/FIXME markers.**

#### `print.mm_drop1(x, ...)`
- **Exported:** yes (`S3method(print, mm_drop1)`)
- **Behaviour:** Prints header and table without row names.

### 1.2 Internal helpers (not exported)

| Helper | Purpose |
|---|---|
| `mm_assert_comparable_lmm(fits)` | Guards: same nobs, same response variable across all fits. |
| `mm_assert_bootstrap_lrt_pair(null, alternative)` | Guards nestedness, model-frame membership, value identity, weight identity. |
| `mm_prepare_comparison_fits(fits, refit_for_comparison)` | Refits REML→ML; tracks which models were refitted. |
| `mm_compare_table(fits, method, refit)` | Calls `mm_compare_models_json`; parses result via `mm_json_parse_model_comparison_table`. |
| `mm_json_parse_model_comparison_table(json)` | Schema validation (`schema_name == "mixedmodels.model_comparison_table"`, `schema_version == "1.0.0"`). |
| `mm_compare_table_from_rust_payload(payload, fits, refit, method)` | Maps Rust rows to R data frame; derives `status` / `method` columns. |
| `mm_comparison_ledger(...)` | Builds provenance data frame with comparison_id, model roles, refit history. |
| `mm_drop1_comparison_ledger(...)` | Provenance data frame specific to drop1 deletes. |
| `mm_comparison_ledger_empty()` | Zero-row template with all required columns. |
| `mm_comparison_reason(status, reason, reason_code)` | Derives human-readable reason strings from status + reason_code. |
| `mm_table_col(table, col, default)` | Safe column accessor with default. |
| `mm_logical_col(table, col, default)` | Like `mm_table_col` but coerces via `isTRUE`. |
| `mm_comparison_id(formulas, target, method)` | Deterministic 8-hex checksum from formula+target+method string. |
| `mm_fit_status_label(fit)` | Reads `fit$fit_status` or falls back to optimizer certificate status. |
| `mm_lrt_stat(null, alternative)` | `pmax(0, deviance(null) - deviance(alternative))`. |
| `mm_drop_fixed_term_formula(fit, term)` | Reconstructs reduced formula string from artifact semantic IR; re-reads random term text verbatim from `artifact$semantic_model$random_terms[*]$source_syntax$text`. |

### 1.3 Rust FFI calls made

| R wrapper | Rust function | Schema |
|---|---|---|
| `mm_compare_models_json(payloads, method, "never")` | `mm_compare_models_json` | Returns `mixedmodels.model_comparison_table` v1.0.0 |
| `mm_bootstrap_lrt_json(...)` | `mm_bootstrap_lrt_json` | Returns `BootstrapLikelihoodRatioTest` JSON; top-level `p_value` field is certified presence |

The Rust `mm_bootstrap_lrt_json` entry point refits both models internally from bridge data (does not use a live handle), runs the parametric-bootstrap loop, accumulates replicate stats, and emits `BootstrapLikelihoodRatioTest` with `successful_replicates`, `completed_replicates`, `boundary_count`, `mcse`, and `seed_record`.

---

## 2. `R/compare-covariance.R`

### 2.1 Exported functions and S3 methods

#### `compare_covariance(spec)`
- **Exported:** yes (`export(compare_covariance)`)
- **Arguments:** `spec` — an `mm_spec` from `compile_model()`. Accepts any object that passes `mm_assert_compiled_spec()` (defined in `R/random-options.R`).
- **Does NOT accept a fitted `mm_lmm` directly** (the docstring says "or, in later phases, an mm_fit" but the body calls `mm_assert_compiled_spec` which enforces a compiled spec class).
- **Mechanism:** Calls `audit_design(spec)` to extract `random_term_cards` and `cross_card_constraints` from the upstream audit JSON. For each card, enumerates the three canonical covariance families (`full`, `diagonal`, `scalar`) using `mm_compare_covariance_card_rows()`.
- **Returns:** `mm_compare_covariance` object with fields `table` (data frame), `cards`, `cross_card_constraints`.
- **Table columns:** `term_id`, `group`, `basis`, `covariance_family`, `theta_parameters`, `assumes_zero`, `design_status`, `current` (logical — TRUE for the family matching the current fit).
- **No Rust FFI call directly** — relies on `audit_design()` which calls into Rust.
- **No TODO/FIXME markers.**

#### `print.mm_compare_covariance(x, ...)`
- **Exported:** yes (`S3method(print, mm_compare_covariance)`)
- **Behaviour:** Iterates rows; marks current family with `" <- current"`. Appends cross-card constraints section if any constraints present.

### 2.2 Internal helpers

| Helper | Purpose |
|---|---|
| `mm_compare_covariance_card_rows(card)` | Builds three data-frame rows (full/diagonal/scalar) for one random-term card. Uses `mm_covariance_label()` (from `R/explain.R`) and `mm_group_label()` (from `R/explain.R`). |
| `mm_covariance_theta_count(family, p)` | Returns `p*(p+1)/2` (full), `p` (diagonal), `1` (scalar, p>0), `0` (scalar, p=0). |
| `mm_covariance_zero_assumption(family, basis)` | Returns `"none"` (p≤1 or full), `"off-diagonal covariances"` (diagonal/scalar), `"unknown"` (unrecognized family). |

### 2.3 Design note
`current` detection has a subtle duplicate: the condition is
```r
identical(family, current) ||
  (identical(family, "diagonal") && identical(current, "diagonal"))
```
The second clause is redundant (covered by the first). No functional impact.

---

## 3. `R/changes.R`

### 3.1 Exported functions and S3 methods

#### `changes(fit, ...)` — generic
- **Exported:** yes (`export(changes)`)
- **Dispatches on class.**

#### `changes.mm_compiled(fit, ...)`
- **Exported:** yes (`S3method(changes, mm_compiled)`)
- **Argument:** `fit` — any `mm_compiled` (or subclass) that `mm_compiled_artifact()` can unwrap; in practice this includes `mm_lmm` and `mm_spec` objects, since both carry the artifact. However `mm_lmm` does NOT inherit `"mm_compiled"` — it inherits `"mm_lmm"`.
- **Gap:** Only `mm_compiled` dispatch is registered in NAMESPACE. There is no `changes.mm_lmm` method. An `mm_lmm` object does NOT inherit `"mm_compiled"`, so `changes(fitted_model)` will fall through to `UseMethod` default and error unless `mm_lmm` objects also carry the `"mm_compiled"` class in their class vector. (Confirmed: the test in `test-audit-verbs.R` line 46 calls `changes(fit)` on a compiled spec, not a fitted `mm_lmm`. The `test-phase2-revive.R` test at line 207 calls `changes(restored)` on a restored object — if that object is an `mm_lmm` without the `mm_compiled` class in its inheritance chain, this will fail.)
- **Mechanism:** Calls `mm_compiled_artifact(fit)` to extract the raw artifact list, then assembles rows from four sub-functions:
  1. `mm_change_formula_rows` — one row: requested formula vs effective formula; status `"unchanged"` or `"canonicalized"`.
  2. `mm_change_reduction_rows` — one row per entry in `artifact$reductions`; stage `"design_time_reduction"`.
  3. `mm_change_covariance_transition_rows` — one row per entry in `artifact$covariance_transitions`; stage `"covariance_transition"`.
  4. `mm_change_effective_covariance_rows` — one row per entry in `artifact$effective_covariance`; stage `"certificate_time"`. Reports requested vs fitted rank.
- **Returns:** `mm_change_log` list with fields `table` (data frame), `reductions`, `covariance_transitions`, `effective_covariance`, `fit_status`.
- **Table columns:** `stage`, `term_id`, `group`, `requested`, `effective`, `fitted`, `status`, `detail`.
- **No Rust FFI call** — reads from the already-materialised artifact JSON structure.
- **No TODO/FIXME markers.**

#### `print.mm_change_log(x, ...)`
- **Exported:** yes (`S3method(print, mm_change_log)`)
- **Behaviour:** Prints `"Model changes:\n"`, or `"  none recorded\n"` for empty table. Otherwise `print()` on the table without row names.

### 3.2 Internal helpers

| Helper | Purpose |
|---|---|
| `mm_change_formula_rows(artifact)` | Single-row: requested_formula vs effective_formula. |
| `mm_change_reduction_rows(artifact)` | Maps `artifact$reductions` list to data-frame rows; uses `from`/`requested` and `to`/`effective` field aliases. |
| `mm_change_covariance_transition_rows(artifact)` | Maps `artifact$covariance_transitions`; uses `from`/`requested_family` and `to`/`effective_family` aliases. |
| `mm_change_effective_covariance_rows(artifact)` | Maps `artifact$effective_covariance`; reads `requested_rank`, `supported_rank` for the detail string. Note: `effective` column is set equal to `requested` (both show `requested_basis`) — the "effective" basis is not separately tracked at certificate time. |
| `mm_change_empty_table()` | Zero-row template. |
| `mm_scalar_text(x, default)` | (from another file) — safe scalar character extractor. |
| `mm_list_text(x)` | (from another file) — safe character list flattener. |

### 3.3 Design note on `mm_change_effective_covariance_rows`
The `effective` column (line 122) is assigned `mm_list_text(summary$requested_basis)` — the same value as `requested`. This means the certificate-time row cannot show a delta in the effective basis, only in the scalar `fitted`/`status` field. This is a known structural constraint (the artifact at certificate time records rank, not a separate effective basis list), not a bug.

---

## 4. Cross-file observations

### 4.1 GLMM coverage
None of the three files contain any `mm_glmm` dispatch. There is no `compare.mm_glmm`, `anova.mm_glmm`, `drop1.mm_glmm`, or `changes.mm_glmm`. The `compare` generic has only an `mm_lmm` method. Users calling `compare(glmm_fit, ...)` will get a dispatch error.

### 4.2 `bootstrap_control()` (defined in `R/inference.R`, used here)
- **Exported:** yes.
- **Arguments:** `nsim` (default 999), `seed` (NULL or non-negative integer), `failed_refit_policy` (`"exclude"` | `"count_extreme"` | `"abort"`; default `"exclude"`).
- **Returns:** `mm_bootstrap_control` list with `requested_replicates`, `seed`, `failed_refit_policy`.

### 4.3 Test coverage summary
- `compare()` / `anova()` multi-model: exercised in `test-phase4.R`, `test-bw-lme-tutorial.R`, `test-brown-2021-lme-tutorial.R`, `test-codingclub-mixed-models.R`.
- `parametric_bootstrap()`: exercised in `test-phase4.R` (nsim=2, seed=101); REML refusal tested.
- `drop1()`: exercised in `test-phase4.R`.
- `anova()` single-model: exercised in `test-inference.R`.
- `compare_covariance()`: exercised in `test-random-options.R`, `test-no-advice.R`, `test-manifest.R`.
- `changes()`: exercised in `test-audit-verbs.R` (on compiled spec), `test-phase2-revive.R` (on restored object), `test-no-advice.R`.

---

## 5. Gaps

| # | Title | Severity | Classification | Evidence |
|---|---|---|---|---|
| G1 | `changes()` has no `mm_lmm` dispatch — `changes(fitted_lmm)` may silently fall through or error | major | in-scope-missing | NAMESPACE registers only `S3method(changes, mm_compiled)`; `mm_lmm` does not inherit `mm_compiled`. Test in `test-audit-verbs.R` line 46 calls on spec, not fitted model. |
| G2 | No `compare.mm_glmm` / `anova.mm_glmm` / `drop1.mm_glmm` | major | in-scope-missing | All three files are LMM-only; `compare(glmm_fit)` will dispatch-error. GLMM comparisons are a core use case not flagged as non-goal in PRD §3. |
| G3 | `compare_covariance()` accepts only `mm_spec` (compiled), not a fitted `mm_lmm` — docstring claims "or, in later phases, an mm_fit" | minor | partial | `mm_assert_compiled_spec` is called unconditionally; no `mm_lmm` pathway is present. Docstring claim is aspirational/deferred. |
| G4 | `drop1()` uses plain asymptotic chi-square for the LRT, with no boundary-aware correction for fixed-effect removal | minor | out-of-scope-by-design | Fixed-effect LRT does not involve a boundary parameter; standard chi-square is correct. Not a gap in practice but undocumented. |
| G5 | `mm_change_effective_covariance_rows` sets `effective = requested_basis` (not the certificate-time effective basis) | cosmetic | partial | Line 122 of `changes.R`; consequence is certificate-time rows never show a basis delta. Structural constraint of the artifact, not a crash risk. |
| G6 | `current` detection in `mm_compare_covariance_card_rows` has a redundant OR clause | cosmetic | works | `identical(family, "diagonal") && identical(current, "diagonal")` is subsumed by the preceding `identical(family, current)` check; no functional impact. |
