# Survey: mixeff-src-9
## Files: R/revive.R, R/as-json.R, R/json.R, R/schema.R

Survey date: 2026-05-31  
Package: mixeff (R wrapper for mixeff-rs Rust crate)

---

## 1. Overview

These four files collectively implement:
- **Persistence/revival**: serialization state machine (`revive.R`, `as-json.R`)
- **JSON artifact parsing**: typed parsers for every artifact that crosses the Rust-R FFI (`json.R`)
- **Schema negotiation**: fast-fail version-gating before any artifact body is consumed (`schema.R`)

No TODO/FIXME/stub/"not yet"/"deferred" markers were found in any of the four files.

---

## 2. R/schema.R

### Exported functions

#### `mm_json_negotiate(header)`
- **Contract**: validates a `list(schema_name, schema_version)` pair against the closed set of known schemas. Raises `mm_schema_error` on any mismatch (missing fields, unknown name, version skew).
- **Returns**: invisibly `TRUE` on success.
- **Rust FFI**: delegates to `mm_json_negotiate_one(name, version)` via `.Call(wrap__mm_json_negotiate_one, ...)`. The Rust side owns the authoritative closed set; the R side re-raises structured errors.
- **Design intent**: fast-fail primitive — any code path that consumes a Rust artifact must call this before parsing, so version skew between crate and wrapper produces a single clean error rather than a field-by-field decode failure.
- **Refuses/NA**: none; always raises on mismatch rather than returning a status value.
- **Stubs/gaps**: none.

#### `mm_json_known_schemas()`
- **Contract**: returns a data frame `(name, version)` of all `(schema_name, schema_version)` pairs the current wrapper build understands.
- **Returns**: a two-column data frame with character columns `name` and `version`.
- **Rust FFI**: `.Call(wrap__mm_json_known_schemas)` — the Rust crate owns the canonical list.
- **Stubs/gaps**: none.

---

## 3. R/json.R

All functions in this file are internal (`@noRd`, not exported). They are the typed parsers for every artifact class that crosses the bridge.

### `mm_json_parse_artifact(json)`
- **Contract**: parses a single-string JSON artifact emitted by the Rust bridge. Uses `jsonlite::fromJSON(simplifyVector = FALSE)` to preserve nested mixed structures as nested lists. Validates the `schema` header via `mm_json_negotiate()`. Attaches the raw JSON string as `attr(parsed, "raw_json")` so downstream code can forward the bytes to other FFI calls without re-encoding.
- **Returns**: a nested list with `attr(.,"raw_json")` set.
- **Rust FFI**: none directly — parses the JSON string returned by fit/compile FFI calls.
- **Refuses/NA**: raises `mm_schema_error` on parse failure, missing `schema` header, or schema mismatch.
- **Stubs/gaps**: none.

### `mm_json_parse_audit_report(json)`
- **Contract**: same parse pattern as `mm_json_parse_artifact` but for audit report JSON, where schema fields are at top level rather than nested under `$schema`. Also validates the schema header on every nested `random_term_cards` entry.
- **Returns**: nested list with `attr(.,"raw_json")`.
- **Stubs/gaps**: none.

### `mm_json_parse_fixed_effect_inference_table(table)`
- **Contract**: parses the `fixed_effect_inference_table` sub-object from an already-parsed artifact. Validates the schema header. Iterates `table$rows` and calls `mm_fixed_effect_inference_row()` for each entry.
- **Returns**: an `mm_fixed_effect_inference_table` object with `$table` (data frame), `$raw`, `$schema_name`, `$schema_version`, `$crate_version`. Returns `NULL` if input is `NULL`.
- **Rust FFI**: none; works on already-parsed R list.
- **Row fields produced** (all honoring the "no fabrication" contract; missing values become `NA`): `term`, `label`, `kind`, `estimate`, `std_error`, `df`, `numerator_df`, `denominator_df`, `statistic`, `statistic_name`, `p_value`, `method`, `status`, `reliability`, `reliability_reason`, `reason`, `reason_code`, `reason_detail`, `estimability`, `details`, `notes`.
- **Stubs/gaps**: none.

### Internal helpers in json.R

| Helper | Purpose |
|---|---|
| `mm_fixed_effect_inference_row(row)` | Maps one Rust inference-table row to a single-row data frame. All fields guarded by `mm_optional_numeric`/`mm_optional_text`/`mm_scalar_text`. |
| `mm_fixed_effect_inference_empty_table()` | Returns a typed zero-row data frame with the full column schema; used when `rows` is empty. |
| `mm_optional_numeric(x)` | Returns first element coerced to numeric, or `NA_real_` if absent/NA/invalid. |
| `mm_optional_text(x)` | Returns first element as character, or `NA_character_` if absent/NA/invalid. |

---

## 4. R/as-json.R

### Exported functions

#### `as_json(x, pretty = FALSE, ...)`  (generic)
- **Contract**: serializes the public R-side state of a compiled spec or fitted object to a JSON string. Declared as the secondary persistence path; `saveRDS`/`readRDS` + `revive()` is primary.
- **Dispatch**: single registered method `as_json.mm_compiled` covering all `mm_compiled` objects (both `mm_fit` and compiled-but-unfitted specs).
- **Returns**: length-1 character string (JSON).
- **Rust FFI**: none; uses `jsonlite::toJSON`.

#### `as_json.mm_compiled(x, pretty = FALSE, ...)`
- **Serializes**:
  - `schema`: `{schema_name: "mixeff.r_object", schema_version: 1}`
  - `object_class`, `object_type` (`"fit"` or `"spec"`)
  - `formula` (deparsed)
  - `vars` (character vector)
  - `artifact_json` (raw JSON string from `attr(artifact, "raw_json")` or re-encoded artifact)
  - If `mm_fit`: `REML`, `beta`, `beta_names`, `theta`, `sigma`, `logLik`, `deviance`, `df_residual`, `fit_status`, `nobs`, `fitted`, `residuals`, plus `fit$schema`
- **What is NOT serialized** (notable omissions for round-trip fidelity):
  - `std_errors` — not included in `payload$fit`
  - `random_effects` / `ranef` — not included
  - `model_frame` — not included (prevents full lazy-cache reconstruction from JSON alone)
  - GLMM-specific fields (`family`, `link`, `n_agq`, `method`) — no branch for `mm_glmm`
- **Deserialization counterpart**: **none**. There is no `from_json()` or `as_json.default` recovery path. A JSON artifact can only be fed back to the Rust bridge; it cannot reconstruct an R `mm_fit` object for use with R-side extractors.
- **Refuses/NA**: none; serializes what it has.
- **Stubs/gaps**:
  - No `as_json.default` error method registered; calling on an unknown class falls through to `UseMethod` stop.
  - No `from_json()` complement — JSON round-trip is one-way only.
  - `std_errors` field omitted from serialized `fit` block.
  - GLMM-specific metadata not captured despite `mm_glmm` inheriting `mm_fit`.

---

## 5. R/revive.R

### Exported functions

#### `revive(fit, ...)` (generic)

#### `revive.mm_fit(fit, ...)`
- **Contract**: restores the process-local parts of an `mm_fit` after `readRDS()` or worker restart. The artifact and flat extractor values are the durable source of truth; the Rust handle is only a cache and is explicitly abandoned (`fit$rust_handle <- NULL`).
- **Steps**:
  1. Parses `fit$artifact` from `fit$fit$artifact_json` if not already a list (handles legacy layout).
  2. Validates that a parsed artifact list is present.
  3. Restores `attr(artifact, "raw_json")` from `fit$fit$artifact_json` if absent.
  4. Rebuilds `fit$schema` via `mm_object_schema(fit$artifact)` if absent.
  5. Sets `fit$rust_handle <- NULL`.
  6. Reinitializes `fit$lazy_cache` as a fresh empty environment.
  7. Ensures `"mm_fit"` and `"mm_compiled"` are in the class vector.
- **Returns**: an `mm_fit` with class `c(..., "mm_fit", "mm_compiled")`, live `lazy_cache` env, null Rust handle.
- **Rust FFI**: none.
- **Refuses/NA**: raises `mm_arg_error` if `fit` is not a list, or if no artifact can be recovered.
- **Stubs/gaps**: none.

#### `revive.default(fit, ...)`
- Raises `mm_arg_error` for non-mixeff inputs.

#### `fit_handle_alive(fit, ...)` (generic)

#### `fit_handle_alive.mm_fit(fit, ...)`
- **Contract**: tests whether the process-local Rust handle (`fit$rust_handle`) is a live `externalptr`.
- **Returns**: length-1 logical.
- **Note**: `FALSE` does not mean the fit is unusable; Phase 2 extractors read from the durable artifact.
- **Stubs/gaps**: none.

#### `fit_handle_alive.default(fit, ...)`
- Returns `FALSE` unconditionally.

#### `getME(object, name, ...)` (generic)

#### `getME.mm_lmm(object, name, ...)`
- **Contract**: lme4-compatible low-level extractor for the named component(s). Components are rebuilt lazily from the serialized R object (no live Rust handle required).
- **Supported names**: `X`, `Z`, `Zt`, `Lambda`, `Lambdat`, `theta`, `beta`, `fixef`, `y`, `mu`, `flist`, `cnms`.
- **Returns**: the single component, or a named list for multiple names.
- **Rust FFI**: none; all components rebuilt from `artifact` + `model_frame` using R-side lazy helpers.
- **Lazy cache keys**: `X` → `mm_fixed_model_matrix`, `Z` → `mm_random_model_matrix`, `Lambda` → `mm_lambda_matrix`, `flist` → `mm_random_flist`, `cnms` → `mm_random_cnms`.
- **Refuses/NA**: raises `mm_arg_error` for unsupported component names.
- **Stubs/gaps**:
  - Only registered for `mm_lmm`; `getME.mm_glmm` is not registered. Calling `getME()` on an `mm_glmm` object dispatches to `getME.default` which raises an error. This is a gap because `mm_glmm` inherits `mm_fit` and shares the same artifact structure.
  - `getME.default` raises `mm_arg_error` (intentionally).

#### `model.matrix.mm_lmm(object, type = c("fixed","random"), ...)`
#### `model.matrix.mm_glmm` (alias for `model.matrix.mm_lmm`)
- **Contract**: returns fixed-effect design matrix (type `"fixed"`) or sparse random-effect design matrix (type `"random"`). Both are rebuilt from the artifact without a Rust handle.
- **Returns**: dense matrix for `"fixed"`, `Matrix::dgCMatrix` for `"random"`.
- **Stubs/gaps**: none; GLMM shares the same implementation via the alias.

#### `vcov.mm_lmm(object, type = c("fixed","theta"), ...)`
#### `vcov.mm_glmm` (alias for `vcov.mm_lmm`)
- **Contract**:
  - `type = "fixed"`: returns the fixed-effect covariance matrix, preferring `object$fixed_effect_vcov`, then decoding `artifact$fixed_effect_covariance_matrix`, then falling back to a diagonal matrix of stored SEs.
  - `type = "theta"`: returns an all-`NA` matrix with `mm_unavailable_reason = "theta_covariance_unavailable"` — theta covariance is explicitly refused.
- **Returns**: named numeric matrix with `mm_*` attributes (`mm_method`, `mm_status`, `mm_reliability`, `mm_reason`, `mm_schema_name`, `mm_schema_version`, etc.).
- **Rust FFI**: none; reads from stored artifact.
- **Stubs/gaps**:
  - Theta covariance is permanently refused (NA matrix with reason attribute). This matches the PRD out-of-scope designation.
  - `mm_fixed_effect_vcov_from_payload` validates the payload strictly, raising `mm_schema_error` on shape, symmetry, finiteness, or name mismatch violations.

#### `random_blocks(fit, ...)` (generic)

#### `random_blocks.mm_compiled(fit, ...)`
- **Contract**: summarizes random-effect block structure from `artifact$design_audit$random_terms`.
- **Returns**: `mm_random_blocks` with `$table` (data frame: `term_id`, `group`, `basis`, `covariance`, `theta_parameters`, `group_levels`, `min_rows_per_group`, `median_rows_per_group`, `status`, `reason`) and `$random_terms`, `$semantic_terms` (raw lists).
- **Rust FFI**: none.
- **Stubs/gaps**: none.

#### `optimizer_certificate(fit, ...)` (generic)

#### `optimizer_certificate.mm_compiled(fit, ...)`
- **Contract**: extracts `artifact$optimizer_certificate` and presents it as a metric/value table.
- **Returns**: `mm_optimizer_certificate` with `$raw` and `$table` (metrics: `status`, `optimizer`, `objective`, `iterations`, `free_gradient_norm`, `projected_gradient_norm`, `hessian_eigen_min`, `hessian_rank`, `information_rank`).
- **Rust FFI**: none.
- **Stubs/gaps**: none.

#### `inference_table(fit, ...)` (generic)

#### `inference_table.mm_lmm(fit, method = c("auto","satterthwaite","kenward_roger","asymptotic","none"), ...)`
- **Contract**:
  - `method = "auto"`: returns the artifact-cached table resolved by the engine at fit time via `mm_json_parse_fixed_effect_inference_table`. Falls back to a table of `NA` inference values with `reason = "fixed_effect_inference_table_unavailable_legacy_object"` for legacy objects missing the field.
  - `method != "auto"`: recomputes by dispatching one `contrast()` per fixed-effect coefficient with the requested method (via `mm_inference_table_recompute`), ensuring user method requests are honored rather than silently replaced.
- **Returns**: `mm_inference_table` with `$table` and `$raw`.
- **Rust FFI** (indirect, via `contrast()`): re-fitting path uses `mm_fixed_effect_contrast_json` for non-auto methods.
- **Refuses/NA**: legacy objects get fully NA inference rows rather than an error, with an honest `reason` string.
- **Stubs/gaps**:
  - Only registered for `mm_lmm`; **no `inference_table.mm_glmm`** is registered. GLMM inference is not available through this extractor.

#### `reproducibility(fit, ...)` (generic)

#### `reproducibility.mm_compiled(fit, ...)`
- **Contract**: exposes `artifact$reproducibility` metadata: `fit_intent`, `random_state_used`, and threshold table.
- **Returns**: `mm_reproducibility` with `$raw` and `$thresholds` (data frame of name/value pairs).
- **Rust FFI**: none.
- **Stubs/gaps**: none.

#### `is_singular(x, tol = 1e-4, ...)` (generic)

#### `is_singular.mm_lmm(x, tol = 1e-4, ...)`
- **Contract**: returns `TRUE` if `fit_status()` indicates a boundary/reduced-rank convergence, or if any `artifact$effective_covariance` entry has `status` in `{boundary, reduced_rank, singular}`. The `tol` argument is accepted for lme4 compatibility but is not used (mixeff classifies singularity at fit time via the Rust engine).
- **Returns**: length-1 logical.
- **Stubs/gaps**:
  - `tol` argument is silently ignored; tolerance-based computation is deferred to the Rust engine.
  - Only registered for `mm_lmm`; `is_singular.mm_glmm` is not registered.

### Internal lazy-cache helpers in revive.R

| Helper | Purpose |
|---|---|
| `mm_empty_lazy_cache()` | Creates a fresh `environment(parent = emptyenv())` for caching. |
| `mm_object_schema(artifact)` | Extracts `schema_name`, `schema_version`, `crate_version`, `package_version` from artifact header. |
| `.mm_lazy(fit, key, producer)` | Cache-aside pattern: checks `fit$lazy_cache` env for `key`, calls `producer(fit)` on miss, stores and returns. Calls `revive(fit)` if `lazy_cache` is not an environment. |
| `mm_fixed_model_matrix(fit)` | Rebuilds fixed-effect design matrix from artifact semantic model + stored model frame. |
| `mm_fixed_formula(fit)` | Reconstructs the fixed-effects-only formula from `artifact$semantic_model$fixed_terms`. |
| `mm_response_name(fit)` | Extracts response variable name from artifact, falling back to formula. |
| `mm_response_vector(fit)` | Extracts response column from `fit$model_frame`. |
| `mm_random_model_matrix(fit)` | Rebuilds sparse random-effect Z matrix by iterating `artifact$semantic_model$random_terms`. |
| `mm_random_term_matrix(fit, term, index)` | Builds one block of Z for a single random term (group × basis Kronecker product). |
| `mm_random_term_group_label(fit, term, index)` | Resolves grouping factor name from `design_audit` then semantic model. |
| `mm_group_factor(frame, group_label)` | Reconstructs grouping factor; supports interaction terms via `:` split. |
| `mm_basis_label(basis)` | Resolves display label for a random-effect basis element. |
| `mm_basis_values(basis, frame)` | Retrieves numeric values for a random-effect basis column from the stored model frame. Raises `mm_inference_unavailable` for non-numeric basis columns. |
| `mm_random_flist(fit)` | Builds `mm_flist` from all random terms. |
| `mm_random_cnms(fit)` | Builds `mm_cnms` (coefficient names per grouping factor) via `ranef()`. |
| `mm_lambda_matrix(fit)` | Rebuilds the relative covariance factor Lambda from stored theta. Handles `full`, `diagonal`/`diag`, and `scalar` families. |
| `mm_random_block_row(term)` | Extracts one row of the random-blocks table from a design-audit term entry. |
| `mm_random_blocks_empty_table()` | Returns typed zero-row data frame for the blocks table schema. |
| `mm_optimizer_certificate_table(cert)` | Flattens certificate fields to metric/value rows. |
| `mm_repro_threshold_table(thresholds)` | Flattens reproducibility threshold list to name/value data frame. |
| `mm_inference_table_recompute(fit, method)` | Dispatches one `contrast()` per coefficient to recompute inference with the requested method. |
| `mm_inference_row_from_contrast(ct, term)` | Reshapes a single-row contrast table to the inference-table row schema. |
| `mm_fixed_effect_vcov_from_payload(payload, beta, std_errors)` | Decodes `fixed_effect_covariance_matrix` payload; handles available/unavailable statuses and validates shape, symmetry, finiteness, and coefficient names. Falls back to diagonal from SEs if no payload present. |
| `mm_validate_fixed_effect_vcov_payload(V, payload_names, coef_names, payload)` | Strict validator for covariance payload; raises `mm_schema_error` on any contract violation. |
| `mm_numeric_matrix_from_rows(rows)` | Converts row-list or matrix to plain numeric matrix. |

---

## 6. Rust FFI surface (only calls made from these four files)

| R wrapper | FFI symbol | Called from |
|---|---|---|
| `mm_json_negotiate_one(name, version)` | `wrap__mm_json_negotiate_one` | `schema.R` (via `mm_json_negotiate`) |
| `mm_json_known_schemas()` | `wrap__mm_json_known_schemas` | `schema.R` (exported) |

All other work in these four files is pure R operating on already-parsed artifact lists. The heavy FFI calls (fit, contrast, compile) live in `fit-lmm.R`, `inference.R`, etc.

---

## 7. Gap inventory

### G1 — `getME` not registered for `mm_glmm` (minor)
`getME.mm_lmm` is registered but `getME.mm_glmm` is not. `mm_glmm` inherits `mm_fit` and carries the same artifact structure; there is no structural reason the extractor cannot work. Calling `getME()` on a GLMM dispatches to `getME.default` and raises an error.

### G2 — `inference_table` not registered for `mm_glmm` (major)
`inference_table.mm_lmm` exists; there is no `inference_table.mm_glmm`. GLMM inference is not accessible through this extractor. The GLMM class uses Wald-z inference (see `summary.mm_glmm`), but `inference_table()` cannot surface it.

### G3 — `is_singular` not registered for `mm_glmm` (minor)
`is_singular.mm_lmm` exists; no `is_singular.mm_glmm`. Dispatches to `is_singular.default` and errors on GLMM objects.

### G4 — `as_json` omits `std_errors` from serialized fit block (minor)
`std_errors` is a first-class field on `mm_fit` objects and is used in the `vcov` fallback path, but it is not written into `payload$fit` by `as_json.mm_compiled`. A JSON-deserialized fit (if a `from_json` path ever existed) would lose SE values.

### G5 — `as_json` omits GLMM-specific fields (minor)
`mm_glmm` inherits `mm_fit` and `mm_compiled`, so `as_json.mm_compiled` fires for GLMM objects. However, GLMM-specific metadata (`family`, `link`, `n_agq`, `method`) is not captured in the payload. The serialized JSON cannot identify whether the fit was LMM or GLMM beyond `object_class`.

### G6 — No `from_json()` / round-trip deserialization (partial)
`as_json()` is documented as a secondary persistence path but there is no `from_json()` or equivalent. A JSON string produced by `as_json()` cannot be deserialized back to a usable R `mm_fit` object. `saveRDS`/`readRDS` + `revive()` is the only working round-trip. This is noted in the docs but represents an incomplete persistence interface.

### G7 — `model_frame` not serialized by `as_json` (minor)
`as_json` does not include `fit$model_frame`. Without the model frame, lazy-cache reconstruction (X, Z, Lambda, flist, cnms, y) cannot proceed from JSON alone, even if a `from_json` were added. This compounds G6.

### G8 — `tol` argument to `is_singular` silently ignored (cosmetic)
Accepted for lme4 API compatibility but has no effect. Not documented in the Rd as ignored.

### G9 — `as_json.default` not registered (cosmetic)
Calling `as_json()` on a non-`mm_compiled` object dispatches to `UseMethod` stop rather than a structured `mm_arg_error`. Minor ergonomic inconsistency compared to other generics in this package.
