# mixeff source survey: fit-lmm.R, parse-formula.R, compile.R, parameterization.R

Survey date: 2026-05-31  
Family: mixeff-src-0  
Files covered:
- `R/fit-lmm.R`
- `R/parse-formula.R`
- `R/compile.R`
- `R/parameterization.R`
Supporting files read: `R/mm-control.R`, `R/data-translate.R`, `R/json.R`, `R/schema.R`,
`R/revive.R` (selected helpers), `R/methods-extract.R` (selected helpers),
`src/rust/src/lib.rs` (FFI implementations for `mm_parse_formula`, `mm_compile_model_json`,
`mm_fit_lmm_json`), `NAMESPACE`.

---

## 1. Exported functions and S3 methods

### 1.1 `mm_parse_formula(formula)` — `R/parse-formula.R`

**Export:** `export(mm_parse_formula)` in NAMESPACE.

**Contract:**  
Accepts a length-1 character string or an R `formula` object. Coerces R formulas to
character via `format()` + `paste(trimws(...), collapse=" ")`. Delegates to Rust FFI
`wrap__mm_parse_formula` → `mm_parse_formula(&str)`, which calls the upstream
`parse_formula()` parser and returns the canonical `Display` rendering. Returns a
length-1 character string (the canonical form). Equivalent formula spellings produce
identical output strings, enabling string-comparison equivalence-class tests.

**FFI call:** `.Call(wrap__mm_parse_formula, input)` → Rust `mm_parse_formula(formula: &str)`
→ `parse_formula(formula).map(|f| format!("{}", f))`.

**Error handling:**  
- Tagged bridge error `mm_formula_error` is caught and re-thrown as typed
  `mm_formula_error` condition with `formula` field.
- Untagged bridge errors become `mm_bridge_error`.

**Returns/refuses/NA:** Returns canonical string. Raises `mm_formula_error` on parse
failure (typed, catchable). Never returns NA.

**Helper (unexported):** `mm_coerce_formula_string(formula)` — shared by
`compile_model()` and `lmm()`. Validates type/length/non-empty; raises
`mm_formula_error` if not a character scalar or R formula.

**TODOs/stubs:** None.

---

### 1.2 `compile_model(formula, data)` — `R/compile.R`

**Export:** `export(compile_model)` in NAMESPACE.

**Contract:**  
Phase 0/1 pre-fit pipeline. Validates that `formula` is a two-sided R formula,
`data` is a data.frame, all formula variables exist in `data`, and no design variable
has NA values (`mm_check_no_na` — hard refusal, no silent row-dropping). Calls
`mm_translate_data()` to decompose columns into the FFI wire format, then calls Rust
`mm_compile_model_json` to run the semantic-IR / design-audit pipeline. Parses the
returned JSON via `mm_json_parse_artifact()` (schema-negotiated). Returns an `mm_spec`
/ `mm_compiled` object.

**FFI call:** `.Call(wrap__mm_compile_model_json, formula_string, column_order,
numeric_columns, categorical_values, categorical_levels)` → Rust
`mm_compile_model_json(formula, column_order, numeric_columns, categorical_values,
categorical_levels)` → `parse_formula` → `compile_formula_ir` →
`CompiledModelArtifact::new` + `attach_design_audit` → `serde_json::to_string`.

**Returns:**  
An S3 object of class `c("mm_spec", "mm_compiled")` with slots:
- `call` — matched call
- `formula` — input R formula
- `vars` — character vector of all formula variables (from `all.vars(formula)`)
- `model_frame` — narrowed `data[, vars]` data.frame
- `artifact` — parsed JSON artifact (list), schema `mixedmodels.compiled_model_artifact`
  v1; raw JSON preserved as `attr(artifact, "raw_json")`

**Error handling (typed conditions, all inherit `mm_condition`):**
- `mm_formula_error` — non-formula input or parse failure
- `mm_data_error` — non-data.frame, missing variables, NAs in design columns, or
  unsupported column types (anything other than numeric/integer/logical/factor/character)
- `mm_schema_error` — Rust artifact JSON fails schema negotiation

**S3 method registered:** `print.mm_spec` (registered via `S3method(print, mm_spec)` in
NAMESPACE; also exported). Prints formula, effective formula (from artifact), fixed term
count, random term count, schema name+version. Directs user to `audit_design()`.

**data-translate wire format (from `mm_translate_data()`):**
- numeric/integer/logical → `numeric_columns` (logicals coerced to 0/1)
- factor → `categorical_values` (as.character), `categorical_levels` (levels(col))
- character → `categorical_values`, `categorical_levels` (unique, first-appearance order)
- Date/POSIXct/list/etc → hard `mm_data_error`

**NA policy:** `mm_check_no_na` enumerates all offending columns and raises one
`mm_data_error` listing all. No silent omission.

**TODOs/stubs:** None found in the R layer. The artifact is produced entirely by Rust;
unsupported formula syntax (e.g., `I()`, `poly()`, `splines::ns()`, `offset()`,
`cbind(y1,y2)~`) would fail at the Rust parser/compiler, not here. The R layer does
not pre-screen for these.

---

### 1.3 `lmm(formula, data, REML, weights, control)` — `R/fit-lmm.R`

**Export:** `export(lmm)` in NAMESPACE.

**Contract:**  
Phase 1.E fit driver. Validates REML (must be length-1 logical non-NA), evaluates
`weights` in the caller's environment, validates them (positive finite numeric, length
== nrow(data)), calls `compile_model()` for the design audit, optionally prints
`explain_model()` (controlled by `control$verbose`), then calls Rust `mm_fit_lmm_json`
for the numerical optimization. Parses the returned fit payload JSON. Returns an
`mm_lmm` / `mm_fit` / `mm_compiled` S3 object.

**Arguments:**
- `formula` — two-sided lme4-style R formula
- `data` — data.frame
- `REML` — logical scalar, default `TRUE`
- `weights` — optional positive numeric case-weights vector (evaluated in `data`
  and caller's frame); `NULL` means unweighted; passed to Rust as empty `numeric()`
  when NULL
- `control` — `mm_control()` list; currently only `verbose` is used at the R layer;
  the full control JSON is passed to Rust but Rust parses it into `_control: Value`
  and does not currently act on any field (see §2.1 below)

**FFI call:** `.Call(wrap__mm_fit_lmm_json, formula_string, isTRUE(REML),
column_order, numeric_columns, categorical_values, categorical_levels, weights_vec,
control_json)` → Rust `mm_fit_lmm_json(...)` → `LinearMixedModel::new` →
`model.fit(reml)` → serializes `lmm_fit_result` JSON (schema
`mixeff.lmm_fit_result` v1).

**Rust fit payload fields extracted by `mm_json_parse_lmm_fit()`:**
`beta`, `beta_names`, `theta`, `sigma`, `log_likelihood`, `deviance`, `aic`, `bic`,
`nobs`, `dof`, `df_residual`, `fit_status`, `std_errors`, `fixed_fitted`, `fitted`,
`residuals`, `ranef`, `varcorr`, `fit_summary`, `artifact_json`, `optimizer`.

**Returns an `mm_lmm` list** (classes `c("mm_lmm", "mm_fit", "mm_compiled")`):

| Field | Type | Source |
|---|---|---|
| `call` | matched call | R |
| `formula` | R formula | R |
| `REML` | logical | R |
| `control` | mm_control list | R |
| `vars` | character vec | compile_model |
| `model_frame` | data.frame | compile_model |
| `weights` | numeric or NULL | R |
| `artifact` | parsed JSON list | Rust (post-fit artifact) |
| `fit` | raw fit JSON list | Rust |
| `fit_summary` | parsed fit_summary | Rust (FitSummaryPayload) |
| `schema` | list(name,version,…) | mm_object_schema(artifact) |
| `rust_handle` | NULL | always NULL at construction |
| `lazy_cache` | environment | mm_empty_lazy_cache() |
| `beta` | named numeric | Rust |
| `theta` | numeric | Rust |
| `sigma` | numeric | Rust |
| `logLik` | numeric | Rust |
| `deviance` | numeric | Rust |
| `AIC` | numeric | Rust |
| `BIC` | numeric | Rust |
| `nobs` | integer | Rust |
| `dof` | integer | Rust |
| `df_residual` | integer | Rust |
| `fit_status` | character | Rust (optimizer_certificate.status) |
| `std_errors` | named numeric | Rust |
| `fixed_effect_vcov` | matrix or NA-matrix | mm_fixed_effect_vcov_from_payload |
| `fixed_fitted` | numeric | Rust |
| `fitted` | numeric | Rust |
| `residuals` | numeric | Rust |
| `random_effects` | mm_ranef list | mm_ranef_from_terms |
| `varcorr` | mm_varcorr list | mm_varcorr_from_result |

**rust_handle is always NULL at construction.** The handle is treated as a rebuildable
cache. Inference verbs that need a live Rust handle call `revive()` lazily.

**Error handling:**
- `mm_fit_error` — REML validation, weights validation, or Rust-side fit/construction
  failure
- `mm_formula_error` — formula parse failure (from compile_model or fit)
- `mm_data_error` — NA in data, unsupported column type, weights mismatch
- `mm_schema_error` — JSON schema negotiation failure, missing fit_summary payload
- `mm_bridge_error` — untagged Rust panic

**TODOs/stubs:** None in R. Control JSON is parsed by Rust into `_control: Value`
(prefixed underscore = intentionally unused). No optimizer hyperparameters (tolerance,
max iterations, algorithm choice) are currently forwarded from R to the Rust engine.

---

### 1.4 `parameterization(fit, ...)` — `R/parameterization.R`

**Export:** `export(parameterization)` in NAMESPACE.  
**S3 dispatch:** `S3method(parameterization, mm_compiled)` in NAMESPACE.

**Contract:**  
Generic + `parameterization.mm_compiled` method. Reads the post-fit (or pre-fit,
for `mm_spec`) compiled artifact via `mm_compiled_artifact()` (defined in
`diagnostics.R`). Extracts `covariance_parameter_traces` and `theta_maps` from the
artifact. Converts each trace record to a data.frame row via
`mm_parameterization_trace_row()`.

**Returns:** An S3 object of class `mm_theta_map` with slots:
- `table` — data.frame with columns: `term_id`, `group`, `source_syntax`,
  `covariance_family`, `user_basis`, `optimizer_basis`, `theta_index`, `theta_name`,
  `theta_status`, `constraint`, `theta_value`, `lambda_row`, `lambda_col`,
  `lambda_row_basis`, `lambda_col_basis`, `lambda_value`, `varcorr_entries`.
  Empty table returned (correct structure, zero rows) when no traces present.
- `traces` — raw trace list from artifact
- `theta_maps` — raw theta_map list from artifact

**S3 method:** `print.mm_theta_map` (registered). Prints visible subset of columns
(term_id, group, source_syntax, covariance_family, theta_name, theta_value,
theta_status, varcorr_entries); reports count of hidden columns.

**Helper functions (unexported):**
- `mm_parameterization_trace_row(trace)` — converts one trace record to data.frame
  row; uses `%||%` and `NA` defaults throughout; gracefully handles absent fields.
- `mm_varcorr_entry_text(entries)` — formats varcorr entry list as a semicolon-
  separated string; returns `""` for empty.
- `mm_parameterization_empty_table()` — canonical empty data.frame with all 17
  columns typed correctly.

**FFI calls:** None. Reads from the artifact JSON that was already parsed at fit time.

**Accepts:** Any object inheriting `mm_compiled` — works for both `mm_spec`
(pre-fit, theta fields will mostly be NA) and `mm_lmm`/`mm_glmm` (post-fit).

**TODOs/stubs:** None.

---

## 2. Supporting infrastructure (shared helpers used by the four files)

### 2.1 `mm_control()` / `mm_validate_control()` — `R/mm-control.R`

`mm_control(verbose = 0L)` returns a list of class `mm_control`. Only `verbose`
is currently exposed. The control list is serialized to JSON and passed to Rust, but
Rust's `mm_fit_lmm_json` binds it into `let _control: Value` (underscore prefix —
intentionally unused). Optimizer hyperparameters (tolerance, maxit, algorithm, BOBYQA
vs Nelder-Mead, etc.) are not yet exposed at the R level.

### 2.2 `mm_translate_data()` — `R/data-translate.R`

Internal wire-format decomposer. Columns: numeric/integer/logical → `numeric_columns`;
factor → `categorical_values` + `categorical_levels` (canonical factor order);
character → same but levels derived first-appearance. All other types raise
`mm_data_error`. No NA tolerance — expects clean data.

### 2.3 `mm_check_no_na()` — `R/data-translate.R`

Hard NA refusal for design variables. Enumerates all offending columns; raises one
`mm_data_error` with `columns` and `na_counts` fields. Called by `compile_model()`.
Not called separately by `lmm()` (which defers to `compile_model()` for that check).

### 2.4 `mm_json_parse_artifact()` — `R/json.R`

Parses Rust artifact JSON. Validates schema via `mm_json_negotiate()` against the
closed schema set from `mm_json_known_schemas()`. Preserves raw JSON as
`attr(parsed, "raw_json")`.

### 2.5 `mm_fixed_effect_vcov_from_payload()` — `R/revive.R`

Constructs the `fixed_effect_vcov` matrix from the Rust-supplied payload.
If status is `"available"` and a matrix is present: builds full matrix, attaches
`mm_method`, `mm_status`, `mm_reliability`, `mm_reason`, `mm_details`, `mm_notes`,
`mm_schema_name`, `mm_schema_version` attributes. If status is `"unavailable"`:
returns an NA-filled matrix of correct dimension with `mm_unavailable_reason` attribute.
Strict schema enforcement — raises `mm_schema_error` for unknown status, available
without matrix, or unavailable with matrix.

### 2.6 `mm_ranef_from_terms()` / `mm_varcorr_from_result()` — `R/methods-extract.R`

`mm_ranef_from_terms` assembles a named list of data.frames (one per grouping factor)
from Rust ranef payload. Returns `mm_ranef` class. Handles multiple random terms per
group by cbind-merging.

`mm_varcorr_from_result` assembles VarCorr table from fit_summary.varcorr payload.
Returns list with `table` (data.frame) and `residual_sd`.

---

## 3. Rust FFI summary

| R-side name | Rust function | Input | Returns |
|---|---|---|---|
| `wrap__mm_parse_formula` | `mm_parse_formula(&str)` | formula string | canonical string |
| `wrap__mm_compile_model_json` | `mm_compile_model_json(...)` | formula + typed column lists | compiled artifact JSON |
| `wrap__mm_fit_lmm_json` | `mm_fit_lmm_json(...)` | formula + columns + REML + weights + control_json | lmm_fit_result JSON |

All three return `Result<String, String>` with tagged error prefixes (`mm_formula_error:`,
`mm_fit_error:`, `mm_schema_error:`, `mm_data_error:`). The R bridge catches the error
condition, splits the tag with `mm_split_tagged_error()`, and re-raises as a typed
`mm_condition`.

The fit Rust payload includes: `beta`, `beta_names`, `theta`, `sigma`,
`log_likelihood`, `deviance`, `aic`, `bic`, `nobs`, `dof`, `df_residual`, `fit_status`,
`std_errors`, `fixed_fitted`, `fitted`, `residuals`, `ranef`, `varcorr`, `fit_summary`,
`artifact_json`, `optimizer` block (backend, algorithm, return_value, feval, fmin, reml).

---

## 4. Gaps and stubs

### GAP-1: control_json not acted on by Rust (partial)

`mm_control()` only exposes `verbose`. Even `verbose` is not forwarded to the Rust
engine — it is used only in R to gate the `explain_model()` print. The full control
JSON is parsed into `_control: Value` in Rust (prefixed underscore, intentionally
unused at lines 355, 471, 1143 of `lib.rs`). No optimizer-tuning knobs
(tolerance, max-iterations, algorithm) are exposed to R users.

**Severity:** minor (optimizer uses Rust defaults; defaults appear to work for
standard cases). Would become major if users need to tune convergence for difficult
models.

### GAP-2: Formula transformation pre-screening absent (partial)

`mm_coerce_formula_string` and `compile_model` do not pre-screen for R-side formula
constructs that the Rust parser does not support: `I(x^2)`, `poly(x, 2)`,
`splines::ns(x, 3)`, `offset(log(n))`, `cbind(y1, y2) ~`, `log(y) ~`. These would
fail with a Rust-side `mm_formula_error` rather than a clear R-level message. The
`mm_formula_manifest` records `transformations` (implicit_intercept,
nested_grouping_expansion, interaction_grouping) but does not enumerate unsupported
patterns. No pre-flight guard at the R level.

**Severity:** minor (error still surfaces as typed `mm_formula_error`; message may be
less user-friendly than an R-level guard could provide).

### GAP-3: print.mm_spec uses effective_formula with fallback (cosmetic)

`print.mm_spec` prints `artifact$effective_formula %||% artifact$requested_formula`.
If neither field is present (schema evolution), the print line would show `NULL`.
No guard or informative fallback message.

### GAP-4: parameterization on pre-fit mm_spec returns NA-filled table (works by design)

When called on an `mm_spec` (pre-fit), theta values and lambda values are NA because
the artifact has no optimizer output. This is correct behavior — the empty/NA table
documents the parameterization structure before fitting. No user-visible warning that
results are pre-fit. Minor UX gap only.

### GAP-5: lmm() does not forward optimizer options to Rust (partial)

`mm_control()` has no parameters for optimizer algorithm, tolerance, or iteration
limits. The upstream Rust `LinearMixedModel::fit()` uses its own defaults. Users
who need BOBYQA vs Nelder-Mead switching, or tolerance tuning for near-singular
models, have no R-level knob. This is consistent with the current Phase 1 scope but
is a known gap relative to `lme4::lmerControl()`.

### GAP-6: No update() method for mm_lmm (out-of-scope-by-design for Phase 1)

The four surveyed files do not implement or reference `update.mm_lmm`. NAMESPACE
confirms no `S3method(update, mm_lmm)`. Standard `update()` would not work. Deferred
per PRD.

---

## 5. What is NOT present in these four files

- No `glmm()` logic (separate `R/glmm.R`).
- No inference/p-value logic (separate `R/inference.R`).
- No `predict()` logic (separate `R/predict.R`).
- No `simulate()` logic (separate `R/simulate.R`).
- No `refit()` logic (separate `R/revive.R`).
- No `bootstrap_control()` / `parametric_bootstrap()` (separate files).
- No `confint()`, `anova()`, `drop1()` implementations (separate files).
