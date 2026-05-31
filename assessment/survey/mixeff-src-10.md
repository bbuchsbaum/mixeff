# mixeff source survey — file group 10: audit.R, explain.R, diagnostics.R, conditions.R

Survey date: 2026-05-31
Surveyor: subagent (claude-sonnet-4-6)

---

## 1. File overview

| File | Lines | Purpose |
|------|-------|---------|
| `R/audit.R` | 149 | `audit_design()` / `audit()` / `print.mm_audit` — delegates full audit rendering to the Rust FFI |
| `R/explain.R` | 289 | `explain_model()` / `print.mm_explanation` and ~15 private helpers that format `RandomTermCard` payloads in R |
| `R/diagnostics.R` | 301 | `diagnostics()` / `fit_status()` / `print.mm_diagnostics` + registries `mm_diagnostic_code_registry` and `mm_response_diagnostic_reason_registry` + guards |
| `R/conditions.R` | 77 | `mm_abort()` / `mm_split_tagged_error()` and documentation for all typed condition classes |

Class hierarchy used across these files:

- `mm_spec` inherits `mm_compiled` (set in `R/compile.R`)
- `mm_lmm` inherits `mm_fit`, `mm_compiled` (set in `R/fit-lmm.R`)

---

## 2. Exported functions and S3 methods

### 2.1 `audit_design(spec)` — `R/audit.R` l.59

**Exported.** Accepts an `mm_spec` or `mm_fit`; aborts with `mm_schema_error` otherwise.

Contract:
1. Extracts `attr(artifact, "raw_json")`.
2. Calls `.Call(wrap__mm_audit_report_text, raw_json)` → Rust `mm_audit_report_text` → `artifact.audit_report().to_text()`. Returns the display-rendered audit report string.
3. Calls `.Call(wrap__mm_audit_report_json, raw_json)` → Rust `mm_audit_report_json` → `serde_json::to_string(&artifact.audit_report())`. Returns the structured JSON.
4. Parses the JSON through `mm_json_parse_audit_report()` (in `R/json.R`), which validates schema headers and per-card schema headers via `mm_json_negotiate()`.
5. Returns `mm_audit` list with: `text`, `design_audit` (from artifact), `report` (parsed), `random_term_cards`, `cross_card_constraints`, `diagnostics` (report-level falling back to artifact-level).

Returns: S3 object `mm_audit`.

FFI calls:
- `wrap__mm_audit_report_text` (Rust: `fn mm_audit_report_text(artifact_json: &str) -> Result<String, String>`)
- `wrap__mm_audit_report_json` (Rust: `fn mm_audit_report_json(artifact_json: &str) -> Result<String, String>`)

What it refuses / NA's:
- No fallback if `raw_json` is missing or empty — raises `mm_schema_error` immediately. This is intentional: the no-silent-surgery contract requires every audit claim to trace to a JSON artifact.
- Rust errors are intercepted, the `mm_<tag>: message` prefix is stripped via `mm_split_tagged_error()`, and re-raised as the tagged class (falling back to `mm_bridge_error`).

Stubs / TODOs: none. Phase note in docstring: "Sections that depend on a fit (Optimizer / Inference) report `not assessed` until Phase 1.E lands `lmm()`" — this is now resolved in practice (lmm() exists), but the docstring was not updated.

---

### 2.2 `audit(fit, ...)` and `audit.mm_fit(fit, ...)` — `R/audit.R` l.132–140

**Exported.** Generic + single S3 method. `audit.mm_fit` delegates directly to `audit_design(fit)`. No `audit.mm_spec` method is registered (it is not needed because `mm_spec` also passes `inherits(spec, "mm_spec")` check in `audit_design()`).

NAMESPACE: `export(audit)`, `S3method(audit, mm_fit)`.

---

### 2.3 `print.mm_audit(x, ...)` — `R/audit.R` l.144

**S3.** Calls `cat(x$text)`, appends a trailing newline if not already present. Returns `x` invisibly.

NAMESPACE: `S3method(print, mm_audit)`.

---

### 2.4 `explain_model(spec)` — `R/explain.R` l.36

**Exported.** Accepts an `mm_spec` or `mm_fit`; aborts with `mm_schema_error` otherwise.

Contract:
1. Calls `audit_design(spec)` internally to obtain `mm_audit`.
2. Builds explanation text entirely in R using the `RandomTermCard` / diagnostic payloads already parsed by `audit_design()`. No additional FFI call.
3. Returns `mm_explanation` list: `text`, `cards`, `cross_card_constraints`, `diagnostics`, `report`.

The explanation text is assembled by these private helpers (all in `R/explain.R`):

| Helper | Role |
|--------|------|
| `mm_explanation_text(spec, audit)` | Top-level assembler: header, blocks, constraints, notes sections |
| `mm_spec_formula_text(spec)` | Pulls `effective_formula` → `requested_formula` → `deparse1(spec$formula)` |
| `mm_group_block_summaries(cards)` | Emits "X has N separate random-effect blocks" lines |
| `mm_card_lines(card)` | Formats one `RandomTermCard`: term_id, original/canonical fragment, per-block named form + scope + covariance + theta count, support |
| `mm_group_label(group)` | Handles `single$name` / `cell$names` / fallback |
| `mm_slopes_label(slopes)` | Flattens slope list to comma string |
| `mm_covariance_label(covariance)` | Collapses covariance tag list |
| `mm_support_label(support)` | Formats group levels, min/median rows per group |
| `mm_variation_label(variation)` | Formats within-group variation map |
| `mm_constraint_lines(constraint, cards)` | "card_A <-> card_B (basis_A <-> basis_B): reason" |
| `mm_design_note_lines(diagnostics)` | Filters bucket="design_note"; delegates to `mm_bucket_advice_lines` |
| `mm_fit_note_lines(diagnostics)` | Filters bucket="fit_note"; delegates to `mm_bucket_advice_lines` |
| `mm_bucket_advice_lines(diagnostics, bucket)` | Formats "code: message" lines, deduplicated |
| `mm_repair_lines(diagnostics)` | Formats bucket="repair" as numbered action list, deduplicated via `\x1f` key |
| `mm_singularity_lines(report)` | Scans report sections for "Effective Covariance", extracts lines whose `detail` mentions rank-deficient/effective rank |

Returns: S3 object `mm_explanation`.

FFI calls: none directly — all via `audit_design()`.

What it refuses / NA's:
- Falls through to "none" if `random_term_cards` is empty.
- `mm_singularity_lines` returns `character()` if the "Effective Covariance" section is absent (pre-fit spec).
- "Fit notes" and "Possible repairs" sections are silently omitted when no matching diagnostics exist.

Stubs / TODOs: none in-file. The Phase 1.C docstring reference is historical.

---

### 2.5 `print.mm_explanation(x, ...)` — `R/explain.R` l.59

**S3.** Same pattern as `print.mm_audit`: `cat(x$text)` + trailing newline guard.

NAMESPACE: `S3method(print, mm_explanation)` — **NOT registered in NAMESPACE.** The method is defined but there is no `S3method(print, mm_explanation)` line in NAMESPACE. `print.mm_explanation` will only be found if the package is attached (via `.S3method` generic dispatch on method names), not via `getS3method`. This is a minor packaging gap.

---

### 2.6 `diagnostics(fit, severity=NULL, stage=NULL, ...)` — `R/diagnostics.R` l.18

**Exported generic.** One registered method:

#### `diagnostics.mm_compiled(fit, severity, stage, ...)`

NAMESPACE: `S3method(diagnostics, mm_compiled)`. Because `mm_spec` and `mm_fit` both inherit `mm_compiled`, this single method covers both compiled specs and fitted models.

Contract:
1. Calls `mm_compiled_artifact(fit)` — validates `inherits(x, "mm_compiled")` and `is.list(x$artifact)`, aborts with `mm_schema_error` if not.
2. Calls `mm_artifact_diagnostics(artifact)` — concatenates three diagnostic lists: `artifact$diagnostics`, `artifact$design_audit$diagnostics`, `artifact$optimizer_certificate$diagnostics`.
3. Calls `mm_diagnostics_table(raw)` — builds a data.frame with columns `code`, `severity`, `stage`, `message`, `affected_terms`; each `affected_terms` is collapsed with `", "`. Passes through `mm_diagnostics_guard()`.
4. Applies `severity` and `stage` filters on the table's corresponding columns.
5. Returns `mm_diagnostics` list: `diagnostics` (filtered raw list), `table` (filtered data.frame), `severity`, `stage`.

Returns: S3 object `mm_diagnostics`.

FFI calls: none — reads from already-parsed artifact.

What it refuses / NA's:
- Empty diagnostics returns a zero-row data.frame (typed columns preserved).
- `mm_diagnostics_guard()` fires a session-scoped warn-once warning for any `code` not registered in `mm_diagnostic_code_registry` and attaches `attr(table, "mm_unrecognized_diagnostic_code")`.

Note: there is no `diagnostics.mm_fit` method — dispatch falls through to `diagnostics.mm_compiled` via the `mm_compiled` superclass.

---

### 2.7 `fit_status(fit, ...)` — `R/diagnostics.R` l.51

**Exported generic.** Two S3 methods:

#### `fit_status.mm_fit(fit, ...)`
Returns `fit$fit_status %||% fit$artifact$optimizer_certificate$status %||% "not_assessed"`.

#### `fit_status.mm_compiled(fit, ...)`
Returns `fit$artifact$optimizer_certificate$status %||% "not_assessed"`.

NAMESPACE: `export(fit_status)`, `S3method(fit_status, mm_fit)`, `S3method(fit_status, mm_compiled)`.

Returns: length-1 character string.

FFI calls: none.

---

### 2.8 `print.mm_diagnostics(x, ...)` — `R/diagnostics.R` l.71

**S3.** Prints "Diagnostics:\n  none" for empty tables. Otherwise:
- Selects and deduplicates columns `code`, `severity`, `stage`, `affected_terms` (intersecting with names present).
- Calls `print(..., row.names = FALSE)`.
- If `message` column exists, follows with a "Messages:" section: for each unique (code, message) pair, wraps the message to width 78 with 4-space exdent.

NAMESPACE: `S3method(print, mm_diagnostics)`.

---

## 3. Internal helpers in conditions.R

### `mm_abort(message, class, ..., call, parent)` — `R/conditions.R` l.52

Internal. Wraps `rlang::abort()`, prepending `class` with `"mm_condition"`. All typed mixeff errors inherit `mm_condition`. Not exported.

### `mm_split_tagged_error(message)` — `R/conditions.R` l.69

Internal. Parses the `"mm_<name>: <message>"` prefix that the Rust bridge attaches. Returns `list(tag, message)` — `tag` is `NA_character_` if no prefix found. Used in `audit.R` and other FFI call sites to route Rust errors to the right R condition class.

### Documented condition classes (by `@name mm-conditions` NULL-doc block):

| Class | Meaning |
|-------|---------|
| `mm_formula_error` | Formula parse / canonicalization failure |
| `mm_schema_error` | Wrong object type / missing required artifact fields / JSON parse failure |
| `mm_bridge_error` | Fallback for untagged Rust errors |
| `mm_data_error` | Data shape/type problems at compile_model() |
| `mm_fit_error` | Fit construction / optimization failure |
| `mm_inference_unavailable` | Inference, extractor, or prediction the engine cannot certify |
| `mm_not_identifiable` | Model not identifiable |
| `mm_design_refusal` | Structural refusal — design problem blocks fit |
| `mm_fit_not_optimized` | Fit was not fully optimized |
| `mm_arg_error` | Caller passed an invalid/malformed argument (distinct from domain refusals) |

All classes inherit `mm_condition`.

---

## 4. Registries in diagnostics.R

### `mm_diagnostic_code_registry` (l.157)

Maps every Rust `DiagnosticCode` variant (snake_case) to a bucket:

| Bucket | Role |
|--------|------|
| `"design_note"` | Rendered under "Design notes:" in `explain_model()` |
| `"repair"` | Rendered under "Possible repairs, not applied automatically:" |
| `"fit_note"` | Rendered under "Fit notes:" |
| `"raw_only"` | Not rendered on advice surfaces; visible via `diagnostics()$table` |

27 variants registered. Notable `raw_only` entries with rationale comments:
- `formula_canonicalized` — duplicative of `syntax_expansion` / `covariance_assumption`
- `optimizer_not_assessed` — pre-fit cosmetic
- `serialization_not_assessed` — pre-fit serialization state

The coverage gate in `tests/testthat/test-diagnostic-formatters.R` parses the vendored Rust source (`src/rust/upstream/mixeff-rs/src/compiler/diagnostics.rs`) and asserts that every enum variant is registered and no stale entries exist.

### `mm_response_diagnostic_reason_registry` (l.222)

Maps 5 `ResponseDiagnosticReason` variants (from `src/model/batch.rs`), all `raw_only`. Forward-compat slot for the GLMM batch path — no R surface currently exposes `ResponseColumnDiagnostic`. Coverage gate in the same test file.

### `mm_diagnostics_guard(table)` (l.259)

Forward-compat guard: warns once per session for any unrecognized `code` not in `mm_diagnostic_code_registry`. Uses session-scoped environment `mm_unknown_diag_state$seen` (initialized at package load) to suppress duplicate warnings.

### `mm_scalar_text(x, default="")` / `mm_list_text(x, default="")` (l.285–300)

Internal text-extraction helpers used when building the diagnostics table from mixed JSON structures (handles NULL, list, vector).

---

## 5. FFI call inventory

| R call site | Rust function | Rust return |
|-------------|---------------|-------------|
| `audit_design()` | `mm_audit_report_text(artifact_json)` | `Result<String, String>` — Display rendering of `ModelAuditReport` |
| `audit_design()` | `mm_audit_report_json(artifact_json)` | `Result<String, String>` — `serde_json::to_string(&artifact.audit_report())` |

Both Rust functions are registered in `src/rust/src/lib.rs` at lines 1585–1602 and 1875–1876. All other operations (explain_model, diagnostics, fit_status, conditions) operate purely on already-parsed R lists.

---

## 6. Test coverage summary

| Test file | What it covers |
|-----------|----------------|
| `test-compile-audit.R` | `audit_design()` round-trip, design_audit field, diagnostic codes exposed |
| `test-audit-verbs.R` | `audit()` post-fit alias, `diagnostics()` filtering, `fit_status()` |
| `test-explain-model.R` | `explain_model()` structure, 8 syntax patterns, repair rendering, R9 forbidden phrases, error on bad input |
| `test-explain-model-snapshots.R` | Snapshot tests for explanation text |
| `test-diagnostic-formatters.R` | Coverage gate: Rust enum vs. R registry (DiagnosticCode + ResponseDiagnosticReason) |
| `test-schema-versioning.R` | Verifies `mm_condition` class on schema errors |
| `test-manifest.R` | Capability manifest includes audit_design, explain_model, diagnostics, fit_status |

No dedicated test file for `conditions.R` internals (`mm_abort`, `mm_split_tagged_error`). These are indirectly exercised through error-path tests in the above files.

---

## 7. Gaps and observations

1. **`print.mm_explanation` not in NAMESPACE** (`minor` / `test-gap`): `print.mm_explanation` is defined in `R/explain.R` l.59 but there is no `S3method(print, mm_explanation)` registration in NAMESPACE. The method is found via informal dispatch when the package is attached but would not appear in `getS3method("print", "mm_explanation")` or be visible to `methods::existsMethod`. A `@method print mm_explanation / @export` roxygen tag is missing.

2. **Phase 1.A docstring not updated** (`cosmetic`): `audit_design()` docstring still says "Sections that depend on a fit (Optimizer / Inference) report `not assessed` until Phase 1.E lands `lmm()`" — `lmm()` has shipped; the note is stale.

3. **No `diagnostics.mm_fit` override** — intentional by design: `mm_fit` inherits `mm_compiled`, so `diagnostics.mm_compiled` covers it. Not a gap.

4. **`mm_response_diagnostic_reason_registry` — no R surface** (`out-of-scope-by-design`): All 5 variants are `raw_only` because `ResponseColumnDiagnostic` from the GLMM batch path has no R-facing surface yet. This is explicitly documented as a forward-compat slot (bd-01KRCKCYZ51H7H2WN8C5D7FNGT).

5. **No conditions.R unit tests** (`test-gap`): `mm_abort` and `mm_split_tagged_error` are exercised indirectly but have no dedicated test. In particular `mm_split_tagged_error`'s regex and the `tag = NA_character_` fallback path are only implicitly covered.

6. **`mm_unknown_diag_state` is process-global** (`minor`): The warn-once state is a package-level environment. If tests mutate it (by triggering an unknown code), subsequent tests in the same session will not re-warn for that code. No reset mechanism is exposed. This can cause test order dependency for the guard path.

7. **`mm_explanation_text` is pure R, not Rust-authored** — this is by design for the explanation surface, but means the R9 "no advice creep" contract relies on the registry buckets and `suggested_actions` fields from Rust rather than full delegation. The audit surface (`audit_design`) fully delegates; `explain_model` post-processes. Both are in scope.
