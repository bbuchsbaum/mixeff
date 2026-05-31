# Survey: R/reporting.R

**Family:** mixeff-src-8
**File:** `/Users/bbuchsbaum/code/mixeff/R/reporting.R`
**Lines:** 843
**Date surveyed:** 2026-05-31

---

## Overview

`reporting.R` implements the publication-oriented reporting surface of mixeff.
Two exported generics — `model_report()` and `reporting_table()` — assemble
structured, provenance-tagged data frames from a fitted model's in-memory
artifact. The design is explicitly audit-first: every cell in every section
carries a `source`, `status`, and `reason` column; a parallel `$unavailable`
ledger records fields that could not be populated and why.

No Rust FFI calls are made directly in this file. All Rust artifact data is
accessed through the R-side wrapper objects already attached to `fit`
(`fit$artifact`, `fit$fit_summary`, `fit$schema`, etc.) or through calls to
other mixeff verbs (`inference_table()`, `VarCorr()`, `audit_design()`,
`random_blocks()`, `changes()`, `optimizer_certificate()`, `reproducibility()`).

---

## Exported API

### `model_report(fit, sections = "all", ...)`

**Generic** dispatched on `fit`.

**Implemented methods:**
- `model_report.mm_fit` — the primary method (described below).
  No method exists for `mm_glmm` beyond inheriting `mm_fit`.

**Contract:**
1. Validates `sections` via `mm_report_sections_arg()` — accepts `"all"` or
   any subset of the 11 known section names; errors with `mm_schema_error` on
   unknown names.
2. Runs all 11 section builders (each `fit → list(table, source, unavailable)`).
3. Collapses `unavailable` rows across all sections.
4. Returns an `mm_model_report` list with four fields:
   - `$metadata` — data frame: `created_at`, `package_version`,
     `crate_version`, `artifact_schema`, `artifact_schema_version`,
     `fit_class`, `sections`.
   - `$sections` — named list of 11 data frames (one per section).
   - `$unavailable` — ledger data frame recording fields that could not be
     populated, with columns `section`, `field`, `status`, `reason`,
     `source`, `action_taken`.
   - `$provenance` — data frame: one row per section, recording which
     upstream source the section came from.

**Known issues / gaps:**
- No `model_report.mm_glmm` method is registered; `mm_glmm` inherits
  `mm_fit` implicitly and receives the LMM-oriented overview label
  (`"LMM"` vs `class(fit)[[1L]]` branch). The `model_class` field in the
  overview will show `"mm_glmm"` (the generic `class(fit)[[1L]]` branch),
  which is correct but undocumented.
- The `sections` argument filtering logic builds *all* sections first and
  then subsets by name — wasteful when only one section is requested — but
  this is a performance issue, not a correctness issue.

---

### `reporting_table(fit, section = "all", view = c("compact", "audit"), ...)`

**Generic** dispatched on `fit`.

**Implemented methods:**

| Method | `fit` class | Notes |
|--------|-------------|-------|
| `reporting_table.mm_fit` | `mm_fit` (LMM/GLMM) | builds full report then extracts section |
| `reporting_table.mm_model_report` | `mm_model_report` | operates on already-built report, avoids rebuild |
| `reporting_table.mm_model_comparison` | `mm_model_comparison` | delegates to `mm_report_comparison_object_table()` |
| `reporting_table.mm_drop1` | `mm_drop1` | same as above |
| `reporting_table.mm_random_effect_test` | `mm_random_effect_test` | direct table extraction, no section builder |

**`view` argument:**
- `"compact"` — retains only human-readable columns per section (defined in
  `mm_report_view_table()`); drops `source`, `reason`, `details`, `notes`.
  Sets `attr(out, "view") = "compact"` and records dropped column names in
  `attr(out, "audit_columns")`.
- `"audit"` — returns the full table unchanged.

**Return values:**
- Single section → data frame.
- `section = "all"` → named list of data frames.

**`mm_report_comparison_object_table()` (internal helper for comparison objects):**
- Accepts `mm_model_comparison` and `mm_drop1`.
- Only exposes the `"comparison_ledger"` section; requesting any other section
  raises `mm_schema_error`.
- Passes `fit$ledger %||% mm_comparison_ledger_empty()` through
  `mm_report_section()`.

**`reporting_table.mm_random_effect_test` specific behaviour:**
- Does not go through the section builder pipeline at all.
- Directly uses `fit$table`.
- `"compact"` view: keeps columns `term`, `group`, `statistic`,
  `statistic_name`, `p_value`, `reference_distribution`, `method`, `status`,
  `reason_code` (intersection of those names with whatever is present).
- `"audit"` view: returns `fit$table` unchanged.

---

## Section Builders (all internal, not exported)

Called from `mm_report_builders()` which returns a named list.

### `mm_report_overview(fit)`
**Data sources:** `mm_compiled_artifact(fit)`, `inference_table(fit)$table`,
`fit$schema`, `fit$formula`, `fit$REML`, `fit$nobs`, `fit$model_frame`.

**Fields produced:** `model_class`, `formula`, `effective_formula`,
`fit_method` (REML/ML), `mode` (from `artifact$reproducibility$fit_intent`),
`nobs`, `fit_status`, `inference` (`"N/M available fixed-effect rows"`),
`artifact_schema`, `crate_version`, `package_version`.

**Fallbacks:** If `artifact$effective_formula` is NULL, falls back to
`artifact$requested_formula`, then to `deparse1(fit$formula)`.
`mode` falls back to `"not_recorded"` if the Rust artifact does not carry
`reproducibility$fit_intent`.

---

### `mm_report_model_specification(fit)`
**Data sources:** `mm_compiled_artifact(fit)`, `changes(fit)$table`,
`fit$call`, `fit$formula`, `artifact$semantic_model`.

**Fields produced:** `call`, `formula`, `requested_formula`,
`effective_formula`, `fixed_terms`, `random_terms`.
If `changes(fit)` returns non-empty rows, appends one row per change as
`change:N` fields.

---

### `mm_report_data_design(fit)`
**Data sources:** `random_blocks(fit)$table`, `audit_design(fit)`,
`fit$model_frame`.

**Fields produced:** `group`, `role`, `group_levels`,
`min_rows_per_group`, `median_rows_per_group`, `max_rows_per_group`,
`status`, `reason`.

`max_rows_per_group` is computed live from the model frame via
`mm_group_factor()` + `tabulate()`. If that call errors, the column
is silently set to `NA_integer_`.

`role` is looked up from `audit_design(fit)$random_term_cards[*]$role_origin$role`.
Falls back to `"not_recorded"` if not present.

If `random_blocks(fit)` returns an empty table, returns an empty-schema data
frame rather than erroring.

---

### `mm_report_random_terms(fit)`
**Data sources:** `audit_design(fit)$random_term_cards`,
`audit_design(fit)$cross_card_constraints`.

**Fields produced** (per block within each card): `term_id`,
`original_fragment`, `canonical_fragment`, `group`, `block`, `basis`,
`intercept`, `slopes`, `covariance`, `theta_parameters`, `english`,
`constraints`, `design_status`.

Cross-card constraints add extra rows with `design_status = "constraint"`
and `source = "cross_card_constraints"`.

If `audit_design(fit)` returns no cards, returns an empty-schema data frame.

---

### `mm_report_random_effects(fit)`
**Data sources (primary path):** `fit$fit_summary$varcorr` (schema
`mixedmodels.fit_summary`), parsed via `mm_varcorr_from_result()`.

**Data sources (fallback path):** `VarCorr(fit)` if
`fit$fit_summary$varcorr` is not a list.

**Fields produced:** `group`, `term_id` (always `NA`), `basis_lhs`,
`basis_rhs`, `kind` (`"variance"` for RE rows, `"residual_variance"` for
the residual row), `variance`, `std_dev`, `correlation`,
`covariance_family`, `status`, `reason`.

**Residual row:** Appended when `vc_obj$residual_sd` is not `NA`.

**Unavailable ledger entry:** When the fallback path is used (i.e., no
`fit_summary$varcorr`), records:
```
section = "random_effects"
field   = "stable_random_effect_variance_covariance_payload"
status  = "schema_gap"
reason  = "using_fit_varcorr_until_rust_report_payload_is_available"
source  = "planning/reporting_artifact_requirements.md"
```
This entry is *absent* when the primary (fit_summary) path is used.
The test `test-reporting.R:140` confirms the primary path is taken in
practice for standard LMM fits.

---

### `mm_report_fixed_effects(fit)`
**Data sources:** `inference_table(fit)$table`.

Passes the table through without transformation other than adding
`source = "fixed_effect_inference_table"`. All per-row `status`,
`reason`, `method`, `reliability`, etc. are Rust-authored or set by
the inference layer.

---

### `mm_report_fit_statistics(fit)`
**Data sources:** `fit$logLik`, `fit$deviance`, `fit$AIC`, `fit$BIC`,
`fit$nobs`, `fit$df_residual`, `fit$sigma` — all R-side slots, no Rust
call.

**Fields produced:** KV table with `logLik`, `deviance`, `AIC`, `BIC`,
`nobs`, `df_residual`, `sigma`. If any slot is NULL/NA it will appear as
`NA` in the value column but no special handling applies.

---

### `mm_report_optimizer(fit)`
**Data sources:** `optimizer_certificate(fit)$table`.

Passes through with `source = "OptimizerCertificate"` added.

---

### `mm_report_comparison_ledger(fit)`
Always returns an *empty* table (zero rows) with the full schema:
`comparison_id`, `model_id`, `formula`, `fit_method`, `refit`,
`comparison_method`, `statistic`, `df`, `p_value`, `status`, `reason`.

Records an unavailable entry:
```
status = "not_applicable"
reason = "no_model_comparison_recorded_on_this_fit"
```
This is by design: comparisons are recorded on `mm_model_comparison`
objects, not on individual fits.

---

### `mm_report_reproducibility(fit)`
**Data sources:** `reproducibility(fit)`, `R.version`, `getOption("contrasts")`.

**Fields produced:** `fit_intent`, `random_state_used`, `mixeff_version`,
`r_version`, `platform`, `contrasts`. If `reproducibility(fit)$thresholds`
is non-empty, appends one row per threshold as `threshold:<name>`.

`fit_intent` and `random_state_used` fall back to `"not_recorded"` if absent.

---

## Internal Infrastructure

### Section normalization (`mm_report_normalize_table`)
Ensures every section table has `source`, `status`, `reason` columns.
- Missing `source`: filled with the section's source string.
- Missing `status`: filled with `"available"`.
- Missing `reason`: filled with `NA_character_`.
- Non-available rows with no reason: backfilled with `"reason_not_recorded"`.

### `mm_report_view_table(table, section, view)`
Maps each section name to a fixed ordered list of compact columns.
Calls `mm_drop_empty_report_columns()` to remove columns where every value
is NA, blank, or `"reason_not_recorded"`.
Attaches `attr(out, "view") = "compact"` and
`attr(out, "audit_columns")` (names dropped from the full table).

### `mm_drop_empty_report_columns(table)`
Drops columns that are entirely empty/NA/blank/"reason_not_recorded".
Handles list columns (checks `lengths > 0`), numeric/logical (checks `!is.na`),
and character (checks `nzchar` and not `"reason_not_recorded"`).

### `mm_report_sections_arg(sections, allow_many = TRUE)`
Validates the `sections`/`section` argument.
- `"all"` expands to all 11 known names when `allow_many = TRUE`, or
  remains as `"all"` when `allow_many = FALSE`.
- Unknown names → `mm_schema_error`.
- Multiple names when `allow_many = FALSE` → `mm_schema_error`.

### `mm_report_kv(field, value, source, status, reason)`
Constructs a key-value data frame. `source`, `status`, `reason` are
recycled to `length(field)`.

### `mm_report_unavailable(section, field, status, reason, source)`
Constructs a single-row unavailable ledger entry (always
`action_taken = "reported"`).

### `mm_report_unavailable_empty()`
Returns a zero-row unavailable data frame with the required schema
(6 character columns).

### `mm_report_metadata(fit, sections)`
Builds the `$metadata` data frame. Reads `fit$schema` (or derives it
via `mm_object_schema(fit$artifact)`); records `created_at` (live
wall-clock time), `package_version`, `crate_version`, `artifact_schema`,
`artifact_schema_version`, `fit_class`, `sections`.

### `mm_report_random_term_card_rows(card)`
Converts one random-term card (from Rust `ModelAuditReport`) to a list
of data frames (one per block). If the card has no `blocks`, creates one
empty block row.

### `mm_report_cross_card_rows(constraints)`
Converts cross-card constraint entries to data frame rows with
`design_status = "constraint"`.

### `mm_report_constraint_text(constraints)`
Serialises implied-constraint list to a semicolon-separated string like
`"reason [between_basis]"`.

### `mm_group_ir_label(group)`
Resolves a Rust group IR value to a string. Tries four path patterns:
`name`, `single$name`, `interaction$names`, `cell$names`; falls back to
`mm_scalar_text(group)`.

### `mm_report_group_roles(fit, groups)` / `mm_report_group_max_rows(fit, group)`
Helper functions for `mm_report_data_design()`.
`mm_report_group_max_rows` wraps `mm_group_factor()+tabulate()` in
`tryCatch`, returning `NA_integer_` on error.

### `mm_report_package_version()`
Wraps `utils::packageVersion("mixeff")` in `tryCatch`; returns
`NA_character_` if the package is not installed.

---

## `print.mm_model_report`

Prints a short header, then an overview sub-table with fields
`formula`, `fit_method`, `nobs`, `fit_status`, `inference` (columns
`field`, `value`, `status`), then a bullet list of section names, then
the unavailable count if non-zero.
Does not print any recommendation language (enforced by test).

---

## Known Sections and Compact Column Sets

| Section | Compact columns |
|---------|----------------|
| `overview` | `field`, `value` |
| `model_specification` | `field`, `value`, `status` |
| `data_design` | `group`, `role`, `group_levels`, `min_rows_per_group`, `median_rows_per_group`, `max_rows_per_group`, `status` |
| `random_terms` | `term_id`, `original_fragment`, `group`, `basis`, `covariance`, `theta_parameters`, `design_status`, `english` |
| `random_effects` | `group`, `basis_lhs`, `kind`, `variance`, `std_dev`, `correlation`, `status` |
| `fixed_effects` | `term`, `estimate`, `std_error`, `df`, `statistic`, `statistic_name`, `p_value`, `method`, `status`, `reliability` |
| `fit_statistics` | `field`, `value` |
| `optimizer` | `metric`, `value`, `status` |
| `reproducibility` | `field`, `value` |
| `unavailable` | `section`, `field`, `status`, `reason` |
| `comparison_ledger` | `comparison_id`, `formula`, `comparison_method`, `statistic`, `p_value`, `status`, `reason` |

All compact sets are further trimmed by `mm_drop_empty_report_columns()`.

---

## Rust FFI calls

None directly. All Rust-authored data is accessed through:
- `fit$artifact` (the `CompiledModelArtifact` JSON-backed object)
- `fit$fit_summary` (parsed from `mixedmodels.fit_summary` schema)
- `fit$schema`
- R-level verbs: `inference_table()`, `VarCorr()`, `audit_design()`,
  `random_blocks()`, `changes()`, `optimizer_certificate()`,
  `reproducibility()`, `fit_status()`

---

## Stubs / Deferred / Known Gaps

1. **`random_effects` section — fallback path still documented in unavailable
   ledger.** When `fit$fit_summary$varcorr` is absent (older artifact schemas
   or non-LMM paths), the section falls back to `VarCorr(fit)` and logs a
   `schema_gap` unavailable entry. The entry references
   `planning/reporting_artifact_requirements.md` as source, suggesting a
   planned migration that is not yet complete for all code paths.

2. **`comparison_ledger` on individual fits is always empty.** This is
   by design, but the section exists structurally to make the schema
   consistent; the non-applicable reason is recorded explicitly.

3. **`model_report.mm_glmm` not explicitly registered.** GLMM fits work
   via inheritance but `model_class` in overview emits `class(fit)[[1L]]`
   — the `"LMM"` label is only applied when `inherits(fit, "mm_lmm")`.
   No test covers GLMM-specific reporting.

4. **`section` filter for `reporting_table.mm_fit` builds all sections
   first.** When `allow_many = FALSE` and a single section is requested,
   `model_report(fit, sections = section)` is called which runs all 11
   builders. This is a performance issue only.

5. **`mm_report_fit_statistics` does not guard against NULL slots.**
   `fit$logLik`, `fit$deviance`, etc. are directly indexed without
   `%||%` protection; if a slot is absent the resulting data frame value
   will be `"NULL"` as character. No known case triggers this but it
   is a latent fragility.

6. **The `"details"` and `"notes"` columns referenced in the audit-view
   test expectation** (`test-reporting.R:95`) are expected to appear in
   the `fixed_effects` audit table, but they are not added by
   `mm_report_fixed_effects()` itself — they must come from
   `inference_table()`. Whether they are always present depends on the
   inference layer, not reporting.R. If `inference_table()` ever stops
   producing them, the test will fail at the reporting layer even though
   the bug is in inference.

---

## Test Coverage Summary

`tests/testthat/test-reporting.R` contains 9 test cases:

| Test | What it checks |
|------|---------------|
| assembles required sections | structure of `mm_model_report` |
| no provenance holes | every section has `source`/`status`/`reason`; non-available rows have reasons |
| unavailable ledger is explicit | all 6 required columns present, no blank reasons |
| compact report fields | `overview` compact = `c("field","value")`; required fields present |
| compact and audit views | audit has `source`, `reason`, `details`, `notes`; compact subsets audit |
| fixed-effect rows preserve inference status | `method`, `status`, `reason` match `inference_table()` |
| random-term rows preserve Rust cards | required columns present; `group == "subject"` exists; english non-empty |
| data-design includes grouping counts | `group_levels`, `min_rows_per_group`, `max_rows_per_group` correct |
| random-effect VarCorr payload | primary path (`fit_summary`) taken; no schema_gap unavailable entry |
| no recommendation language | print output lacks forbidden phrases |
| durable comparison ledgers | `compare()` and `drop1()` results accessible via `reporting_table()` |
| saved fits preserve report | `saveRDS`/`readRDS` round-trip yields identical tables |
