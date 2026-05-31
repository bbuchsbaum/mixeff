# Test Survey: tests-6 — Reporting, Phase 4, Explain-Model, Diagnostic Formatters

Survey date: 2026-05-31  
Files surveyed:
- `tests/testthat/test-reporting.R`
- `tests/testthat/test-phase4.R`
- `tests/testthat/test-explain-model.R`
- `tests/testthat/test-diagnostic-formatters.R`

Companion snapshot file: `tests/testthat/test-explain-model-snapshots.R` (read for context).

---

## test-reporting.R

### What is tested

| Test | Assertions |
|---|---|
| `model_report()` assembles required sections | S3 class `mm_model_report`; exact section names; presence of `metadata`, `sections`, `unavailable`, `provenance`; non-zero metadata rows; provenance row count == section count |
| No provenance holes | Every section data-frame has `source`, `status`, `reason` columns; no blank `source`; non-available rows must have a non-blank `reason` |
| Unavailable ledger is explicit | Required columns present; no blank `reason` or `source`; all `action_taken == "reported"` |
| `reporting_table()` overview | Correct column names `c("field","value")`; expected `field` values present; idempotent when called on `mm_model_report` vs `mm_fit` |
| Compact vs audit views | Compact hides `source`/`details`; audit exposes `source`, `reason`, `details`, `notes`; compact columns are a subset of audit |
| Fixed-effect rows preserve Rust inference status | `method`, `status`, `reason` match `inference_table()`; `source` present; at least one `"available"` row |
| Random-term rows preserve Rust cards | Required column names; `group == "subject"` present; non-empty `english`; `source` in expected set |
| Data-design grouping unit counts | `subject` row present; `group_levels == 10`; `min/max_rows_per_group == 5`; `status` column present |
| Random-effect report consumes VarCorr payload | Schema name/version check; `kind == "variance"` row; `group == "Residual"` row; all `status == "available"`; all `source == "mixedmodels.fit_summary.varcorr"`; stable payload not in unavailable ledger |
| No recommendation language | `print(model_report(fit))` does not match five forbidden patterns |
| Comparison ledger via `reporting_table(cmp, ...)` | Compact/audit columns; `source` absent from compact; `all_sections` matches; `dropped_ledger` row count equals `dropped$ledger`; `reference_formula` identity; requesting non-comparison section raises `mm_schema_error` |
| Saved fits preserve report sections | `reporting_table()` of revived fit is identical to original for `fixed_effects` and `random_effects` |

**Fixture:** `mk_reporting_fit()` — 10 subjects × 5 obs, `y ~ x + (1 | subject)`, seed 41. Single-slope random intercept model. `compare(reduced, full)` also exercised.

**Tolerances:** Identity comparisons only (`expect_identical`, `expect_equal`). No numeric tolerances.

**Skips:** None.

### What is NOT tested (gaps)

- `model_report()` on an **mm_glmm** fit — all tests use `mm_lmm`. The GLMM report path (especially `fixed_effects` / `random_effects` for binomial/Poisson) is untested.
- `model_report()` with an explicit **subset of sections** (`sections = c("overview", "fixed_effects")`) rather than `"all"`.
- `reporting_table()` with **`section = "all"`** on an `mm_fit` — only individual sections are extracted.
- `reporting_table()` on an **`mm_random_effect_test`** object — the method exists (`S3method(reporting_table, mm_random_effect_test)`) but is not exercised anywhere in these files.
- **`model_specification`** section values are never inspected; changes rows are not verified when the compiler canonicalizes a formula.
- **`fit_statistics`** section contents (`logLik`, `AIC`, `BIC`, `sigma` values) are not numerically verified.
- **`optimizer`** section content is not checked beyond its presence in `section_names`.
- **`reproducibility`** section field values are not verified.
- `reporting_table()` called on a **dropped object (`mm_drop1`)** — `dropped_ledger` is created and `nrow`/`reference_formula` are checked, but compact vs. audit column content is not compared.
- **Print output** of `model_report` is only tested for forbidden phrases; the actual rendered content (section listing, nobs value, formula string) has no assertion.
- **Revive path** only checks `fixed_effects` and `random_effects`; `random_terms`, `data_design`, and `fit_statistics` sections are not verified post-revive.
- **Invalid `section` names** passed to `reporting_table()` — error class `mm_schema_error` is only demonstrated for the comparison-object case.
- `mm_report_group_max_rows` falling back to `NA_integer_` when the group factor is not in `model_frame` is not exercised.

---

## test-phase4.R

### What is tested

| Test | Assertions |
|---|---|
| `glmm()` cbpp binomial smoke (profiled PIRLS) | S3 classes; `family`/`link` fields; `method`; `nobs`; finite `fixef`, `theta`, `logLik`, `AIC`; `summary`; finite `vcov` diagonal; `model.matrix` row count; `summary(tests="coefficients")` column names and `vcov_status` |
| `glmm()` family/link surface matches contract | All 6 supported family/link combos fit without error; `logLik` and `fixef` finite |
| `glmm()` refuses off-contract family/link with stable reason code | `mm_inference_unavailable` raised; `reason_code == "unsupported_glmm_family_link"`; `family`/`link` echoed; `supported` contains expected columns |
| `mm_glmm` revive | S3 class; handle dead; lazy cache is env; `fixef`/`ranef`/`VarCorr` equal to original; `predict` raises `mm_inference_unavailable` |
| `glmm()` validates family and reports unavailable joint backend | `mm_fit_error` for `joint_laplace`; message contains `"estimation_method_unavailable"`; metadata echoes family/link/method; gaussian raises `mm_inference_unavailable` with correct reason code |
| `simulate.mm_lmm` reproducible | S3 class `data.frame`; `dim` correct; two seeds equal; `attr(,"mm_method")`; `refit()` returns `mm_lmm` with same `nobs` and `formula` |
| `compare()` and `anova()` refit REML→ML | `mm_model_comparison`; table/ledger row counts; required column names; `fit_method == "ML"` throughout; at least one `refit`; `REML == FALSE`; finite LRT; `lrt_available`; ledger/table identity checks; error for `refit_for_comparison = "error"` |
| `compare()` non-nested validity rows | Finite AIC/BIC; `lrt_available == FALSE`; NA LRT/p-value; `reason_code == "non_nested_models_lrt_invalid"`; `comparison_class == "non_nested_fixed_effects"`; error on `method = "lrt"` |
| `drop1()` preserves random effects | S3 class `mm_drop1`; required table/ledger columns; `x` and `z` present in `dropped`; random term present in all formulas; `reference_formula` identity |
| Parametric bootstrap tiny nsim | `mm_parametric_bootstrap`; `status` in `c("available","not_assessed")`; finite `p_value` or NA with non-NA reason; `simulated` is numeric; `successful_replicates` present; `compare()` bootstrap slot; method string `"parametric_bootstrap_lrt"` |
| Bootstrap refuses non-nested / mismatched data | `mm_arg_error` with `reason_code == "bootstrap_lrt_requires_same_observations"` and `"bootstrap_lrt_requires_nested_models"` |
| `compare(method="bootstrap")` validates nsim | Error class `mm_arg_error` for `nsim = NA` and `nsim = -1` |
| Bootstrap refuses REML fits | `mm_inference_unavailable`; `reason_code == "bootstrap_lrt_requires_ml"` |
| Manifest advertises capabilities | `cap$simulate`, `cap$inference`, `cap$fit_glmm` all `TRUE` |

**Fixture:** `mk_phase4_fit()` — 8 subjects × 5 obs, `y ~ x + z + (1 | subject)` (slope variant) or `y ~ x + (1 | subject)` (reduced). Also `mk_cbpp_glmm_fit()` using lme4's cbpp expanded to binary. `mk_glmm_contract_data()` for each family.

**Tolerances:** Identity/equality checks; no numeric tolerances.

**Skips:** `skip_if_not_installed("lme4")` in `mk_cbpp_glmm_fit()`.

### What is NOT tested (gaps)

- **Poisson and Gamma GLMM** `summary(tests="coefficients")` — Wald-z table is only verified for the binomial cbpp fit; GLMM coefficient table for other families is untested.
- **`glmm()` with random slopes** — all GLMM tests use `(1 | group)` intercept-only; `(1 + x | group)` GLMM is not exercised.
- **AGQ method** (`method = "agq"`) — only `"pirls_profiled"` and `"joint_laplace"` (to test rejection) appear; AGQ convergence and coefficient recovery are not tested.
- **`drop1()` on a GLMM** — `drop1.mm_lmm` is exported but there is no `mm_glmm` method; the absence of the method is not tested.
- **`compare()` with three or more models** — only two-model comparisons are exercised.
- **`compare()` REML fits with `refit_for_comparison = FALSE`** — the warning/error path for this edge case is not covered.
- **`simulate()` with `nsim = 1`** — only `nsim = 2` is tested; edge case of a single simulation is not verified.
- **`simulate()` method attribute** is only checked for the LMM Gaussian path; a "non-Gaussian" or GLMM simulation path is not tested.
- **`refit()` with a formula that differs from the stored formula** — only same-formula refit is tested.
- **`parametric_bootstrap()` with large `nsim`** and `seed` reproducibility — only `nsim = 2` is used; no check that two calls with the same seed produce identical results.
- **`bootstrap_control()` arguments** are not passed to `compare(method="bootstrap")`; `bootstrap_control` constructor is not exercised in these tests.
- **`mm_glmm` `predict()`** on a live (non-revived) fit — the test only checks that predict raises on a revived (dead-handle) fit; the live case is not checked.
- **`anova()` on a single model** or error handling for equal-formula models in `compare()`.
- **Manifest `mm_formula_manifest()`** — only `$capabilities` is checked; formula manifest rows, capability flags beyond the three tested, and unsupported patterns are not inspected.

---

## test-explain-model.R

### What is tested

| Test | Assertions |
|---|---|
| `explain_model()` returns printable `mm_explanation` | S3 class; `text` type/length; `cards` and `report` are lists; printed output contains `"Random effects explanation"` |
| Eight Phase 1.C syntax patterns | Each of the eight patterns produces text matching expected named-form strings, `scope_note:` tag, English sentences, `theta parameters: 3`, `double-bar syntax`, `separate random-effect blocks`, `syntax_expansion:` tag, `r0 <-> r1` constraint, nested canonical expansion, interaction `s:i`, crossed two-group output |
| Structural refusal renders as "Possible repairs" | `structural_refusal` in text; `"Possible repairs, not applied automatically:"` header; `"\`between\` does not vary within \`g\`"` phrase |
| R9 forbidden advice phrases | Five forbidden patterns do not appear in explain output for four formula patterns |
| `explain_model()` refuses non-compiled inputs | `mm_schema_error` raised for `list()` input |

**Fixture:** `mk_explain_design()` — `expand.grid(s=6, i=4, b=2)` plus `a`, `t`, `y`. Also inline `df` for the between-subject refusal case.

**Tolerances:** `expect_match` / `expect_false(grepl(...))` only.

**Skips:** None.

### What is NOT tested (gaps)

- **`explain_model()` on a fitted `mm_lmm`** — all tests call `explain_model(compile_model(...))` (i.e., an `mm_spec`). The `mm_fit` dispatch path is not exercised in this file (test-explain-model-snapshots.R also only uses `compile_model`).
- **`explain_model()` on a fitted `mm_glmm`** — GLMM explain output is entirely untested.
- **Fit notes section** (`"Fit notes:"`) — `mm_fit_note_lines()` output is not tested via `explain_model()`; it only appears through the `diagnostics()` function indirectly.
- **Singularity section** (`"Fitted covariance state:"`) — `mm_singularity_lines()` is only exercised if the Rust artifact returns a rank-deficient covariance; no test creates this condition.
- **Multiple random-effect terms with the same group** (`(1 | s) + (0 + t | s)`) — the `"s has 2 separate random-effect blocks."` summary line is checked (split pattern), but the constraint cross-card line rendering (`r0 <-> r1 (Intercept <-> t): separate random-effect blocks...`) is only checked via `expect_match` without verifying it does not duplicate.
- **`print.mm_explanation`** trailing newline behaviour (`if (!grepl("\n$", x$text)) cat("\n")`) is not tested.
- **`explain_model()` when `audit$random_term_cards` is empty** — the `"none"` fallback path is not triggered by any test.
- **`random_options()` output** is not tested here (it is tested in the snapshot file, but not in the main unit test file).
- **`compare_covariance()` output** — only snapshot-tested; no assertion on column names, row count, or absence of "recommended" column in the unit-test file.
- **Formula with `I()` or `poly()` fixed effects** passed to `compile_model()` — outside full v1 scope but partial-support cases are not covered.

---

## test-diagnostic-formatters.R

### What is tested

| Test | Assertions |
|---|---|
| `mm_diagnostic_code_registry` covers every Rust `DiagnosticCode` | `setdiff(snake, registered) == character()` (no missing); `setdiff(registered, snake) == character()` (no stale) |
| `mm_response_diagnostic_reason_registry` covers every `ResponseDiagnosticReason` | Same bidirectional coverage check |
| Every registered code is bound to a valid bucket | All `bucket` values in `c("design_note","repair","fit_note","raw_only")`; `raw_only` entries all have non-empty `rationale` |
| Formattable buckets have non-empty advice renderers | `mm_design_note_lines`, `mm_fit_note_lines`, `mm_repair_lines` each produce ≥ 1 line for a fixed list of hand-crafted diagnostics |
| `mm_repair_lines` falls back to message when `suggested_actions` is empty | Output length 1; matches `"not_identifiable: MLE does not exist"` |
| `mm_diagnostics_table` warns once and attaches attribute on unknown code | Warning text matches `"unrecognized DiagnosticCode"`; attribute `mm_unrecognized_diagnostic_code` set to `"fictional_future_code"`; no second warning in same session; `mm_diagnostic_bucket("fictional_future_code")` is `NA` |
| Registry/bucket helper alignment | `mm_diagnostic_bucket(code)` matches `reg[[code]]$bucket` for every registered code |

**Fixture:** Hand-crafted diagnostic lists inline; paths to vendored Rust source resolved relative to test directory.

**Tolerances:** Identity checks only.

**Skips:** `testthat::skip()` inside helpers if vendored Rust source file or enum not found (graceful skip, not `skip_on_cran`).

### What is NOT tested (gaps)

- **`diagnostics.mm_compiled()` public API** — the formatter registry tests call internal helpers directly; no test checks that `diagnostics(compile_model(...))` produces a correctly structured `mm_diagnostics` object with a populated `$table`.
- **`diagnostics()` `severity` and `stage` filter arguments** — the filter paths in `diagnostics.mm_compiled()` are not exercised.
- **`print.mm_diagnostics`** for a non-empty table — only the empty-diagnostics path would be exercised by other tests; the branch that prints the message block (`if ("message" %in% names(x$table))`) is not directly tested.
- **`fit_status()` on an `mm_fit`** and on an `mm_compiled` — `fit_status` is exported and used in `mm_report_overview()` but has no dedicated test; the fallback chain (`$fit_status`, `$artifact$optimizer_certificate$status`, `"not_assessed"`) is not tested.
- **Multiple unknown codes** in a single `mm_diagnostics_table()` call — the warn-once test uses one unknown code; two simultaneous unknowns are not tested.
- **Warn-once state reset between test sessions** — the state env `mm_unknown_diag_state$seen` is reset manually in the test, but there is no test verifying that a new R session starts with an empty `seen` set.
- **`mm_design_note_lines` / `mm_repair_lines` / `mm_fit_note_lines` with an empty diagnostics list** — the `character()` early-return paths are not directly tested (though they are implied to work by tests that produce no output).
- **`mm_repair_lines` deduplication** — when multiple `repair`-bucket diagnostics share the same code+action pair, the deduplication logic is not exercised.
- **`mm_diagnostic_code_registry` `raw_only` forward-compat** — no test verifies that codes in the `raw_only` bucket do not appear in any of the three advice-line renderer outputs.
- **`ResponseDiagnosticReason` R surface** — all five registered variants have `bucket = "raw_only"` with no R-side formatter; there is no test that the batch-engine path that would emit these codes does not accidentally surface them through `explain_model()` or `print()`.

---

## Cross-cutting gaps

- **`model_report()` on a GLMM** — neither `test-reporting.R` nor `test-phase4.R` exercises `model_report(glmm_fit)`. The GLMM-specific sections (`fixed_effects` via Wald-z, `random_effects` for non-Gaussian) are untested in the reporting contract.
- **`reporting_table(mm_random_effect_test, ...)` compact vs audit** — the S3 method exists but is not tested.
- **`diagnostics()` filter arguments** — `severity` and `stage` filtering are untested across all four files.
- **`fit_status()` fallback chain** — untested in isolation.
- **`explain_model()` on an `mm_fit`** — only `mm_spec` inputs are tested.
