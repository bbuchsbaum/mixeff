# Survey: mixeff-src-3
## Files: R/inference-options.R, R/random-options.R, R/mm-control.R

Survey date: 2026-05-31
Surveyor: automated subagent (claude-sonnet-4-6)

---

## 1. R/mm-control.R

### Exported symbols

#### `mm_control(verbose = 0L)`
- **Class returned:** `mm_control` (a plain list with one element: `verbose`)
- **Contract:** Validates that `verbose` is a single non-NA numeric, coerces to integer. Aborts with class `mm_arg_error` on bad input.
- **Semantics of `verbose`:**
  - `>= 0` — prints the pre-fit `explain_model()` output once before optimization (handled in `lmm()` R-side, not Rust-side).
  - `-1` — suppresses the automatic model explanation printout.
- **Rust FFI:** The entire `mm_control` list is serialized to JSON via `jsonlite::toJSON(unclass(control), ...)` and passed as `control_json` to Rust FFI entry points (`mm_fit_lmm`, `mm_fit_glmm`, `mm_loglik_lmm`, etc.). However, on the Rust side the parsed value is bound to `_control` (leading underscore = intentionally unused). The control JSON is validated/parsed for schema correctness but **no Rust-side behavior is currently gated on its content.** The `optimizer_certificate` in the artifact is set by the Rust optimizer internally, not by anything in the control struct.
- **Callers:** `lmm()`, `glmm()`, `simulate.mm_lmm()`, `predict.mm_lmm()`, `logLik.mm_lmm()`, plus internal refit paths in `inference.R` and `compare.R` (all pass `mm_control(verbose = -1)` to suppress output during automated refits).

#### `mm_validate_control(control)` (internal, not exported)
- Accepts `NULL` or missing → returns `mm_control()` defaults.
- Accepts a list → extracts `verbose` (defaulting to `0L`) and returns a validated `mm_control`.
- Rejects non-list → `mm_arg_error`.

### Gaps relative to lme4::lmerControl
lme4 exposes: `optimizer`, `optCtrl` (optimizer-specific options), `calc.derivs`, `use.last.params`, `sparseX`, `check.*` family of validators, `checkControl`, `nAGQ0initStep`. **None of these are currently surfaced by `mm_control`.** This is by design: the Rust optimizer is self-contained and not user-configurable from the R side yet.

---

## 2. R/inference-options.R

### Exported symbols

#### `inference_options(fit, term = NULL, nsim = 1000L, ...)` — generic
- **Dispatches** via `UseMethod`.
- `term`: reserved for future per-term refinement; validated against `mm_fixed_effect_terms(fit)` but otherwise currently unused in row construction (table is fit-level, not term-level).
- `...`: reserved for future methods, currently absorbed silently.

#### `inference_options.mm_lmm(fit, term = NULL, nsim = 1000L, ...)`
- **Returns:** S3 object of class `mm_inference_options`:
  - `$table` — data frame, one row per method, columns: `method`, `expected_status`, `expected_reliability_reason`, `r_verb`, `approx_cost`, `notes`, `current`, `display_status`, `display_reason`, `what_to_do_next`.
  - `$fit_status`, `$is_reml`, `$n_groups_max`, `$term`.
- **Rust FFI:** None. This is a pure R-side prediction from fit metadata (`fit$fit_status`, `fit$REML`, `ranef(fit)`, `inference_table(fit)`). Does not refit or call any Rust entry point.
- **Seven rows produced** (one per method):

| method | status logic |
|---|---|
| `asymptotic_wald_z` | always `"available"` |
| `satterthwaite` | `"not_assessed"` if boundary/reduced-rank, else `"available"` |
| `kenward_roger` | `"not_assessed"` if boundary/reduced-rank, else `"available"` |
| `bootstrap` | always `"available"` |
| `bootstrap_lrt` | `"not_assessed"` if REML fit, else `"available"` |
| `cluster_bootstrap` | always `"not_assessed"` (p-values not certified) |
| `profile_ci` | `"not_assessed"` if boundary or REML; else `"available"` |

- **`approx_cost` column:** order-of-magnitude estimate only; uses a heuristic of 0.02 s/replicate for n<200 and 0.05 s/replicate otherwise. Explicitly documented as "out of scope" for calibration.
- **`current` column:** marks the method that `inference_table(fit)$table$method[[1]]` resolved to (i.e., what `auto` picked).
- **`not_yet_wired` status code** is defined in `mm_inference_options_display_status()` and `mm_inference_options_next_step()` but is never set by any of the seven row-building helpers — it exists as a forward slot for methods wired upstream but not yet bridged to R.
- **No recommendation column.** The table explicitly tests `expect_false("recommended" %in% tolower(names(opt$table)))`.

#### `print.mm_inference_options(x, ...)`
- Prints: fit_status header + subset of table columns (`method`, `display_status`, `display_reason`, `what_to_do_next`, `approx_cost`, `current`) using human-readable display columns (not raw enum strings).
- Raw enum columns are accessible via `<obj>$table`.

### Internal helpers (all unexported)

| Helper | Purpose |
|---|---|
| `mm_inference_options_n_groups_max(fit)` | Max grouping-level count across random effects; used as cluster-bootstrap viability proxy |
| `mm_inference_options_row_wald/satterthwaite/kenward_roger/bootstrap/bootstrap_lrt/cluster_bootstrap/profile_ci` | Build one named list (one row) per method |
| `mm_inference_options_bootstrap_lrt_reliability_reason(nsim)` | `"bootstrap_monte_carlo_replicates"` iff nsim >= 999, else `"bootstrap_insufficient_replicates"` |
| `mm_inference_options_add_display(tab)` | Adds `display_status`, `display_reason`, `what_to_do_next` via `mapply` over enum columns |
| `mm_inference_options_display_status(status)` | Enum → human-readable status string |
| `mm_inference_options_display_reason(reason, status, method)` | Enum → human-readable reason; falls back to `gsub("_", " ", reason)` for unknown keys |
| `mm_inference_options_next_step(method, status, reason, r_verb)` | Produces the `what_to_do_next` guidance string |
| `mm_inference_options_format_cost(fit, nsim, factor)` | Heuristic seconds/minutes estimate |

---

## 3. R/random-options.R

### Exported symbols

#### `random_options(spec, group, slope = NULL)` — generic (no UseMethod, direct function)
- **Input:** `spec` must be an `mm_spec` (from `compile_model()`) or `mm_fit`; validated by `mm_assert_compiled_spec()`.
- **`group`:** bare or string; validated against random-effect cards for that group via `audit_design(spec)`.
- **`slope`:** bare or string; if missing/NULL, resolved by `mm_default_slope()` (priority: current random slopes → scope-note fixed effects → first non-intercept fixed effect). Errors if no slope variable is resolvable or if the resolved variable is not in `spec$model_frame`.
- **Rust FFI:** Indirect — calls `compile_model()` and `audit_design()` per candidate, which in turn call the Rust compilation/audit FFI. Does not call fit-side Rust entry points.
- **Five candidates always evaluated:**
  1. `(1 | group)` — intercept-only ("punt")
  2. `(0 + slope | group)` — slope only
  3. `(1 | group) + (0 + slope | group)` — split uncorrelated
  4. `(1 + slope || group)` — double-bar synonym
  5. `(1 + slope | group)` — full correlated

- **Returns:** S3 object of class `mm_random_options`:
  - `$group`, `$slope` — labels
  - `$options` — data frame: `formula`, `varying_coefficients`, `covariance_family`, `theta_parameters`, `design_status`, `plain_meaning`, `note`, `current`
  - `$cards` — named list (by candidate key) of upstream `RandomTermCard` records
  - `$constraints` — named list of cross-card constraint lists per candidate
  - `$reports` — named list of full audit reports per candidate
- **`current` column:** set by exact-match of candidate fragment against `mm_current_random_fragment(current_cards)` (reconstructs original fragment from card `original_fragment` fields).
- **No recommendation column.** Tested explicitly (`expect_false("recommended" %in% ...)`).
- **Errors thrown:** `mm_arg_error` if spec invalid; `mm_schema_error` if no cards for group or no slope resolvable; `mm_data_error` if slope variable not in model frame.

#### `print.mm_random_options(x, ...)`
- Prints: group header, current model block (formula + plain_meaning), then all nearby options with formula, varying_coefficients, covariance_family, theta_parameters, design_status, plain_meaning.

### Internal helpers (all unexported)

| Helper | Purpose |
|---|---|
| `mm_random_option_candidates(group, slope)` | Builds the fixed 5-candidate list |
| `mm_random_option_row(candidate, cards, constraints, current)` | Produces one data.frame row from card data |
| `mm_option_plain_meaning(cards, constraints)` | Assembles English text from block `english` fields + constraint reasons |
| `mm_option_covariance_family(key, covariances)` | Overrides covariance label for split/double-bar candidates to "diagonal via separate blocks" |
| `mm_current_random_fragment(cards)` | Reconstructs original fragment string from card `original_fragment` fields |
| `mm_default_slope(spec, audit, group)` | Three-tier slope resolution (card slopes → scope notes → fixed effects) |
| `mm_candidate_formula(spec, random_fragment)` | Rebuilds a full formula from spec semantic model + new random fragment |
| `mm_expr_label(expr, allow_null)` | Converts bare name or string to a label string |
| `mm_assert_compiled_spec(spec)` | Guards mm_spec/mm_fit inheritance |
| `mm_spec_model_frame(spec)` | Extracts model frame or errors |
| `mm_cards_for_group(cards, group)` | Filters card list by group label |

---

## 4. Cross-cutting observations

### Rust FFI contract
- `mm_control` is serialized to JSON and threaded through every Rust entry point, but the Rust side currently ignores all fields (parsed to `_control`, discarded). The only optimizer metadata flowing back is the `optimizer_certificate` in the fit artifact, which is Rust-internal.
- `inference_options` and `random_options` make no direct Rust FFI calls; they operate on metadata from prior calls to `compile_model()`, `audit_design()`, and `inference_table()`.

### Refusal / NA patterns
- `inference_options`: never returns NA for method rows; instead uses structured `expected_status = "not_assessed"` with a stable `expected_reliability_reason` enum.
- `random_options`: errors loudly (never silently produces an empty or NA table) if preconditions fail.
- `mm_control`: errors on invalid `verbose`; never silently coerces to a default.

### Forward slots / stubs
- `inference-options.R` line 276: `not_yet_wired = "available upstream; R bridge pending"` — a display-status enum value defined and handled (line 335 produces a guidance string) but never assigned by any current row builder. This is a deliberate forward slot for methods that exist in the Rust engine but lack an R bridge.
- `inference_options(fit, term = ...)`: `term` is validated but unused in row construction — per-term refinement is reserved for a future phase.
- `inference_options(..., ...)`: `...` is reserved for future methods dispatch.
- `mm_control(verbose)`: `verbose` is the only current field. Optimizer selection, tolerance controls, and convergence checking (all offered by lme4's `lmerControl`) are absent and unsurfaced.

### Test coverage
- `test-inference-options.R`: 9 tests covering boundary/non-boundary status, method enumeration, bootstrap/bootstrap_lrt/cluster_bootstrap refusals, REML gating, multi-df F-test, and term validation.
- `test-inference-options-display.R`: 4 tests covering display column population, no-enum-leak invariant, print output, and route-table refusal reason alignment.
- `test-random-options.R`: 5 tests covering nearby map rendering, current-marker accuracy for split/double-bar, card JSON round-trip, forbidden advice phrases, and input validation. (Also covers `compare_covariance()` incidentally.)
- No dedicated test file for `mm_control` in isolation; it is exercised indirectly throughout all fitting tests via `mm_control(verbose = -1)`.
