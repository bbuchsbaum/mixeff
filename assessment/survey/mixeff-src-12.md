# mixeff source survey — batch 12

Files surveyed: `R/manifest.R`, `R/roles.R`, `R/data-translate.R`,
`R/simulate.R`, `R/zzz.R`, `R/mixeff-package.R`

Survey date: 2026-05-31

---

## R/manifest.R

### `mm_formula_manifest()` — exported

**Purpose.** Machine-readable capability document for the current build.
Answers "what does this version of mixeff know how to do?" Every `mm_fit`
object is supposed to snapshot this at construction time so `audit()` can
answer the question for cold/revived fits without a live Rust handle.

**Signature.** `mm_formula_manifest()` — no arguments.

**FFI call.** `.Call(wrap__mm_formula_manifest)` — single Rust entry point;
the entire response is produced in Rust and returned as a named R list.

**Returns.** Named list with five keys:

| key | type | content |
|-----|------|---------|
| `mixeff_rust_version` | `character(1)` | version of the bundled extendr crate |
| `crate_version` | `character(1)` | version of the upstream `mixedmodels` crate |
| `schema_versions` | named list | one entry per artifact schema (`formula`, `mixedmodels.compiled_model_artifact`, `mixedmodels.model_audit_report`, `mixedmodels.random_term_card`, `mixedmodels.fixed_effect_inference_table`, `mixedmodels.model_comparison_table`, `mixedmodels.fit_summary`, `mixedmodels.marginal_quantity_table`) |
| `formula_features` | named list | `operators`, `intercept_forms`, `random_term_forms`, `transformations` |
| `capabilities` | named list of logicals | 25 flags (see below) |

**Capability flags documented in the function's `@return` section:**
`parse_formula`, `compile_model`, `audit_design`, `explain_model`,
`random_options`, `compare_covariance`, `fit_lmm`, `fit_glmm`, `audit`,
`changes`, `diagnostics`, `fit_status`, `parameterization`, `roles`,
`as_json`, `simulate`, `inference`, `model_comparison_table`,
`fit_summary_payload`, `marginal_quantity_table`, `marginal_quantities`.

The test suite (`test-manifest.R`) expects 25 flags including four not listed
in the docstring: `fixed_effect_inference_table`, `satterthwaite`,
`kenward_roger_explicit`, `bootstrap_fixed_effect_payload`. This is a minor
docstring discrepancy — the Rust side emits all 25.

**Refuses / NA behaviour.** None — all validation is inside Rust.

**TODO/FIXME markers.** None in this file.

**Test coverage.** Full: `tests/testthat/test-manifest.R` (7 tests),
`tests/testthat/test-schema-versioning.R`, and a capability-flag check in
`test-phase4.R`.

---

## R/roles.R

### `roles(...)` — exported

**Purpose.** Phase 1.F verb. Two operating modes:

1. **Declaration mode** — named string arguments, e.g.
   `roles(subject = "sampled_unit", x = "fixed_condition")`. Builds an
   `mm_roles` object whose `source` is `"declared"`.
2. **Observation mode** — single unnamed `mm_compiled` object (spec or fit).
   Extracts grouping factors and fixed terms from the compiled artifact's
   `semantic_model` and classifies each as `sampled_unit`,
   `observed_within_group`, `observed_between_group`, or
   `observed_fixed_effect`. Sets `source = "observed_from_data"`.

**Returns.** `mm_roles` S3 object: list with `$table` (data.frame with columns
`variable`, `role`, `origin`, `group`, `evidence`) and `$source` (character).

**No FFI calls.** Pure R; reads from `mm_compiled_artifact(x)` (defined in
`diagnostics.R`) when in observation mode.

**Zero-argument call.** Returns an empty `mm_roles` with a zero-row table.

**Refuses.** Mixed unnamed/named arguments raise `mm_arg_error`. Non-string
role values raise `mm_arg_error`.

### `print.mm_roles(x, ...)` — exported S3 method

Prints `"Design roles:\n"` header then calls `print()` on `x$table`, or
`"  none declared\n"` if the table is empty. Returns `invisible(x)`.

### Internal helpers

- `mm_roles_from_compiled(x)` — observation-mode logic; iterates
  `random_terms` for grouping factors, then `fixed_terms` for fixed effects,
  calling `mm_observed_fixed_role()` for each.
- `mm_observed_fixed_role(variable, groups, frame)` — classifies one fixed
  variable against all grouping factors by splitting the data and counting
  unique values per group. Returns `within_group`, `between_group`, or plain
  `fixed_effect`.
- `mm_roles_empty_table()` — constructs the canonical empty data.frame
  schema.

**TODO/FIXME markers.** None.

**Test coverage.** `test-audit-verbs.R` (lines 74–92): declared mode and
observation-mode fallback via `mk_audit_fit()`. `test-no-advice.R` line 38:
`print(roles(fit))` does not emit advice text.

---

## R/data-translate.R

### `mm_translate_data(data)` — internal (`@noRd`)

**Purpose.** Converts an R `data.frame` into the three-list wire format
(`numeric_columns`, `categorical_values`, `categorical_levels`) consumed by
the Rust FFI at every compilation and fitting entry point.

**Called from.** `compile.R:108`, `fit-lmm.R:56`, `glmm.R:68`,
`predict.R:191+195`, `inference.R:1076`, `methods-extract.R:149`.

**Type mapping:**

| R type | FFI classification | Notes |
|--------|-------------------|-------|
| `numeric` / `integer` / `logical` | `numeric` | logicals coerced via `as.numeric()` (0/1) |
| `factor` | `categorical` | levels from `levels(col)` — canonical factor order |
| `character` | `categorical` | levels from `unique(col)` — first-appearance order |
| `Date`, `POSIXct`, list, other | error | raises `mm_data_error` |

**Returns.** Named list with four entries:
- `column_order` — `character` vector of original column names
- `numeric_columns` — named list of `numeric` vectors
- `categorical_values` — named list of `character` vectors (observed values)
- `categorical_levels` — named list of `character` vectors (canonical levels)

**Validation contract.** Refuses (raises `mm_data_error`) if:
- `data` is not a `data.frame`
- `data` has zero columns
- column names are NULL, contain empty strings, or are duplicated
- any column has an unsupported type

**NA behaviour.** Does not check for NAs — that is delegated to
`mm_check_no_na()` which is called before `mm_translate_data()` at the
`compile.R` entry point.

### `mm_check_no_na(data, vars, .call)` — internal (`@noRd`)

**Purpose.** "No silent surgery" enforcement. Enumerates all design variables
with any NA values and raises one `mm_data_error` listing all offenders rather
than failing on the first.

**Called from.** `compile.R:105` only (before `mm_translate_data`).

**Returns.** `invisible(TRUE)` if no NAs found.

**Refuses.** Raises `mm_data_error` with columns `columns`, `na_counts`, and
a message instructing the user to call `na.omit(data)` explicitly.

**TODO/FIXME markers.** None. No deferred paths.

**Test coverage.** Tested indirectly through compile/fit pipelines. No
dedicated unit test file for `data-translate.R` was found.

---

## R/simulate.R

### `refit(object, newresp, ...)` — exported generic

**Purpose.** S3 generic for refitting a model with a new response vector.
Mirrors `lme4::refit()`.

**Returns.** Dispatches to `refit.mm_lmm()`.

### `refit.mm_lmm(object, newresp, ...)` — exported S3 method

**Mechanism.** Copies `object$model_frame`, replaces the response column
(via `mm_response_name()`), and calls `lmm()` with the stored formula,
`REML` flag, and weights. Attaches `$refit` metadata to the result recording
`source = "refit"` and `original_fit_status`.

**Validates.** `newresp` must be numeric, length `nobs(object)`, no NAs.
Raises `mm_arg_error` otherwise.

**No FFI calls.** Delegates entirely to `lmm()`, which makes its own FFI
calls.

**TODO/FIXME.** None.

### `simulate.mm_lmm(object, nsim, seed, re.form, ...)` — exported S3 method

**Purpose.** Parametric Gaussian simulation from stored fixed effects,
random-effect covariances, and residual sigma. Mirrors `lme4::simulate()`.

**Arguments:**
- `nsim` — positive integer (default 1)
- `seed` — optional integer; if provided, sets `set.seed()` and restores the
  original `.Random.seed` on exit
- `re.form` — `NULL` (conditional; new random effects drawn per simulation)
  or `NA` / `~0` (population-level only; RE contribution omitted)
- `...` — reserved, ignored

**Returns.** `data.frame` with `nsim` columns named `sim_1`, ..., `sim_<nsim>`.
Row names match `object$model_frame`. Attributes:
- `attr(..., "seed")` — the seed argument (may be `NULL`)
- `attr(..., "mm_method")` — always `"r_side_gaussian_parametric"`

**Refuses.** Any `re.form` other than `NULL`/`NA`/`~0` raises
`mm_inference_unavailable` with message "re.form requests beyond NULL and NA
are not available for simulation."

**No GLMM support.** Only `simulate.mm_lmm` is defined; there is no
`simulate.mm_glmm`. Simulating from a GLMM fit would fall through to the
default `stats::simulate()` method, which will not understand the object —
effectively missing.

**Covariance reconstruction path (`mm_random_term_covariance`).** Reads
`artifact$covariance_parameter_traces` from the stored JSON artifact, matching
entries by `term_id`. Two passes: first fills diagonal (standard deviations →
variances), then fills off-diagonal (correlations → covariances). Fallback:
if any diagonal entry is still zero or negative after the trace pass, looks up
in `fit$varcorr$table`. If the fallback also fails (row not found), that
diagonal stays at 0 and the draw proceeds with a rank-deficient matrix.

**Degenerate Sigma handling (`mm_rmvnorm`).** Attempts `chol(Sigma)`; if
Cholesky fails (non-PD), falls back to spectral decomposition with
`pmax(eigenvalues, 0)` — positive semi-definite projection. No warning is
emitted on fallback.

**Basis reconstruction.** Helpers `mm_basis_label()`, `mm_basis_values()`,
`mm_group_factor()`, `mm_random_term_group_label()` are all defined in
`revive.R`. They reconstruct random-effect design columns from the stored
model frame and artifact JSON without re-fitting. If a basis column cannot be
found in the model frame, `mm_basis_values()` raises `mm_data_error`.

**TODO/FIXME markers.** None in the file.

**Test coverage.** `test-phase4.R` lines 189–203: reproducibility (same seed
→ same draws), dimension check, `mm_method` attribute, and that `refit()`
returns a valid `mm_lmm` of the same size.

---

## R/zzz.R

### `.onLoad(libname, pkgname)` — package hook

Calls `mm_register_external_s3()` then returns `invisible()`.

### `.onUnload(libpath)` — package hook

Calls `library.dynam.unload("mixeff", libpath)` to unload the compiled Rust
shared library.

### `` `%||%`(x, y) `` — package-internal null-coalescing operator

Defined here so every `R/` file has access without importing from rlang.
`if (is.null(x)) y else x`. Not exported.

### `mm_register_external_s3()` — internal

Registers S3 methods for lme4 and emmeans generics using `setHook` on
`packageEvent` plus an immediate registration path if either namespace is
already loaded at `onLoad` time. The two conditional paths cover:

1. **lme4 not yet loaded** — hook fires when lme4 loads.
2. **lme4 already loaded** — registers immediately.

Same logic for emmeans.

### `mm_register_lme4_s3()` — internal

Registers the following methods into the lme4 namespace:

| generic | class |
|---------|-------|
| `fixef` | `mm_lmm`, `mm_glmm` |
| `ranef` | `mm_lmm`, `mm_glmm` |
| `VarCorr` | `mm_lmm`, `mm_glmm` |
| `getME` | `mm_lmm` |
| `refit` | `mm_lmm` |

### `mm_register_emmeans_s3()` — internal

Calls `emmeans::.emm_register("mm_lmm", "mixeff")`. This registers
`recover_data.mm_lmm` and `emm_basis.mm_lmm` (defined in `R/emmeans.R`) via
emmeans' internal registration mechanism. Only `mm_lmm` is registered; there
is no `mm_glmm` registration here (though `recover_data.mm_glmm` and
`emm_basis.mm_glmm` exist in `emmeans.R` and would need a separate
`.emm_register` call to be discoverable).

**TODO/FIXME markers.** None.

---

## R/mixeff-package.R

### `"_PACKAGE"` docstring — package-level Roxygen sentinel

Standard `"_PACKAGE"` declaration. Provides the package-level `@description`
for the generated `.Rd` file and declares `@importFrom` entries for base
stats generics used throughout:

`AIC`, `BIC`, `coef`, `deviance`, `df.residual`, `fitted`, `formula`,
`logLik`, `model.frame`, `nobs`, `predict`, `residuals`, `setNames`, `sigma`,
`update`.

No functions defined. No exports. No FFI.

---

## Cross-cutting observations

### Wiring correctness

- `mm_translate_data` is called at every FFI boundary (compile, fit, predict,
  inference, ranef). The type contract is consistent and documented.
- `mm_check_no_na` is only called at `compile.R`; the `predict.R` and
  `inference.R` paths that call `mm_translate_data` on `newdata` do **not**
  call `mm_check_no_na` first. If `newdata` contains NAs the Rust side will
  receive `NaN` numeric values or `"NA"` strings in categoricals. This is a
  silent partial path.
- The `%||%` operator in `zzz.R` is the only source; it is not re-exported
  and rlang's `%||%` is not used directly by these files.

### Simulation — partial paths

- Only `mm_lmm` has a `simulate` method. GLMM simulation is absent
  (`simulate.mm_glmm` does not exist).
- The `mm_rmvnorm` degenerate-Sigma fallback (spectral projection) proceeds
  silently — no warning is emitted when Cholesky fails.
- `re.form` support is binary: `NULL` or `NA`/`~0`. Partial RE formulas
  (subset of grouping factors) raise `mm_inference_unavailable` with a clear
  message.

### emmeans registration gap

`mm_register_emmeans_s3()` calls `.emm_register("mm_lmm", "mixeff")` but
does not register `"mm_glmm"`. The `recover_data.mm_glmm` and
`emm_basis.mm_glmm` methods in `emmeans.R` are therefore not discoverable
via the standard emmeans dispatch pathway without a separate
`.emm_register("mm_glmm", "mixeff")` call.

### No deferred/TODO markers

None of the six files contain TODO, FIXME, stub, "not yet", "deferred", or
"NOT IMPLEMENTED" comments.
