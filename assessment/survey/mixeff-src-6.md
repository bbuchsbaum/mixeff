# mixeff source survey: R/marginal.R and R/emmeans.R

**Survey date:** 2026-05-31
**Files:** `R/marginal.R` (641 lines), `R/emmeans.R` (191 lines)
**Family:** mixeff-src-6

---

## 1. R/marginal.R — Native Marginal Surface

### 1.1 Exported generics and S3 methods

All four public entry points are exported as S3 generics with `.mm_lmm` dispatch
only. No `.mm_glmm` dispatch is registered anywhere in the package.

#### `mm_grid()` / `mm_grid.mm_lmm()`

**Contract:** Build a reference grid (Cartesian product of fixed-predictor
levels/values) from a fitted `mm_lmm`.

**Arguments:**

| Arg | Default | Notes |
|-----|---------|-------|
| `fit` | — | `mm_lmm` object |
| `specs` | — | character vector or one-sided formula `~ x` or `~ x \| by` |
| `by` | `NULL` | grouping variable(s); also parsed from `~ x \| by` syntax |
| `at` | `list()` | named list overriding predictor values in the grid |
| `cov.reduce` | `mean` | function reducing undisplayed numeric predictors to a scalar |
| `...` | — | reserved |

**What it returns:** An `mm_grid` S3 object with slots `grid` (data.frame),
`X` (fixed-effect design matrix), `specs`, `by`, `at`, `grid_id`,
`factor_vars`, `numeric_vars`.

**Rust FFI calls:** None directly. Calls `stats::model.matrix()` and
`fixef(fit)` (which delegates to the Rust handle).

**Refusals / NAs:**
- Unknown predictor names in `specs`/`by`/`at` → `mm_abort(class="mm_arg_error")`.
- `at` not a named list → `mm_abort`.
- Design-matrix column mismatch vs. `fixef(fit)` → `mm_abort(class="mm_inference_unavailable")`.
- Unsupported predictor type (not factor/character/logical/numeric) → `mm_abort`.
- Numeric `at[[var]]` with unknown factor levels → `mm_abort`.
- Numeric covariates not explicitly listed in `specs`/`by` and having >10 unique
  values are reduced to a single scalar via `cov.reduce` (default: `mean`);
  those with ≤10 unique values are kept as-is.

**Stubs / deferred:** None. Fully implemented for `mm_lmm`.

---

#### `mm_predictions()` / `mm_predictions.mm_lmm()`

**Contract:** Population-level fixed-effect predictions at each reference-grid
cell, with inference routed through `contrast()`.

**Arguments:** `fit`, `grid` (optional `mm_grid`), `specs`, `by`, `at`,
`method`, `level`, `target`, `scale`, `...`

**`method` values (matched):** `"auto"`, `"satterthwaite"`, `"kenward_roger"`,
`"bootstrap"`, `"asymptotic"`, `"none"`.

**`target`:** Only `"population"` accepted (single-element `match.arg`).

**`scale`:** `"response"` or `"link"` accepted; for Gaussian LMMs these are
identical. The value is stored in the output table but no transformation is
applied.

**What it returns:** `mm_marginal_quantity` S3 object with slots `table`
(data.frame conforming to `mixedmodels.marginal_quantity_table` schema v1.0.0),
`grid`, `L`, `requested_method`.

**Row columns in `table`:** `quantity`, `label`, `estimate`, `rhs`,
`std_error`, `df`, `statistic`, `statistic_name`, `p_value`, `conf_low`,
`conf_high`, `method`, `requested_method`, `status`, `reliability`,
`estimability`, `reason`, `target`, `scale`, `weights`, `comparison`, `by`,
`specs`, `grid_id`, `details`, `notes`.

**Confidence intervals:** Computed from `std_error` + `df` via `stats::qt` (or
`stats::qnorm` when `df` is non-finite), at `level`.

**Rust FFI:** Indirect — via `contrast(fit, L, method=method)`.

---

#### `mm_means()` / `mm_means.mm_lmm()`

**Contract:** Marginal (least-squares) means for cells defined by `specs` (and
optionally `by`), averaging over nuisance fixed-factor levels via a weighted
linear combination of grid rows.

**Additional argument vs. `mm_predictions()`:** `weights = c("equal",
"proportional")`.

- `"equal"`: each reference-grid cell contributes equally (1/n per group).
- `"proportional"`: cells weighted by observed counts of factor levels in
  `fit$model_frame`.

**What it returns:** Same `mm_marginal_quantity` shape as `mm_predictions()`,
but `quantity = "mean"` and `comparison = "identity"`.

**Stubs / deferred:** None. Fully implemented.

---

#### `mm_comparisons()` / `mm_comparisons.mm_lmm()`

**Contract:** Pairwise differences of marginal means, grouped within `by`
levels (all C(n,2) pairs per by-group).

**Additional argument:** `comparison = c("difference", "ratio", "odds_ratio")`.

**`"ratio"` and `"odds_ratio"` are explicitly refused** via
`mm_match_marginal_comparison()` → `mm_abort(class="mm_inference_unavailable")`.
Only `"difference"` is implemented.

**What it returns:** `mm_marginal_quantity` with `quantity = "comparison"` and
`comparison = "difference"`. Labels are `"A - B"` strings.

**Stubs / deferred:**
- `comparison = "ratio"` and `comparison = "odds_ratio"` listed in the argument
  signature but immediately refused — stubs present in the API surface.

---

### 1.2 Print methods

- `print.mm_grid`: prints `grid$grid` without row names.
- `print.mm_marginal_quantity`: prints `x$table` without row names.

Both are exported S3 methods.

---

### 1.3 Internal helpers (not exported)

| Function | Purpose |
|----------|---------|
| `mm_resolve_grid()` | Coerce or build `mm_grid`; validates class |
| `mm_parse_marginal_specs()` | Parse formula or character `specs` + `by`, including `~ x \| by` syntax |
| `mm_fixed_predictor_vars()` | Extract RHS predictor names from fixed formula |
| `mm_fixed_factor_vars()` | Subset to factor/character predictors |
| `mm_fixed_numeric_vars()` | Subset to numeric predictors |
| `mm_reference_grid_vars()` | Determine which predictors go into the grid |
| `mm_grid_values_like()` | Validate and cast `at[[var]]` values |
| `mm_default_grid_values()` | Default grid values per type (levels / unique / mean) |
| `mm_restore_grid_classes()` | Re-attach factor levels after `expand.grid` |
| `mm_fixed_basis()` | Compute fixed-effect design matrix `X` for a grid, respecting training contrasts |
| `mm_grid_labels()` | Human-readable `"var=val, ..."` labels per grid row |
| `mm_group_basis()` | Average design rows into per-group `L` matrix |
| `mm_rows_match()` | Boolean row selector for matching a group |
| `mm_cell_weights()` | Compute equal or proportional cell weights |
| `mm_marginal_rows_from_contrast()` | Assemble `mm_marginal_quantity` table from `contrast()` output |
| `mm_marginal_intervals()` | Compute `conf_low`/`conf_high` from SE + df |
| `mm_new_marginal_quantity()` | Construct and class the return object |
| `mm_validate_marginal_level()` | Validate `level` argument |
| `mm_marginal_group_frame()` | Parse `label` strings back into a group data frame for pairwise logic |
| `mm_pairwise_rows()` | Generate all C(n,2) index pairs within each by-group |
| `mm_match_marginal_comparison()` | Validate `comparison`; refuse non-`"difference"` values |

---

### 1.4 What marginal.R does NOT do

- No `.mm_glmm` dispatch for any of the four generics — `mm_grid`, `mm_means`,
  `mm_comparisons`, `mm_predictions` are LMM-only.
- No `comparison = "ratio"` or `"odds_ratio"` (signature accepts them, body
  refuses them).
- No non-estimability detection: `nbasis` in the basis is not propagated; the
  function always uses `estimability::all.estble` (set in the emmeans bridge,
  not in the native surface).
- No p-value adjustment (no `adjust=` argument, unlike emmeans `contrast()`).
- No `joint_tests()` equivalent.
- No `emtrends()` equivalent (slope estimates over a continuous variable).
- No `eff_size()` equivalent.
- No `plot.mm_marginal_quantity()`.
- No `confint.mm_marginal_quantity()` / `test.mm_marginal_quantity()` wrappers.

---

## 2. R/emmeans.R — emmeans Bridge

### 2.1 Exported functions

#### `recover_data.mm_lmm()`

**Contract:** emmeans plumbing — reconstruct the model-frame data for an
`mm_lmm` so emmeans can build its reference grid.

**Arguments:** `object` (mm_lmm), `data` (optional override), `...`

**What it does:** Extracts the fixed-effect terms, delegates to
`emmeans::recover_data(object$call, trms, na.action, frame=frame)`.

**Requires:** `emmeans` package (soft dependency, checked via
`requireNamespace`).

**Stubs / deferred:** None.

---

#### `emm_basis.mm_lmm()`

**Contract:** emmeans plumbing — supply the linear basis (`X`, `bhat`, `V`,
`dffun`) for an `mm_lmm` to emmeans.

**Arguments:** `object`, `trms`, `xlev`, `grid`, `method` (default `"auto"`),
`...`

**`method` values:** `"auto"`, `"satterthwaite"`, `"kenward_roger"`,
`"asymptotic"`, `"none"` — note: **no `"bootstrap"`** in this bridge (unlike
the native marginal surface).

**What it returns:** Named list with elements:
- `X`: fixed-effect design matrix for the emmeans grid, columns aligned with
  `fixef(object)`.
- `bhat`: `as.numeric(fixef(object))`.
- `nbasis`: `estimability::all.estble` — always full rank (no
  non-estimability detection).
- `V`: from `vcov(object, type="fixed")`, i.e.,
  `mixedmodels.fixed_effect_covariance_matrix` payload when available.
- `dffun`: closure calling `df_for_contrast()` from the mixeff namespace;
  returns `Inf` if non-finite.
- `dfargs`: `list(object=object, method=method)`.
- `misc`: `list(initMesg=...)` carrying covariance status message.

**Requires:** `emmeans` and `estimability` packages.

**Mode/naming note:** emmeans expects `mode="kenward-roger"` (hyphen) in its
own option system; `emm_basis.mm_lmm()` accepts `method="kenward_roger"`
(underscore) and maps it into `df_for_contrast()`. The `mode=` argument from
`emmeans::emm_basis.merMod()` is not forwarded — callers must pass `method=`
explicitly when calling `emm_basis.mm_lmm()` directly; going through
`emmeans::emmeans(fit, ..., method="asymptotic")` passes `method` via `...`.

**Stubs / deferred:**
- `nbasis` is hardcoded to `estimability::all.estble`; rank-deficient models
  will silently produce numeric results rather than surfacing `nonEst` cells.
- No size-limit guards analogous to lme4's `pbkrtest.limit` / `lmerTest.limit`.

---

#### `recover_data.mm_glmm()`

**Contract:** Same shape as `recover_data.mm_lmm()` but for `mm_glmm`.

**Behavior:** Identical pattern — extracts fixed terms, delegates to
`emmeans::recover_data()`.

---

#### `emm_basis.mm_glmm()`

**Contract:** emmeans plumbing for `mm_glmm`; asymptotic z (df = Inf) only.

**Arguments:** `object`, `trms`, `xlev`, `grid`, `...` — no `method` argument.

**What it returns:** Same seven-element list as `emm_basis.mm_lmm()` except:
- `dffun` is `function(k, dfargs) Inf` (hardcoded asymptotic).
- `dfargs` is `list()`.
- `misc` gets family link info via `emmeans::.std.link.labels(fam, misc)` when
  `object$family` has a `$link` field, enabling `type="response"` back-transform.

**Stubs / deferred:**
- `nbasis` is again hardcoded to `estimability::all.estble`.
- Profile-LL CIs for GLMM via emmeans are out of scope (PRD §3).

---

### 2.2 Internal helpers

| Function | Purpose |
|----------|---------|
| `mm_emmeans_vcov()` | Thin wrapper: `stats::vcov(object, type="fixed")` |
| `mm_emmeans_init_messages()` | Constructs `initMesg` string from `mm_status` / `mm_method` attributes on V; distinguishes available vs. unavailable covariance |

---

### 2.3 S3 Registration

`recover_data.mm_lmm` and `emm_basis.mm_lmm` are:
1. Exported from the package NAMESPACE (`export()` directives).
2. Registered with emmeans at load time via `emmeans::.emm_register("mm_lmm", "mixeff")` in `.onLoad()` (`R/zzz.R` line 50).

`recover_data.mm_glmm` and `emm_basis.mm_glmm` are exported but **not
registered** via `.emm_register()`. They are directly exported so callers can
invoke them, but they do not auto-wire into `emmeans::emmeans(fit, ...)` when
`fit` is `mm_glmm`.

---

### 2.4 What emmeans.R does NOT do

- `emm_basis.mm_glmm` is not registered via `.emm_register("mm_glmm", "mixeff")`.
- `emm_basis.mm_lmm` does not implement the `mode="kenward-roger"` (hyphen)
  argument that lme4's `emm_basis.merMod` accepts; uses underscore form only.
- No size-limit guards (`pbkrtest.limit` / `lmerTest.limit` equivalents).
- `nbasis` always `estimability::all.estble`; no non-estimability detection.
- No `bootstrap` method in `emm_basis.mm_lmm`.

---

## 3. Gap Summary

| # | Gap | Severity | Classification |
|---|-----|----------|----------------|
| 1 | `mm_grid`/`mm_means`/`mm_predictions`/`mm_comparisons` have no `.mm_glmm` dispatch | major | in-scope-missing |
| 2 | `comparison="ratio"` and `comparison="odds_ratio"` listed in signature but immediately refused | minor | partial |
| 3 | `emm_basis.mm_glmm` not registered via `.emm_register()`; `emmeans(glmm_fit, ~x)` does not auto-dispatch | major | in-scope-missing |
| 4 | `nbasis` hardcoded to `estimability::all.estble`; rank-deficient models produce numeric garbage instead of `nonEst` cells | major | in-scope-missing |
| 5 | No size-limit guards for Satterthwaite/KR in emmeans bridge (lme4 has `pbkrtest.limit=3000`) | minor | out-of-scope-by-design |
| 6 | `emm_basis.mm_lmm` accepts `method="kenward_roger"` (underscore) but emmeans ecosystem passes `mode="kenward-roger"` (hyphen) via `lmer.df` option | minor | partial |
| 7 | No `plot.mm_marginal_quantity()` | minor | out-of-scope-by-design |
| 8 | No `p.adjust`/`adjust=` on `mm_comparisons()` | minor | out-of-scope-by-design |
| 9 | No `emtrends()` equivalent | minor | out-of-scope-by-design |
| 10 | No `joint_tests()` equivalent (omnibus F through emmeans machinery) | minor | out-of-scope-by-design |
| 11 | No `eff_size()` equivalent | minor | out-of-scope-by-design |
| 12 | Tests cover all four native marginal functions and both GLMM emmeans methods | — | works |
