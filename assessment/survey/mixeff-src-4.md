# mixeff Source Survey: methods-extract.R, methods-summary.R, methods-print.R

Survey date: 2026-05-31
Source files: R/methods-extract.R, R/methods-summary.R, R/methods-print.R
Adjacent files consulted: R/revive.R (vcov, model.matrix, getME, fitted, residuals, is_singular),
  R/predict.R (fitted.mm_lmm, residuals.mm_lmm), NAMESPACE

---

## 1. methods-extract.R

### 1.1 `fixef` / `fixef.mm_lmm` / `fixef.mm_glmm`

**Exported:** yes (`export(fixef)`, `S3method(fixef, mm_lmm)`, `S3method(fixef, mm_glmm)`)

**Contract:** Returns the fixed-effects vector stored at `object$beta`. Named
numeric vector. No recomputation; reads directly from the R object without a
Rust FFI call. Works after `saveRDS()`/`readRDS()` with no live Rust handle.

**What it refuses/NAs:** Nothing; always returns whatever is in `$beta`.

**lme4 gap:** lme4's `fixef()` accepts `add.dropped` (reinstate aliased
coefficients as NA) and `noScale` (undo autoscaling). Neither argument is
implemented here; the signatures accept `...` but ignore all extras.

**Stubs/TODO:** None.

---

### 1.2 `ranef` / `ranef.mm_lmm` / `ranef.mm_glmm`

**Exported:** yes

**Contract (condVar = FALSE):** Returns `object$random_effects` directly — a
`mm_ranef`-classed named list of data frames (one per grouping factor, rows =
levels, columns = RE terms). No FFI call.

**Contract (condVar = TRUE, LMM):** Attempts to compute conditional
variances via `mm_cond_var_postvars(object)`, which calls the Rust FFI
`wrap__mm_lmm_cond_var_json` (schema `mixeff.lmm_cond_var` v1). On success,
each data frame in the returned list receives a `postVar` attribute: a
`p × p × n` numeric array (dimnames = slope-names × slope-names × level-names).
On failure, the array is filled with `NA_real_` and the data frames receive
`mm_unavailable_reason = "random_effect_conditional_variance_unavailable"` and
`mm_cond_var_error = <message>`.

**Contract (condVar = TRUE, GLMM):** Always takes the typed-refusal path.
Returns NA-filled postVar arrays with
`mm_unavailable_reason = "random_effect_conditional_variance_unavailable_for_glmm"`.
No FFI call.

**cond_var FFI detail:** `mm_compute_cond_var_postvars()` re-runs a full Rust
LMM fit via `.Call(wrap__mm_lmm_cond_var_json, ...)`, passing the formula,
REML flag, model-frame columns (translated via `mm_translate_data`), weights,
and control JSON. The result is cached on `fit$lazy_cache` via `.mm_lazy()` so
repeated calls do not re-invoke Rust.

**Block-diagonal merge:** When the same grouping factor appears in multiple RE
terms (e.g., `(1|g) + (0+t|g)`), `mm_merge_block_diag_postvar()` assembles a
block-diagonal postVar treating the two blocks as independent. Off-diagonal
blocks are explicitly set to zero.

**lme4 gaps:**
- `drop` argument (simplify single-column data frames to vectors): not implemented.
- `whichel` argument (select subset of grouping factors): not implemented.
- `postVar` alias (deprecated synonym for `condVar`): not implemented.
- `as.data.frame.ranef.mer` (`grpvar`/`term`/`grp`/`condval`/`condsd`
  long-format): no `as.data.frame.mm_ranef` method exists.
- `dotplot` / `qqmath` visualization methods: not implemented (out-of-scope
  for CRAN).

**Stubs/TODO:** None explicitly, but GLMM `condVar` is a permanent typed
refusal (by design; GLMM conditional variance not yet supported in Rust).

---

### 1.3 `coef.mm_lmm` / `coef.mm_glmm`

**Exported:** yes (`S3method(coef, mm_lmm)`, `S3method(coef, mm_glmm)`)

**Contract:** Returns subject-specific (BLUP) coefficients as an
`mm_coef`-classed named list of data frames. Each column that matches a
fixed-effect name has `ranef + fixef` added. Columns in the RE data frame that
do not appear in `fixef` are returned as pure random modes (no fixed
contribution added).

**Rust FFI:** None. Arithmetic over `ranef(object)` and `fixef(object)`.

**lme4 gaps:**
- `plot.coef.mer` visualization: no `plot.mm_coef` defined.
- The identity `coef(fit)$G[i,j] == fixef(fit)[j] + ranef(fit)$G[i,j]` holds
  for matched columns only; no warning is emitted for unmatched columns.

**Stubs/TODO:** None.

---

### 1.4 `VarCorr` / `VarCorr.mm_lmm` / `VarCorr.mm_glmm`

**Exported:** yes (`export(VarCorr)`, `S3method(VarCorr, mm_lmm)`,
`S3method(VarCorr, mm_glmm)`)

**Contract:** Returns `object$varcorr` — a pre-built `mm_varcorr` object
constructed at fit time by `mm_varcorr_from_result()`. Structure: a list with
`$table` (data frame with columns `group`, `name`, `variance`, `std_dev`,
`correlation`, `boundary`) and `$residual_sd` (scalar numeric).

**Boundary flag:** `mm_varcorr_boundary_flag()` marks a component as on-boundary
when `std_dev <= max(1e-8, 1e-6 * residual_sd)`. This is model state, not a
warning (PRD §9.4).

**Correlation text:** Lower-triangle entries only; formatted as `"+0.xx"` strings
stored in the `correlation` column. First row always `""`.

**lme4 gaps:**
- lme4 returns a `VarCorr.merMod` object: named list of variance-covariance
  matrices with `"stddev"`, `"correlation"`, `"theta"`, `"rho"`, `"profpar"`,
  and `"sc"` attributes. mixeff returns a flat data-frame table instead. Code
  that calls `as.data.frame(VarCorr(fit))` or indexes by group name
  (`VarCorr(fit)$Subject`) will fail.
- `as.data.frame.VarCorr.merMod` columns (`grp`, `var1`, `var2`, `vcov`,
  `sdcor`) are not produced; mixeff produces `group`/`name`/`variance`/`std_dev`
  instead.
- `sigma` argument (scale multiplier): not implemented.
- `print.VarCorr.merMod` `comp` argument (`"Variance"` vs `"Std.Dev."` display):
  not applicable; mixeff prints both columns unconditionally.

**Stubs/TODO:** None.

---

### 1.5 `sigma.mm_lmm` / `sigma.mm_glmm`

**Exported:** yes (`S3method(sigma, mm_lmm)`, `S3method(sigma, mm_glmm)`)

**Contract:** Returns `object$sigma` — scalar numeric. No FFI call.

**Stubs/TODO:** None.

---

### 1.6 `logLik.mm_lmm` / `logLik.mm_glmm`

**Exported:** yes

**Contract:** Returns a `"logLik"` S3 object wrapping `object$logLik` with
`df = object$dof` and `nobs = object$nobs` attributes.

**lme4 gap:** lme4's `logLik.merMod` accepts a `REML` argument that switches
between the REML and ML log-likelihood at evaluation time. mixeff returns
whichever was computed at fit time; no on-the-fly switching.

**Stubs/TODO:** None.

---

### 1.7 `deviance.mm_lmm` / `deviance.mm_glmm`

**Exported:** yes

**Contract:** Returns `object$deviance` — scalar numeric. No FFI call.

**lme4 gap:** lme4's `deviance.merMod` also accepts a `REML` argument.
Not implemented in mixeff.

**Stubs/TODO:** None.

---

### 1.8 `AIC.mm_lmm` / `AIC.mm_glmm`

**Exported:** yes

**Contract:** Returns `-2 * object$logLik + k * object$dof`. Accepts `k`
(default `2`). Deliberately **refuses** multi-model `AIC(fit1, fit2, ...)` with
an `mm_inference_unavailable` condition, directing users to `compare()`.

**Stubs/TODO:** None.

---

### 1.9 `BIC.mm_lmm` / `BIC.mm_glmm`

**Exported:** yes

**Contract:** Returns `object$BIC` (pre-computed at fit time). Refuses
multi-model `BIC(...)` with an `mm_inference_unavailable` condition.

**Note:** Does not recompute from `logLik`/`nobs`/`dof`; uses the stored
scalar. If the fit was run with REML, the BIC is REML-based unless Rust
returns ML-BIC; callers should check `$REML`.

**Stubs/TODO:** None.

---

### 1.10 `nobs.mm_lmm` / `nobs.mm_glmm`

**Exported:** yes. Returns `object$nobs` — integer scalar.

---

### 1.11 `df.residual.mm_lmm` / `df.residual.mm_glmm`

**Exported:** yes. Returns `object$df_residual` — numeric scalar.

---

### 1.12 `formula.mm_lmm` / `formula.mm_glmm`

**Exported:** yes. Returns `x$formula` — the original R formula object stored
on the fit. No `fixed.only` / `random.only` filtering (lme4 supports both).

---

### 1.13 `model.frame.mm_lmm` / `model.frame.mm_glmm`

**Exported:** yes. Returns `formula$model_frame` — the model frame data frame
stored at fit time. No `fixed.only` argument (lme4 supports this).

---

### 1.14 Internal helpers (not exported)

| Helper | Role |
|---|---|
| `mm_attach_ranef_postvar_unavailable(df, reason)` | Attaches NA-filled `p×p×n` postVar array with reason attributes |
| `mm_attach_ranef_postvar(df)` | Back-compat alias for above (old name) |
| `mm_attach_ranef_postvars(ranef_list, postvars)` | Aligns and attaches real postVar arrays from Rust payload |
| `mm_cond_var_postvars(fit)` | Lazy-cached dispatcher: calls `mm_compute_cond_var_postvars` via `.mm_lazy()` |
| `mm_compute_cond_var_postvars(fit)` | Full FFI round-trip to `wrap__mm_lmm_cond_var_json`; parses JSON schema `mixeff.lmm_cond_var` v1 |
| `mm_merge_block_diag_postvar(existing, incoming, group)` | Block-diagonal postVar merge for `(1\|g) + (0+t\|g)` style |
| `mm_ranef_from_terms(terms)` | Parses Rust JSON `terms` list into `mm_ranef` list of data frames; handles same-group merging via `cbind` |
| `mm_varcorr_from_result(varcorr)` | Parses Rust `varcorr` JSON into `mm_varcorr` (table + residual_sd) |
| `mm_varcorr_boundary_flag(std_dev, residual_sd)` | Returns logical vector; threshold `max(1e-8, 1e-6 * scale)` |
| `mm_varcorr_correlation_text(correlations, row_index)` | Formats lower-triangle correlation entries as `"+0.xx"` strings |

---

## 2. methods-summary.R

### 2.1 `summary.mm_lmm`

**Exported:** yes (`S3method(summary, mm_lmm)`)

**Signature:**
```r
summary.mm_lmm(object,
               tests  = c("coefficients", "none"),
               method = c("auto", "satterthwaite", "kenward_roger",
                          "bootstrap", "asymptotic", "none"), ...)
```

**Contract:**

- `tests = "coefficients"` (default): calls `inference_table(object, method = inf_method)` (defined in R/inference.R) to build a per-coefficient inference table. The `method` argument is forwarded; `"bootstrap"` is accepted syntactically but falls through to `"auto"` in the method guard (it is not a recognized `inf_method` for coefficients).
- `tests = "none"`: skips inference; all statistic/p-value columns in the coefficient table are `NA` with `method = "not_computed"`.

**Returns:** S3 object of class `"summary.mm_lmm"` with fields:
`call`, `formula`, `REML`, `coefficients` (data frame), `sigma`, `logLik`,
`AIC`, `BIC`, `nobs`, `df_residual`, `fit_status`, `varcorr`, `tests`,
`inference`, `requested_method`.

**Coefficient table (`mm_summary_coefficients`):**
Columns: `Estimate`, `Std. Error`, `df`, `<stat_col>`, `<p_col>`, `method`.
The statistic and p-value column names are dynamically selected:
`"t value"` / `"Pr(>|t|)"` for Satterthwaite/KR, `"z value"` / `"Pr(>|z|)"`
for asymptotic, `"statistic"` / `"p.value"` when mixed or unknown.

**Verbosity gate:** The inference-rows footer (reliability_reason column
printout) is suppressed by default. Opt-in via `print(summary(fit), verbose = TRUE)`
or `options(mixeff.verbose = 1L)`. Implemented in `mm_summary_verbose()`.

**Stubs/TODO:**
- `method = "bootstrap"` is listed in the signature but is not routed to
  a bootstrap inference path; it falls through the `inf_method` guard to
  `"auto"`. This is a partial/aspirational argument.

---

### 2.2 `print.summary.mm_lmm`

**Exported:** yes (`S3method(print, summary.mm_lmm)`)

**Contract:** Prints:
1. Fit method (REML/ML) and formula.
2. Fit status string.
3. Variance components (`print(x$varcorr)`).
4. Fixed-effects coefficient table.
5. Inference status table (term/method/status/reliability/reliability_reason
   columns) when `tests = "coefficients"` was used.
6. Per-term reason lines when any `$reason` column entry is non-NA/non-empty.
7. Long footer about reliability_reason only when `verbose = TRUE` or
   `options(mixeff.verbose >= 1L)`.

Returns `invisible(x)`.

---

### 2.3 `summary.mm_glmm`

**Exported:** yes (`S3method(summary, mm_glmm)`)

**Signature:**
```r
summary.mm_glmm(object, tests = c("none", "coefficients"), ...)
```

**Note on default:** For GLMM the default is `tests = "none"` (unlike LMM where
`"coefficients"` is the default). Users must explicitly request `tests = "coefficients"`.

**Contract:**

- `tests = "coefficients"`: calls `mm_glmm_wald_z_inference(object)`, which
  calls `stats::vcov(object)` to retrieve the fixed-effect covariance matrix,
  then computes `z = beta / se` and `p = 2 * pnorm(|z|)`. Uses stored
  `$std_errors` when the vcov is unavailable (diagonal fallback). When
  `mm_vcov_status$status` is not `"available"`, all z/p are set to `NA_real_`
  and `method = "not_computed"`.
- `tests = "none"`: no inference computed; all statistic columns NA.

**Returns:** S3 object of class `"summary.mm_glmm"` with fields:
`call`, `formula`, `family`, `method`, `nAGQ`, `coefficients`, `dispersion`,
`logLik`, `AIC`, `BIC`, `nobs`, `df_residual`, `fit_status`, `varcorr`,
`tests`, `inference`, `vcov_status`.

**Stubs/TODO:**
- Profile-LL CIs for GLMM are out-of-scope (PRD §3).
- Satterthwaite/Kenward-Roger df for GLMM not available; only asymptotic
  Wald-z is implemented.
- `denominator_df` column is always `NA_real_` for GLMM (no df method).

---

### 2.4 `print.summary.mm_glmm`

**Exported:** yes (`S3method(print, summary.mm_glmm)`)

**Contract:** Prints:
1. Header line, formula, family/link, method/nAGQ, fit status.
2. Variance components.
3. Fixed-effects coefficient table.
4. Wald-z reliability note when `vcov_status$reliability` is non-NA,
   non-empty, and not `"available"` (i.e., degraded reliability).

Returns `invisible(x)`.

---

### 2.5 Internal helpers

| Helper | Role |
|---|---|
| `mm_glmm_wald_z_inference(object)` | Builds Wald-z `mm_inference_table` from stored vcov; typed-refusal path returns NA z/p with status attributes |
| `mm_summary_verbose(...)` | Returns TRUE when `verbose=TRUE` arg or `options(mixeff.verbose >= 1)` |
| `mm_summary_coefficients(object, inference)` | Merges inference table rows with stored beta/std_errors; builds display data frame |
| `mm_summary_statistic_column(statistic_name)` | Picks display column name (`"t value"`, `"z value"`, etc.) from unique statistic types |
| `mm_summary_p_value_column(statistic_name)` | Picks display column name (`"Pr(>|t|)"`, `"Pr(>|z|)"`, etc.) |

---

## 3. methods-print.R

### 3.1 `print.mm_lmm`

**Exported:** yes (`S3method(print, mm_lmm)`)

**Contract:** Prints a compact model summary:
1. Fit method (REML/ML), formula, fit status.
2. Optimizer line: name, iteration count, objective value — from
   `x$artifact$optimizer_certificate`. Suppressed silently when the
   certificate is absent (returns `""`).
3. Artifact schema name/version and crate version (from `x$schema` or
   `mm_object_schema(x$artifact)`).
4. nobs, sigma, logLik.
5. Fixed effects (6 significant digits via `signif(fixef(x), 6)`).
6. Singular-fit summary when `is_singular(fit)` is TRUE (see §3.5 below).
7. Audit verb reminder: `audit()`, `diagnostics()`, `inference_table()`,
   `model_report()`.

Returns `invisible(x)`.

---

### 3.2 `print.mm_glmm`

**Exported:** yes (`S3method(print, mm_glmm)`)

**Contract:** Similar to `print.mm_lmm` but:
1. No REML/ML distinction; header is `"Generalized linear mixed model fit"`.
2. Prints family/link and method/nAGQ lines.
3. Prints `dispersion` instead of `sigma`.
4. No singular-fit summary block (not called for GLMM; `is_singular` only
   dispatches on `mm_lmm`).
5. Audit verbs line omits `inference_table()`.

Returns `invisible(x)`.

---

### 3.3 `print.mm_varcorr`

**Exported:** yes (`S3method(print, mm_varcorr)`)

**Contract:** Prints variance components table with `variance` and `std_dev`
columns (6 sig figs), a `note` column showing `"[boundary]"` for boundary
components, and a footer `"[boundary]: ..."` when any component is flagged.
Prints `"Residual std. dev.:"` when `residual_sd` is finite.
Prints `"none"` when the table is empty.

**lme4 gap:** lme4's `print.VarCorr.merMod` has a `comp` argument to choose
between `"Variance"` and `"Std.Dev."` display. mixeff always prints both.
lme4 formats correlations inline in the table; mixeff stores them as formatted
text in the `correlation` column.

---

### 3.4 `print.mm_ranef`

**Exported:** yes (`S3method(print, mm_ranef)`)

**Contract:** Iterates over the named list and prints each group's data frame
with a `$groupname` header. Prints `"none"` for an empty list.

---

### 3.5 `print.mm_coef`

**Exported:** yes (`S3method(print, mm_coef)`)

**Contract:** Iterates over the named list with `"Conditional coefficients:"`
header and `$groupname` subheaders. Prints `"none"` for an empty list.

---

### 3.6 Internal helpers

| Helper | Role |
|---|---|
| `mm_singular_render_lines(fit)` | Returns character vector of singular-fit message lines; empty when not singular. Sources `fit$artifact$effective_covariance` for rank info. Points to `changes(fit)` and `random_options(spec, group=...)`. |
| `mm_singular_first_group(fit)` | Finds first grouping factor from `artifact$design_audit$random_term_cards`; falls back to `names(fit$random_effects)[[1]]`; then `"<group>"`. |
| `mm_print_optimizer_line(x)` | Formats optimizer certificate (name/iterations/objective) from `x$artifact$optimizer_certificate`; returns `""` when absent. |

---

## 4. Dependency map: FFI calls in this surface

| Function | Rust FFI | Schema |
|---|---|---|
| `ranef.mm_lmm(condVar=TRUE)` | `.Call(wrap__mm_lmm_cond_var_json, ...)` | `mixeff.lmm_cond_var` v1 |
| `ranef.mm_glmm(condVar=TRUE)` | none — typed refusal | n/a |
| All other extractors | none | n/a |

All other extractors (fixef, coef, VarCorr, sigma, logLik, deviance, AIC, BIC,
nobs, df.residual, formula, model.frame, vcov, fitted, residuals,
model.matrix) read directly from the R object; they do not require a live Rust
handle and work after `saveRDS()`/`readRDS()`.

---

## 5. Gap summary

### 5.1 `fixef`: `add.dropped` / `noScale` not implemented

lme4 supports `add.dropped = TRUE` to reinstate aliased columns as NA. mixeff
does not expose rank-deficient column handling at the R level (the Rust
compiler drops them silently and records it in the design audit). Severity:
minor for most users; a gap for users with rank-deficient fixed designs.

### 5.2 `ranef`: `drop`, `whichel`, `postVar` arguments not implemented

`drop = TRUE` (simplify single-column RE frames to vectors) and `whichel`
(subset grouping factors) are standard lme4 usage patterns (e.g., used in
`as.data.frame.ranef.mer`). `postVar` is a deprecated lme4 synonym for
`condVar`. Not implemented. Minor severity; cosmetic/convenience.

### 5.3 `ranef`: no `as.data.frame.mm_ranef`

lme4's `as.data.frame(ranef(fit))` produces a long-format data frame
(`grpvar`, `term`, `grp`, `condval`, `condsd`) used by `ggplot2` caterpillar
plots and `broom::tidy`. No `as.data.frame.mm_ranef` exists. Moderate severity
for users building caterpillar plots without lattice.

### 5.4 GLMM `condVar` permanently NA

`ranef.mm_glmm(condVar = TRUE)` always returns NA-filled postVar arrays.
The GLMM conditional variance path is not yet implemented in the Rust crate.
Classified upstream-blocked.

### 5.5 `VarCorr` returns `mm_varcorr` not `VarCorr.merMod`

The return type is a flat table, not the lme4-style named list of matrices.
Downstream code that calls `as.data.frame(VarCorr(fit))` expecting `grp`/`var1`/`var2`/`vcov`/`sdcor` columns, or that indexes `VarCorr(fit)$Subject` to get a variance-covariance matrix, will fail. The `sigma` argument is also absent. Severity: major for any lme4-compatible downstream (broom, performance, see).

### 5.6 `summary.mm_glmm` default is `tests = "none"`

lme4/lmerTest default is to show coefficient tests in summary. mixeff's GLMM
summary defaults to no tests. By design (Wald-z inference requires the vcov
to be available, which is not guaranteed for revived objects). Minor severity
as a documentation/discoverability issue.

### 5.7 `summary.mm_lmm` method `"bootstrap"` is listed but not routed

`method = "bootstrap"` is accepted in the signature but falls through the
`inf_method` guard (it is not in the recognized set for `inference_table`) and
resolves to `"auto"`. This is aspirational API surface, not a real bootstrap
implementation. Classified partial/in-scope-missing.

### 5.8 `logLik`/`deviance`: no `REML` switching argument

lme4 accepts `REML = NULL` on `logLik`/`deviance` to retrieve the ML or REML
criterion on demand. mixeff returns whichever was computed at fit time. Minor.

### 5.9 `ngrps` not implemented

lme4's `ngrps(object)` returns the number of levels per grouping factor —
commonly used in model summaries and sample-size checks. No `ngrps` method
exists for `mm_lmm` or `mm_glmm`. Minor severity.

### 5.10 `getME` subset is smaller than lme4's

lme4 exposes 38+ named components including `"b"`, `"u"`, `"sigma"`,
`"devcomp"`, `"L"`, and `"ALL"`. mixeff's `getME.mm_lmm` supports:
`X`, `Z`, `Zt`, `Lambda`, `Lambdat`, `theta`, `beta`/`fixef`, `y`, `mu`,
`flist`, `cnms`. Any unsupported name raises `mm_arg_error`. Severity: minor
for most workflows; moderate if downstream packages call `getME(fit, "ALL")`.

### 5.11 `residuals.mm_lmm` supports only `type = "response"`

lme4 supports `"response"`, `"pearson"`, `"deviance"`, `"working"` and a
`scaled` argument. mixeff only implements `"response"` (stored
`object$residuals`). Severity: minor for standard use; moderate for diagnostics
workflows relying on Pearson or scaled residuals.

### 5.12 `formula`/`model.frame`: no `fixed.only`/`random.only` arguments

lme4's formula and model.frame methods accept `fixed.only` and `random.only`
flags. Not implemented. Minor.

### 5.13 No `REMLcrit` or `isREML`

lme4 exports `REMLcrit(object)` (REML criterion at optimum) and
`isREML(object)` (logical). mixeff exposes `$REML` directly on the object but
does not register these generics. Minor.
