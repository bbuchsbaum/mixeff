# mixeff-src-1 — Survey of `R/glmm.R` and the `mm_glmm` surface

Generated: 2026-05-31  
Source files surveyed: `R/glmm.R` (primary), plus all files that define or dispatch
on `mm_glmm`: `R/methods-extract.R`, `R/methods-print.R`, `R/methods-summary.R`,
`R/emmeans.R`, `R/predict.R`, `R/inference.R`, `R/revive.R`, `R/extendr-wrappers.R`,
`NAMESPACE`, `tests/testthat/test-glmm.R`.

---

## 1. Exported API (`NAMESPACE` registrations for `mm_glmm`)

| Symbol | Kind | Source file |
|---|---|---|
| `glmm()` | exported function | `R/glmm.R` |
| `print.mm_glmm` | S3 method | `R/methods-print.R` |
| `summary.mm_glmm` | S3 method | `R/methods-summary.R` |
| `print.summary.mm_glmm` | S3 method | `R/methods-summary.R` |
| `fixef.mm_glmm` | S3 method (alias of `fixef.mm_lmm`) | `R/methods-extract.R` |
| `ranef.mm_glmm` | S3 method | `R/methods-extract.R` |
| `coef.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `VarCorr.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `sigma.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `logLik.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `deviance.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `AIC.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `BIC.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `nobs.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `df.residual.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `formula.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `model.frame.mm_glmm` | S3 method (alias) | `R/methods-extract.R` |
| `model.matrix.mm_glmm` | S3 method (alias) | `R/revive.R` |
| `vcov.mm_glmm` | S3 method (alias of `vcov.mm_lmm`) | `R/revive.R` |
| `fitted.mm_glmm` | S3 method (alias) | `R/predict.R` |
| `residuals.mm_glmm` | S3 method | `R/predict.R` |
| `predict.mm_glmm` | S3 method — typed refusal | `R/predict.R` |
| `recover_data.mm_glmm` | S3 method (emmeans bridge) | `R/emmeans.R` |
| `emm_basis.mm_glmm` | S3 method (emmeans bridge) | `R/emmeans.R` |
| `mm_lincomb.mm_glmm` | S3 method | `R/inference.R` |

**Not registered for `mm_glmm`** (only `mm_lmm` dispatches exist):
`contrast`, `test_effect`, `confint`, `df_for_contrast`, `inference_table`,
`estimability`, `test_random_effect`.

---

## 2. `glmm()` — fit function

### Signature
```r
glmm(formula, data, family,
     random = NULL, weights = NULL, subset = NULL,
     na.action = na.omit, contrasts = NULL,
     method = c("pirls_profiled", "joint_laplace"),
     nAGQ = 1L,
     inference = c("auto", "none", "asymptotic", "bootstrap"),
     control = mm_control(), ...)
```

### Contract and flow
1. `match.arg(method)` — `"joint_laplace"` is declared in the choice vector but
   **refused at fit time** via the Rust engine (the `nlopt` feature is not compiled
   for CRAN; the error crosses the bridge as a typed `mm_fit_error`).
2. `mm_glmm_family_info(family)` validates that `family/link` is in the certified
   set: `binomial/{logit,probit,cloglog}`, `poisson/{log,sqrt}`, `Gamma/log`.
   Unsupported family/link raises `mm_inference_unavailable` (reason code
   `unsupported_glmm_family_link`) before any Rust call.
3. `random`, `weights`, `subset`, custom `na.action`, `contrasts` are reserved and
   raise `mm_fit_error` immediately if non-null/non-default.
4. `mm_glmm_validate_nagq(nAGQ, method)` validates nAGQ ≥ 1;
   `joint_laplace + nAGQ > 1` raises `mm_arg_error`.
5. `compile_model(formula, data)` — shared LMM compile path.
6. `explain_model(spec)` is printed if `control$verbose >= 0`.
7. FFI call: `.Call(wrap__mm_fit_glmm_json, formula_string, family, link, method,
   nAGQ, column_order, numeric_columns, categorical_values, categorical_levels,
   control_json)` — returns a JSON string.
8. Bridge errors are routed through `mm_abort_from_bridge()`.
9. JSON is parsed by `mm_json_parse_glmm_fit()`, which validates schema header
   `mixeff.glmm_fit_result v1`.
10. The returned `mm_glmm` list is assembled from parsed fields.

### What is stored on the returned object
`call`, `formula`, `family` (list with `$family` and `$link`), `method`,
`nAGQ`, `inference_request`, `control`, `vars`, `model_frame`, `weights=NULL`,
`artifact`, `fit`, `fit_summary`, `schema`, `rust_handle=NULL`,
`lazy_cache` (env), `beta` (named numeric), `theta`, `sigma`, `dispersion`,
`logLik`, `deviance`, `AIC`, `BIC`, `nobs`, `dof`, `df_residual`,
`fit_status`, `std_errors`, `fixed_effect_vcov`, `fixed_fitted=NULL`,
`fitted`, `residuals`, `random_effects` (mm_ranef), `varcorr` (mm_varcorr).

### Explicit stubs / deferred fields
- `fixed_fitted = NULL` — population-level fitted values are not computed for GLMM
  (contrast: LMM sets this from `fit_result$fixed_fitted`).
- `rust_handle = NULL` — no persistent Rust object; all subsequent operations
  re-run the fit from the stored formula + data.
- `inference` argument is stored as `inference_request` but is **not used** to
  trigger inference during fit; the `auto`/`asymptotic`/`bootstrap` choices have
  no effect on the returned object.

---

## 3. Family/link validation helpers (unexported)

| Function | Purpose |
|---|---|
| `mm_glmm_family_info(family)` | Accepts R family object or constructor; maps `Gamma` → `"gamma"` for Rust |
| `mm_glmm_supported_family_links()` | Returns the certified list |
| `mm_glmm_supported_family_link_table()` | Formats it as a data.frame for error messages |
| `mm_abort_glmm_unsupported_family_link()` | Raises `mm_inference_unavailable` with structured payload |
| `mm_glmm_validate_nagq()` | Validates nAGQ; enforces `joint_laplace → nAGQ ≤ 1` |

---

## 4. JSON parsing helpers (unexported)

`mm_json_parse_glmm_fit(json)` — validates schema header
`mixeff.glmm_fit_result v1`, returns parsed list. Raises `mm_schema_error` on
bad JSON or unknown schema.

---

## 5. FFI primitive

`mm_fit_glmm_json(formula, family, link, method, n_agq, column_order,
numeric_columns, categorical_values, categorical_levels, control_json)`
→ `.Call(wrap__mm_fit_glmm_json, ...)` (generated binding in
`R/extendr-wrappers.R`).

`method = "pirls_profiled"` maps to upstream `fast = true`. `"joint_laplace"`
requires the `nlopt` feature (not enabled); any attempt reaches Rust and comes
back as a bridge error.

---

## 6. S3 methods provided for `mm_glmm`

### Extractors (aliases of `mm_lmm` implementations)
`fixef`, `coef`, `VarCorr`, `sigma`, `logLik`, `deviance`, `AIC`, `BIC`,
`nobs`, `df.residual`, `formula`, `model.frame`, `model.matrix`, `vcov`.
All delegate directly to the `mm_lmm` implementation via `= fixef.mm_lmm` etc.

### `ranef.mm_glmm`
Returns `object$random_effects`. When `condVar = TRUE`, immediately returns NA
postVar arrays with `reason = "random_effect_conditional_variance_unavailable_for_glmm"`
(no Rust bridge call attempted — typed refusal, not a code stub).

### `vcov.mm_glmm` (= `vcov.mm_lmm`)
`type = "fixed"` (default): returns `object$fixed_effect_vcov` if present,
otherwise calls `mm_fixed_effect_vcov_from_payload()` which reads the artifact's
`fixed_effect_covariance_matrix` field. Falls back to a diagonal matrix from
stored `std_errors` with `mm_status = "unavailable"`.
`type = "theta"`: always returns NA matrix with reason `theta_covariance_unavailable`.

### `predict.mm_glmm`
**Hard refusal**: raises `mm_inference_unavailable` with message
`"GLMM prediction is not certified by the current Rust contract."` regardless of
arguments. `newdata`, `re.form`, etc. are not inspected.

### `fitted.mm_glmm` / `residuals.mm_glmm`
Both are aliases of the LMM implementations and return the in-sample values stored
on the fit object.

### `print.mm_glmm`
Prints family/link, method, nAGQ, fit status, optimizer line, artifact schema,
nobs/dispersion/logLik, and fixed effects. Does **not** print singularity lines
(those are LMM-only). Audit verbs listed: `audit()`, `diagnostics()`,
`model_report()` — note **no** `inference_table()` in the GLMM print footer.

### `summary.mm_glmm`
`tests = c("none", "coefficients")`. Default is `"none"` (note: reversed from
LMM where default is `"coefficients"`). When `tests = "coefficients"`,
calls `mm_glmm_wald_z_inference()` to build a Wald-z table from
`vcov(object)`. Returns `summary.mm_glmm` object; `print.summary.mm_glmm`
displays it with optional vcov reliability warning.

### `inference_table` for `mm_glmm`
**Not registered** — `NAMESPACE` shows only `S3method(inference_table, mm_lmm)`.
Calling `inference_table(fit)` on an `mm_glmm` will dispatch to S3 default
(likely error or wrong method).

### `contrast` / `test_effect` / `confint` / `df_for_contrast` / `estimability` for `mm_glmm`
**Not registered** — all these generic S3 dispatches are registered only for
`mm_lmm`. Calling them on an `mm_glmm` will not dispatch through the
dedicated LMM Rust inference bridge (which re-fits an LMM internally).

### `mm_lincomb.mm_glmm`
Exported and registered. Accepts `method = "asymptotic"` or `"auto"` only.
Builds the Wald-z from `vcov(fit)` R-side (no Rust bridge). Returns a
single-row data.frame with `estimate`, `std_error`, `statistic` (z),
`df = NA`, `p_value`, `lower`, `upper`, `method = "asymptotic"`.

### emmeans bridge (`recover_data.mm_glmm` / `emm_basis.mm_glmm`)
Both exported. `emm_basis.mm_glmm` uses `dffun = function(k, dfargs) Inf`
(asymptotic z, no df). Passes the family/link to `emmeans::.std.link.labels`
so `type = "response"` applies the inverse link. Fixed-effect covariance
comes from `vcov(object)`. Note: the `family` field on the fit is
`list(family, link)` not a full R family object with `linkfun`/`linkinv`;
the emmeans bridge checks `!is.null(fam$link)` before calling `.std.link.labels`.

---

## 7. Deferred/stub/refused paths — summary

| Path | Classification | Evidence |
|---|---|---|
| `method = "joint_laplace"` | out-of-scope-by-design | Doc: "requires optional upstream `nlopt` backend … intentionally does not enable" |
| `predict.mm_glmm(newdata = ...)` | in-scope-missing (deferred) | Hard `mm_inference_unavailable` refusal; no Rust path exists for GLMM `predict_new` |
| `ranef(condVar = TRUE)` for GLMM | out-of-scope-by-design | Typed refusal, reason `random_effect_conditional_variance_unavailable_for_glmm` |
| `inference_table.mm_glmm` | in-scope-missing | Not registered; no S3 dispatch |
| `contrast.mm_glmm` | in-scope-missing | Not registered; inference bridge is LMM-only |
| `test_effect.mm_glmm` | in-scope-missing | Not registered |
| `confint.mm_glmm` | in-scope-missing | Not registered; only `confint.mm_lmm` exists |
| `df_for_contrast.mm_glmm` | in-scope-missing | Not registered |
| `estimability.mm_glmm` | in-scope-missing | Not registered |
| `test_random_effect.mm_glmm` | out-of-scope-by-design (LMM-only boundary-LRT) | Only `mm_lmm` dispatch registered |
| `profile confint for GLMM` | out-of-scope-by-design (PRD §3) | Explicitly deferred v2 |
| `fixed_fitted` (population fitted) | in-scope-missing | `fixed_fitted = NULL` hard-coded at line 145 |
| `inference` argument effect | partial | Stored as `inference_request` but not acted on during or after fit |
| AGQ with `nAGQ > 1` | partial | Accepted as argument, passed to Rust, but no R-side AGQ-specific processing; outcome depends entirely on Rust response |
| `weights` / `subset` / `contrasts` | out-of-scope-by-design | Reserved; raise `mm_fit_error` when supplied |
| `sigma` for GLMM | partial | Returns `object$sigma` which holds the dispersion parameter; for `binomial`/`poisson` families this may be 1.0 (not NA) — no guard against misleading output |
| `vcov(type = "theta")` | out-of-scope-by-design (typed refusal) | Always NA with reason `theta_covariance_unavailable` |
| `emm_basis.mm_glmm` family object | partial | `object$family` is `list(family=, link=)` not a real R family; `.std.link.labels` gets `list` not a proper family — may silently fail if emmeans expects `$linkinv` |

---

## 8. Test coverage

`tests/testthat/test-glmm.R` (139 lines):
- Two parity cases: `cbpp_binomial_logit_profiled_pirls`, `grouseticks_poisson_log_profiled_pirls`
- Both declared `expected_mismatch` (not strict lme4 parity)
- Tests: `nobs`, `family`, `link`, `method`, `fixef`, `theta`, `logLik`, `deviance`
- **Not tested**: `summary`, `print`, `vcov`, `ranef`, `predict`, `emmeans`,
  `mm_lincomb`, `sigma`, `AIC`, `BIC`, `VarCorr`, `fitted`, `residuals`,
  `Gamma` family, `probit`/`cloglog` binomial, `sqrt` Poisson link

---

## 9. Cross-cutting notes

- **`rust_handle = NULL`**: every inference operation that reaches Rust
  (for LMM) re-runs the full fit from formula + data. For GLMM this is
  moot because no post-fit Rust inference is wired; all inference is R-side Wald-z.
- **No `REML` field**: `mm_glmm` has no `$REML` slot; the inference bridge
  functions (`mm_rust_contrast_table`, etc.) call `isTRUE(fit$REML)` which
  evaluates to `FALSE` on a GLMM. This is safe for the functions that are
  not dispatched for GLMM, but would be a latent bug if any shared path
  were ever called on a GLMM.
- **`is_singular` not defined for GLMM**: `is_singular.mm_lmm` inspects
  `fit$artifact$effective_covariance`; there is no `is_singular.mm_glmm`
  registered, so calling `is_singular(glmm_fit)` would fall through to
  `is_singular.default` and raise `mm_arg_error`.
