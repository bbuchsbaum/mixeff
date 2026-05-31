# mixeff source survey — predict.R / simulate.R
**Family:** mixeff-src-5
**Date:** 2026-05-31
**Files surveyed:** `R/predict.R`, `R/simulate.R`
**lme4 reference:** `assessment/survey/lme4-prediction.md`

---

## 1. Exported API inventory

### 1.1 `predict.mm_lmm` (S3 — `stats::predict`)

**File:** `R/predict.R` lines 36–97

**Signature:**
```r
predict.mm_lmm(object, newdata = NULL, re.form = NULL,
               allow.new.levels = FALSE,
               type = c("response", "link"),
               se.fit = FALSE,
               interval = c("none", "confidence", "prediction"),
               level = 0.95, ...)
```

**Returns:** Named numeric vector (length = nrow prediction target). When `se.fit = TRUE`, returns a list `list(fit = <numeric>, se.fit = <NA_real_ vector>)` with `mm_unavailable_reason = "prediction_se_unavailable_phase_2"` on both the inner `se.fit` vector and the outer list.

**Rust FFI used:**
- Conditional newdata path: `.Call(wrap__mm_lmm_predict_new_json, ...)` — dispatches through `mm_lmm_predict_new_json` (extendr wrapper for `mm_lmm_predict_new_json` in `lib.rs` line 1713). Rebuilds the model from training data + formula + control on the Rust side; returns JSON with schema `mixeff.lmm_predict_new` version 1.
- Population/fixed-only path (`re.form = NA` or `~0`): **no Rust FFI**; pure R matrix multiply `mm_new %*% beta` using `stats::model.matrix` against training contrasts/xlevels.
- In-sample path (`newdata = NULL`): **no Rust FFI**; reads cached `object$fitted` (conditional) or `object$fixed_fitted` (population).

**What it does:**
| Branch | Trigger | Behavior |
|---|---|---|
| In-sample conditional | `newdata=NULL, re.form=NULL` | Returns `object$fitted` with rownames from `object$model_frame` |
| In-sample population | `newdata=NULL, re.form=NA` or `~0` | Returns `object$fixed_fitted` with rownames |
| Out-of-sample conditional | `newdata=<df>, re.form=NULL` | Calls Rust `predict_new`; allow.new.levels maps to `"population"` vs `"error"` policy |
| Out-of-sample population | `newdata=<df>, re.form=NA` | R-side: model.matrix × beta; no Rust call |

**What it refuses / NA-returns:**
- `interval != "none"`: hard abort (`mm_inference_unavailable`), reason: "not certified by current Rust inference contract"; tagged `prediction_se_unavailable_phase_2`.
- `re.form` other than `NULL`, `NA`, or `~0`: hard abort (`mm_inference_unavailable`), reason: subset RE formulas unsupported by current Rust prediction contract.
- `allow.new.levels` not scalar logical: hard abort (`mm_arg_error`).
- `newdata` missing required formula variables: hard abort (`mm_data_error`).
- `se.fit = TRUE`: returns `NA` se.fit vector (does not abort); the absence reason is surfaced as an attribute rather than a thrown condition. This is a soft partial — caller gets the fit but se.fit is all-NA.

**Missing lme4 arguments (not in signature):**
- `newparams` — counterfactual parameter substitution
- `random.only` — RE-only contribution without FE
- `na.action` — NA-row propagation policy in newdata
- Offset evaluation against newdata in the fixed-only path

**Positional-alignment note (line 293–298):** column alignment of `model.matrix` to `beta` is positional, not by name, because upstream encodes coefficient names in mixeff-specific format (e.g. `"recipe: B"` vs. `"recipeB"`). This means the fixed-only prediction path is brittle if the newdata produces a different column order from training.

---

### 1.2 `fitted.mm_lmm` (S3 — `stats::fitted`)

**File:** `R/predict.R` lines 101–105

**Signature:** `fitted.mm_lmm(object, ...)`

**Returns:** Named numeric vector from `object$fitted`; names = `rownames(object$model_frame)`.

**Rust FFI:** none.

**Notes:** Simple alias over cached field. No additional arguments accepted. Equivalent to `predict(object)` with defaults.

---

### 1.3 `residuals.mm_lmm` (S3 — `stats::residuals`)

**File:** `R/predict.R` lines 109–114

**Signature:** `residuals.mm_lmm(object, type = c("response"), ...)`

**Returns:** Named numeric vector from `object$residuals`; names = `rownames(object$model_frame)`.

**Rust FFI:** none.

**What it refuses:** `match.arg` accepts only `"response"`. Any other type string is an error.

**Missing lme4 arguments:**
- `type = "pearson"` / `"deviance"` / `"working"` — not implemented
- `scaled` argument — absent
- `na.action` / `naresid` re-insertion — absent

---

### 1.4 `fitted.mm_glmm` / `residuals.mm_glmm` (S3)

**File:** `R/predict.R` lines 118–122

Both are direct aliases (`fitted.mm_glmm <- fitted.mm_lmm`; `residuals.mm_glmm <- residuals.mm_lmm`). Inherit all limitations of the LMM versions.

---

### 1.5 `predict.mm_glmm` (S3 — `stats::predict`)

**File:** `R/predict.R` lines 126–132

**Signature:** `predict.mm_glmm(object, ...)`

**Always aborts** with `mm_inference_unavailable`: "GLMM prediction is not certified by the current Rust contract." No computation is attempted. This is an explicit stub/deferred marker.

---

### 1.6 `refit` / `refit.mm_lmm` (generic + S3)

**File:** `R/simulate.R` lines 13–38

**Signature:**
```r
refit(object, newresp, ...)
refit.mm_lmm(object, newresp, ...)
```

**Returns:** A new `mm_lmm` object fitted to the replacement response, with an extra `$refit` list recording `source = "refit"` and `original_fit_status`.

**Rust FFI:** none directly — calls `lmm()` (the top-level fit function) which internally calls Rust. The response column in `object$model_frame` is replaced with `newresp` before refitting.

**Validation:** `newresp` must be numeric, length == `nobs(object)`, no NAs; otherwise `mm_arg_error`.

**Control:** accepts `control` via `...`; defaults to `mm_control(verbose = -1)` if not supplied.

**Notes:** This is the parametric bootstrap building block. It is exported and documented. No `simulate.mm_glmm` / `refit.mm_glmm` exists.

---

### 1.7 `simulate.mm_lmm` (S3 — `stats::simulate`)

**File:** `R/simulate.R` lines 57–84

**Signature:**
```r
simulate.mm_lmm(object, nsim = 1, seed = NULL, re.form = NULL, ...)
```

**Returns:** A `data.frame` with `nrow(object$model_frame)` rows and `nsim` columns named `sim_1`, …, `sim_<nsim>`. Rownames = `rownames(object$model_frame)`. Attributes: `seed` (the supplied seed value), `mm_method = "r_side_gaussian_parametric"`.

**Rust FFI:** none. Simulation is entirely R-side.

**What it does:**
- `re.form = NULL` (default): draws new random effects from the fitted covariance using `mm_simulate_random_mean()`, then adds Gaussian noise at `sigma = fit$sigma`.
- `re.form = NA` or `~0`: uses `fit$fixed_fitted` as the linear predictor, then adds Gaussian noise only. No RE draws.
- `seed`: saved/restored correctly via `mm_with_seed()` — saves `.Random.seed`, calls `set.seed(seed)`, restores on exit.

**Random-effect simulation path (`mm_simulate_random_mean`):**
1. Reads `fit$artifact$semantic_model$random_terms` — list of RE term descriptors from the JSON artifact.
2. For each term, reads the group factor from `fit$model_frame`, extracts levels and basis labels.
3. Builds the covariance matrix from `fit$artifact$covariance_parameter_traces` — reads `standard_deviation` entries first (diagonal), then `correlation` entries (off-diagonal).
4. Falls back to `fit$varcorr$table` for any missing diagonal entries.
5. Draws from the multivariate normal via `mm_rmvnorm` (uses Cholesky; falls back to eigendecomposition + `pmax(vals, 0)` if Cholesky fails — i.e. handles near-singular/boundary covariance gracefully).
6. Accumulates RE contribution into `eta` for each basis.

**What it refuses:**
- `nsim` not a positive scalar numeric: `mm_arg_error`.
- `re.form` anything other than `NULL`, `NA`, or `~0`: `mm_inference_unavailable`.

**Missing lme4 arguments:**
- `use.u` — deprecated lme4 arg mapping to `re.form`; absent
- `newdata` — simulate for a new covariate grid; absent
- `newparams` — counterfactual parameters; absent
- `formula` — formula override; absent
- `family` — family override; absent
- `cluster.rand` — custom RE generator; absent
- `weights` / `offset` — absent
- `allow.new.levels` — absent
- `cond.sim` — absent

**No `simulate.mm_glmm` exists.**

---

## 2. Internal helpers (unexported)

| Helper | Location | Purpose |
|---|---|---|
| `mm_prediction_target(re.form)` | predict.R:134 | Maps `re.form` to `"conditional"`, `"population"`, or `"unsupported"` |
| `mm_predict_newdata(fit, newdata, target, allow_new_levels)` | predict.R:156 | Validates newdata columns, routes to fixed-only or conditional path |
| `mm_predict_conditional_newdata(fit, newdata, allow_new_levels)` | predict.R:188 | Builds Rust FFI call; parses JSON response; validates schema |
| `mm_predict_fixed_only(fit, newdata)` | predict.R:263 | R-side model.matrix × beta; positional alignment |
| `mm_training_xlevels(fit)` | predict.R:302 | Extracts factor levels from training model_frame for model.frame consistency |
| `mm_training_contrasts(fit)` | predict.R:311 | Returns `NULL`; defers to R defaults — contrast argument effectively unused |
| `mm_simulate_once(fit, target)` | simulate.R:86 | Single simulation draw; calls `mm_simulate_random_mean` or uses `fixed_fitted` |
| `mm_simulate_random_mean(fit)` | simulate.R:95 | Accumulates RE draws into linear predictor from artifact traces |
| `mm_random_term_covariance(fit, term_id, basis_labels)` | simulate.R:121 | Reconstructs covariance matrix from artifact traces + varcorr fallback |
| `mm_rmvnorm(n, Sigma)` | simulate.R:171 | MVN draw; Cholesky with eigendecomposition fallback for boundary Sigma |
| `mm_with_seed(seed, expr)` | simulate.R:189 | Save/restore `.Random.seed` around seeded expression |

---

## 3. Stubs and deferred markers

| Location | Marker text / pattern | Nature |
|---|---|---|
| `predict.R:91` | `"prediction_se_unavailable_phase_2"` | `se.fit=TRUE` returns NA; phase-2 deferred |
| `predict.R:55–64` | `interval` != "none" aborts with "not certified by current Rust inference contract" | Interval prediction deferred |
| `predict.R:127–131` | `predict.mm_glmm` always aborts | GLMM prediction is a declared stub |
| `predict.R:293–298` | Positional alignment comment | Fragile alignment; acknowledged as workaround |
| `simulate.R:65–73` | `re.form` subset formula aborts | Partial RE simulation not implemented |
| No `simulate.mm_glmm` | — | GLMM simulation entirely absent |
| No `bootMer` equivalent | — | Parametric bootstrap infrastructure absent (refit exists but no orchestration) |

---

## 4. Parity gap summary vs. lme4

### predict.mm_lmm

| lme4 feature | mixeff status | Classification |
|---|---|---|
| In-sample conditional `predict(m)` | works | works |
| `predict(m, newdata=, re.form=NULL)` | works | works |
| `predict(m, newdata=, re.form=NA)` | works | works |
| `re.form = ~(1\|grp)` subset | refused with structured error | partial |
| `random.only = TRUE` | absent | in-scope-missing |
| `newparams` | absent | in-scope-missing |
| `na.action` propagation in newdata | absent | in-scope-missing |
| `se.fit = TRUE` (actual SEs) | returns NA + reason attr | partial |
| Offset in newdata (fixed-only path) | absent | in-scope-missing |
| GLMM predict | stub abort | partial |

### residuals.mm_lmm

| lme4 feature | mixeff status | Classification |
|---|---|---|
| `type="response"` | works | works |
| `type="pearson"` | absent | in-scope-missing |
| `type="deviance"` | absent | in-scope-missing |
| `type="working"` | absent | in-scope-missing |
| `scaled=TRUE` | absent | in-scope-missing |
| NA re-insertion via `naresid` | absent | in-scope-missing |

### simulate.mm_lmm

| lme4 feature | mixeff status | Classification |
|---|---|---|
| `nsim`, `seed` | works | works |
| `re.form=NULL` (new RE draws) | works | works |
| `re.form=NA` (marginal) | works | works |
| `newdata` | absent | in-scope-missing |
| `newparams` | absent | in-scope-missing |
| `use.u` | absent | out-of-scope-by-design (deprecated lme4 arg) |
| `formula` / `family` override | absent | out-of-scope-by-design |
| `cluster.rand` | absent | out-of-scope-by-design |
| `weights` / `offset` | absent | in-scope-missing |
| `allow.new.levels` | absent | in-scope-missing |
| `cond.sim=FALSE` | absent | out-of-scope-by-design |
| GLMM simulate | absent | partial (GLMM phase) |

### Bootstrap / bootMer

| lme4 feature | mixeff status | Classification |
|---|---|---|
| `refit()` building block | works | works |
| `bootMer()` orchestration | absent | in-scope-missing |
| Parallel bootstrap | absent | in-scope-missing |
| Prediction-interval CI infrastructure | absent | in-scope-missing |
