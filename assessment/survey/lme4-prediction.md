# lme4 / lmerTest — Prediction & Residuals Surface Survey

**Date:** 2026-05-31
**lme4 version surveyed:** 2.0.1 (installed)
**lmerTest version:** 3.2.1
**Purpose:** Exhaustive reference of every user-facing function, argument, and
behavior in the Prediction & Residuals family, against which mixeff is judged.

---

## 1. `predict.merMod`

Full signature (confirmed via `args(lme4:::predict.merMod)`):

```r
predict(object, newdata = NULL, newparams = NULL,
        re.form = NULL, random.only = FALSE,
        terms = NULL,
        type = c("link", "response"),
        allow.new.levels = FALSE,
        na.action = na.pass,
        se.fit = FALSE,
        ...)
```

### 1.1 Arguments

| Argument | What it does | Why users rely on it |
|---|---|---|
| `newdata` | New data frame to predict into. `NULL` returns in-sample fitted values. | Core use: out-of-sample prediction, partial-effects grids, cross-validation. |
| `newparams` | Named list (or vector) of replacement fixed/random parameters; passed to internal `setParams()`. | Counterfactual prediction at hypothetical parameter values without refitting. |
| `re.form` | Formula controlling which random effects to include. `NULL` = all RE from the model; `NA` or `~0` = no RE (population-level/marginal); a one-sided formula like `~(1|group)` = only that grouping factor. | Controls conditional vs. marginal vs. partial-RE prediction; critical for plotting typical subjects vs. average population. |
| `random.only` | `TRUE` = return only the random-effects contribution (i.e., BLUPs × design), dropping the fixed-effects contribution. | Decomposing fitted values into FE and RE parts; diagnostic plots of RE shrinkage. |
| `terms` | Reserved; throws `"terms functionality for predict not yet implemented"` in lme4 2.0.1. | Not currently usable; future partial-effects/terms-only predictions. |
| `type` | `"link"` (default, linear predictor η) or `"response"` (μ = g⁻¹(η)). For LMM the two are identical. For GLMM they differ by the inverse link. | Response-scale vs. link-scale output for GLMM; required for probabilities/counts rather than log-odds/log-rates. |
| `allow.new.levels` | `FALSE` (default) raises an error when `newdata` contains a grouping level not seen in training. `TRUE` silently treats the unseen level's RE as zero (population mean). | Out-of-sample / leave-one-out prediction into held-out groups. |
| `na.action` | Function applied to rows with `NA` in `newdata`. Default `na.pass` returns `NA` predictions for missing rows (preserving row count). | Controls whether NA rows are removed (`na.omit`) or kept (`na.pass`) in the output vector. |
| `se.fit` | `FALSE` = return a plain numeric vector. `TRUE` = return a list `list(fit=..., se.fit=...)` where `se.fit` is the standard error of the linear predictor (or, for GLMM response scale, delta-method transformed). For non-LMM an approximation warning is issued. | Confidence-band plots; downstream use with `qnorm`/`qt` for Wald intervals. |

### 1.2 Return values

- `se.fit = FALSE` (default): named numeric vector, length = number of prediction rows. Names are row names of the model frame (in-sample) or `newdata`.
- `se.fit = TRUE`: named list with elements `$fit` (numeric vector) and `$se.fit` (numeric vector of standard errors). No `df` element; users must supply degrees of freedom themselves for t-based CIs.

### 1.3 Behavioral details

- **In-sample shortcut:** when `newdata = NULL` and `re.form` and `newparams` are all unset, lme4 returns `fitted(object)` for LMM/NLMM or the stored `eta`/`mu` slots for GLMM directly, bypassing matrix arithmetic. This is fast and exact.
- **NA propagation:** `na.action = na.pass` inserts `NA` at any row of `newdata` that contains `NA`; `na.omit` drops those rows and the return vector is shorter.
- **Offset handling:** when `newdata` is supplied, offsets declared via `offset()` in the formula are evaluated against `newdata` and added to the fixed-effects linear predictor.
- **`se.fit` for GLMM:** uses delta method — `se(μ) ≈ se(η) × |g'⁻¹(η)|`, where `g'⁻¹` is the derivative of the inverse link. A warning is issued: "se.fit computation uses an approximation to estimate the sampling distribution of the parameters".
- **Deprecated synonyms:** `ReForm`, `REForm`, `REform` are hard errors since lme4 1.1.

### 1.4 lmerTest interaction

`lmerTest` does not override `predict.merMod`; it uses lme4's method unchanged. Predictions from `lmerTest::lmer()` objects behave identically.

---

## 2. `fitted.merMod`

```r
fitted(object, ...)
```

Returns the in-sample fitted values (conditional on estimated BLUPs) as a named numeric vector. Equivalent to `predict(object)` with all defaults. The names come from `rownames(model.frame(object))`.

---

## 3. `residuals.merMod`

```r
residuals(object,
          type = if (isGLMM(object)) "deviance" else "response",
          scaled = FALSE,
          ...)
```

### 3.1 `type` options (confirmed by testing both LMM and GLMM)

| `type` | LMM | GLMM | Definition |
|---|---|---|---|
| `"response"` | works | works | y − μ (raw residual on response scale) |
| `"pearson"` | works | works | (y − μ) / √Var(y) using the variance function; for LMM equals response/sigma |
| `"deviance"` | works | works (default for GLMM) | signed square-root of the per-observation deviance contribution |
| `"working"` | works | works | working response residual from IRLS; for LMM equals response residuals |
| `"partial"` | **error** | **error** | "partial residuals are not implemented yet" — throws error in lme4 2.0.1 |

### 3.2 `scaled`

`TRUE` divides residuals by `sigma(object)`. Only meaningful for `"response"` type in LMM (standardized residuals). For GLMM, `sigma` = 1 so the division is a no-op.

### 3.3 NA handling

If `na.action` was `na.omit` at fit time, `residuals()` re-inserts `NA` at the omitted positions via `naresid()`, so the returned vector always matches the original data length.

---

## 4. `simulate.merMod`

Full signature (confirmed via `args(lme4:::.simulateFun)`):

```r
simulate(object, nsim = 1, seed = NULL,
         use.u = FALSE,
         re.form = NA,
         newdata = NULL, newparams = NULL,
         formula = NULL, family = NULL,
         cluster.rand = rnorm,
         weights = NULL, offset = NULL,
         allow.new.levels = FALSE,
         na.action = na.pass,
         cond.sim = TRUE,
         ...)
```

Note: the S3 generic `simulate.merMod` delegates to `lme4:::.simulateFun` with
`match.call()`, so all `.simulateFun` arguments are accessible via `simulate()`.

### 4.1 Key arguments

| Argument | What it does |
|---|---|
| `nsim` | Number of simulation draws. Returns a data frame with `nsim` columns named `sim_1`, `sim_2`, …. |
| `seed` | Integer passed to `set.seed()` before simulation; reproducibility. |
| `use.u` | **Deprecated in favour of `re.form`.** `TRUE` = condition on the estimated BLUPs (like `re.form = NULL`). When `FALSE` (default), new RE draws are used. |
| `re.form` | Same semantics as in `predict()`: `NA` = simulate from the marginal/population model (no RE); `NULL` = simulate new RE draws from the fitted covariance. |
| `newdata` | Simulate responses for a new covariate grid rather than the training data. |
| `newparams` | Simulate at counterfactual parameter values. |
| `formula` | Override the model formula for simulation (e.g., a different RE structure). |
| `family` | Override the response family (e.g., to simulate from a different link). |
| `cluster.rand` | Function for generating random effects; default `rnorm` (Gaussian). |
| `weights` | Observation weights for simulation. |
| `offset` | Offset vector for new predictions. |
| `allow.new.levels` | As in `predict()`: allow unseen grouping levels via zero-RE fallback. |
| `cond.sim` | `TRUE` (default) conditions on estimated fixed + random parameters; `FALSE` re-draws the fixed parameters from the asymptotic Gaussian posterior (experimental). |

### 4.2 Return value

A `data.frame` with `nrow(model.frame(object))` rows and `nsim` columns (`sim_1`, …, `sim_nsim`). The class includes `"data.frame"` (used by `refit()` in lme4's parametric bootstrap workflow).

---

## 5. `bootMer` (parametric & semiparametric bootstrap)

```r
bootMer(x, FUN, nsim = 1, seed = NULL,
        use.u = FALSE,
        re.form = NA,
        type = c("parametric", "semiparametric"),
        verbose = FALSE,
        .progress = "none", PBargs = list(),
        parallel = c("no", "multicore", "snow"),
        ncpus = getOption("boot.ncpus", 1L),
        cl = NULL)
```

### 5.1 Purpose and use

`bootMer` generates `nsim` bootstrap replicates by:
- **parametric** (`type = "parametric"`): simulating new responses from the fitted model via `simulate.merMod`, refitting via `refit()`, then applying `FUN` to each refit.
- **semiparametric** (`type = "semiparametric"`): resampling residuals (LMM only, experimental).

`FUN` is any function of a fitted `merMod`; typical uses:
- `FUN = fixef` → bootstrap confidence intervals on fixed effects
- `FUN = function(m) predict(m, newdata=grid)` → simultaneous prediction intervals

### 5.2 Return value

Returns an object of class `c("bootMer", "boot")` compatible with `boot::boot.ci()` and `boot::confint.boot`. Key slots: `$t0` (observed statistic), `$t` (nsim × p matrix of bootstrap statistics), `$R` (= nsim).

### 5.3 Parallelism

`parallel = "multicore"` uses `parallel::mclapply` (Unix only). `"snow"` uses a PSOCK cluster passed via `cl`.

---

## 6. Related extraction functions (used in prediction workflows)

These are not strictly "prediction" functions but are routinely combined with the above:

| Function | Signature | Purpose |
|---|---|---|
| `ranef(object, condVar=FALSE, drop=FALSE, whichel=NULL, postVar=FALSE, ...)` | `lme4:::ranef.merMod` | Returns BLUPs (conditional modes of the random effects). `condVar=TRUE` attaches posterior variances (diagonal of conditional variance matrix) as an attribute — used for caterpillar/dotplot of RE uncertainty. |
| `fixef(object, add.dropped=FALSE)` | `lme4:::fixef.merMod` | Fixed-effect coefficient vector β. `add.dropped=TRUE` re-inserts dropped aliased columns as NA. |
| `getME(object, name)` | generic | Returns internal model matrices. For prediction workflows, key names are `"X"` (FE design matrix), `"Z"` (RE design matrix), `"beta"` (=`fixef`), `"u"` (conditional modes of u = Λᵀu̅), `"b"` (= Λu, BLUPs on original scale). |
| `sigma(object)` | `lme4:::sigma.merMod` | Residual standard deviation σ. Used in standardised residual computation. |
| `vcov(object, correlation=FALSE, use.hessian=NULL)` | `lme4:::vcov.merMod` | Variance-covariance matrix of β. Used internally by `predict(..., se.fit=TRUE)` and downstream by emmeans. |

---

## 7. Simulate-based prediction intervals (bootstrap workflow)

lme4 provides no single function for prediction intervals; the standard pattern is:

```r
# 1. Generate bootstrap predictions
boot_out <- bootMer(fit,
                    FUN = function(m) predict(m, newdata = newgrid),
                    nsim = 999, seed = 1)

# 2. Compute quantile CIs via boot package
ci <- t(apply(boot_out$t, 2, quantile, probs = c(0.025, 0.975)))
```

This yields prediction intervals that account for:
- Uncertainty in fixed effects (β)
- Uncertainty in variance components (θ, σ)
- RE-level uncertainty when `re.form = NULL`

---

## 8. lmerTest extensions (relevant to prediction)

lmerTest does **not** add new prediction, residual, or simulation methods. Its
contributions to inference (Satterthwaite df, `summary()` with p-values,
`anova()` with F-tests) are orthogonal to the prediction surface. The `predict`
method for `lmerModLmerTest` objects dispatches to `lme4:::predict.merMod`
unchanged.

---

## 9. Summary table: full lme4 prediction surface

| Feature | Function/arg | Notes |
|---|---|---|
| In-sample fitted values | `fitted(m)` / `predict(m)` | Identical for LMM with defaults |
| Out-of-sample conditional | `predict(m, newdata=, re.form=NULL)` | Requires known group levels by default |
| Population-level (marginal) | `predict(m, newdata=, re.form=NA)` | Fixed-effects-only; works for any newdata |
| RE subset | `predict(m, re.form=~(1|grp))` | Only that grouping factor's BLUP included |
| RE-only contribution | `predict(m, random.only=TRUE)` | Returns BLUP × Z; no fixed-effect contribution |
| Counterfactual params | `predict(m, newparams=list(...))` | Sensitivity analysis without refitting |
| New factor levels | `predict(m, allow.new.levels=TRUE)` | Zero RE for unseen levels |
| NA handling in newdata | `na.action=na.pass/na.omit` | Controls output vector length |
| GLMM response scale | `predict(m, type="response")` | Applies inverse link; default is link scale |
| Standard errors (Wald) | `predict(m, se.fit=TRUE)` | Delta-method for GLMM; exact for LMM |
| Offset in newdata | automatic via formula | Offset terms evaluated against newdata |
| Response residuals | `residuals(m, type="response")` | y − μ |
| Pearson residuals | `residuals(m, type="pearson")` | (y−μ)/√Var(y) |
| Deviance residuals | `residuals(m, type="deviance")` | Default for GLMM; sqrt of deviance contrib |
| Working residuals | `residuals(m, type="working")` | IRLS working response |
| Partial residuals | `residuals(m, type="partial")` | **Not implemented** in lme4 2.0.1 |
| Scaled residuals | `residuals(m, scaled=TRUE)` | Divides by σ |
| Parametric simulation | `simulate(m, nsim=, seed=, re.form=)` | Returns data.frame; used by bootMer/refit |
| Marginal simulation | `simulate(m, re.form=NA)` | No RE draws; from population mean |
| New-data simulation | `simulate(m, newdata=)` | Simulate for a different covariate grid |
| Parametric bootstrap | `bootMer(m, FUN, nsim=, type="parametric")` | CI for any function of the model |
| Semiparametric bootstrap | `bootMer(m, type="semiparametric")` | Residual resampling; experimental |
| Parallel bootstrap | `bootMer(..., parallel="multicore")` | Speedup via `parallel` package |
| Posterior RE variances | `ranef(m, condVar=TRUE)` | Uncertainty in BLUPs; for caterpillar plots |
| BLUPs | `ranef(m)` | Conditional modes of u |

---

## 10. mixeff coverage assessment

This section documents each lme4 surface item against what mixeff currently
implements (per `R/predict.R`, `R/simulate.R`, `NAMESPACE`, and
`tests/testthat/test-predict-newdata.R`).

### predict / fitted

| Item | mixeff status | Notes |
|---|---|---|
| `predict(m)` in-sample conditional | **works** | Returns `object$fitted` |
| `predict(m, newdata=, re.form=NULL)` | **works** | Via Rust `predict_new` FFI |
| `predict(m, newdata=, re.form=NA)` | **works** | Fixed-only via R-side matrix multiply |
| `predict(m, re.form=~(1|grp))` | **partial** | Raises `mm_inference_unavailable`; documented as unsupported |
| `predict(m, random.only=TRUE)` | **missing** | Argument not in signature |
| `predict(m, newparams=)` | **missing** | Argument not in signature |
| `predict(m, type="link"/"response")` | **works** | Accepted; LMM-only (link=response for LMM) |
| `predict(m, allow.new.levels=TRUE)` | **works** | Routes to `NewReLevels::population` policy |
| `predict(m, allow.new.levels=FALSE)` | **works** | Raises `mm_inference_unavailable` on unseen level |
| `predict(m, se.fit=TRUE)` | **partial** | Returns `NA` se.fit with `mm_unavailable_reason` attribute; not computed |
| `predict(m, na.action=)` | **missing** | Argument not in signature; no NA propagation in newdata |
| Offset in newdata | **missing** | `mm_predict_fixed_only` does not evaluate `offset()` terms against newdata |
| `terms` argument | **missing** (same as lme4; lme4 itself doesn't implement it) | Not applicable |
| `fitted(m)` | **works** | Returns `object$fitted` with rownames |
| GLMM `predict` | **missing** | `predict.mm_glmm` always raises `mm_inference_unavailable` |

### residuals

| Item | mixeff status | Notes |
|---|---|---|
| `residuals(m, type="response")` | **works** | Returns `object$residuals` (single accepted type) |
| `residuals(m, type="pearson")` | **missing** | Only "response" is in `match.arg` |
| `residuals(m, type="deviance")` | **missing** | Not implemented |
| `residuals(m, type="working")` | **missing** | Not implemented |
| `residuals(m, scaled=TRUE)` | **missing** | `scaled` argument absent |
| NA re-insertion via `naresid` | **missing** | No `na.action` tracking in residuals |
| GLMM `residuals` | **partial** | `residuals.mm_glmm` delegates to `residuals.mm_lmm` — same single-type limitation |

### simulate

| Item | mixeff status | Notes |
|---|---|---|
| `simulate(m, nsim=, seed=)` | **works** | Returns correctly shaped `data.frame`; seed is restored |
| `simulate(m, re.form=NULL)` conditional | **works** | Draws new RE from fitted covariance |
| `simulate(m, re.form=NA)` marginal | **works** | Uses `fixed_fitted` only |
| `simulate(m, use.u=TRUE)` | **missing** | `use.u` not in signature |
| `simulate(m, newdata=)` | **missing** | `newdata` argument absent |
| `simulate(m, newparams=)` | **missing** | `newparams` argument absent |
| `simulate(m, re.form=~(...))` formula subset | **missing** | Raises `mm_inference_unavailable` |
| `simulate(m, cond.sim=FALSE)` | **missing** | Not implemented |
| `simulate(m, cluster.rand=)` | **missing** | Not implemented |
| `simulate(m, weights=, offset=)` | **missing** | Not implemented |
| GLMM simulate | **missing** | No `simulate.mm_glmm` |

### bootMer / prediction intervals

| Item | mixeff status | Notes |
|---|---|---|
| `bootMer()` parametric bootstrap | **missing** | No equivalent function |
| `bootMer()` semiparametric bootstrap | **missing** | No equivalent function |
| Simulate-based prediction intervals | **missing** | No CI infrastructure for intervals |
| `se.fit` = actual computed SE | **missing** | `se.fit=TRUE` returns `NA` with a reason attribute |

### ranef / BLUPs

| Item | mixeff status | Notes |
|---|---|---|
| `ranef(m)` BLUPs | **works** | Via `ranef.mm_lmm` / `ranef.mm_glmm` |
| `ranef(m, condVar=TRUE)` posterior variances | **works** | Implemented via Rust `cond_var()` FFI in `methods-extract.R`; results cached; parity test against lme4 within 1e-3 tolerance |
