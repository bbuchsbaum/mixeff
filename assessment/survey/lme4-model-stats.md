# lme4/lmerTest Model-Statistics Generics — Reference Survey

**Surveyed:** lme4 2.0.1, lmerTest 3.2.1  
**Date:** 2026-05-31  
**Purpose:** Canonical lme4 surface for the "Model-statistics generics" capability family, against which mixeff is assessed. Does **not** assess mixeff.

---

## 1. Core Likelihood Generics

### `logLik(object, REML = NULL, ...)`

Registered as `logLik.merMod`. Extracts the log-likelihood from a fitted `lmerMod` or `glmerMod` object.

**`REML` argument semantics (LMM only):**
- `REML = NULL` (default): returns the criterion the model was optimized under — REML log-likelihood if the model was fitted with `REML = TRUE`, ML log-likelihood otherwise.
- `REML = TRUE`: forces return of the REML log-likelihood regardless of how the model was fitted (re-evaluates at the stored parameters).
- `REML = FALSE`: forces return of the ML log-likelihood regardless of how the model was fitted (re-evaluates at the stored parameters).

**Return value:** An object of class `"logLik"` (S3, inherits from `numeric`). Attributes:
- `df`: integer — total number of estimated parameters (fixed effects + variance components + residual variance for LMM). For a random-slopes LMM with 2 fixed effects and 3 variance components, `df = 6`.
- `nobs`: integer — number of observations used in fitting (equals `N`; equal to `nall` attribute as well, since lme4 does not use case weights by default).
- `nall`: integer — same as `nobs` in the unweighted case; present for compatibility with the `logLik` generic contract.
- `class`: `"logLik"`.

**GLMMs:** `REML` argument is ignored; always returns the Laplace-approximated marginal log-likelihood. `df` counts estimated parameters (fixed effects + variance components; no residual variance for canonical GLMMs with known dispersion).

**Users rely on this for:** likelihood-ratio tests via `anova()`, manual LRT computation, AIC/BIC derivation, `bbmle`/`AICcmodavg` workflows.

---

### `deviance(object, REML = FALSE, ...)`

Registered as `deviance.merMod`. Returns `-2 * logLik(object, REML = FALSE)` — the ML deviance — regardless of how the model was fitted.

**`REML` argument:** The argument is accepted but **does not change the output** for a REML-fitted model in lme4 2.0.1; both `deviance(m_reml, REML=TRUE)` and `deviance(m_reml, REML=FALSE)` return the ML deviance `-2*logLik(m, REML=FALSE)`. The REML criterion is accessed via `REMLcrit()`, not `deviance()`.

**GLMMs:** Returns the deviance from the Laplace approximation (not the saturated-model deviance). Equal to `-2 * logLik(object)`.

**Relationship to `REMLcrit`:** `deviance()` ≠ `REMLcrit()`. Specifically:
- `deviance(m_reml)` = `-2 * logLik(m_reml, REML=FALSE)` ≈ 1751.99 (sleepstudy example)
- `REMLcrit(m_reml)` = `-2 * logLik(m_reml, REML=TRUE)` ≈ 1743.63 (sleepstudy example)

**Users rely on this for:** model comparison in likelihood-ratio test workflows; extracting the ML objective for GLMM diagnostics.

---

### `REMLcrit(object)`

lme4-specific (not a base-R generic). Returns the REML criterion: `-2 * logLik(object, REML = TRUE)`. Defined only for `merMod` objects. Errors or gives meaningless output on ML-fitted models.

**Return value:** A scalar numeric (not a `logLik` object).

**Users rely on this for:** directly accessing the REML optimization criterion; confirming that the optimizer minimized the REML criterion to the expected value; sometimes used in custom convergence checks.

---

## 2. Information Criteria

### `AIC(object, ..., k = 2)`

Uses base-R `AIC` generic dispatched via `logLik.merMod`. Formula: `-2 * logLik + k * df`.

**Multi-model form:** `AIC(m1, m2, ...)` returns a data frame with columns `df` and `AIC`, one row per model. lme4 issues a **warning** when comparing a REML-fitted model to an ML-fitted model (or REML to REML with different fixed effects), because REML likelihoods are not comparable across models with different fixed-effect specifications.

**`k` argument:** Default 2 (standard AIC). Setting `k = log(n)` gives BIC-equivalent penalty.

**Return value (single model):** A named numeric scalar. **Return value (multi-model):** A data frame with `df` and `AIC` columns, row names = deparsed call expressions.

**For REML fits:** `AIC` is computed from the REML log-likelihood (i.e. `REML = NULL` default), which is the criterion value under which the model was optimized. This is sometimes considered incorrect for model selection (ML is preferred for fixed-effect comparison); lme4 warns but does not error.

**Users rely on this for:** model selection via AIC, `bbmle::AICtab()`, `MuMIn::model.sel()`.

---

### `BIC(object, ...)`

Uses base-R `BIC` generic dispatched via `logLik.merMod`. Formula: `-2 * logLik + log(nobs) * df`.

**Return value:** Same structure as `AIC()`. Same multi-model form and REML warning behavior.

**Users rely on this for:** model selection with stronger complexity penalty than AIC; BIC is preferred when the goal is identifying the true model rather than minimizing prediction error.

---

### `extractAIC(fit, scale = 0, k = 2, ...)`

Registered as `extractAIC.merMod`. Returns a length-2 numeric vector: `c(df, AIC)`, where `df` is the number of estimated parameters (same as `attr(logLik(fit), "df")`) and `AIC = -2 * logLik + k * df`.

**Difference from `AIC()`:** `extractAIC()` always returns a 2-element vector; does not accept multiple models. It is used internally by `step()` and related model-selection utilities in base R.

**`scale` argument:** Ignored for `merMod` objects (scale is always derived from the model).

**Users rely on this for:** compatibility with `step()`, `drop1()`, and other base-R stepwise selection tools that call `extractAIC` internally.

---

### `llikAIC(object, devianceFUN = devCrit, chkREML = TRUE, devcomp = object@devcomp)`

lme4-internal (exported but primarily used internally). Returns a list with:
- `logLik`: the `logLik` object
- `AICtab`: a named numeric vector — for REML fits, contains `REML` = the REML criterion; for ML fits, contains `AIC`, `BIC`, `logLik`, `deviance`, `df.resid`.

**Users rely on this for:** computing the AIC table that appears in `summary.merMod`; also used by `anova.merMod` to build comparison tables. Not typically called directly by users.

---

## 3. Observation and Degrees-of-Freedom Counts

### `nobs(object, ...)`

Registered as `nobs.merMod`. Returns the number of observations as an integer scalar. Equal to `getME(object, "N")` for unweighted fits. For LMMs and GLMMs this is the number of rows in the data after NA removal.

**Users rely on this for:** BIC computation (uses `log(nobs)`), sample-size checks, reporting.

---

### `df.residual(object, ...)`

Registered as `df.residual.merMod`. Returns `nobs(object) - length(fixef(object))` — i.e., `N - p` where `p` is the number of fixed-effect parameters.

**Important:** This does **not** subtract random-effect parameters or variance components. It is a conventional approximation; lme4 docs explicitly caution that the "correct" residual df for mixed models is not well-defined. For a 180-observation sleepstudy model with 2 fixed effects, `df.residual = 178` (not 174 as the REML denominator would suggest under Satterthwaite).

**GLMMs:** Same formula `N - p`.

**Users rely on this for:** conventional residual df in downstream tools (e.g., `confint` with `method="Wald"` uses `df.residual` from the underlying object in some paths); t-distribution reference in non-lmerTest workflows.

---

## 4. REML/ML State Inspection and Refitting

### `isREML(object, ...)`

Registered as `isREML.merMod`. Returns `TRUE` if the model was fitted with REML, `FALSE` otherwise. Internally reads `getME(object, "REML")` (an integer: `2` = REML, `0` = ML).

**Users rely on this for:** programmatic checks before likelihood comparisons; guards in custom functions that require ML fits.

---

### `refitML(object, optimizer = "bobyqa", ...)`

Registered as `refitML.merMod`. Refits a REML-fitted model by ML (equivalent to calling `lmer(..., REML=FALSE)` with the same formula and data). Returns a new `lmerMod` object with `isREML() == FALSE`.

**Users rely on this for:** converting REML fits to ML before likelihood-ratio tests (`anova(m1, m2)` with `refit=TRUE` calls this internally); also used when `AIC`/`BIC` comparisons are desired under ML.

---

### `refit(object, newresp, ...)`

Registered as `refit.merMod`. Refits the model with a new response vector `newresp` (same formula, random-effect structure, data). Used in parametric bootstrap loops and simulation-based inference.

**Return value:** A new fitted `merMod` object of the same class.

**Users rely on this for:** parametric bootstrap (e.g., `bootMer`); sensitivity analyses where the response is permuted or simulated.

---

## 5. Internal Component Extraction via `getME`

### `getME(object, name, ...)`

Registered as `getME.merMod`. Extracts named low-level components from a fitted `merMod`. `name` can be a single string or a character vector (returns a named list for multiple names).

**Full list of valid `name` values** (confirmed in lme4 2.0.1):

| Name | Type | Description |
|------|------|-------------|
| `"X"` | `matrix` | Fixed-effect design matrix (dense, N×p) |
| `"Z"` | `dgCMatrix` | Random-effect design matrix (sparse, N×q) |
| `"Zt"` | `dgCMatrix` | Transpose of Z (q×N) |
| `"Ztlist"` | `list` | Z split by random-effect term |
| `"mmList"` | `list` | Model-matrix list by random-effect term |
| `"y"` | `numeric` | Response vector |
| `"mu"` | `numeric` | Fitted (conditional mean) values |
| `"u"` | `numeric` | Conditional modes of random effects (scaled) |
| `"b"` | `dgeMatrix` | Conditional modes of random effects (unscaled) |
| `"Gp"` | `integer` | Group pointer vector for random-effect blocks |
| `"Tp"` | `numeric` | Theta pointer |
| `"L"` | `dCHMsimpl` | Cholesky factor of the downdated cross-product |
| `"Lambda"` | `dgCMatrix` | Relative covariance factor (lower triangular) |
| `"Lambdat"` | `dgCMatrix` | Transpose of Lambda |
| `"Lind"` | `integer` | Index vector mapping theta to Lambda entries |
| `"Tlist"` | `list` | List of lower-triangular theta matrices |
| `"A"` | `dgCMatrix` | `Lambdat %*% Zt` (used in penalized LS) |
| `"RX"` | `matrix` | Cholesky factor of the fixed-effect information |
| `"RZX"` | `matrix` | Cross-product factor |
| `"sigma"` | `numeric` | Residual standard deviation |
| `"flist"` | `list` | Named list of grouping factors (as factors) |
| `"fixef"` | `numeric` | Fixed-effect coefficient vector (alias for `beta`) |
| `"beta"` | `numeric` | Fixed-effect coefficient vector |
| `"theta"` | `numeric` | Variance-component parameter vector (Cholesky factors) |
| `"ST"` | `list` | S and T matrices for `theta` reparameterization |
| `"par"` | `numeric` | Full parameter vector passed to optimizer |
| `"REML"` | `integer` | `2L` if REML, `0L` if ML |
| `"is_REML"` | `logical` | `TRUE` if REML, `FALSE` if ML |
| `"n_rtrms"` | `integer` | Number of random-effect terms |
| `"n_rfacs"` | `integer` | Number of grouping factors |
| `"N"` | `integer` | Total number of observations |
| `"n"` | `integer` | Number of observations (= N for LMM) |
| `"p"` | `integer` | Number of fixed-effect parameters |
| `"q"` | `integer` | Total number of random-effect coefficients |
| `"p_i"` | `integer` | Number of columns per random-effect term |
| `"l_i"` | `integer` | Number of levels per grouping factor |
| `"q_i"` | `integer` | Number of random coefficients per term |
| `"k"` | `integer` | Number of random-effect terms |
| `"m_i"` | `integer` | Number of theta parameters per term |
| `"m"` | `integer` | Total number of theta parameters |
| `"cnms"` | `list` | Column names of random-effect model matrices per term |
| `"devcomp"` | `list` | Full deviance components list (see `devcomp()` below) |
| `"offset"` | `numeric` | Model offset vector |
| `"lower"` | `numeric` | Lower bounds for theta in optimizer |
| `"devfun"` | `function` | The optimizer's objective (deviance) function |
| `"devarg"` | `list` | Arguments passed to `devfun` |
| `"glmer.nb.theta"` | `numeric` | Negative-binomial theta (glmer.nb only) |

**Users rely on this for:** accessing internal Cholesky factors (`Lambda`, `L`), grouping structure (`flist`, `cnms`), optimizer internals (`devfun`, `lower`), and matrix components for downstream linear algebra (e.g., hat matrix computation, custom standard errors).

---

## 6. `devcomp` — Deviance Components

### `devcomp(object, ...)`

lme4-specific (exported). Returns a list with two elements:

**`$cmp` — numeric named vector of deviance components:**
| Name | Description |
|------|-------------|
| `ldL2` | log-determinant of L (the Cholesky factor), times 2 |
| `ldRX2` | log-determinant of RX (fixed-effect Cholesky), times 2 |
| `wrss` | weighted residual sum of squares |
| `ussq` | sum of squares of scaled random effects (`u'u`) |
| `pwrss` | penalized weighted residual sum of squares (`wrss + ussq`) |
| `drsum` | deviance residual sum (for GLMMs) |
| `REML` | REML criterion value (if fitted by REML; `NA` for ML) |
| `dev` | ML deviance (`-2 * logLik_ML`; `NA` for REML fits) |
| `sigmaML` | ML estimate of residual standard deviation |
| `sigmaREML` | REML estimate of residual standard deviation |

**`$dims` — integer named vector of model dimensions:**
| Name | Description |
|------|-------------|
| `N` | total observations |
| `n` | observations (= N for LMM) |
| `p` | fixed-effect parameters |
| `nmp` | `n - p` |
| `q` | total random-effect coefficients |
| `nth` | number of theta parameters |
| `useSc` | `1` if model uses a scale parameter (LMM/Gaussian), `0` otherwise |
| `reTrms` | number of random-effect terms |
| `spFe` | `1` if sparse fixed-effect design, `0` otherwise |
| `REML` | `2` if REML, `0` if ML |
| `GLMM` | `1` if GLMM, `0` if LMM |
| `NLMM` | `1` if NLMM, `0` otherwise |

**Users rely on this for:** debugging optimizer internals; computing the PWRSS (penalized weighted RSS) for custom variance-component estimates; accessing log-determinants for custom information criteria.

---

## 7. Ancillary Model-Statistics Utilities

### `devfun2(fm, useSc = TRUE, scale = c("sdcor", "varcov"), ...)`

Returns a reparameterized deviance function (class `"devfun"`) over variance-component parameters in either standard-deviation/correlation (`"sdcor"`) or variance/covariance (`"varcov"`) scale. Used to profile the likelihood over interpretable parameterizations.

**Attributes of return value:** `useSc`, `scale`, and optionally `lower`/`upper` bounds if constrained.

**Users rely on this for:** profiling the likelihood for confidence intervals (`profile.merMod` calls this); custom optimization over interpretable variance-component scales.

---

### `varianceProf(x, ranef = TRUE)`

Converts a `profile.merMod` object (`x`) to a variance-parameterization profile. Used for converting standard-deviation profiles to variance-scale profiles.

**Users rely on this for:** post-processing profile objects when reporting variance (rather than SD) confidence intervals.

---

## 8. lmerTest Additions

### `anova.lmerModLmerTest(object, ..., type = c("III","II","I","3","2","1"), ddf = c("Satterthwaite","Kenward-Roger","lme4"))`

lmerTest overrides `anova()` for `lmerMod` objects promoted to `lmerModLmerTest` class. The key additions to the model-statistics surface:
- **`type` argument:** Type I (sequential), II, or III sums of squares. lme4's `anova()` only does sequential (Type I).
- **`ddf` argument:** Denominator degrees of freedom method — `"Satterthwaite"` (default), `"Kenward-Roger"`, or `"lme4"` (infinite df / z-test). This determines how `df` and `p-value` columns in the anova table are computed for the model terms.

**Users rely on this for:** F-tests with Satterthwaite or Kenward-Roger df in multi-term anova tables; this is lmerTest's primary user-facing contribution.

### `as_lmerModLmerTest(model, tol = 1e-08)`

Promotes a fitted `lmerMod` to `lmerModLmerTest` class, computing Satterthwaite gradient information. Required before lmerTest's `anova()`, `summary()` with df, and `contest()` can be used.

**Users rely on this for:** applying lmerTest inference methods to models originally fitted with `lme4::lmer()` rather than `lmerTest::lmer()`.

---

## 9. Behavioral Notes Critical for Parity Assessment

1. **`logLik.merMod` has a REML override argument.** `logLik(m, REML=FALSE)` on a REML-fitted model re-evaluates the ML log-likelihood at the stored REML parameters — it does **not** refit. The returned value is therefore the ML likelihood at REML-optimal theta, not the ML optimum. This is different from `refitML()` which optimizes.

2. **`deviance()` always returns ML deviance.** The `REML` argument on `deviance.merMod` is accepted but its effect in lme4 2.0.1 is to return the ML deviance in both cases (not the REML criterion). `REMLcrit()` is the correct accessor for the REML objective.

3. **`AIC()` uses whichever logLik the model was fitted under (REML or ML).** This is intentional — AIC computed on a REML log-likelihood is not directly comparable to AIC on ML for models with different fixed effects, but lme4 emits a warning rather than an error, leaving it to the user to decide.

4. **`AIC(m1, m2)` multi-model call uses `logLik` dispatch via base R's `AIC` generic.** lme4 does not define a separate `AIC.merMod`; the behavior for single-model and multi-model calls comes from the base R `AIC` generic applied to `logLik.merMod`.

5. **`df.residual` uses `N - p` (fixed-effects rank only).** It does not account for random-effect parameters. This is widely known to be an approximation for mixed models.

6. **`logLik` has an `nall` attribute** (in addition to `nobs` and `df`). In the unweighted case `nall == nobs`. Downstream tools like `bbmle` inspect `nall`.

7. **`getME` is an lme4-specific generic, not a base-R generic.** Packages that want to use it on non-lme4 objects must either re-export it from lme4 or define their own generic. When mixeff defines its own `getME`, it shadows lme4's.

8. **`refit()` is an lme4 generic.** It accepts a `newresp` argument (new response vector), refitting the full model at the new y. Distinct from `refitML()` which changes the estimation criterion.

9. **`REMLcrit()` is the only lme4 function that directly exposes the REML criterion as a scalar.** `logLik(m, REML=TRUE)` returns a `logLik` object (with `df` and `nobs` attributes); `REMLcrit(m)` returns `numeric(1)` equal to `-2 * as.numeric(logLik(m, REML=TRUE))`.
