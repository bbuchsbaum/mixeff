# lme4 / lmerTest Public Surface: `simulate` / `refit` / `update`

Survey date: 2026-05-31  
lme4 version: 2.0.1 (installed)  
lmerTest version: 3.2.1 (installed)  
Scope: every user-facing function, argument, and behaviour in this family.

---

## 1. `simulate.merMod` — draw responses from a fitted model

### Signature

```r
simulate(object, nsim = 1, seed = NULL,
         use.u = FALSE, re.form = NA,
         newdata = NULL, newparams = NULL, family = NULL,
         cluster.rand = rnorm,
         allow.new.levels = FALSE, na.action = na.pass, ...)
```

### What it does

Draws `nsim` response vectors from the conditional distribution of a fitted `lmerMod` or `glmerMod` object, returning a data frame with one column per simulated dataset (columns named `sim_1`, `sim_2`, …).

### Arguments

| Argument | Role / why users rely on it |
|---|---|
| `object` | The fitted `merMod` model to simulate from. |
| `nsim` | Number of independent simulated response vectors. Vectorised simulation in one call is far more efficient than looping. |
| `seed` | Passed to `set.seed()` before simulation; ensures reproducible samples. The actual `.Random.seed` state is stored as an attribute on the returned data frame so the exact RNG state can be recovered. |
| `use.u` | Legacy alias for `re.form`. `TRUE` ↔ `re.form = NULL` (condition on fitted random-effect estimates); `FALSE` ↔ `re.form = ~0` (simulate new random effects). Kept for backward compatibility; `re.form` is preferred. |
| `re.form` | Controls which random effects are conditioned on vs. simulated from scratch. `NULL` = condition on all (conditional simulation); `NA` or `~0` = marginalise over all (population-level simulation); a partial formula = condition on specified random effects only and simulate the rest. Default is `NA`. |
| `newdata` | A data frame of predictor values at which to generate predictions/simulations, enabling out-of-sample simulation from the fitted parameter estimates. |
| `newparams` | Named list with components `theta`, `beta`, and optionally `sigma` — allows simulation from a *hypothetical* (user-specified) parameter vector rather than the fitted estimates. Key for power analysis and model checking. |
| `family` | Override the response family (e.g., swap negative-binomial theta). |
| `cluster.rand` | A function of `n` that generates standardised random cluster effects (default `rnorm`). Supports testing robustness to non-normal random effects (e.g., heavy-tailed, truncated normal, mixture). |
| `allow.new.levels` | If `TRUE`, new grouping-factor levels in `newdata` receive population-level (zero) random effects instead of triggering an error. Required for out-of-sample simulation. |
| `na.action` | How to handle `NA` in `newdata`; default `na.pass` passes them through. |

### Return value

A `data.frame` with `nrow(model.frame(object))` rows and `nsim` columns. The `.Random.seed` state at call time is stored as `attr(result, "seed")`. Each column can be passed directly to `refit()`.

### Key behaviours

- **Reproducibility contract**: seeding is done by storing and restoring `.Random.seed` so calling code is not permanently affected.
- **NA handling in original response**: if `newresp` (from a simulate call) carries an `na.action` attribute, `refit()` knows NAs have been removed and does not require length-matching.
- **GLMM support**: generates binomial proportions (or success/failure matrices), Poisson counts, negative-binomial counts, Gamma values, etc., matching the fitted family.
- **`simulate.formula` variant** (see §4): can simulate from a formula + `newdata` + `newparams` without a fitted model at all.

---

## 2. `simulate.formula` — simulate without a fitted model

### Signature (stats generic, lme4 method)

```r
simulate(object, nsim = 1, seed = NULL, ..., basis, newdata, data)
```

### What it does

Dispatches to `.simulateFun()` when `object` is a one-sided mixed-model formula. Generates responses from scratch given user-specified parameters, design data, and family. No fitted model is required.

### `.simulateFun()` — the workhorse

```r
.simulateFun(object, nsim = 1, seed = NULL, use.u = FALSE,
             re.form = NA,
             newdata = NULL, newparams = NULL,
             formula = NULL, family = NULL,
             cluster.rand = rnorm,
             weights = NULL, offset = NULL,
             allow.new.levels = FALSE, na.action = na.pass,
             cond.sim = TRUE, ...)
```

Additional arguments beyond `simulate.merMod`:

| Argument | Role |
|---|---|
| `formula` | One-sided mixed-model formula; allows simulation when no fitted object exists. |
| `weights` | Prior weights (as in `lmer`/`glmer`). Required for binomial proportions with non-1 denominators. |
| `offset` | Offset term (as in `glmer`). |
| `cond.sim` | (Experimental) if `FALSE`, return only the group-level random-effect draws without adding conditional residuals. Useful for inspecting just the random-effects realisation. |

### Usage pattern

```r
params <- list(theta = 0.5, beta = c(2, -1, -2, -3))
simdat <- ...
form   <- formula(gm1)[-2]   # drop LHS to get one-sided formula
simulate(form, newdata = simdat, family = binomial, newparams = params)
```

This is the primary path for power analysis and simulation studies where the user constructs the data-generating process manually.

---

## 3. `refit` — refit the same model to a new response

### Signature

```r
refit(object, newresp, ...)

## S3 method for class 'merMod'
refit(object, newresp = NULL, newweights = NULL,
      rename.response = FALSE,
      maxit = 100, ...)
```

### What it does

Refits a `merMod` model with a new response vector (and optionally new weights), bypassing model-frame construction and formula parsing to go directly to the optimisation step. Intended as the fast inner loop of parametric bootstrap or simulation-study workflows.

### Arguments

| Argument | Role / why users rely on it |
|---|---|
| `object` | The reference fitted model whose structure (formula, random-effects grouping, design matrices) is reused. |
| `newresp` | New numeric response vector (same length as original) **or** a single-column data frame as produced by `simulate(object)`. May carry an `na.action` attribute (see NA handling below). |
| `newweights` | Optional replacement prior-weights vector. Supports refit with changed observation weights. |
| `rename.response` | When `TRUE`, replaces the response variable name in the formula/model-frame with the name of `newresp`. Needed if downstream operations (e.g., a further `update()`) must see the new name. Note: `terms` component is not updated, so the result is slightly different from `update()`. |
| `maxit` | For GLMMs only: maximum number of PWRSS update iterations in the inner loop. |
| `control` | Passed via `...`; allows changing optimizer/convergence settings for the refit. |

### NA handling

If `newresp` has an `na.action` attribute, lme4 assumes NAs have already been removed (the vector is shorter than the original). This lets `simulate(object)[, i]` (which inherits the model's `na.action`) be passed directly even when the original response had `NA` values.

### Return value

An object of the same class as `object`, fit to the new response. Computationally cheaper than `update()` for the common case of same formula + same data structure.

### See also

- `refitML()` — specialised refit to change REML→ML criterion (§4).
- `bootMer()` — wraps the simulate + refit cycle into a bootstrap harness (§5).

---

## 4. `refitML` — refit switching from REML to ML

### Signature

```r
refitML(x, ...)

## S3 method for class 'merMod'
refitML(x, optimizer = "bobyqa", ...)
```

### What it does

Refits an LMM that was estimated by REML using the maximum-likelihood criterion. Primarily used internally by `anova.merMod` (which by default sets `refit = TRUE`) so that LMMs with different fixed effects can be compared on a valid likelihood scale.

### Arguments

| Argument | Role |
|---|---|
| `x` | A REML-fitted `lmerMod`. |
| `optimizer` | String naming the optimizer; default `"bobyqa"`. |

---

## 5. `update.merMod` — refit with modified formula or arguments

### Signature

```r
## S3 method for class 'merMod'
update(object, formula., ..., evaluate = TRUE)
```

### What it does

Reconstructs and re-evaluates the original model call with modifications to the formula and/or any other arguments (`data`, `control`, `REML`, `weights`, `offset`, `subset`, etc.). More flexible than `refit()` but slower because it rebuilds the full model representation.

### Arguments

| Argument | Role / why users rely on it |
|---|---|
| `object` | A fitted `merMod`. |
| `formula.` | Formula update specification (passed to `update.formula()`). Standard R update notation: `. ~ . + x` adds a term, `. ~ . - x` removes one. Can change both fixed and random parts. |
| `...` | Any named argument accepted by the original fitting function (`lmer`/`glmer`). Common uses: change `data`, `control`, `REML`, `weights`, `offset`, `subset`, `na.action`. |
| `evaluate` | If `FALSE`, return the modified unevaluated call rather than fitting the model. Useful for inspection or deferred evaluation. |

### Key behaviours

- Modifies the stored call, then evaluates it in the formula's environment (falling back to the calling frame).
- Can change random-effects structure, add/remove fixed-effects terms, switch REML on/off, or swap the data source.
- Used internally by `allFit()` to re-run with different optimizers.

### lmerTest override: `update.lmerModLmerTest`

lmerTest provides its own `update.lmerModLmerTest` method with the same signature. It delegates to the same call-manipulation logic, then calls `as_lmerModLmerTest()` on the result if needed, ensuring the updated model is upgraded to `lmerModLmerTest` class (so Satterthwaite df, `ranova()`, etc., are available). lmerTest does **not** add `simulate` or `refit` methods for `lmerModLmerTest`; those dispatch to the lme4 `merMod` methods.

---

## 6. `bootMer` — model-based (semi-)parametric bootstrap

### Signature

```r
bootMer(x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA,
        type = c("parametric", "semiparametric"),
        verbose = FALSE, .progress = "none", PBargs = list(),
        parallel = c("no", "multicore", "snow"),
        ncpus = getOption("boot.ncpus", 1L), cl = NULL)
```

### What it does

The "simulate-then-estimate" harness: for each of `nsim` replicates it (i) simulates a response from the fitted model via `simulate.merMod`, (ii) refits the model via `refit()`, and (iii) applies `FUN` to the refitted object to extract statistics of interest. Returns a `boot`-class object compatible with `boot::boot.ci()`.

### Arguments

| Argument | Role / why users rely on it |
|---|---|
| `x` | Fitted `merMod`. |
| `FUN` | User-supplied summary function; takes a fitted `merMod` and returns a named numeric vector (e.g., fixed effects, variance components, custom contrasts). |
| `nsim` | Bootstrap replications. |
| `seed` | Reproducibility. |
| `use.u` / `re.form` | Same semantics as `simulate.merMod`: whether to condition on estimated random effects (`use.u = TRUE` / `re.form = NULL`) or simulate new ones (`use.u = FALSE` / `re.form = NA`). |
| `type` | `"parametric"`: simulate from fitted Gaussian/GLMM distribution. `"semiparametric"`: resample residuals (LMM only; experimental). |
| `verbose` | Print progress. |
| `.progress` / `PBargs` | Progress bar type and options. |
| `parallel` | `"no"`, `"multicore"` (Unix fork), or `"snow"` (cluster). Enables parallel bootstrap. |
| `ncpus` | Number of CPU cores for parallel execution. |
| `cl` | Pre-existing `snow`/`parallel` cluster object. |

### Return value

An S3 object of class `"boot"` (from the `boot` package). Error/warning metadata during bootstrap is stored as attributes `bootFail`, `boot.fail.msgs`, `boot.all.msgs`.

### Primary use cases

- Parametric-bootstrap confidence intervals for variance components (where profile-likelihood may be unreliable near boundaries).
- Bootstrap p-values for fixed effects (vs. Satterthwaite/Kenward-Roger).
- Any user-defined statistic drawn from refitted models.

---

## 7. `allFit` — convergence check across optimizers

### Signature

```r
allFit(object, meth.tab = NULL, data = NULL,
       verbose = TRUE, show.meth.tab = FALSE,
       maxfun = 1e5,
       parallel = c("no", "multicore", "snow"),
       ncpus = getOption("allFit.ncpus", 1L), cl = NULL,
       catch.errs = TRUE,
       start_from_mle = TRUE)
```

### What it does

Calls `update()` repeatedly on a fitted model, substituting each available optimizer. Returns an `allFit`-class list of fitted models so that users can verify that parameter estimates are stable across optimizers (convergence robustness check).

### Arguments

| Argument | Role |
|---|---|
| `object` | A fitted `merMod`. |
| `meth.tab` | Custom method table with columns `method` and `optimizer`. |
| `data` | Data to bundle with result for debugging. |
| `maxfun` | Maximum function evaluations; auto-translated to each optimizer's convention. |
| `parallel` / `ncpus` / `cl` | Parallel execution options (same as `bootMer`). |
| `catch.errs` | If `TRUE` (default), failed optimizer runs are silently skipped; if `FALSE`, errors propagate. |
| `start_from_mle` | Initialize each refit from the fitted model's MLE. |

### Return value

List of class `allFit` with one entry per optimizer that succeeded. `summary(allFit_result)` extracts `$which.OK`, `$llik`, `$fixef`, `$sdcor`, `$theta`, `$icors` (information criteria) across optimizers.

---

## 8. lmerTest-specific additions in this family

lmerTest adds **no** new simulate, refit, or bootMer methods. Its contribution to this family is:

| Function | What it adds |
|---|---|
| `update.lmerModLmerTest` | Same update logic as `update.merMod` but ensures the returned object is upgraded to `lmerModLmerTest` class so Satterthwaite machinery is available post-update. |
| `as_lmerModLmerTest(x)` | Converts a plain `lmerMod` to `lmerModLmerTest`; called internally by `update.lmerModLmerTest` and `lmerTest::lmer`. Not directly part of the simulate/refit surface but is the recovery path after any refit that loses the lmerTest class. |

lmerTest's `simulate` and `refit` dispatch falls through to lme4's `merMod` methods because `lmerModLmerTest` inherits from `lmerMod` → `merMod`.

---

## 9. Cross-cutting behaviours and idioms

### The canonical bootstrap loop

```r
# lme4 idiomatic parametric bootstrap
boo <- bootMer(fit, FUN = fixef, nsim = 500, seed = 42)
boot::boot.ci(boo, type = "perc")
```

`refit()` is the fast path inside `bootMer`. Users who want custom outer loops call `simulate()` + `refit()` directly:

```r
sims <- simulate(fit, nsim = 200, seed = 1)
results <- lapply(sims, refit, object = fit)
sapply(results, fixef)
```

### NA round-trip

```r
# Original data has NAs; simulate produces a shorter vector with na.action attr
fm0 <- lmer(y ~ x + (1|g), data_with_na, na.action = na.exclude)
ss  <- simulate(fm0, nsim = 1)   # length = nobs(fm0), carries na.action attr
refit(fm0, ss)                    # works: lme4 detects the na.action attr
```

### `evaluate = FALSE` for call inspection

```r
update(fit, . ~ . + z, evaluate = FALSE)
# returns the unevaluated call; useful for checking formula before fitting
```

### Partial `re.form` conditioning

```r
# Condition on subject random effect but simulate new item random effects
simulate(fit, nsim = 10, re.form = ~ (1 | subject))
```

### Simulation from hypothetical parameters (power analysis)

```r
params  <- list(beta = c(1, 0.5), theta = 0.3, sigma = 1.2)
design  <- expand.grid(x = 0:1, id = 1:50)
design$y <- 0
form    <- y ~ x + (1 | id)
sims    <- simulate(form, newdata = design, newparams = params,
                    nsim = 1000, seed = 99)
# sims is a 100×1000 data frame; each column is one power-analysis replicate
```

---

## 10. Summary table

| Function | Primary purpose | Key args not in base R generic |
|---|---|---|
| `simulate.merMod` | Draw nsim responses from fitted model | `use.u`, `re.form`, `newdata`, `newparams`, `cluster.rand`, `allow.new.levels` |
| `simulate.formula` / `.simulateFun` | Draw responses from formula + params (no fit needed) | `formula`, `weights`, `offset`, `cond.sim` |
| `refit.merMod` | Fast refit to new response, same structure | `newweights`, `rename.response`, `maxit` |
| `refitML` | Refit REML model as ML | `optimizer` |
| `update.merMod` | Modify formula/args and refit | `formula.`, `evaluate` |
| `update.lmerModLmerTest` | Same + preserves lmerTest class | same as `update.merMod` |
| `bootMer` | Parametric/semiparametric bootstrap | `FUN`, `type`, `parallel`, `ncpus`, `cl` |
| `allFit` | Convergence check across all optimizers | `meth.tab`, `catch.errs`, `start_from_mle` |
