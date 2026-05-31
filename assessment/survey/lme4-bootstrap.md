# lme4 / lmerTest / pbkrtest — Bootstrap Capability Surface

Survey date: 2026-05-31  
Packages confirmed installed: lme4 2.0.1, lmerTest 3.2.1, pbkrtest (version on system), boot (base R)

---

## Overview

lme4 provides a complete model-based bootstrap infrastructure for mixed models. The
primary entry point is `bootMer()`, which generates bootstrap replicates of any scalar
statistic derived from a refitted model. `confint.merMod()` wraps `bootMer()` for the
common case of parameter confidence intervals. The companion package `pbkrtest`
provides `PBmodcomp()` / `PBrefdist()` for parametric bootstrap model comparison (LRT
p-values). `simulate.merMod()` is the lower-level data-generation engine that
`bootMer()` calls internally.

---

## 1. `bootMer()` — core bootstrap engine

**Package:** lme4  
**Signature:**
```r
bootMer(x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA,
        type = c("parametric", "semiparametric"),
        verbose = FALSE, .progress = "none", PBargs = list(),
        parallel = c("no", "multicore", "snow"),
        ncpus = getOption("boot.ncpus", 1L), cl = NULL)
```

### Arguments

| Argument | Default | Purpose |
|---|---|---|
| `x` | — | Fitted `merMod` object (lmer/glmer/nlmer) to bootstrap from |
| `FUN` | — | User-supplied function `f(merMod) → named numeric vector`; the statistic of interest; **this is what gives bootMer its generality** |
| `nsim` | `1` | Number of bootstrap replicates (the bootstrap B / R) |
| `seed` | `NULL` | Passed to `set.seed()` for reproducibility |
| `use.u` | `FALSE` | If `TRUE`, hold spherical random effects u fixed; all inference is conditional on observed u. If `FALSE`, new u are drawn each replicate |
| `re.form` | `NA` | Formula specifying which random effects to incorporate; `NA` ≡ `use.u=FALSE`, `NULL` ≡ `use.u=TRUE`; finer-grained than `use.u` |
| `type` | `"parametric"` | Bootstrap variant: `"parametric"` resamples from fitted parametric distributions; `"semiparametric"` resamples residuals (partially implemented; lmer/glmer only; experimental) |
| `verbose` | `FALSE` | Print iteration progress to console |
| `.progress` | `"none"` | Progress bar type: `"none"`, `"txt"`, `"tk"` (requires tcltk), `"win"` (Windows). Disabled automatically during parallel runs |
| `PBargs` | `list()` | Extra arguments forwarded to the progress bar constructor (e.g., `list(style=3)`) |
| `parallel` | `"no"` | Parallelisation backend: `"no"`, `"multicore"` (fork, not Windows), `"snow"` (socket cluster) |
| `ncpus` | `getOption("boot.ncpus", 1L)` | Number of parallel worker processes |
| `cl` | `NULL` | Pre-existing `parallel`/`snow` cluster object; if not supplied and `parallel="snow"`, a local cluster is created for the duration of the call |

### Bootstrap modes (type × use.u combinations)

| `type` | `use.u` | Behaviour |
|---|---|---|
| `"parametric"` | `FALSE` | Re-draw both u (new normal deviates) and ε; fully marginal parametric bootstrap (the standard/recommended mode) |
| `"parametric"` | `TRUE` | Re-draw only ε (or GLMM response values); u fixed at estimates — conditional parametric bootstrap |
| `"semiparametric"` | `TRUE` | Re-sample response residuals; u fixed. Experimental; warning issued for GLMMs |
| `"semiparametric"` | `FALSE` | **Not implemented** (Morris 2002 argues against resampling from u estimates) |

### Return value

An object of S3 class `"boot"` compatible with the `boot` package's `boot()` output.
Users can therefore apply all of `boot::boot.ci()`, `plot.boot()`, `confint.boot()`,
`as.data.frame()`, etc. directly to the result.

Diagnostic attributes on the returned object:
- `attr(., "bootFail")` — integer count of replicate failures (errors)
- `attr(., "boot.fail.msgs")` — character vector of error messages from failed replicates
- `attr(., "boot.all.msgs")` — all messages, warnings, and errors across replicates

### Typical FUN patterns

```r
# Fixed-effects only
FUN <- function(.) fixef(.)

# Fixed + variance components
FUN <- function(.) {
  s <- sigma(.)
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta")))
}

# Alternative using VarCorr
FUN <- function(.) c(beta = fixef(.), sigma = sigma(.), sig01 = sqrt(unlist(VarCorr(.))))
```

---

## 2. `confint.merMod()` — bootstrap confidence intervals via wrapper

**Package:** lme4  
**Signature:**
```r
## S3 method for class 'merMod'
confint(object, parm, level = 0.95,
        method = c("profile", "Wald", "boot"), zeta,
        nsim = 500,
        boot.type = c("perc", "basic", "norm"),
        FUN = NULL, quiet = FALSE,
        oldNames, signames = TRUE,
        boot.scale = c("sdcor", "vcov"),
        ...)
```

### Arguments relevant to bootstrap

| Argument | Default | Purpose |
|---|---|---|
| `object` | — | Fitted `merMod` model |
| `parm` | — | Parameters to compute CIs for: integer positions, character names, `"theta_"` (variance-covariance), or `"beta_"` (fixed effects) |
| `level` | `0.95` | Confidence level |
| `method` | `"profile"` | Use `"boot"` to invoke `bootMer()` internally |
| `nsim` | `500` | Number of bootstrap replicates (passed to `bootMer()`) |
| `boot.type` | `"perc"` | CI type from `boot::boot.ci()`: `"perc"` (percentile), `"basic"`, `"norm"` (normal approximation). Note: `"stud"` and `"bca"` are **unavailable** because they require extra computed components |
| `FUN` | `NULL` | Custom extraction function passed to `bootMer()`; if `NULL`, the default returns fixed effects and variance-covariance parameters on the sd/correlation scale |
| `quiet` | `FALSE` | Suppress computation-intensive messages |
| `signames` | `TRUE` | Use `.sigNN` abbreviated names for variance-covariance parameters; `FALSE` gives `sd_(Intercept)|Subject` style names |
| `boot.scale` | `"sdcor"` | Scale for bootstrap CIs: `"sdcor"` (standard deviation / correlation) or `"vcov"` (variance-covariance) |
| `...` | — | Additional arguments forwarded to `bootMer()` (including `.progress`, `PBargs`, `parallel`, `ncpus`, `cl`) |

### `confint.thpr()` — on profile objects

```r
## S3 method for class 'thpr'
confint(object, parm, level = 0.95, zeta, non.mono.tol = 1e-2, ...)
```
Computes CIs from a pre-computed likelihood profile object (not bootstrap, but the sibling method users compare against).

---

## 3. `simulate.merMod()` — data generation engine

**Package:** lme4  
**Signature:**
```r
## S3 method for class 'merMod'
simulate(object, nsim = 1, seed = NULL,
         use.u = FALSE, re.form = NA,
         newdata = NULL, newparams = NULL, family = NULL,
         cluster.rand = rnorm,
         allow.new.levels = FALSE, na.action = na.pass, ...)
```

`bootMer()` calls `simulate.merMod()` internally to generate each replicate dataset.
Users can call it directly to generate predictive simulations without refitting.

### Arguments

| Argument | Default | Purpose |
|---|---|---|
| `object` | — | Fitted `merMod` |
| `nsim` | `1` | Number of response vectors to simulate |
| `seed` | `NULL` | Seed for reproducibility |
| `use.u` | `FALSE` | Hold u fixed (see `bootMer()` description) |
| `re.form` | `NA` | Formula controlling which random effects are conditioned on vs. re-drawn |
| `newdata` | `NULL` | New data frame for out-of-sample simulation |
| `newparams` | `NULL` | List with `theta`, `beta`, optionally `sigma` — simulate from non-fitted parameter values |
| `family` | `NULL` | GLM family override |
| `cluster.rand` | `rnorm` | Function generating standardised random cluster effects; swap for non-normal random effects (e.g., heavy-tailed, truncated) |
| `allow.new.levels` | `FALSE` | Allow new grouping levels in `newdata` (uses population-level u=0) |
| `na.action` | `na.pass` | Handling of NAs in new data |

### `.simulateFun()` — internal workhorse

Also exported; extends `simulate.merMod` with additional `formula`, `weights`, `offset`, `cond.sim` arguments for programmatic use. `cond.sim = FALSE` simulates only random effects without drawing from the conditional distribution.

---

## 4. `PBmodcomp()` — parametric bootstrap model comparison

**Package:** pbkrtest  
**Signature:**
```r
PBmodcomp(largeModel, smallModel, nsim = 1000, ref = NULL,
          seed = NULL, cl = NULL, details = 0)
```

Users rely on this when comparing nested mixed models via LRT and need valid p-values
without assuming chi-squared null distributions (which are often incorrect when testing
variance components near the boundary).

### Arguments

| Argument | Default | Purpose |
|---|---|---|
| `largeModel` | — | The full (larger) fitted model. S3 dispatch supports `merMod` (lmer/glmer) and `lm` |
| `smallModel` | — | The restricted (smaller) nested model. Can be a fitted model object, a character string of the term to drop, a formula `~. - term`, or a contrast matrix L |
| `nsim` | `1000` | Bootstrap replicates to form the reference LRT distribution |
| `ref` | `NULL` | Pre-computed reference distribution (numeric vector from `PBrefdist()`); avoids re-simulation when trying different comparisons on the same model pair |
| `seed` | `NULL` | Seed for reproducibility |
| `cl` | `NULL` | Cluster specification: integer (number of local cores) or a `parallel` cluster object |
| `details` | `0` | Verbosity level for debugging output |

### Returned p-value columns

`PBmodcomp()` returns a data frame with one row per test type:

| Column | Description |
|---|---|
| `LRT` | Standard likelihood ratio test p-value (chi-squared assumption) |
| `PBtest` | Fraction of simulated LRT values ≥ observed LRT |
| `Bartlett` | Bartlett-corrected LRT p-value (mean of simulated LRTs used as correction) |
| `Gamma` | P-value from gamma distribution fitted to simulated reference LRTs |
| `F` | F-approximation p-value (LRT / df, denominator df matched by moment) |

### REML handling

If either model is fitted with `REML=TRUE`, `PBmodcomp()` automatically refits with
`REML=FALSE` before computing p-values. Users do not need to worry about this.

### Negative LRT values

Some bootstrap replicates may yield numerically negative LRT statistics (overfitting
artefact). The function reports the count of usable (positive) samples and excludes
the rest.

---

## 5. `PBrefdist()` — pre-compute reference distribution

**Package:** pbkrtest  
**Signature:**
```r
PBrefdist(largeModel, smallModel, nsim = 1000, seed = NULL,
          cl = NULL, details = 0)
```

Separates the expensive simulation step from the test step. Users simulate once,
store the result, and supply it as `ref` to multiple `PBmodcomp()` calls (e.g., to
try different `nsim` values or inspect the distribution). Supports the same `cl`
parallelism as `PBmodcomp()`.

---

## 6. `seqPBmodcomp()` — sequential parametric bootstrap

**Package:** pbkrtest  
**Note:** This function appears in the help page for `PBmodcomp()` but is **not
exported** in the installed version of pbkrtest on this system.

Described interface:
```r
seqPBmodcomp(largeModel, smallModel, h = 20, nsim = 1000, cl = 1)
```
Stops sampling once `h` extreme cases (LRT ≥ observed) have been collected, giving an
adaptive early-stopping procedure.

---

## 7. `refit()` — refit with new response

**Package:** lme4  
**Signature:**
```r
refit(object, newresp = NULL, newoffset = NULL, ...)
```

Refits a `merMod` with a new response vector, reusing the existing model structure and
starting from current parameter estimates. This is the mechanism `bootMer()` uses
internally for each replicate. Users can call it directly for custom bootstrap loops
that manage simulation and refitting separately.

---

## 8. Cross-cutting behaviours users rely on

### Progress bars (`.progress` / `PBargs`)

`bootMer()` and `confint(..., method="boot")` both accept `.progress` (forwarded to
`bootMer()`). Supported values:

| Value | Requires |
|---|---|
| `"none"` | default; no output |
| `"txt"` | `utils::txtProgressBar`; works everywhere |
| `"tk"` | `tcltk::tkProgressBar`; requires tcltk package |
| `"win"` | `utils::winProgressBar`; Windows only |

Progress bars are silently disabled during parallel operation.

`PBargs` is a list forwarded to the progress bar constructor; `list(style=3)` gives a
cleaner style with the package authors' preferred rendering.

### Parallel computation

Both `bootMer()` (via `parallel`, `ncpus`, `cl`) and `PBmodcomp()` / `PBrefdist()`
(via `cl`) support parallel execution. For `"snow"` clusters, users must pre-load lme4
on all workers:

```r
clusterEvalQ(cl, library("lme4"))
```

and export any environment objects needed by `FUN` via `clusterExport()`.

### Seed control

All bootstrap functions accept a `seed` argument passed to `set.seed()` for
reproducible results.

### Integration with the `boot` package

`bootMer()` returns a `"boot"`-class object. Users apply:
- `boot::boot.ci(boo, index=i, type=c("norm","basic","perc"))` — per-parameter CIs with optional transformations (`h`, `hdot`, `hinv`)
- `plot.boot(boo, index=i)` — density/QQ plots of bootstrap distribution
- `confint.boot(boo)` — shorthand CIs across all parameters
- `as.data.frame(boo)` — extract replicate matrix

### Semiparametric bootstrap (partial implementation)

The semiparametric variant (`type="semiparametric"`, `use.u=TRUE`) resamples response
residuals. Caveats:
- Only implemented for lmer / glmer results
- Marked as **experimental**
- A warning is generated for GLMMs (resampled data no longer has the same properties)
- The `use.u=FALSE` + semiparametric combination is not implemented (Morris 2002)

### `cluster.rand` — non-normal random effects simulation

`simulate.merMod()` accepts a `cluster.rand` argument (default `rnorm`) allowing users
to substitute any function generating standardised random cluster effects. This enables
sensitivity analyses for departures from the normality assumption on random effects
(e.g., heavy-tailed, truncated-normal, or mixture distributions).

---

## 9. Summary table of user-facing functions

| Function | Package | What it does | Why users need it |
|---|---|---|---|
| `bootMer(x, FUN, nsim, ...)` | lme4 | General model-based (semi-)parametric bootstrap; refit model nsim times, collect FUN statistics | Core bootstrap engine; handles any scalar statistic |
| `confint(merMod, method="boot", ...)` | lme4 | Parameter CIs via `bootMer()` + `boot::boot.ci()` | Convenient single call for standard parameter CIs |
| `simulate.merMod(object, ...)` | lme4 | Generate response replicates from fitted model | Lower-level data generation; also powers bootMer internally |
| `.simulateFun(...)` | lme4 | Extended internal simulation engine with formula/weights/offset/cond.sim | Programmatic simulation from non-fitted models |
| `refit(object, newresp)` | lme4 | Refit model with new response, reusing structure | Used inside custom bootstrap loops |
| `PBmodcomp(large, small, nsim, ...)` | pbkrtest | Parametric bootstrap LRT p-values for nested model comparison | Valid p-values for variance component tests near boundary |
| `PBrefdist(large, small, nsim, ...)` | pbkrtest | Pre-compute LRT reference distribution separately | Reuse expensive simulation across multiple comparisons |
| `seqPBmodcomp(...)` | pbkrtest | Adaptive sequential bootstrap (stop after h extremes) | Efficient p-value estimation; not exported in current install |
