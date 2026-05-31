# lme4 / lmerTest / pbkrtest — Inference: confint / anova / drop1 / profile

**Survey date:** 2026-05-31  
**Packages surveyed:** lme4 2.0.1, lmerTest 3.2.1, pbkrtest (installed)  
**Purpose:** Exhaustive reference catalogue of the public inference surface for this
capability family. This document records what lme4 and its inference companions offer;
mixeff parity assessment is a separate step.

---

## 1. `confint()` — Confidence Intervals

### 1.1 `confint.merMod` (lme4)

```r
confint(object, parm, level = 0.95,
        method = c("profile", "Wald", "boot"), zeta,
        nsim = 500,
        boot.type = c("perc", "basic", "norm"),
        FUN = NULL, quiet = FALSE,
        oldNames, signames = TRUE,
        boot.scale = c("sdcor", "vcov"), ...)
```

**What it does:** Computes CIs for all or a subset of parameters in an `[ng]lmer` fit.
The method argument selects the statistical approach:

| `method` | What it computes | Coverage characteristics |
|----------|-----------------|--------------------------|
| `"profile"` (default) | Likelihood profile CIs: finds the parameter value where the profile log-likelihood drops by the chi-square cutoff corresponding to `level`. Calls `profile.merMod` internally. | Gold standard; respects non-normality of VC sampling distributions near boundary. |
| `"Wald"` | Normal-theory CI: `estimate ± z_{α/2} * SE`. Valid only for fixed-effect (beta) parameters; returns `NA` for all variance-covariance parameters. | Fast; assumes large-sample normality. |
| `"boot"` | Parametric bootstrap: simulates `nsim` data sets from the fit, refits, and applies `boot.ci`. | Computationally expensive; more robust than Wald. |

**Key arguments:**

- `parm` — integer positions, character names, `"theta_"` (VC parameters only), or
  `"beta_"` (fixed effects only). Determines which parameters are profiled or
  bootstrapped.
- `level` — confidence level (default 0.95).
- `zeta` — manual likelihood cutoff override (profile method only).
- `nsim` — bootstrap replicates (boot method only, default 500).
- `boot.type` — bootstrap CI flavour: `"perc"` (percentile), `"basic"`, or `"norm"`
  (normal approximation). `"stud"` (studentised) and `"bca"` are unavailable because
  they require extra calculations.
- `FUN` — user-supplied bootstrap function; default returns fixed effects + VC
  parameters on SD/correlation scale.
- `quiet` — suppress profiling progress messages.
- `signames` — if `FALSE`, use descriptive names like `sd_(Intercept)|Subject` instead
  of `.sig01` etc.
- `boot.scale` — scale for bootstrap CIs: `"sdcor"` (standard deviations and
  correlations, default) or `"vcov"` (variances).

**Return value:** Numeric matrix with row names = parameter names, column names =
probability bounds (e.g. `2.5 %` and `97.5 %`). CIs are on the standard-deviation
scale for VC parameters.

**Typical workflow note:** Because profiling is expensive, users normally call
`pp <- profile(object)` once and then reuse `confint(pp, level=0.90)`,
`confint(pp, level=0.99)`, etc.

### 1.2 `confint.thpr` (lme4) — CI from a pre-computed profile

```r
confint(object, parm, level = 0.95, zeta, non.mono.tol = 1e-2, ...)
```

Takes an already-computed profile object (`thpr` class, output of `profile.merMod`)
and extracts the CIs at a specified level without re-profiling. `non.mono.tol` controls
tolerance for detecting and warning about non-monotonic profiles before falling back to
linear interpolation.

---

## 2. `profile()` — Likelihood Profile

### 2.1 `profile.merMod` (lme4)

```r
profile(fitted, which = NULL, alphamax = 0.01,
        maxpts = 100, delta = NULL, delta.cutoff = 1/8,
        verbose = 0, devtol = 1e-9, devmatchtol = 1e-5,
        maxmult = 10, startmethod = "prev",
        optimizer = NULL, control = NULL, signames = TRUE,
        parallel = c("no", "multicore", "snow"),
        ncpus = getOption("profile.ncpus", 1L), cl = NULL,
        prof.scale = c("sdcor", "varcov"), ...)
```

**What it does:** Computes a likelihood profile for each parameter in `which` (default:
all) by stepping the parameter away from its MLE in both directions and recording the
deviance at each step. Returns a `thpr` object (a list of data frames, one per
parameter).

**Key arguments:**

- `which` — parameters to profile. `NULL` = all. Integer indexing follows: (1) theta
  (VC) parameters ordered as the lower-triangle of the VC matrix (SDs on diagonal,
  correlations off-diagonal), (2) residual SD / scale, (3) fixed-effect (beta)
  parameters. Character values `"beta_"` or `"theta_"` select entire groups; individual
  names like `".sigma"` or `"(Intercept)"` select single parameters.
- `alphamax` — controls the profiling range: the profile is extended until the profile
  deviance exceeds the 1-alpha cutoff for `alpha = alphamax`.
- `maxpts` — maximum profile evaluation points per parameter per direction (default 100;
  total ≤ 200 per parameter).
- `delta` / `delta.cutoff` — stepping scale. `delta.cutoff = 1/n` gives approximately
  2n points per parameter.
- `parallel` / `ncpus` / `cl` — parallelisation across parameters via `multicore` or
  `snow`.
- `prof.scale` — whether VC parameters are reported on `"sdcor"` (SD + correlation,
  default) or `"varcov"` (variance/covariance) scale.
- `signames` — same as in `confint.merMod`.

**Profile transformation utilities:**

| Function | What it does |
|----------|--------------|
| `as.data.frame(x)` | Converts `thpr` to a data frame for custom plotting. |
| `log(x)` / `logProf(x, ranef=, sigIni=)` | Transforms the profile to use log(SD) for random-effect SD parameters. Useful for visualising skewed profiles near boundary. |
| `varianceProf(x, ranef=)` | Converts from SD scale to variance scale. |
| `xyplot.thpr(x, levels=, conf=, absVal=)` | Lattice plot of profile zeta (signed square-root-deviance) traces for each parameter. Standard confidence levels can be overlaid. |
| `splom.thpr(x, levels=, conf=, draw.lower=, draw.upper=)` | Lattice scatter-plot matrix of pairwise profile joint contours. |

---

## 3. `anova()` — ANOVA Tables and Model Comparison

### 3.1 `anova.merMod` (lme4) — single-model and multi-model

```r
anova(object, ..., refit = TRUE, model.names = NULL)
```

**What it does (single-model call):** Returns a sequential (Type I) ANOVA table of
deviance differences for fixed-effect terms. No p-values are attached; only chi-square
statistics and degrees of freedom. Users should treat these as approximate.

**What it does (multi-model call, e.g. `anova(m1, m2)`):** Compares nested models via
likelihood-ratio tests. Produces a table with log-likelihoods, AIC, BIC, deviance
difference, df difference, and chi-square p-values.

- `refit = TRUE` (default) — lmerMod fits are refitted from REML to ML before
  comparing. Prevents the common mistake of comparing REML models with different fixed
  effects (whose REML log-likelihoods are not comparable).
- `model.names` — optional character vector to label models in the table.

**Limitations:** No denominator df, no F statistics, no p-values for individual terms
in the single-model call. Users needing p-values must use lmerTest's override (§3.2).

### 3.2 `anova.lmerModLmerTest` (lmerTest) — F-test ANOVA with denominator df

```r
anova(object, ...,
      type = c("III", "II", "I", "3", "2", "1"),
      ddf = c("Satterthwaite", "Kenward-Roger", "lme4"))
```

**What it does:** When lmerTest is loaded, `lmer()` returns an `lmerModLmerTest`
object and `anova()` dispatches here instead of to `anova.merMod`. Adds F statistics
and p-values to the ANOVA table using approximate denominator degrees-of-freedom.

- `type` — sum-of-squares decomposition type:
  - `"I"` — sequential (Type I): each term adjusted for all terms above it in the model
    formula.
  - `"II"` — marginal (Type II): each term adjusted for all other main effects of equal
    or lower order; preferred when interaction terms are absent.
  - `"III"` (default) — marginal, adjusted for all other terms including interactions;
    matches SAS PROC MIXED default. Required for balanced interpretation when
    interactions are present.
- `ddf` — denominator df method:
  - `"Satterthwaite"` (default) — Welch-Satterthwaite approximation to the denominator
    df; computationally cheap, good performance in practice.
  - `"Kenward-Roger"` — calls `pbkrtest::KRmodcomp` internally; uses the
    Kenward-Roger bias-corrected covariance matrix and a more precise df approximation;
    better in small samples but requires REML fit.
  - `"lme4"` — falls back to `anova.merMod` output (no p-values, ignores `type`).

**Multi-model call:** When two or more models are passed, `type` and `ddf` are ignored
and the call delegates to sequential LRT comparison (same as `anova.merMod`).

---

## 4. `drop1()` — Single-Term Deletion

### 4.1 `drop1.merMod` (lme4)

```r
drop1(object, scope, scale = 0,
      test = c("none", "Chisq", "user"),
      k = 2, trace = FALSE, sumFun, ...)
```

**What it does:** Fits a separate reduced model for each fixed-effect term in `scope`
(default: all removable terms), compares each to the full model, and returns an ANOVA-
like table. More principled than the sequential Type I ANOVA table from `anova()`.

- `scope` — formula giving terms to consider; default uses `drop.scope(object)`.
- `test`:
  - `"none"` — reports AIC only.
  - `"Chisq"` — appends an asymptotic LRT chi-square statistic and p-value.
  - `"user"` — uses a caller-supplied `sumFun` function; allows plugging in
    `KRmodcomp` or `PBmodcomp` per term (see `sumFun` argument).
- `k` — AIC penalty constant (default 2; set to `log(n)` for BIC-like penalty).
- `sumFun` — for `test="user"`: a function `(object, objectDrop, ...)` returning a
  named numeric vector with an `"method"` attribute. Used to plug in Kenward-Roger or
  parametric-bootstrap tests per term.

**Formula environment caveat:** `drop1` must be able to locate the data within the
model formula environment. Creating formulas in separate functions can break this;
workaround is to set `environment(formula)` explicitly.

### 4.2 `drop1.lmerModLmerTest` (lmerTest)

```r
drop1(object, scope,
      ddf = c("Satterthwaite", "Kenward-Roger", "lme4"),
      force_get_contrasts = FALSE, ...)
```

**What it does:** lmerTest's override of `drop1` for `lmerModLmerTest` objects. Returns
an F-test table with Satterthwaite or Kenward-Roger denominator df and p-values for
each dropped term, rather than the asymptotic chi-square test lme4 provides.

- `ddf` — same three options as `anova.lmerModLmerTest` (§3.2).
- `force_get_contrasts` — internal flag for contrast extraction.

---

## 5. `KRmodcomp` (pbkrtest) — Kenward-Roger F-test

```r
KRmodcomp(largeModel, smallModel, betaH = 0, details = 0)
```

**What it does:** Computes an approximate F-test using the Kenward-Roger (1997)
adjustment for small-sample bias in the fixed-effect covariance matrix. Applicable only
to linear mixed models; models must have the same covariance structure.

- `largeModel` — the full lmer model (must be REML; re-fitted automatically if ML).
- `smallModel` — one of:
  1. A nested lmer model with the same VC structure.
  2. A restriction matrix `L` (k × p) specifying `L·β = L·βH`.
  3. A formula or text string to derive the reduced model from `largeModel`.
- `betaH` — the hypothesis value vector (`L·β = L·βH`); default 0.
- `details` — verbosity level.

**Return value:** An object with a `$stats` data frame containing: `ndf`, `ddf`
(Kenward-Roger denominator df), `Fstat`, `p.value`, `F.scaling`.

**Integration with lmerTest:** `anova(..., ddf="Kenward-Roger")` calls `KRmodcomp`
internally for each term.

**Companion function `getKR()`:** Extracts components (F statistic, ddf, p-value, etc.)
from a `KRmodcomp` result object.

---

## 6. `PBmodcomp` (pbkrtest) — Parametric Bootstrap LRT

```r
PBmodcomp(largeModel, smallModel,
          nsim = 1000, ref = NULL, seed = NULL, cl = NULL, details = 0)
```

**What it does:** Compares nested mixed models via a parametric bootstrap. Simulates
`nsim` data sets from the fitted small model, refits both models to each, and uses the
distribution of the LRT statistic as the reference distribution.

- Models are automatically refitted to ML if fitted with REML.
- `ref` — pre-computed reference distribution (output of `PBrefdist()`); reuse to avoid
  re-simulation.
- `seed` — RNG seed passed to simulation.
- `cl` — cluster object for parallel computation.

**Return value:** An object with a `$test` data frame containing multiple p-value rows:

| Row | Method | What it computes |
|-----|--------|-----------------|
| `LRT` | Chi-square | Assumes LRT ~ chi-square (same as `anova()`) |
| `PBtest` | Parametric bootstrap | Fraction of simulated LRT values ≥ observed |
| `Bartlett` | Bartlett correction | Scales observed LRT by mean of simulated distribution |
| `Gamma` | Gamma reference | Fits a gamma distribution to simulated LRT values |
| `F` | F approximation | Treats LRT/df as F-distributed; denominator df from moment matching |

**Companion `seqPBmodcomp()`:** Sequential parametric bootstrap — stops early once `h`
extreme cases are accumulated; more efficient when p-values are known to be large.

**Integration with `drop1`:** Can be plugged in via `test="user"` + a `sumFun` wrapper
(shown in lme4 `drop1` examples).

---

## 7. `SATmodcomp` (pbkrtest) — Satterthwaite Model Comparison

```r
SATmodcomp(largeModel, smallModel, betaH = 0, details = 0,
           eps = sqrt(.Machine$double.eps))
```

**What it does:** Computes a Satterthwaite-approximated F-test for comparing nested
lmer models. Lighter than Kenward-Roger (no bias correction to the covariance matrix)
but shares the same general interface.

---

## 8. `X2modcomp` / `x2_modcomp` (pbkrtest) — Chi-square Model Comparison

```r
X2modcomp(largeModel, smallModel, betaH = 0, details = 0, ...)
```

**What it does:** Asymptotic chi-square (LRT) comparison between nested models;
essentially the same result as `anova()` but in the pbkrtest framework, allowing
consistent output format alongside `KRmodcomp` / `PBmodcomp`.

---

## 9. `ranova` (lmerTest) — Random-Effect ANOVA

```r
ranova(model, reduce.terms = TRUE, ...)
```

**What it does:** Tests each random-effect term in the model by dropping it (or
reducing it) and computing a likelihood-ratio test with the correct boundary correction.
Returns a table analogous to `anova()` but for random-effect variance components.

- `reduce.terms = TRUE` — correlated random slopes `(x | g)` are first reduced to
  uncorrelated `(1 | g) + (0 + x | g)` before being dropped entirely, yielding two
  separate tests; `FALSE` tests full term removal only.
- Uses Self-Liang (1987) boundary-mixture chi-square distributions when appropriate.
- Analogous to lme4's `anova(m1, m0)` for VC testing but with proper boundary
  correction built in.

---

## 10. `contest` / `contest1D` / `contestMD` (lmerTest) — Contrast Tests

```r
contest(model, L, ...)
contest1D(model, L, ...)
contestMD(model, L, ...)
```

**What they do:** Test user-specified linear contrasts in fixed effects with
Satterthwaite or Kenward-Roger denominator df.

- `contest1D` — one-dimensional (single df) contrast: returns t-statistic with
  Satterthwaite df approximation.
- `contestMD` — multi-dimensional (k df) contrast matrix: returns F-statistic.
- `contest` — dispatches to `contest1D` or `contestMD` based on rank of `L`.
- `L` — contrast vector (for `contest1D`) or matrix (for `contestMD`); must conform to
  the fixed-effects parameter vector.

---

## 11. `show_tests` (lmerTest) — Inspect Test Contrasts

```r
show_tests(object, ...)
```

**What it does:** For an `lmerModLmerTest` or `anova` result, shows the contrast
matrices `L` that lmerTest constructs internally for each fixed-effect term in the ANOVA
table. Allows users to verify exactly which linear hypotheses are being tested for each
term (critical when Type II or Type III decompositions are non-trivial).

---

## 12. `bootMer` (lme4) — Parametric Bootstrap Engine

```r
bootMer(x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA,
        type = c("parametric", "semiparametric"), verbose = FALSE,
        .progress = "none", PBargs = list(),
        parallel = c("no", "multicore", "snow"),
        ncpus = getOption("boot.ncpus", 1L), cl = NULL)
```

**What it does:** The lower-level bootstrap engine used by `confint(method="boot")`.
Simulates `nsim` response vectors from the fitted model and evaluates `FUN` on each
refitted model. Returns a `boot` object (from the `boot` package) for downstream CI
computation.

- `FUN` — a function from a fitted model to a numeric vector; default in `confint` is
  the combined fixed-effects + SD/correlation VC parameters.
- `use.u = FALSE` — simulate new random effects (parametric bootstrap); `TRUE` fixes
  the estimated random effects (conditional bootstrap / semiparametric).
- `re.form` — formula for which random effects to condition on during simulation.
- `type` — `"parametric"` (default) simulates from the conditional model;
  `"semiparametric"` resamples residuals.
- `parallel` / `ncpus` / `cl` — parallelisation options.

---

## 13. `vcovAdj` (pbkrtest) — KR-Adjusted Covariance Matrix

```r
vcovAdj(object, details = 0)
```

**What it does:** Computes the Kenward-Roger bias-corrected (inflated) covariance matrix
for fixed effects. This is the building block consumed by `KRmodcomp` and lmerTest's
`ddf="Kenward-Roger"`. Exposed publicly so users can extract the adjusted SEs directly.

---

## 14. `PBrefdist` (pbkrtest) — Pre-compute Bootstrap Reference Distribution

```r
PBrefdist(largeModel, smallModel, nsim = 1000, seed = NULL, cl = NULL, details = 0)
```

**What it does:** Generates the parametric bootstrap reference distribution (vector of
simulated LRT values under the null) without computing p-values. The result can be
passed as `ref` to `PBmodcomp()` to avoid re-simulation when trying multiple
comparisons against the same null.

---

## 15. Summary Table — All User-Facing Calls

| Function | Package | Purpose | Parameters of note |
|----------|---------|---------|-------------------|
| `confint(method="profile")` | lme4 | Profile-likelihood CIs for all/selected parameters | `parm`, `level`, `zeta`, `signames`, `prof.scale` |
| `confint(method="Wald")` | lme4 | Normal-theory CIs for fixed effects only | `parm`, `level` |
| `confint(method="boot")` | lme4 | Bootstrap CIs for all/selected parameters | `parm`, `nsim`, `boot.type`, `FUN`, `boot.scale` |
| `profile()` | lme4 | Likelihood profile object (prerequisite for reuse) | `which`, `alphamax`, `maxpts`, `parallel`, `prof.scale` |
| `confint.thpr()` | lme4 | CIs from pre-computed profile | `level`, `zeta`, `non.mono.tol` |
| `log()` / `logProf()` | lme4 | Log-transform profile for plotting | `ranef`, `sigIni` |
| `varianceProf()` | lme4 | Variance-scale profile | `ranef` |
| `xyplot.thpr()` | lme4 | Lattice plot of profile traces | `levels`, `conf`, `absVal` |
| `splom.thpr()` | lme4 | Lattice scatter-plot matrix of joint profile contours | `levels`, `conf` |
| `anova()` (single-model) | lme4 | Sequential deviance table; no p-values | `refit`, `model.names` |
| `anova()` (multi-model) | lme4 | LRT model comparison table | `refit`, `model.names` |
| `anova(ddf="Satterthwaite")` | lmerTest | F-test ANOVA with Satterthwaite df | `type`, `ddf` |
| `anova(ddf="Kenward-Roger")` | lmerTest | F-test ANOVA with KR df | `type`, `ddf` |
| `drop1(test="none")` | lme4 | AIC-only single-term deletion table | `scope`, `k` |
| `drop1(test="Chisq")` | lme4 | LRT chi-square per dropped term | `scope` |
| `drop1(test="user", sumFun=)` | lme4 | Custom test per dropped term (KR or PB hookable) | `sumFun` |
| `drop1(ddf="Satterthwaite")` | lmerTest | F-test drop1 with Satterthwaite df | `ddf` |
| `drop1(ddf="Kenward-Roger")` | lmerTest | F-test drop1 with KR df | `ddf` |
| `KRmodcomp()` | pbkrtest | KR F-test; large vs. small model or restriction matrix | `largeModel`, `smallModel`, `betaH` |
| `PBmodcomp()` | pbkrtest | Parametric bootstrap LRT comparison | `largeModel`, `smallModel`, `nsim`, `ref`, `seed`, `cl` |
| `seqPBmodcomp()` | pbkrtest | Sequential (early-stop) parametric bootstrap | `h`, `nsim` |
| `SATmodcomp()` | pbkrtest | Satterthwaite F-test; same interface as KRmodcomp | `betaH`, `eps` |
| `X2modcomp()` | pbkrtest | Chi-square LRT in pbkrtest framework | `betaH` |
| `ranova()` | lmerTest | LRT table for random-effect terms with boundary correction | `reduce.terms` |
| `contest()` / `contest1D()` / `contestMD()` | lmerTest | Contrast tests with Satterthwaite df | `L` |
| `show_tests()` | lmerTest | Inspect internal contrast matrices used in anova() | — |
| `bootMer()` | lme4 | Low-level parametric bootstrap engine | `FUN`, `nsim`, `use.u`, `re.form`, `type`, `parallel` |
| `vcovAdj()` | pbkrtest | KR bias-corrected covariance matrix | `details` |
| `PBrefdist()` | pbkrtest | Pre-compute bootstrap null distribution for reuse | `nsim`, `seed`, `cl` |
| `getKR()` | pbkrtest | Extract KR components from KRmodcomp result | — |

---

## 16. Behavioural Details Users Rely On

1. **REML auto-refit.** Both `anova.merMod(refit=TRUE)` and `KRmodcomp`/`PBmodcomp`
   automatically refit REML models to ML (or REML for KR) before comparing. Users
   expect not to think about this.

2. **Profile parameter ordering.** `profile(..., which=)` uses a specific parameter
   index order: theta (VC lower triangle), sigma (residual SD), then beta (fixed
   effects). Users rely on this ordering when selecting `parm` in `confint`.

3. **Boundary parameters.** Profile CIs for VC parameters near zero have non-symmetric,
   possibly truncated distributions. `xyplot.thpr` and `logProf` expose this
   visually. Users near singularity expect the profile to correctly reflect the
   boundary.

4. **Type II vs. Type III.** lmerTest `anova(type="II")` vs `"III"` produce different
   F-values when predictors are correlated or interactions are present. Users must
   choose consciously; lmerTest default is `"III"`.

5. **`drop1` formula environment.** `drop1.merMod` requires the formula environment
   to contain the data. Users who construct formulas programmatically in helper
   functions must set `environment(f) <- parent.frame()`.

6. **Bootstrap failure fraction.** `PBmodcomp` reports the number of usable replicates
   (some simulated LRT values may be negative due to numerical issues). Users inspect
   this to assess reliability.

7. **KR requires linear VC.** `KRmodcomp` is only valid when the model's covariance
   structure is a linear combination of known matrices (standard random-intercept and
   random-slope models satisfy this; non-standard structures may not).

8. **`ranova` boundary correction.** Unlike `anova(m, m0)`, `ranova` uses the
   Self-Liang 50:50 mixture chi-square for variance-component tests at the boundary.
   This halves the naive p-value for standard one-parameter VC tests.

9. **`confint` on stored profile.** `confint(pp, level=0.90)` is fast; recomputing
   from scratch each time is very slow for complex models. The separation of
   `profile()` and `confint()` is intentional and users exploit it heavily.

10. **`show_tests` transparency.** lmerTest's `show_tests()` reveals the exact `L`
    matrix used for each term so users can audit whether the ANOVA table tests the
    hypothesis they intend.
