# Error-message quality probe: single-obs-per-group

**Scenario:** Every group has exactly one observation — the classic degenerate
random-effects case where variance components are unidentifiable.

**Date:** 2026-05-31  
**mixeff commit:** d030be6  
**Probe script:** `/tmp/probe-single-obs-v3.R`

---

## 1. Data setup

```r
set.seed(42)
n_groups <- 20
df_single <- data.frame(
  y       = rnorm(n_groups),
  y_bin   = rbinom(n_groups, 1, 0.5),
  x       = rnorm(n_groups),
  subject = factor(seq_len(n_groups))   # 20 groups, 1 obs each
)
```

---

## 2. lme4 / glmer exact messages

### 2a. `lme4::lmer` — **hard error (good)**

```
ERROR class: simpleError, error, condition
ERROR message: number of levels of each grouping factor must be < number of observations (problems: subject)
```

lme4 refuses to fit.  The message is terse but unambiguous: it names the
offending grouping factor (`subject`) and states the violated constraint.

### 2b. `lme4::glmer` (binomial) — **silent fit (concerning)**

```
Result class: glmerMod
isSingular: FALSE
theta: c("subject.(Intercept)" = 0.152072...)
VarCorr:
      grp        var1 var2       vcov     sdcor
1 subject (Intercept) <NA> 0.02312593 0.1520721
```

glmer returns a fit without any error or warning — the random-effect variance
shrinks toward zero (sd ≈ 0.15) but `isSingular` is FALSE.  This is lme4's own
silent-wrong-answer behaviour for GLMM; it is not a useful reference point.

---

## 3. mixeff exact messages

### 3a. `mixeff::lmm` random-intercept — **silent fit, no diagnostic**

```
Result class: mm_lmm, mm_fit, mm_compiled
theta: 1
boundary: NULL
VarCorr: subject (Intercept)  variance=0.8298  std_dev=0.9109  boundary=FALSE
logLik: -33.056
nobs: 20
fit_status: converged_interior
is_singular: FALSE
```

No error. No warning. No condition of any kind.

The `summary()` output:

```
Linear mixed model fit by REML
Formula: y ~ x + (1 | subject)
Fit status: converged_interior

Variance components:
   group        name variance  std_dev correlation note
 subject (Intercept) 0.829753 0.910908

Fixed effects:
              Estimate Std. Error df    z value  Pr(>|z|)            method
(Intercept)  0.1162800  0.2937491 NA  0.3958481 0.6922171 asymptotic_wald_z
x           -0.3938359  0.2996993 NA -1.3141033 0.1888114 asymptotic_wald_z

Inference status:
        term            method    status reliability reliability_reason
 (Intercept) asymptotic_wald_z available         low      not_available
           x asymptotic_wald_z available         low      not_available
```

Observations:
- `df = NA` for all fixed effects — degrees of freedom cannot be computed.
- `reliability = low` with `reliability_reason = not_available` — the package
  detects something is wrong but does not say *what*.
- `boundary = FALSE` and `is_singular = FALSE` despite the problem being
  structural unidentifiability.
- The residual SD equals the random-effect SD (both 0.9109), which is the
  hallmark of a degenerate fit where sigma absorbs all variance.

### 3b. `mixeff::lmm` random-slope — **silent fit, no diagnostic**

```
Result class: mm_lmm, mm_fit, mm_compiled
theta: c(5.226, 0.644, 1.036)
boundary: NULL
VarCorr: subject (Intercept) variance=1.586; subject x variance=0.086; corr=+0.53
```

Even a random-slope model with `n_groups` = `n_obs` (far fewer observations
than random-effect parameters) fits silently.

### 3c. `mixeff::glmm` (binomial) — **silent fit, inference not computed**

```
Result class: mm_glmm, mm_fit, mm_compiled
theta: 0.108016
VarCorr: subject (Intercept) variance=0.01167  std_dev=0.108016  boundary=FALSE
is_singular: FALSE
optimizer_status: FTOL_REACHED

Fixed effects summary:
(Intercept)  -0.2166  SE=0.547  statistic=NA  p.value=NA  method=not_computed
x            -1.3781  SE=0.709  statistic=NA  p.value=NA  method=not_computed
```

Inference is silently not computed (method = `not_computed`, p-values = NA) but
no condition is raised and no explanation is given for why.

---

## 4. Comparison table

| Aspect | lme4::lmer | lme4::glmer | mixeff::lmm | mixeff::glmm |
|---|---|---|---|---|
| Raises error/warning? | YES — hard error | NO (silent) | NO (silent) | NO (silent) |
| Names the problem? | Yes ("# levels must be < n_obs") | n/a | No | No |
| Returns a fit? | No | Yes (questionable) | Yes | Yes |
| `boundary` flag? | n/a | isSingular=FALSE | FALSE | FALSE |
| Inference suppressed? | n/a | No | Low-reliability (no reason) | not_computed (no reason) |
| Typed mm_condition raised? | n/a | n/a | None | None |

---

## 5. Assessment

### Is this a bug?

**Yes — classified as `bug` (silent-wrong-answer).**

The core promise of mixeff is "less inscrutable errors than lme4."  For
`lmm()`, lme4 hard-errors with a named, actionable message.  mixeff silently
returns a fit that is structurally meaningless: the model has as many groups as
observations, so all variance is confounded between the random effect and
residual error.  The residual SD equalling the random-effect SD is the
numerical fingerprint of this degeneracy, but no diagnostic is surfaced.

The `reliability_reason = not_available` in the inference table is a symptom
(it cannot compute df) but does not identify the cause (degenerate data
structure).  `boundary = FALSE` and `is_singular = FALSE` are actively
misleading: the model is not boundary-singular in the theta-on-boundary sense,
but it is structurally unidentifiable, which is worse.

### What a correct mixeff response would look like

Per PRD §8.1 and the `mm_not_identifiable` / `mm_data_error` condition
hierarchy, this case should raise a typed `mm_data_error` (or
`mm_not_identifiable`) **before** entering the Rust engine, at the
`compile_model()` / data-validation stage.  The message should:

1. State the violated structural constraint (n_groups ≥ n_obs for the grouping
   factor).
2. Name the offending grouping factor.
3. Explain *why* this makes the model unidentifiable (variance components
   cannot be separated from residual error when each group has one observation).
4. Be a typed condition (`mm_not_identifiable` or `mm_data_error`) so
   `tryCatch(..., mm_condition = ...)` works.

Example of an acceptable message:

```
Error in lmm() : grouping factor 'subject' has 20 levels but the data contain
only 20 observations — one observation per group makes random-effect variance
unidentifiable (cannot separate group variance from residual error).
[mm_not_identifiable]
```

### Scope classification

**in-scope-missing** — This is a data-validation check that belongs in the
pre-fit audit path (`compile_model()`), is consistent with the audit-first
design, and is not listed in PRD §3 non-goals.  The related condition class
(`mm_not_identifiable`) already exists in the conditions registry.

---

## 6. Probe script location

`/tmp/probe-single-obs-v3.R` (not committed; reproduce with the data setup
above).
