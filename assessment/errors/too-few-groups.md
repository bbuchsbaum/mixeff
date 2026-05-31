# Error-message quality probe: "too-few-groups"

**Scenario:** random effect grouping factor has only 2 levels.

**Date:** 2026-05-31
**mixeff version:** installed from /Users/bbuchsbaum/code/mixeff
**lme4 version:** 2.0.1

---

## Setup

```r
set.seed(42)
df_lmm <- data.frame(
  y       = rnorm(20),
  x       = rnorm(20),
  subject = factor(rep(c("A", "B"), each = 10))   # 2 levels
)
df_glmm <- data.frame(
  y     = rbinom(20, 1, 0.5),
  x     = rnorm(20),
  group = factor(rep(c("A", "B"), each = 10))      # 2 levels
)
```

---

## lme4 behavior

### lme4::lmer — random intercept `(1 | subject)`, 2 groups

**Result:** fits silently; no warning, no error.

```
[lmer SUCCESS]: fit object returned (class: lmerMod)
  VarCorr: subject intercept variance = 0.04606879
  sigma = 1.202528
```

lme4 returns a fit with no signal whatsoever that only 2 groups is unusual.
The variance estimate (0.046) is almost certainly unreliable, but the user
receives no indication of this.

### lme4::lmer — random slope `(1 + x | subject)`, 2 groups

**Result:** fits with a singular-fit message only.

```
boundary (singular) fit: see help('isSingular')
```

This message (issued as a `message()` not a `warning()`) names singularity but
does not explain *why* the model is singular (too few groups for a 3-parameter
covariance). The user is pointed to `help('isSingular')` for more context.

### lme4::glmer — random intercept `(1 | group)`, binomial, 2 groups

**Result:** fits with a singular-fit message only.

```
boundary (singular) fit: see help('isSingular')
```

Same terse message as the slope case above.

---

## mixeff behavior

### mixeff::lmm — random intercept `(1 | subject)`, 2 groups

**Result:** fit succeeds; structured diagnostic table emitted.

```
[lmm SUCCESS]: fit returned (class: mm_lmm, mm_fit, mm_compiled)
[lmm fit_status]: converged_interior

DIAGNOSTICS:
code                       severity  message
random_effect_few_levels   warning   2 levels are fit-eligible for a scalar random
                                     intercept but below the v0 reliability threshold 5
support_note               info      the requested covariance structure is
                                     information-hungry relative to the observed
                                     grouping levels
scope_note                 info      `x` varies within `subject`, so a `subject`-level
                                     slope is structurally possible
```

The `random_effect_few_levels` diagnostic is also visible at compile time via
`compile_model()` / `explain_model()` **before** fitting:

```
Design notes:
  random_effect_few_levels: 2 levels are fit-eligible for a scalar random intercept
    but below the v0 reliability threshold 5
  support_note: the requested covariance structure is information-hungry relative
    to the observed grouping levels
  scope_note: `x` varies within `subject`, so a `subject`-level slope is
    structurally possible
```

### mixeff::lmm — random slope `(1 + x | subject)`, 2 groups

**Result:** fit completes but `fit_status = not_optimized`; richer diagnostics.

```
[lmm slope fit_status]: not_optimized

DIAGNOSTICS:
code                      severity  message
covariance_too_rich       warning   2 levels are below the v0 full-covariance
                                    threshold 15 for 3 covariance parameters
optimizer_nonconvergence  warning   optimizer stopped before an acceptable
                                    convergence criterion with return code
                                    'MAXEVAL_REACHED'
boundary_parameter        info      standard deviation for x in (1 + x | subject)
                                    is on its lower bound
```

### mixeff::glmm — random intercept `(1 | group)`, binomial, 2 groups

**Result:** fit returns; `fit_status = converged_boundary`.

```
[glmm fit_status]: converged_boundary

DIAGNOSTICS:
code                       severity  message
random_effect_few_levels   warning   2 levels are fit-eligible for a scalar random
                                     intercept but below the v0 reliability threshold 5
support_note               info      the requested covariance structure is
                                     information-hungry relative to the observed
                                     grouping levels
scope_note                 info      `x` varies within `group`, so a `group`-level
                                     slope is structurally possible
boundary_parameter         info      GLMM covariance state classified as
                                     ValidZeroVariance
```

---

## Direct comparison

| Dimension              | lme4                                              | mixeff                                                                              |
|------------------------|---------------------------------------------------|-------------------------------------------------------------------------------------|
| Random intercept, 2 grp | **Silent success** — no warning, no diagnostic  | `random_effect_few_levels` warning at compile AND fit time, names threshold (5)     |
| Random slope, 2 grp    | Terse `boundary (singular) fit: see help(...)` message | `covariance_too_rich` + `optimizer_nonconvergence` + `boundary_parameter`; `fit_status=not_optimized` |
| GLMM, 2 grp            | Same terse `boundary (singular)` message         | Same structured diagnostics as LMM case + `boundary_parameter` info                |
| Pre-fit visibility     | None (lme4 has no pre-fit audit path)            | `compile_model()` + `explain_model()` expose the `few_levels` warning before fitting |
| Actionability          | `isSingular` help page is general                | Message names the *specific count* (2) and the *specific threshold* (5); links to the affected term |
| Condition type         | `message()` (easy to miss in scripts)            | Typed `mm_condition` subclass, catchable by `tryCatch(..., mm_condition=)`          |

---

## Issues found

### Minor: duplicated diagnostic rows

In every case the `random_effect_few_levels` / `support_note` / `scope_note`
rows appear **twice** in the diagnostics table. The rows are identical; the
duplication is cosmetic but noisy. This appears to be an upstream engine issue
(the compiler emits the diagnostic from two code paths — design audit and
artifact-level diagnostics). Classification: **needs-work** (minor cosmetic
bug; does not affect correctness or message clarity).

---

## Verdict

**good** — mixeff is materially clearer than lme4 on this scenario.

For the random-intercept case where lme4 is completely silent, mixeff emits a
named, typed, severity-tagged diagnostic (`random_effect_few_levels: warning`)
that identifies the count (2), the reliability threshold (5), and the affected
grouping factor, and does so at pre-fit compile time as well as post-fit. For
the random-slope case, mixeff provides three distinct typed diagnostics
(`covariance_too_rich`, `optimizer_nonconvergence`, `boundary_parameter`) that
together explain both *why* the model is hard to fit and *what* the optimizer
did, whereas lme4 produces only the generic `boundary (singular) fit` message.

The one blemish — duplicated rows in the diagnostics table — is cosmetic and
does not obscure any information.
