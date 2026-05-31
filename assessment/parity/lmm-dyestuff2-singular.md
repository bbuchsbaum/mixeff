# Parity Assessment: lmm-dyestuff2-singular

**Cell:** lmm-dyestuff2-singular  
**Date:** 2026-05-31  
**Dataset:** `Dyestuff2` (lme4 built-in; 30 obs, 6 batches)  
**Formula:** `Yield ~ 1 + (1|Batch)`  
**Focus:** near-zero variance / boundary case  
**REML:** TRUE  

## Environment

| Package  | Version |
|----------|---------|
| lme4     | 2.0.1   |
| lmerTest | 3.2.1   |
| mixeff   | 0.1.0   |

---

## 1. Raw Script Output

### lme4 fit

```
boundary (singular) fit: see help('isSingular')
lme4 wall time: 0.027 s

Random effects:
 Groups   Name        Variance Std.Dev.
 Batch    (Intercept)  0.00    0.000
 Residual             13.81    3.716
Number of obs: 30, groups: Batch, 6

Fixed effects:
            Estimate Std. Error      df t value Pr(>|t|)
(Intercept)   5.6656     0.6784 29.0000   8.352 3.32e-09 ***

isSingular: TRUE
theta (Batch.(Intercept)): 0
```

### mixeff fit

```
mixeff wall time: 0.010 s
fit_status: converged_reduced_rank

Variance components:
 group        name variance std_dev correlation       note
 Batch (Intercept)        0       0             [boundary]
[boundary]: variance component is at the boundary of the parameter space.
Residual std. dev.: 3.71568

Fixed effects:
            Estimate Std. Error df  z value Pr(>|z|)            method
(Intercept)   5.6656   0.678388 NA 8.351562        0 asymptotic_wald_z

Inference status:
        term            method    status reliability reliability_reason
 (Intercept) asymptotic_wald_z available         low      not_available
```

---

## 2. Numeric Comparison Table

Tolerances: fixef 1e-4, theta 1e-3, logLik 1e-3, sigma 1e-4.

| Quantity          | lme4 value      | mixeff value    | Max abs diff | Tol  | Status     |
|-------------------|-----------------|-----------------|-------------|------|------------|
| fixef (Intercept) | 5.6656          | 5.6656          | 0.00e+00    | 1e-4 | WITHIN-TOL |
| SE (Intercept)    | 0.678388        | 0.678388        | 2.22e-16    | 1e-4 | WITHIN-TOL |
| vcov[1,1]         | 0.46021         | 0.46021         | 3.33e-16    | 1e-4 | WITHIN-TOL |
| theta             | 0               | 0               | 0.00e+00    | 1e-3 | WITHIN-TOL |
| sigma             | 3.715684        | 3.715684        | 1.33e-15    | 1e-4 | WITHIN-TOL |
| logLik            | -80.914139      | -80.914139      | 1.42e-14    | 1e-3 | WITHIN-TOL |
| AIC               | 167.828278      | 167.828278      | 2.84e-14    | 1e-3 | WITHIN-TOL |
| BIC               | 172.031870      | 172.031870      | 2.84e-14    | 1e-3 | WITHIN-TOL |
| VarCorr Batch var | 0               | 0               | 0.00e+00    | —    | WITHIN-TOL* |
| ranef (all 6)     | 0,0,0,0,0,0     | 0,0,0,0,0,0     | 0.00e+00    | 0.05 | WITHIN-TOL |
| fitted (max diff) | —               | —               | 0.000000    | —    | WITHIN-TOL |

*The probe script's extraction of `VarCorr Batch` initially showed `MISSING` due to `as.numeric()` being applied to a data.frame (mixeff's `VarCorr` returns a named list with a `$table` data.frame, not a numeric matrix). Manual inspection confirms `vc$table$variance == 0` — numerically correct. This is a probe-script extraction issue, not a mixeff defect.

---

## 3. Boundary / Singular Behaviour

| Aspect | lme4 | mixeff |
|--------|------|--------|
| Boundary detection | `isSingular = TRUE`; message: "boundary (singular) fit: see help('isSingular')" | `fit_status = "converged_reduced_rank"`; prints `[boundary]` annotation in VarCorr table; `reliability = "low"` in inference table |
| theta at boundary | 0 | 0 |
| Inference at boundary | lmerTest provides df=29, t=8.352, p=3.32e-09 (Satterthwaite) | mixeff falls back to `asymptotic_wald_z`; df=NA; reliability="low"; `reliability_reason="not_available"` |

### Assessment of boundary handling

mixeff's behaviour is honest and by design. At a singular boundary:

- It correctly identifies `fit_status = "converged_reduced_rank"`.
- It annotates the boundary variance component in the summary with `[boundary]`.
- It downgrades inference reliability to `"low"` and uses `asymptotic_wald_z` (z-test, df=NA) rather than Satterthwaite df.

lme4/lmerTest reports df=29 (Satterthwaite) at the boundary. This is technically the correct Satterthwaite answer when the random effect variance is pinned at zero. mixeff chooses not to provide a finite df at the boundary and flags reliability as low — this is a conservative, honest refusal rather than a capability gap. Per PRD §3 (non-goals), mixeff is an audit-first system: every inference claim must trace to a reliable artifact. A boundary singular fit makes the Satterthwaite df unreliable, so degrading reliability is the correct behaviour.

**Classification:** `works` (by design; the typed diagnostic is clear)

---

## 4. Speed

| Engine  | Wall time |
|---------|-----------|
| lme4    | 0.027 s   |
| mixeff  | 0.010 s   |
| ratio   | 0.37× (mixeff is ~2.7× faster on this tiny dataset) |

Note: Dyestuff2 is very small (30 obs, 6 groups); startup costs dominate. The ratio is not meaningful at this scale. For a speed cell, a scaled dataset with repeated runs would be needed.

---

## 5. Findings Summary

**All numeric quantities are within tolerance** (differences at floating-point machine epsilon or zero). Specifically:

- fixef, SE, vcov: exact match (differences ≤ 2.22e-16, tol 1e-4)
- theta: exact match (both 0)
- sigma: exact match (diff 1.33e-15, tol 1e-4)
- logLik: exact match (diff 1.42e-14, tol 1e-3)
- AIC, BIC: exact match (diff 2.84e-14, tol 1e-3)
- ranef: all zeros, exact match
- fitted: all equal to grand mean (5.6656), exact match

**One cosmetic/classification issue:**

The `VarCorr` extractor surface: mixeff returns a structured `mm_varcorr` list (with `$table` data.frame and `$residual_sd`) rather than a numeric matrix like lme4's `VarCorr`. Callers using `as.data.frame(VarCorr(lme4_fit))$vcov` will need a different extraction idiom for mixeff. This is a **cosmetic/API-surface difference**, not a numeric discrepancy, and is consistent with mixeff's design (structured, serializable artifacts over bare numerics).

**Inference degradation at boundary:**

mixeff falls back to `asymptotic_wald_z` with `reliability="low"` instead of Satterthwaite df=29. The point estimate and SE are identical, but df=NA in mixeff vs df=29 in lme4. This is a deliberate, typed diagnostic (not a silent failure) and aligns with PRD §3 (audit-first, no fabricated inference).

**Severity of findings:**

| Finding | Severity |
|---------|----------|
| All numeric quantities within tolerance | none |
| VarCorr API shape differs from lme4 | cosmetic |
| Boundary inference: df=NA + reliability=low vs lme4 df=29 | cosmetic (by design; honest refusal) |

**Overall outcome: `within-tol`** — all parity quantities match within tolerance; boundary behaviour is honest and by design.

---

## 6. Script

Script: `/Users/bbuchsbaum/code/mixeff/assessment/parity/lmm-dyestuff2-singular-probe.R`
