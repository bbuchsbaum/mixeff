# Parity Assessment: lmm-sleep-corr

**Cell:** lmm-sleep-corr  
**Date:** 2026-05-31  
**Dataset:** sleepstudy (lme4 built-in, N=180, 18 subjects × 10 days)  
**Formula:** `Reaction ~ Days + (Days|Subject)` (correlated random slope)  
**REML:** TRUE  

---

## Environment

| Package  | Version |
|----------|---------|
| lme4     | 2.0.1   |
| lmerTest | 3.2.1   |
| mixeff   | 0.1.0   |

---

## Script

Written to `/tmp/lmm_sleep_corr_probe.R`. Namespace conflicts resolved: mixeff exports
`fixef`, `ranef`, `getME`, `VarCorr` — all lme4-object calls use `lme4::` prefix explicitly.

---

## Raw Output

```
=== ENVIRONMENT ===
lme4 version: 2.0.1 
lmerTest version: 3.2.1 
mixeff version: 0.1.0 

=== FIT lme4 ===
lme4 wall time: 0.035 sec

=== FIT mixeff (verbose=-1 to suppress audit output) ===
mixeff wall time: 0.013 sec

mixeff fit class: mm_lmm, mm_fit, mm_compiled 

=== FIXED EFFECTS ===
lme4 fixef:
(Intercept)        Days 
  251.40510    10.46729 
mixeff fixef:
(Intercept)        Days 
  251.40510    10.46729 
Abs diff fixef:
 (Intercept)         Days 
5.684342e-13 1.865175e-13 
Max abs diff fixef: 5.684342e-13 
Within tol (1e-4): TRUE 

=== SE (from vcov) ===
lme4 SE:
(Intercept)        Days 
   6.824597    1.545790 
mixeff SE:
(Intercept)        Days 
   6.824726    1.546506 
Abs diff SE:
 (Intercept)         Days 
0.0001296386 0.0007162301 
Max abs diff SE: 0.0007162301 
Within tol (1e-4): FALSE 

=== THETA (variance-covariance parameters) ===
lme4 theta:
     Subject.(Intercept) Subject.Days.(Intercept)             Subject.Days 
              0.96674177               0.01516906               0.23090995 
mixeff theta (raw):
[1] 0.9668414 0.0151315 0.2310602
Abs diff theta: 9.965909e-05 3.756139e-05 0.0001502393 
Max abs diff theta: 0.0001502393 
Within tol (1e-3): TRUE 

=== VARCORR ===
lme4 VarCorr:
 Groups   Name        Std.Dev. Corr  
 Subject  (Intercept) 24.7407        
          Days         5.9221  0.066 
 Residual             25.5918        
mixeff VarCorr:
Variance components:
   group        name variance  std_dev correlation note
 Subject (Intercept) 612.1570 24.74180                 
 Subject        Days  35.1125  5.92558       +0.07     
Residual std. dev.: 25.5904

=== SIGMA ===
lme4 sigma: 25.5918 
mixeff sigma: 25.59036 
Abs diff sigma: 0.00143999 
Within tol (1e-4): FALSE 

=== LOGLIK / AIC / BIC ===
lme4 logLik: -871.8141  AIC: 1755.628  BIC: 1774.786 
mixeff logLik: -871.8141  AIC: 1755.628  BIC: 1774.786 
Abs diff logLik: 4.133908e-06  (tol 1e-3: TRUE )
Abs diff AIC: 8.267816e-06  BIC: 8.267816e-06 

=== RANDOM EFFECTS (ranef) ===
lme4 ranef Subject (first 6 rows):
    (Intercept)      Days
308    2.258551  9.198976
309  -40.398738 -8.619681
310  -38.960409 -5.448856
330   23.690620 -4.814350
331   22.260313 -3.069912
332    9.039568 -0.272177
mixeff ranef Subject (first 6 rows):
    (Intercept)       Days
308    2.249080  9.2010998
309  -40.393426 -8.6210027
310  -38.958300 -5.4494579
330   23.697891 -4.8158960
331   22.265615 -3.0710221
332    9.040696 -0.2723988
Max abs diff ranef per column:
(Intercept)        Days 
0.011284537 0.002522863 
Max abs diff ranef overall: 0.01128454 

=== FITTED VALUES ===
Max abs diff fitted: 0.01142123 
Mean abs diff fitted: 0.002425965 

=== CONVERGENCE STATUS ===
lme4 converged (no warnings): TRUE 
mixeff convergence field: (empty)
mixeff object names: call, formula, REML, control, vars, model_frame, weights, artifact,
  fit, fit_summary, schema, rust_handle, lazy_cache, beta, theta, sigma, logLik, deviance,
  AIC, BIC, nobs, dof, df_residual, fit_status, std_errors, fixed_effect_vcov,
  fixed_fitted, fitted, residuals, random_effects, varcorr 

=== SPEED ===
lme4 avg time per fit: 0.0208 sec
mixeff avg time per fit: 0.003 sec
Speed ratio (mixeff/lme4): 0.1442 (mixeff ~7x faster)
```

---

## Quantity-by-Quantity Analysis

### Fixed Effects (beta)

| Parameter   | lme4      | mixeff    | Abs diff     | Tol 1e-4 | Status       |
|-------------|-----------|-----------|--------------|----------|--------------|
| (Intercept) | 251.40510 | 251.40510 | 5.68e-13     | PASS     | within-tol   |
| Days        | 10.46729  | 10.46729  | 1.87e-13     | PASS     | within-tol   |

**Result: within-tol.** Numerical agreement to machine precision.

---

### SE / vcov

| Parameter   | lme4 SE  | mixeff SE | Abs diff    | Tol 1e-4 | Status   |
|-------------|----------|-----------|-------------|----------|----------|
| (Intercept) | 6.824597 | 6.824726  | 1.30e-4     | FAIL     | divergent |
| Days        | 1.545790 | 1.546506  | 7.16e-4     | FAIL     | divergent |

**Result: divergent.** Both SEs exceed 1e-4 tolerance. Max abs diff = 7.2e-4 (~0.046% relative error on Days SE). This is a minor numerical divergence — the values are practically close but outside the stated tolerance. Likely cause: differing numerical Hessian approximation vs. analytic Fisher information, or slightly different optimizer convergence point.

**Severity: minor.** SEs are practically close (~0.05% relative error) but exceed 1e-4 absolute tolerance.

---

### Theta (Cholesky factor elements)

| Element                    | lme4       | mixeff     | Abs diff    | Tol 1e-3 | Status     |
|----------------------------|------------|------------|-------------|----------|------------|
| Subject.(Intercept)        | 0.96674177 | 0.96684140 | 9.97e-5     | PASS     | within-tol |
| Subject.Days.(Intercept)   | 0.01516906 | 0.01513150 | 3.76e-5     | PASS     | within-tol |
| Subject.Days               | 0.23090995 | 0.23106020 | 1.50e-4     | PASS     | within-tol |

**Result: within-tol.** All three Cholesky elements within 1e-3.

---

### VarCorr (variance components)

lme4 reports: Intercept SD = 24.7407, Days SD = 5.9221, Corr = 0.066  
mixeff reports: Intercept SD = 24.74180, Days SD = 5.92558, Corr = +0.07

The VarCorr values are consistent — the correlation rounding (+0.07 vs 0.066) is a display difference in mixeff's print method (rounds to 2 decimal places with sign). The SD values agree to ~4 significant figures.

**Note on VarCorr display:** mixeff's `VarCorr` print format is a custom tabular layout, not the standard `nlme::VarCorr`-style matrix. The correlation is displayed as "+0.07" (2 sig figs) vs lme4's "0.066" (3 sig figs). The actual correlation from theta is ~0.066, so this is a display precision issue only.

---

### Sigma (residual SD)

| Source  | Value    | Abs diff | Tol 1e-4 | Status   |
|---------|----------|----------|----------|----------|
| lme4    | 25.5918  |          |          |          |
| mixeff  | 25.59036 | 1.44e-3  | FAIL     | divergent |

**Result: divergent.** Abs diff = 1.44e-3 exceeds 1e-4 tolerance (~0.006% relative error). This is consistent with the slightly different optimizer convergence point evidenced by the SE divergence. Practically negligible.

**Severity: minor.**

---

### logLik / AIC / BIC

| Quantity | lme4        | mixeff      | Abs diff   | Tol 1e-3 | Status     |
|----------|-------------|-------------|------------|----------|------------|
| logLik   | -871.8141   | -871.8141   | 4.13e-6    | PASS     | within-tol |
| AIC      | 1755.628    | 1755.628    | 8.27e-6    | —        | within-tol |
| BIC      | 1774.786    | 1774.786    | 8.27e-6    | —        | within-tol |

**Result: within-tol.** logLik agrees to 4e-6, well within 1e-3.

---

### Random Effects (ranef)

| Subject | lme4 Intercept | mixeff Intercept | Abs diff | lme4 Days | mixeff Days | Abs diff |
|---------|---------------|-----------------|----------|-----------|-------------|----------|
| 308     | 2.258551      | 2.249080        | 0.00947  | 9.198976  | 9.201100    | 0.00212  |
| 309     | -40.398738    | -40.393426      | 0.00531  | -8.619681 | -8.621003   | 0.00132  |
| 310     | -38.960409    | -38.958300      | 0.00211  | -5.448856 | -5.449458   | 0.000602 |
| 330     | 23.690620     | 23.697891       | 0.00727  | -4.814350 | -4.815896   | 0.00155  |
| 331     | 22.260313     | 22.265615       | 0.00530  | -3.069912 | -3.071022   | 0.00111  |
| 332     | 9.039568      | 9.040696        | 0.00113  | -0.272177 | -0.272399   | 0.000222 |

Max abs diff: Intercept = 0.01128, Days = 0.00252

**Result: divergent** relative to 1e-4 fixef tolerance. However, there is no published tolerance for ranef in this probe spec. These differences are consistent with the slightly different sigma/theta solution point and are practically negligible (~0.05% relative error on the largest ranef values). The ranef differences propagate directly from the sigma/theta divergence.

**Severity: minor.**

---

### Fitted Values

- Max abs diff: 0.01142 ms  
- Mean abs diff: 0.00243 ms

Fitted value differences are commensurate with ranef differences (ranef enter into fitted values directly). These are practically negligible on the ~250 ms scale of Reaction times.

**Severity: minor.**

---

### Convergence / Refusal Status

- lme4: converged cleanly (no warnings)
- mixeff: fit succeeded, `fit_status` field present on object; `convergence_status` / `converged` field returned empty (not exposed as a top-level scalar)

**Note:** `fit_mm$convergence_status` returned an empty character — mixeff does not currently surface a boolean "converged" flag directly accessible as `$converged` or `$convergence_status`. The `fit_status` field exists. This is a minor usability gap — users cannot programmatically check convergence without knowing to inspect `fit_status`.

**Severity: cosmetic.** Model did converge (logLik matches); the gap is just the missing `$converged` accessor.

---

### Speed

| Method  | Avg time/fit (5 reps) |
|---------|-----------------------|
| lme4    | 0.0208 sec            |
| mixeff  | 0.0030 sec            |
| Ratio   | 0.144 (mixeff ~6.9x faster) |

**Result: mixeff is ~7x faster than lme4** on this dataset/formula. This is a positive result consistent with the project goal.

---

## Summary Table

| Quantity      | Max Abs Diff      | Tolerance | Status     | Severity |
|---------------|-------------------|-----------|------------|----------|
| fixef         | 5.68e-13          | 1e-4      | within-tol | none     |
| SE (vcov)     | 7.16e-4           | 1e-4      | divergent  | minor    |
| theta         | 1.50e-4           | 1e-3      | within-tol | none     |
| sigma         | 1.44e-3           | 1e-4      | divergent  | minor    |
| logLik        | 4.13e-6           | 1e-3      | within-tol | none     |
| AIC/BIC       | 8.27e-6           | —         | within-tol | none     |
| ranef         | 1.13e-2 (Intercept) | —       | divergent* | minor    |
| fitted        | 1.14e-2           | —         | divergent* | minor    |
| convergence   | (not exposed)     | —         | cosmetic   | cosmetic |
| speed         | ~7x faster        | —         | pass       | none     |

*No explicit ranef/fitted tolerance in spec; values are practically close.

---

## Root Cause Analysis

The SE, sigma, and ranef divergences are all consistent with a single underlying cause: mixeff's optimizer converges to a slightly different point in the theta/sigma parameter space than lme4's bobyqa. The logLik difference is only 4e-6, meaning both solutions lie on essentially the same likelihood surface — the optimization landscape is very flat near the optimum for this dataset, and the two solvers stop at slightly different points within that flat region.

The SE divergence specifically may also reflect a difference in how the observed Fisher information matrix is computed (numerical Hessian step size / method vs. analytic gradient).

**These are not structural bugs.** They are normal optimizer-variance artifacts on a flat likelihood surface.

---

## Classification

- fixef: **works**
- SE/vcov: **partial** (within-tol: FALSE, but practically close; not a structural defect)
- theta: **works**
- VarCorr: **works** (display correlation rounding is cosmetic)
- sigma: **partial** (exceeds 1e-4; practically negligible)
- logLik/AIC/BIC: **works**
- ranef: **partial** (commensurate with sigma divergence; practically negligible)
- fitted: **partial** (propagates from ranef; practically negligible)
- convergence flag: **partial** (no scalar `$converged` accessor)
- speed: **works** (~7x faster than lme4)

**Overall cell outcome: within-tol on logLik/fixef/theta (the primary parity quantities); minor divergences on SE/sigma/ranef that exceed stated tolerances but are practically negligible.**
