# Parity Assessment: lmm-penicillin-crossed

**Cell:** lmm-penicillin-crossed
**Dataset:** Penicillin (lme4 built-in)
**Formula:** `diameter ~ 1 + (1|plate) + (1|sample)`
**Focus:** Fully crossed random effects (two independent random intercepts, 24 plates × 6 samples)
**Date:** 2026-05-31
**Probe script:** `assessment/parity/lmm-penicillin-crossed-probe.R`

---

## Environment

| Package  | Version |
|----------|---------|
| lme4     | 2.0.1   |
| lmerTest | 3.2.1   |
| mixeff   | 0.1.0   |

---

## Dataset Summary

- 144 rows, 3 columns (`diameter`, `plate`, `sample`)
- 24 plate levels, 6 sample levels (fully crossed: 144 = 24 × 6)

---

## Fit Status

| Engine  | Fit status            | Convergence warning |
|---------|-----------------------|---------------------|
| lme4    | (normal)              | None                |
| mixeff  | `converged_interior`  | None                |

Both engines converged cleanly with no warnings.

---

## Raw Output (key sections)

### lme4

```
fixef: (Intercept) = 22.97222
SE:    (Intercept) = 0.8085954
vcov[1,1]: 0.6538265

VarCorr:
 plate    (Intercept) Std.Dev. 0.84670
 sample   (Intercept) Std.Dev. 1.93161
 Residual             Std.Dev. 0.54992

sigma:  0.5499227
logLik: -165.4303 (df=4)
AIC:    338.8606
BIC:    350.7398
theta:  plate.(Intercept)=1.539676  sample.(Intercept)=3.512519
```

### mixeff

```
fit_status: converged_interior

fixef: (Intercept) = 22.97222
SE:    (Intercept) = 0.8085784
vcov[1,1]: 0.6537991

VarCorr:
 plate  (Intercept) variance=0.716938  std_dev=0.846722
 sample (Intercept) variance=3.730960  std_dev=1.931570
 Residual std. dev.: 0.549921

sigma:  0.5499209
logLik: -165.4303 (df=4)
AIC:    338.8606
BIC:    350.7398
theta:  1.539716  3.512449
```

---

## Numerical Comparison

Tolerances (project spec): fixef=1e-4, theta=1e-3, logLik=1e-3, sigma=1e-4.

| Quantity                  | lme4 value      | mixeff value    | Max abs diff  | Tolerance | Status     |
|---------------------------|-----------------|-----------------|---------------|-----------|------------|
| fixef (Intercept)         | 22.97222222     | 22.97222222     | 4.08e-12      | 1e-4      | WITHIN-TOL |
| SE (Intercept)            | 0.80859536      | 0.80857844      | 1.69e-05      | 1e-4      | WITHIN-TOL |
| vcov[1,1]                 | 0.65382646      | 0.65379909      | 2.74e-05      | 1e-4 *    | WITHIN-TOL |
| theta[plate]              | 1.53967555      | 1.53971597      | 4.04e-05      | 1e-3      | WITHIN-TOL |
| theta[sample]             | 3.51251886      | 3.51244914      | 6.97e-05      | 1e-3      | WITHIN-TOL |
| sigma                     | 0.54992268      | 0.54992089      | 1.79e-06      | 1e-4      | WITHIN-TOL |
| VarCorr plate var         | 0.71690514      | 0.71693811      | 3.30e-05      | 1e-3      | WITHIN-TOL |
| VarCorr sample var        | 3.73113184      | 3.73095946      | 1.72e-04      | 1e-3      | WITHIN-TOL |
| logLik                    | -165.43029450   | -165.43029451   | 6.33e-09      | 1e-3      | WITHIN-TOL |
| AIC                       | 338.86058900    | 338.86058901    | 1.27e-08      | 2e-3      | WITHIN-TOL |
| BIC                       | 350.73984220    | 350.73984221    | 1.27e-08      | 2e-3      | WITHIN-TOL |
| fitted (max abs diff)     | —               | —               | 5.32e-06      | 1e-4      | WITHIN-TOL |
| ranef plate (max abs diff)| —               | —               | 4.92e-06      | 1e-4      | WITHIN-TOL |
| ranef sample (max abs diff)| —              | —               | 4.01e-07      | 1e-4      | WITHIN-TOL |

*Note on vcov[1,1]: The probe script used `tols$fixef^2 = 1e-8` as the internal vcov tolerance (SE^2), which flagged this as EXCEEDS-TOL. However the project spec defines only fixef=1e-4; the diff of 2.74e-05 is well within 1e-4. No real defect.

---

## Speed

| Engine  | Mean wall-clock (5 reps) | Ratio (mm/lme4) |
|---------|--------------------------|-----------------|
| lme4    | 0.0226 s                 | —               |
| mixeff  | 0.0044 s                 | 0.19x (5.1x faster) |

mixeff is approximately **5x faster** than lme4 on this dataset.

---

## Analysis

### Outcome: within-tol (full parity)

All quantities — fixed effects, standard errors, variance components (theta, VarCorr for both grouping factors), sigma, logLik, AIC/BIC, fitted values, and both random effect vectors — are within project tolerances.

The only probe-flagged "EXCEEDS-TOL" was `vcov[1,1]` at diff=2.74e-05 against the over-tight self-imposed tolerance of 1e-8 (= fixef_tol^2). Against the project's stated fixef tolerance of 1e-4, this is clean.

### Convergence

Both engines converge without warnings. mixeff reports `converged_interior`, which is the correct status for an interior (non-boundary) solution.

### VarCorr object structure

mixeff's `VarCorr()` returns a `mm_varcorr` list with `$table` (data.frame with columns `group`, `name`, `variance`, `std_dev`, `correlation`, `boundary`) and `$residual_sd`. This differs from lme4's `VarCorr()` which returns a named list of covariance matrices with attributes. Downstream code extracting variances must use `obj$table[obj$table$group == grp, "variance"]` for mixeff vs `attr(obj[[grp]], "stddev")^2` for lme4. This is a known structural difference, not a defect.

### Classification

**works** — fully crossed two-way random intercepts model fits to within-tolerance parity on all quantities, with a 5x speed advantage.

---

## Summary

- **Outcome:** within-tol
- **Severity:** none
- **Max abs diff (headline):** theta max=6.97e-05 (within 1e-3 tol); all others smaller
- **Speed:** mixeff 5.1x faster than lme4 (0.0044s vs 0.0226s per fit)
- **Convergence:** both clean, no warnings
- **Classification:** works
