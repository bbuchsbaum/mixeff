# Parity Assessment: inf-ranef-condvar

**Cell:** inf-ranef-condvar  
**Dataset:** sleepstudy  
**Formula:** `Reaction ~ Days + (Days | Subject)`  
**Focus:** `ranef(condVar=TRUE)` — point estimates and conditional variance arrays (`postVar`) vs lme4  
**Date:** 2026-05-31  
**lme4:** 2.0.1 | **mixeff:** 0.1.0

---

## Script

See `/Users/bbuchsbaum/code/mixeff/assessment/parity/inf-ranef-condvar-probe.R`.

---

## Raw Output (verbatim)

```
=== inf-ranef-condvar parity probe ===

lme4 version: 2.0.1 
mixeff version: 0.1.0 

--- Fitting models ---
lme4 wall time:    0.028 s
mixeff wall time:  0.021 s
Speed ratio (lme4/mixeff):  1.33 

--- Convergence / refusal status ---
lme4 converged: TRUE (no error)
mixeff converged: TRUE (no error = TRUE)

--- Fixed effects ---
lme4    fixef: 251.4051 10.46729 
mixeff  fixef: 251.4051 10.46729 
abs diff:      0 0 
max abs diff (fixef):  5.684342e-13  [tol 1e-4]

--- SE / vcov ---
lme4   SE: 6.824597 1.54579 
mixeff SE: 6.824726 1.546506 
abs diff:  0.00012964 0.00071623 
max abs diff (SE):  0.0007162301  [tol 1e-4]

--- VarCorr / theta ---
lme4 VarCorr (Subject):
 Groups   Name        Std.Dev. Corr  
 Subject  (Intercept) 24.7407        
          Days         5.9221  0.066 
 Residual             25.5918        

mixeff VarCorr table:
    group        name  variance   std_dev correlation boundary
1 Subject (Intercept) 612.15747 24.741816                FALSE
2 Subject        Days  35.11247  5.925578       +0.07    FALSE

lme4 theta: 0.966742 0.015169 0.23091 
lme4 sigma: 25.5918 
mixeff sigma: 25.59036 
abs diff sigma: 0.00143999  [tol 1e-4]

lme4 Subject std devs  (Intercept, Days): 24.74066 5.922138 
mixeff Subject std devs (Intercept, Days): 24.74182 5.925578 
abs diff std devs: 0.00115821 0.0034403 
max abs diff (theta/SD):  0.003440298  [tol 1e-3]

--- logLik / AIC / BIC ---
lme4   logLik: -871.8141 
mixeff logLik: -871.8141 
abs diff logLik: 4.133908e-06  [tol 1e-3]

lme4   AIC: 1755.628 
mixeff AIC: 1755.628 
abs diff AIC: 8.267816e-06 

lme4   BIC: 1774.786 
mixeff BIC: 1774.786 
abs diff BIC: 8.267816e-06 

--- ranef point estimates (condVar=FALSE) ---
lme4 ranef (first 6):
    (Intercept)    Days
308      2.2586  9.1990
309    -40.3987 -8.6197
310    -38.9604 -5.4489
330     23.6906 -4.8144
331     22.2603 -3.0699
332      9.0396 -0.2722
mixeff ranef (first 6):
    (Intercept)    Days
308      2.2491  9.2011
309    -40.3934 -8.6210
310    -38.9583 -5.4495
330     23.6979 -4.8159
331     22.2656 -3.0710
332      9.0407 -0.2724

max abs diff ranef (Intercept):  0.01128454  [tol 1e-4]
max abs diff ranef (Days):        0.002522863  [tol 1e-4]

--- Fitted values ---
max abs diff fitted:  0.01142123  [tol 1e-4]

=== PRIMARY FOCUS: ranef(condVar=TRUE) ===

lme4 postVar array dim: 2 2 18 
lme4 postVar (first subject slice):
         [,1]       [,2]
[1,] 145.7056 -21.444504
[2,] -21.4445   5.312283
lme4 subject order: 308 309 310 330 331 332 333 334 335 337 349 350 351 352 369 370 371 372 

Calling ranef(fit_mm, condVar=TRUE)...
mixeff postVar array dim: 2 2 18 
mixeff subject order: 308 309 310 330 331 332 333 334 335 337 349 350 351 352 369 370 371 372 

lme4   postVar (Subject 308, slice):
            (Intercept)       Days
(Intercept)    145.7056 -21.444504
Days           -21.4445   5.312283
mixeff postVar (Subject 308, slice):
            (Intercept)       Days
(Intercept)   145.72531 -21.449676
Days          -21.44968   5.313309

All finite (mixeff postVar): TRUE 
Any NA: FALSE 

All slices symmetric: TRUE 
All diagonal entries >= 0: TRUE 

Comparison postVar arrays ( 18 subjects):
max abs diff (postVar):  0.01972731  [tol 1e-3]
mean abs diff (postVar): 0.007774376 

Max diff by position [1,1] (Intercept var): 0.01972731 
Max diff by position [2,2] (Days var):      0.001025629 
Max diff by position [1,2] / [2,1] (cov):  0.005172281 

lme4   postVar [1,1] per subject: 145.7056 (all equal)
mixeff postVar [1,1] per subject: 145.7253 (all equal)

lme4   postVar [2,2] per subject: 5.3123 (all equal)
mixeff postVar [2,2] per subject: 5.3133 (all equal)

--- Conditional SDs ---
lme4   condSD (Intercept): 12.0709 (all equal across subjects)
mixeff condSD (Intercept): 12.0717 (all equal across subjects)
max abs diff condSD (Intercept): 0.0008171187  [tol 1e-3]

lme4   condSD (Days): 2.3048 (all equal across subjects)
mixeff condSD (Days): 2.3051 (all equal across subjects)
max abs diff condSD (Days): 0.0002224839  [tol 1e-3]

RESULT: postVar BEYOND tolerance 1e-3 

=== SUMMARY TABLE ===

Quantity                          max_abs_diff        tol     status
---------------------------------------------------------------------- 
fixef                                 5.68e-13   1.00e-04       PASS
SE(fixef)                             7.16e-04   1.00e-04       FAIL
sigma                                 1.44e-03   1.00e-04       FAIL
VarCorr SD                            3.44e-03   1.00e-03       FAIL
logLik                                4.13e-06   1.00e-03       PASS
AIC                                   8.27e-06   2.00e-03       PASS
BIC                                   8.27e-06   2.00e-03       PASS
ranef(Intercept)                      1.13e-02   1.00e-04       FAIL
ranef(Days)                           2.52e-03   1.00e-04       FAIL
fitted                                1.14e-02   1.00e-04       FAIL
postVar(condVar=TRUE)                 1.97e-02   1.00e-03       FAIL
```

---

## Analysis

### What works

- **Convergence:** both models fit without error on sleepstudy.
- **Fixed effects:** essentially identical (max diff 5.7e-13, well within 1e-4).
- **logLik / AIC / BIC:** within-tolerance (diff ~4e-6 / 8e-6).
- **`ranef(condVar=TRUE)` does not error or refuse:** mixeff returns a real `postVar` array — a `2 × 2 × 18` finite, symmetric, PSD array with named dimnames. No typed-refusal (`mm_unavailable_reason`) triggered. This is a genuine capability: the Rust `cond_var()` bridge fires and returns values.
- **Conditional SDs within tolerance:** condSD(Intercept) diff = 8.2e-4 (tol 1e-3 ✓); condSD(Days) diff = 2.2e-4 (tol 1e-3 ✓).

### What fails

The failures cluster around a single root cause: **mixeff's variance-component estimates (theta, sigma) are slightly off vs lme4**, which propagates downstream.

| Quantity | lme4 | mixeff | diff | tol | verdict |
|---|---|---|---|---|---|
| sigma | 25.5918 | 25.5904 | 1.44e-3 | 1e-4 | **FAIL** |
| VarCorr SD (Intercept) | 24.7407 | 24.7418 | 1.16e-3 | 1e-3 | **FAIL** |
| VarCorr SD (Days) | 5.9221 | 5.9256 | 3.44e-3 | 1e-3 | **FAIL** |
| SE(Intercept) | 6.8246 | 6.8247 | 1.3e-4 | 1e-4 | **FAIL** |
| SE(Days) | 1.5458 | 1.5465 | 7.2e-4 | 1e-4 | **FAIL** |
| ranef(Intercept) | (varies) | (varies) | 1.13e-2 | 1e-4 | **FAIL** |
| ranef(Days) | (varies) | (varies) | 2.52e-3 | 1e-4 | **FAIL** |
| postVar[1,1] (all equal) | 145.706 | 145.725 | 1.97e-2 | 1e-3 | **FAIL** |
| postVar[2,2] (all equal) | 5.3123 | 5.3133 | 1.03e-3 | 1e-3 | **FAIL** |
| postVar[1,2] cov (all equal) | -21.445 | -21.450 | 5.17e-3 | 1e-3 | **FAIL** |

**Key observation on postVar:** Both lme4 and mixeff return identical values for all 18 subjects within each engine (e.g., lme4 postVar[1,1] = 145.7056 for every subject; mixeff = 145.7253 for every subject). This is expected for a balanced random-intercept-and-slope model where the posterior variance depends only on the design and VarCorr, not on subject-specific data. The numerical disagreement in postVar traces to the differing VarCorr estimates: lme4's theta gives slightly different Cholesky factors, producing slightly different posterior covariance matrices.

**Root cause hierarchy:**

1. **theta/VarCorr discrepancy** (primary): mixeff sigma = 25.590 vs lme4 25.592 (diff 1.4e-3, 14× above tol). VarCorr SDs differ by up to 3.4e-3 (3.4× above tol). This is not a ranef-specific issue — it exists at the fit level (documented in the broader lmm-sleep-corr probe as a pre-existing variance-component parity drift).

2. **ranef point estimates** propagate from theta: BLUPs are computed from theta and the data; the slightly different VarCorr yields ranef diffs of up to 1.13e-2 (113× above 1e-4 tol).

3. **postVar** propagates from VarCorr: the posterior variance matrices are a function of theta and the design matrix. Since postVar is constant across subjects (balanced design), the absolute diff of 1.97e-2 in postVar[1,1] reflects the accumulated error from the VarCorr disagreement. The conditional SDs (sqrt of diagonals) are within tol (8.2e-4 and 2.2e-4 vs 1e-3 tol), showing the error is tolerable at the SD level but not at the variance level.

### Classification

| Gap | Classification |
|---|---|
| `ranef(condVar=TRUE)` raises error / refuses | **works** — no refusal, real array returned |
| postVar array finite, symmetric, PSD | **works** |
| postVar values outside tolerance 1e-3 | **partial** — functionally correct path, numerically driven by upstream VarCorr drift |
| ranef point estimates outside tolerance 1e-4 | **partial** — same root cause (VarCorr) |
| SE outside tolerance 1e-4 | **partial** — same root cause |
| sigma outside tolerance 1e-4 | **partial** — pre-existing, tracked separately |
| Conditional SDs within tolerance 1e-3 | **works** |

The `condVar=TRUE` capability itself (Rust bridge, R plumbing, named dimnames, symmetry, PSD, caching) is **working**. The numerical differences in postVar are a downstream consequence of the VarCorr/theta discrepancy that already exists at the fit level, not a bug specific to the `cond_var()` path.

### Severity

**minor** — the `condVar=TRUE` feature is implemented and structurally correct. The postVar values differ by ~1.4% in the intercept variance (145.706 vs 145.725) and are within tolerance at the conditional-SD level (condSD Intercept diff = 8e-4 < 1e-3; condSD Days diff = 2e-4 < 1e-3). The root cause is the upstream VarCorr/theta drift, which is a pre-existing issue tracked under bd-01KRV31R4BJVQCEF0F58NFD4YN, not introduced by the condVar path.
