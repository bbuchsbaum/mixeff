# Parity Assessment: inf-getME

**Cell:** inf-getME
**Date:** 2026-05-31
**Dataset:** sleepstudy (lme4 built-in, N=180, 18 subjects × 10 days)
**Formula:** `Reaction ~ Days + (Days|Subject)` (correlated random slope)
**REML:** TRUE
**Focus:** `getME()` names parity: X, Z, theta, beta, Lambda dimensions/values

---

## Environment

| Package | Version |
|---------|---------|
| lme4    | 2.0.1   |
| mixeff  | 0.1.0   |

---

## Script

Written to `/tmp/inf_getME_probe.R`. All lme4 extractor calls use explicit `lme4::getME()` to avoid dispatch to mixeff's `getME` S3 method (which masks lme4's via `zzz.R` registration). All mixeff calls use explicit `mixeff::getME()`.

---

## Raw Output

```
=== ENVIRONMENT ===
lme4 version: 2.0.1 
mixeff version: 0.1.0 

=== FIT lme4 ===
lme4 wall time: 0.025 sec

=== FIT mixeff ===
Random effects explanation:
  formula: Reaction ~ 1 + Days + (1 + Days | Subject)

Random effects:
  r0:
    wrote:      (Days | Subject)
    canonical:  (1 + Days | Subject)
    named form: re(group = Subject, intercept = TRUE, slopes = Days, cov = "full")
    scope:      `Subject` units differ in baseline and `Days` slope; the model estimates whether
                these are associated.
    covariance: full; theta parameters: 3
    support:    sufficient; group levels: 18; min rows/group: 10; median rows/group: 10
    variation:  Days=present; intercept=not_assessed
mixeff wall time: 0.013 sec
mixeff fit class: mm_lmm, mm_fit, mm_compiled 

=== getME: X (fixed-effect design matrix) ===
lme4 X dim: 180 x 2 
mixeff X dim: 180 x 2 
dim match: TRUE 
Max abs diff X: 0 
Within tol (1e-10): TRUE 
lme4 X colnames: (Intercept), Days 
mixeff X colnames: (Intercept), Days 
colnames match: TRUE 

=== getME: Z (random-effect design matrix) ===
lme4 Z dim: 180 x 36 
mixeff Z dim: 180 x 36 
dim match: TRUE 
Max abs diff Z (direct column order): 0 
Within tol (1e-10): TRUE 
lme4 Z nnz: 342 
mixeff Z nnz: 342 
nnz match: TRUE 
Max abs diff sorted Z col-norms: 0 
Within tol (1e-10): TRUE 

=== getME: theta ===
lme4 theta (length 3 ): 0.9667418 0.01516906 0.23091 
lme4 theta names: Subject.(Intercept), Subject.Days.(Intercept), Subject.Days 
mixeff theta (length 3 ): 0.9668414 0.0151315 0.2310602 
mixeff theta names:  
length match: TRUE 
Abs diff theta: 9.965909e-05 3.756139e-05 0.0001502393 
Max abs diff theta: 0.0001502393 
Within tol (1e-3): TRUE 

=== getME: beta / fixef ===
lme4 getME(beta): 251.4051 10.46729 
mixeff getME(beta): 251.4051 10.46729 
names match: FALSE 
lme4 beta names: 
mixeff beta names: (Intercept) Days 
Abs diff beta: 5.684342e-13 1.865175e-13 
Max abs diff beta: 5.684342e-13 
Within tol (1e-4): TRUE 
getME('beta') == getME('fixef') lme4: FALSE 
getME('beta') == getME('fixef') mixeff: TRUE 

=== getME: Lambda (relative covariance factor) ===
lme4 Lambda dim: 36 x 36 
mixeff Lambda dim: 36 x 36 
dim match: TRUE 
lme4 Lambda nnz: 54 
mixeff Lambda nnz: 54 
nnz match: TRUE 
Max abs diff Lambda (direct): 0.0001502393 
Within tol (1e-3): TRUE 
Max abs diff Lambda'Lambda eigenvalues (sorted): 0.0001915083 
Within tol (1e-3): TRUE 

=== getME: Lambdat (transposed Lambda) ===
lme4 Lambdat dim: 36 x 36 
mixeff Lambdat ERROR: argument is not a matrix 

=== getME: Zt (transposed Z) ===
lme4 Zt dim: 36 x 180 
mixeff Zt ERROR: argument is not a matrix 

=== getME: flist (grouping factors) ===
lme4 flist names: Subject 
mixeff flist names: Subject 
names match: TRUE 
lme4 flist$Subject levels: 308, 309, 310, 330, 331, 332, 333, 334, 335, 337, 349, 350,
                            351, 352, 369, 370, 371, 372 
mixeff flist$Subject levels: 308, 309, 310, 330, 331, 332, 333, 334, 335, 337, 349, 350,
                              351, 352, 369, 370, 371, 372 
levels match: TRUE 
values match: TRUE 

=== getME: cnms (random coefficient names) ===
lme4 cnms: Subject -> (Intercept), Days 
mixeff cnms: Subject -> (Intercept), Days 
group names match: TRUE 
Subject cnms match: TRUE 

=== getME: y (response vector) ===
lme4 y length: 180 
mixeff y length: 180 
Max abs diff y: 0 
Within tol (1e-10): TRUE 

=== getME: mu (fitted values) ===
lme4 mu length: 180 
mixeff mu length: 180 
Max abs diff mu: 0.01142123 
Within tol (1e-4): FALSE 

=== getME: unavailable names (error handling) ===
lme4 unknown name response: ERROR: 'arg' should be one of "X", "Z", "Zt", "Ztlist",
  "mmList", "y", "mu", "u", "b", "Gp", "Tp", "L", "Lambda", "Lambdat", "Lind", "Tlist",
  "A", "RX", "RZX", "sigma", "flist", "fixef", "beta", "theta", "ST", "par", "REML",
  "is_REML", "n_rtrms", "n_rfacs", "N", "n", "p", "q", "p_i", "l_i", "q_i", "k",
  "m_i", "m", "cnms", "devcomp", "offset", "lower", "devfun", "devarg",
  "glmer.nb.theta" 
mixeff unknown name response: ERROR: `getME()` component `NONEXISTENT` is not available. 

=== lme4 getME names supported by lme4 but tested on mixeff ===
  'sigma       ': lme4=TRUE     mixeff=FALSE    gap=YES
  'delta       ': lme4=FALSE    mixeff=FALSE    gap=both_missing
  'offset      ': lme4=TRUE     mixeff=FALSE    gap=YES
  'Lind        ': lme4=TRUE     mixeff=FALSE    gap=YES
  'is_REML     ': lme4=TRUE     mixeff=FALSE    gap=YES
  'ST          ': lme4=TRUE     mixeff=FALSE    gap=YES
  'n_rtrms     ': lme4=TRUE     mixeff=FALSE    gap=YES

=== CORE SUMMARY TABLE ===
Component            lme4 dim/len   mm dim/len     dim_match    max_abs_diff   in_tol    
------------------------------------------------------------------------------------------
X                    180x2          180x2          TRUE         0e+00          TRUE      
Z (direct)           180x36         180x36         TRUE         0e+00          TRUE      
theta                3              3              TRUE         1.5e-04        TRUE      
beta/fixef           2              2              TRUE         5.68e-13       TRUE      
Lambda               36x36          36x36          TRUE         1.5e-04        TRUE      
Lambdat              36x36          ERROR          ERROR        ERROR          ERROR     
Zt                   36x180         ERROR          ERROR        ERROR          ERROR     
y                    180            180            TRUE         0e+00          TRUE      
mu (fitted)          180            180            TRUE         1.14e-02       FALSE     
------------------------------------------------------------------------------------------
```

---

## Quantity-by-Quantity Analysis

### X — Fixed-Effect Design Matrix

- lme4: 180×2, columns `(Intercept)`, `Days`
- mixeff: 180×2, columns `(Intercept)`, `Days`
- Max abs diff: **0** (exact match)
- Column names: **identical**

**Result: works.** Exact match to machine precision.

---

### Z — Random-Effect Design Matrix

- lme4: 180×36 sparse (18 subjects × 2 RE terms)
- mixeff: 180×36 sparse
- Max abs diff (direct column order): **0** (exact match)
- nnz: both 342 (identical sparsity pattern)
- Sorted column norms: max diff **0**

**Result: works.** Exact match including column ordering. Both packages produce the same block-level assignment matrix.

---

### theta — Cholesky Factor Elements

- lme4 theta: `[0.9667418, 0.01516906, 0.23091]`, names: `Subject.(Intercept)`, `Subject.Days.(Intercept)`, `Subject.Days`
- mixeff theta: `[0.9668414, 0.0151315, 0.2310602]`, **no names**
- Length: both 3
- Max abs diff: **1.50e-4** — within tolerance (1e-3)

**Gap — missing names:** lme4's theta vector carries element names (`Subject.(Intercept)`, etc.); mixeff's theta is an unnamed numeric vector. This is a usability gap — programmatic code that relies on `names(getME(fit, "theta"))` to identify which Cholesky element corresponds to which RE term will receive `NULL` from mixeff.

**Numerical result: within-tol.** Values agree within 1e-3.

**Severity: minor** (missing theta names; values OK).

---

### beta / fixef

- lme4 `getME(fit, "beta")`: `[251.4051, 10.46729]`, **unnamed** (lme4 returns beta without names via getME)
- mixeff `getME(fit, "beta")`: `[251.4051, 10.46729]`, names: `(Intercept)`, `Days`
- Max abs diff: **5.68e-13** — machine precision
- `getME("beta") == getME("fixef")`: FALSE for lme4 (lme4 returns different objects), TRUE for mixeff (both return identical named vector)

**Note:** lme4's `getME(fit, "beta")` returns an unnamed numeric vector while `getME(fit, "fixef")` returns a named vector. mixeff's `getME(fit, "beta")` and `getME(fit, "fixef")` both return the same named vector — this is actually *better* than lme4's behavior. The `names match: FALSE` flag reflects that lme4 strips names from `beta` whereas mixeff keeps them.

**Result: works** (values identical; mixeff is more informative than lme4 here).

---

### Lambda — Relative Covariance Factor

- lme4: 36×36 sparse, nnz=54
- mixeff: 36×36 sparse, nnz=54
- Max abs diff (direct): **1.50e-4** — within tolerance (1e-3)
- Lambda'Lambda eigenvalues max diff: **1.92e-4** — within tolerance (1e-3)

The 1.50e-4 diff is the same magnitude as the theta diff, confirming Lambda is correctly rebuilt from stored theta.

**Result: within-tol.** Lambda dimensions, sparsity pattern, and values all correct within 1e-3.

---

### Lambdat — Transposed Lambda

**Error in mixeff:**
```
Error in t.default(.mm_lazy(object, "Lambda", mm_lambda_matrix)):
  argument is not a matrix
```

**Root cause:** In `revive.R` line 130, `getME.mm_lmm` implements `Lambdat` as:
```r
Lambdat = t(.mm_lazy(object, "Lambda", mm_lambda_matrix)),
```
`mm_lambda_matrix()` returns a `Matrix::sparseMatrix` object. The bare `t()` call dispatches to `base::t.default` instead of `Matrix::t` because there is no explicit `Matrix::` qualification and `Matrix` is not re-exported into the function's call environment at that point. `t.default` rejects non-matrix objects with `"argument is not a matrix"`.

**Result: in-scope bug.** `getME("Lambdat")` crashes with a method dispatch error. lme4 succeeds. Fix: change `t(...)` to `Matrix::t(...)` in `getME.mm_lmm`.

**Severity: major** — crashes with an error on a documented supported name.

---

### Zt — Transposed Z

**Error in mixeff:**
```
Error in t.default(stats::model.matrix(object, type = "random")):
  argument is not a matrix
```

**Root cause:** Same bug as Lambdat. In `revive.R` line 128:
```r
Zt = t(stats::model.matrix(object, type = "random")),
```
`model.matrix.mm_lmm` returns a sparse `Matrix` object. Again `t()` dispatches to `t.default`. Fix: `Matrix::t(stats::model.matrix(object, type = "random"))`.

**Result: in-scope bug.** `getME("Zt")` crashes. lme4 succeeds. Same root cause as Lambdat — a single character fix in `revive.R` handles both.

**Severity: major** — crashes on a documented supported name.

---

### flist — Grouping Factors

- lme4: `$Subject` factor, 18 levels
- mixeff: `$Subject` factor, 18 levels
- Group names: **identical** (`Subject`)
- Levels: **identical** (308, 309, ..., 372)
- Values: **identical** (factor assignments match observation-by-observation)

**Result: works.** Perfect parity.

---

### cnms — Random Coefficient Names

- lme4: `Subject -> (Intercept), Days`
- mixeff: `Subject -> (Intercept), Days`
- Group names: **identical**
- Subject cnms: **identical** (both have `(Intercept)` and `Days`)

**Result: works.** Perfect parity.

---

### y — Response Vector

- Both: length 180
- Max abs diff: **0** (exact match)

**Result: works.**

---

### mu — Fitted Values (via getME)

- Both: length 180
- Max abs diff: **0.01142** — exceeds 1e-4 tolerance

This is consistent with the sigma/ranef divergence documented in `lmm-sleep-corr.md` — mixeff's optimizer converges to a slightly different theta/sigma point. The fitted-value divergence propagates from the ranef differences, which in turn propagate from the optimizer stopping at a slightly different point on a flat likelihood surface.

**Result: divergent** relative to 1e-4. Practically negligible on the ~250 ms scale. Same underlying cause as the main lmm-sleep-corr findings.

**Severity: minor** (pre-existing optimizer convergence issue, not a getME-specific defect).

---

### Unsupported lme4 getME Names

| Name      | lme4 | mixeff | Gap classification |
|-----------|------|--------|--------------------|
| sigma     | OK   | ERROR  | in-scope-missing   |
| delta     | FAIL | FAIL   | both_missing (lme4 also doesn't have it) |
| offset    | OK   | ERROR  | out-of-scope-by-design (offset not yet in mixeff) |
| Lind      | OK   | ERROR  | in-scope-missing (internal parameterization detail) |
| is_REML   | OK   | ERROR  | in-scope-missing (useful flag) |
| ST        | OK   | ERROR  | out-of-scope-by-design (internal lme4 representation) |
| n_rtrms   | OK   | ERROR  | in-scope-missing (number of random terms) |

These are gap names not in the core probe spec (X/Z/theta/beta/Lambda). They are noted for completeness. `sigma`, `is_REML`, and `n_rtrms` are the most commonly used by downstream packages (e.g. merTools, emmeans) and represent in-scope-missing gaps.

---

### Error Handling for Unknown Names

- lme4: emits a `match.arg` error listing all valid names
- mixeff: emits `getME() component 'NONEXISTENT' is not available.` — clear, structured, honest

**Result: works.** mixeff's error message is arguably cleaner than lme4's for unknown names.

---

## Root Cause Summary: Lambdat and Zt Bug

Both failures share one root cause: bare `t()` on a `Matrix` sparse object.

**Location:** `/Users/bbuchsbaum/code/mixeff/R/revive.R`, function `getME.mm_lmm`, lines 128–130:

```r
Zt     = t(stats::model.matrix(object, type = "random")),   # line 128
...
Lambdat = t(.mm_lazy(object, "Lambda", mm_lambda_matrix)),  # line 130
```

**Fix:** Replace both with `Matrix::t(...)`. This is a two-character fix per line. The `Matrix` package is already imported and used throughout the package.

---

## Summary Table

| getME Component | lme4 dim/len | mixeff dim/len | dim_match | max_abs_diff | tol    | Status      | Severity |
|-----------------|-------------|----------------|-----------|--------------|--------|-------------|----------|
| X               | 180×2       | 180×2          | TRUE      | 0            | 1e-10  | within-tol  | none     |
| Z               | 180×36      | 180×36         | TRUE      | 0            | 1e-10  | within-tol  | none     |
| theta           | 3           | 3              | TRUE      | 1.50e-4      | 1e-3   | within-tol  | minor*   |
| beta/fixef      | 2           | 2              | TRUE      | 5.68e-13     | 1e-4   | within-tol  | none     |
| Lambda          | 36×36       | 36×36          | TRUE      | 1.50e-4      | 1e-3   | within-tol  | none     |
| Lambdat         | 36×36       | ERROR          | ERROR     | ERROR        | —      | divergent   | major    |
| Zt              | 36×180      | ERROR          | ERROR     | ERROR        | —      | divergent   | major    |
| flist           | {Subject}   | {Subject}      | TRUE      | exact        | —      | within-tol  | none     |
| cnms            | {Subject}   | {Subject}      | TRUE      | exact        | —      | within-tol  | none     |
| y               | 180         | 180            | TRUE      | 0            | 1e-10  | within-tol  | none     |
| mu (fitted)     | 180         | 180            | TRUE      | 1.14e-2      | 1e-4   | divergent   | minor†   |

*theta values within-tol but **unnamed** — lme4 provides named theta, mixeff does not.
†mu divergence is pre-existing optimizer convergence issue, not a getME-specific defect.

---

## Classification

- X: **works**
- Z: **works**
- theta values: **within-tol** / theta names: **partial** (missing names)
- beta/fixef: **works** (mixeff actually more informative than lme4 here)
- Lambda: **within-tol**
- Lambdat: **in-scope bug** — crashes with `t.default` dispatch error; fix is `Matrix::t()`
- Zt: **in-scope bug** — crashes with `t.default` dispatch error; same fix
- flist: **works**
- cnms: **works**
- y: **works**
- mu: **partial** (pre-existing optimizer convergence divergence)
- sigma/is_REML/n_rtrms/Lind: **in-scope-missing** (not in core probe spec)

**Overall cell outcome: mixed.** Core primary quantities (X, Z, beta, Lambda) are within tolerance. Two documented getME names (`Zt`, `Lambdat`) crash with a method dispatch bug that has a trivial two-line fix (`t()` → `Matrix::t()`). The bug is in `/Users/bbuchsbaum/code/mixeff/R/revive.R` lines 128 and 130.
