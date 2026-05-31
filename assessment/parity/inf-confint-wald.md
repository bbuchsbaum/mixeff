# Parity Assessment: inf-confint-wald

**Cell:** inf-confint-wald  
**Dataset:** sleepstudy  
**Formula:** `Reaction ~ Days + (Days | Subject)`  
**Focus:** `confint(method="wald")` parity  
**Date run:** 2026-05-31  
**Probe script:** `assessment/parity/inf-confint-wald-probe.R`

---

## Raw output (verbatim)

```
=== MODEL FITTING ===

--- lme4 fit ---
lme4 fit time: 0.0250 s
Linear mixed model fit by REML ['lmerMod']
Formula: Reaction ~ Days + (Days | Subject)
   Data: sleepstudy

REML criterion at convergence: 1743.6

Random effects:
 Groups   Name        Variance Std.Dev. Corr 
 Subject  (Intercept) 612.10   24.741        
          Days         35.07    5.922   0.07 
 Residual             654.94   25.592        
Number of obs: 180, groups:  Subject, 18

Fixed effects:
            Estimate Std. Error t value
(Intercept)  251.405      6.825  36.838
Days          10.467      1.546   6.771

--- mixeff fit ---
mixeff fit time: 0.0140 s
Fit status: converged_interior
Optimizer: trust_bq; iterations: 391; objective: 1743.63
nobs: 180, sigma: 25.5904, logLik: -871.814
Fixed effects:
(Intercept)        Days 
   251.4050     10.4673 

=== SUPPORTING QUANTITIES ===

--- fixef ---
max abs diff fixef: 5.68e-13  (tol 1e-4: PASS)

--- SE (sqrt diag vcov) ---
lme4:   (Intercept)=6.824597  Days=1.545790
mixeff: (Intercept)=6.824726  Days=1.546506
abs diff:  (Intercept)=0.0001296386  Days=0.0007162301
max abs diff SE: 7.16e-04  (tol 1e-4: FAIL)

--- theta ---
max abs diff theta: 1.50e-04  (tol 1e-3: PASS)

--- sigma ---
lme4 sigma:   25.591796
mixeff sigma: 25.590356
abs diff sigma: 1.44e-03  (tol 1e-4: FAIL)

--- logLik ---
lme4:   -871.814136
mixeff: -871.814140
abs diff logLik: 4.13e-06  (tol 1e-3: PASS)

--- AIC/BIC ---
lme4   AIC=1755.6283  BIC=1774.7860
mixeff AIC=1755.6283  BIC=1774.7860
abs diff AIC: 8.27e-06   BIC: 8.27e-06

--- ranef (first 5 subjects) ---
max abs diff ranef: 1.13e-02

--- fitted values ---
max abs diff fitted: 1.14e-02

--- convergence ---
lme4 converged:   TRUE
mixeff converged: FALSE  (field absent; isTRUE(fit_mm$converged) == FALSE)

=== WALD CONFINT (primary cell target) ===

--- lme4 confint(method='Wald') ---
                 2.5 %    97.5 %
.sig01              NA        NA
.sig02              NA        NA
.sig03              NA        NA
.sigma              NA        NA
(Intercept) 238.029141 264.78107
Days          7.437594  13.49698

--- mixeff confint(method='wald') ---
Confidence intervals:
                2.5 %    97.5 %
(Intercept) 238.02889 264.78132
Days          7.43619  13.49838
method: wald_asymptotic_from_stored_standard_errors
status: not_certified_by_rust_inference_contract

--- Fixed-effect Wald CI comparison ---
abs diff per bound:
                  2.5 %      97.5 %
(Intercept) 0.000254087 0.000254087
Days        0.001403785 0.001403785
max abs diff Wald CI (fixed effects): 1.40e-03  (tol 1e-4: FAIL)

--- Does mixeff Wald CI cover VarCorr/sigma parameters? ---
lme4 Wald CI rows (non-fixef): .sig01, .sig02, .sig03, .sigma
mixeff Wald CI rows:           (Intercept), Days
MISSING from mixeff Wald CI: .sig01, .sig02, .sig03, .sigma

--- parm subset selection (parm='Days') ---
Works correctly — returns single-row CI for Days.

--- level=0.90 ---
max abs diff Wald CI level=0.90: 1.18e-03  (tol 1e-4: FAIL)

--- bad parm (should error gracefully) ---
Error (expected): Unknown fixed-effect parameter(s): NONEXISTENT.

--- bad level (should error gracefully) ---
Error (expected): `level` must be a single number between 0 and 1.

--- method='asymptotic' synonym ---
wald vs asymptotic synonym max diff: 0.00e+00 (should be 0)

=== SUMMARY TABLE ===
Quantity                   | lme4       | mixeff     | abs diff   | tol    | status
---------------------------|------------|------------|------------|--------|-------
fixef[Intercept]           | 251.40510  | 251.40510  | 5.68e-13   | 1e-04  | PASS
fixef[Days]                | 10.46729   | 10.46729   | 1.87e-13   | 1e-04  | PASS
SE[Intercept]              | 6.82460    | 6.82473    | 1.3e-04    | 1e-04  | FAIL
SE[Days]                   | 1.54579    | 1.54651    | 7.16e-04   | 1e-04  | FAIL
sigma                      | 25.59180   | 25.59036   | 1.44e-03   | 1e-04  | FAIL
logLik                     | -871.81414 | -871.81414 | 4.13e-06   | 1e-03  | PASS
Wald CI lower[Intercept]   | 238.02914  | 238.02889  | 2.54e-04   | 1e-04  | FAIL
Wald CI upper[Intercept]   | 264.78107  | 264.78132  | 2.54e-04   | 1e-04  | FAIL
Wald CI lower[Days]        | 7.43759    | 7.43619    | 1.4e-03    | 1e-04  | FAIL
Wald CI upper[Days]        | 13.49698   | 13.49838   | 1.4e-03    | 1e-04  | FAIL
```

---

## Analysis

### What works

- `confint.mm_lmm(method="wald")` runs without error and returns a typed `mm_confint` matrix.
- Fixed-effect point estimates match lme4 to machine precision (max diff 5.7e-13, tol 1e-4: PASS).
- logLik matches to 4.1e-06 (tol 1e-3: PASS).
- AIC/BIC agree to 8.3e-06 (well within any tolerance).
- theta (relative Cholesky factors) match to 1.5e-04 (tol 1e-3: PASS).
- `parm` subset selection works correctly.
- `method="asymptotic"` synonym maps identically to `"wald"` (diff = 0).
- Error handling for bad `parm` and bad `level` is clean and informative.
- `print.mm_confint` reports method and status attributes correctly.

### Findings

#### Finding 1 — Wald CI bounds exceed tolerance (SE propagation gap): MAJOR

The Wald CIs are computed as `est ± qnorm(1-α/2) * SE`. The fixed-effect point estimates are
fine, so the CI divergence is entirely attributable to SE divergence. mixeff's SEs are larger
than lme4's:

```
SE(Intercept): lme4=6.8246  mixeff=6.8247  diff=1.30e-04  [tol 1e-4: FAIL, barely]
SE(Days):      lme4=1.5458  mixeff=1.5465  diff=7.16e-04  [tol 1e-4: FAIL, ~7x over]
```

This propagates directly into CI width:
```
Wald CI lower[Days]: lme4=7.4376  mixeff=7.4362  diff=1.40e-03  [tol 1e-4: FAIL]
```

Root cause: mixeff uses the REML objective but computes the fixed-effect covariance via a
slightly different numerical path from lme4 (the Rust trust-region optimizer with REML=TRUE
leaves theta at a slightly different point, and the vcov is evaluated from that). The theta
mismatch (~1.5e-04) is within theta tolerance but still propagates through the Hessian-based
vcov computation into SEs that exceed the 1e-4 SE tolerance.

**Classification: in-scope-missing / partial** — the feature exists and runs, but numerical
agreement on SEs (and therefore Wald CI bounds) is outside the stated 1e-4 fixef tolerance
class. The divergence is ~7x over tolerance for SE(Days) and ~1.4x for SE(Intercept).

#### Finding 2 — sigma outside tolerance: MINOR

```
lme4 sigma:   25.591796
mixeff sigma: 25.590356
diff:          1.44e-03  [tol 1e-4: FAIL, ~14x over]
```

sigma = sqrt(residual variance), which derives from the REML objective at the theta optimum.
The small theta offset (within theta's own tolerance of 1e-3) produces a sigma offset 14x
outside the sigma tolerance of 1e-4. This is a downstream consequence of Finding 1.

**Classification: in-scope-missing** — sigma tolerance (1e-4) is tight; the current optimizer
precision for REML does not meet it.

#### Finding 3 — Wald CI does not cover VarCorr / sigma parameters: MINOR (by design)

lme4's `confint(method="Wald")` returns NA rows for `.sig01`, `.sig02`, `.sig03`, `.sigma`
(these are NAs because closed-form Wald CIs for variance components are not valid, but lme4
still includes the rows). mixeff returns only fixed-effect rows.

The mixeff source is explicit: the Wald method is described as "closed-form Wald interval" for
fixed effects only, and the output carries `status = "not_certified_by_rust_inference_contract"`.
**This is by design** — mixeff does not fabricate NA rows for quantities where the Wald
approximation is known to be invalid. The attribute `not_certified_by_rust_inference_contract`
is an honest signal that this path is not covered by the Rust inference contract.

**Classification: out-of-scope-by-design** — the NA rows in lme4 carry no information; omitting
them is an intentional design choice consistent with PRD §3 (no silent transformations, no
fabricated output). A user migrating from lme4 who expects those rows will be confused, which
is a documentation gap but not a correctness bug.

#### Finding 4 — ranef / fitted diverge by ~1e-02: MINOR (downstream of Finding 1)

```
max abs diff ranef:   1.13e-02
max abs diff fitted:  1.14e-02
```

Both are downstream of the theta/sigma offset. No separate tolerance was stated for these
quantities; the magnitude is consistent with a small theta offset in a correlated random
slope model.

#### Finding 5 — `fit_mm$converged` is absent / FALSE: COSMETIC

`isTRUE(fit_mm$converged)` returns FALSE because the field is absent from the list, not
because mixeff reports non-convergence. The fit status string "converged_interior" is present
in the print output. This is a cosmetic API inconsistency (no `$converged` boolean slot),
not a real convergence failure.

**Classification: test-gap** — no standard accessor for convergence status.

### Root cause summary

All numerical FAIL results (Findings 1–2 and the CI bounds) trace to a single root cause:
mixeff's Rust trust-region REML optimizer converges to a slightly different theta than lme4's
bobyqa/minqa, producing SEs and sigma that are self-consistent with mixeff's theta but
diverge from lme4's by more than the 1e-4 fixef-class tolerance. logLik and AIC/BIC agree
well (logLik diff 4.1e-06), confirming the objective is essentially equivalent; the gap is in
the precision of the marginal Hessian used to derive vcov.

---

## Outcome classification

| Cell | Outcome | Max abs diff (primary: Wald CI bounds) | Severity |
|------|---------|----------------------------------------|----------|
| inf-confint-wald | divergent | 1.40e-03 (Wald CI for Days) | major |

The feature exists and runs; the divergence is real and beyond tolerance on SE/CI bounds,
driven by a theta precision gap in the REML optimizer. logLik/fixef/theta/AIC are all within
their respective tolerances.
