# Parity probe: lmm-contrasts

**Cell:** lmm-contrasts  
**Dataset:** Simulated, n=120 (12 groups × 10 obs/group), 4-level factor `fac`, random intercept `(1|g)`  
**Formula:** `y ~ fac + (1|g)`  
**REML:** TRUE  
**Date run:** 2026-05-31  
**lme4:** 2.0.1 | **lmerTest:** 3.2.1 | **mixeff:** 0.1.0

---

## Script

See `lmm-contrasts-probe.R` in this directory.

---

## Raw output (abridged)

```
=== SCENARIO 1: treatment contrasts (default) ===

lme4 coef names: (Intercept), facB, facC, facD
mm   coef names: (Intercept), fac: B, fac: C, fac: D

fixef[1] (Intercept)  lme4=6.17633711  mm=6.17633711  maxAbsDiff=2.576e-14  [WITHIN-TOL]
fixef[2] facB         lme4=1.84090979  mm=1.84090979  maxAbsDiff=4.663e-15  [WITHIN-TOL]
fixef[3] facC         lme4=-1.63019887 mm=-1.63019887 maxAbsDiff=3.109e-15  [WITHIN-TOL]
fixef[4] facD         lme4=3.09490182  mm=3.09490182  maxAbsDiff=6.217e-15  [WITHIN-TOL]
SE[1]-(4)            all WITHIN-TOL (max diff 7.3e-7)
vcov max abs diff     maxAbsDiff=6.593e-07  [WITHIN-TOL]
theta                 maxAbsDiff=3.075e-06  [WITHIN-TOL]
sigma                 maxAbsDiff=1.828e-07  [WITHIN-TOL]
VarCorr g variance    lme4=2.14923826  mm=2.14924629  maxAbsDiff=8.032e-06  [WITHIN-TOL]
logLik                maxAbsDiff=3.965e-11  [WITHIN-TOL]
AIC                   maxAbsDiff=7.930e-11  [WITHIN-TOL]
BIC                   maxAbsDiff=7.930e-11  [WITHIN-TOL]
fitted max abs diff   maxAbsDiff=4.128e-07  [WITHIN-TOL]
ranef g max abs diff  maxAbsDiff=4.128e-07  [WITHIN-TOL]
wall-clock  lme4=0.029s  mm=0.011s  ratio=0.38x  (mixeff ~2.6x faster)

=== SCENARIO 2: sum contrasts (contr.sum) ===

lme4 (contr.sum via factor attribute) fixef:
  (Intercept)=7.0027403  fac1=-0.8264032  fac2=1.0145066  fac3=-2.4566021
  logLik=-187.7136

lme4 (contr.sum via contrasts= arg) fixef: SAME as above (correct)

mixeff fit (contr.sum via factor attribute) fixef:
  (Intercept)=6.176337  fac: B=1.840910  fac: C=-1.630199  fac: D=3.094902
  logLik=-186.3273   <-- WRONG: same as treatment-contrast fit

mixeff::lmm has 'contrasts' formal argument: FALSE
mixeff + contrasts= arg: ERROR — unused argument (contrasts = list(fac = "contr.sum"))

fixef(sum)[1] (Intercept)  maxAbsDiff=8.264e-01  [EXCEEDS-TOL]
fixef(sum)[2] fac1         maxAbsDiff=2.667e+00  [EXCEEDS-TOL]
fixef(sum)[3] fac2         maxAbsDiff=2.645e+00  [EXCEEDS-TOL]
fixef(sum)[4] fac3         maxAbsDiff=5.552e+00  [EXCEEDS-TOL]
SE(sum) max abs diff       maxAbsDiff=1.110e-01  [EXCEEDS-TOL]
logLik(sum)                maxAbsDiff=1.386e+00  [EXCEEDS-TOL]

lme4 logLik invariant check: treatment=-186.327276  sum=-187.713570  diff=1.386e+00
  (lme4 logLik differs between coding schemes — expected: different log-likelihood
   due to unequal cell sizes in the simulated data; REML penalises differently)
mm  logLik: treatment=-186.327276  sum=-186.327276  diff=0.000e+00
  (mixeff fit is IDENTICAL for sum and treatment — confirms it ignored the
   contrasts attribute)

=== SCENARIO 3: Helmert contrasts (contr.helmert) ===

lme4 fixef: (Intercept)=7.003  fac1=0.920  fac2=-0.850  fac3=0.756
            logLik=-189.505

mixeff fixef: (Intercept)=6.177  fac: B=1.841  fac: C=-1.630  fac: D=3.095
              logLik=-186.327   <-- WRONG: same as treatment-contrast fit

fixef(helm)[1]-(4): all EXCEEDS-TOL (diffs 0.83 to 2.34)
logLik(helm): maxAbsDiff=3.178  [EXCEEDS-TOL]

=== SCENARIO 4: 5-level factor (treatment, default) ===

All fixef WITHIN-TOL. logLik WITHIN-TOL. sigma WITHIN-TOL.
```

---

## Analysis

### Quantity-by-quantity

#### Scenario 1: Treatment contrasts (default, 4 levels)

| Quantity | lme4 | mixeff | maxAbsDiff | Tol | Status |
|---|---|---|---|---|---|
| fixef (Intercept) | 6.17634 | 6.17634 | 2.6e-14 | 1e-4 | within-tol |
| fixef facB | 1.84091 | 1.84091 | 4.7e-15 | 1e-4 | within-tol |
| fixef facC | -1.63020 | -1.63020 | 3.1e-15 | 1e-4 | within-tol |
| fixef facD | 3.09490 | 3.09490 | 6.2e-15 | 1e-4 | within-tol |
| SE (all 4) | — | — | 7.3e-7 max | 1e-4 | within-tol |
| vcov (all entries) | — | — | 6.6e-7 | 1e-4 | within-tol |
| theta | 1.496268 | 1.496272 | 3.1e-6 | 1e-3 | within-tol |
| sigma | 0.979789 | 0.979789 | 1.8e-7 | 1e-4 | within-tol |
| VarCorr g var | 2.149238 | 2.149246 | 8.0e-6 | 1e-3 | within-tol |
| logLik | -186.3273 | -186.3273 | 4.0e-11 | 1e-3 | within-tol |
| AIC | 384.6546 | 384.6546 | 7.9e-11 | 2e-3 | within-tol |
| BIC | 401.3795 | 401.3795 | 7.9e-11 | 2e-3 | within-tol |
| fitted (max) | — | — | 4.1e-7 | 1e-4 | within-tol |
| ranef g (max) | — | — | 4.1e-7 | 1e-4 | within-tol |
| convergence | no warnings | converged_interior | — | — | works |
| speed | 0.029s | 0.011s | ratio=0.38x | — | ~2.6x faster |

#### Scenario 2: Sum contrasts (contr.sum)

**Root cause:** mixeff ignores the `contrasts` attribute on the factor in the data frame. It uses treatment (dummy) coding regardless of the `contrasts()` attribute set on the factor column. The fit is numerically identical to the treatment-contrast fit.

lme4 correctly uses the sum-contrast design matrix when the attribute is pre-set (or passed via `contrasts=`). The fitted model is a different parameterization with different fixef and different logLik (note: logLik differs between coding schemes here because of unequal cell sizes and REML, so this is a genuine parameterization difference, not just a re-labeling).

Additionally, `mixeff::lmm()` does **not** accept a `contrasts=` argument. Calling `lmm(..., contrasts = list(fac = "contr.sum"))` produces an R error: `unused argument (contrasts = list(fac = "contr.sum"))`.

| Quantity | lme4 | mixeff | maxAbsDiff | Tol | Status |
|---|---|---|---|---|---|
| fixef (Intercept) | 7.003 | 6.177 | 8.3e-1 | 1e-4 | EXCEEDS-TOL |
| fixef fac1/fac:B | -0.826 | 1.841 | 2.7e+0 | 1e-4 | EXCEEDS-TOL |
| fixef fac2/fac:C | 1.015 | -1.630 | 2.6e+0 | 1e-4 | EXCEEDS-TOL |
| fixef fac3/fac:D | -2.457 | 3.095 | 5.6e+0 | 1e-4 | EXCEEDS-TOL |
| SE (max) | 0.433 | 0.454 | 1.1e-1 | 1e-4 | EXCEEDS-TOL |
| logLik | -187.714 | -186.327 | 1.4e+0 | 1e-3 | EXCEEDS-TOL |
| sigma | within-tol | — | 1.8e-7 | 1e-4 | within-tol |
| contrasts= arg | accepted | error | — | — | **MISSING** |

#### Scenario 3: Helmert contrasts (contr.helmert)

Same root cause as Scenario 2. mixeff ignores the Helmert contrasts attribute and fits treatment coding silently. No error or diagnostic is emitted.

| Quantity | lme4 | mixeff | maxAbsDiff | Tol | Status |
|---|---|---|---|---|---|
| fixef[1] (Intercept) | 7.003 | 6.177 | 8.3e-1 | 1e-4 | EXCEEDS-TOL |
| fixef[2]-[4] | various | treatment vals | 0.78–2.34 | 1e-4 | EXCEEDS-TOL |
| logLik | -189.505 | -186.327 | 3.2e+0 | 1e-3 | EXCEEDS-TOL |
| sigma | within-tol | — | 1.8e-7 | 1e-4 | within-tol |

#### Scenario 4: 5-level treatment contrasts

All quantities within tolerance. Confirms mixeff handles multi-level factors correctly when treatment coding is used (the default).

---

## Coefficient name mismatch (cosmetic)

lme4 names factor coefficients `"facB"`, `"facC"`, `"facD"` (variable name glued to level). mixeff names them `"fac: B"`, `"fac: C"`, `"fac: D"` (variable name + colon-space + level). This is a cosmetic difference in the coefficient name format. It does not affect numerical values but will break any downstream code that looks up coefficients by name using lme4-style names (e.g., a custom contrast matrix built with `c("facB" = 1)`), and it will mismatch with `emmeans` or `multcomp` workflows that rely on coefficient name conventions.

---

## Classification

| Gap | Classification | Severity |
|---|---|---|
| `contrasts` attribute on factor is silently ignored — mixeff always uses treatment coding | **in-scope-missing** | **major** |
| No `contrasts=` argument in `lmm()` (lme4 has it) | **in-scope-missing** | **major** |
| No diagnostic/warning when a non-default contrasts attribute is found on a factor | **in-scope-missing** | **major** |
| Coefficient names use `"fac: B"` vs lme4's `"facB"` | **in-scope-missing** | **minor** |
| Treatment-contrast fit: all quantities within tolerance | works | none |
| 5-level factor (treatment): all quantities within tolerance | works | none |

---

## Overall verdict

**Outcome: divergent**

Scenario 1 (treatment coding, default) and Scenario 4 (5-level treatment) are fully within tolerance and show no defects. However, the contrasts-parity cell specifically tests whether non-default factor coding is respected. It is not.

When `contrasts()` is pre-set on a factor column (the standard R idiom for sum or Helmert coding), mixeff silently discards the attribute and fits treatment contrasts. The resulting fixef, SEs, and logLik are all wrong for the requested parameterization. No error, warning, or diagnostic is emitted — this violates the "no silent surgery" principle stated in CLAUDE.md.

Additionally, `lmm()` lacks a `contrasts=` argument, so the lme4-style `contrasts = list(fac = "contr.sum")` call-site pattern is an immediate error.

**Severity: major** — Silent wrong answer for a very common workflow (sum contrasts for ANOVA-style models). The user gets a converged fit with plausible-looking numbers that are a completely different parameterization than requested.
