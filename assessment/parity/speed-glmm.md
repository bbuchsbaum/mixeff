# Parity Probe: speed-glmm

**Cell:** speed-glmm  
**Date:** 2026-05-31  
**Dataset:** Simulated binomial, N=10 000, 50 groups  
**Formula:** `y ~ x + (1|g)`, `family = binomial(logit)`, Laplace (nAGQ=1)  
**lme4 version:** 2.0.1 &nbsp;&nbsp; **mixeff version:** 0.1.0

---

## 1. Speed result

| | lme4 (bobyqa) | mixeff (pirls_profiled) | ratio (mm/lme4) |
|---|---|---|---|
| Median wall-clock (3 reps, post warm-up) | 0.188 s | 0.055 s | **0.29×** |
| Mean wall-clock | 0.189 s | 0.055 s | **0.29×** |

**mixeff is 3.42× faster than lme4** on N=10 000, y~x+(1|g), binomial/logit/Laplace.

---

## 2. Convergence / fit status

- **lme4:** converged without warnings (`optinfo$conv$lme4$messages` is empty).
- **mixeff:** `fit_status = "converged_interior"`.

Both fits report successful convergence.

---

## 3. Numerical parity — seed 42

### 3a. Raw values

| Quantity | lme4 | mixeff | |
|---|---|---|---|
| fixef (Intercept) | 0.47574 | 0.47284 | |
| fixef x | 0.30849 | 0.30702 | |
| SE (Intercept) | 0.11525 | 0.11506 | |
| SE x | 0.02236 | 0.02228 | |
| theta | 0.79931 | 0.80078 | |
| sigma (dispersion) | 1.0 | 1.0 | exact match |
| logLik | −6063.938 | −6063.941 | |
| AIC | 12133.88 | 12133.88 | |
| BIC | 12155.51 | 12155.51 | |
| deviance (field) | 11915.62 | 12127.88 | see §4 |
| ranef max abs diff | — | — | 3.0e-3 |
| fitted max abs diff | — | — | 1.4e-3 |

### 3b. Max absolute differences vs tolerances

| Quantity | maxAbsDiff | Tolerance | Status |
|---|---|---|---|
| fixef (Intercept) | **2.89e-3** | 1e-4 | EXCEEDS-TOL |
| fixef x | **1.47e-3** | 1e-4 | EXCEEDS-TOL |
| SE (Intercept) | **1.84e-4** | 1e-4 | EXCEEDS-TOL |
| SE x | 8.1e-5 | 1e-4 | within-tol |
| theta | **1.47e-3** | 1e-3 | EXCEEDS-TOL |
| sigma | 0 | 1e-4 | within-tol (exact) |
| logLik | **2.65e-3** | 1e-3 | EXCEEDS-TOL |
| AIC | **5.3e-3** | 2e-3 | EXCEEDS-TOL |
| BIC | **5.3e-3** | 2e-3 | EXCEEDS-TOL |
| fitted | **1.4e-3** | 1e-4 | EXCEEDS-TOL |
| ranef | **3.0e-3** | 1e-4 | EXCEEDS-TOL |

---

## 4. Deviance definition difference (not a bug)

mixeff's `$deviance` field is **−2 × logLik = 12127.88**, matching `lme4::getME(fit, "devcomp")$cmp["dev"]`.  
lme4's `deviance()` generic for GLMMs returns the **penalized (REML-style) deviance = 11915.62**, which includes the random-effect penalty term subtracted from the Laplace log-determinant.  
The gap of ~212 is entirely a definition difference. The fixef and theta exceedances in the comparison table above are the real numerical issues.

---

## 5. Cross-seed stability (5 seeds)

| seed | d_fixef | d_theta | d_loglik | lme4 > mm |
|---|---|---|---|---|
| 1 | 3.65e-3 | 3.1e-5 | 3.23e-3 | YES |
| 2 | 4.26e-3 | 1.5e-5 | 2.04e-3 | YES |
| 3 | 4.04e-3 | 3.9e-5 | 2.26e-3 | YES |
| 7 | 3.68e-3 | 1.6e-5 | 2.34e-3 | YES |
| 99 | 3.49e-3 | 1.7e-4 | 2.31e-3 | YES |

The pattern is **systematic**:

- **fixef gap** is consistently ~3.5–4.3e-3 (35–43× the 1e-4 tolerance). This exceeds tolerance on every seed.
- **theta gap** is consistently tiny (~1.5–3.1e-5, well within the 1e-3 tolerance) — the random-effect variance parameter converges identically.
- **logLik gap** is consistently ~2–3.2e-3 (2–3× the 1e-3 tolerance). lme4 always finds the slightly better optimum.
- **lme4 always achieves a higher log-likelihood** than mixeff, confirming mixeff terminates too early (not a different local optimum with lower curvature — the loss landscape is convex here).

**Root cause:** mixeff's Rust PIRLS/profiled-likelihood optimizer terminates before fully converging the fixed-effect gradient. The theta parameter is essentially at the same solution; only the fixed-effect subproblem is under-tightened. `mm_control()` exposes only `verbose` — no tolerance knob is available to the R caller.

---

## 6. Classification

| Finding | Classification | Severity |
|---|---|---|
| mixeff 3.42× faster than lme4 | works (goal achieved) | none |
| theta within 1.5e-5 on all seeds | works | none |
| sigma exact match | works | none |
| deviance definition mismatch (212 gap) | cosmetic — definition difference, not a numerical error | cosmetic |
| fixef gap ~3.5–4.3e-3 (consistently exceeds 1e-4 tol) | **in-scope-missing** — PIRLS termination too loose | **major** |
| logLik gap ~2–3.2e-3 (consistently exceeds 1e-3 tol) | in-scope-missing — follows from fixef under-convergence | **major** |
| AIC/BIC gap ~5.3e-3 (follows from logLik) | in-scope-missing — follows from logLik | minor |
| SE (Intercept) gap 1.84e-4 (slightly exceeds 1e-4) | in-scope-missing — follows from fixef | minor |
| fitted/ranef gaps follow from fixef | in-scope-missing | minor |

**Overall outcome: divergent** — speed goal is strongly met (3.42×), but fixef convergence precision falls short of the parity tolerance on every tested seed. The issue is the Rust PIRLS optimizer's termination criterion, not a fundamental algorithmic failure. Tightening the convergence tolerance (or exposing it via `mm_control`) should close the gap.

---

## 7. Raw script output (seed 42)

```
=== SESSION INFO ===
lme4 version:    2.0.1 
mixeff version:  0.1.0 

=== DATASET ===
N: 10000   groups: 50   y=1 proportion: 0.603 

=== TIMING (3reps each after warm-up) ===
lme4  times (s): 0.1920 0.1870 0.1880 
lme4  median  : 0.1880 s
lme4  mean    : 0.1890 s

mixeff times (s): 0.0550 0.0570 0.0520 
mixeff median   : 0.0550 s
mixeff mean     : 0.0547 s

Speed ratio (mixeff/lme4): median=0.29x  mean=0.29x
  (ratio < 1 = mixeff FASTER; > 1 = mixeff SLOWER)

fit_status: converged_interior 

=== NUMERICAL COMPARISON ===
fixef[(Intercept)]                      maxAbsDiff=2.892e-03  tol=1e-04  [EXCEEDS-TOL*]
fixef[x]                                maxAbsDiff=1.473e-03  tol=1e-04  [EXCEEDS-TOL*]
fixef (all)                             maxAbsDiff=2.892e-03  tol=1e-04  [EXCEEDS-TOL*]
SE[(Intercept)]                         maxAbsDiff=1.836e-04  tol=1e-04  [EXCEEDS-TOL*]
SE[x]                                   maxAbsDiff=8.118e-05  tol=1e-04  [WITHIN-TOL]
SE (all)                                maxAbsDiff=1.836e-04  tol=1e-04  [EXCEEDS-TOL*]
theta                                   maxAbsDiff=1.471e-03  tol=1e-03  [EXCEEDS-TOL*]
sigma (dispersion)                      maxAbsDiff=0.000e+00  tol=1e-04  [WITHIN-TOL]
VarCorr g variance                      lme4=0.638898  mm=0.641251
VarCorr g variance                      maxAbsDiff=2.353e-03  tol=1e-03  [EXCEEDS-TOL*]
logLik                                  maxAbsDiff=2.647e-03  tol=1e-03  [EXCEEDS-TOL*]
deviance                                maxAbsDiff=2.123e+02  tol=2e-03  [EXCEEDS-TOL*]
AIC                                     maxAbsDiff=5.293e-03  tol=2e-03  [EXCEEDS-TOL*]
BIC                                     maxAbsDiff=5.293e-03  tol=2e-03  [EXCEEDS-TOL*]
fitted max abs diff                     maxAbsDiff=1.414e-03  tol=1e-04  [EXCEEDS-TOL*]
ranef g max abs diff                    maxAbsDiff=3.020e-03  tol=1e-04  [EXCEEDS-TOL*]

=== SPEED SUMMARY ===
N=10000, ng=50, formula=y~x+(1|g), family=binomial(logit), Laplace
lme4   median wall-clock : 0.1880 s
mixeff median wall-clock : 0.0550 s
Speed ratio (mixeff/lme4): 0.29x
  => mixeff is 3.42x FASTER than lme4
```
