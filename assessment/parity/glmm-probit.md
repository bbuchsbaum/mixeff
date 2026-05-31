# Parity Assessment: glmm-probit

**Cell:** glmm-probit  
**Model:** `y ~ x + (1|g)`, `family = binomial(link = "probit")`  
**Dataset:** Simulated Bernoulli, N=200, 20 groups of 10 (seed 42, true b0=0.3, b1=0.8, RE SD=0.6)  
**Date:** 2026-05-31  
**mixeff method:** `pirls_profiled` (default)

---

## Raw Script Output

```
=== mixeff fit class ===
[1] "mm_glmm"     "mm_fit"      "mm_compiled"

=== fixef lme4 ===
(Intercept)           x
 0.08396309  0.81562899
=== fixef mixeff ===
(Intercept)           x
 0.07875249  0.77150495
=== fixef |diff| ===
(Intercept)           x
0.005210596 0.044124041
MAX fixef |diff|: 0.04412404

=== SE lme4 ===
(Intercept)           x
  0.1764690   0.1435798
=== SE mixeff ===
(Intercept)           x
  0.1572643   0.1135006
=== SE |diff| ===
(Intercept)           x
 0.01920469  0.03007916
MAX SE |diff|: 0.03007916

=== theta lme4 ===
g.(Intercept)
    0.6403111
=== theta mixeff ===
[1] 0.625
=== theta |diff| ===
g.(Intercept)
   0.01531115
MAX theta |diff|: 0.01531115

=== VarCorr lme4 ===
 Groups Name        Std.Dev.
 g      (Intercept) 0.64031

=== VarCorr mixeff ===
Variance components:
 group        name variance std_dev correlation note
     g (Intercept) 0.390625   0.625

=== sigma lme4: 1  ===
=== sigma mixeff: 1 ===
sigma |diff|: 0

=== logLik lme4: -114.7077 ===
=== logLik mixeff: -114.7568 ===
logLik |diff|: 0.04905932

=== AIC lme4: 235.4154  mixeff: 235.5135  |diff|: 0.09811863 ===
=== BIC lme4: 245.3104  mixeff: 245.4085  |diff|: 0.09811863 ===

=== ranef (first 5) lme4 ===
[1] -0.5436232  0.2617379  0.2812325  0.5949613 -0.6615291 -0.7991445
=== ranef (first 5) mixeff ===
[1] -0.5178878  0.2539347  0.2704184  0.5723032 -0.6446737 -0.7691475
MAX ranef |diff|: 0.03062494

=== fitted (first 5) lme4 ===
        1         2         3         4         5         6
0.7449023 0.1787226 0.4350695 0.5225366 0.4483120 0.2924579
=== fitted (first 5) mixeff ===
[1] 0.7318988 0.1908406 0.4368423 0.5195887 0.4493751 0.3011796
MAX fitted |diff|: 0.01767683

=== lme4 convergence warnings ===
n messages: 0
=== mixeff fit_status ===
converged_interior

=== Timing ===
lme4   elapsed: 0.047 s
mixeff elapsed: 0.015 s
speed ratio (lme4/mixeff): 3.133333
```

---

## Quantity-by-Quantity Analysis

| Quantity        | lme4 value          | mixeff value       | Max \|diff\|   | Tolerance | Status          |
|----------------|---------------------|--------------------|----------------|-----------|-----------------|
| fixef intercept | 0.08396             | 0.07875            | 0.00521        | 1e-4      | **BEYOND TOL**  |
| fixef x         | 0.81563             | 0.77150            | 0.04412        | 1e-4      | **BEYOND TOL**  |
| SE intercept    | 0.17647             | 0.15726            | 0.01920        | (1e-4)    | **BEYOND TOL**  |
| SE x            | 0.14358             | 0.11350            | 0.03008        | (1e-4)    | **BEYOND TOL**  |
| theta           | 0.64031             | 0.62500            | 0.01531        | 1e-3      | **BEYOND TOL**  |
| sigma           | 1.0                 | 1.0                | 0.0            | 1e-4      | within-tol      |
| logLik          | -114.7077           | -114.7568          | 0.04906        | 1e-3      | **BEYOND TOL**  |
| AIC             | 235.4154            | 235.5135           | 0.09812        | (derived) | **BEYOND TOL**  |
| BIC             | 245.3104            | 245.4085           | 0.09812        | (derived) | **BEYOND TOL**  |
| ranef max diff  | —                   | —                  | 0.03062        | —         | **BEYOND TOL**  |
| fitted max diff | —                   | —                  | 0.01768        | —         | **BEYOND TOL**  |
| convergence     | 0 warnings          | converged_interior | —              | —         | OK (both converged) |

---

## Interpretation

### What worked
- mixeff accepted `binomial(link = "probit")` without error — the certified surface is correctly plumbed.
- `fit_status` is `converged_interior`; lme4 reports 0 convergence warnings. Both optimizers agree the model converged.
- `sigma` is exactly 1.0 for both (correct for binomial GLMM — dispersion is fixed).
- mixeff is **3.1× faster** than lme4 on this N=200 / 20-group dataset (0.015 s vs 0.047 s).
- All accessor methods dispatched without error: `fixef()`, `vcov()`, `VarCorr()`, `ranef()`, `$theta`, `$logLik`, `$AIC`, `$BIC`, `$fitted`.

### What diverged
All numerical quantities exceed the stated tolerances. The divergence is **systematic**, not random:

- **fixef**: mixeff intercept and slope are both shrunk toward zero relative to lme4 (intercept −6%, slope −5.4%). This suggests the Laplace approximation in the `pirls_profiled` path is less accurate than lme4's Laplace approximation for this probit model, or a subtle difference in the probit link implementation in the Rust engine.
- **SE / vcov**: mixeff SEs are smaller than lme4's (intercept −11%, slope −21%). Underestimated SEs are concerning — if the Hessian-based covariance computation in the Rust engine is not accounting for the full uncertainty correctly, inference (Wald z-tests, CIs) will be anti-conservative.
- **theta**: 0.625 vs 0.640. Small but above the 1e-3 tolerance.
- **logLik**: −114.757 vs −114.708, diff = 0.049. This is the key diagnostic: lme4 finds a slightly higher (better) likelihood at its optimum, which is consistent with the fixef and theta differences — mixeff converges to a nearby but suboptimal point.
- **ranef / fitted**: Downstream from the fixef/theta differences; within ~3% in absolute terms but systematically offset.

### Root-cause hypothesis
The `logLik` gap strongly suggests mixeff is **not finding the same optimum** as lme4. Possible causes (in order of likelihood):
1. Different optimizer convergence criteria or step tolerances in the Rust `pirls_profiled` path causing premature termination.
2. Probit-specific numerical issue in the PIRLS working weights (the probit inverse-link derivative `dnorm(eta)` is flatter than logistic's, making PIRLS more sensitive to step size).
3. Different parameterization of theta (Cholesky vs SD) causing the penalty surface to be traversed differently.

The SE underestimation relative to lme4 is a separate concern — it may indicate the Hessian at the mixeff optimum is incorrectly computed, or the observed-information vs expected-information choice differs.

### Classification
- **fixef divergence**: `in-scope-missing` — probit fit accuracy is below the 1e-4 tolerance; the optimizer is converging to a suboptimal point.
- **SE divergence**: `in-scope-missing` — Hessian-based SE computation for probit link appears less accurate than lme4.
- **logLik divergence**: `in-scope-missing` — flows from the optimizer gap.
- **Sigma, convergence, API surface**: `works`.
- **Speed**: `works` (mixeff is faster).

### Severity
**major** — fixef and SE are the primary quantities users rely on for inference. A ~5% fixef difference and ~20% SE underestimation for a probit model means Wald-z statistics and p-values will be noticeably wrong. This is above the stated 1e-4 fixef tolerance. The issue is specific to the probit link (logit should be checked separately); it does not affect model acceptance or crash the package.

---

## Script location
`assessment/parity/glmm-probit-probe.R`
