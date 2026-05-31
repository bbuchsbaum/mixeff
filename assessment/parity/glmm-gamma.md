# Parity probe: glmm-gamma

**Cell:** glmm-gamma  
**Dataset:** Simulated Gamma (n=300, 20 groups x 15 obs, seed=42)  
**Formula:** `y ~ x + (1|g)`  
**Family:** `Gamma(link="log")`  
**Date run:** 2026-05-31  
**lme4:** 2.0.1 | **mixeff:** 0.1.0

---

## Script

See `glmm-gamma-probe.R` in this directory.

---

## Raw output (n=300, seed=42)

```
=== SESSION INFO ===
lme4 version: 2.0.1
mixeff version: 0.1.0

=== DATASET ===
nrow: 300   ncol: 3
y summary: min= 0.2909  mean= 4.7822  max= 27.6751

=== lme4 FIT ===
Warning: Model failed to converge with max|grad| = 0.00983153 (tol = 0.002, component 1)
lme4 wall-clock (seconds): 0.046

-- fixef --
(Intercept)           x
  1.4179466   0.3978094
-- SE --
(Intercept)           x
0.006776343 0.006570537
-- VarCorr --
 Groups Name        Std.Dev.
 g      (Intercept) 0.49575
-- sigma (dispersion) --
[1] 0.5474072
-- logLik --
'log Lik.' -645.2302 (df=4)
-- AIC --
[1] 1298.46
-- BIC --
[1] 1313.276

=== mixeff FIT ===
mixeff wall-clock (seconds): 0.013
fit_status: converged_interior
estimation_method: fast_pirls_profiled
objective_definition: profiled_glmm_deviance
response_constants: dropped

-- fixef --
(Intercept)           x
  1.4388065   0.3996769
-- SE --
(Intercept)           x
 0.04919685  0.03455463
-- VarCorr --
 group        name  variance  std_dev
     g (Intercept) 0.0231347 0.152101
Residual std. dev.: 0.549358
-- sigma (dispersion) --
[1] 0.5493575
-- logLik --
'log Lik.' -651.4268 (df=4)
-- AIC --
[1] 1310.854
-- BIC --
[1] 1325.669

=== NUMERICAL COMPARISON ===
fixef (Intercept)    lme4=1.41794661  mm=1.43880653  maxAbsDiff=2.086e-02  tol=1e-04  [EXCEEDS-TOL]
fixef x              lme4=0.39780940  mm=0.39967691  maxAbsDiff=1.868e-03  tol=1e-04  [EXCEEDS-TOL]
SE (Intercept)       lme4=0.00677634  mm=0.04919685  maxAbsDiff=4.242e-02  tol=1e-04  [EXCEEDS-TOL]
SE x                 lme4=0.00657054  mm=0.03455463  maxAbsDiff=2.798e-02  tol=1e-04  [EXCEEDS-TOL]
theta                lme4=0.49574739  mm=0.27687073  maxAbsDiff=2.189e-01  tol=1e-03  [EXCEEDS-TOL]
sigma (dispersion)   lme4=0.54740721  mm=0.54935753  maxAbsDiff=1.950e-03  tol=1e-04  [EXCEEDS-TOL]
VarCorr g variance   lme4=0.24576547  mm=0.02313472  maxAbsDiff=2.226e-01  tol=1e-03  [EXCEEDS-TOL]
logLik               lme4=-645.23025  mm=-651.42679  maxAbsDiff=6.197e+00  tol=1e-03  [EXCEEDS-TOL]
AIC                  lme4=1298.4605   mm=1310.8536   maxAbsDiff=1.239e+01  tol=2e-03  [EXCEEDS-TOL]
BIC                  lme4=1313.2756   mm=1325.6687   maxAbsDiff=1.239e+01  tol=2e-03  [EXCEEDS-TOL]
fitted max abs diff  maxAbsDiff=2.755e+00  tol=1e-04  [EXCEEDS-TOL]
ranef g max abs diff maxAbsDiff=1.779e-01  tol=1e-04  [EXCEEDS-TOL]

wall-clock elapsed   lme4=0.0460s  mm=0.0130s  ratio(mm/lme4)=0.28x  (mixeff ~3.5x faster)
```

---

## Analysis

### Root cause: objective convention difference (profiled vs. joint Laplace)

This is a **documented divergence by design**, not an optimizer bug. The contract is explicit
(`docs/glmm_support_contract.md`, "Approximation Semantics" and "Parity Claim Classes"):

- mixeff uses `fast_pirls_profiled` / `profiled_glmm_deviance` / `response_constants = dropped`.
  This is the **MixedModels.jl `fast=true` behavior** — profiled fast-PIRLS, covariance parameters
  optimized on the profiled objective with response normalising constants dropped.
- lme4 uses **joint Laplace estimation**, where `[β; θ]` are estimated on a joint deviance
  objective that includes response constants.

These are different approximations. The contract explicitly classifies such rows as
`documented_divergence`, not `release_blocking_parity`. The divergences below are entirely
consistent with this approximation-semantics difference.

### Quantity-by-quantity table

| Quantity | lme4 (joint Laplace) | mixeff (fast PIRLS) | maxAbsDiff | Tol | Status | Classification |
|---|---|---|---|---|---|---|
| fixef (Intercept) | 1.41795 | 1.43881 | 2.09e-2 | 1e-4 | exceeds-tol | documented_divergence |
| fixef x | 0.39781 | 0.39968 | 1.87e-3 | 1e-4 | exceeds-tol | documented_divergence |
| SE (Intercept) | 0.00678 | 0.04920 | 4.24e-2 | 1e-4 | exceeds-tol | documented_divergence |
| SE x | 0.00657 | 0.03455 | 2.80e-2 | 1e-4 | exceeds-tol | documented_divergence |
| theta | 0.49575 | 0.27687 | 2.19e-1 | 1e-3 | exceeds-tol | documented_divergence |
| sigma (dispersion) | 0.54741 | 0.54936 | 1.95e-3 | 1e-4 | exceeds-tol | documented_divergence |
| VarCorr g variance | 0.24577 | 0.02313 | 2.23e-1 | 1e-3 | exceeds-tol | documented_divergence |
| logLik | -645.230 | -651.427 | 6.20e+0 | 1e-3 | exceeds-tol | response-constant convention difference |
| AIC | 1298.46 | 1310.85 | 12.39 | 2e-3 | exceeds-tol | response-constant convention difference |
| BIC | 1313.28 | 1325.67 | 12.39 | 2e-3 | exceeds-tol | response-constant convention difference |
| fitted (max) | — | — | 2.76 | 1e-4 | exceeds-tol | follows from RE gap |
| ranef g (max) | — | — | 1.78e-1 | 1e-4 | exceeds-tol | follows from theta gap |
| convergence | FAILED (grad=0.0098) | converged_interior | — | — | lme4 not clean | see note |
| speed | 0.046s | 0.013s | ratio=0.28x | — | ~3.5x faster | passes |

### Notes on individual quantities

**fixef and SE divergences:** The profiled-PIRLS approximation integrates out θ differently than
joint Laplace. On this dataset the random-effect variance is substantially different between the
two approximations (0.246 vs 0.023), which pulls the intercept and propagates into standard errors.
The SE ratio (~7x larger in mixeff) reflects that mixeff's SEs are computed from the working
Hessian on the profiled objective with the much-shrunk random-effect variance, while lme4's SEs
come from the joint Hessian at a different θ. This is an approximation-gap artefact, not an SE
calculation bug.

**theta / VarCorr gap:** The profiled fast-PIRLS path finds a very different random-effect variance
(0.023) vs. lme4's joint Laplace (0.246). This is the central expression of the objective-convention
difference. On this simulated dataset (true RE sd=0.5, shape=4) the profiled path is substantially
under-estimating the RE variance. This is a known limitation of the fast-PIRLS approximation: it is
expected to be less accurate for inference when the profiled approximation is stressed. The contract
explicitly notes "it can be less accurate for inference when the profiled approximation is stressed,
especially overdispersed Poisson/binomial models." Gamma with overdispersion falls in this category.

**logLik / AIC / BIC:** The ~12.4 unit AIC/BIC gap is almost entirely the response-constant
convention (`response_constants = dropped` vs. lme4's `included`). The contract mandates that
`response_constants` is always reported in the artifact, which it is. Objectives are not
directly comparable across conventions; the contract labels this a convention difference,
not a fit failure. The remaining ~6 unit logLik gap above the constant difference reflects the
different θ estimates.

**lme4 convergence warning:** lme4 itself did not fully converge on this dataset (bobyqa
max|grad|=0.00983, tol=0.002). `nlminbwrap` converges and gives nearly identical estimates
(fixef≈1.4175/0.3977, theta≈0.4964, logLik≈-645.2304), confirming the bobyqa point is
essentially correct. lme4's gradient check is noisy for Gamma GLMM; the underlying estimates
are stable across optimizers.

**Speed:** mixeff is ~3.5x faster (0.013s vs 0.046s) despite the fundamental algorithmic
difference. The profiled fast-PIRLS path is designed to be faster; this is a design win.

### Scale-up check (n=1500, seed=123)

On a larger dataset (50 groups × 30 obs/group), the fixef gap narrows somewhat (Intercept
diff ~0.006, x diff ~0.004) but the theta gap persists (lme4=0.524, mm=0.384, diff=0.140).
VarCorr g: lme4=0.274, mm=0.036. lme4 also fails to converge on this dataset (max|grad|=0.013).
The pattern is consistent: profiled fast-PIRLS systematically under-estimates the random-effect
variance for Gamma GLMM relative to joint Laplace, and the gap does not close with more data.

---

## Classification summary

**Outcome: divergent**

All numerical divergences are **documented_divergence** (approximation-semantics difference:
profiled fast-PIRLS vs. joint Laplace), or **response-constant convention differences** (logLik/AIC/BIC).
Neither class is a `release_blocking_parity` row per `docs/glmm_support_contract.md`. The Gamma
family with log link is in the certified support list; mixeff fits it without error.

However, two concerns deserve mention:

1. **Severity: major** — The RE variance estimate is ~10x smaller than lme4's on these datasets.
   While the profiled objective is by-design different, a user comparing mixeff output to lme4 for
   a Gamma GLMM will see very different random-effect structure. The current `mm_glmm` print/summary
   does not proactively warn that `fast_pirls_profiled` may produce substantially different RE
   variance estimates than lme4 for Gamma models. An honest diagnostic note at fit time would
   satisfy the "no silent surgery" principle.

2. **Missing transparency token:** The artifact correctly records `estimation_method`,
   `objective_definition`, and `response_constants`. But the R-side `glmm()` caller receives no
   visible warning or message that the profiled approximation may be less accurate for Gamma —
   the only signal is `fit_status: converged_interior`. The contract says wrappers should not need
   to guess model semantics. Adding a `mm_fit_note` when `estimation_method == "fast_pirls_profiled"`
   and `family == "gamma"` would close this gap without blocking the fit.

**Severity: major** — not a bug in the optimizer or SE calculation, but the absence of a
user-visible warning for a known approximation gap that produces ~10x RE variance discrepancy
relative to lme4 on Gamma models. The fit does not fail; it just silently diverges from lme4
in ways that matter for inference.
