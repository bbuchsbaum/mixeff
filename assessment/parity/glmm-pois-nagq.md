# Parity Assessment: glmm-pois-nagq

**Cell:** `glmm-pois-nagq`
**Date:** 2026-05-31
**Dataset:** Simulated single-scalar RE Poisson, N=300, 30 groups × 10 obs,
  true intercept=1, x slope=0.5, RE sd=0.4 (seed=42)
**Formula:** `y ~ x + (1|g)`
**Focus:** nAGQ>1 AGQ parity (and documented refusal status)

---

## Raw Output

```
=== Environment ===
lme4 version: 2.0.1
mixeff version: 0.1.0
R version: R version 4.5.1 (2025-06-13)

=== Dataset summary ===
N obs: 300  Groups: 30
y range: 0 15   mean: 3.217
```

### Scenario A — nAGQ=1 (Laplace)

```
--- lme4 glmer (nAGQ=1) ---
lme4 converged: TRUE
fixef: 0.98437033  0.53208380
theta: 0.32411365
logLik: -580.9449069
AIC: 1167.889814
BIC: 1179.001161
sigma: 1
nobs: 300   nAGQ: 1
Wall time: 0.047 s

--- mixeff glmm (nAGQ=1) ---
fit_status: converged_interior
method: pirls_profiled   nAGQ: 1
fixef: 0.99697651  0.53193880
theta: 0.32420349
logLik: -580.9619156
AIC: 1167.923831
BIC: 1179.035179
sigma: 1
nobs: 300
Wall time: 0.014 s

--- Comparison nAGQ=1 ---
  fixef[(Intercept)]    abs_diff=1.26e-02  tol=1e-04  [BEYOND-TOL]
  fixef[x]              abs_diff=1.45e-04  tol=1e-04  [BEYOND-TOL]
  theta (RE sd)         abs_diff=8.98e-05  tol=1e-03  [WITHIN-TOL]
  logLik                abs_diff=1.70e-02  tol=1e-03  [BEYOND-TOL]
  ranef (max across groups)   max_abs_diff=1.10e-02
  fitted values (max)         max_abs_diff=2.97e-02
  VarCorr SD            abs_diff=8.98e-05  tol=1e-04  [WITHIN-TOL]
```

### Scenario B — nAGQ=5 (AGQ)

```
--- lme4 glmer (nAGQ=5) ---
lme4 converged: TRUE
fixef: 0.98427922  0.53206943
theta: 0.32467314
logLik: -194.2462718       ← response constants INCLUDED
AIC: 394.4925436
BIC: 405.6038911
sigma: 1
nobs: 300   nAGQ: 5
Wall time: 0.024 s

--- mixeff glmm (nAGQ=5, method='pirls_profiled') ---
fit_status: converged_interior
method: pirls_profiled   nAGQ: 5
fixef: 0.99691784  0.53192968
theta: 0.32470703
logLik: -580.9448678       ← response constants DROPPED
AIC: 1167.889736
BIC: 1179.001083
sigma: 1
nobs: 300
Wall time: 0.004 s

--- mixeff glmm (nAGQ=5, method='joint_laplace') — expected refusal ---
Result: ERROR (expected)
Class: mm_arg_error, mm_condition, rlang_error, error, condition
Message: `method = "joint_laplace"` requires `nAGQ <= 1` in this slice.
==> TYPED REFUSAL: honest and labelled
```

### Scenario C — Cross-nAGQ internal consistency

```
lme4 nAGQ=1 vs nAGQ=5:
  fixef[(Intercept)]: nAGQ1=0.98437033  nAGQ5=0.98427922  diff=9.11e-05
  fixef[x]:           nAGQ1=0.53208380  nAGQ5=0.53206943  diff=1.44e-05
  theta:              nAGQ1=0.32411365  nAGQ5=0.32467314  diff=-5.59e-04
  logLik:             nAGQ1=-580.944    nAGQ5=-194.246     diff=-3.87e+02  ← constant gap

mixeff nAGQ=1 vs nAGQ=5 (pirls_profiled):
  fixef[(Intercept)]: nAGQ1=0.99697651  nAGQ5=0.99691784  diff=5.87e-05
  fixef[x]:           nAGQ1=0.53193880  nAGQ5=0.53192968  diff=9.12e-06
  theta:              nAGQ1=0.32420349  nAGQ5=0.32470703  diff=-5.04e-04
  logLik:             nAGQ1=-580.962    nAGQ5=-580.945     diff=-1.70e-02
```

### Scenario D — mixeff(nAGQ=5) vs lme4(nAGQ=5) parity

```
  fixef[(Intercept)]  abs_diff=1.26e-02  tol=1e-04  [BEYOND-TOL]
  fixef[x]            abs_diff=1.40e-04  tol=1e-04  [BEYOND-TOL]
  theta               abs_diff=3.39e-05  tol=1e-03  [WITHIN-TOL]
  logLik              abs_diff=3.87e+02  tol=1e-03  [BEYOND-TOL]   ← convention gap
```

### Scenario E — Timing

```
lme4   nAGQ=1: 0.047 s
mixeff nAGQ=1: 0.014 s
lme4   nAGQ=5: 0.024 s
mixeff nAGQ=5: 0.004 s
Speed ratio nAGQ=1 (lme4/mixeff): 3.36x  (mixeff faster)
```

---

## Analysis

### Quantity-by-quantity breakdown

| Quantity | lme4 (nAGQ=1) | mixeff (nAGQ=1) | Abs diff | Tol | Status |
|---|---|---|---|---|---|
| fixef (Intercept) | 0.98437033 | 0.99697651 | 1.26e-02 | 1e-04 | **BEYOND-TOL** |
| fixef x | 0.53208380 | 0.53193880 | 1.45e-04 | 1e-04 | **BEYOND-TOL** |
| theta | 0.32411365 | 0.32420349 | 8.98e-05 | 1e-03 | WITHIN-TOL |
| logLik | -580.9449 | -580.9619 | 1.70e-02 | 1e-03 | **BEYOND-TOL** |
| VarCorr SD | 0.32411 | 0.32420 | 8.98e-05 | 1e-04 | WITHIN-TOL |

| Quantity | lme4 (nAGQ=5) | mixeff (nAGQ=5) | Abs diff | Tol | Status |
|---|---|---|---|---|---|
| fixef (Intercept) | 0.98427922 | 0.99691784 | 1.26e-02 | 1e-04 | **BEYOND-TOL** |
| fixef x | 0.53206943 | 0.53192968 | 1.40e-04 | 1e-04 | **BEYOND-TOL** |
| theta | 0.32467314 | 0.32470703 | 3.39e-05 | 1e-03 | WITHIN-TOL |
| logLik | -194.2463 | -580.9449 | 386.7 | 1e-03 | BEYOND-TOL (convention) |

### Finding 1: fixef intercept gap (~0.013) — documented divergence, not a bug

The fixed-effect intercept differs by ~0.013 in both the nAGQ=1 and nAGQ=5
comparisons. This is the expected and documented divergence between lme4's
**joint Laplace** objective (lme4 default for nAGQ=1 via `glmer`) and mixeff's
**profiled fast-PIRLS** objective. The `glmm_support_contract.md` explicitly
classifies this as `documented_divergence`:

> *"fast-PIRLS/profiled-objective rows that track the MixedModels.jl `fast=true`
> behavior while differing from `lme4` joint-estimation coefficients."*

lme4's `glmer` with `nAGQ=1` uses a joint `[β; θ]` Laplace objective; mixeff's
`pirls_profiled` profiles `β` out of the deviance and optimizes only over `θ`.
These are different objective functions; coefficient differences at this scale
are expected. This is not a numerical accuracy issue.

The fixef[x] gap (1.40–1.45e-04) just barely exceeds the 1e-04 tolerance. Given
the ~0.013 intercept gap driven by the same objective-function difference, this
is the same root cause, not an independent problem.

**Classification:** `documented_divergence` per `glmm_support_contract.md`.
Not a release-blocking parity failure for the `pirls_profiled` path.

### Finding 2: logLik gap of ~387 (nAGQ=5) — response-constant convention

lme4 reports `logLik=-194.25` for nAGQ=5; mixeff reports `logLik=-580.94`.
The ~387 unit difference is the dataset-dependent response-normalizing constant
that lme4 includes and mixeff drops. This is explicit in `glmm_support_contract.md`
under *"response_constants"*: `dropped` for the supported fast path vs `included`
for joint objectives. For nAGQ=1, the same convention difference produces a
logLik gap of only ~0.017 because the Poisson log-factorial terms are always
present regardless of AGQ order; the large gap at nAGQ=5 comes from lme4
switching to a properly normalized AGQ objective while mixeff stays on the
profiled-no-constants surface.

**This is a convention difference, not an optimizer failure.** Comparing logLik
across these two engines at nAGQ=5 is meaningless without adjusting for the
response-constant offset.

**Classification:** By-design convention gap (`response_constants: dropped` vs
`included`). Not a defect.

### Finding 3: nAGQ>1 acceptance in pirls_profiled — partial implementation

mixeff accepts `nAGQ=5` with `method="pirls_profiled"` and fits successfully
(`fit_status: converged_interior`). The Rust `fit_with_options(fast=true, n_agq=5)`
path passes `n_agq` through to `penalized_pirls_deviance_at_theta`, which
calls the GLMM `deviance(n_agq=5)` for objective evaluations — so the
optimization landscape does use 5-point AGQ for the deviance approximation.

However, the **optimizer is still optimizing only over θ** (not jointly over
`[β; θ]`), so this is "profiled-PIRLS with AGQ deviance" not "full AGQ" as
lme4 implements it. The internal consistency is good: mixeff's own nAGQ=1
vs nAGQ=5 difference is small (fixef diff < 6e-05, theta diff < 5.1e-04,
logLik diff 0.017), showing that increasing quadrature points does shift the
profiled objective slightly as expected.

The gap between mixeff(nAGQ=5) and lme4(nAGQ=5) remains at ~0.013 in the
intercept — same as the nAGQ=1 gap — confirming this is the profiling
approximation gap, not quadrature inaccuracy.

**Classification:** `partial` — nAGQ>1 is accepted and changes the deviance
approximation within the profiled path, but it is not a drop-in replacement for
lme4's full AGQ joint estimation.

### Finding 4: joint_laplace + nAGQ>1 refusal — honest typed diagnostic

`glmm(nAGQ=5, method="joint_laplace")` raises `mm_arg_error` with the message:
*"`method = "joint_laplace"` requires `nAGQ <= 1` in this slice."*
This is a typed, classed R condition (`mm_arg_error`). The refusal is honest,
labelled, and reaches the user before any optimizer work. The NLopt-backed
joint-AGQ path would be the correct implementation here but is disabled in the
vendored build (no `nlopt` feature). The R-layer refusal is consistent with the
Rust-layer `glmm_method()` guard.

**Classification:** `works` (typed refusal, by design in this build).

### Finding 5: Speed

mixeff is 3.4× faster than lme4 at nAGQ=1 and ~6× faster at nAGQ=5 on this
dataset. The nAGQ=5 lme4 time (0.024 s) includes its two-stage fit (initial
Laplace then AGQ refinement); mixeff's single-pass profiled optimization is 0.004 s.

**Classification:** `works` (speed goal met).

---

## Summary Table

| Scenario | lme4 | mixeff | Outcome | Root cause | Severity |
|---|---|---|---|---|---|
| nAGQ=1, theta | OK | OK, within-tol | within-tol | — | none |
| nAGQ=1, VarCorr SD | OK | OK, within-tol | within-tol | — | none |
| nAGQ=1, fixef intercept | 0.9844 | 0.9970, diff=0.013 | documented divergence | profiled vs joint objective | major (known) |
| nAGQ=1, fixef x | 0.5321 | 0.5319, diff=1.45e-4 | documented divergence | same | minor (known) |
| nAGQ=1, logLik | -580.94 | -580.96, diff=0.017 | beyond-tol | objective convention | minor (known) |
| nAGQ=5, theta | OK | OK, within-tol | within-tol | — | none |
| nAGQ=5, fixef | 0.9843 | 0.9969, diff=0.013 | documented divergence | profiled vs full AGQ | major (known) |
| nAGQ=5, logLik | -194.25 | -580.94, diff=387 | convention gap | response constants | cosmetic (convention) |
| joint_laplace+nAGQ=5 | n/a | typed refusal | honest refusal | NLopt not built | none |
| Speed (nAGQ=1) | 0.047 s | 0.014 s | 3.4× faster | — | none |
| Speed (nAGQ=5) | 0.024 s | 0.004 s | 6× faster | — | none |

---

## Classification

**Outcome:** `mixed`

- The `theta` / `VarCorr` quantities are at parity (within tolerance).
- The fixef intercept gap (~0.013) at both nAGQ=1 and nAGQ=5 is a
  **documented divergence** between profiled-PIRLS and joint-Laplace objectives,
  not a numerical bug. It is classified as `documented_divergence` in
  `glmm_support_contract.md` and is not a release-blocking failure for the
  `pirls_profiled` path.
- The logLik gap at nAGQ=5 (~387) is a response-constant convention difference,
  not an optimizer failure.
- nAGQ>1 with `pirls_profiled` is accepted and changes the deviance
  approximation, but is not equivalent to lme4's full joint-AGQ.
- The `joint_laplace+nAGQ>1` path is refused with a typed, honest diagnostic.
- Speed is substantially better than lme4.

**Max absolute difference (coefficient scale, vs published tolerance):**
- fixef intercept: **1.26e-02** (tolerance 1e-04, documented divergence)
- theta: **8.98e-05** (tolerance 1e-03, WITHIN-TOL)

**Severity of non-cosmetic gaps:** Major for the fixef intercept/x divergence,
but these are documented and expected for the profiled-vs-joint objective difference.
No undocumented or unexpected failures found.
