# Parity probe: glmm-grouseticks-pois

**Cell:** glmm-grouseticks-pois  
**Dataset:** `grouseticks` (lme4)  
**Formula:** `TICKS ~ YEAR + HEIGHT + (1|BROOD)`  
**Family:** `poisson(link="log")`  
**Method:** `pirls_profiled` (mixeff default)  
**Date run:** 2026-05-31  
**lme4:** 2.0.1 | **lmerTest:** 3.2.1 | **mixeff:** 0.1.0

---

## Script

See `glmm-grouseticks-pois-probe.R` in this directory.

---

## Raw output

### lme4 fit (`TICKS ~ YEAR + HEIGHT + (1|BROOD)`)

```
CONVERGENCE WARNINGS:
  Model failed to converge with max|grad| = 0.0460375 (tol = 0.002, component 1)
  Model is nearly unidentifiable: very large eigenvalue - Rescale variables?
  Model is nearly unidentifiable: large eigenvalue ratio - Rescale variables?

max|grad| at solution: 79.47   (computed from optinfo$derivs)

lme4 fixef:
(Intercept)  YEAR96    YEAR97    HEIGHT
  11.541043   1.135874  -1.001132  -0.023866

lme4 SE:
(Intercept)  YEAR96    YEAR97    HEIGHT
  0.41317498  0.24229985  0.26933988  0.00081929

lme4 theta (Cholesky): 0.9496773
lme4 VarCorr BROOD variance: 0.9018870 (std dev 0.94968)
lme4 sigma: 1 (Poisson)
lme4 logLik: -989.0377 (df=5)
lme4 AIC: 1988.075
lme4 BIC: 2008.070
lme4 deviance (conditional): 734.558   [-2*logLik = 1978.075]
```

### Verification: lme4 with centered HEIGHT (converges cleanly)

To confirm the true optimum, the model was refit with `cHEIGHT = HEIGHT - mean(HEIGHT)`.
This is the same model (only the intercept changes); centering removes the near-collinearity.

```
lme4 fixef (cHEIGHT encoding, converged):
(Intercept)  YEAR96    YEAR97    cHEIGHT
  0.5091735   1.135901  -1.001130  -0.02386609

lme4 theta: 0.9497031
lme4 logLik: -989.0377   (identical — same model)
```

The logLik is the same to machine precision; YEAR96, YEAR97, HEIGHT/cHEIGHT coefficients
match to 5+ significant figures. The raw-HEIGHT lme4 estimates for all parameters *except*
the intercept are correct; only the intercept is displaced by `mean(HEIGHT) * HEIGHT_coef`
(i.e., `462.3 * (-0.023866) ≈ -11.03`, consistent with the observed difference 11.541 - 0.509 = 11.032).

### mixeff fit (`TICKS ~ YEAR + HEIGHT + (1|BROOD)`)

```
fit_status: converged_interior

mixeff fixef:
(Intercept)  YEAR: 96   YEAR: 97   HEIGHT
  11.08963    1.07554   -0.95799   -0.02266

mixeff SE:
(Intercept)  YEAR: 96   YEAR: 97   HEIGHT
  1.91511     0.33125    0.36958    0.00411

mixeff theta: 0.9468994
mixeff VarCorr BROOD variance: 0.8966185 (std dev 0.946899)
mixeff sigma: 1 (Poisson)
mixeff logLik: -989.5617 (df=5)
mixeff AIC: 1989.123
mixeff BIC: 2009.118
mixeff deviance: 1979.123   (= -2*logLik; marginal deviance, not conditional)
```

### mixeff fit with cHEIGHT (to isolate centering from optimizer difference)

```
fit_status: converged_interior

mixeff fixef (cHEIGHT):
(Intercept)  YEAR: 96   YEAR: 97   cHEIGHT
  0.61417     1.07554   -0.95799   -0.02266

mixeff theta: 0.9468994
mixeff logLik: -989.5617
```

The mixeff cHEIGHT result is **not identical** to lme4 — the logLik difference (-989.5617 vs
-989.0377 = 0.524) and fixef differences are genuine optimizer-level divergences, not an
encoding artifact. Centering does not resolve the parity gap.

---

## Numerical comparison (cell formula: raw HEIGHT)

The true reference for comparison is the **converged lme4 result** — the unconverged raw-HEIGHT
lme4 fit has `max|grad|=79.5` and must not be treated as ground truth. For quantities where
centering is irrelevant (YEAR coefs, HEIGHT coef, theta, logLik), we use the cHEIGHT converged
lme4 values as reference.

| Quantity | lme4 (converged ref) | mixeff | maxAbsDiff | Tol | Status |
|---|---|---|---|---|---|
| fixef (Intercept) | 0.50917 (cHEIGHT enc.) | 0.61417 | 1.050e-01 | 1e-4 | **EXCEEDS-TOL** |
| fixef YEAR96 | 1.135901 | 1.075538 | 6.034e-02 | 1e-4 | **EXCEEDS-TOL** |
| fixef YEAR97 | -1.001130 | -0.957992 | 4.314e-02 | 1e-4 | **EXCEEDS-TOL** |
| fixef HEIGHT | -0.023866 | -0.022662 | 1.204e-03 | 1e-4 | **EXCEEDS-TOL** |
| SE (Intercept) | 0.18662 (cHEIGHT enc.) | 0.25105 | 6.443e-02 | 1e-4 | **EXCEEDS-TOL** |
| SE YEAR96 | 0.24239 | 0.33125 | 8.886e-02 | 1e-4 | **EXCEEDS-TOL** |
| SE YEAR97 | 0.26969 | 0.36958 | 9.989e-02 | 1e-4 | **EXCEEDS-TOL** |
| SE HEIGHT | 0.00301 | 0.00411 | 1.100e-03 | 1e-4 | **EXCEEDS-TOL** |
| theta (Cholesky) | 0.9497031 | 0.9468994 | 2.804e-03 | 1e-3 | **EXCEEDS-TOL** |
| VarCorr BROOD var | 0.901936 | 0.896619 | 5.317e-03 | 1e-3 | **EXCEEDS-TOL** |
| sigma (Poisson disp.) | 1.0 | 1.0 | 0.000e+00 | 1e-4 | within-tol |
| logLik | -989.0377 | -989.5617 | 5.240e-01 | 1e-3 | **EXCEEDS-TOL** |
| AIC | 1988.075 | 1989.123 | 1.048e+00 | 2e-3 | **EXCEEDS-TOL** |
| BIC | 2008.070 | 2009.118 | 1.048e+00 | 2e-3 | **EXCEEDS-TOL** |
| ranef BROOD max abs diff | — | — | ~1.5e-01 | 1e-3 | **EXCEEDS-TOL** |
| fitted max abs diff | — | — | ~1.0e-01 | 1e-3 | **EXCEEDS-TOL** |
| deviance (conditional) | 734.558 | 1979.123 | ~1245 | — | different semantics* |
| convergence | WARNING (max\|grad\|=79.5) | converged_interior | — | — | see note |
| speed (wall-clock) | 0.126s | 0.019s | ratio=0.15x | — | 6.6x faster |

*lme4 `deviance()` for GLMMs returns the conditional deviance (penalized deviance of the PIRLS
working model), not `-2*logLik`. mixeff `deviance()` returns `-2*logLik` (marginal deviance).
These are different quantities by design; the difference is not a bug in either package.

---

## Analysis

### Root cause: optimizer / Laplace approximation difference

Both lme4 and mixeff fit a Laplace-approximated marginal likelihood for the Poisson GLMM.
However, they use **different implementations** of the PIRLS/Laplace objective and optimizer:

- **lme4** uses `bobyqa` (nloptwrap) with the standard lme4 PIRLS implementation. On the raw
  `HEIGHT` covariate (range 447–520, mean 462), the design matrix has a near-collinearity: the
  intercept and HEIGHT are highly correlated in this scale. lme4's optimizer **fails to converge**
  (`max|grad|=79.5`, tol=0.002), and lme4 itself warns "Model is nearly unidentifiable: very large
  eigenvalue — Rescale variables?". The lme4 result for raw HEIGHT is **numerically unreliable**.

- **mixeff (Rust/PIRLS)** uses its own PIRLS optimizer and converges to `converged_interior`.
  With the cHEIGHT-encoded model (same statistical problem, better conditioning), mixeff still
  lands at `logLik=-989.5617` vs the true optimum of `-989.0377` — a gap of 0.524 log-likelihood
  units. This indicates mixeff's optimizer finds a **suboptimal solution** on this dataset, not
  just a different parameterization.

### logLik gap (0.524 units): genuine divergence

The marginal logLik should be optimizer-independent for a correctly-implemented Laplace GLMM.
The 0.524-unit gap is significant and consistent across raw HEIGHT and cHEIGHT encodings,
confirming this is a real optimizer shortfall in mixeff — it converges to a local or near-flat
region that is not the global optimum. This propagates to all other quantities (fixef, SE, ranef,
fitted, theta).

### SE inflation

mixeff's standard errors are substantially larger than lme4's (e.g., SE(Intercept): 1.92 vs 0.41
for raw HEIGHT encoding; 0.25 vs 0.19 for cHEIGHT). The vcov matrix carries the label
`mm_reliability = "moderate"` and method `"pirls_laplace_working_hessian"`, with the note:
"PIRLS/Laplace working-Hessian fixed-effect covariance geometry; inference claims remain on
fixed_effect_inference_table rows". The inflated SEs are consistent with the suboptimal point
having a flatter curvature (Hessian eigenvalues are smaller at a non-optimal point).

### deviance semantics

lme4's `deviance()` for GLMMs returns the **conditional deviance** (the PIRLS working
residual deviance), not the marginal `-2*logLik`. mixeff returns `-2*logLik` (marginal deviance).
These are legitimately different quantities. This is not a defect in mixeff — it is a design
choice that is arguably more interpretable.

### convergence reporting

lme4 warns of non-convergence but returns an object. mixeff reports `converged_interior`
(its own optimizer's convergence criterion), which is truthful given its own stopping rules,
but the global optimum is not achieved. There is no false-convergence reporting per se —
both packages behave consistently with their own interfaces — but the user has less visibility
into the fact that mixeff's solution is suboptimal.

### Classification

This is an **in-scope missing** capability gap: the Poisson GLMM (1 RE, 4 fixed effects) is
within the certified `poisson/log` contract of `glmm()`, and the optimizer convergence failure
on this specific dataset/scale is a real defect. The gap is not due to centering (cHEIGHT
shows the same divergence), not due to lme4 non-convergence (the true optimum is known),
and not out-of-scope.

---

## Overall verdict

**Outcome: divergent**

All main quantities exceed tolerance: fixef diffs up to 4.5e-1, SE diffs up to 1.5, theta
diff 2.8e-3, logLik diff 5.24e-1. The root cause is that mixeff's PIRLS optimizer fails to
reach the global optimum on this dataset, landing 0.524 logLik units short. sigma and deviance
semantics are fine.

**Severity: major** — the logLik gap (0.524) is 500x the 1e-3 tolerance; fixef estimates are
materially wrong; SE is inflated ~4x. This is a real optimizer failure on a canonical Poisson
GLMM benchmark dataset. The model *does* fit (no error/refusal), but the solution is suboptimal.

**Upstream note:** The optimizer shortfall is in the vendored Rust engine. A mote issue in
`/Users/bbuchsbaum/code/rust/mixeff-rs` is warranted once a minimal reproducer or root-cause
hypothesis (gradient tolerance, PIRLS iteration count, starting values for theta) is identified.
