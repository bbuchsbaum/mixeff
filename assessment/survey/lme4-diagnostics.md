# lme4 / lmerTest Diagnostics & Convergence — Surface Survey

**Scope:** All user-facing functions, arguments, slots, and behaviours in the
`lme4` (2.0.1) and `lmerTest` (3.2.1) public surface that relate to model
diagnostics, convergence assessment, and influence/leverage analysis.  This
document is the reference against which `mixeff` parity is judged; it does
**not** assess mixeff.

---

## 1. Singularity Detection

### `isSingular(object, tol = getSingTol())`

Tests whether a fitted `merMod` is (near-)singular — i.e., any
random-effects variance-covariance matrix is estimated at or near the
boundary of the feasible parameter space (variance → 0, or |correlation| → 1).
Returns a single logical.  Users call this after fitting to decide whether the
random-effects structure is over-parameterised.

- `object`: a fitted `merMod` or `Covariance` object.
- `tol`: threshold below which a scaled eigenvalue is declared zero
  (default `1e-4`, changeable via `options(lme4.singular.tolerance = ...)`).

### `getSingTol()`

Helper that returns the active singularity tolerance — either the package
default (`1e-4`) or the value set via `options(lme4.singular.tolerance)`.
Used inside `isSingular` and as the default for `check.conv.singular`.

---

## 2. PCA of Random-Effects Covariance

### `rePCA(x)`

Performs a PCA of the random-effects variance-covariance estimates of a
fitted `merMod`.  Returns a `prcomplist` — a list of `prcomp` results, one
per random-effects grouping factor.  Each component shows the standard
deviations of orthogonal variance dimensions and the rotation (eigenvector)
matrix.  Users use this to identify which dimensions are effectively zero
(the "near-singular" pattern) even when individual variances and correlations
look plausible, which `isSingular` alone may miss for matrices of dimension
≥ 3.

---

## 3. Extended Convergence Checking

### `checkConv(derivs, coefs, ctrl, lbound, ubound, debug, nobs, ndim)`

Primary internal function (exported but rarely called directly) that
implements the post-optimisation convergence checks.  Inspects the gradient
and Hessian of the deviance at the reported optimum and emits warnings or
messages if convergence criteria are violated.

Key arguments:
- `derivs`: list with `$gradient` (numeric vector) and `$Hessian` (matrix)
  from finite-difference approximation at the optimum.
- `coefs`: current coefficient estimates (for scaling).
- `ctrl`: a list of lists produced by `lmerControl()$checkConv`, each
  specifying `action` and `tol` for one check.
- `lbound`, `ubound`: bounds for random-effects parameters (length used to
  count RE params).
- `debug`: logical; print triggering details even for "ignore"-level checks.
- `nobs`, `ndim`: thresholds for skipping derivative-based checks on large
  problems.

Returns a list with:
- `$code`: integer convergence code.
- `$messages`: character vector of triggered warning/message strings.

### `convergence` (help topic)

A narrative help page documenting the recommended workflow for diagnosing and
resolving convergence warnings:

1. Double-check model specification and data.
2. Tighten optimizer tolerances (`optCtrl = list(xtol_abs=1e-8, ftol_abs=1e-8)`).
3. Center and scale predictors.
4. Recompute gradient/Hessian via Richardson extrapolation (`numDeriv::hessian`,
   `numDeriv::grad`) and compare with `fm@optinfo$derivs`.
5. Restart from reported optimum or a slightly perturbed point.
6. Use `allFit()` as the "gold standard" multi-optimizer check.

---

## 4. Control Parameters Governing Convergence Checks

### `lmerControl(...)` / `glmerControl(...)`

The convergence-relevant arguments (all other args omitted here):

| Argument | Default | Purpose |
|---|---|---|
| `optimizer` | `"nloptwrap"` (lmer), `c("bobyqa","Nelder_Mead")` (glmer) | Optimizer(s) to use |
| `restart_edge` | `TRUE` (lmer), `FALSE` (glmer) | Restart when solution is at boundary |
| `boundary.tol` | `1e-5` | Distance from boundary that triggers boundary check |
| `calc.derivs` | `NULL` (auto) | Compute gradient/Hessian post-fit? Auto-skips for large problems |
| `use.last.params` | `FALSE` | Return last params evaluated rather than minimum-deviance params |
| `check.conv.nobsmax` | `1e4` | Skip derivative checks when `nobs` exceeds this |
| `check.conv.nparmax` | `10` (lmer), `20` (glmer) | Skip derivative checks when RE dim exceeds this |
| `check.conv.grad` | `.makeCC("warning", tol=2e-3)` | Rule for scaled-gradient check |
| `check.conv.singular` | `.makeCC("message", tol=getSingTol())` | Rule for singularity check |
| `check.conv.hess` | `.makeCC("warning", tol=1e-6)` | Rule for Hessian positive-definiteness check |
| `check.nobs.vs.rankZ` | `"ignore"` | Check nobs > rank(Z) |
| `check.nobs.vs.nlev` | `"stop"` | Check nobs < nlevels (identifiability) |
| `check.nlev.gtreq.5` | `"ignore"` | Warn if any RE grouping has < 5 levels |
| `check.nlev.gtr.1` | `"stop"` | Error if any RE grouping has ≤ 1 level |
| `check.nobs.vs.nRE` | `"stop"` | Check nobs vs RE levels |
| `check.rankX` | `"message+drop.cols"` | Check rank of fixed-effects design matrix |
| `check.scaleX` | `"warning"` | Check for badly scaled fixed-effects columns |
| `check.formula.LHS` | `"stop"` | Check formula has LHS |
| `optCtrl` | `list()` | Optimizer-specific control (maxfun, xtol_abs, ftol_abs, etc.) |

#### Optimizer-specific `optCtrl` keys

**nloptwrap** (lmer default):
- `ftol_abs` (default 1e-6) — stop on small absolute deviance change
- `ftol_rel` (default 0) — stop on small relative deviance change
- `xtol_abs` (default 1e-6) — stop on small parameter change
- `xtol_rel` (default 0) — stop on small relative parameter change
- `maxeval` (default 1000) — maximum function evaluations

**bobyqa** (glmer first stage):
- `rhobeg` (default 2e-3) — initial trust-region radius
- `rhoend` (default 2e-7) — final trust-region radius (~= xtol_abs)
- `maxfun` (default 10000)

**Nelder_Mead** (glmer second stage):
- `FtolAbs` (default 1e-5)
- `FtolRel` (default 1e-15)
- `XtolRel` (default 1e-7)
- `maxfun` (default 10000)

### `.makeCC(action, tol, relTol, ...)`

Utility that constructs a convergence-check control list.  `action` is one of
`"ignore"`, `"message"`, `"warning"`, `"stop"`.  `tol` is the absolute
tolerance; `relTol` is an optional relative tolerance.  Users call this when
customising check thresholds in `lmerControl()`.

---

## 5. `optinfo` Slot — Internal Convergence State

After fitting, `fm@optinfo` is a list with:

| Field | Contents |
|---|---|
| `$optimizer` | character: name of optimizer used |
| `$control` | list: optimizer control arguments that were passed |
| `$derivs$gradient` | numeric vector: finite-difference gradient at optimum |
| `$derivs$Hessian` | numeric matrix: finite-difference Hessian at optimum |
| `$conv$opt` | integer: optimizer return code (0 = success) |
| `$conv$lme4` | list: lme4-level convergence check results (messages) |
| `$feval` | integer: number of function evaluations |
| `$message` | character: optimizer-reported convergence message |
| `$warnings` | list: any optimizer warnings |
| `$val` | numeric vector: optimal theta parameter values |

Users inspect `fm@optinfo$derivs` directly (or via `getME(fm,"devfun")`) to
recheck convergence with external tools like `numDeriv`.  The reciprocal
condition number of the Hessian (`Matrix::rcond(fm@optinfo$derivs$Hessian)`)
is a standard diagnostic for near-singular optimisation.

---

## 6. Multi-Optimizer Convergence Check

### `allFit(object, meth.tab, data, verbose, show.meth.tab, maxfun, parallel, ncpus, cl, catch.errs, start_from_mle)`

Refits a model with all available optimizers and compares results.  Considered
the "gold standard" for determining whether a convergence warning is a false
positive.  Returns an `allFit` object (a named list of `merMod` fits).

Key arguments:
- `object`: a fitted `merMod`.
- `meth.tab`: optional data frame with `method` and `optimizer` columns to
  restrict which optimizers are tried.
- `maxfun`: maximum function evaluations (auto-translated per optimizer).
- `verbose`: logical; print progress.
- `show.meth.tab`: if `TRUE`, return the methods table instead of fits.
- `catch.errs`: wrap fits in `tryCatch` (default `TRUE`).
- `start_from_mle`: initialise from the MLE of the original fit.
- `parallel`, `ncpus`, `cl`: parallel execution support.

**Built-in optimizers tried:** `bobyqa`, `Nelder_Mead`, `nlminbwrap`,
`nloptwrap.NLOPT_LN_NELDERMEAD`, `nloptwrap.NLOPT_LN_BOBYQA`.
Optional (require packages): `optimx`-wrapped L-BFGS-B and nlminb,
`dfoptim::nmkb` (appears as `nmkbw`).

`summary(allFit_object)` returns a list with:
- `$which.OK` — logical vector: which optimizers converged without error.
- `$msgs` — messages/warnings per optimizer.
- `$fixef` — matrix of fixed-effect estimates across optimizers.
- `$llik` — vector of log-likelihoods.
- `$sdcor` — matrix of random-effect SDs and correlations.
- `$theta` — matrix of Cholesky-scale RE parameters.
- `$times` — timing per optimizer.
- `$feval` — function evaluation counts.

---

## 7. `getME()` — Low-Level Parameter Access

### `getME(object, name, ...)`

Extracts internal model components from a `merMod`.  Diagnostically relevant
names:

| Name | Contents |
|---|---|
| `"theta"` | Cholesky factor entries of relative RE covariance matrices |
| `"ST"` | S and T factors in the TSST' Cholesky factorisation |
| `"Lambda"` | Relative covariance factor Λ |
| `"Lambdat"` | Transpose of Λ |
| `"L"` | Sparse Cholesky factor of the penalised RE model |
| `"lower"` | Lower bounds on theta (0 for diagonal elements, -Inf for off-diagonal) |
| `"devfun"` | The REML or ML deviance function closure (for external gradient/Hessian checks) |
| `"devarg"` | Named list of current argument values to devfun |
| `"par"` | RE covariance parameter estimates (same as theta for simple models) |

Users extract `theta` and `lower` to identify which parameters are on the
boundary (`theta == lower` indicates a boundary constraint is active).

---

## 8. Deviance Function Access

### `devfun2(fm, useSc, scale)`

Returns the deviance function reparameterised on the SD/correlation scale
(rather than the Cholesky / theta scale).  Useful for manual profiling and
for checking the Hessian in a more interpretable parameterisation.

- `useSc`: whether a scale parameter (sigma) is included.
- `scale`: `"sdcor"` (default) or `"varcov"`.

Returns a function whose attributes carry:
- `$optimum` — named parameter vector at the MLE.
- `$basedev` — deviance at the MLE.
- `$thopt` — optimal theta.
- `$stderr` — SEs of fixed-effect parameters.

---

## 9. Profile Likelihood

### `profile(fitted, which, alphamax, maxpts, delta, delta.cutoff, verbose, devtol, devmatchtol, maxmult, startmethod, optimizer, control, signames, parallel, ncpus, cl, prof.scale, ...)`

Computes the profile likelihood for each parameter (or a subset specified by
`which`).  Returns a `thpr` object (inherits `data.frame`) with columns
`.zeta` (signed square-root likelihood ratio), the profiled parameter value,
and all other parameters at their conditional optima.

Key arguments:
- `which`: `NULL` (all), integer indices, `"theta_"`, `"beta_"`, or
  specific parameter names such as `".sigma"`.
- `alphamax`: maximum alpha for the profiled range (default 0.01, giving
  99% coverage).
- `maxpts`: max profile points per parameter in each direction.
- `prof.scale`: `"sdcor"` (SD/correlation) or `"varcov"` (variance/covariance).

### `varianceProf(x, ranef)` / `logProf(x, base, ranef, sigIni)`

Transform a profile object: `varianceProf` converts from SD to variance
scale; `logProf` converts to log-SD scale.  Both return modified `thpr`
objects for plotting or further analysis.

### `confint(profile_object, ...)` (method for `thpr`)

Extracts profile likelihood confidence intervals from a `thpr` object.
Also: `confint.merMod(..., method="profile")` which calls `profile()` then
extracts CIs directly.

---

## 10. Influence Diagnostics

### `influence(model, groups, data, maxfun, do.coef, start, parallel, ncpus, cl, ...)`
(S3 method for `merMod`)

Computes deletion influence by refitting the model with each group omitted in
turn.  Returns an `influence.merMod` object.

Key arguments:
- `groups`: character vector of grouping factor name(s) to delete by group;
  if omitted, each individual observation row is deleted.
- `maxfun`: max function evaluations after deletion (default large; set
  `maxfun=20` for fast approximation for LMMs, `100` for GLMMs).
- `do.coef`: if `FALSE`, skip fixed-effect coefficient collection
  (returns only hat values; faster).
- `start`: starting values for refits.
- `parallel`, `ncpus`, `cl`: parallel execution.

The `influence.merMod` object contains:
- `"fixed.effects"` — full-data fixed-effect estimates.
- `"fixed.effects[-groups]"` — matrix of leave-one-group-out fixed effects.
- `"var.cov.comps"` — full-data variance-covariance parameters.
- `"var.cov.comps[-groups]"` — matrix of leave-one-out VC parameters.
- `"vcov"` — full-data fixed-effects covariance matrix.
- `"vcov[-groups]"` — list of leave-one-out fixed-effects covariance matrices.
- `"groups"` — grouping factor name(s).
- `"deleted"` — the composite deletion factor.
- `"converged"` — logical vector: did each deletion refit converge?
- `"function.evals"` — vector of function evaluation counts per deletion.

### `cooks.distance(model, ...)`
(S3 method for `influence.merMod`)

Computes Cook's distance for each deleted group from an `influence.merMod`
object.  Measures overall influence on the fixed-effect estimates.

### `dfbeta(model, which, ...)`
(S3 method for `influence.merMod`)

Raw change in coefficient estimates when each group is deleted.

- `which`: `"fixed"` (default) — influence on fixed effects; `"var.cov"` —
  influence on variance-covariance components.

### `dfbetas(model, ...)`
(S3 method for `influence.merMod`)

Standardised `dfbeta` — divided by the standard error of each coefficient.

---

## 11. Hat Values / Leverage

### `hatvalues(model, fullHatMatrix, ...)`
(S3 method for `merMod`)

Returns the diagonal elements of the hat matrix H = Z(Z'Z)^{-1}Z' + X(X'X)^{-1}X'
(schematically; the actual computation uses the internal Cholesky factors).
High leverage points have hat values near 1.

- `fullHatMatrix`: if `TRUE`, return the full n × n hat matrix rather than
  just the diagonal.

**Note:** Only interpretable for LMMs; concept is not well-defined for GLMMs.

---

## 12. Residual Diagnostics

### `residuals(object, type, scaled, ...)`
(S3 method for `merMod`)

Extracts model residuals.

- `type`: `"response"` (default for LMMs) or `"deviance"` (default for
  GLMMs) — also `"pearson"`, `"working"`.
- `scaled`: if `TRUE`, scale by the residual standard deviation.

### `plot(x, form, abline, id, idLabels, grid, ...)`
(S3 method for `merMod`)

Residual diagnostic plot (Pearson residuals vs. fitted values by default).
- `form`: a formula specifying what to plot (default
  `resid(., type="pearson") ~ fitted(.)`).
- `id`: proportion of extreme residuals to label (passed to `lattice`).
- `idLabels`: custom labels for extreme points.

### `fortify.merMod(model, data, ...)`

Adds fitted values, residuals (response, Pearson, deviance, working, partial),
hat values, and Cook's distances to the model's data frame.  Used by
`ggplot2::ggplot()`.  Deprecated in favour of `broom.mixed`.

### `dotplot.ranef.mer(x, data, main, ...)`

Dotplots of random effects with 95% prediction intervals (conditional modes
± 1.96 × conditional standard deviations).  Useful for identifying groups
with extreme random effects.

---

## 13. Random Effects — Conditional Variances

### `ranef(object, condVar, ...)`

Extracts conditional modes (BLUPs) of random effects.  When `condVar=TRUE`
(the default), attaches the conditional variance-covariance matrices as a
`"postVar"` attribute on each random-effects data frame.  These conditional
variances are used to construct:

- Prediction intervals on random effects (dotplots).
- The `arm::se.ranef()` function and similar tools.

---

## 14. Rank / Design Diagnostics

### `check.rankX` (lmerControl argument)

Computes `Matrix::rankMatrix(X)` and compares it to `ncol(X)`.  If X is
rank-deficient, drops columns (various options: `"message+drop.cols"`,
`"silent.drop.cols"`, `"warn+drop.cols"`, `"stop.deficient"`, `"ignore"`).
Users see this as a message or warning at fit time.

### `check.scaleX` (lmerControl argument)

Detects badly scaled fixed-effects columns (e.g., predictors on very different
scales).  Options: `"warning"`, `"stop"`, `"silent.rescale"`,
`"message+rescale"`, `"warn+rescale"`, `"ignore"`.

### `check.nobs.vs.rankZ` (lmerControl argument)

Checks `nobs > rank(Z)` (identifiability of random effects variance).

### `check.nobs.vs.nlev` / `check.nobs.vs.nRE` (lmerControl arguments)

Check that there are more observations than random-effects levels — necessary
for variance identifiability.

---

## 15. Optimizer Utilities

### `Nelder_Mead(fn, par, lower, upper, control, ...)`

Built-in Nelder-Mead optimizer.  Exported so users can call it directly or
specify it as an optimizer in `lmerControl`.

### `nloptwrap(par, fn, lower, upper, control, ...)`

Wrapper around `nloptr` optimizers; default for `lmer`.  `control` accepts
`algorithm` (default `"NLOPT_LN_BOBYQA"`), `xtol_rel`, `ftol_rel`,
`xtol_abs`, `ftol_abs`, `maxeval`.

### `nlminbwrap(par, fn, lower, upper, control, ...)`

Wrapper around `stats::nlminb`; accepts standard `nlminb` control arguments.

### `allFit(show.meth.tab=TRUE)`

When called with `show.meth.tab=TRUE`, returns a data frame of all
optimizer/method combinations available in the current session (including
optional `optimx` and `dfoptim` entries).

---

## 16. isLMM / isGLMM / isNLMM / isREML — Model Type Predicates

These predicates are frequently used in diagnostic workflows to condition on
model type before applying type-specific checks:

- `isLMM(object)` — TRUE for `lmer` fits.
- `isGLMM(object)` — TRUE for `glmer` fits.
- `isNLMM(object)` — TRUE for `nlmer` fits.
- `isREML(object)` — TRUE if fitted by REML (relevant: gradient/Hessian
  checks are on the REML criterion for LMMs).

---

## 17. `VarCorr` — Variance-Covariance of Random Effects

### `VarCorr(x, sigma, ...)`

Extracts the random-effects variance-covariance matrices and residual
variance.  Returns a `VarCorr.merMod` object whose `print` method shows
standard deviations and correlations.  Diagnostically used to check:
- Whether any variance is estimated as exactly zero (singularity).
- Whether any correlation is exactly ±1.

---

## 18. lmerTest Diagnostics

`lmerTest` wraps `lmer` to add Satterthwaite or Kenward-Roger denominator
degrees of freedom.  Its diagnostic additions are minimal — the main
convergence/singularity infrastructure is inherited from `lme4`.

### `lmerTest::lmer(...)` convergence behaviour

The `lmerTest::lmer` wrapper passes all `lme4::lmerControl` arguments
through unchanged; convergence checks, `isSingular`, and `allFit` behaviour
are identical to `lme4`.

### `anova(model, ddf, ...)` (lmerTest method)

Satterthwaite/KR ddf computation will warn (or produce `NaN` ddf) if the
model is singular or the Hessian is ill-conditioned, which is a diagnostic
signal for the user.

---

## 19. Summary of User Workflow

The typical lme4 diagnostics workflow is:

```
1. Fit:          fm <- lmer(y ~ x + (x|g), data)
2. Singularity:  isSingular(fm)          # quick check
3. RE structure: rePCA(fm)               # detailed decomposition
4. Convergence:  fm@optinfo$conv$lme4    # check for messages/warnings
                 fm@optinfo$derivs       # gradient, Hessian
5. Recheck:      numDeriv::hessian(devfun, getME(fm,"theta"))
6. Multi-opt:    allFit(fm)              # gold standard
7. Leverage:     hatvalues(fm)
8. Influence:    inf <- influence(fm, "g")
                 cooks.distance(inf)
                 dfbeta(inf)
                 dfbetas(inf)
9. Residuals:    plot(fm)
                 residuals(fm, type="pearson", scaled=TRUE)
10. Profile CI:  profile(fm, which="theta_")
                 confint(fm, method="profile")
```

---

## 20. Cross-Reference Index

| Function / Argument | Section |
|---|---|
| `isSingular` | §1 |
| `getSingTol` | §1 |
| `rePCA` | §2 |
| `checkConv` | §3 |
| `convergence` (help topic) | §3 |
| `lmerControl` / `glmerControl` | §4 |
| `.makeCC` | §4 |
| `optCtrl` keys (ftol_abs, xtol_abs, …) | §4 |
| `fm@optinfo` (gradient, Hessian, rcond) | §5 |
| `allFit` | §6 |
| `summary.allFit` | §6 |
| `getME("theta")`, `getME("devfun")` | §7 |
| `devfun2` | §8 |
| `profile.merMod` | §9 |
| `varianceProf`, `logProf` | §9 |
| `influence.merMod` | §10 |
| `cooks.distance.influence.merMod` | §10 |
| `dfbeta.influence.merMod` | §10 |
| `dfbetas.influence.merMod` | §10 |
| `hatvalues.merMod` | §11 |
| `residuals.merMod` | §12 |
| `plot.merMod` | §12 |
| `fortify.merMod` | §12 |
| `dotplot.ranef.mer` | §12 |
| `ranef(..., condVar=TRUE)` | §13 |
| `check.rankX`, `check.scaleX` | §14 |
| `check.nobs.vs.*` | §14 |
| `Nelder_Mead`, `nloptwrap`, `nlminbwrap` | §15 |
| `isLMM`, `isGLMM`, `isREML` | §16 |
| `VarCorr` | §17 |
| lmerTest diagnostic behaviour | §18 |
