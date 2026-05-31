# lme4 / lmerTest — LMM Fitting & Control: Public Surface Reference

Survey date: 2026-05-31  
lme4 version: 2.0.1 | lmerTest version: 3.2.1  
Purpose: exhaustive reference of every user-facing function, argument, and
behavior in the "LMM fitting & control" capability family.  This document is
the reference against which mixeff is assessed; it does **not** evaluate
mixeff.

---

## 1. Primary fit function

### `lmer()` (lme4 and lmerTest)

```
lmer(formula, data = NULL, REML = TRUE, control = lmerControl(),
     start = NULL, verbose = 0L, subset, weights, na.action,
     offset, contrasts = NULL, devFunOnly = FALSE)
```

lmerTest ships its own `lmer()` that wraps `lme4::lmer` and adds slots for
Satterthwaite denominator-df computation.  The argument signature is
identical; the return class changes from `lmerMod` to `lmerModLmerTest`.

| Argument | Type | What it does and why users rely on it |
|---|---|---|
| `formula` | `formula` | Two-sided lme4-style formula. Left of `~` is the response; right combines fixed-effect terms and random-effect terms separated by `\|`. Users rely on it for full formula-based model specification including complex random structures (multiple grouping factors, slopes, nested terms). Covariance-structure tags (`us`, `diag`, `cs`, `ar1`) new in lme4 2.0. |
| `data` | `data.frame` or `NULL` | Data frame containing all variables. Strongly recommended; omitting it breaks `update()`, `drop1()`, and other post-fit methods. |
| `REML` | `logical` | `TRUE` (default) fits by restricted maximum likelihood; `FALSE` fits by maximum likelihood. ML is required for likelihood-ratio tests comparing fixed-effect structures; REML gives less biased variance estimates. Users switch this deliberately for model comparison via `anova()`. |
| `control` | `merControl` list | Control object from `lmerControl()`. Governs optimizer choice, convergence tolerances, pre-fit data checks, and post-fit convergence checks. See §2. |
| `start` | numeric vector or named list | Starting values for the covariance (theta) parameters. Numeric input treated as `list(par = start)`. Values are on the Cholesky-factor scale of the relative covariance matrix. Named lists may have component `par` or `theta`. Users supply this to avoid convergence problems or to reproduce results from an earlier fit. |
| `verbose` | `integer` | `0` (silent), `1` (print optimizer trace), `2` (also print PIRLS steps). Indispensable for diagnosing convergence. |
| `subset` | expression | Expression selecting a subset of rows. Can be logical, numeric index, or character rownames. Avoids pre-filtering the data frame; interacts correctly with `na.action`. |
| `weights` | numeric vector | Prior case weights (not normalized). Diagonal of residual covariance is `sigma^2 / weights_i`. Used for heterogeneous measurement precision (e.g., known measurement errors, aggregated observations). |
| `na.action` | function | Governs missing-data handling. Default `na.omit` drops rows with any NA. `na.exclude` preserves indices so `fitted()` / `residuals()` return vectors of original length with `NA` fills. `na.pass` / `na.fail` are also valid. |
| `offset` | numeric vector | Fixed a-priori offset added to the linear predictor before fitting. Length must equal `nrow(data)`. Can also appear in the formula as `offset(var)`. Multiple offset terms are summed. Used for exposure adjustment, known baselines, or offsets derived from design. |
| `contrasts` | named list or `NULL` | Passed to `model.matrix.default` as `contrasts.arg`. Overrides the default contrast coding (treatment coding for unordered factors, polynomial for ordered). Users rely on this to specify Helmert, sum-to-zero, or custom contrast matrices without altering global options. |
| `devFunOnly` | `logical` | If `TRUE`, returns the deviance function instead of a fitted model. The function takes a single numeric vector (the theta/par vector) and returns the (restricted) deviance. Used for manual optimization, profiling, gradient checks with `numDeriv`, and convergence diagnostics. |

---

## 2. `lmerControl()` — fit control object

```
lmerControl(optimizer = "nloptwrap",
            restart_edge = TRUE, boundary.tol = 1e-5,
            calc.derivs = NULL, use.last.params = FALSE,
            sparseX = FALSE, standardize.X = FALSE, autoscale = NULL,
            check.nobs.vs.rankZ = "ignore",
            check.nobs.vs.nlev  = "stop",
            check.nlev.gtreq.5  = "ignore",
            check.nlev.gtr.1    = "stop",
            check.nobs.vs.nRE   = "stop",
            check.rankX = c("message+drop.cols", ...),
            check.scaleX = c("warning", ...),
            check.formula.LHS = "stop",
            check.conv.nobsmax = 1e4,
            check.conv.nparmax = 10,
            check.conv.grad     = .makeCC("warning", tol=2e-3),
            check.conv.singular = .makeCC("message", tol=getSingTol()),
            check.conv.hess     = .makeCC("warning", tol=1e-6),
            optCtrl = list(),
            mod.type = "lmer")
```

Returns a list of class `merControl` consumed by `lmer()`.  Some `check.*`
(pre-fit checking) parameters can also be set globally via `options()`.

### 2a. General / optimizer parameters

| Parameter | Default | What it does |
|---|---|---|
| `optimizer` | `"nloptwrap"` | Name of the nonlinear optimizer. Built-ins: `"nloptwrap"` (NLopt BOBYQA via `nloptr`), `"bobyqa"` (minqa BOBYQA), `"Nelder_Mead"` (lme4's internal), `"nlminbwrap"` (base `nlminb`). Also accepts any function satisfying the box-constrained optimizer interface (takes `fn`, `par`, `lower`, `upper`, `control`; returns `par`, `fval`, `conv`, `message`). `optimx`-wrapped optimizers (L-BFGS-B, nlminb, etc.) are supported by passing the method name to `optCtrl`. `NULL` sets up all model structures but skips optimization (all parameters returned as `NA`). |
| `restart_edge` | `TRUE` | Whether to restart optimization when a solution lands on the boundary (zero variance or ±1 correlation). Only implemented for LMM. Helps escape degenerate saddle-point solutions at boundaries. |
| `boundary.tol` | `1e-5` | Distance from a parameter boundary within which a boundary check is triggered. Set to 0 to disable. |
| `calc.derivs` | `NULL` | Whether to compute gradient and Hessian at the optimum for convergence checking. `NULL` means compute only when both `nobs < check.conv.nobsmax` and `ndim < check.conv.nparmax`. Gradient/Hessian are estimated by finite differences. |
| `use.last.params` | `FALSE` | Return the last evaluated parameter values rather than the values at the minimum deviance. A backward-compatibility escape hatch; not recommended for new code. |
| `sparseX` | `FALSE` | Use a sparse fixed-effects design matrix. Currently inactive in lme4 (reserved). |
| `standardize.X` | `FALSE` | Scale columns of the X matrix before fitting. Not yet implemented. |
| `autoscale` | `NULL` | Automatically scale continuous covariates and back-transform estimates. LMM-only; incompatible with `glmer`. Prevents numeric stability issues without user intervention. |
| `optCtrl` | `list()` | Named list of optimizer-specific control parameters passed through to the optimizer. For `nloptwrap`: `algorithm` (default `"NLOPT_LN_BOBYQA"`), `ftol_abs` (1e-6), `xtol_abs` (1e-6), `ftol_rel` (0), `xtol_rel` (0), `maxeval` (1000). For `bobyqa`/`Nelder_Mead`: `maxfun`. For `nlminbwrap`: `nlminb` defaults. For `optimx`-based: `method`, `kkt`, `maxit`. |

### 2b. Pre-fit data and formula checking options

Each `check.*` parameter accepts one of `"stop"`, `"warning"`, `"message"`,
`"ignore"` (plus variant forms for some).  The pre-fit checks can be
over-ridden globally via `options(lmerControl = list(check.* = ...))`.

| Parameter | Default | What it checks |
|---|---|---|
| `check.nobs.vs.rankZ` | `"ignore"` | Whether `nobs > rank(Z)` (random-effects design matrix). When violated, random-effect variances are typically not identifiable. |
| `check.nobs.vs.nlev` | `"stop"` | Whether `nobs < nlevels(grouping_factor)`. Fires when there are more grouping levels than observations, making variance estimation impossible. |
| `check.nlev.gtreq.5` | `"ignore"` | Whether all random effects have ≥ 5 levels. Fewer levels make variance estimation unreliable. |
| `check.nlev.gtr.1` | `"stop"` | Whether all random effects have > 1 level. Exactly 1 level means the random effect is unidentifiable. |
| `check.nobs.vs.nRE` | `"stop"` | Whether `nobs > number_of_RE_parameters`. Related to identifiability of random-effect covariances. |
| `check.rankX` | `"message+drop.cols"` | Rank of the fixed-effects design matrix X. If `X` is rank-deficient, columns are dropped with a message. Options control whether to drop silently, warn, or error. |
| `check.scaleX` | `"warning"` | Whether columns of X are on very different scales. Poor scaling causes numerical instability in the optimizer. |
| `check.formula.LHS` | `"stop"` | Whether the formula has a left-hand side. Primarily for internal use with `simulate.merMod`. |

### 2c. Post-fit convergence checking options

| Parameter | Default | What it checks |
|---|---|---|
| `check.conv.nobsmax` | `1e4` | Skip derivative-based convergence checks when `nobs` exceeds this threshold (finite-difference Hessians are unreliable and slow for large N). Set to `Inf` to always check. |
| `check.conv.nparmax` | `10` | Skip derivative-based checks when the number of variance/covariance parameters exceeds this. Set to `Inf` to always check. |
| `check.conv.grad` | `warning, tol=2e-3` | Check scaled gradient of deviance at the optimum. Scaled relative to inverse Cholesky of the Hessian (Wald SD scale). Action + tolerance configured via `.makeCC()`. |
| `check.conv.singular` | `message, tol=getSingTol()` | Check for singular fits (parameters on the feasible space boundary). Default tolerance from `getSingTol()`. |
| `check.conv.hess` | `warning, tol=1e-6` | Check reciprocal condition number of the Hessian at the optimum. Near-singular Hessian indicates flat likelihood ridges. |

### 2d. Helper: `.makeCC(action, tol, relTol, ...)`

Constructs a convergence-check specification list with fields `action`,
`tol`, and optionally `relTol`.  Used as the default value of
`check.conv.grad`, `check.conv.singular`, `check.conv.hess`.  Users call it
when they want to customize the tolerance and action together.

---

## 3. Optimizer wrappers

### `nloptwrap(par, fn, lower, upper, control = list(), ...)`

Default LMM optimizer wrapper.  Routes to NLopt's BOBYQA
(`"NLOPT_LN_BOBYQA"`) by default.  Accepts any algorithm string from
`nloptr::nloptr.print.options()` (e.g., `"NLOPT_LN_NELDERMEAD"`,
`"NLOPT_LN_COBYLA"`, `"NLOPT_LN_SBPLX"`).  Provides finer tolerance control
(`ftol_abs`, `ftol_rel`, `xtol_abs`, `xtol_rel`, `maxeval`) compared to
`minqa::bobyqa`.

### `nlminbwrap(par, fn, lower, upper, control = list(), ...)`

Wrapper for base R `nlminb()`.  Provides `nlminb` access without requiring
`optimx`.  Gabor Grothendieck's implementation.

### Built-in `"Nelder_Mead"` / `"bobyqa"`

`"Nelder_Mead"` uses lme4's internal implementation (box-constrained).
`"bobyqa"` uses `minqa::bobyqa`.  Both accept `maxfun` in `optCtrl`.

### External: `optimx` package

Setting `optimizer = "optimx"` and supplying `method` in `optCtrl` unlocks
L-BFGS-B, nlminb, and other optimizers from the `optimx` package.  Users
rely on this for gradient-based optimizers when the default derivative-free
methods fail or are too slow.

---

## 4. `allFit()` — multi-optimizer robustness check

```
allFit(object, meth.tab = NULL, data = NULL, verbose = TRUE,
       show.meth.tab = FALSE, maxfun = 1e5,
       parallel = c("no", "multicore", "snow"),
       ncpus = 1L, cl = NULL,
       catch.errs = TRUE, start_from_mle = TRUE)
```

Re-fits a fitted `[g]lmerMod` with all available optimizers.  The gold
standard for deciding whether a convergence warning is a false positive.
If all optimizers agree to within practical tolerance, the warning is
spurious.

| Argument | What it does |
|---|---|
| `object` | A fitted `merMod`. |
| `meth.tab` | Custom table of `(method, optimizer)` pairs to override defaults. |
| `data` | Optional data to attach to the result for later debugging. |
| `verbose` | Print progress per optimizer. |
| `show.meth.tab` | Return the method table rather than fits. |
| `maxfun` | Max function evaluations; auto-converted per optimizer convention. |
| `parallel` | `"no"`, `"multicore"`, or `"snow"` parallelization. |
| `ncpus` / `cl` | Number of CPUs or a `snow` cluster object. |
| `catch.errs` | Wrap individual fits in `tryCatch` to skip optimizer errors. |
| `start_from_mle` | Initialize refits from the MLE of the original fit. |

`summary(allFit(...))` returns a list with `$which.OK`, `$llik`, `$fixef`,
`$sdcor`, and `$theta` slots for comparison.

---

## 5. Modular fit pipeline

lme4 exposes the four internal stages of an LMM fit as public functions.
Users rely on these for custom optimizers, profiling, deviance-function
inspection, and bootstrapping.

### `lFormula(formula, data, REML, subset, weights, na.action, offset, contrasts, control, ...)`

Parses formula and data, validates inputs, returns a list with components:
`fr` (model frame), `X` (fixed-effects design matrix), `reTrms` (random
terms info), `REML`.

### `mkLmerDevfun(fr, X, reTrms, REML, start, verbose, control, ...)`

Constructs the deviance function to be minimized.  Returns a closure whose
environment holds `merPredD` and `lmResp` reference-class objects.  **Deep
copies are needed** if the user wants to isolate environments.

### `optimizeLmer(devfun, optimizer, restart_edge, boundary.tol, start, verbose, control, ...)`

Runs the nonlinear optimizer over the theta parameter vector.  Returns an
optimization result list (`par`, `fval`, `feval`, `conv`, `message`).
Accepts `calc.derivs` and `use.last.params` via `...`.

### `mkMerMod(rho, opt, reTrms, fr, mc, lme4conv)`

Packages environment `rho` of the deviance function plus optimization
results into a `[g]lmerMod` object.

### `mkReTrms(bars, fr, drop.unused.levels, reorder.terms, reorder.vars)`

(Internal but documented) Constructs the random-effects term list from the
bars-list parsed from the formula and the model frame.  Users call this when
building completely custom fitting pipelines.

---

## 6. Convergence diagnostics

### `isSingular(object, tol = getSingTol())`

Returns `TRUE` if any random-effect variance is at or near zero (or
correlation at ±1), i.e. a boundary/singular fit.  Users call this to
classify boundary fits programmatically before inference.

### `getSingTol()` / `setSingTol()`

Get/set the global singularity tolerance used by `isSingular()` and
`check.conv.singular`.

### `checkConv(derivs, coefs, ctrl, lbound, debug)`

Internal function called post-optimization.  Checks gradient, Hessian, and
singularity.  Not typically called directly but its output is stored in
`fit@optinfo$conv`.

### Convergence help page (`?convergence`)

Documents the full recommended workflow: tighten tolerances, center/scale
predictors, Richardson-extrapolation Hessian check, restart from optimum,
use `allFit`.  Key `nloptwrap` tolerances are `ftol_abs` (1e-6) and
`xtol_abs` (1e-6); tightening to 1e-8 is the first recommended step.

---

## 7. lmerTest-specific fitting behaviors

lmerTest's `lmer()` returns an `lmerModLmerTest` object (S4, inherits from
`lmerMod`).  Additional slots store the Jacobians and other quantities needed
for Satterthwaite denominator-df computation.  All standard lme4 post-fit
methods work on `lmerModLmerTest`.

**Fallback**: if Satterthwaite df computation fails after a valid
`lmerMod` is produced, lmerTest returns a plain `lmerMod` (not
`lmerModLmerTest`), with a message.

No new fitting arguments are introduced; the Satterthwaite machinery is
computed automatically at fit time.

---

## 8. Formula syntax surface (within `formula` argument)

| Syntax | Meaning | Notes |
|---|---|---|
| `y ~ x + (1 \| g)` | Random intercept per group `g` | Standard random intercept |
| `y ~ x + (1 + x \| g)` | Correlated random intercept + slope | Unstructured (general PSD) covariance |
| `y ~ x + (x \| g)` | Correlated slope (implicit intercept) | Equivalent to above with intercept |
| `y ~ x + (0 + x \| g)` | Random slope, no random intercept | Suppresses intercept in RE term |
| `y ~ x + (1 \|\| g)` | Uncorrelated intercept (diagonal) | Only for continuous predictors by default |
| `y ~ x + diag(x \| g)` | Diagonal RE covariance | Explicit diagonal structure (lme4 ≥ 2.0) |
| `y ~ x + cs(x \| g)` | Compound-symmetric RE covariance | All pairwise correlations equal (lme4 ≥ 2.0) |
| `y ~ x + ar1(x \| g)` | AR(1) RE covariance | Autoregressive order-1, homogeneous by default (lme4 ≥ 2.0) |
| `y ~ x + us(x \| g)` | Unstructured RE covariance | Explicit; same as default `(x \| g)` |
| `y ~ x + (1 \| g1/g2)` | Nested grouping | Expands to `(1\|g1) + (1\|g1:g2)` |
| `y ~ x + (1 \| g1) + (1 \| g2)` | Crossed random effects | Two independent grouping factors |
| `offset(z)` | Offset in formula | Equivalent to `offset` argument |
| `I(x^2)` / `poly(x, 2)` | Fixed-effect transformations | Standard R formula machinery |

Global option `lme4.doublevert.default = "diag_special"` extends `||` to
categorical predictors.

---

## 9. Update and refit methods

### `update.merMod(object, formula., ..., evaluate = TRUE)`

Re-fits the model with modified arguments.  Respects the stored `call`;
reliable only when the original fit included an explicit `data` argument.
Used pervasively in multi-optimizer robustness checks and REML ↔ ML switches.

### `refit(object, newresp, newweights, ...)`

Re-fits the same model structure with a new response vector (and optionally
new weights).  Used in parametric bootstrap loops.  More efficient than
`update()` when only the response changes.

---

## 10. Key post-fit control and diagnostic accessors

These are all ways users interact with the fit's control/convergence state.

| Function | What it returns |
|---|---|
| `fit@optinfo` | List with `optimizer`, `control`, `derivs` (gradient and Hessian), `warnings`, `conv` (convergence codes and messages), `feval` (function evaluations) |
| `fit@optinfo$derivs` | List with `gradient` and `Hessian` at the optimum (if computed) |
| `fit@optinfo$conv` | Convergence check output: `lme4` (optimizer code + messages) and `opt` (named list from `checkConv`) |
| `getME(fit, "theta")` | Estimated theta vector (Cholesky factor elements) |
| `getME(fit, "par")` | Alias for `theta` |
| `getME(fit, "X")` | Fixed-effects design matrix (with `col.dropped` attribute if rank-deficient) |
| `getME(fit, "Z")` | Random-effects design matrix |
| `getME(fit, "Lambda")` | Relative covariance factor |
| `isSingular(fit)` | `TRUE` if boundary/singular |
| `deviance(fit)` | (Restricted) deviance at optimum |

---

## 11. Convergence controls summary (nloptwrap defaults)

| Control | Default | Effect |
|---|---|---|
| `ftol_abs` | 1e-6 | Stop on absolute change in deviance |
| `ftol_rel` | 0 | Stop on relative change in deviance (disabled) |
| `xtol_abs` | 1e-6 | Stop on absolute change in parameters |
| `xtol_rel` | 0 | Stop on relative change in parameters (disabled) |
| `maxeval` | 1000 | Max function evaluations |

For `bobyqa`: `rhobeg` (2e-3), `rhoend` (2e-7), `maxfun` (10000).  
For `Nelder_Mead`: `FtolAbs` (1e-5), `FtolRel` (1e-15), `XtolRel` (1e-7),
`maxfun` (10000).

---

## 12. Global options affecting fit behavior

| Option | Effect |
|---|---|
| `lmerControl` | Named list; sets defaults for `check.*` pre-fit options globally |
| `lme4.doublevert.default` | `"diag_special"` enables `||` for categorical predictors |
| `na.action` | Default na-action function if not passed explicitly |

---

*End of lme4/lmerTest LMM fitting & control surface reference.*
