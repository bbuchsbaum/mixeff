# lme4 GLMM Fitting — Public Surface Reference

**Surveyed packages:** lme4 2.0.1 (primary), lmerTest 3.2.1 (noted where relevant)
**Survey date:** 2026-05-31
**Purpose:** Exhaustive reference of every user-facing function, argument, and behavior
in the lme4 GLMM fitting family. This document is the reference against which mixeff
is judged; it does NOT assess mixeff.

---

## 1. Top-level fitting functions

### 1.1 `glmer()`

The primary GLMM fitting function. Fits a generalized linear mixed-effects model via
maximum likelihood. The linear predictor is related to the conditional mean of the
response through the inverse link function defined in the GLM family.

```r
glmer(
  formula,
  data       = NULL,
  family     = gaussian,
  control    = glmerControl(),
  start      = NULL,
  verbose    = 0L,
  nAGQ       = 1L,
  subset,
  weights,
  na.action,
  offset,
  contrasts  = NULL,
  mustart,
  etastart,
  devFunOnly = FALSE
)
```

| Argument | Type | Purpose |
|---|---|---|
| `formula` | two-sided formula | lme4-style: `response ~ fixed + (re | group)`. Random-effects terms separated by `|`. |
| `data` | data.frame | Source data. Strongly recommended; required for `update()` / `drop1()` reliability. |
| `family` | family object or constructor | GLM family. See §2 for full family/link matrix. |
| `control` | `glmerControl()` list | Optimizer, convergence checks, tolerances. See §4. |
| `start` | numeric vector or named list | Starting values. Named list may have components `par`/`theta` (covariance parameters, Cholesky factor of relative covariance) and `fixef`/`beta` (fixed-effect coefficients). When both stages are run, stage-2 start extends stage-1 optimum with fixed effects. |
| `verbose` | integer (0/1/2) | 0 = silent; 1 = outer optimizer trace; 2 = PIRLS step trace. |
| `nAGQ` | integer ≥ 0 | Likelihood approximation: 0 = fast PQL-like (nAGQ=0 optimizes RE+FE jointly in PIRLS); 1 = Laplace (default); >1 = adaptive Gauss-Hermite quadrature. **Constraint:** nAGQ > 1 requires a model with exactly one scalar (intercept-only) random-effects term. Max value is 100. |
| `subset` | expression or vector | Row subset; same semantics as `lm()`. |
| `weights` | numeric vector | Prior weights (observation-level). For binomial CBind responses, encodes the denominator instead of `cbind()`. |
| `na.action` | function | Default `na.omit`; strips observations with any NA. |
| `offset` | numeric vector | A priori known additive component to the linear predictor. Can also appear as `offset(x)` in the formula. Multiple `offset()` terms are summed. |
| `contrasts` | named list | Contrast specifications for factor columns; passed to `model.matrix`. |
| `mustart` | numeric vector | Starting values on the scale of the conditional mean (passed to PIRLS initialization, like `glm`). |
| `etastart` | numeric vector | Starting values on the linear predictor scale. |
| `devFunOnly` | logical | If `TRUE`, return only the deviance evaluation function rather than a fitted model. The deviance function closes over its environment; subsequent calls are not guaranteed identical but are within machine tolerance. |

**Returns:** Object of class `glmerMod` (subclass of `merMod`).

**Two-stage optimization:** By default, `glmer()` runs two optimization stages. Stage 1 uses `nAGQ = 0` (fast PQL-like), and Stage 2 refines with the requested `nAGQ` (default 1 = Laplace). This two-stage approach can be disabled via `control = glmerControl(nAGQ0initStep = FALSE)`.

---

### 1.2 `glmer.nb()`

Fits a negative-binomial GLMM with the overdispersion parameter `theta` estimated
jointly with the other model parameters, using a profile likelihood (golden-section
search) over `theta`.

```r
glmer.nb(
  ...,
  interval   = log(th) + c(-3, 3),
  tol        = 5e-05,
  verbose    = FALSE,
  nb.control = NULL,
  initCtrl   = list(limit = 20, eps = 2 * tol, trace = verbose, theta = NULL)
)
```

| Argument | Purpose |
|---|---|
| `...` | All `glmer()` arguments (formula, data, control, start, nAGQ, subset, weights, na.action, offset, contrasts, mustart, etastart, devFunOnly). |
| `interval` | Search interval for `log(theta)` in golden-section search. Default is centered on an initial estimate. |
| `tol` | Convergence tolerance for `theta` estimation. |
| `verbose` | Trace the `theta` search. |
| `nb.control` | A `glmerControl()` list used for the inner fits. |
| `initCtrl` | Controls for the initial `theta` estimation step. `theta = NULL` means estimate from data; supply a numeric to use a fixed starting value. |

**Note:** The negative-binomial dispersion parameter `theta` is accessible via
`getME(fit, "glmer.nb.theta")`. The variance function is `mu + mu^2/theta`.

---

## 2. Family/link matrix

lme4 accepts any base-R GLM family object or constructor, plus its own
`negative.binomial()`. The certified supported combinations (those tested and used
in lme4's own examples and tests) are:

| Family | Default link | Additional valid links | Notes |
|---|---|---|---|
| `binomial` | `logit` | `probit`, `cauchit`, `log`, `cloglog` | Binary 0/1 or cbind(successes, failures). |
| `poisson` | `log` | `identity`, `sqrt` | Count data. |
| `Gamma` | `inverse` | `identity`, `log` | Positive continuous. |
| `gaussian` | `identity` | `log`, `inverse` | Continuous; same as LMM when using identity link. |
| `inverse.gaussian` | `1/mu^2` | `inverse`, `identity`, `log` | Rare; positive continuous. |
| `negative.binomial(theta)` | `log` | Custom links passed to `MASS::negative.binomial` | Overdispersed counts; `theta` is the dispersion parameter. |

**Quasi-families (`quasi`, `quasibinomial`, `quasipoisson`) are NOT supported** in
`glmer()` — lme4 explicitly throws an error: `"quasi" families cannot be used in glmer`.

**Custom link objects:** Any valid `stats::make.link()` link object can be passed as
the `link` argument to a family constructor. The GLM family mechanism passes it
through to lme4's PIRLS machinery without additional validation.

---

## 3. `cbind()` binomial response

For grouped binomial data, the response is specified as a two-column matrix:

```r
glmer(cbind(successes, failures) ~ x + (1 | group), data, family = binomial)
```

- Left column: number of successes.
- Right column: number of failures (not total trials).
- Total trials = successes + failures.
- Equivalent to using `weights = n` with a proportion response, but the cbind form
  is preferred for clarity.
- nAGQ > 1 is permitted when there is a single scalar RE term.

---

## 4. `glmerControl()` — optimizer and convergence control

```r
glmerControl(
  optimizer              = c("bobyqa", "Nelder_Mead"),
  restart_edge           = FALSE,
  boundary.tol           = 1e-05,
  calc.derivs            = NULL,
  use.last.params        = FALSE,
  sparseX                = FALSE,
  standardize.X          = FALSE,
  autoscale              = NULL,
  check.nobs.vs.rankZ    = "ignore",
  check.nobs.vs.nlev     = "stop",
  check.nlev.gtreq.5     = "ignore",
  check.nlev.gtr.1       = "stop",
  check.nobs.vs.nRE      = "stop",
  check.rankX            = c("message+drop.cols", "silent.drop.cols",
                              "warn+drop.cols", "stop.deficient", "ignore"),
  check.scaleX           = c("warning", "stop", "silent.rescale",
                              "message+rescale", "warn+rescale", "ignore"),
  check.formula.LHS      = "stop",
  check.conv.nobsmax     = 10000,
  check.conv.nparmax     = 20,
  check.conv.grad        = .makeCC("warning", tol = 0.002, relTol = NULL),
  check.conv.singular    = .makeCC(action = "message", tol = getSingTol()),
  check.conv.hess        = .makeCC(action = "warning", tol = 1e-06),
  optCtrl                = list(),
  mod.type               = "glmer",
  tolPwrss               = 1e-07,
  compDev                = TRUE,
  nAGQ0initStep          = TRUE,
  check.response.not.const = "stop"
)
```

Key parameters:

| Parameter | Purpose |
|---|---|
| `optimizer` | Default is `c("bobyqa", "Nelder_Mead")`. For nAGQ > 0 stage, `Nelder_Mead` is used; for nAGQ = 0 stage, `bobyqa`. Can be set to a single optimizer name or a function. |
| `restart_edge` | Whether to restart from boundary when optimizer hits a boundary. |
| `boundary.tol` | Tolerance for declaring a parameter at a boundary. |
| `calc.derivs` | Whether to compute gradient/Hessian at convergence. `NULL` = use package default. |
| `use.last.params` | Use the parameter values at the last PIRLS step if convergence fails. |
| `tolPwrss` | Convergence tolerance for the penalized weighted residual sum of squares in PIRLS. |
| `compDev` | Whether to use compiled deviance function (C++ vs R). |
| `nAGQ0initStep` | Whether to run the nAGQ=0 stage first before the main nAGQ stage. |
| `check.conv.*` | Per-criterion convergence check actions: `"stop"`, `"warning"`, `"message"`, `"ignore"`. |
| `check.nobs.vs.*` | Data adequacy checks: observations vs rank, levels, random effects. |
| `optCtrl` | Named list passed directly to the optimizer (e.g., `list(maxfun = 2e5)`). |
| `sparseX` | Use sparse fixed-effect design matrix. |
| `standardize.X` | Standardize predictors before fitting. |

---

## 5. Modular fitting interface

lme4 exposes a low-level modular API for fine-grained control over GLMM fitting:

### 5.1 `glFormula()` — formula processing

```r
glFormula(
  formula, data = NULL, family = gaussian,
  subset, weights, na.action, offset, contrasts = NULL,
  start, mustart, etastart, control = glmerControl(), ...
)
```

Parses the formula, builds the fixed-effect design matrix (`X`), the random-effects
terms (`reTrms`), and the response object. Returns a named list with components
`fr`, `X`, `reTrms`, `family`.

### 5.2 `mkGlmerDevfun()` — deviance function construction

```r
mkGlmerDevfun(
  fr, X, reTrms, family,
  nAGQ    = if (control$nAGQ0initStep) 0L else 1L,
  verbose = 0L,
  maxit   = 100L,
  control = glmerControl(), ...
)
```

Constructs the GLMM deviance function. The `nAGQ` argument here controls which
approximation the deviance function evaluates.

### 5.3 `updateGlmerDevfun()` — switch approximation after stage 1

```r
updateGlmerDevfun(devfun, reTrms, nAGQ = 1L)
```

Takes a stage-1 deviance function (typically fitted with nAGQ=0) and updates it to
use the requested `nAGQ` (Laplace or AGQ) for stage 2.

### 5.4 `optimizeGlmer()` — optimization

```r
optimizeGlmer(
  devfun,
  optimizer    = if (nAGQ > 0L) "Nelder_Mead" else "bobyqa",
  restart_edge = FALSE,
  boundary.tol = formals(glmerControl)$boundary.tol,
  verbose      = 0L,
  control      = list(),
  nAGQ         = if (missing(stage) || stage == 1L) 0L else 1L,
  stage,
  start        = NULL, ...
)
```

Runs the optimizer on the deviance function. The `stage` argument distinguishes
stage-1 (covariance + FE) from stage-2 (covariance only at fixed FE) optimization.

### 5.5 `mkMerMod()` — model object assembly

```r
mkMerMod(rho, opt, reTrms, fr, mc, lme4conv = NULL)
```

Assembles a `merMod` object from the optimization result and model components.

---

## 6. Quadrature: `GHrule()` and `GQN` / `GQdk`

### 6.1 `GHrule()`

```r
GHrule(ord, asMatrix = TRUE)
```

Returns the Gauss-Hermite quadrature nodes and weights for `ord` points.
`ord` must be between 1 and 100. Used internally by `updateGlmerDevfun()` to
pre-compute the AGQ grid when `nAGQ > 1`.

Users can call `GHrule()` directly to inspect or pre-compute quadrature grids.

### 6.2 `GQN`, `GQdk`

Pre-computed tables of quadrature nodes (`GQN`) and weights (`GQdk`) exported from
lme4 for reference. These are the same grids used by `glmer()` internally.

---

## 7. Post-fit extraction methods on `glmerMod`/`merMod`

All methods below apply to objects returned by `glmer()` and `glmer.nb()`.

### 7.1 Fixed effects

| Function | Signature | Returns |
|---|---|---|
| `fixef(object, add.dropped=FALSE, noScale=NULL, ...)` | Extract fixed-effect coefficient vector. `add.dropped=TRUE` re-inserts dropped collinear columns as NA. | Named numeric vector. |
| `coef(object, ...)` | Combined fixed + conditional RE coefficients per group level. | List of data frames. |
| `vcov(object, ...)` | Variance-covariance matrix of fixed effects. | Matrix. |
| `summary(object, correlation=..., use.hessian=NULL, ...)` | Full model summary including fixed-effect table (Wald z-values, p-values for GLMM), random-effect table, deviance, AIC, BIC. | `summary.merMod` object. |

### 7.2 Random effects

| Function | Signature | Returns |
|---|---|---|
| `ranef(object, condVar=TRUE, drop=FALSE, whichel=..., postVar=FALSE, ...)` | Extract random-effect BLUPs. `condVar=TRUE` attaches conditional variances as attributes (used by `dotplot`). `postVar` is a deprecated alias for `condVar`. | Named list of data frames. |
| `VarCorr(x, sigma=1, ...)` | Extract variance-covariance components of the random effects. | `VarCorr.merMod` object (printable, plottable). |
| `sigma(object, ...)` | Residual standard deviation. For GLMMs with canonical link, typically 1 (fixed dispersion). For Gamma/inverse.gaussian, returns the estimated dispersion. | Scalar numeric. |

### 7.3 Fitted values and predictions

| Function | Signature | Returns |
|---|---|---|
| `fitted(object, ...)` | In-sample fitted values on the response scale. | Named numeric vector. |
| `residuals(object, type=c("deviance","pearson","working","response","partial"), scaled=FALSE, ...)` | Residuals. Default `type="deviance"` for GLMMs. `scaled=TRUE` divides by `sigma`. | Named numeric vector. |
| `predict(object, newdata=NULL, newparams=NULL, re.form=NULL, random.only=FALSE, terms=NULL, type=c("link","response"), allow.new.levels=FALSE, na.action=na.pass, se.fit=FALSE, ...)` | Predictions. `re.form=NULL` uses all RE; `re.form=NA` or `~0` gives population-level (FE-only); partial RE via one-sided formula. `type="link"` returns linear predictor; `type="response"` applies inverse link. `allow.new.levels=TRUE` uses zero for unseen group levels. `se.fit=TRUE` returns SE of prediction (via delta method). | Numeric vector or list. |

### 7.4 Model structure

| Function | Signature | Returns |
|---|---|---|
| `formula(object, ...)` | Extract model formula. | Formula. |
| `terms(object, ...)` | Terms object. | `terms`. |
| `model.frame(object, ...)` | Recover the model frame. | Data frame. |
| `model.matrix(object, ...)` | Fixed-effect design matrix. | Matrix. |
| `getData(object, ...)` | Alias for `model.frame`. | Data frame. |
| `nobs(object, ...)` | Number of observations used in the fit. | Integer. |
| `ngrps(object, ...)` | Number of levels per grouping factor. | Named integer vector. |
| `df.residual(object, ...)` | Residual degrees of freedom (n - p). | Integer. |
| `deviance(object, REML=NULL, ...)` | Model deviance (−2 log-likelihood). | Scalar. |
| `logLik(object, REML=NULL, ...)` | Log-likelihood. For GLMMs always ML (REML not applicable). | `logLik` object. |
| `AIC(object, ...)` / `BIC(object, ...)` | Information criteria. | Scalar. |
| `family(object, ...)` | Extract the GLM family object used in the fit. | `family` object. |

### 7.5 Lower-level matrix/parameter access: `getME()`

```r
getME(object, name, ...)
```

Extracts named internal components. Valid `name` values for GLMMs:

| Name | Description |
|---|---|
| `"X"` | Fixed-effect design matrix |
| `"Z"` | Random-effect design matrix (transposed Zt, transposed) |
| `"Zt"` | Transposed random-effect design matrix |
| `"Ztlist"` | List of per-term transposed RE design matrices |
| `"mmList"` | Per-term model matrices |
| `"y"` | Response vector |
| `"mu"` | Conditional mean vector |
| `"u"` | Random-effect modes (spherical) |
| `"b"` | Random-effect modes (original scale) |
| `"L"` | Cholesky factor of the system matrix |
| `"Lambda"` | RE relative covariance factor (lower triangular) |
| `"Lambdat"` | Transpose of Lambda |
| `"A"` | Lambda %*% Zt |
| `"Lind"` | Index vector mapping theta to Lambda entries |
| `"RX"` | Cholesky factor of FE precision matrix |
| `"RZX"` | Cross-factor matrix |
| `"sigma"` | Residual standard deviation |
| `"Gp"` | Group pointer vector |
| `"Tp"` | Theta pointer vector |
| `"flist"` | List of grouping factors |
| `"fixef"` | Fixed-effect coefficients |
| `"beta"` | Fixed-effect coefficients (alias) |
| `"theta"` | Covariance parameter vector (Cholesky parameterization) |
| `"par"` | Full parameter vector as used in optimization |
| `"ST"` | List of Cholesky factor blocks per term |
| `"Tlist"` | List of Lambda submatrices per term |
| `"n_rtrms"` | Number of RE terms |
| `"n_rfacs"` | Number of RE grouping factors |
| `"N"` | Total number of observations (before subsetting) |
| `"n"` | Number of rows used in fit |
| `"p"` | Number of fixed-effect parameters |
| `"q"` | Total number of random-effect parameters |
| `"p_i"` | Number of FE params per RE term |
| `"l_i"` | Number of levels per grouping factor |
| `"q_i"` | Random effects per term (p_i * l_i) |
| `"k"` | Number of RE terms |
| `"m_i"` | Number of covariance parameters per term |
| `"m"` | Total number of covariance parameters (nth) |
| `"cnms"` | Column names of random-effects terms |
| `"devcomp"` | Deviance components list |
| `"offset"` | Offset vector |
| `"lower"` | Lower bounds for covariance parameters |
| `"devfun"` | Deviance function (reconstructed) |
| `"devarg"` | Parameter vector at optimum |
| `"glmer.nb.theta"` | Negative-binomial dispersion `theta` (NA if not NB model) |
| `"ALL"` | Returns all extractable components as a list |

---

## 8. Model comparison and testing

### 8.1 `anova()` for model comparison

```r
anova(object, ..., refit = TRUE, model.names = NULL)
```

- Compares nested or non-nested `merMod` objects.
- `refit = TRUE` (default): refits REML models as ML before comparison (for LMMs). For GLMMs, which are always ML, this is effectively a no-op.
- Returns a table with Df, AIC, BIC, logLik, deviance, Chisq, Df difference, p-value.
- For single-model `anova(fit)`, returns a type-III-style Wald table.

### 8.2 `drop1()` for term deletion

```r
drop1(object, scope, test = "Chisq", ...)
```

Drops each eligible fixed-effect term and performs an LRT. `test = "Chisq"` uses chi-squared test.

### 8.3 `isSingular()` — boundary detection

```r
isSingular(object, tol = getSingTol())
```

Tests whether the fitted model is on the boundary (singular fit): a random-effect
variance is effectively zero or the correlation is ±1. `getSingTol()` returns the
default tolerance (1e-4).

---

## 9. Confidence intervals: `confint()`

```r
confint(
  object, parm, level = 0.95,
  method    = c("profile", "Wald", "boot"),
  zeta, nsim = 500,
  boot.type  = c("perc", "basic", "norm"),
  FUN        = NULL, quiet = FALSE,
  oldNames, signames = TRUE,
  boot.scale = c("sdcor", "vcov"), ...
)
```

| Method | Description |
|---|---|
| `"profile"` | Profile-likelihood CIs. Computationally expensive; requires `profile.merMod`. For GLMMs, this is the most accurate method. |
| `"Wald"` | Wald CIs using `±z * SE`. Fast but approximate; based on asymptotic normality. |
| `"boot"` | Parametric or semiparametric bootstrap CIs via `bootMer`. |

Parameters:
- `parm`: subset of parameters to compute CIs for.
- `zeta`: zeta scale value for profile CIs.
- `nsim`: number of bootstrap replicates.
- `boot.type`: bootstrap CI type (`"perc"` = percentile, `"basic"`, `"norm"`).
- `signames`: use sigma-parameterization names (`.sig01`, `.sigma`) vs verbose names.
- `boot.scale`: parameterization for boot CIs (`"sdcor"` = SD/correlation, `"vcov"` = variance-covariance).

---

## 10. Profile likelihood: `profile()`

```r
profile(
  fitted, which = NULL, alphamax = 0.01, maxpts = 100,
  delta = NULL, delta.cutoff = 1/8,
  verbose = 0, devtol = 1e-09, devmatchtol = 1e-05,
  maxmult = 10, startmethod = "prev",
  optimizer = NULL, control = NULL,
  signames = TRUE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("profile.ncpus", 1L),
  cl = NULL,
  prof.scale = c("sdcor", "varcov"), ...
)
```

Computes the profile likelihood over each parameter, used for profile CIs and
zeta plots. Key parameters:
- `which`: which parameters to profile (by index or name).
- `alphamax`: maximum alpha for the CI extent.
- `maxpts`: maximum profile points per parameter.
- `parallel`: parallel computation backend.
- `prof.scale`: whether to profile on SD/correlation or variance/covariance scale.

**Note for GLMMs:** Profile CIs for GLMMs are computationally feasible but slow
because each profile point requires a GLMM refit.

---

## 11. Bootstrap: `bootMer()`

```r
bootMer(
  x, FUN, nsim = 1, seed = NULL,
  use.u   = FALSE,
  re.form = NA,
  type    = c("parametric", "semiparametric"),
  verbose = FALSE,
  .progress = "none",
  PBargs = list(),
  parallel = c("no", "multicore", "snow"),
  ncpus   = getOption("boot.ncpus", 1L),
  cl      = NULL
)
```

| Argument | Purpose |
|---|---|
| `x` | Fitted `merMod` object. |
| `FUN` | Function applied to each bootstrap replicate; returns a numeric vector. |
| `nsim` | Number of bootstrap replicates. |
| `seed` | Random seed for reproducibility. |
| `use.u` | Whether to condition on observed random effects (`TRUE`) or marginalize (`FALSE`, default). |
| `re.form` | Random-effects formula for simulation. `NA` = marginalize over RE; `NULL` = include all RE. |
| `type` | `"parametric"` simulates from the fitted model; `"semiparametric"` resamples residuals. |
| `parallel` | Parallelism backend. |

`bootMer()` is the engine behind `confint(method="boot")` and `simulate.merMod()`.

---

## 12. Simulation: `simulate()`

```r
simulate(
  object, nsim = 1, seed = NULL,
  use.u   = FALSE,
  re.form = NA,
  newdata = NULL,
  newparams = NULL,
  family  = NULL,
  cluster.rand = rnorm,
  allow.new.levels = FALSE,
  na.action = na.pass, ...
)
```

Simulates responses from the fitted model.
- `use.u = FALSE` marginalizes over RE; `TRUE` conditions on estimated RE BLUPs.
- `re.form` controls which RE to include.
- `newdata` / `newparams` allow simulation at new covariate values or parameter values.
- `family` can override the family for simulation (e.g., simulate overdispersed data from a Poisson fit).
- `cluster.rand`: function to draw cluster-level random effects (default `rnorm`).

---

## 13. Optimizer utilities

### 13.1 Built-in optimizers

| Optimizer | Function | Notes |
|---|---|---|
| `bobyqa` | `nlminbwrap` wrapping `minqa::bobyqa` | Default for nAGQ=0 stage; handles box constraints. |
| `Nelder_Mead` | `Nelder_Mead()` (C++ implementation) | Default for nAGQ>0 stage. |
| `nloptwrap` | `nloptwrap()` | Wraps NLopt algorithms (requires `nloptr` package). |
| `nlminbwrap` | `nlminbwrap()` | Wraps `stats::nlminb`. |
| `golden` | Internal | Golden-section search used in `glmer.nb`. |

### 13.2 `allFit()`

```r
allFit(
  object, meth.tab = NULL, data = NULL,
  verbose = TRUE, show.meth.tab = FALSE,
  maxfun = 1e5,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("allFit.ncpus", 1L),
  cl = NULL,
  catch.errs = TRUE,
  start_from_mle = TRUE
)
```

Re-fits the model with all available optimizers and returns a list of fits.
Useful for convergence diagnostics: if all optimizers agree, the fit is reliable.
`catch.errs = TRUE` (default) prevents failed optimizer fits from aborting the loop.

---

## 14. Diagnostics and model utilities

### 14.1 `rePCA()` — random-effects PCA

```r
rePCA(x)
```

Performs PCA on the random-effects covariance structure to diagnose singularity.
Near-zero eigenvalues indicate redundant RE dimensions. Complements `isSingular()`.

### 14.2 `devfun2()`

```r
devfun2(fm, useSc = if (isLMM(fm)) TRUE else NA, scale = c("sdcor", "varcov"), ...)
```

Returns a version of the deviance function reparameterized on SD/correlation
or variance/covariance scale, for profiling and inspection.

### 14.3 `refit()` / `refitML()`

```r
refit(object, newresp, ...)
refitML(x, ...)
```

- `refit(object, newresp)`: refit the model with a new response vector (useful in
  bootstrap loops and simulation studies).
- `refitML(x)`: refit a REML-fitted LMM using ML. For GLMMs this is effectively
  a no-op since GLMMs always use ML.

### 14.4 `checkConv()` — convergence check

```r
checkConv(derivs, coefs, ctrl, lbound, debug = FALSE)
```

Checks gradient, Hessian, and boundary conditions at the optimum.
Returns a list of messages and action codes. Called internally at the end of fitting.

---

## 15. Weights and dispersion

### 15.1 Prior weights

`weights` in `glmer()` takes a numeric vector of prior weights. For standard GLMMs:
- **Binomial:** weights specify the denominator when the response is a proportion
  (alternative to `cbind()`). When using `cbind()`, the denominator is extracted
  automatically and no additional `weights` argument is needed.
- **Poisson/Gamma:** weights act as the standard GLM prior weights, scaling the
  contribution of each observation to the log-likelihood.
- `weights(fitted_model, type = c("prior", "working"))`: extracts prior or PIRLS
  working weights from a fitted model.

### 15.2 Dispersion

For families with fixed dispersion (binomial, Poisson), `sigma(fit)` returns 1.
For Gamma and inverse.gaussian, `sigma(fit)` returns the estimated dispersion
parameter (reciprocal of the GLM shape parameter). For negative binomial,
the extra dispersion is captured by `getME(fit, "glmer.nb.theta")`.

There is **no free overdispersion parameter** in `glmer()` for binomial or Poisson
families. Overdispersion is modeled by adding an observation-level random effect:
```r
data$obs <- seq_len(nrow(data))
glmer(y ~ x + (1|group) + (1|obs), family = poisson)
```

---

## 16. Type-checking predicates

| Function | Purpose |
|---|---|
| `isGLMM(object)` | Returns `TRUE` if object is a GLMM (i.e., glmerMod). |
| `isLMM(object)` | Returns `TRUE` if object is an LMM. |
| `isNLMM(object)` | Returns `TRUE` if object is a nonlinear MM. |
| `isREML(object)` | Returns `TRUE` if fitted by REML (always `FALSE` for GLMMs). |
| `isSingular(object, tol)` | Returns `TRUE` if fit is on the boundary. |

---

## 17. lmerTest relevance to GLMM

**lmerTest does not extend `glmer()`.** Its functionality is entirely restricted to
LMMs fitted by `lmer()`. Specifically:
- `lmer()` in lmerTest returns a `lmerModLmerTest` object (subclass of `lmerMod`).
- `step()`, `ranova()`, `contest()`, `contestMD()`, `contest1D()`, `ls_means()`,
  `calcSatterth()`, `show_tests()` — all operate only on `lmerMod`/`lmerModLmerTest`.
- `glmer()` results remain plain `glmerMod`; lmerTest adds no methods for them.

For GLMM inference, users rely on:
- Wald z-tests in `summary(glmer_fit)`.
- LRT via `anova(fit1, fit2)`.
- Profile CIs via `confint(fit, method="profile")`.
- Bootstrap via `confint(fit, method="boot")` or `bootMer()`.

---

## 18. `negative.binomial()` family constructor

```r
negative.binomial(theta = stop("'theta' must be specified"), link = "log")
```

This is the lme4-provided (re-exported from MASS) negative binomial family
constructor. It creates a `family` object where variance = `mu + mu^2/theta`.
- `theta` must be specified explicitly when calling `glmer()` with a fixed `theta`.
- When `theta` is to be estimated, use `glmer.nb()` instead.

---

## 19. Offset specification

Offsets can be specified in two equivalent ways:
1. Via the `offset` argument: `glmer(..., offset = log(exposure))`.
2. Via `offset()` in the formula: `glmer(y ~ x + offset(log(exposure)) + (1|g), ...)`.
3. Multiple `offset()` terms in the formula are summed.

Offsets are stored in the model frame and accessible via `getME(fit, "offset")`.
They shift the linear predictor by a fixed amount and are commonly used in Poisson
rate models (e.g., modeling events per person-year).

---

## 20. Summary of exported symbols in the GLMM family

The following lme4 symbols are directly relevant to GLMM fitting:

**Primary fitting:**
`glmer`, `glmer.nb`, `negative.binomial`, `glmerControl`

**Modular fitting:**
`glFormula`, `mkGlmerDevfun`, `updateGlmerDevfun`, `optimizeGlmer`, `mkMerMod`

**Quadrature:**
`GHrule`, `GQN`, `GQdk`

**Post-fit extraction:**
`fixef`, `ranef`, `VarCorr`, `sigma`, `coef`, `vcov`,
`fitted`, `residuals`, `predict`, `weights`,
`logLik`, `deviance`, `AIC`, `BIC`, `nobs`, `ngrps`, `df.residual`,
`formula`, `terms`, `model.frame`, `model.matrix`, `getData`, `family`,
`getME`, `summary`, `print`

**Comparison/testing:**
`anova`, `drop1`, `isSingular`, `rePCA`

**Uncertainty:**
`confint`, `profile`, `bootMer`, `simulate`, `devfun2`

**Refitting:**
`refit`, `refitML`, `update`, `allFit`

**Type predicates:**
`isGLMM`, `isLMM`, `isREML`, `isNLMM`, `checkConv`

**Optimizer utilities:**
`nloptwrap`, `nlminbwrap`, `NelderMead`, `Nelder_Mead`, `golden`

**Family utilities:**
`glmFamily` (reference class wrapping a GLM family for PIRLS use)

---

*End of lme4 GLMM Fitting Surface Reference*
