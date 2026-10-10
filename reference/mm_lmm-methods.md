# Extract components from a fitted mixeff LMM

These methods provide the common lme4-style extractor surface for
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) fits.
The required values are stored directly on the R object or rebuilt
lazily from the serialized artifact, so these methods do not require a
live Rust handle after
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
[`readRDS()`](https://rdrr.io/r/base/readRDS.html).

`ngrps()` returns a named integer vector giving the number of levels of
each random-effect grouping factor, mirroring
[`lme4::ngrps()`](https://rdrr.io/pkg/lme4/man/ngrps.html).

Produces the long form returned by `as.data.frame(lme4::VarCorr(.))`:
one row per variance (`var2 = NA`) and one row per covariance (`var1`,
`var2` both set), with a final `Residual` row for LMMs. `vcov` holds the
(co)variance and `sdcor` the standard deviation (diagonal) or
correlation (off-diagonal). This is the shape
[`broom.mixed::tidy()`](https://generics.r-lib.org/reference/tidy.html)
expects.

Produces the long form returned by `as.data.frame(lme4::ranef(.))`:
columns `grpvar`, `term`, `grp`, `condval`, and `condsd`. `condsd` is
the conditional standard deviation, taken from the `postVar` attribute
when the modes were extracted with `condVar = TRUE`, and `NA` otherwise.

## Usage

``` r
fixef(object, ...)

# Default S3 method
fixef(object, ...)

# S3 method for class 'mm_lmm'
fixef(object, add.dropped = FALSE, ...)

# S3 method for class 'mm_glmm'
fixef(object, add.dropped = FALSE, ...)

ranef(object, ...)

# Default S3 method
ranef(object, ...)

# S3 method for class 'mm_lmm'
ranef(object, condVar = FALSE, ...)

# S3 method for class 'mm_glmm'
ranef(object, condVar = FALSE, ...)

# S3 method for class 'mm_lmm'
coef(object, ...)

# S3 method for class 'mm_glmm'
coef(object, ...)

VarCorr(x, ...)

# Default S3 method
VarCorr(x, ...)

# S3 method for class 'mm_lmm'
VarCorr(x, ...)

# S3 method for class 'mm_glmm'
VarCorr(x, ...)

# S3 method for class 'mm_lmm'
sigma(object, ...)

# S3 method for class 'mm_glmm'
sigma(object, ...)

# S3 method for class 'mm_lmm'
logLik(object, REML = NULL, ...)

# S3 method for class 'mm_glmm'
logLik(object, REML = NULL, ...)

# S3 method for class 'mm_lmm'
deviance(object, REML = NULL, ...)

# S3 method for class 'mm_glmm'
deviance(object, ...)

# S3 method for class 'mm_lmm'
AIC(object, ..., k = 2)

# S3 method for class 'mm_glmm'
AIC(object, ..., k = 2)

# S3 method for class 'mm_lmm'
BIC(object, ...)

# S3 method for class 'mm_glmm'
BIC(object, ...)

# S3 method for class 'mm_lmm'
nobs(object, ...)

# S3 method for class 'mm_glmm'
nobs(object, ...)

# S3 method for class 'mm_lmm'
df.residual(object, ...)

# S3 method for class 'mm_glmm'
df.residual(object, ...)

# S3 method for class 'mm_lmm'
formula(x, ...)

# S3 method for class 'mm_glmm'
formula(x, ...)

# S3 method for class 'mm_lmm'
model.frame(formula, fixed.only = FALSE, ...)

# S3 method for class 'mm_glmm'
model.frame(formula, fixed.only = FALSE, ...)

ngrps(object, ...)

# Default S3 method
ngrps(object, ...)

# S3 method for class 'mm_lmm'
ngrps(object, ...)

# S3 method for class 'mm_glmm'
ngrps(object, ...)

# S3 method for class 'mm_lmm'
extractAIC(fit, scale, k = 2, ...)

# S3 method for class 'mm_glmm'
extractAIC(fit, scale, k = 2, ...)

# S3 method for class 'mm_lmm'
terms(x, ...)

# S3 method for class 'mm_glmm'
terms(x, ...)

# S3 method for class 'mm_varcorr'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'mm_ranef'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'mm_lmm'
model.matrix(object, type = c("fixed", "random"), ...)

# S3 method for class 'mm_glmm'
model.matrix(object, type = c("fixed", "random"), ...)

# S3 method for class 'mm_lmm'
vcov(object, type = c("fixed", "theta"), correlation = FALSE, ...)

# S3 method for class 'mm_glmm'
vcov(object, type = c("fixed", "theta"), correlation = FALSE, ...)
```

## Arguments

- object, x, formula, fit:

  A fitted `mm_lmm` or `mm_glmm` object.

- ...:

  Reserved for generic compatibility. For
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html), further fitted models:
  with several models the result is the data frame
  [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) returns (`df` and
  `AIC`/`BIC`, rows named after the arguments).

- add.dropped:

  For `fixef()`: if `TRUE`, coefficients dropped for rank deficiency are
  included as `NA` (lme4's `add.dropped`); by default they are omitted,
  as in lme4.

- condVar:

  Logical; when `TRUE`, attach a `postVar` array of conditional
  variances to each random-effects table. Gaussian
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) fits
  receive the conditional covariances computed by the engine.
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits
  receive lme4's Laplace-approximation conditional covariances
  `sigma^2 * Lambda (Lambda'Z'WZ Lambda + I)^{-1} Lambda'` at the fitted
  modes, with `W` the final PIRLS working weights. A fit whose
  conditional-variance computation refuses receives an all-`NA`
  `postVar` plus an `mm_unavailable_reason` attribute rather than
  fabricated conditional variances.

- REML:

  Ignored; included for S3 compatibility with likelihood and deviance
  generics.

- k:

  Penalty per parameter for [`AIC()`](https://rdrr.io/r/stats/AIC.html).

- fixed.only:

  For [`model.frame()`](https://rdrr.io/r/stats/model.frame.html):
  return only the fixed-effect variables, as
  [`lme4::model.frame.merMod()`](https://rdrr.io/pkg/lme4/man/merMod-class.html).

- scale:

  Ignored; included for S3 compatibility with
  [`extractAIC()`](https://rdrr.io/r/stats/extractAIC.html).

- row.names, optional:

  Ignored; present for S3 consistency.

- type:

  For [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
  `"fixed"` returns the fixed-effect design matrix and `"random"`
  returns the sparse random-effect design matrix. For
  [`vcov()`](https://rdrr.io/r/stats/vcov.html), `"fixed"` returns the
  fixed-effect covariance surface and `"theta"` returns an unavailable
  theta-covariance matrix with a reason attribute.

- correlation:

  Logical; accepted for S3 compatibility with
  [`vcov()`](https://rdrr.io/r/stats/vcov.html).

## Value

A named integer vector of group counts.

## Details

These methods follow lme4's shapes.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns, per grouping
factor, every fixed-effect column (fixed-only columns repeated) plus the
conditional modes, in lme4's column order. `VarCorr()` returns lme4's
`VarCorr.merMod` structure: a named list of covariance matrices (one per
random-effect term, so `VarCorr(fit)$Subject` is a matrix) with
`"stddev"` and `"correlation"` attributes and list attributes `"sc"`
([`sigma()`](https://rdrr.io/r/stats/sigma.html)) and `"useSc"` (`TRUE`
for LMMs; `FALSE` for every GLMM, as in lme4 \>= 2.1-0, so GLMM
random-effect SDs are absolute and no Residual row is printed); it
prints like lme4 and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives
lme4's long form. mixeff's full-precision long table remains available
as `VarCorr(fit)$table` (with `$residual_sd`).
[`sigma()`](https://rdrr.io/r/stats/sigma.html) of a binomial, Poisson
or negative-binomial GLMM is 1, as in lme4 (the negative-binomial theta
is `getME(fit, "glmer.nb.theta")`). For GLMMs with a free dispersion
parameter (Gamma, inverse Gaussian, and Gaussian with a non-identity
link) it is `sqrt(phi)` with
`phi = sum(deviance residuals^2) / (n - rank([X, Z]))`, the dispersion
profiled during the fit, as in lme4 \>= 2.1-0 (`glmerControl()` defaults
`disp_method = "moment"`, `disp_dof_correction = TRUE`; lme4 \< 2.1-0
used a different, biased estimate).
`mm_control(disp_dof_correction = FALSE)` divides by `n` instead, and
`mm_control(max_phi_iter = )` caps the `phi` iterations. With
`mm_control(disp_method = "old/buggy")` the fit reproduces lme4 2.1-0's
`disp_method = "old/buggy"`: PIRLS runs with `phi = 1`, so
`weights(fit, "working")` does not carry `1/phi`, and `VarCorr()`
reports theta itself as the random-effect SD, as lme4 2.1-0 does. Under
lme4 \< 2.1-0 that SD was `sigma() * theta`.
[`deviance()`](https://rdrr.io/r/stats/deviance.html) of a GLMM is the
sum of squared deviance residuals, as
[`lme4::deviance.merMod()`](https://rdrr.io/pkg/lme4/man/deviance.html);
`-2 * logLik()` is the Laplace objective reported in
[`anova()`](https://rdrr.io/r/stats/anova.html).
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) returns
lme4's frame (transformed columns such as `log(y)`, a `terms` attribute,
`(weights)` and `(offset)` columns); the raw variables the engine used
are in `fit$model_frame`.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) of a joint-Laplace GLMM is
the inverse of the finite-difference Hessian over the fixed effects and
covariance parameters, glmer's default (`use.hessian = TRUE`). When that
Hessian is not positive definite (or unavailable) the engine falls back
to the fixed-effect block RX of the Laplace penalized least-squares
factorization at the optimum, conditional on the covariance parameters –
glmer's `vcov(fit, use.hessian = FALSE)`, which glmer also falls back to
with a warning. mixeff does the same: the matrix carries
`attr(, "mm_method") == "laplace_rx_conditional_on_theta"`, reliability
`"low"` and the engine's notes (`attr(, "mm_notes")`), and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) warns (class
`mm_vcov_rx_fallback`). Such standard errors ignore the uncertainty in
the covariance parameters.

Known lme4 2.1-0 discrepancy: for GLMMs with a free dispersion parameter
(Gamma, inverse Gaussian, Gaussian with a non-identity link) the
penalized least-squares factor already carries the `1 / phi` working
weights, so mixeff's RX-based
[`vcov()`](https://rdrr.io/r/stats/vcov.html) (`nAGQ = 0` fits, and the
fallback above) is `unsc()` with no further `sigma^2` factor. lme4 \>=
2.1-0 still multiplies by `sigma()^2` there (and in
`ranef(condVar = TRUE)`), so those values differ from lme4 2.1's by the
factor `phi = sigma()^2`; Hessian-based
[`vcov()`](https://rdrr.io/r/stats/vcov.html) values agree.

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rnorm(60), x = rnorm(60),
  g = factor(rep(seq_len(10), each = 6))
)
fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
fixef(fit)
#> (Intercept)           x 
#>   0.1122739  -0.0416571 
VarCorr(fit)
#>  Groups   Name        Std.Dev.           
#>  g        (Intercept) 0.00000  [boundary]
#>  Residual             0.86165            
#> [boundary]: variance component is at the boundary of the parameter space.
head(ranef(fit)$g)
#>   (Intercept)
#> 1           0
#> 2           0
#> 3           0
#> 4           0
#> 5           0
#> 6           0
sigma(fit)
#> [1] 0.8616536
logLik(fit)
#> 'log Lik.' -77.65847 (df=4)
nobs(fit)
#> [1] 60
```
