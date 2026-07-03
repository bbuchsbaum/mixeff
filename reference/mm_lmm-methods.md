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

# S3 method for class 'mm_lmm'
fixef(object, ...)

# S3 method for class 'mm_glmm'
fixef(object, ...)

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
deviance(object, REML = NULL, ...)

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
model.frame(formula, ...)

# S3 method for class 'mm_glmm'
model.frame(formula, ...)

ngrps(object, ...)

# Default S3 method
ngrps(object, ...)

# S3 method for class 'mm_lmm'
ngrps(object, ...)

# S3 method for class 'mm_glmm'
ngrps(object, ...)

# S3 method for class 'mm_lmm'
weights(object, ...)

# S3 method for class 'mm_glmm'
weights(object, ...)

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

  Reserved for generic compatibility.

- condVar:

  Logical; when `TRUE`, Phase 2 returns the random-effects tables with
  an `NA` `postVar` array and an `mm_unavailable_reason` attribute
  rather than fabricating conditional variances.

- REML:

  Ignored; included for S3 compatibility with likelihood and deviance
  generics.

- k:

  Penalty per parameter for [`AIC()`](https://rdrr.io/r/stats/AIC.html).

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
#> Variance components:
#>  group        name variance std_dev correlation       note
#>      g (Intercept)        0       0             [boundary]
#> [boundary]: variance component is at the boundary of the parameter space.
#> Residual std. dev.: 0.861654
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
