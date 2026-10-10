# Refit a mixeff model with a new response

`refit()` fits the same model to a new response, keeping every other
setting of the original fit, like
[`lme4::refit()`](https://rdrr.io/pkg/lme4/man/refit.html). For an
`mm_lmm` it calls
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) with the
stored model frame, `REML` setting and prior weights; for an `mm_glmm`
it calls
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) with
the stored family (negative-binomial fits re-estimate theta when the
original did), prior weights, offset, `method`, `nAGQ` and control
settings (verbosity silenced).

## Usage

``` r
refit(object, newresp, ...)

# Default S3 method
refit(object, newresp, ...)

# S3 method for class 'mm_lmm'
refit(object, newresp, ...)

# S3 method for class 'mm_glmm'
refit(object, newresp, ...)
```

## Arguments

- object:

  A fitted `mm_lmm` or `mm_glmm`.

- newresp:

  The new response (see Details).

- ...:

  `control =` to override the stored control settings.

## Value

A new fit of the same class as `object`.

## Details

`newresp` may be a numeric vector, or a one-column data frame such as
the output of [`simulate()`](https://rdrr.io/r/stats/simulate.html) with
`nsim = 1`. For binomial GLMMs it may also be a two-column
`(successes, failures)` matrix (the shape
[`simulate()`](https://rdrr.io/r/stats/simulate.html) returns for a
[`cbind()`](https://rdrr.io/r/base/cbind.html) response; the trial
counts become the new prior weights), a two-level factor, or a logical
vector.

As in lme4, `newresp` is on the scale of the model's response: for a
transformed response such as `log(y) ~ ...` it replaces `log(y)`, and
the refit's formula names it `.mm_refit_response`.

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rpois(60, 3), x = rnorm(60),
  g = factor(rep(seq_len(10), each = 6))
)
fit <- glmm(y ~ x + (1 | g), df, family = poisson(),
            control = mm_control(verbose = -1))
ystar <- simulate(fit, nsim = 1, seed = 2)
fixef(refit(fit, ystar))
#> (Intercept)           x 
#>  1.06657255  0.08442761 
```
