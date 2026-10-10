# REML criterion of a linear mixed model

`REMLcrit()` mirrors
[`lme4::REMLcrit()`](https://rdrr.io/pkg/lme4/man/merMod-class.html):
the REML criterion (`-2` times the restricted log-likelihood) at the
fitted parameters. For a REML fit this is `-2 * logLik(fit)`. For an ML
fit it is the REML criterion evaluated at the ML estimates of the
covariance parameters, computed from the penalized least-squares
decomposition (`getME(fit, "devcomp")`), exactly as lme4 does.
[`lme4::REMLcrit()`](https://rdrr.io/pkg/lme4/man/merMod-class.html) is
a plain function rather than a generic, so mixeff exports its own; it
forwards `merMod` objects to lme4.

## Usage

``` r
REMLcrit(object)
```

## Arguments

- object:

  A fitted `mm_lmm` (or an lme4 `merMod`).

## Value

A single number.

## Examples

``` r
set.seed(1)
df <- data.frame(y = rnorm(60), x = rnorm(60),
                 g = factor(rep(seq_len(10), each = 6)))
fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
REMLcrit(fit)
#> [1] 155.3169
```
