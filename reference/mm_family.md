# Family objects and prior weights for mixeff fits

[`family()`](https://rdrr.io/r/stats/family.html) returns the R
[`family()`](https://rdrr.io/r/stats/family.html) object of a fit:
[`gaussian()`](https://rdrr.io/r/stats/family.html) for
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) fits and
the GLMM family (with its link) for
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits.
For negative-binomial fits it is `MASS::negative.binomial(theta)` at the
fitted (or fixed) theta, like
[`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html); when
MASS is not installed an equivalent family object with the same
variance, deviance and link functions is returned.

## Usage

``` r
# S3 method for class 'mm_lmm'
family(object, ...)

# S3 method for class 'mm_glmm'
family(object, ...)

# S3 method for class 'mm_lmm'
weights(object, type = c("prior", "working"), ...)

# S3 method for class 'mm_glmm'
weights(object, type = c("prior", "working"), ...)
```

## Arguments

- object:

  A fitted `mm_lmm` or `mm_glmm` object.

- ...:

  Unused.

- type:

  `"prior"` or `"working"` (GLMMs only).

## Value

A `family` object, or a numeric vector of weights.

## Details

[`weights()`](https://rdrr.io/r/stats/weights.html) mirrors
`lme4:::weights.merMod()`: `type = "prior"` returns the prior weights, a
vector of ones when the model was fitted without weights;
`type = "working"` returns the final PIRLS working weights of a GLMM.

## Examples

``` r
set.seed(1)
df <- data.frame(y = rnorm(60), x = rnorm(60),
                 g = factor(rep(seq_len(10), each = 6)))
fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
family(fit)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
head(weights(fit))
#> [1] 1 1 1 1 1 1
```
