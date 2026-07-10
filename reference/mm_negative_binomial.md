# Negative-binomial family for `glmm()`

NB2 family (variance `mu + mu^2/theta`, log link). With `theta = NULL`
(the default) the size parameter is estimated alongside the model,
matching
[`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html).
Supplying a positive `theta` fits conditional on that value, matching
`glmer(family = MASS::negative.binomial(theta))` — which
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) also
accepts directly.

## Usage

``` r
mm_negative_binomial(theta = NULL)
```

## Arguments

- theta:

  Optional positive NB2 size (dispersion) parameter. `NULL` estimates
  theta from the data.

## Value

A `family` object accepted by
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md).

## Examples

``` r
fam <- mm_negative_binomial()      # glmer.nb-style: theta estimated
fam_fixed <- mm_negative_binomial(theta = 2.5)
```
