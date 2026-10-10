# Per-group linear (or generalized linear) model fits

`mm_lmlist()` is mixeff's counterpart of
[`lme4::lmList()`](https://rdrr.io/pkg/lme4/man/lmList.html): it fits
[`lm()`](https://rdrr.io/r/stats/lm.html) (or
[`glm()`](https://rdrr.io/r/stats/glm.html) when `family` is given)
separately within each level of a grouping factor, using lme4's
`response ~ predictors | group` formula. It is named `mm_lmlist()` so
attaching mixeff does not mask
[`lme4::lmList()`](https://rdrr.io/pkg/lme4/man/lmList.html) /
[`nlme::lmList()`](https://rdrr.io/pkg/nlme/man/lmList.html).

## Usage

``` r
mm_lmlist(
  formula,
  data,
  family = NULL,
  subset = NULL,
  weights = NULL,
  na.action = stats::na.omit,
  pool = is.null(family)
)
```

## Arguments

- formula:

  A formula `y ~ x | g`.

- data:

  A data frame.

- family:

  Optional GLM family; `NULL` fits
  [`lm()`](https://rdrr.io/r/stats/lm.html).

- subset, weights, na.action:

  Passed to
  [`lm()`](https://rdrr.io/r/stats/lm.html)/[`glm()`](https://rdrr.io/r/stats/glm.html)
  within each group (`subset` is applied before splitting). `weights`
  may name a column of `data`.

- pool:

  Pool the residual standard deviation across groups for inference.
  Defaults to `TRUE` for `lm` fits and `FALSE` for GLMs.

## Value

An object of class `mm_lmlist`: a named list of `lm`/`glm` fits (one per
group) with attributes `call`, `pool`, `groups`, and `failed`.

## Details

[`coef()`](https://rdrr.io/r/stats/coef.html) returns one row per group
(lme4's shape); [`confint()`](https://rdrr.io/r/stats/confint.html)
returns a `group x (lower, upper) x coefficient` array. With
`pool = TRUE` (the default for `lm` fits, as in lme4)
[`confint()`](https://rdrr.io/r/stats/confint.html) uses the residual
standard deviation pooled across groups, with the pooled residual
degrees of freedom, and [`sigma()`](https://rdrr.io/r/stats/sigma.html)
returns that pooled value; `confint(x, pool = FALSE)` gives each group's
own `lm` t interval (Wald z for GLMs), which is what
[`confint()`](https://rdrr.io/r/stats/confint.html) on an
[`lme4::lmList()`](https://rdrr.io/pkg/lme4/man/lmList.html) fit
returns. A group whose fit fails is kept as `NULL` and reported, never
silently dropped.

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:6, each = 8)), x = rep(0:7, 6))
d$y <- 2 + rnorm(6)[d$g] + (0.5 + rnorm(6, sd = 0.2)[d$g]) * d$x +
  rnorm(48)
ml <- mm_lmlist(y ~ x | g, d)
coef(ml)
#>   (Intercept)         x
#> 1   0.4814386 0.8733317
#> 2   2.8241058 0.4190855
#> 3   1.6432213 0.4689668
#> 4   3.6488567 0.5038056
#> 5   2.0853955 0.8825206
#> 6   1.5779637 0.5230055
confint(ml)[, , "x"]
#>       2.5 %    97.5 %
#> 1 0.5990784 1.1475851
#> 2 0.1448322 0.6933388
#> 3 0.1947135 0.7432202
#> 4 0.2295523 0.7780589
#> 5 0.6082673 1.1567739
#> 6 0.2487522 0.7972588
```
