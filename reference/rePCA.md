# PCA of the random-effects covariance (lme4's `rePCA()`)

For each grouping factor, the singular values and left singular vectors
of the relative covariance factor (the block-diagonal Cholesky factor of
the random-effect covariance divided by the residual variance, lme4's
`Lambda`), returned as `prcomp` objects. Near-zero standard deviations
flag a random-effect structure that is over-parameterised (singular).
[`summary()`](https://rdrr.io/r/base/summary.html) gives each
component's proportion of variance.

## Usage

``` r
rePCA(x)

# Default S3 method
rePCA(x)

# S3 method for class 'mm_lmm'
rePCA(x)

# S3 method for class 'mm_glmm'
rePCA(x)
```

## Arguments

- x:

  A fitted `mm_lmm` or `mm_glmm` (for GLMM families without a dispersion
  parameter the covariance itself is used, as in lme4).

## Value

An object of class `c("mm_prcomplist", "prcomplist")`: a named list (one
entry per grouping factor) of `prcomp` objects.

## Details

Singular values match
[`lme4::rePCA()`](https://rdrr.io/pkg/lme4/man/rePCA.html); the rotation
columns are eigenvectors of the same matrix and may differ from lme4's
in sign.

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:10, each = 6)), x = rep(0:5, 10))
d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(10, sd = 0.3)[d$g] * d$x +
  rnorm(60)
fit <- lmm(y ~ x + (x | g), d, control = mm_control(verbose = -1))
summary(rePCA(fit))
#> $g
#> Importance of components:
#>                          [,1]    [,2]
#> Standard deviation     1.0550 0.29653
#> Proportion of Variance 0.9268 0.07322
#> Cumulative Proportion  0.9268 1.00000
#> 
```
