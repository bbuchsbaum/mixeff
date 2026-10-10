# Caterpillar and Q-Q plots of random effects

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the result
of
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
draws a caterpillar plot (dot chart) of the conditional modes for each
grouping factor and coefficient, with levels sorted by their value and,
when the object carries conditional variances
(`ranef(fit, condVar = TRUE)`), `level` prediction intervals, mirroring
[`lattice::dotplot()`](https://rdrr.io/pkg/lattice/man/xyplot.html) on
`lme4`'s
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
output. [`qqnorm()`](https://rdrr.io/r/stats/qqnorm.html) draws normal
Q-Q plots of the conditional modes, one panel per grouping factor and
coefficient.

## Usage

``` r
# S3 method for class 'mm_ranef'
plot(x, level = 0.95, type = c("caterpillar", "qq"), ...)

# S3 method for class 'mm_ranef'
qqnorm(y, ...)
```

## Arguments

- x, y:

  An `mm_ranef` object from
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md).

- level:

  Coverage of the intervals (default 0.95).

- type:

  `"caterpillar"` (default) or `"qq"`.

- ...:

  Further graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x` (or `y`), invisibly.

## Details

When the lattice package is loaded,
[`lattice::dotplot()`](https://rdrr.io/pkg/lattice/man/xyplot.html) and
[`lattice::qqmath()`](https://rdrr.io/pkg/lattice/man/qqmath.html) also
work on
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
output and return one trellis object per grouping factor, as they do for
lme4.

Intervals are omitted for coefficients whose conditional variances are
not available (for example GLMM fits where `ranef(condVar = TRUE)`
returns `NA` variances).

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rnorm(60), x = rnorm(60),
  g = factor(rep(seq_len(10), each = 6))
)
fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
re <- ranef(fit, condVar = TRUE)
plot(re)

qqnorm(re)

```
