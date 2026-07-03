# Update and re-fit a mixeff model

[`update()`](https://rdrr.io/r/stats/update.html) re-fits an
[`mm_lmm`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) or
[`mm_glmm`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) with
a modified model specification, mirroring
[`stats::update()`](https://rdrr.io/r/stats/update.html) and `lme4`'s
[`update()`](https://rdrr.io/r/stats/update.html) for the common cases:
changing the formula (`. ~ . - x`), toggling `REML`, swapping
`weights`/`offset`/`family`/`control`, or supplying new `data`.

## Usage

``` r
# S3 method for class 'mm_lmm'
update(object, formula., ..., evaluate = TRUE)

# S3 method for class 'mm_glmm'
update(object, formula., ..., evaluate = TRUE)
```

## Arguments

- object:

  A fitted `mm_lmm` or `mm_glmm`.

- formula.:

  A formula-change applied with
  [`stats::update.formula()`](https://rdrr.io/r/stats/update.formula.html);
  omit to keep the current formula. Random-effect terms (`(x | g)`,
  `(x || g)`) are preserved across `. ~ .` edits.

- ...:

  Arguments to override on the re-fit. For `mm_lmm`: `data`, `REML`,
  `weights`, `control`. For `mm_glmm`: additionally `family`, `offset`,
  `method`, `nAGQ`, `inference`.

- evaluate:

  If `TRUE` (default) re-fit and return the new model; if `FALSE` return
  the unevaluated call.

## Value

A new fitted model of the same class as `object`, or an unevaluated call
when `evaluate = FALSE`.

## Details

The re-fit reuses the fitted model frame
([`model.frame()`](https://rdrr.io/r/stats/model.frame.html)) as the
default data source, so formula edits that *remove* terms or change
estimation options work without re-supplying data. A formula edit that
introduces a **new** variable absent from the original model frame
requires an explicit `data =` argument.

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rnorm(80), x = rnorm(80), z = rnorm(80),
  g = factor(rep(seq_len(10), each = 8))
)
fit <- lmm(y ~ x + z + (1 | g), df, control = mm_control(verbose = -1))
# drop a fixed term
fit2 <- update(fit, . ~ . - z)
# refit by ML for a likelihood-ratio comparison
fit_ml <- update(fit, REML = FALSE)
fixef(fit2)
#> (Intercept)           x 
#>  0.07858531 -0.28479350 
```
