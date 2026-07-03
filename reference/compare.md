# Compare fitted mixeff models

`compare()` is the namespace-qualified model-comparison front door. For
LMMs it reports likelihood, information criteria, and asymptotic
likelihood-ratio comparisons. REML fits are refit by ML when
`refit_for_comparison = "auto"` or `"ml"`; `"error"` refuses that
comparison.

## Usage

``` r
compare(object, ...)

# S3 method for class 'mm_lmm'
compare(
  object,
  ...,
  target = c("fixed_effects", "random_effects", "prediction"),
  method = c("auto", "lrt", "bootstrap", "aic"),
  refit_for_comparison = c("auto", "error", "ml"),
  nsim = 0L,
  seed = NULL
)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- ...:

  Additional fitted `mm_lmm` objects.

- target:

  Comparison target label.

- method:

  `"auto"` / `"lrt"` for asymptotic likelihood-ratio rows, `"aic"` for
  information criteria only, or `"bootstrap"` for a small
  parametric-bootstrap LRT when `nsim > 0`.

- refit_for_comparison:

  How to handle REML fits.

- nsim:

  Number of bootstrap simulations for `method = "bootstrap"`.

- seed:

  Optional bootstrap seed.

## Value

An `mm_model_comparison` object with a data-frame `table`.
