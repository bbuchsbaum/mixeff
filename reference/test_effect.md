# Test a fixed-effect term

`test_effect()` asks Rust to construct fixed-effect term hypotheses and
returns the corresponding fixed-effect inference rows.

## Usage

``` r
test_effect(
  fit,
  term,
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "bootstrap_lrt",
    "cluster_bootstrap", "asymptotic", "boundary_lrt", "none"),
  bootstrap = NULL,
  group = NULL,
  ...
)

# S3 method for class 'mm_lmm'
test_effect(
  fit,
  term,
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "bootstrap_lrt",
    "cluster_bootstrap", "asymptotic", "boundary_lrt", "none"),
  bootstrap = NULL,
  group = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted `mm_lmm`.

- term:

  A fixed-effect term label.

- method:

  Requested inference method.

- bootstrap:

  Optional
  [`bootstrap_control()`](https://bbuchsbaum.github.io/mixeff/reference/bootstrap_control.md)
  object for bootstrap-backed methods.

- group:

  Optional grouping factor for `method = "cluster_bootstrap"`. Required
  for crossed or multi-grouping-factor models. In schema 1.0.0, cluster
  resampling is an estimator-distribution target and term-level p-values
  return `not_assessed` with a stable reason code.

- ...:

  Reserved for future methods.

## Value

An `mm_effect_test` object.
