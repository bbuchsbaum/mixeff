# Confidence intervals for fixed effects of a mixeff GLMM

Asymptotic Wald intervals (`estimate +/- z * SE`) built from the Rust
fixed-effect inference table. Profile and bootstrap intervals are not
certified for GLMMs by the upstream contract and are refused with a
typed reason rather than approximated.

## Usage

``` r
# S3 method for class 'mm_glmm'
confint(
  object,
  parm,
  level = 0.95,
  method = c("asymptotic", "wald", "profile", "bootstrap"),
  ...
)
```

## Arguments

- object:

  A fitted `mm_glmm`.

- parm:

  Optional fixed-effect names or indices; defaults to all.

- level:

  Confidence level.

- method:

  `"asymptotic"` (the default; the package-wide name for the closed-form
  Wald interval) or its synonym `"wald"`.

- ...:

  Unused.

## Value

An `mm_confint` matrix of lower/upper bounds.
