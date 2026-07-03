# Assess contrast estimability

Routes each requested contrast row through the Rust fixed-effect
inference bridge and reports the upstream estimability assessment
verbatim. Returned rows carry `status` (the closed enum from upstream:
`estimable`, `not_estimable`, `aliased`, ...), a boolean `estimable`
convenience flag, the contrast `rank` and `requested_rank`, and a stable
`reason` populated only when the engine refuses the contrast.

## Usage

``` r
estimability(fit, L = NULL, ...)

# S3 method for class 'mm_lmm'
estimability(fit, L = NULL, ...)
```

## Arguments

- fit:

  A fitted `mm_lmm`.

- L:

  Optional contrast vector or matrix. Defaults to the fixed-effect
  coefficient basis.

- ...:

  Reserved for future methods.

## Value

An `mm_estimability` object.
