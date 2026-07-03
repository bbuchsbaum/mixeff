# Degrees of freedom for a contrast

Degrees of freedom for a contrast

## Usage

``` r
df_for_contrast(
  fit,
  L,
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  ...
)

# S3 method for class 'mm_lmm'
df_for_contrast(
  fit,
  L,
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  ...
)
```

## Arguments

- fit:

  A fitted `mm_lmm`.

- L:

  A contrast vector or matrix.

- method:

  Requested degrees-of-freedom method.

- ...:

  Reserved for future methods.

## Value

A numeric vector of `NA` degrees of freedom with an
`mm_unavailable_reason` attribute.
