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

An `mm_df_for_contrast` object with `$table` (one row per contrast:
`contrast`, `df`, `method`, `requested_method`, `reason`), `$df` (the
named numeric vector), and `$method`. When the method is `"none"` or the
engine refuses, `df` is `NA` and `reason` records why.
