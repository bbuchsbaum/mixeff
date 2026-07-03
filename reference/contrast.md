# Contrast fixed effects

Note: this is not
[`emmeans::contrast`](https://rvlenth.github.io/emmeans/reference/contrast.html).
`contrast()` is mixeff's fixed-effect contrast front door. R validates
the contrast matrix shape, then asks Rust to evaluate estimability,
method prerequisites, standard errors, degrees of freedom, statistics,
p-values, reliability, and unavailable reasons.

## Usage

``` r
# S3 method for class 'mm_glmm'
contrast(fit, L, rhs = 0, method = c("asymptotic", "wald"), ...)

contrast(
  fit,
  L,
  rhs = 0,
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic",
    "boundary_lrt", "none"),
  bootstrap = NULL,
  ...
)

# S3 method for class 'mm_lmm'
contrast(
  fit,
  L,
  rhs = 0,
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic",
    "boundary_lrt", "none"),
  bootstrap = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted `mm_lmm`.

- L:

  A numeric contrast vector or matrix with one column per fixed effect.

- rhs:

  Numeric right-hand side, recycled to the number of contrasts.

- method:

  Requested inference method.

- ...:

  Reserved for future methods.

- bootstrap:

  Optional
  [`bootstrap_control()`](https://bbuchsbaum.github.io/mixeff/reference/bootstrap_control.md)
  object for `method = "bootstrap"`.

## Value

An `mm_contrast` object with a data-frame `table`. The `estimate` column
is the tested difference, `L beta_hat - rhs`.

## Details

For `mm_glmm` fits, contrasts use an asymptotic Wald z-test built from
the stored fixed-effect covariance (the GLMM contract does not provide
finite-sample df), so `method` accepts only `"asymptotic"` (alias
`"wald"`).
