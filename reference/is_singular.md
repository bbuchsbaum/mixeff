# Test whether a fit is singular or reduced-rank

Test whether a fit is singular or reduced-rank

## Usage

``` r
is_singular(x, tol = 1e-04, ...)

# S3 method for class 'mm_lmm'
is_singular(x, tol = 1e-04, ...)
```

## Arguments

- x:

  A fitted `mm_lmm`.

- tol:

  Reserved for compatibility with lme4's `isSingular()`.

- ...:

  Reserved for future methods.

## Value

A length-one logical value.
