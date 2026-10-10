# Test whether a fit is singular or reduced-rank

Follows
[`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html): a
fit is singular when any covariance parameter with a lower bound of zero
(a diagonal element of a relative covariance factor,
`getME(x, "lower") == 0`) is below `tol`. Works for
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits.
[`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html) is
a plain function, so it cannot dispatch on mixeff fits; call
`is_singular()` instead. Fits without a stored theta map (very old
serialized objects) fall back to the engine's fit-status and
effective-covariance labels.

## Usage

``` r
is_singular(x, tol = 1e-04, ...)

# S3 method for class 'mm_lmm'
is_singular(x, tol = 1e-04, ...)

# S3 method for class 'mm_glmm'
is_singular(x, tol = 1e-04, ...)
```

## Arguments

- x:

  A fitted `mm_lmm` or `mm_glmm`.

- tol:

  Tolerance on the zero-bounded theta values (default `1e-4`, as in
  lme4).

- ...:

  Reserved for future methods.

## Value

A length-one logical value.
