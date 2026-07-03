# Refit a mixeff LMM with a new response

`refit()` fits the same model formula to a new response by calling
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) with the
stored model frame and `REML` setting.

## Usage

``` r
refit(object, newresp, ...)

# S3 method for class 'mm_lmm'
refit(object, newresp, ...)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- newresp:

  Numeric response for `refit()`.

- ...:

  Reserved for future methods.

## Value

A new `mm_lmm`.
