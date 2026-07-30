# Refit a mixeff LMM with a new response

`refit()` fits the same model formula to a new response by calling
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) with the
stored model frame and `REML` setting. Refitting is implemented for
linear mixed models only; calling `refit()` on a GLMM fit signals a
typed `mm_inference_unavailable` error.

## Usage

``` r
refit(object, newresp, ...)

# Default S3 method
refit(object, newresp, ...)

# S3 method for class 'mm_glmm'
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
