# Test whether a mixeff fit has a live native handle

The native handle is a process-local cache. A `FALSE` result does not
mean the fit is unusable: Phase 2 extractors read from the durable
artifact and flat R-side payload, and
[`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
recreates the lazy cache after serialization.

## Usage

``` r
fit_handle_alive(fit, ...)

# S3 method for class 'mm_fit'
fit_handle_alive(fit, ...)
```

## Arguments

- fit:

  A fitted `mm_fit` object.

- ...:

  Reserved for future methods.

## Value

A length-one logical value.
