# Test whether a mixeff fit has a live native handle

The native handle is a process-local cache of the fitted engine model,
created by
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md). While
it is alive, computations that need the engine model (contrasts,
[`summary()`](https://rdrr.io/r/base/summary.html) tests, predictions
with intervals, `ranef(condVar = TRUE)`, profiles, bootstraps, model
comparison,
[`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md))
reuse it instead of refitting. A `FALSE` result does not mean the fit is
unusable: the handle does not survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) or a new R process,
and those computations then refit the model from the stored model frame,
with identical results; extractors read from the durable artifact and
flat R-side payload, and
[`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
recreates the lazy cache after serialization. Set
`options(mixeff.keep_handle = FALSE)` before fitting to not keep handles
(each one holds a copy of the model data for as long as the fit object
lives).

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
