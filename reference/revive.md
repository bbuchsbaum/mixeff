# Revive a serialized mixeff object

`revive()` restores the process-local parts of a `mixeff` object after
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) or a worker restart.
The fitted artifact and flat extractor values are the durable source of
truth; the Rust handle is only a cache and may be absent. Revival
recreates the lazy R-side cache, keeps a live native handle, and clears
a dead one (an external pointer restored by
[`readRDS()`](https://rdrr.io/r/base/readRDS.html)); it does not refit
the model. Without a live handle, computations that need the engine
model (contrasts, [`summary()`](https://rdrr.io/r/base/summary.html)
tests, predictions with intervals, profiles, bootstraps) refit it from
the stored model frame, with identical results.

## Usage

``` r
revive(fit, ...)

# S3 method for class 'mm_fit'
revive(fit, ...)
```

## Arguments

- fit:

  A fitted `mm_fit` object.

- ...:

  Reserved for future methods.

## Value

A revived `mm_fit` object.
