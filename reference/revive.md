# Revive a serialized mixeff object

`revive()` restores the process-local parts of a `mixeff` object after
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) or a worker restart.
The fitted artifact and flat extractor values are the durable source of
truth; the Rust handle is only a cache and may be absent. In the current
bridge, revival recreates the lazy R-side cache and explicitly leaves
`rust_handle = NULL`.

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
