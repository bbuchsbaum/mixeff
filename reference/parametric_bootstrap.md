# Parametric bootstrap likelihood-ratio comparison

Runs the engine-certified parametric-bootstrap likelihood-ratio test
between two nested ML-fitted LMMs through the Rust
`mm_bootstrap_lrt_json` entry point. The smaller model (fewer estimated
parameters) is the reduced model; the larger is the alternative. The
returned object carries the engine's replicate accounting (successful
and completed replicates, boundary count, Monte-Carlo standard error,
seed) rather than a bare [`mean()`](https://rdrr.io/r/base/mean.html)
p-value, so every reported number traces back to a versioned Rust
payload.

## Usage

``` r
parametric_bootstrap(null, alternative, nsim = 100L, seed = NULL, ...)
```

## Arguments

- null, alternative:

  Fitted `mm_lmm` objects. Order is irrelevant; the model with fewer
  parameters is treated as the reduced model.

- nsim:

  Number of bootstrap replicates.

- seed:

  Optional bootstrap seed.

- ...:

  Reserved for future methods.

## Value

An `mm_parametric_bootstrap` object.

## Details

The engine refuses REML fits: refit with `lmm(..., REML = FALSE)` before
calling. (`compare(method = "bootstrap")` refits REML to ML
automatically.)
