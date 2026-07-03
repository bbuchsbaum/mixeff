# Inspect inference methods available for this fit

`inference_options()` is the audit verb for fixed-effect inference. It
does not run any test; it predicts, from the fit's metadata, which
inference methods will succeed on this fit and at what approximate cost.
The goal is to remove trial-and-error: a user reading the table can see
which routes are immediately available, which will refuse and why, and
which require a bootstrap.

## Usage

``` r
inference_options(fit, term = NULL, nsim = 1000L, ...)

# S3 method for class 'mm_lmm'
inference_options(fit, term = NULL, nsim = 1000L, ...)
```

## Arguments

- fit:

  A fitted `mm_lmm`.

- term:

  Optional fixed-effect term name. Reserved for future per-term
  refinement; currently unused (the table is fit-level).

- nsim:

  Bootstrap replicate count to use when estimating cost. Used only to
  format the `approx_cost` column.

- ...:

  Reserved for future methods.

## Value

An `mm_inference_options` object with a `table` data frame of one row
per candidate method.

## Details

Like
[`random_options()`](https://bbuchsbaum.github.io/mixeff/reference/random_options.md),
this function does not rank or recommend. There is no "best method" row.
