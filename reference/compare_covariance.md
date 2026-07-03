# Compare covariance parameterizations for current random terms

`compare_covariance()` is a compact alternate view of the same upstream
random-term cards used by
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
and
[`random_options()`](https://bbuchsbaum.github.io/mixeff/reference/random_options.md).
For each current random-term card, it lays out the full, diagonal, and
scalar covariance families without ranking them.

## Usage

``` r
compare_covariance(spec)
```

## Arguments

- spec:

  An `mm_spec` from
  [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  or, in later phases, an `mm_fit`.

## Value

An object of class `mm_compare_covariance` with a `table` data frame and
the upstream cards it was derived from.
