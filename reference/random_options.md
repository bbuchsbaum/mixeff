# Inspect nearby random-effect spellings for one grouping factor

`random_options()` is an opt-in map over nearby random-effect structures
for a grouping factor. It recompiles each displayed spelling through the
same upstream audit path as
[`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md),
so support facts and block meanings come from Rust-authored
`RandomTermCard` records.

## Usage

``` r
random_options(spec, group, slope = NULL)
```

## Arguments

- spec:

  An `mm_spec` from
  [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  or, in later phases, an `mm_fit`.

- group:

  Grouping factor to inspect. May be supplied bare (`group = subject`)
  or as a string.

- slope:

  Optional slope variable to use for nearby slope-bearing spellings.
  When omitted, the function uses the first current random slope for
  `group`, then any scope-note fixed effect for `group`, then the first
  non-intercept fixed effect.

## Value

An object of class `mm_random_options` with an `options` data frame, the
upstream candidate `cards`, and the candidate audit reports.
