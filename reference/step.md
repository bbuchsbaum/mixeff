# Stepwise model selection

[`stats::step()`](https://rdrr.io/r/stats/step.html) is not a generic,
so mixeff (like lmerTest) exports a `step()` generic whose default
method is [`stats::step()`](https://rdrr.io/r/stats/step.html); the
`mm_lmm` method is the lmerTest-style backward elimination documented in
[`step.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/step.mm_lmm.md).

## Usage

``` r
step(object, ...)

# Default S3 method
step(object, ...)
```

## Arguments

- object:

  A fitted model.

- ...:

  Passed to the method (for the default, to
  [`stats::step()`](https://rdrr.io/r/stats/step.html)).

## Value

See the method.
