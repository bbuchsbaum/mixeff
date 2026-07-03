# Parse and canonicalize an lme4-style formula

`mm_parse_formula()` parses a formula string through the Rust formula
parser and returns its canonical `Display` rendering. This is the Phase
0 round-trip primitive: equivalent formula spellings produce identical
canonical strings, so equivalence-class testing in R is just string
comparison on the canonical form.

## Usage

``` r
mm_parse_formula(formula)
```

## Arguments

- formula:

  A single character string (length 1, non-NA, non-empty), or a
  one-sided / two-sided R `formula` object. R `formula` objects are
  coerced to character via
  [`format()`](https://rdrr.io/r/base/format.html) before parsing.

## Value

A single character string: the canonical rendering of the parsed
formula.

## Errors

Parse failures are signalled as a typed `mm_formula_error` condition
(also inheriting from `mm_condition` and `error`). The condition object
carries the original input string in its `formula` field. Catch with
`tryCatch(..., mm_formula_error = handler)`.

## Examples

``` r
mm_parse_formula("y ~ x + (1 | g)")
#> [1] "y ~ 1 + x + (1 | g)"
mm_parse_formula(y ~ x + (1 | g))
#> [1] "y ~ 1 + x + (1 | g)"
```
