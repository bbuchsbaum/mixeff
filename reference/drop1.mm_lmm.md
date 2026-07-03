# Drop one fixed-effect term at a time

`drop1.mm_lmm()` refits reduced fixed-effect models and compares them to
the original fit. It is conservative: random-effect terms are preserved
exactly, and the reduced formulas are reported in the result table.

## Usage

``` r
# S3 method for class 'mm_lmm'
drop1(
  object,
  scope = NULL,
  test = c("none", "Chisq"),
  refit_for_comparison = c("auto", "error", "ml"),
  ...
)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- scope:

  Optional character vector of fixed-effect terms to drop.

- test:

  Comparison test label. `"Chisq"` reports asymptotic LRT rows; `"none"`
  reports information criteria only.

- refit_for_comparison:

  How to handle REML fits.

- ...:

  Reserved for future methods.

## Value

An `mm_drop1` object.
