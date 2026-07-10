# Produce reporting tables for a fitted mixeff model

`model_report()` assembles a structured, publication-oriented report
from the Rust artifact fields carried by a fitted model plus R-owned
provenance such as the call and session metadata. `reporting_table()`
extracts one section as a data-frame-compatible object.

## Usage

``` r
model_report(fit, sections = "all", ...)

# S3 method for class 'mm_fit'
model_report(fit, sections = "all", ...)

reporting_table(object, section = "all", view = c("compact", "audit"), ...)

# S3 method for class 'mm_fit'
reporting_table(object, section = "all", view = c("compact", "audit"), ...)

# S3 method for class 'mm_model_comparison'
reporting_table(
  object,
  section = "comparison_ledger",
  view = c("compact", "audit"),
  ...
)

# S3 method for class 'mm_drop1'
reporting_table(
  object,
  section = "comparison_ledger",
  view = c("compact", "audit"),
  ...
)

# S3 method for class 'mm_random_effect_test'
reporting_table(object, section = "all", view = c("compact", "audit"), ...)
```

## Arguments

- fit:

  A fitted `mm_fit`, usually from
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md).

- sections:

  Character vector of report sections, or `"all"`.

- ...:

  Reserved for future methods.

- object:

  For `reporting_table()`: a fitted `mm_fit`, an `mm_model_report`, or a
  comparison/test object with a durable ledger.

- section:

  One section name, or `"all"`.

- view:

  `"compact"` for reader-facing columns, or `"audit"` for the full
  provenance table with `source`, `reason`, `details`, and related audit
  columns.

## Value

`model_report()` returns an `mm_model_report`. `reporting_table()`
returns an `mm_reporting_table` object: `$table` holds the section's
data frame (or `$sections` the named list when `section = "all"`).
