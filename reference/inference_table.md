# Fixed-effect inference table

Fitted artifacts may carry Rust-owned fixed-effect inference rows. When
present, those rows are the source of truth for estimates, standard
errors, degrees of freedom, statistics, p-values, methods, status,
reliability, and unavailable reasons. Legacy objects without this
artifact field fall back to an unavailable table.

## Usage

``` r
inference_table(fit, ...)

# S3 method for class 'mm_lmm'
inference_table(
  fit,
  method = c("auto", "satterthwaite", "kenward_roger", "asymptotic", "none"),
  ...
)
```

## Arguments

- fit:

  A fitted `mm_lmm`.

- ...:

  Reserved for future methods.

- method:

  Inference method. `"auto"` (the default) returns the artifact-cached
  table that the engine resolved at fit time. Any other value
  (`"satterthwaite"`, `"kenward_roger"`, `"asymptotic"`, `"none"`)
  recomputes the table by dispatching one
  [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
  per fixed-effect term with the requested method, so refusals and
  reasons are surfaced honestly rather than silently swapped for the
  auto-resolved row.

## Value

An `mm_inference_table`.
