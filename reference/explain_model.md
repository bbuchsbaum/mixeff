# Explain the random-effects structure of a compiled model

`explain_model()` renders the random-effects guidance surface for an
`mm_spec` returned by
[`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
or, in later phases, an `mm_fit`. It formats the upstream
`RandomTermCard` and diagnostic payloads; Rust remains the source of
truth for per-block English wording and design facts.

## Usage

``` r
explain_model(spec)
```

## Arguments

- spec:

  An `mm_spec` produced by
  [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  or an `mm_fit`.

## Value

An object of class `mm_explanation` carrying:

- `text`:

  the rendered explanation text

- `cards`:

  the upstream random-term cards

- `cross_card_constraints`:

  report-level constraints between cards

- `diagnostics`:

  the upstream diagnostics used for design notes

- `report`:

  the parsed upstream `ModelAuditReport`

## Errors

Raises an `mm_schema_error` if `spec` is not an `mm_spec`/`mm_fit` or
does not carry a valid compiled artifact.

## See also

[`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md),
[`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
  y = rnorm(20),
  t = rep(0:3, 5),
  s = factor(rep(1:5, each = 4))
)
explain_model(compile_model(y ~ t + (1 | s), df))
} # }
```
