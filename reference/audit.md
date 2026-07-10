# Audit a compiled model spec or fitted model

`audit()` returns the user-facing audit report attached to an `mm_spec`
(pre-fit) or an `mm_fit` (post-fit). The text is rendered by the
upstream Rust crate (the `mixedmodels.model_audit_report` schema's
`Display` impl) — Rust authors the wording, R formats nothing. Routing
every printed audit line through the upstream renderer is what enforces
the R9 "no advice creep" contract: drift in scope notes / tone is
visible in one place rather than scattered across R formatters.

## Usage

``` r
audit(object, ...)

# Default S3 method
audit(object, ...)

# S3 method for class 'mm_compiled'
audit(object, ...)
```

## Arguments

- object:

  An `mm_spec` produced by
  [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  or an `mm_fit` from
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md).

- ...:

  Reserved for future methods.

## Value

An object of class `mm_audit` carrying:

- `text`:

  the rendered report text (a single character string,
  newline-separated)

- `summary_text`:

  the compact report rendered by the upstream
  `ModelAuditReport::render_summary` (Audit Summary plus the Requested
  Model section)

- `design_audit`:

  the parsed `design_audit` field from the `CompiledModelArtifact`
  (random-term audits, fixed-effect rank, covariance kernel graph, ...)
  — `NULL` on uncompilable formulas

- `report`:

  the parsed upstream `ModelAuditReport` v2, including Rust-authored
  `random_term_cards` for downstream explanation verbs

- `random_term_cards`:

  the report's per-random-term cards, copied to the top level for
  convenient inspection

- `cross_card_constraints`:

  report-level constraints between random-term cards

- `diagnostics`:

  the parsed report diagnostics, falling back to artifact diagnostics
  when needed

`print.mm_audit` defaults to the compact upstream-rendered summary in
`summary_text`. Use `print(x, full = TRUE)` for the complete upstream
report stored in `text`.

## Details

`audit()` accepts an `mm_spec` from
[`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
or an `mm_fit` from
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
and emits the report sections (Requested Model, Model State,
Fixed/Random Effects, Information Budget, Dependence Paths,
Parameterization Trace, Effective Covariance, Policy Recommendations,
Optimizer, Inference, Diagnostics). Sections that depend on a fit
(Optimizer / Inference) report `not applicable before fitting` on a
pre-fit spec.

## Errors

Raises an `mm_schema_error` if the supplied object does not carry a
parsed artifact with the expected schema header.

## See also

[`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md).

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
  y       = rnorm(20),
  x       = rnorm(20),
  subject = factor(rep(letters[1:5], each = 4))
)
audit(compile_model(y ~ x + (1 + x | subject), df))
} # }
```
