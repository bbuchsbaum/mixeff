# Compile a mixed-effects model spec without fitting

`compile_model()` parses the formula, runs the upstream semantic-IR /
design-audit pipeline against the supplied data, and returns an
`mm_spec` object — the audit-first analogue of the design-only step in
base [`lm()`](https://rdrr.io/r/stats/lm.html)'s
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) /
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) chain.
Nothing is optimized; nothing is fitted.
[`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md),
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md),
[`random_options()`](https://bbuchsbaum.github.io/mixeff/reference/random_options.md),
and (in Phase 1.E)
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) all
consume the same artifact.

## Usage

``` r
compile_model(formula, data)
```

## Arguments

- formula:

  A two-sided lme4-style formula, e.g. `y ~ x + (1 + x | subject)`.

- data:

  A `data.frame` whose columns include every variable named in
  `formula`. Variables with missing values raise an `mm_data_error`;
  pass `na.omit(data)` explicitly if that is what you want.

## Value

An object inheriting from `mm_spec` and containing:

- `call`:

  the matched call

- `formula`:

  the input formula

- `vars`:

  character vector of variables read from `data`

- `model_frame`:

  the data columns used to compile the artifact, retained so prefit
  audit views can evaluate nearby formula spellings

- `artifact`:

  parsed JSON artifact (the `mixedmodels.compiled_model_artifact` v1
  schema)

The raw artifact JSON is attached as `attr(spec$artifact, "raw_json")`
so the post-compile FFI calls (e.g., the internal `mm_audit_report_text`
primitive) can round-trip without re-encoding.

## Details

The compiled artifact is the structured truth: every print, summary, and
audit verb in mixeff reads back from it rather than re-deriving meaning
from formula text. R formats; Rust authors wording (PRD §9.6).

Phase 1 compile scope: returns a populated `mm_spec` with the JSON
artifact attached.
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md),
[`random_options()`](https://bbuchsbaum.github.io/mixeff/reference/random_options.md),
and
[`compare_covariance()`](https://bbuchsbaum.github.io/mixeff/reference/compare_covariance.md)
render random-effects guidance from upstream random-term cards; the fit
driver ([`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md))
lands in 1.E.

## Errors

Raises typed conditions (all inheriting from `mm_condition`):

- `mm_formula_error` — formula is not a two-sided R formula or fails
  parsing.

- `mm_data_error` — `data` is not a data.frame, refers to unknown
  variables, contains NAs in design columns, or has an unsupported
  column type.

- `mm_schema_error` — the artifact JSON returned by Rust does not match
  the wrapper's known schema set.

## See also

[`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md) for
the printed audit report.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
  y       = rnorm(20),
  x       = rnorm(20),
  subject = factor(rep(letters[1:5], each = 4))
)
spec <- compile_model(y ~ x + (1 + x | subject), df)
audit(spec)
} # }
```
