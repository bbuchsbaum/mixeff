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
and [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) all
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

Compiling returns a populated `mm_spec` with the JSON artifact attached.
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md),
[`random_options()`](https://bbuchsbaum.github.io/mixeff/reference/random_options.md),
and
[`compare_covariance()`](https://bbuchsbaum.github.io/mixeff/reference/compare_covariance.md)
render random-effects guidance from upstream random-term cards;
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fit
the model itself.

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
set.seed(1)
df <- data.frame(
  y       = rnorm(20),
  x       = rnorm(20),
  subject = factor(rep(letters[1:5], each = 4))
)
spec <- compile_model(y ~ x + (1 + x | subject), df)
audit(spec)
#> Audit Summary:
#>   overall [WARNING]: 6 warning(s); review attention lines before treating inference as routine
#>   attention [WARNING]: Model State / changes: Recommended:DesignTime:r0 -> 5 levels are below the v0 full-covariance threshold 15 for 3 covariance parameters
#>   attention [WARNING]: Random Effects / subject: group=subject, rows=20, levels=5, obs_per_level=4..4, basis=2, covariance=full, params=3, budget=too_rich; reason=5 levels are below the v0 full-covariance threshold 15 for 3 covariance parameters
#>   attention [WARNING]: Random-Effect Information Budget / subject: levels=5, rows=20, obs_per_level=4..4, basis=2, cov_params=3, levels/basis=2.50, levels/param=1.67, rows/param=6.67; total rows can be misleading for covariance support; risk=maximal covariance structure is too rich for the grouping-level budget; recommendation=full covariance asks for 3 parameter(s); options include using design_compiled with diagonal or reduced-rank covariance, treating the grouping factor as fixed when it is designed or directly compared, or collecting at least 15 grouping levels; explanation=20 rows are clustered into 5 grouping levels; covariance support is limited by grouping levels, not by total rows
#>   attention [WARNING]: Random Term Cards / subject: `subject` units differ in baseline and `x` slope; the model estimates whether these are associated. It uses 3 covariance parameters for intercept, x. The data contain 5 `subject` levels, with at least 4 rows per level (median 4). Within-group support is present for `x`. Within-group support was not assessed for `intercept`. The requested covariance is too rich for the available grouping-level information. Formula detail: `(1 + x | subject)`.
#>   attention [WARNING]: Policy Recommendations / subject: reduce_covariance: 5 levels are below the v0 full-covariance threshold 15 for 3 covariance parameters; lower-dimensional covariance=diagonal; inference=correlation parameters would not be interpreted as confirmatory
#>   attention [WARNING]: Diagnostics / covariance_too_rich: 5 levels are below the v0 full-covariance threshold 15 for 3 covariance parameters; affected=(1 + x | subject)
#> 
#> Requested Model:
#>   formula [INFO]: y ~ 1 + x + (1 + x | subject)
#>   model kind [INFO]: linear_mixed_model
#>   distribution/link [INFO]: gaussian/identity
#>   objective [INFO]: exact_gaussian
#>   convergence certificate [INFO]: exact_objective
#>   fixed terms [INFO]: 1, x
#>   random terms [INFO]: 1
#>   covariance parameter maps [INFO]: 1 map(s)
```
