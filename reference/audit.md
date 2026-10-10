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
set.seed(1)
df <- data.frame(
  y       = rnorm(20),
  x       = rnorm(20),
  subject = factor(rep(letters[1:5], each = 4))
)
audit(compile_model(y ~ x + (1 + x | subject), df))
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
