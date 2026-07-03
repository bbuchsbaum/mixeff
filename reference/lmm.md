# Fit a linear mixed-effects model

`lmm()` is mixeff's Phase 1 linear mixed-model fit driver. It compiles
the requested lme4-style formula, emits the same
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
view that pre-fit audit users see as a message (silence it with
[`suppressMessages()`](https://rdrr.io/r/base/message.html) or
`mm_control(verbose = -1)`), then delegates the numerical fit to the
upstream Rust `LinearMixedModel`.

## Usage

``` r
lmm(
  formula,
  data,
  REML = TRUE,
  weights = NULL,
  subset = NULL,
  na.action = NULL,
  contrasts = NULL,
  control = mm_control()
)
```

## Arguments

- formula:

  A two-sided lme4-style formula, e.g. `y ~ x + (1 + x | subject)`.

- data:

  A `data.frame` containing all variables in `formula`.

- REML:

  Logical; fit by restricted maximum likelihood when `TRUE`.

- weights:

  Optional positive numeric case weights, either a vector with one value
  per row or an expression evaluated in `data`.

- subset:

  Optional expression selecting rows of `data`, evaluated in `data` (as
  in [`stats::lm()`](https://rdrr.io/r/stats/lm.html)).

- na.action:

  Optional function controlling missing-value handling, applied to the
  model variables before fitting (e.g.
  [stats::na.omit](https://rdrr.io/r/stats/na.fail.html)). The default
  (`NULL`) refuses any `NA` in a model variable with a typed
  `mm_data_error` (audit-first: missing-data dropping must be opt-in).
  Pass `na.action = na.omit` for lme4's complete-case behaviour.

- contrasts:

  Optional named list of factor contrasts. The engine codes all factors
  with treatment contrasts; a request for any other coding is refused
  (recode the factor instead).

- control:

  A list from
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md).

## Value

An object of class `mm_lmm`, also inheriting from `mm_fit` and
`mm_compiled`.

## Details

The returned object is deliberately serializable: fixed effects, theta,
sigma, likelihood summaries, fitted values, residuals, random effects,
and the post-fit compiler artifact are all stored directly on the R
object. The native Rust handle is treated as a rebuildable cache, not as
the source of truth.

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rnorm(80),
  x = rnorm(80),
  subject = factor(rep(seq_len(20), each = 4))
)
fit <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))
fixef(fit)
#> (Intercept)           x 
#>  0.07858531 -0.28479350 
VarCorr(fit)
#> Variance components:
#>    group        name variance std_dev correlation       note
#>  subject (Intercept)        0       0             [boundary]
#> [boundary]: variance component is at the boundary of the parameter space.
#> Residual std. dev.: 0.867021
summary(fit)
#> Linear mixed model fit by REML
#> Formula: y ~ x + (1 | subject)
#> Fit status: converged_reduced_rank
#> 
#> Variance components:
#>    group        name variance std_dev correlation       note
#>  subject (Intercept)        0       0             [boundary]
#> [boundary]: variance component is at the boundary of the parameter space.
#> Residual std. dev.: 0.867021
#> 
#> Fixed effects:
#>                Estimate Std. Error   z value Pr(>|z|)            method
#> (Intercept)  0.07858531  0.0974727  0.806229 0.420111 asymptotic_wald_z
#> x           -0.28479350  0.1055607 -2.697913 0.006978 asymptotic_wald_z
#> 
#> Inference status:
#>         term            method    status reliability         reliability_reason
#>  (Intercept) asymptotic_wald_z available         low asymptotic_wald_z_fallback
#>            x asymptotic_wald_z available         low asymptotic_wald_z_fallback
#> 
#> Notes:
#>   asymptotic Wald z is a labeled fallback, not a finite-sample correction
```
