# Fit a linear mixed-effects model

`lmm()` is mixeff's linear mixed-model fit driver. It compiles the
requested lme4-style formula, emits the same
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
  na.action = getOption("na.action"),
  contrasts = NULL,
  control = mm_control(),
  offset = NULL
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

  A function (or function name) controlling missing-value handling, as
  in [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) and
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html). The default is
  `getOption("na.action")` (`na.omit` unless changed). It is applied
  after `subset` to the variables the model uses only – the formula
  variables (fixed and random parts, including evaluated transforms such
  as `log(x)`) plus `weights` and `offset` – so an `NA` in an unused
  column of `data` never drops a row, exactly like
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html). Dropped
  rows are never silent: a typed `mm_rows_dropped` message reports how
  many rows were dropped and which variables had missing values (silence
  it with `mm_control(verbose = -1)`). The dropped rows are recorded as
  `na.action(fit)` and `attr(model.frame(fit), "na.action")`, and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) counts the rows used.
  [stats::na.omit](https://rdrr.io/r/stats/na.fail.html) keeps the
  compact shape;
  [stats::na.exclude](https://rdrr.io/r/stats/na.fail.html) pads
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and in-sample
  [`predict()`](https://rdrr.io/r/stats/predict.html) back to the rows
  of `data` with `NA`;
  [stats::na.fail](https://rdrr.io/r/stats/na.fail.html) refuses missing
  values with a typed `mm_data_error`;
  [stats::na.pass](https://rdrr.io/r/stats/na.fail.html) (or no action
  configured) is refused with a typed `mm_data_error`, since a mixed
  model cannot be fitted to missing values. Factor levels left unused
  after `subset`/NA removal are dropped, as
  `model.frame(drop.unused.levels = TRUE)` does in lme4.

- contrasts:

  Optional named list of factor contrasts. The engine codes unordered
  factors with treatment contrasts (`contr.treatment`) and ordered
  factors with orthonormal polynomial contrasts (`contr.poly`), matching
  lme4/R defaults. A request for any other coding is refused (recode the
  factor instead).

- control:

  A list from
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md).

- offset:

  Optional known component of the linear predictor: a numeric vector
  with one value per row of `data`, or an expression evaluated in `data`
  (as in `lme4::lmer(offset = )`).
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the formula
  are supported too and are added to it. The engine fits `y - offset`;
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) add the offset
  back, so they match lmer.

## Value

An object of class `mm_lmm`, also inheriting from `mm_fit` and
`mm_compiled`.

## Details

The returned object is deliberately serializable: fixed effects, theta,
sigma, likelihood summaries, fitted values, residuals, random effects,
and the post-fit compiler artifact are all stored directly on the R
object. The native Rust handle is treated as a rebuildable cache, not as
the source of truth.

Optimization runs inside a single native call with no progress output:
the pre-fit explanation block (when `verbose >= 0`) is the last thing
printed before the fitted result returns. The fit checks for a user
interrupt (Ctrl-C / Esc) between optimizer evaluations and stops with a
typed `mm_interrupted` error (the check happens at evaluation
boundaries, so a single very expensive evaluation finishes first).
Evaluation budgets are bounded (a bounded budget caps optimizer
iterations; it does not prove every native evaluation terminates);
runtime on large problems is governed by `mm_control(max_feval = )`.

## Formula language

The fixed-effects part accepts R's formula language as
[`lm()`](https://rdrr.io/r/stats/lm.html)/`lmer()` do:
[`factor()`](https://rdrr.io/r/base/factor.html),
[`relevel()`](https://rdrr.io/r/stats/relevel.html),
[`cut()`](https://rdrr.io/r/base/cut.html),
[`scale()`](https://rdrr.io/r/base/scale.html),
[`poly()`](https://rdrr.io/r/stats/poly.html),
[`splines::ns()`](https://rdrr.io/r/splines/ns.html)/`bs()`,
[`as.numeric()`](https://rdrr.io/r/base/numeric.html), comparisons such
as `I(x > 0)`, [`offset()`](https://rdrr.io/r/stats/offset.html), and
the operators `*`, `:`, `^`, `%in%`, `-` and parenthesised groups such
as `(a + b)^2` or `a * (b + c)`. Terms the engine cannot evaluate itself
are evaluated in R (with
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
semantics: on the full `data`, before `subset`/`na.action`) and sent to
the engine as synthetic columns; coefficient names match lme4's
(`poly(x, 2)1`, `factor(cyl)6`), and
[`predict()`](https://rdrr.io/r/stats/predict.html) on new data reuses
the training basis (`predvars`), as
[`lm()`](https://rdrr.io/r/stats/lm.html) does. A logical predictor is
coded like a factor with levels `FALSE`/`TRUE`, as
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) codes it.
Such transforms are not supported inside random-effect terms, where the
engine refuses them with a typed error.

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
#>  Groups   Name        Std.Dev.           
#>  subject  (Intercept) 0.00000  [boundary]
#>  Residual             0.86702            
#> [boundary]: variance component is at the boundary of the parameter space.
summary(fit)
#> Linear mixed model fit by REML
#> Formula: y ~ x + (1 | subject)
#> Fit status: converged_reduced_rank
#> 
#>  Groups   Name        Std.Dev.           
#>  subject  (Intercept) 0.00000  [boundary]
#>  Residual             0.86702            
#> [boundary]: variance component is at the boundary of the parameter space.
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
