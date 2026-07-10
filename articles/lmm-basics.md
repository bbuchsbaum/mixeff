# Fitting Linear Mixed Models

``` r

library(mixeff)
```

You use a linear mixed model when observations are not independent:
visits within clinics, trials within people, students within classrooms.
The standard R answer is
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).
[`mixeff::lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
accepts the same formula language — `(1 | group)`, `(x | group)`,
`(1 + x || group)`, `(1 | a/b)` — and the fitted object answers to
[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and the rest of the
generics you already use.

What it adds, and the reason to use it for a walkthrough like this one,
is that the fitted object also carries the design, the convergence
status, and the inference labels with it. You can audit, summarize,
save, and reload without recomputing.

This vignette fits one clinic-visit model and reads the pieces an
analyst usually reaches for first: fixed effects, p-values, variance
components, fitted values, residuals, and a compact design report.

## What data are we fitting?

``` r

head(clinic_visits)
#>      score week treatment clinic
#> 1 7.330830    0     usual      1
#> 2 7.439355    1     usual      1
#> 3 7.210269    2     usual      1
#> 4 7.199514    3     usual      1
#> 5 6.597443    4     usual      1
#> 6 5.487172    5     usual      1
```

`score` is the response, `week` is a numeric predictor, `treatment` is a
fixed-effect factor, and `clinic` identifies repeated visits from the
same clinic.

## What happens when you call `lmm()`?

The model estimates average week and treatment effects while allowing
clinics to have different baselines.

``` r

fit <- lmm(
  score ~ week + treatment + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)
```

Printing the fit gives the formula, convergence status, likelihood
summary, residual scale, and fixed effects.

``` r

fit
#> Linear mixed model fit by REML
#> Formula: score ~ week + treatment + (1 | clinic)
#> Fit status: converged_interior
#> Optimizer: pattern_search; iterations: 23; objective: 95.0585
#> nobs: 72, sigma: 0.376063, logLik: -47.5293
#> Fixed effects:
#>      (Intercept)             week treatmentcoached 
#>         7.682880        -0.278399        -0.899475 
#> Audit verbs: audit(), diagnostics(), inference_table(), model_report()
```

## How do you read the coefficient table?

[`summary()`](https://rdrr.io/r/base/summary.html) gives estimates,
standard errors, degrees of freedom, test statistics, p-values, and
method labels when those inference rows are available.

``` r

coef_table <- summary(fit, method = "auto")$coefficients
knitr::kable(coef_table, digits = 4)
```

|  | Estimate | Std. Error | df | t value | Pr(\>\|t\|) | method |
|:---|---:|---:|---:|---:|---:|:---|
| (Intercept) | 7.6829 | 0.1965 | 12.5650 | 39.1065 | 0.0000 | satterthwaite |
| week | -0.2784 | 0.0260 | 58.9997 | -10.7280 | 0.0000 | satterthwaite |
| treatmentcoached | -0.8995 | 0.2623 | 9.9993 | -3.4298 | 0.0064 | satterthwaite |

For a focused term-level test, use
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md).

``` r

test_effect(fit, "treatment", method = "kenward_roger")
#> Effect tests:
#>       term den_df statistic statistic_name     p_value        method    status
#>  treatment     10 -3.429835              t 0.006440269 kenward_roger available
#> Full audit columns available in `x$table` (9 hidden).
```

## Which familiar extractors work?

The usual fixed-effect and fit-statistic extractors are available.

``` r

fixef(fit)
#>      (Intercept)             week treatmentcoached 
#>        7.6828778       -0.2783994       -0.8994747
sigma(fit)
#> [1] 0.3760633
logLik(fit)
#> 'log Lik.' -47.52925 (df=5)
```

[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
reports fitted variance components, and
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
returns conditional random effects by grouping factor.

``` r

VarCorr(fit)
#> Variance components:
#>   group        name variance  std_dev correlation
#>  clinic (Intercept) 0.182755 0.427498            
#> Residual std. dev.: 0.376063
head(ranef(fit)$clinic)
#>   (Intercept)
#> 1 -0.09694543
#> 2 -0.45868391
#> 3 -0.32452643
#> 4  0.14989203
#> 5 -0.48550430
#> 6 -0.17802563
```

## How do prediction and residuals line up?

For fitted data, [`predict()`](https://rdrr.io/r/stats/predict.html)
returns in-sample fitted values. Use `re.form = NA` for the fixed-effect
part only.

``` r

prediction_check <- data.frame(
  score = clinic_visits$score,
  fitted = predict(fit),
  fixed_only = predict(fit, re.form = NA),
  residual = residuals(fit)
)
head(prediction_check)
#>      score   fitted fixed_only   residual
#> 1 7.330830 7.585932   7.682878 -0.2551023
#> 2 7.439355 7.307533   7.404478  0.1318223
#> 3 7.210269 7.029134   7.126079  0.1811353
#> 4 7.199514 6.750734   6.847680  0.4487797
#> 5 6.597443 6.472335   6.569280  0.1251080
#> 6 5.487172 6.193935   6.290881 -0.7067635
```

## Where is the design audit?

Use reporting tables when you want a compact, data-frame result for a
report or review. The data-design table is often the first one to
inspect.

``` r

reporting_table(fit, "data_design")
#>   group    role group_levels min_rows_per_group median_rows_per_group
#>  clinic unknown           12                  6                     6
#>  max_rows_per_group     status
#>                   6 sufficient
```

The random-term table translates the random-effects part of the formula
into rows.

``` r

reporting_table(fit, "random_terms")
#>  term_id original_fragment  group     basis covariance theta_parameters
#>       r0      (1 | clinic) clinic intercept     scalar                1
#>  design_status                                       english
#>     sufficient `clinic` units may differ in average outcome.
```

For lower-level checks, use
[`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md),
[`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md),
and
[`parameterization()`](https://bbuchsbaum.github.io/mixeff/reference/parameterization.md).

``` r

fit_status(fit)
#> [1] "converged_interior"
diagnostics(fit)
#> Diagnostics:
#>        code severity        stage affected_terms
#>  scope_note     info design_audit             r0
#> 
#> Messages:
#>   scope_note: `week` varies within `clinic`, so a `clinic`-level slope is structurally
#>         possible
parameterization(fit)
#> Covariance parameterization:
#>  term_id  group source_syntax covariance_family                   theta_name
#>       r0 clinic  (1 | clinic)            scalar theta[0:intercept,intercept]
#>  theta_value theta_status                        varcorr_entries
#>     1.136772         free standard_deviation[intercept]=0.427498
#> Full theta/Lambda columns available in `x$table` (9 hidden).
```

## What should you read next?

Use
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
for p-values, contrasts, term tests, and model comparisons. Use
[`vignette("demystifying-formulas", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/demystifying-formulas.md)
when you want to understand how `(1 | clinic)`, `(week | clinic)`, split
blocks, and `||` change the random-effects structure.
