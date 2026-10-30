# Reporting Linear Mixed Models

``` r

library(mixeff)
```

A mixed-model report should describe the observations and grouping
units, model specification, estimates, and inference methods. It should
also record convergence problems and quantities that could not be
computed.

`mixeff` stores this information in the fitted object.
[`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md)
extracts individual sections, and
[`model_report()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md)
assembles them with software-version information. This article shows how
to use those tables and turn the results into methods and results prose.

## Example model

The example uses repeated weekly observations from clinics, with clinic
as the grouping factor.

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

The analysis estimates average changes by week and treatment while
allowing clinics to have different baselines.

## Checking the specification

Before fitting, compile the model and read the design audit. This
separates formula interpretation from optimization and lets you check
the requested random-effects structure.

``` r

spec <- compile_model(score ~ week + treatment + (1 | clinic), clinic_visits)
audit(spec)
#> Audit Summary:
#>   overall [OK]: clean: no warnings or attention items
#>   attention [OK]: no warnings or unchecked inference-critical items
#> 
#> Requested Model:
#>   formula [INFO]: score ~ 1 + week + treatment + (1 | clinic)
#>   model kind [INFO]: linear_mixed_model
#>   distribution/link [INFO]: gaussian/identity
#>   objective [INFO]: exact_gaussian
#>   convergence certificate [INFO]: exact_objective
#>   fixed terms [INFO]: 1, week, treatment
#>   random terms [INFO]: 1
#>   covariance parameter maps [INFO]: 1 map(s)
```

## Fitting and summarizing the model

Fit the model and inspect its summary:

``` r

fit <- lmm(
  score ~ week + treatment + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)
```

The overview table is a compact first pass over the fitted model:
formula, fitting mode, number of observations, fit status, inference
availability, and versioned artifact information.

``` r

reporting_table(fit, "overview")
#>              field                                       value
#>        model_class                                         LMM
#>            formula     score ~ week + treatment + (1 | clinic)
#>  effective_formula score ~ 1 + week + treatment + (1 | clinic)
#>         fit_method                                        REML
#>               mode                   confirmatory_as_specified
#>               nobs                                          72
#>         fit_status                          converged_interior
#>          inference             3/3 available fixed-effect rows
#>    artifact_schema       mixedmodels.compiled_model_artifact 1
#>      crate_version                                  1.0.0-rc.5
#>    package_version                                       0.2.0
```

## Data and random-effects structure

The data-design section gives the grouping factors, their number of
levels, and rows per group.

``` r

reporting_table(fit, "data_design")
#>   group    role group_levels min_rows_per_group median_rows_per_group
#>  clinic unknown           12                  6                     6
#>  max_rows_per_group     status
#>                   6 sufficient
```

The random-term section records the grouping factor, basis, covariance
family, parameter count, and a description of each term.

``` r

reporting_table(fit, "random_terms")
#>  term_id original_fragment  group     basis covariance theta_parameters
#>       r0      (1 | clinic) clinic intercept     scalar                1
#>  design_status                                       english
#>     sufficient `clinic` units may differ in average outcome.
```

## Fixed effects and inference

[`summary()`](https://rdrr.io/r/base/summary.html) gives a familiar
coefficient table. Use it for console review.

``` r

coef_table <- summary(fit, method = "auto")$coef_table
knitr::kable(coef_table, digits = 4)
```

|  | Estimate | Std. Error | df | t value | Pr(\>\|t\|) | method |
|:---|---:|---:|---:|---:|---:|:---|
| (Intercept) | 7.6829 | 0.1965 | 12.5655 | 39.1065 | 0.0000 | satterthwaite |
| week | -0.2784 | 0.0260 | 59.0002 | -10.7280 | 0.0000 | satterthwaite |
| treatmentcoached | -0.8995 | 0.2623 | 9.9997 | -3.4298 | 0.0064 | satterthwaite |

For report assembly, use
[`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md).
The fixed-effect section keeps the estimate, uncertainty, statistic,
p-value, method, row status, and reliability label together.

``` r

reporting_table(fit, "fixed_effects")
#>              term   estimate  std_error        df  statistic statistic_name
#>       (Intercept)  7.6828778 0.19646018 12.565500  39.106539              t
#>              week -0.2783994 0.02595083 59.000150 -10.727955              t
#>  treatmentcoached -0.8994747 0.26225014  9.999671  -3.429835              t
#>       p_value        method    status reliability
#>  1.665335e-14 satterthwaite available    moderate
#>  1.776357e-15 satterthwaite available    moderate
#>  6.440574e-03 satterthwaite available    moderate
```

When you need to audit where those rows came from, request the audit
view.

``` r

fixed_audit <- reporting_table(fit, "fixed_effects", view = "audit")$table
fixed_audit[, c("term", "method", "status", "reliability", "source")]
#>               term        method    status reliability
#> 1      (Intercept) satterthwaite available    moderate
#> 2             week satterthwaite available    moderate
#> 3 treatmentcoached satterthwaite available    moderate
#>                         source
#> 1 fixed_effect_inference_table
#> 2 fixed_effect_inference_table
#> 3 fixed_effect_inference_table
```

## Turning the tables into methods and results prose

For this simulated example, the following prose uses the evaluated fit
and coefficient table. An applied report would also describe
recruitment, exclusions, and the study design; those facts cannot be
inferred from a model object.

> We fitted a Gaussian linear mixed model by restricted maximum
> likelihood to 72 observations from 12 clinics. The model included
> fixed effects of week and treatment and a random intercept for clinic.
> Fixed-effect tests used the satterthwaite method. The optimizer
> reported an interior solution (converged_interior).
>
> Adjusting for treatment, the estimated weekly change in score was
> -0.278 points (SE 0.026; p \< .001). This estimate describes the
> common weekly slope assumed by the model. The simulated example does
> not establish a treatment effect in a clinical population.

## Fitted variance components

The random-effect table reports fitted variance components, together
with their source and availability status.

``` r

reporting_table(fit, "random_effects")
#>     group   basis_lhs              kind  variance   std_dev    status
#>    clinic (Intercept)          variance 0.1827548 0.4274983 available
#>  Residual    Residual residual_variance 0.1414236 0.3760633 available
```

## Unavailable results and other qualifications

Build the full model report when you want the section map, software
provenance, and list of unavailable or inapplicable fields.

``` r

report <- model_report(fit)
report
#> mixeff model report
#>       field                                   value    status
#>     formula score ~ week + treatment + (1 | clinic) available
#>  fit_method                                    REML available
#>        nobs                                      72 available
#>  fit_status                      converged_interior available
#>   inference         3/3 available fixed-effect rows available
#> 
#> Sections:
#>   overview
#>   model_specification
#>   data_design
#>   random_terms
#>   random_effects
#>   fixed_effects
#>   fit_statistics
#>   optimizer
#>   comparison_ledger
#>   reproducibility
#>   unavailable
#> 
#> Unavailable/caveated fields: 1
```

The `unavailable` section identifies missing or inapplicable fields and
other qualifications, with reason codes and a record of the action
taken.

``` r

reporting_table(report, "unavailable")
#>            section             field         status
#>  comparison_ledger comparison_ledger not_applicable
#>                                    reason
#>  no_model_comparison_recorded_on_this_fit
```

## Assembling a report

Use the report sections as the source material for prose:

- `overview` records the formula, fit method, observations, fit status,
  and versioned artifact information.
- `data_design` records grouping-unit counts and rows per group.
- `random_terms` records the random-effects specification and its design
  support.
- `fixed_effects` records estimates and inference labels.
- `random_effects` records fitted variance components.
- `unavailable` records unavailable quantities and their reasons.

Use
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
for term tests, contrasts, and model comparisons. Use the introduction
vignette’s save-and-revive section
([`vignette("mixeff", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/mixeff.md))
for reports based on saved fits.
