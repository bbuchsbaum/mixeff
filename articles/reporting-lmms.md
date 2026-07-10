# Reporting Linear Mixed Models

``` r

library(mixeff)
```

A mixed-model report needs more than coefficient estimates. The reader
needs to know which observations and grouping units were used, what
random-effects structure was fitted, which inference method produced
each row, whether the optimizer reached an interior solution, and which
quantities were not available. Most of those are facts about the fit
that `lme4` does not store in the fitted object; the analyst has to
remember them, or reconstruct them from console output weeks after the
fact.

`mixeff` stores every one of them. This vignette follows a structured
mixed-model reporting checklist — describe the data and design, report
the model specification, show fixed and random effects, label inference
methods, preserve software provenance, and make caveats explicit instead
of hiding them in prose — drawing each section from the same fitted
object.

## What model will we report?

Use one dataset and one fitted model all the way through. Here, clinics
are the grouping units and each clinic contributes repeated weekly
observations.

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

## What did the formula request?

Before fitting, compile the model and read the design audit. This
separates formula interpretation from optimization and gives you a
stable place to check the requested random-effects structure.

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

## What fit was used?

Fit the model once. In reporting work, the fitted object is the source
for summary output, report sections, and provenance.

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
#>      crate_version                                  1.0.0-rc.1
#>    package_version                                       0.2.0
```

## Which design facts belong in the report?

A useful mixed-model report names the grouping units and their
information budget. The data-design section gives the number of levels
and rows per group.

``` r

reporting_table(fit, "data_design")
#>   group    role group_levels min_rows_per_group median_rows_per_group
#>  clinic unknown           12                  6                     6
#>  max_rows_per_group     status
#>                   6 sufficient
```

The random-term section translates the random-effects formula into
auditable rows. This is where the report records the grouping factor,
basis, covariance family, parameter count, and Rust-authored
plain-language description.

``` r

reporting_table(fit, "random_terms")
#>  term_id original_fragment  group     basis covariance theta_parameters
#>       r0      (1 | clinic) clinic intercept     scalar                1
#>  design_status                                       english
#>     sufficient `clinic` units may differ in average outcome.
```

## How are estimates and p-values labelled?

[`summary()`](https://rdrr.io/r/base/summary.html) gives a familiar
coefficient table. Use it for console review.

``` r

coef_table <- summary(fit, method = "auto")$coefficients
knitr::kable(coef_table, digits = 4)
```

|  | Estimate | Std. Error | df | t value | Pr(\>\|t\|) | method |
|:---|---:|---:|---:|---:|---:|:---|
| (Intercept) | 7.6829 | 0.1965 | 12.5650 | 39.1065 | 0.0000 | satterthwaite |
| week | -0.2784 | 0.0260 | 58.9997 | -10.7280 | 0.0000 | satterthwaite |
| treatmentcoached | -0.8995 | 0.2623 | 9.9993 | -3.4298 | 0.0064 | satterthwaite |

For report assembly, use
[`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md).
The fixed-effect section keeps the estimate, uncertainty, statistic,
p-value, method, row status, and reliability label together.

``` r

reporting_table(fit, "fixed_effects")
#>                term   estimate  std_error  statistic statistic_name
#>         (Intercept)  7.6828778 0.19646018  39.106539              z
#>                week -0.2783994 0.02595083 -10.727955              z
#>  treatment: coached -0.8994747 0.26225014  -3.429835              z
#>       p_value            method    status reliability
#>  0.0000000000 asymptotic_wald_z available         low
#>  0.0000000000 asymptotic_wald_z available         low
#>  0.0006039485 asymptotic_wald_z available         low
```

When you need to audit where those rows came from, request the audit
view.

``` r

fixed_audit <- reporting_table(fit, "fixed_effects", view = "audit")$table
fixed_audit[, c("term", "method", "status", "reliability", "source")]
#>                 term            method    status reliability
#> 1        (Intercept) asymptotic_wald_z available         low
#> 2               week asymptotic_wald_z available         low
#> 3 treatment: coached asymptotic_wald_z available         low
#>                         source
#> 1 fixed_effect_inference_table
#> 2 fixed_effect_inference_table
#> 3 fixed_effect_inference_table
```

## How are random effects reported?

The random-effect table reports variance components on the fitted scale.
In the current contract, those rows come from the Rust-authored
`mixedmodels.fit_summary` payload, so the source and availability status
travel with the report.

``` r

reporting_table(fit, "random_effects")
#>     group   basis_lhs              kind  variance   std_dev    status
#>    clinic (Intercept)          variance 0.1827548 0.4274983 available
#>  Residual    Residual residual_variance 0.1414236 0.3760633 available
```

## What is unavailable or caveated?

Build the full model report when you want the section map, software
provenance, and ledger of unavailable or not-applicable fields.

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

The unavailable ledger is part of the report, not an error condition. It
is where `mixeff` records schema gaps, not-applicable sections, and
other caveats with stable reasons and an action taken.

``` r

reporting_table(report, "unavailable")
#>            section             field         status
#>  comparison_ledger comparison_ledger not_applicable
#>                                    reason
#>  no_model_comparison_recorded_on_this_fit
```

## What should go into your written report?

Use the report sections as the source material for prose:

- `overview` records the formula, fit method, observations, fit status,
  and versioned artifact information.
- `data_design` records grouping-unit counts and rows per group.
- `random_terms` records the random-effects specification and its design
  support.
- `fixed_effects` records estimates and inference labels.
- `random_effects` records fitted variance components.
- `unavailable` records caveats that should not disappear from the
  analysis record.

The important habit is to report from the fitted object rather than from
memory. Use
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
for term tests, contrasts, and model comparisons. Use
[`vignette("saving-and-reviving", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/saving-and-reviving.md)
when the report needs to survive an RDS round trip.
