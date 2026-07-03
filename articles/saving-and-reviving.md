# Saving and Reloading Fits

``` r

library(mixeff)
```

A mixed-model fit often outlives the R session that produced it. You fit
a model today, save it with the analysis, hand the script to a
collaborator, and reopen the fit six months later for a contrast, a
revision, or a referee response.

`mixeff` stores the fitted values, the random-effects design, the
convergence record, and the inference labels inside the R object — so
each of those tasks works after
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) without recomputing
the fit, and without depending on the original Rust handle.

## Fit a model

``` r

fit <- lmm(
  score ~ week + treatment + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)
```

Before saving, the ordinary extractors work as expected.

``` r

fixef(fit)
#>        (Intercept)               week treatment: coached 
#>          7.6828778         -0.2783994         -0.8994747
reporting_table(fit, "fixed_effects")
#>                 term   estimate  std_error  statistic statistic_name
#> 1        (Intercept)  7.6828778 0.19646018  39.106539              z
#> 2               week -0.2783994 0.02595083 -10.727955              z
#> 3 treatment: coached -0.8994747 0.26225014  -3.429835              z
#>        p_value            method    status reliability
#> 1 0.0000000000 asymptotic_wald_z available         low
#> 2 0.0000000000 asymptotic_wald_z available         low
#> 3 0.0006039485 asymptotic_wald_z available         low
```

## Round trip through RDS

``` r

path <- tempfile(fileext = ".rds")
saveRDS(fit, path)

restored <- readRDS(path)
restored <- revive(restored)
```

The restored object still answers the same fitted-model questions.

``` r

fixef(restored)
#>        (Intercept)               week treatment: coached 
#>          7.6828778         -0.2783994         -0.8994747
head(predict(restored))
#>        1        2        3        4        5        6 
#> 7.585932 7.307533 7.029134 6.750734 6.472335 6.193935
reporting_table(restored, "fixed_effects")
#>                 term   estimate  std_error  statistic statistic_name
#> 1        (Intercept)  7.6828778 0.19646018  39.106539              z
#> 2               week -0.2783994 0.02595083 -10.727955              z
#> 3 treatment: coached -0.8994747 0.26225014  -3.429835              z
#>        p_value            method    status reliability
#> 1 0.0000000000 asymptotic_wald_z available         low
#> 2 0.0000000000 asymptotic_wald_z available         low
#> 3 0.0006039485 asymptotic_wald_z available         low
```

## Rebuild design matrices when needed

Design extractors can be rebuilt from the stored formula and model
frame.

``` r

X <- model.matrix(restored, type = "fixed")
Z <- model.matrix(restored, type = "random")

dim(X)
#> [1] 72  3
dim(Z)
#> [1] 72 12
class(Z)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```

[`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md)
provides a small familiar subset for code that expects lme4-style names.

``` r

getME(restored, c("theta", "beta", "cnms"))
#> $theta
#> [1] 1.136772
#> 
#> $beta
#>        (Intercept)               week treatment: coached 
#>          7.6828778         -0.2783994         -0.8994747 
#> 
#> $cnms
#> $clinic
#> [1] "(Intercept)"
#> 
#> attr(,"class")
#> [1] "mm_cnms" "list"
```

## What stays explicit?

Quantities that the Rust inference contract cannot certify are marked
rather than fabricated. For full-rank fits,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the model-based
fixed-effect covariance from the upstream
`fixed_effect_covariance_matrix` payload. For rank-deficient or
otherwise uncertified fits, the matrix carries an
`mm_unavailable_reason` attribute and the values are `NA`.

``` r

V <- vcov(restored)
attr(V, "mm_unavailable_reason")
#> NULL
V
#>                     (Intercept)          week treatment: coached
#> (Intercept)         0.038596604 -1.683614e-03      -3.438757e-02
#> week               -0.001683614  6.734457e-04      -4.654064e-18
#> treatment: coached -0.034387568 -4.654064e-18       6.877514e-02
#> attr(,"mm_method")
#> [1] "model_based"
#> attr(,"mm_status")
#> [1] "available"
#> attr(,"mm_reliability")
#> [1] "high"
#> attr(,"mm_reason")
#> [1] NA
#> attr(,"mm_details")
#> attr(,"mm_details")$rank
#> [1] 3
#> 
#> attr(,"mm_details")$expected_rank
#> [1] 3
#> 
#> attr(,"mm_details")$aliased
#> list()
#> 
#> attr(,"mm_details")$matrix_rows
#> [1] 3
#> 
#> attr(,"mm_details")$matrix_cols
#> [1] 3
#> 
#> attr(,"mm_details")$finite
#> [1] TRUE
#> 
#> attr(,"mm_details")$symmetric
#> [1] TRUE
#> 
#> attr(,"mm_notes")
#> [1] "model-based fixed-effect covariance geometry; inference claims remain on fixed_effect_inference_table rows"
#> attr(,"mm_schema_name")
#> [1] "mixedmodels.fixed_effect_covariance_matrix"
#> attr(,"mm_schema_version")
#> [1] "1.0.0"
```

Conditional variances for random effects also survive the round trip.
With `condVar = TRUE`, each grouping table carries a finite `postVar`
array.

``` r

re <- ranef(restored, condVar = TRUE)
attr(re, "mm_unavailable_reason")
#> NULL
dim(attr(re$clinic, "postVar"))
#> [1]  1  1 12
```

## What should you report?

Use
[`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md)
for durable tables and
[`summary()`](https://rdrr.io/r/base/summary.html) for console output.
Both continue to work after a save/load cycle.

``` r

coef_table <- summary(restored, method = "auto")$coefficients
knitr::kable(coef_table, digits = 4)
```

|  | Estimate | Std. Error | df | t value | Pr(\>\|t\|) | method |
|:---|---:|---:|---:|---:|---:|:---|
| (Intercept) | 7.6829 | 0.1965 | 12.5650 | 39.1065 | 0.0000 | satterthwaite |
| week | -0.2784 | 0.0260 | 58.9997 | -10.7280 | 0.0000 | satterthwaite |
| treatment: coached | -0.8995 | 0.2623 | 9.9993 | -3.4298 | 0.0064 | satterthwaite |

``` r

reporting_table(restored, "overview")
#>                field                                       value
#> 1        model_class                                         LMM
#> 2            formula     score ~ week + treatment + (1 | clinic)
#> 3  effective_formula score ~ 1 + week + treatment + (1 | clinic)
#> 4         fit_method                                        REML
#> 5               mode                   confirmatory_as_specified
#> 6               nobs                                          72
#> 7         fit_status                          converged_interior
#> 8          inference             3/3 available fixed-effect rows
#> 9    artifact_schema       mixedmodels.compiled_model_artifact 1
#> 10     crate_version                                  1.0.0-rc.1
#> 11   package_version                                       0.1.0
```
