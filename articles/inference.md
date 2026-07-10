# Inference You Can Report

``` r

library(mixeff)
```

Mixed-model users want p-values. They also want to know what those
p-values mean: which method produced them, and how trustworthy that
method is for this particular fit. `lme4` reports the number; `mixeff`
reports the number alongside its provenance.

When a coefficient, contrast, or term test has an available method, the
p-value is printed with the method name. When the requested method is
unavailable on this fit, the row says so, with a stable
[reason](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.md)
rather than an apologetic substitute.

## What model are we fitting?

``` r

fit <- lmm(
  score ~ week + treatment + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)
```

The fixed effects ask whether scores change across weeks and whether the
coached program differs from usual care, while clinic baselines are
allowed to vary.

## Coefficient p-values

Start with the ordinary coefficient table.

``` r

coef_table <- summary(fit, tests = "coefficients", method = "auto")$coefficients
knitr::kable(coef_table, digits = 4)
```

|  | Estimate | Std. Error | df | t value | Pr(\>\|t\|) | method |
|:---|---:|---:|---:|---:|---:|:---|
| (Intercept) | 7.6829 | 0.1965 | 12.5650 | 39.1065 | 0.0000 | satterthwaite |
| week | -0.2784 | 0.0260 | 58.9997 | -10.7280 | 0.0000 | satterthwaite |
| treatmentcoached | -0.8995 | 0.2623 | 9.9993 | -3.4298 | 0.0064 | satterthwaite |

The last column tells you the method used for each available p-value.

``` r

inference_table(fit)
#> Inference table:
#>                term            label        kind   estimate  std_error df
#>         (Intercept)      (Intercept) coefficient  7.6828778 0.19646018 NA
#>                week             week coefficient -0.2783994 0.02595083 NA
#>  treatment: coached treatmentcoached coefficient -0.8994747 0.26225014 NA
#>  numerator_df denominator_df  statistic statistic_name      p_value
#>            NA             NA  39.106539              z 0.0000000000
#>            NA             NA -10.727955              z 0.0000000000
#>            NA             NA  -3.429835              z 0.0006039485
#>             method    status reliability         reliability_reason reason
#>  asymptotic_wald_z available         low asymptotic_wald_z_fallback   <NA>
#>  asymptotic_wald_z available         low asymptotic_wald_z_fallback   <NA>
#>  asymptotic_wald_z available         low asymptotic_wald_z_fallback   <NA>
#>  reason_code reason_detail estimability details        notes
#>         <NA>          <NA> fixed_co....         asymptot....
#>         <NA>          <NA> fixed_co....         asymptot....
#>         <NA>          <NA> fixed_co....         asymptot....
```

## Contrasts

[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
is the direct route when you know the fixed-effect comparison you want.
This contrast asks whether the coached program differs from usual care.

``` r

L <- c(0, 0, 1)
names(L) <- names(fixef(fit))

contrast(fit, L, method = "satterthwaite")
#> Fixed-effect contrasts:
#>  contrast   estimate rhs std_error       df statistic statistic_name
#>        c1 -0.8994747   0 0.2622501 9.999273 -3.429835              t
#>      p_value        method requested_method    status reliability
#>  0.006440943 satterthwaite    satterthwaite available    moderate
#>                             reliability_reason estimability reason reason_code
#>  satterthwaite_finite_difference_approximation fixed_co....   <NA>        <NA>
#>  reason_detail      details        notes
#>           <NA> list(fam.... Satterth....
```

## Term tests

Use
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md)
for a named fixed-effect term.

``` r

test_effect(fit, "treatment", method = "kenward_roger")
#> Effect tests:
#>       term den_df statistic statistic_name     p_value        method    status
#>  treatment     10 -3.429835              t 0.006440269 kenward_roger available
#> Full audit columns available in `x$table` (9 hidden).
```

Single-model [`anova()`](https://rdrr.io/r/stats/anova.html) gives the
same kind of term-level table.

``` r

anova(fit, type = "III", method = "kenward_roger")
#> Type III analysis of fixed effects (method: kenward_roger):
#>       term num_df den_df statistic      p_value        method
#>       week      1     59 115.08903 1.776357e-15 kenward_roger
#>  treatment      1     10  11.76377 6.440269e-03 kenward_roger
#> Full provenance columns available in `$table` (type, statistic_name, requested_method, status, reliability, reason, details, notes).
```

## Model comparisons

For nested fixed-effect comparisons, fit the reduced model and compare
it with the full model.

``` r

reduced <- lmm(
  score ~ week + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)

compare(reduced, fit)
#> Model comparison:
#>  model                                     formula nobs df    logLik deviance
#>     m1             score ~ 1 + week + (1 | clinic)   72  4 -47.83234 95.66467
#>     m2 score ~ 1 + week + treatment + (1 | clinic)   72  5 -43.16629 86.33259
#>        AIC      BIC delta_aic delta_bic  REML refit         fit_status delta_df
#>  103.66467 112.7713  7.332086   5.05542 FALSE  TRUE converged_interior       NA
#>   96.33259 107.7159  0.000000   0.00000 FALSE  TRUE converged_interior        1
#>       LRT     p_value         method          status reason reason_code
#>        NA          NA asymptotic_lrt reference_model               <NA>
#>  9.332086 0.002251758 asymptotic_lrt       available               <NA>
#>      comparison_class lrt_available information_criteria_available
#>                  <NA>         FALSE                           TRUE
#>  nested_fixed_effects          TRUE                           TRUE
#>  requires_ml_refit loglik_within_optimizer_tol rust_method rust_refit_policy
#>              FALSE                       FALSE        auto             never
#>              FALSE                       FALSE        auto             never
```

[`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)
records that likelihood-ratio p-values are asymptotic. If you want a
simulation-based check for a small example, use the bootstrap path.

``` r

compare(reduced, fit, method = "bootstrap", nsim = 10, seed = 7)
#> Model comparison:
#>  model                                     formula nobs df    logLik deviance
#>     m1             score ~ 1 + week + (1 | clinic)   72  4 -47.83234 95.66467
#>     m2 score ~ 1 + week + treatment + (1 | clinic)   72  5 -43.16629 86.33259
#>        AIC      BIC delta_aic delta_bic  REML refit         fit_status delta_df
#>  103.66467 112.7713  7.332086   5.05542 FALSE  TRUE converged_interior       NA
#>   96.33259 107.7159  0.000000   0.00000 FALSE  TRUE converged_interior        1
#>       LRT p_value                   method          status
#>        NA      NA           asymptotic_lrt reference_model
#>  9.332086       0 parametric_bootstrap_lrt       available
#>                                               reason reason_code
#>                                                             <NA>
#>  parametric bootstrap LRT (10/10 replicates, MCSE=0)        <NA>
#>      comparison_class lrt_available information_criteria_available
#>                  <NA>         FALSE                           TRUE
#>  nested_fixed_effects          TRUE                           TRUE
#>  requires_ml_refit loglik_within_optimizer_tol rust_method rust_refit_policy
#>              FALSE                       FALSE        auto             never
#>              FALSE                       FALSE        auto             never
```

## Unavailable is still useful information

Population-level prediction standard errors and intervals are available
via `re.form = NA` (the Wald SE of the fixed-effect linear predictor):

``` r

pop <- predict(fit, re.form = NA, se.fit = TRUE)
head(pop$se.fit)
#>         1         2         3         4         5         6 
#> 0.1964602 0.1894804 0.1858923 0.1858923 0.1894804 0.1964602
head(predict(fit, re.form = NA, interval = "confidence"))
#>        fit      lwr      upr
#> 1 7.682878 7.297823 8.067933
#> 2 7.404478 7.033104 7.775853
#> 3 7.126079 6.761737 7.490421
#> 4 6.847680 6.483338 7.212022
#> 5 6.569280 6.197906 6.940655
#> 6 6.290881 5.905826 6.675936
```

*Conditional* prediction standard errors (the default, `re.form = NULL`)
come from the engine’s prediction-variance payload, which adds the
random-effect (BLUP) variance and the fixed/random covariance to the
fixed-effect Wald variance — a surface `lme4::predict()` does not offer
at all. Conditional confidence and prediction intervals come from the
same payload.

``` r

pred <- predict(fit, se.fit = TRUE)
head(pred$fit)
#>        1        2        3        4        5        6 
#> 7.585932 7.307533 7.029134 6.750734 6.472335 6.193935
head(pred$se.fit)
#>         1         2         3         4         5         6 
#> 0.1597990 0.1511355 0.1466119 0.1466119 0.1511355 0.1597990
head(predict(fit, interval = "confidence"))
#>        fit      lwr      upr
#> 1 7.585932 7.272732 7.899133
#> 2 7.307533 7.011313 7.603753
#> 3 7.029134 6.741780 7.316488
#> 4 6.750734 6.463380 7.038088
#> 5 6.472335 6.176115 6.768555
#> 6 6.193935 5.880735 6.507136
```

The boundary has moved, not vanished: rows the engine cannot certify —
for example an unseen grouping level predicted with
`allow.new.levels = TRUE`, where there is no posterior variance for the
missing level — still return `NA` with the engine’s reason in the
`mm_reason` attribute rather than a fabricated number.
