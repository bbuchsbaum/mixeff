# Marginal Means and Comparisons

``` r

library(mixeff)
```

A fixed-effects coefficient is not a group mean. For a model with
multiple predictors and interactions, `fixef(fit)["trtactive"]` is the
treatment effect *at the reference level of every other predictor* — not
the average treatment effect across the population. Marginal means give
you the latter: averages over a reference grid of predictor
combinations, with standard errors from the fixed-effect covariance.

`mixeff` computes these with four functions —
[`mm_grid()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
— that route all inference through the same machinery as
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md).
Each row in the returned table carries `method`, `status`,
`reliability`, and `reason` columns naming how it was computed and how
far to trust it.

## The study

A rehabilitation trial assigns patients to coached or usual-care
treatment. Each patient is assessed before and after the intervention.
Patients are grouped within clinics.

``` r

head(rehab)
#>   subj clinic     trt time    score
#> 1   S1     C1 control  pre 50.83695
#> 2   S1     C1 control post 46.17501
#> 3   S2     C1 control  pre 50.35095
#> 4   S2     C1 control post 46.31640
#> 5   S3     C1 control  pre 51.03367
#> 6   S3     C1 control post 47.24603
```

Clinics are not balanced between treatment arms: clinic C1 has five
control patients and one active, clinic C3 the reverse. So raw group
means are confounded with clinic effects — each arm’s mean leans on the
clinics that arm over-samples, and the clinics differ — and a naive
`tapply(rehab$score, rehab$trt, mean)` will not answer the treatment
question.

## Fit the model

``` r

fit <- lmm(score ~ trt * time + (1 | clinic) + (1 | subj), rehab)
summary(fit, tests = "coefficients")
#> Linear mixed model fit by REML
#> Formula: score ~ trt * time + (1 | clinic) + (1 | subj)
#> Fit status: converged_interior
#> 
#> Variance components:
#>   group        name variance  std_dev correlation
#>    subj (Intercept) 0.328791 0.573403            
#>  clinic (Intercept) 1.267670 1.125910            
#> Residual std. dev.: 0.399888
#> 
#> Fixed effects:
#>                     Estimate Std. Error        df    t value  Pr(>|t|)
#> (Intercept)        51.871482  0.6018249  3.425072  86.190327 7.755e-07
#> trtactive          -1.971661  0.3156838 25.438945  -6.245684 1.443e-06
#> timepost           -3.962867  0.1632535 22.009671 -24.274316   < 1e-16
#> trtactive:timepost -2.818065  0.2308753 22.009671 -12.206006 2.848e-11
#>                           method
#> (Intercept)        satterthwaite
#> trtactive          satterthwaite
#> timepost           satterthwaite
#> trtactive:timepost satterthwaite
#> 
#> Inference status:
#>                term        method    status reliability
#>         (Intercept) satterthwaite available    moderate
#>           trtactive satterthwaite available    moderate
#>            timepost satterthwaite available    moderate
#>  trtactive:timepost satterthwaite available    moderate
#>                             reliability_reason
#>  satterthwaite_finite_difference_approximation
#>  satterthwaite_finite_difference_approximation
#>  satterthwaite_finite_difference_approximation
#>  satterthwaite_finite_difference_approximation
#> 
#> Notes:
#>   Satterthwaite denominator df computed from finite-difference vcov_beta Jacobian and deviance Hessian over varpar
```

The interaction coefficient `trtactive:timepost` tells you the
*additional* post-treatment change for the active arm relative to
control. It is not the average treatment effect. For that you need
marginal means.

## Reference grids with `mm_grid()`

[`mm_grid()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
constructs the cross-product of all fixed-predictor levels. By default,
factor predictors expand to all their levels; numeric predictors
collapse to their mean.

``` r

g <- mm_grid(fit, specs = ~ trt * time)
g
#> Marginal grid:
#>      trt time
#>  control  pre
#>   active  pre
#>  control post
#>   active post
```

The grid has four rows — one for each treatment × timepoint cell — and
retains the model matrix needed for inference.

## Cell predictions with `mm_predictions()`

[`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
evaluates the fixed-effect prediction at each grid row, with a
confidence interval from the certified covariance.

``` r

preds <- mm_predictions(fit, specs = ~ trt * time)
preds$table[, c("label", "estimate", "conf_low", "conf_high", "method")]
#>                    label estimate conf_low conf_high        method
#> 1  trt=control, time=pre 51.87148 50.08390  53.65906 satterthwaite
#> 2   trt=active, time=pre 49.89982 48.11224  51.68740 satterthwaite
#> 3 trt=control, time=post 47.90862 46.12103  49.69620 satterthwaite
#> 4  trt=active, time=post 43.11889 41.33131  44.90647 satterthwaite
```

These are the four population-level cell means. Each row carries its
inference method so the provenance is visible without digging into model
objects.

## Marginal means with `mm_means()`

Marginal means average the reference grid over dimensions you want to
collapse. Here: average over timepoints to get the *overall* treatment
effect.

``` r

mt <- mm_means(fit, specs = ~ trt)
mt$table[, c("label", "estimate", "conf_low", "conf_high", "method")]
#>         label estimate conf_low conf_high        method
#> 1 trt=control 49.89005 48.08653  51.69356 satterthwaite
#> 2  trt=active 46.50935 44.70584  48.31287 satterthwaite
```

Compare these to the raw means:

``` r

tapply(rehab$score, rehab$trt, mean)
#>  control   active 
#> 49.45125 46.94815
```

The clinic imbalance pulls each raw mean about 0.44 points toward the
clinics its arm over-samples; the marginal means average over the
model’s fixed-effect grid instead and are free of that shift. The size
of the gap depends on how unbalanced the design is and how large the
group effects are.

## Pairwise comparisons with `mm_comparisons()`

[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
takes all pairwise differences among the marginal means and applies the
same inference method.

``` r

ct <- mm_comparisons(fit, specs = ~ trt)
ct$table[, c("label", "estimate", "conf_low", "conf_high", "p_value", "method")]
#>                      label  estimate conf_low conf_high     p_value
#> 1 trt=active - trt=control -3.380694 -3.99461 -2.766778 3.86692e-10
#>          method
#> 1 satterthwaite
```

The `active - control` row is the average treatment effect across both
timepoints, tested with a Satterthwaite *t* and graded `moderate`
reliability in the same row.

## Conditional comparisons with `by =`

The `by` argument splits comparisons within levels of another variable —
the analogue of simple effects in a factorial design.

``` r

ct_by <- mm_comparisons(fit, specs = ~ trt | time)
ct_by$table[, c("label", "estimate", "conf_low", "conf_high", "p_value", "method")]
#>                                            label  estimate  conf_low conf_high
#> 1 trt=active, time=post - trt=control, time=post -4.789727 -5.439321 -4.140132
#> 2   trt=active, time=pre - trt=control, time=pre -1.971661 -2.621256 -1.322067
#>        p_value        method
#> 1 2.953193e-14 satterthwaite
#> 2 1.442769e-06 satterthwaite
```

Two rows: the treatment difference *at pre-intervention* and the
treatment difference *at post-intervention*. The post-intervention gap
is larger because the interaction drives additional improvement in the
active arm.

## Constraining the grid with `at =`

`at` pins predictors to chosen values instead of expanding or averaging
them: a factor to a subset of its levels, a numeric predictor to
specific points rather than its mean.

``` r

mt_time <- mm_means(fit, specs = ~ time, at = list(trt = "active"))
mt_time$table[, c("label", "estimate", "conf_low", "conf_high")]
#>       label estimate conf_low conf_high
#> 1  time=pre 49.89982 48.11224  51.68740
#> 2 time=post 43.11889 41.33131  44.90647
```

This gives the pre/post means *within the active arm* only, holding
`trt` constant at `"active"`.

## Custom contrasts with `mm_lincomb()`

For hypotheses that are not pairwise differences of marginal means,
build the contrast weights directly with
[`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md).
The interaction effect expressed as a contrast: (active post − active
pre) − (control post − control pre).

``` r

beta <- fixef(fit)
names(beta)
#> [1] "(Intercept)"        "trtactive"          "timepost"          
#> [4] "trtactive:timepost"

# Interaction row: coefficient named "trtactive:timepost" (lme4-identical)
w <- setNames(numeric(length(beta)), names(beta))
w["trtactive:timepost"] <- 1

lc <- mm_lincomb(fit, weights = w)
lc[, c("estimate", "lower", "upper", "p_value", "method")]
#>    estimate     lower     upper      p_value        method
#> 1 -2.818065 -3.296859 -2.339271 2.848244e-11 satterthwaite
```

[`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md)
applies the same contract-preserving inference as
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
so the method and status fields are populated identically.

## The `emmeans` bridge

When `emmeans` is installed, `mixeff` registers a basis so you can call
[`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
directly on `mm_lmm` objects.

``` r

if (requireNamespace("emmeans", quietly = TRUE)) {
  em <- emmeans::emmeans(fit, ~ trt)
  print(em)
  print(pairs(em))
}
#>  trt     emmean    SE  df lower.CL upper.CL
#>  control   49.9 0.596 3.3     48.1     51.7
#>  active    46.5 0.596 3.3     44.7     48.3
#> 
#> Results are averaged over the levels of: time 
#> mixeff emmeans bridge: fixed-effect covariance from mixedmodels.fixed_effect_covariance_matrix (model_based); prefer mm_means()/mm_comparisons() when row-level status and reasons are needed. 
#> Degrees-of-freedom method: mixeff auto 
#> Confidence level used: 0.95 
#>  contrast         estimate    SE   df t.ratio p.value
#>  control - active     3.38 0.294 19.5  11.506 <0.0001
#> 
#> Results are averaged over the levels of: time 
#> Degrees-of-freedom method: mixeff auto
```

`emmeans` uses the same `mixedmodels.fixed_effect_covariance_matrix`
payload as
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
so the point estimates and standard errors agree. The bridge prints an
informational message noting which covariance it used.

**When to prefer the native verbs over `emmeans`:**
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
and
[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
carry the full `status`, `reliability`, and `reason` row-level audit
fields from the underlying
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
call. `emmeans` does not propagate these fields; if a row is
`"unavailable"` for a documented reason, that information disappears in
the `emmeans` output. Use `emmeans` when you need what it does and
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
does not — Tukey correction, back-transformation, custom contrast
grammars. Use the native verbs when you need each row to say how it was
computed and why, because those columns do not survive the trip through
`emmeans`.

## Reading `status` and `reason`

Every table returned by
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
and
[`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
carries these columns:

| Column | Values |
|----|----|
| `status` | `"available"`, `"p_value_unavailable"`, `"not_estimable"`, `"not_assessed"`, or `"unsupported"` |
| `reliability` | `"high"`, `"moderate"`, `"low"`, or `"not_available"` |
| `reason` | stable code when a row is degraded or refused, `NA` otherwise (see [`vignette("inference-method-glossary")`](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.md)) |
| `method` | the inference method that was applied |

On the fits in this vignette every row reads `status = "available"`,
`reliability = "moderate"`, `method = "satterthwaite"`. A row whose
status is anything other than `"available"` contains `NA` where the
refused numbers would go, and `reason` says why (rank deficiency,
missing covariance payload, and so on). These codes come from the
engine’s fixed vocabulary, so scripts can test for them by string.

## Where to read next

- [`vignette("inference")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
  — how the underlying
  [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
  machinery works and what inference methods are available.
- [`vignette("inference-method-glossary")`](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.md)
  — reference table of every `reason` code.
- [`vignette("reporting-lmms")`](https://bbuchsbaum.github.io/mixeff/articles/reporting-lmms.md)
  — building a full results section from a fitted object.
