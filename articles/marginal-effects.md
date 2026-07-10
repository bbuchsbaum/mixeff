# Marginal Means and Comparisons

``` r

library(mixeff)
```

A fixed-effects coefficient is not a group mean. For a model with
multiple predictors and interactions, `fixef(fit)["trtactive"]` is the
treatment effect *at the reference level of every other predictor* — not
the average treatment effect across the population. Marginal means give
you the latter: population-level averages at each combination of
interest, properly accounting for the reference grid and the
fixed-effect covariance.

`mixeff` provides a native marginal-quantities surface —
[`mm_grid()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
— that routes all inference through the same contract machinery as
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md).
Each row in the returned table carries `method`, `status`,
`reliability`, and `reason` fields so you know exactly what you are
reporting and when to be cautious.

## The study

A rehabilitation trial assigns patients to coached or usual-care
treatment. Each patient is assessed before and after the intervention.
Patients are grouped within clinics.

``` r

head(rehab)
#>   subj clinic     trt time    score
#> 1   S1     C1 control  pre 50.83695
#> 2   S1     C1 control post 46.17501
#> 3   S2     C1  active  pre 48.25095
#> 4   S2     C1  active post 41.31640
#> 5   S3     C1 control  pre 51.03367
#> 6   S3     C1 control post 47.24603
```

Clinics are not perfectly balanced between treatment arms — clinic C1
has more control patients than average, and clinic C3 has more
active-arm patients. This imbalance means raw group means are confounded
by clinic effects: a naive `tapply(rehab$score, rehab$trt, mean)` will
not answer the treatment question you actually want.

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
#>    subj (Intercept) 0.319551 0.565288            
#>  clinic (Intercept) 1.369530 1.170270            
#> Residual std. dev.: 0.396581
#> 
#> Fixed effects:
#>                     Estimate Std. Error       df    t value  Pr(>|t|)
#> (Intercept)        51.875075  0.6181573  3.39845  83.918896 9.313e-07
#> trtactive          -1.978847  0.2819065 26.32727  -7.019516 1.747e-07
#> timepost           -4.002565  0.1619036 21.96772 -24.721911   < 1e-16
#> trtactive:timepost -2.738668  0.2289662 21.96772 -11.961016 4.314e-11
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

The interaction coefficient `trt: active:time: post` tells you the
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
#> 1  trt=control, time=pre 51.87508 50.03208  53.71807 satterthwaite
#> 2   trt=active, time=pre 49.89623 48.05324  51.73922 satterthwaite
#> 3 trt=control, time=post 47.87251 46.02952  49.71550 satterthwaite
#> 4  trt=active, time=post 43.15499 41.31200  44.99799 satterthwaite
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
#> 1 trt=control 49.87379 48.01524  51.73235 satterthwaite
#> 2  trt=active 46.52561 44.66706  48.38417 satterthwaite
```

Compare these to the raw means:

``` r

tapply(rehab$score, rehab$trt, mean)
#>  control   active 
#> 49.87379 46.52561
```

The raw means are shifted by the clinic imbalance; the marginal means
are not. This difference is small in a balanced simulation but can be
substantial in real data.

## Pairwise comparisons with `mm_comparisons()`

[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
takes all pairwise differences among the marginal means and applies the
same inference method.

``` r

ct <- mm_comparisons(fit, specs = ~ trt)
ct$table[, c("label", "estimate", "conf_low", "conf_high", "p_value", "method")]
#>                      label  estimate  conf_low conf_high      p_value
#> 1 trt=active - trt=control -3.348181 -3.887437 -2.808926 6.814682e-11
#>          method
#> 1 satterthwaite
```

The `active - control` row is the average treatment effect across both
timepoints, with a Satterthwaite *t* test and its certified provenance.

## Conditional comparisons with `by =`

The `by` argument splits comparisons within levels of another variable —
the analogue of simple effects in a factorial design.

``` r

ct_by <- mm_comparisons(fit, specs = ~ trt | time)
ct_by$table[, c("label", "estimate", "conf_low", "conf_high", "p_value", "method")]
#>                                            label  estimate  conf_low conf_high
#> 1 trt=active, time=post - trt=control, time=post -4.717515 -5.296632 -4.138399
#> 2   trt=active, time=pre - trt=control, time=pre -1.978847 -2.557964 -1.399730
#>        p_value        method
#> 1 1.554312e-15 satterthwaite
#> 2 1.746510e-07 satterthwaite
```

Two rows: the treatment difference *at pre-intervention* and the
treatment difference *at post-intervention*. The post-intervention gap
is larger because the interaction drives additional improvement in the
active arm.

## Constraining the grid with `at =`

For numeric predictors, `at` pins specific values rather than collapsing
to the mean.

``` r

mt_time <- mm_means(fit, specs = ~ time, at = list(trt = "active"))
mt_time$table[, c("label", "estimate", "conf_low", "conf_high")]
#>       label estimate conf_low conf_high
#> 1  time=pre 49.89623 48.05324  51.73922
#> 2 time=post 43.15499 41.31200  44.99799
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
#>    estimate     lower     upper     p_value        method
#> 1 -2.738668 -3.213556 -2.263781 4.31376e-11 satterthwaite
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
#>  trt     emmean    SE   df lower.CL upper.CL
#>  control   49.9 0.613 3.28     48.0     51.7
#>  active    46.5 0.613 3.28     44.7     48.4
#> 
#> Results are averaged over the levels of: time 
#> mixeff emmeans bridge: fixed-effect covariance from mixedmodels.fixed_effect_covariance_matrix (model_based); prefer mm_means()/mm_comparisons() when row-level status and reasons are needed. 
#> Degrees-of-freedom method: mixeff auto 
#> Confidence level used: 0.95 
#>  contrast         estimate    SE df t.ratio p.value
#>  control - active     3.35 0.258 19  12.997 <0.0001
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
the `emmeans` output. Use `emmeans` for its richer contrasts grammar
(Tukey correction, back-transformation, custom correction methods);
prefer the native verbs when auditability and report-ready provenance
matter.

## Reading `status` and `reason`

Every table returned by the marginal-quantities surface has these
columns:

| Column | Meaning |
|----|----|
| `status` | `"available"` / `"unavailable"` |
| `reliability` | `"certified"` / `"indicative"` / `"unavailable"` |
| `reason` | stable code (see [`vignette("inference-method-glossary")`](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.md)) |
| `method` | the inference method that was applied |

A row with `status = "unavailable"` contains `NA` for standard error and
p-value — the package refuses to invent them. The `reason` code tells
you why (rank deficiency, missing covariance payload, etc.) and is
stable across package versions so you can guard on it in reproducible
scripts.

## Where to read next

- [`vignette("inference")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
  — how the underlying
  [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
  machinery works and what inference methods are available.
- [`vignette("inference-method-glossary")`](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.md)
  — reference table of every `reason` code.
- [`vignette("reporting-lmms")`](https://bbuchsbaum.github.io/mixeff/articles/reporting-lmms.md)
  — building a full results section from a fitted object.
