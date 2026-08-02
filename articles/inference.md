# Inference

``` r

library(mixeff)
```

Mixed-model users want p-values. They also want to know what those
p-values mean: which method produced them, and how trustworthy that
method is for this particular fit.
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) famously
declines to print LMM p-values at all — the reason `lmerTest` exists.
`mixeff` prints them, and prints the method and a reliability grade
beside each one.

When a coefficient, contrast, or term test has an available method, the
p-value is printed with the method name. When the requested method is
unavailable on this fit, the row says so, with a stable reason code
rather than an apologetic substitute. The full vocabulary of reason
codes, and the grading rubric behind the `reliability` column, live in
the inference-method glossary on the package site
(<https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.html>);
this vignette points there rather than restating them.

This vignette covers the whole inference surface: coefficient tests,
contrasts, term tests, and model comparisons on a well-behaved fit; the
full route map — including what is refused at a boundary fit and why;
the GLMM three-route contract; and the marginal-means verbs, starting
with the question every marginal summary should answer first: *which
estimand are you averaging toward?*

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
Here `"auto"` resolved to Satterthwaite, because the fit is interior and
the finite-sample route is feasible.

[`inference_table()`](https://bbuchsbaum.github.io/mixeff/reference/inference_table.md)
shows a different, lower-level view: the inference rows the engine
cached at fit time. Those rows use asymptotic Wald *z* — the engine’s
always-computable baseline — so their p-values are smaller than the
Satterthwaite ones above, and their `reliability` column says `low`.
This is the raw artifact, useful for audit; for reporting, use
[`summary()`](https://rdrr.io/r/base/summary.html) or
[`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md),
which resolve to the best feasible method and agree with each other row
for row.

``` r

inference_table(fit)
#> Inference table:
#>              term            label        kind   estimate  std_error df
#>       (Intercept)      (Intercept) coefficient  7.6828778 0.19646018 NA
#>              week             week coefficient -0.2783994 0.02595083 NA
#>  treatmentcoached treatmentcoached coefficient -0.8994747 0.26225014 NA
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

`type` selects the term hypothesis. `"III"` (the default) is marginal
Type III: a term’s own contrast columns plus the equally-weighted
average over the levels of every term that contains it. That
construction makes the test invariant to which factor level is the
reference — the defining property of Type III, and what SAS, `car`, and
`lmerTest` compute. `"II"` respects marginality; `"I"` is sequential in
[`terms()`](https://rdrr.io/r/stats/terms.html) order.

There is a fourth value, `"block"`, for the raw coefficient block of
each term. Under treatment coding that is the *simple effect at the
other factors’ reference levels* — a legitimate quantity, but not Type
III. On an unbalanced design the two can disagree sharply: on the
package’s test fixture the same fit gives F = 0.29 for a main effect
under `"block"` and F = 49.95 under `"III"`. If you want simple effects,
ask for them by name; if you want Type III, the default is Type III.

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
simulation-based check for a small example, use the bootstrap path. The
bootstrap p-value is `(b + 1) / (successful + 1)`, where `b` counts
replicates at or above the observed statistic and `successful` counts
the refits that completed with a finite statistic — equal to `nsim` only
when every replicate succeeds, as it does in these seeded runs. So at
`nsim = 10` the smallest reportable value is 1/11 — use hundreds of
replicates before reading the p-value as a number rather than a bound.

``` r

compare(reduced, fit, method = "bootstrap", nsim = 10, seed = 7)
#> Model comparison:
#>  model                                     formula nobs df    logLik deviance
#>     m1             score ~ 1 + week + (1 | clinic)   72  4 -47.83234 95.66467
#>     m2 score ~ 1 + week + treatment + (1 | clinic)   72  5 -43.16629 86.33259
#>        AIC      BIC delta_aic delta_bic  REML refit         fit_status delta_df
#>  103.66467 112.7713  7.332086   5.05542 FALSE  TRUE converged_interior       NA
#>   96.33259 107.7159  0.000000   0.00000 FALSE  TRUE converged_interior        1
#>       LRT    p_value                   method          status
#>        NA         NA           asymptotic_lrt reference_model
#>  9.332086 0.09090909 parametric_bootstrap_lrt       available
#>                                                     reason reason_code
#>                                                                   <NA>
#>  parametric bootstrap LRT (10/10 replicates, MCSE=0.09091)        <NA>
#>      comparison_class lrt_available information_criteria_available
#>                  <NA>         FALSE                           TRUE
#>  nested_fixed_effects          TRUE                           TRUE
#>  requires_ml_refit loglik_within_optimizer_tol rust_method rust_refit_policy
#>              FALSE                       FALSE        auto             never
#>              FALSE                       FALSE        auto             never
```

## When standard p-values are not well justified

Mixed models sometimes converge but still leave us in a difficult
inference situation. The fit exists, the estimates are usable, and the
model may be scientifically reasonable, but some familiar p-values are
no longer well justified.

The common case is a boundary or singular fit: a variance component
estimated at or near zero, or random-effect terms estimated as perfectly
correlated — equivalently, a fitted random-effect covariance with lower
effective rank than the model requested. What those labels mean, how
they are recorded on the object, and how singularity is decided are
covered in
[`vignette("convergence")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md);
here we take the labels as given and ask what they do to inference.

The statistical background is well known. Self and Liang (1987) showed
that likelihood-ratio statistics do not necessarily follow their usual
chi-square reference distributions when a parameter lies on the boundary
of its parameter space. Stram and Lee (1994) specialized this issue to
variance components in linear mixed models. Kenward and Roger (1997)
derived their fixed-effect adjustment from a Taylor expansion of the
variance estimator; that derivation assumes a regular interior solution
for the variance parameters, and at a boundary those derivatives may not
exist or may not behave regularly. Parametric bootstrap methods have a
long history as the better-supported route in this setting; Halekoh and
Højsgaard (2014) implemented them in `pbkrtest` for small-sample and
boundary cases.

`mixeff` handles this by making the inference rules explicit. Asymptotic
rows are labeled with a reliability reason from a fixed vocabulary, such
as `asymptotic_wald_z_fallback` or `*_unavailable_at_boundary`, and
bootstrap inference is exposed as a labeled method rather than an
informal workaround. The goal is not to say that a singular fit is
automatically useless; it is to keep the covariance state attached to
the inference result, so that available numbers and unavailable numbers
are both documented.

This small repeated-measures example has subject-specific intercepts and
slopes that are nearly collinear by construction — enough to produce a
reduced-rank random-effect covariance estimate.

``` r

sfit <- lmm(
  rt ~ days + (1 + days | subj),
  sleep_like,
  control = mm_control(verbose = -1)
)

fit_status(sfit)
#> [1] "converged_reduced_rank"
is_singular(sfit)
#> [1] TRUE
```

## The route map

Before testing a term, ask which inference routes are defined for this
fit.
[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
is an audit map, not a recommender. Each row describes one method that
`mixeff` knows about. The row says whether the method can run, why it
can or cannot run, and what the rough computational cost is.

The printed object is the reader-facing map. It keeps the raw enum
columns in `routes$table` for scripts, but explains the route status
with display columns. The complete contract vocabulary is in the
[inference-method
glossary](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.html).

``` r

routes <- inference_options(sfit, "days", nsim = 200)
routes
#> Inference options (fit_status: converged_reduced_rank, REML: TRUE):
#>             method      display_status
#>  asymptotic_wald_z            runs now
#>      satterthwaite refused on this fit
#>      kenward_roger refused on this fit
#>          bootstrap            runs now
#>      bootstrap_lrt refused on this fit
#>  cluster_bootstrap refused on this fit
#>         profile_ci refused on this fit
#>                                                                 display_reason
#>                                                     asymptotic wald z fallback
#>                            variance-parameter derivative undefined at boundary
#>                            variance-parameter derivative undefined at boundary
#>                                       calibrated by nsim and Monte Carlo error
#>                                    requires an ML fit; refit with REML = FALSE
#>  cluster resampling reports estimator distributions, not fixed-effect p-values
#>                            profile intervals are not certified at the boundary
#>                                                                                                                          what_to_do_next
#>                                                                                                                             summary(fit)
#>                                  Use asymptotic_wald_z or bootstrap; simplify the random-effects structure if the boundary is unintended
#>                                  Use asymptotic_wald_z or bootstrap; simplify the random-effects structure if the boundary is unintended
#>                                                  test_effect(fit, term, method = 'bootstrap', bootstrap = bootstrap_control(nsim = 200))
#>  Refit with lmm(..., REML = FALSE), then run test_effect(fit, term, method = 'bootstrap_lrt', bootstrap = bootstrap_control(nsim = 200))
#>                                                                                 Use bootstrap or bootstrap_lrt for fixed-effect p-values
#>                                  Use asymptotic_wald_z or bootstrap; simplify the random-effects structure if the boundary is unintended
#>                    approx_cost current
#>                      immediate    TRUE
#>                      immediate   FALSE
#>                      immediate   FALSE
#>  200 model refits (nsim = 200)   FALSE
#>  400 model refits (nsim = 200)   FALSE
#>                              -   FALSE
#>            slow profile refits   FALSE
#> 
#> Use `<obj>$table` for raw enum columns (`expected_status`, `expected_reliability_reason`) and notes.
```

Read this table as follows. The Wald route can run immediately, but it
is labeled as a low-reliability asymptotic fallback. Satterthwaite and
Kenward–Roger are refused on this fit by mixeff’s validated-support
policy: both approximations are derived under asymptotics that assume
the variance parameters lie in the interior of their parameter space,
and this fit has a variance parameter on the boundary. Other packages
will compute degrees of freedom here (`lmerTest` returns a number on
this same fit); `mixeff` declines to, and states the refusal as a stable
reason code — `satterthwaite_unavailable_at_boundary`,
`kenward_roger_unavailable_at_boundary` — you can test for. That is a
policy about what the package will certify, not a theorem that the
numbers other packages print are wrong.

The parametric bootstrap route can run and is the documented
fixed-effect testing route for this fit. The bootstrap LRT route is
refused here because this model was fit by REML, while that route
requires an ML fit (`bootstrap_lrt_requires_ml`). The cluster bootstrap
route is refused for fixed-effect p-values because its current target is
an estimator distribution, not a null hypothesis test. Profile
confidence intervals are not certified for this boundary fit
(`profile_ci_unavailable_at_boundary`).

``` r

inference_table(sfit)$table[, c("term", "method", "status",
                                "reliability", "reliability_reason")]
#>          term            method    status reliability
#> 1 (Intercept) asymptotic_wald_z available         low
#> 2        days asymptotic_wald_z available         low
#>           reliability_reason
#> 1 asymptotic_wald_z_fallback
#> 2 asymptotic_wald_z_fallback
```

`asymptotic_wald_z_fallback` on these rows is the engine’s reason: a *t*
reference distribution was the target, the degrees of freedom could not
be computed at this boundary fit, and a standard normal was substituted
— which is why the row is graded `low` rather than hidden. The *z* test
drops the finite-sample correction, so treat its p-values as optimistic
rather than cautious.

The important point is that `mixeff` distinguishes three cases:

1.  A method can run and is considered available.
2.  A method can produce a number, but the number is labeled as low
    reliability.
3.  A method is refused because the method’s assumptions or target are
    not defined for this fit.

That distinction is especially useful in scripts and reports, because
unavailable inference is recorded with stable reason codes rather than
disappearing as an error message.

## How do I test the term anyway?

For the fixed-effect term `days`, request the parametric bootstrap route
explicitly. Here the bootstrap target is a fixed-effect null model.
`mixeff` simulates data under the null hypothesis for the tested term,
refits the model to each simulated data set, and compares the observed
test statistic with the bootstrap reference distribution.

The small `nsim` value below keeps the vignette fast, and it earns a
`low` reliability grade: per the package’s grade rubric (see the
inference-method glossary linked above), a certified bootstrap reaches
`moderate` at 999 or more successful replicates with a finite Monte
Carlo SE, and stays `low` below that floor. The reason to use bootstrap
here is not the grade — it is that the p-value rests on simulation
rather than on asymptotics that mixeff refuses to certify at the
boundary.

``` r

term_boot <- test_effect(
  sfit,
  "days",
  method = "bootstrap",
  bootstrap = bootstrap_control(nsim = 50, seed = 1)
)

term_boot$table[, c("term", "statistic_name", "p_value", "method",
                    "status", "reliability", "reliability_reason")]
#>   term statistic_name    p_value    method    status reliability
#> 1 days              t 0.01960784 bootstrap available         low
#>                 reliability_reason
#> 1 parametric_bootstrap_monte_carlo
```

The p-value is returned as part of the labeled result row. It is not
reconstructed by hand in R.

The result also carries a bootstrap payload in `details`, including the
requested number of replicates, the number of successful refits, the
boundary rate among bootstrap refits, and the Monte Carlo standard
error.

``` r

run <- term_boot$table$details[[1]]$bootstrap
data.frame(
  requested_replicates = run$requested_replicates,
  successful_replicates = run$successful_replicates,
  boundary_rate = round(run$boundary_rate, 3),
  mcse = round(run$mcse, 4)
)
#>   requested_replicates successful_replicates boundary_rate   mcse
#> 1                   50                    50          0.56 0.0196
```

The boundary rate is not a failure count. It reports how often the
bootstrap refits also ended on a covariance boundary. That is useful
diagnostic information: it tells you that the boundary behavior is not
just a one-off feature of the observed data set.

With `nsim = 50`, the bootstrap p-value is also necessarily coarse. For
example, a p-value of about 0.0196 corresponds to the smallest non-zero
value possible under a common plus-one bootstrap adjustment,
`1 / (successful + 1)` — here all 50 replicates succeeded, as the
accounting above shows, so 1/51. Use more replicates when the p-value
itself will be reported or compared with a decision threshold.

## What about confidence intervals?

For fixed-effect confidence intervals, use
`confint(method = "bootstrap")`. This bootstrap has a different target
from the fixed-effect test above. It simulates from the full fitted
model and summarizes the resulting estimator distribution. That target
is appropriate for percentile-style confidence intervals and
diagnostics, but it is not the same as a null distribution for a
fixed-effect p-value.

``` r

ci <- confint(
  sfit,
  parm = "days",
  method = "bootstrap",
  bootstrap = bootstrap_control(nsim = 50, seed = 2)
)

ci
#> Confidence intervals:
#>         2.5 %   97.5 %
#> days 7.588973 10.68778
#> method: bootstrap_full_model_distribution
#> interval: percentile
#> status: available
#> 
#> Bootstrap run:
#>  parameter requested successful failed_refits boundary_rate seed
#>       days        50         50             0          0.58    2
#> notes:
#>   - full-model bootstrap distributions do not certify fixed-effect
#>         hypothesis-test p-values
#>   - 29 successful bootstrap refit(s) ended on a covariance boundary
#> Full bootstrap payload available in `attr(x, "bootstrap")`.
```

The attached bootstrap payload records what was simulated and why the
interval should not be reinterpreted as a fixed-effect hypothesis test.

``` r

payload <- attr(ci, "bootstrap")[[1]]
payload$metadata$target$kind
#> [1] "full_model_distribution"
payload$metadata$notes
#> [[1]]
#> [1] "full-model bootstrap distributions do not certify fixed-effect hypothesis-test p-values"
#> 
#> [[2]]
#> [1] "29 successful bootstrap refit(s) ended on a covariance boundary"
```

This distinction is central:

- `test_effect(..., method = "bootstrap")` uses a fixed-effect null
  target and returns a p-value.
- `confint(..., method = "bootstrap")` uses a full-model
  estimator-distribution target and returns an interval.
- The second result should not be used to reverse-engineer the first.

## Where does cluster bootstrap stand?

Cluster resampling is available as an estimator-distribution target. It
resamples grouping levels and summarizes how estimates vary across those
resamples.

In the current schema, that target does not certify fixed-effect
p-values. Therefore a request for a fixed-effect cluster-bootstrap
p-value refuses cleanly instead of inventing a null distribution.

``` r

cluster_row <- test_effect(
  sfit,
  "days",
  method = "cluster_bootstrap",
  group = "subj"
)

cluster_row$table[, c("term", "method", "status", "p_value", "reason_code")]
#>   term            method       status p_value
#> 1 days cluster_bootstrap not_assessed      NA
#>                                      reason_code
#> 1 bootstrap_cluster_resample_p_value_unavailable
```

That refusal is part of the same inference contract as the available
bootstrap result. Available numbers are labeled with their method and
reliability reason. Unavailable numbers are labeled with stable reason
codes explaining why they are unavailable.

## Variance-component boundary tests

The boundary likelihood-ratio route is for random-effect variance
components, not ordinary fixed effects.

This is a different testing problem. For a variance component, the null
hypothesis can put the parameter exactly on the boundary of the
parameter space, because variances cannot be negative. In the simplest
one-component case, the *asymptotic* reference distribution is the
Self-Liang 50:50 mixture, `0.5 * chi-square(0) + 0.5 * chi-square(1)`.
In finite samples the true null distribution puts more mass at zero than
the mixture, so the test runs conservative (Crainiceanu and Ruppert
2004); with 18 subjects here, read a marginal mixture p-value with that
in mind.

Fit a simpler random-intercept model by ML:

``` r

ri_fit <- lmm(
  rt ~ days + (1 | subj),
  sleep_like,
  REML = FALSE,
  control = mm_control(verbose = -1)
)

re_lrt <- test_random_effect(ri_fit, "subj", method = "boundary_lrt")
re_lrt
#> Random-effect variance-component test:
#>        term group statistic statistic_name p_value
#>  (1 | subj)  subj  178.8188 chi_bar_square       0
#>                     reference_distribution    status reason_code
#>  0.5 * chi-square(0) + 0.5 * chi-square(1) available        <NA>
```

``` r

reporting_table(re_lrt)$table[, c("term", "statistic", "p_value",
                                  "reference_distribution", "status")]
#>         term statistic p_value                    reference_distribution
#> 1 (1 | subj)  178.8188       0 0.5 * chi-square(0) + 0.5 * chi-square(1)
#>      status
#> 1 available
```

The reporting table keeps the reference distribution attached to the
result.

Asking for the same method on a fixed effect returns a typed refusal.

``` r

fixed_boundary <- test_effect(ri_fit, "days", method = "boundary_lrt")
fixed_boundary$table[, c("term", "method", "status", "reason_code")]
#>   term         method      status                                  reason_code
#> 1 days not_applicable unsupported boundary_lrt_not_applicable_to_fixed_effects
```

Boundary likelihood-ratio theory applies to variance-component tests,
where the null hypothesis places a variance on the edge of the parameter
space. A fixed-effect null does not put any parameter on a boundary, so
the mixture distribution has nothing to correct there — which is why the
route refuses rather than running.

## Profile confidence intervals

Profile-likelihood intervals are available for variance parameters under
both fitting criteria, but fixed-effect (beta) profile intervals require
an ML fit: under the REML default, the beta rows are recorded as
unavailable with the reason code `profile_beta_unavailable_under_reml`,
not silently approximated.

``` r

prof <- confint(fit, method = "profile")
prof
#> Confidence intervals:
#>                      2.5 %    97.5 %
#> theta1           0.6712410 2.0292158
#> sigma            0.3172088 0.4556499
#> (Intercept)             NA        NA
#> week                    NA        NA
#> treatmentcoached        NA        NA
#> method: profile_likelihood
#> status: available
attr(prof, "mm_profile")$table[
  ,
  c("parameter", "parameter_kind", "lower", "upper", "reason_code")
]
#>          parameter parameter_kind     lower     upper
#> 1           theta1          theta 0.6712410 2.0292158
#> 2            sigma          sigma 0.3172088 0.4556499
#> 3      (Intercept)           beta        NA        NA
#> 4             week           beta        NA        NA
#> 5 treatmentcoached           beta        NA        NA
#>                           reason_code
#> 1                                <NA>
#> 2                                <NA>
#> 3 profile_beta_unavailable_under_reml
#> 4 profile_beta_unavailable_under_reml
#> 5 profile_beta_unavailable_under_reml
```

Refit with `REML = FALSE` when you want beta profile intervals. On
boundary fits, profile intervals are refused outright
(`profile_ci_unavailable_at_boundary`), as the route map above showed.

## GLMM inference: the three-route contract

For generalized linear mixed models the contract is stricter, because
the default estimator is not certified for Wald inference.
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)’s
default `method = "pirls_profiled"` is a fast profiled PIRLS fitter: it
returns point estimates, fitted values, and variance components, but
**refuses Wald standard errors, z statistics, p-values, and confidence
intervals on every route**. Its covariance payload is a working Hessian
the engine does not certify, and the engine’s own documentation warns
that the profiled path is a different statistical approximation than
`glmer()`’s joint Laplace fit — it can be less accurate for
overdispersed models and models with observation-level random effects.

Exactly three routes lead to reportable GLMM inference, and every
refusal names them:

1.  **Certified joint Wald.** Refit with `method = "joint_laplace"`;
    standard errors, z tests, and Wald intervals come from a covariance
    the engine certifies.
    [`vignette("glmm")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md)
    is the operational quickstart for this route, including its parity
    with `glmer()` and its runtime notice.
2.  **Parametric bootstrap on profiled fits.**
    `confint(fit, method = "bootstrap")` simulates from the fitted model
    and refits the same profiled estimator. It works for every supported
    family, including negative binomial — which has no joint-Laplace
    route, so the bootstrap is its only interval path — and costs `nsim`
    refits. On `joint_laplace` fits the bootstrap is refused at this
    engine pin; those fits have certified Wald intervals instead.
3.  **Working-Hessian approximation, explicit opt-in.** For exploration
    and screening only. Every unlocked row is labeled
    `wald_z_working_hessian` with status `available_noninferential` and
    grade `moderate`; `inference_options(fit)` states the caveat and the
    invocation, and the aphantasia reproduction article
    (<https://bbuchsbaum.github.io/mixeff/articles/reproducing-aphantasia.html>)
    is the worked example of why this route is opt-in.

The practical rule: explore with the default, report with
`joint_laplace` or the bootstrap.

Here is the contract on a small simulated binomial GLMM. On the default
fit, [`summary()`](https://rdrr.io/r/base/summary.html) prints the
estimates and withholds the inference columns:

``` r

gfit <- glmm(
  y ~ x + (1 | g),
  trials,
  family = binomial(),
  control = mm_control(verbose = -1)
)

gtab <- summary(gfit, tests = "coefficients")$coefficients
gtab
#>               Estimate Std. Error statistic p.value       method
#> (Intercept) -0.3993006         NA        NA      NA not_computed
#> x            0.6448091         NA        NA      NA not_computed
```

A Wald [`confint()`](https://rdrr.io/r/stats/confint.html) on the same
fit is a typed refusal, not an error to guess about:

``` r

err <- tryCatch(confint(gfit), error = identity)
class(err)[[1L]]
#> [1] "mm_inference_unavailable"
writeLines(strwrap(conditionMessage(err), width = 72))
#> confint() Wald intervals are withheld for this GLMM estimator: the
#> engine does not certify the profiled working-Hessian covariance for
#> Wald inference. Use confint(fit, method = "bootstrap") for
#> parametric-bootstrap intervals (all supported families), or refit with
#> method = "joint_laplace" for certified Wald inference (not available
#> for negative binomial). inference_options(fit) lists every route.
```

The bootstrap route works on the same fit, and the interval carries its
accounting: requested, successful, and failed replicate counts,
replicate-SD standard errors, the Monte Carlo SE, and a reliability
grade. The replicate floor for a `moderate` grade is the rubric’s
999-successful-replicates convention (see the glossary linked above);
the small, seeded run below stays `low` by that rubric.

``` r

gci <- confint(gfit, method = "bootstrap", nsim = 30, seed = 42)
gci
#>                  2.5 %    97.5 %
#> (Intercept) -0.9729923 0.3359735
#> x            0.1704129 1.2433597
#> attr(,"mm_method")
#> [1] "glmm_parametric_bootstrap_percentile"
#> attr(,"mm_estimator")
#> [1] "pirls_profiled"
#> attr(,"mm_bootstrap")
#> attr(,"mm_bootstrap")$requested
#> [1] 30
#> 
#> attr(,"mm_bootstrap")$successful
#> [1] 30
#> 
#> attr(,"mm_bootstrap")$failed
#> [1] 0
#> 
#> attr(,"mm_bootstrap")$seed
#> [1] 42
#> 
#> attr(,"mm_bootstrap")$std_errors
#> (Intercept)           x 
#>   0.3686198   0.2927409 
#> 
#> attr(,"mm_bootstrap")$mcse
#> [1] 0.06730046
#> 
#> attr(,"mm_bootstrap")$reliability
#> attr(,"mm_bootstrap")$reliability$reliability
#> [1] "low"
#> 
#> attr(,"mm_bootstrap")$reliability$reason
#> [1] "bootstrap_insufficient_replicates"

acct <- attr(gci, "mm_bootstrap")
data.frame(
  requested = acct$requested,
  successful = acct$successful,
  failed = acct$failed,
  mcse = round(acct$mcse, 4),
  reliability = acct$reliability$reliability
)
#>   requested successful failed   mcse reliability
#> 1        30         30      0 0.0673         low
attr(gci, "mm_method")
#> [1] "glmm_parametric_bootstrap_percentile"
```

The estimator mechanics — how the two estimators differ, the `nAGQ`
sensitivity check, and the joint-Laplace parity with `glmer()` — are
covered in
[`vignette("glmm")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md).
If a `joint_laplace` request does not certify and the engine returns a
labelled fallback fit, that substitution is a typed record with its own
notices; see
[`vignette("convergence")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md).

## Estimands: what are you averaging?

A fixed-effect coefficient is not a group mean. For a model with
multiple predictors and interactions, `fixef(fit)["trtactive"]` is the
treatment effect *at the reference level of every other predictor* — one
cell of the design, not a summary over it. Before averaging cells into
something more interpretable, name the estimand. Four quantities are
easy to conflate:

1.  **A prediction at a covariate cell** — the model’s fixed-effect
    prediction at one row of the reference grid. This is what
    [`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
    returns.
2.  **An equal-weight reference-grid average** — the average of cell
    predictions over the grid, weighting every cell equally regardless
    of how many observations fall in it. This is what
    [`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
    returns, and what
    [`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
    differences.
3.  **An empirical-distribution average** — the average of predictions
    over the observed covariate rows, weighted as the data are weighted.
    `mixeff` has no dedicated verb for this; it is what you get by
    averaging [`predict()`](https://rdrr.io/r/stats/predict.html) over
    the model frame.
4.  **A causal average treatment effect** — the expected difference in
    potential outcomes under treatment versus control. Only this one is
    a causal quantity, and no mixed-model contrast supplies it by
    itself: it requires identification assumptions (randomization, or
    conditional ignorability given the model’s covariates) that are
    properties of the study design, not of the fit.

The first three are model-based summaries that differ only in how they
weight covariate cells; the marginal verbs below compute (1) and (2)
with standard errors from the fixed-effect covariance. Whether any of
them deserves a causal reading is a question about the design, and the
fitted model does not certify it.

`mixeff` computes these with four functions —
[`mm_grid()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
— that route all inference through the same machinery as
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md).
Each row in the returned table carries `method`, `status`,
`reliability`, and `reason` columns naming how it was computed and how
far to trust it (the column vocabularies are in the glossary linked
above).

### The study

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

``` r

rfit <- lmm(
  score ~ trt * time + (1 | clinic) + (1 | subj),
  rehab,
  control = mm_control(verbose = -1)
)
summary(rfit, tests = "coefficients")
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
control — a cell-specific quantity, estimand (1) territory. Summaries
over cells are the marginal verbs’ job.

### Reference grids with `mm_grid()`

[`mm_grid()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
constructs the cross-product of all fixed-predictor levels. By default,
factor predictors expand to all their levels; numeric predictors
collapse to their mean.

``` r

grd <- mm_grid(rfit, specs = ~ trt * time)
grd
#> Marginal grid:
#>      trt time
#>  control  pre
#>   active  pre
#>  control post
#>   active post
```

The grid has four rows — one for each treatment × timepoint cell — and
retains the model matrix needed for inference.

### Cell predictions with `mm_predictions()`

[`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
evaluates the fixed-effect prediction at each grid row, with a
confidence interval from the certified covariance. These are estimand
(1): the four population-level cell means.

``` r

preds <- mm_predictions(rfit, specs = ~ trt * time)
preds$table[, c("label", "estimate", "conf_low", "conf_high", "method")]
#>                    label estimate conf_low conf_high        method
#> 1  trt=control, time=pre 51.87148 50.08390  53.65906 satterthwaite
#> 2   trt=active, time=pre 49.89982 48.11224  51.68740 satterthwaite
#> 3 trt=control, time=post 47.90862 46.12103  49.69620 satterthwaite
#> 4  trt=active, time=post 43.11889 41.33131  44.90647 satterthwaite
```

Each row carries its inference method so the provenance is visible
without digging into model objects.

### Marginal means with `mm_means()`

Marginal means average the reference grid over dimensions you want to
collapse — estimand (2), the equal-weight reference-grid average. Here:
average over timepoints to get a treatment-level summary.

``` r

mt <- mm_means(rfit, specs = ~ trt)
mt$table[, c("label", "estimate", "conf_low", "conf_high", "method")]
#>         label estimate conf_low conf_high        method
#> 1 trt=control 49.89005 48.08653  51.69356 satterthwaite
#> 2  trt=active 46.50935 44.70584  48.31287 satterthwaite
```

Compare these to the raw means:

``` r

raw_means <- tapply(rehab$score, rehab$trt, mean)
raw_means
#>  control   active 
#> 49.45125 46.94815

marg_means <- setNames(mt$table$estimate, mt$table$label)
imbalance_pull <- max(abs(
  raw_means[c("control", "active")] -
    marg_means[c("trt=control", "trt=active")]
))
round(imbalance_pull, 2)
#> [1] 0.44
```

The clinic imbalance pulls each raw mean 0.44 points toward the clinics
its arm over-samples; the marginal means average over the model’s
fixed-effect grid instead and are free of that shift. The size of the
gap depends on how unbalanced the design is and how large the group
effects are.

### Pairwise comparisons with `mm_comparisons()`

[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
takes all pairwise differences among the marginal means and applies the
same inference method.

``` r

ctr <- mm_comparisons(rfit, specs = ~ trt)
ctr$table[, c("label", "estimate", "conf_low", "conf_high", "p_value", "method")]
#>                      label  estimate conf_low conf_high     p_value
#> 1 trt=active - trt=control -3.380694 -3.99461 -2.766778 3.86692e-10
#>          method
#> 1 satterthwaite
```

The `active - control` row is the equal-weight reference-grid contrast —
estimand (2) in the list above, averaged over both timepoints — tested
with a Satterthwaite *t* and graded `moderate` reliability in the same
row. In a randomized design you may argue it also estimates the causal
quantity (4); that argument rests on the design’s identification
assumptions, and it is yours to make — the fitted model does not certify
it.

### Conditional comparisons with `by =`

The `by` argument splits comparisons within levels of another variable —
the analogue of simple effects in a factorial design.

``` r

ctr_by <- mm_comparisons(rfit, specs = ~ trt | time)
ctr_by$table[, c("label", "estimate", "conf_low", "conf_high", "p_value", "method")]
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

### Constraining the grid with `at =`

`at` pins predictors to chosen values instead of expanding or averaging
them: a factor to a subset of its levels, a numeric predictor to
specific points rather than its mean.

``` r

mt_time <- mm_means(rfit, specs = ~ time, at = list(trt = "active"))
mt_time$table[, c("label", "estimate", "conf_low", "conf_high")]
#>       label estimate conf_low conf_high
#> 1  time=pre 49.89982 48.11224  51.68740
#> 2 time=post 43.11889 41.33131  44.90647
```

This gives the pre/post means *within the active arm* only, holding
`trt` constant at `"active"`.

### Custom contrasts with `mm_lincomb()`

For hypotheses that are not pairwise differences of marginal means,
build the contrast weights directly with
[`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md).
The interaction effect expressed as a contrast: (active post − active
pre) − (control post − control pre).

``` r

beta <- fixef(rfit)
names(beta)
#> [1] "(Intercept)"        "trtactive"          "timepost"          
#> [4] "trtactive:timepost"

# Interaction row: coefficient named "trtactive:timepost" (lme4-identical)
w <- setNames(numeric(length(beta)), names(beta))
w["trtactive:timepost"] <- 1

lc <- mm_lincomb(rfit, weights = w)
lc[, c("estimate", "lower", "upper", "p_value", "method")]
#>    estimate     lower     upper      p_value        method
#> 1 -2.818065 -3.296859 -2.339271 2.848244e-11 satterthwaite
```

[`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md)
applies the same contract-preserving inference as
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
so the method and status fields are populated identically. A row whose
status is anything other than `"available"` contains `NA` where the
refused numbers would go, with a stable `reason` code saying why; the
code vocabulary is in the glossary linked above. For the `emmeans`
bridge — and the audit columns that do *not* survive the trip through
`emmeans` — see
[`vignette("lme4-migration")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md).

## The covariance contract

The quantities above all lean on the fixed-effect covariance, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) is explicit about what it
returns. For full-rank fits it is the model-based covariance from the
engine’s `fixed_effect_covariance_matrix` payload, with method and
status attributes attached:

``` r

V <- vcov(fit)
V
#>                   (Intercept)          week treatmentcoached
#> (Intercept)       0.038596604 -1.683614e-03    -3.438757e-02
#> week             -0.001683614  6.734457e-04    -4.654064e-18
#> treatmentcoached -0.034387568 -4.654064e-18     6.877514e-02
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
#> [1] "1.1.0"
attr(V, "mm_method")
#> [1] "model_based"
attr(V, "mm_status")
#> [1] "available"
```

For rank-deficient or otherwise uncertified fits,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns an `NA` matrix
carrying an `mm_unavailable_reason` attribute instead of fabricated
numbers — the same available/unavailable/typed-error contract as every
other surface in this vignette. Prediction standard errors and intervals
(population versus conditional, and why the conditional SE can be
*smaller*) are documented in
[`vignette("lmm")`](https://bbuchsbaum.github.io/mixeff/articles/lmm.md).

## Summary

On an interior LMM fit,
[`summary()`](https://rdrr.io/r/base/summary.html),
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md),
[`anova()`](https://rdrr.io/r/stats/anova.html), and
[`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)
return labeled finite-sample inference, with `"auto"` resolving to
Satterthwaite.

On the boundary fit in this vignette:

- the model converges, but the random-effect covariance is reduced rank
  ([`vignette("convergence")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md)
  for what that label establishes);
- Wald inference is available but labeled as low reliability;
- Satterthwaite and Kenward–Roger routes are refused by mixeff’s
  validated-support policy: their degrees-of-freedom approximations
  assume interior variance parameters, which this fit does not have;
- parametric bootstrap testing is available for the fixed effect;
- bootstrap confidence intervals are available, but they use a
  full-model estimator-distribution target, not a fixed-effect null
  target;
- cluster bootstrap currently reports estimator distributions, not
  fixed-effect p-values;
- boundary LRTs are available for variance-component tests, not
  fixed-effect tests.

For GLMMs, three routes lead to reportable inference: certified joint
Wald (`method = "joint_laplace"`), the parametric bootstrap on profiled
fits, and the labeled working-Hessian opt-in for screening.

Every route either returns a labeled result or refuses with a stable
reason.

## References

- Crainiceanu, C. M., and Ruppert, D. (2004). Likelihood ratio tests in
  linear mixed models with one variance component. *Journal of the Royal
  Statistical Society: Series B*, 66(1), 165–185.
  <https://doi.org/10.1111/j.1467-9868.2004.00438.x>
- Halekoh, U., and Højsgaard, S. (2014). A Kenward-Roger approximation
  and parametric bootstrap methods for tests in linear mixed models —
  the `R` package `pbkrtest`. *Journal of Statistical Software*, 59(9),
  1–32. <https://doi.org/10.18637/jss.v059.i09>
- Kenward, M. G., and Roger, J. H. (1997). Small sample inference for
  fixed effects from restricted maximum likelihood. *Biometrics*, 53(3),
  983–997. <https://doi.org/10.2307/2533558>
- Self, S. G., and Liang, K.-Y. (1987). Asymptotic properties of maximum
  likelihood estimators and likelihood ratio tests under nonstandard
  conditions. *Journal of the American Statistical Association*,
  82(398), 605–610. <https://doi.org/10.1080/01621459.1987.10478472>
- Stram, D. O., and Lee, J. W. (1994). Variance components testing in
  the longitudinal mixed effects model. *Biometrics*, 50(4), 1171–1177.
  <https://doi.org/10.2307/2533455>
