# Inference Where Standard Mixed-Model p-values Break Down

``` r

library(mixeff)
```

Mixed models sometimes converge but still leave us in a difficult
inference situation. The fit exists, the estimates are usable, and the
model may be scientifically reasonable, but some familiar p-values are
no longer well justified.

This vignette focuses on one common case: boundary or singular fits. In
a mixed model, this can happen when a variance component is estimated
very close to zero, or when random-effect terms are estimated as
perfectly or nearly perfectly correlated. Equivalently, the fitted
random-effect covariance matrix has lower effective rank than the model
requested.

These are exactly the situations where ordinary large-sample tests can
become least trustworthy. Wald tests, Satterthwaite approximations, and
Kenward-Roger adjustments all rely on regularity conditions that are
strained or broken when a variance parameter is on the boundary of its
parameter space.

The statistical reason is well known. Self and Liang (1987) showed that
likelihood-ratio statistics do not necessarily follow their usual
chi-square reference distributions when a parameter lies on the
boundary. Stram and Lee (1994) specialized this issue to variance
components in linear mixed models. Kenward and Roger (1997) derived
their fixed-effect adjustment from a Taylor expansion of the variance
estimator; that derivation assumes a regular interior solution for the
variance parameters. At a boundary, those derivatives may not exist or
may not behave regularly.

`mixeff` handles this by making the inference contract explicit. First,
it does not silently present every number as equally reliable.
Asymptotic rows are labelled with a [closed-enum reliability
reason](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.md),
such as `asymptotic_wald_z_fallback` or `*_unavailable_at_boundary`.
Second, when bootstrap inference is the better-supported route, `mixeff`
exposes it as a labelled method rather than an informal workaround.
Parametric bootstrap methods have a long history in this setting; for
example, Halekoh and Højsgaard (2014) implemented parametric bootstrap
tests in `pbkrtest` for small-sample and boundary cases where
Kenward-Roger methods may need support.

The goal is not to say that a singular fit is automatically useless. The
goal is to keep the covariance state attached to the inference result,
so that available numbers and unavailable numbers are both documented.

## What fit are we worried about?

This small repeated-measures example has subject-specific intercepts and
slopes. In this data set, those two random effects are nearly collinear,
which is enough to produce a reduced-rank random-effect covariance
estimate.

``` r

fit <- lmm(
  rt ~ days + (1 + days | subj),
  sleep_like,
  control = mm_control(verbose = -1)
)

fit_status(fit)
#> [1] "converged_reduced_rank"
is_singular(fit)
#> [1] TRUE
```

The fit converged, but the covariance estimate is singular. That matters
for inference, so `mixeff` records it in the model diagnostics.

``` r

diagnostics(fit)$table[, c("code", "severity", "stage", "message")]
#>                 code severity         stage
#> 1 boundary_parameter     info certification
#> 2 covariance_reduced     info certification
#>                                                                            message
#> 1           standard deviation for days in (1 + days | subj) is on its lower bound
#> 2 fitted covariance for (1 + days | subj) has effective rank 1 of requested rank 2
```

The requested random-effect covariance for `(1 + days | subj)` has rank
2: one dimension for the subject intercept and one for the subject
slope. The fitted covariance has effective rank 1. In plain language,
the model asked for two random-effect dimensions, but the data support
only one effective dimension.

That does not automatically invalidate the model. It does mean that
inference methods relying on a regular full-rank covariance estimate
need to be treated carefully.

## Which inference routes are available?

Before testing a term, ask which inference routes are defined for this
fit.
[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
is an audit map, not a recommender. Each row describes one method that
`mixeff` knows about. The row says whether the method can run, why it
can or cannot run, and what the rough computational cost is.

The printed object is the reader-facing map. It keeps the raw enum
columns in `routes$table` for scripts, but explains the route status
with display columns. The complete contract vocabulary is in the
[inference method
glossary](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.md).

``` r

routes <- inference_options(fit, "days", nsim = 200)
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
#>          approx_cost current
#>            immediate    TRUE
#>            immediate   FALSE
#>            immediate   FALSE
#>       ~4s @ nsim=200   FALSE
#>       ~8s @ nsim=200   FALSE
#>                    -   FALSE
#>  slow profile refits   FALSE
#> 
#> Use `<obj>$table` for raw enum columns (`expected_status`, `expected_reliability_reason`) and notes.
```

Read this table as follows. The Wald route can run immediately, but it
is labelled as a low-reliability asymptotic fallback. Satterthwaite and
Kenward-Roger are refused because the variance-parameter derivatives
needed by those approximations are not available at the boundary. The
parametric bootstrap route can run and is the documented fixed-effect
testing route for this fit. The bootstrap LRT route is refused here
because this model was fit by REML, while that route requires an ML fit.
The cluster bootstrap route is refused for fixed-effect p-values because
its current target is an estimator distribution, not a null hypothesis
test. Profile confidence intervals are not certified for this boundary
fit.

``` r

inference_table(fit)$table[, c("term", "method", "status",
                               "reliability", "reliability_reason")]
#>          term            method    status reliability
#> 1 (Intercept) asymptotic_wald_z available         low
#> 2        days asymptotic_wald_z available         low
#>           reliability_reason
#> 1 asymptotic_wald_z_fallback
#> 2 asymptotic_wald_z_fallback
```

The important point is that `mixeff` distinguishes three cases:

1.  A method can run and is considered available.
2.  A method can produce a number, but the number is labelled as low
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

The small `nsim` value below keeps the vignette fast. For a real
analysis, use more bootstrap replicates.

``` r

term_boot <- test_effect(
  fit,
  "days",
  method = "bootstrap",
  bootstrap = bootstrap_control(nsim = 50, seed = 1)
)

term_boot$table[, c("term", "statistic_name", "p_value",
                    "method", "status", "reliability_reason")]
#>   term statistic_name    p_value    method    status
#> 1 days              t 0.01960784 bootstrap available
#>                 reliability_reason
#> 1 parametric_bootstrap_monte_carlo
```

The p-value is returned as part of the labelled result row. It is not
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
`1 / (50 + 1)`. Use more replicates when the p-value itself will be
reported or compared with a decision threshold.

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
  fit,
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
  fit,
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
bootstrap result. Available numbers are labelled with their method and
reliability reason. Unavailable numbers are labelled with stable reason
codes explaining why they are unavailable.

## Variance-component boundary tests

The boundary likelihood-ratio route is for random-effect variance
components, not ordinary fixed effects.

This is a different testing problem. For a variance component, the null
hypothesis can put the parameter exactly on the boundary of the
parameter space, because variances cannot be negative. In the simplest
one-component case, the reference distribution is the Self-Liang 50:50
mixture: `0.5 * chi-square(0) + 0.5 * chi-square(1)`.

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

This is not a software gap. It is a methodological distinction. Boundary
likelihood-ratio theory applies naturally to variance-component tests,
where the null hypothesis places a variance on the edge of the parameter
space. It is not the reference distribution for an ordinary fixed-effect
Wald row.

## Summary

A singular mixed-model fit can still be useful, but it changes the
inference problem.

For this fit:

- the model converges, but the random-effect covariance is reduced rank;
- Wald inference is available but labelled as low reliability;
- Satterthwaite and Kenward-Roger routes are refused because the
  variance-parameter derivatives are not defined at the boundary;
- parametric bootstrap testing is available for the fixed effect;
- bootstrap confidence intervals are available, but they use a
  full-model estimator-distribution target, not a fixed-effect null
  target;
- cluster bootstrap currently reports estimator distributions, not
  fixed-effect p-values;
- boundary LRTs are available for variance-component tests, not
  fixed-effect tests.

The main design principle is that `mixeff` does not make the analyst
guess. Every route either returns a labelled result or refuses with a
stable reason.

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
