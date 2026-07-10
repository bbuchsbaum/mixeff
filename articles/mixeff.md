# Why mixeff?

``` r

library(mixeff)
```

If you fit mixed-effects models in R, you most likely use
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).
[`mixeff::lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
aims to be *functionally equivalent*: the formula language is the same,
and the extractors —
[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`anova()`](https://rdrr.io/r/stats/anova.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`update()`](https://rdrr.io/r/stats/update.html) — do what you expect.
Statistical answers agree with `lme4` within documented tolerances on
the parity datasets shipped with the package. It is not a literal
*drop-in*: you call
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) rather
than `lmer()`, results are not bit-exact, and the package is audit-first
— it reports or refuses rather than silently transforming a model.

The reason to switch is what `mixeff` does around the fit. It is faster
on most designs — typically by a factor of two to five — and it makes
four things explicit that `lme4` leaves implicit:

1.  **The formula stays familiar.** Anything you would hand to `lmer()`
    you can hand to
    [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md).
2.  **Singular fits become labelled facts.** A reduced-rank
    random-effect covariance is reported with codes, severity, and
    effective rank, instead of a single warning that scrolls off the
    screen.
3.  **Inference labels what asymptotics can and cannot do.** Each
    p-value carries the method that produced it; where Wald,
    Satterthwaite, and Kenward-Roger are unavailable in principle, a
    parametric bootstrap is offered as a labelled first-class peer.
4.  **The fitted object is a record.** A model saved today reopens six
    months later with the same coefficients, the same audit trail, and
    the same method labels — without depending on the original Rust
    handle.

This page demonstrates each of the four on one small dataset.

## The dataset

A small repeated-measures study: 18 subjects, 10 daily reaction-time
measurements each, intercepts and slopes that are nearly perfectly
correlated by construction. This is the kind of design that produces a
singular fit in any modern engine.

``` r

set.seed(3)
n_subj <- 18L
days   <- 0:9
b0     <- rnorm(n_subj, sd = 30)
b1     <- 0.5 * b0 / 30 * 10 + rnorm(n_subj, sd = 0.5)

sleep_like <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
  data.frame(
    subj = factor(i),
    days = days,
    rt   = 250 + b0[i] + (10 + b1[i]) * days +
           rnorm(length(days), sd = 20)
  )
}))

head(sleep_like, 4)
#>   subj days       rt
#> 1    1    0 247.1492
#> 2    1    1 227.7095
#> 3    1    2 213.1613
#> 4    1    3 254.4247
```

## A. The formula stays familiar

`mixeff` keeps the lme4 random-effects syntax. If you can read
`(x | g)`, you can read `mixeff`.

``` r

fit <- lmm(
  rt ~ days + (1 + days | subj),
  sleep_like,
  control = mm_control(verbose = -1)
)
fit
#> Linear mixed model fit by REML
#> Formula: rt ~ days + (1 + days | subj)
#> Fit status: converged_reduced_rank
#> Optimizer: trust_bq; iterations: 314; objective: 1647.7
#> nobs: 180, sigma: 20.065, logLik: -823.849
#> Fixed effects:
#> (Intercept)        days 
#>   239.87900     9.25095 
#> 
#> Fitted covariance state:
#> The fitted covariance matrix is rank-deficient.
#>   r0: requested rank 2; fitted effective rank 1.
#> Use changes(fit) to see which dimension was unsupported.
#> Use random_options(spec, group = subj) to inspect lower-dimensional covariance choices.
#> Audit verbs: audit(), diagnostics(), inference_table(), model_report()
```

[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`anova()`](https://rdrr.io/r/stats/anova.html), and
[`summary()`](https://rdrr.io/r/base/summary.html) all do what you
expect.

## B. When a fit is degenerate, you find out *which* part

[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) fits this model
and reports the situation in one parsimonious line:
`boundary (singular) fit`. The fact is correct. What it leaves implicit
is *which* variance component reached the boundary, what the effective
rank of the random-effect covariance is, and which downstream inference
methods are no longer defined.

``` r

m <- suppressMessages(lme4::lmer(
  rt ~ days + (1 + days | subj),
  data = sleep_like
))
m
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: rt ~ days + (1 + days | subj)
#>    Data: sleep_like
#> REML criterion at convergence: 1647.697
#> Random effects:
#>  Groups   Name        Std.Dev. Corr 
#>  subj     (Intercept) 23.682        
#>           days         3.221   1.00 
#>  Residual             20.065        
#> Number of obs: 180, groups:  subj, 18
#> Fixed Effects:
#> (Intercept)         days  
#>     239.879        9.251  
#> optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings
lme4::isSingular(m)
#> [1] TRUE
```

`mixeff` reports the same fact, and then unpacks it.
[`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
names the convergence outcome,
[`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
shows the requested-to-effective transition, and
[`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
returns stable codes.

``` r

fit_status(fit)
#> [1] "converged_reduced_rank"
is_singular(fit)
#> [1] TRUE
changes(fit)
#> Model changes:
#>   Fitted covariance for (1 + days | subj): requested rank 2, fitted rank 1 [reduced_rank].
#> Stage-by-stage records available via $table.
diagnostics(fit)$table[, c("code", "severity", "stage", "message")]
#>                 code severity         stage
#> 1 boundary_parameter     info certification
#> 2 covariance_reduced     info certification
#>                                                                            message
#> 1           standard deviation for days in (1 + days | subj) is on its lower bound
#> 2 fitted covariance for (1 + days | subj) has effective rank 1 of requested rank 2
```

A reduced-rank covariance is now a *labelled* fact about the fit, not a
warning that scrolls off the screen.

## C. Expose when asymptotics are weak, then offer bootstrap

`mixeff` has two jobs on the inference side. First, it tells you when
the asymptotic methods (Wald z, Satterthwaite, Kenward-Roger) are
unreliable on your fit. Then, when bootstrap is the defensible route, it
offers it as a labelled, first-class peer to the asymptotic methods —
not a fallback trick.

[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
is the audit verb for that judgment. It enumerates the inference methods
available on the current fit, gives each one a [closed-enum
*reason*](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.md)
for its status, and names the verb you would call to invoke it.

``` r

opt <- inference_options(fit, "days", nsim = 200)
opt$table[, c("method", "expected_status",
              "expected_reliability_reason", "current")]
#>              method expected_status
#> 1 asymptotic_wald_z       available
#> 2     satterthwaite    not_assessed
#> 3     kenward_roger    not_assessed
#> 4         bootstrap       available
#> 5     bootstrap_lrt    not_assessed
#> 6 cluster_bootstrap    not_assessed
#> 7        profile_ci    not_assessed
#>                      expected_reliability_reason current
#> 1                     asymptotic_wald_z_fallback    TRUE
#> 2          satterthwaite_unavailable_at_boundary   FALSE
#> 3          kenward_roger_unavailable_at_boundary   FALSE
#> 4               bootstrap_monte_carlo_replicates   FALSE
#> 5                      bootstrap_lrt_requires_ml   FALSE
#> 6 bootstrap_cluster_resample_p_value_unavailable   FALSE
#> 7             profile_ci_unavailable_at_boundary   FALSE
```

Two routes are available on this fit: asymptotic Wald z (immediate, but
labelled `low` reliability), and bootstrap (~seconds, labelled by
replicate count and Monte-Carlo SE). Satterthwaite and Kenward-Roger
refuse with a stable reason — `*_unavailable_at_boundary` — because at a
boundary fit the variance-parameter derivative they need is not defined.
That is a fact about the math, not a bug.

The asymptotic Wald row carries its own warrant.
[`summary()`](https://rdrr.io/r/base/summary.html) now prints
`reliability_reason` next to `reliability`:

``` r

inf <- inference_table(fit)$table
inf[, c("term", "method", "status", "reliability", "reliability_reason")]
#>          term            method    status reliability
#> 1 (Intercept) asymptotic_wald_z available         low
#> 2        days asymptotic_wald_z available         low
#>           reliability_reason
#> 1 asymptotic_wald_z_fallback
#> 2 asymptotic_wald_z_fallback
```

`degrees_of_freedom_unavailable_so_z_substituted` is the closed-enum
warrant: a t reference distribution was the requested target but the df
could not be computed at this boundary fit, so a standard normal was
substituted. The number is real; the *grade* is calibrated.

For a defensible p-value on this same fit, route through
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
with `method = "bootstrap"`. The Rust engine simulates from the
constrained null, refits each replicate, and returns a labelled
inference row plus a run payload (boundary rate, MCSE, replicate count)
for audit.

``` r

ct <- contrast(fit, c(0, 1), method = "bootstrap",
               bootstrap = bootstrap_control(nsim = 200, seed = 1))
ct$table[, c("contrast", "estimate", "p_value",
             "method", "status", "reliability")]
#>   contrast estimate     p_value    method    status reliability
#> 1       c1 9.250949 0.004975124 bootstrap available         low

run <- ct$table$details[[1]]$bootstrap
data.frame(
  successful_replicates = run$successful_replicates,
  boundary_rate         = round(run$boundary_rate, 3),
  mcse                  = round(run$mcse, 4)
)
#>   successful_replicates boundary_rate  mcse
#> 1                   200          0.46 0.005
```

The bootstrap p-value is `available`, the method is named, the run
payload makes the simulation provenance explicit. The boundary rate is
visible because singular fits propagate boundary behaviour into their
own bootstrap replicates — that is honest, not a bug. `mcse` quantifies
the Monte-Carlo uncertainty of the p-value estimate; raise `nsim` for a
tighter MCSE.

There are three states for any reported quantity: *available with a
named method and a closed-enum warrant*, *unavailable with a stable
reason code*, or *typed error*. There is no fourth state where the
package guesses.

## D. The fit is the record

The fitted object is a serialisable record.
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) followed by
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) and
[`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
reproduces the audit trail and the extractors without depending on the
original Rust handle.

``` r

path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
restored <- revive(readRDS(path))

identical(fixef(restored),               fixef(fit))
#> [1] TRUE
identical(changes(restored)$table,       changes(fit)$table)
#> [1] TRUE
identical(diagnostics(restored)$table,   diagnostics(fit)$table)
#> [1] TRUE
```

A reviewer reading the `.rds` six months from now sees the same
convergence status, the same reduced-rank diagnostic, the same method
labels on the same coefficients.

## What this page did not show

`mixeff` does not (yet, by design):

- match `lme4` numerics bit-for-bit; statistical equivalence within
  documented tolerances on parity datasets is the bar.
- provide the joint-Laplace / AGQ GLMM backend;
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  currently ships the labelled profiled-PIRLS path and refuses
  unavailable joint methods explicitly
  ([`vignette("glmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md)).
- ship Kenward-Roger or profile-likelihood confidence intervals in v0.

Each of those is a stable boundary with a name, not a missing feature
hidden behind a fallback.

## Where to read next

- [`vignette("lmm-basics", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm-basics.md)
  — fitting and the standard extractors at a slower pace.
- [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
  — coefficient tests, contrasts, term tests, and model comparisons.
- [`vignette("demystifying-formulas", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/demystifying-formulas.md)
  — what `(1 | g)`, `(x | g)`, split blocks, and `||` actually mean.
- [`vignette("saving-and-reviving", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/saving-and-reviving.md)
  — the round-trip story in detail.
