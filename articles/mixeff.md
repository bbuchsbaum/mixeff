# Introduction and the Supported Model Contract

``` r

library(mixeff)
```

If you fit mixed-effects models in R, you most likely use
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).
[`mixeff::lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
speaks the same formula language, and the extractors —
[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`anova()`](https://rdrr.io/r/stats/anova.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`update()`](https://rdrr.io/r/stats/update.html) — do what you expect.
On the package’s supported envelope, statistical answers agree with
`lme4` within documented tolerances on the parity datasets shipped with
the package; the envelope itself is documented, closed, and
machine-readable
([`supported_models()`](https://bbuchsbaum.github.io/mixeff/reference/supported_models.md)).
It is not a literal *drop-in*: you call
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) rather
than `lmer()`, results are not bit-exact, and when the package changes
or declines a request it says so, in a diagnostic you can read from the
fitted object.

The reason to switch is what `mixeff` does around the fit. On the
scaling benchmark shipped with the package it ran faster than `lme4` at
the largest tested scale of each of five common designs — all small fits
in absolute terms; the benchmarking article on the package site,
<https://bbuchsbaum.github.io/mixeff/articles/benchmarking.html>,
reports the numbers and what they do and do not establish. And it makes
four things explicit that `lme4` leaves implicit:

1.  **The formula stays familiar.** Anything you would hand to `lmer()`
    you can hand to
    [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md).
2.  **Singular fits become labeled facts.** A reduced-rank random-effect
    covariance is reported on the fitted object, not as a single console
    warning.
3.  **Inference labels its method.** Each p-value carries the method
    that produced it. Where the package declines to compute
    Satterthwaite or Kenward–Roger degrees of freedom — at a boundary
    fit, for example — a parametric bootstrap is offered under the same
    labeling rules.
4.  **The fitted object is a record.** A saved model reopens with the
    same coefficients, the same diagnostics, and the same method labels,
    without depending on the original Rust handle.

Behind the four is one design habit. In ordinary mixed-model work, three
questions come up again and again: what does this random-effects formula
actually mean; are the p-values or tests available, and by what method;
can the model state be reconstructed after the fit has been saved.
`mixeff` makes those questions part of the fitted object.

This page demonstrates each of the four on one small dataset, then
states the supported-model contract and the persistence guarantees. The
other vignettes go deeper on every step.

## The dataset

A small repeated-measures study: 18 subjects, 10 daily reaction-time
measurements each, intercepts and slopes that are nearly perfectly
correlated by construction. This design produces a singular fit in
`lme4` and in `mixeff` alike.

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

## B. When a fit is degenerate, the object says so

By construction, this design pushes the random-effect covariance to
reduced rank. `mixeff` fits the model and records the outcome on the
object, where you can query it later:

``` r

fit_status(fit)
#> [1] "converged_reduced_rank"
is_singular(fit)
#> [1] TRUE
```

What this status means, which variance component reached the boundary,
the diagnostic codes behind the label
([`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)),
the requested-to-effective-to-fitted record
([`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)),
and what a boundary does to downstream inference are the subject of
[`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md).

## C. Inference labels its method

On the inference side, `mixeff` first reports which methods are
available on your fit, and then lets you run the one you choose under
the same labeling rules.

[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
enumerates the inference methods for the current fit, gives each one a
fixed-vocabulary *reason* for its status (the vocabularies are defined
in the [inference-method
glossary](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.html)
on the package site), and names the function you would call to invoke
it.

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

Some routes are immediately available on this fit; others decline with a
stable reason code rather than printing a number without comment. Why a
boundary fit takes some methods off the table is explained in
[`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md);
the routes themselves — including the parametric bootstrap that yields a
defensible p-value on a fit like this one — are walked through in
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).

Every reported quantity is in one of three states: available, with a
named method and reason; unavailable, with a stable reason code; or a
typed error.

## The supported model contract

The envelope is machine-readable.
[`supported_models()`](https://bbuchsbaum.github.io/mixeff/reference/supported_models.md)
returns the same registry the fitting code enforces
(`inst/support/model-support.json`); the test suite holds the fitting
code, the registry, and the documentation to the same envelope.

``` r

sm <- supported_models()
knitr::kable(sm$glmm_families[, c("family", "links", "estimators")],
             caption = "GLMM family/link/estimator cells, from the registry.")
```

| family            | links                  | estimators                    |
|:------------------|:-----------------------|:------------------------------|
| binomial          | logit, probit, cloglog | pirls_profiled, joint_laplace |
| poisson           | log, sqrt              | pirls_profiled, joint_laplace |
| Gamma             | log                    | pirls_profiled, joint_laplace |
| negative_binomial | log                    | pirls_profiled                |

GLMM family/link/estimator cells, from the registry. {.table}

``` r

knitr::kable(sm$features[, c("feature", "lmm", "glmm")],
             caption = "Feature support, from the registry.")
```

| feature | lmm | glmm |
|:---|:---|:---|
| case weights | supported | supported |
| offset | unsupported | supported |
| simulate | supported | refused |
| refit | supported | refused |
| subset / random / custom na.action / custom contrasts (GLMM) | see notes | refused |
| marginal means verbs (mm_means / mm_comparisons / mm_grid / mm_predictions / test_effect) | supported | refused |
| double-bar \|\| with factor terms | mixeff semantics | mixeff semantics |

Feature support, from the registry. {.table}

Explicit non-goals, from the same registry: nonlinear mixed models;
arbitrary GLM family objects or custom links; adaptive Gaussian
quadrature as a general estimator contract (nAGQ \> 1 is a profiled-path
sensitivity mode).

These are the short forms.
[`vignette("lme4-migration", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md)
renders the full registry tables with their notes columns and is the
verb-for-verb, argument-for-argument map for porting an lme4 script.

Two boundary conditions of the contract are worth stating up front:

- Estimates are not bit-for-bit identical to `lme4`; the bar is
  statistical agreement within documented tolerances on the parity
  datasets shipped with the package. Where results are *expected* to
  differ, the differences are classified in
  `inst/extdata/expected-mismatches.json`, with tolerances enforced by
  the test suite.
- Generalized fits carry their own estimator and inference contract; see
  [`vignette("glmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md)
  for the estimators and
  [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
  for the inference routes.

Each limit is reported by the package itself, as a reason code or a
typed error, when you hit it.

## D. The fit is the record

The fitted object is a serializable record.
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) followed by
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) and
[`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
reproduces the audit trail and the extractors without depending on the
original Rust handle.

``` r

path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
restored_fit <- revive(readRDS(path))

identical(fixef(restored_fit),             fixef(fit))
#> [1] TRUE
identical(changes(restored_fit)$table,     changes(fit)$table)
#> [1] TRUE
identical(diagnostics(restored_fit)$table, diagnostics(fit)$table)
#> [1] TRUE
```

A reviewer reading the `.rds` later sees the same convergence status,
the same reduced-rank diagnostic, the same method labels on the same
coefficients.

### Saving, reloading, and reviving in practice

A mixed-model fit often outlives the R session that produced it: you fit
a model today, save it with the analysis, hand the script to a
collaborator, and reopen the fit six months later for a contrast, a
revision, or a referee response. `mixeff` stores the fitted values, the
random-effects design, the convergence record, and the inference labels
inside the R object, so each of those tasks works after
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) without recomputing
the fit.

The round trip above used the singular fit; the same guarantees hold for
an ordinary interior fit.

``` r

fit2 <- lmm(
  score ~ week + treatment + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)

path2 <- tempfile(fileext = ".rds")
saveRDS(fit2, path2)

restored <- readRDS(path2)
restored <- revive(restored)
```

[`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
resets the object’s process-local cache. It is cheap and safe to call
after every [`readRDS()`](https://rdrr.io/r/base/readRDS.html), but the
extractors and inference functions shown below also work on the plain
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) result — the durable
values live in the object itself, and the Rust handle is deliberately
absent after reload.

``` r

fixef(restored)
#>      (Intercept)             week treatmentcoached 
#>        7.6828778       -0.2783994       -0.8994747
head(predict(restored))
#>        1        2        3        4        5        6 
#> 7.585932 7.307533 7.029134 6.750734 6.472335 6.193935
reporting_table(restored, "fixed_effects")
#>              term   estimate  std_error        df  statistic statistic_name
#>       (Intercept)  7.6828778 0.19646018 12.565022  39.106539              t
#>              week -0.2783994 0.02595083 58.999740 -10.727955              t
#>  treatmentcoached -0.8994747 0.26225014  9.999273  -3.429835              t
#>       p_value        method    status reliability
#>  1.665335e-14 satterthwaite available    moderate
#>  1.776357e-15 satterthwaite available    moderate
#>  6.440943e-03 satterthwaite available    moderate
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

Two related surfaces are documented elsewhere: the honest-covariance
contract behind [`vcov()`](https://rdrr.io/r/stats/vcov.html) is part of
the inference contract
([`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)),
and rebuilding design matrices
([`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
[`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md))
from the stored formula and model frame is shown in
[`vignette("lmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm.md).

### The artifact contract

The structured JSON artifacts stored in the fitted object are the source
of truth; external pointers such as the Rust handle are caches.
[`mm_json_known_schemas()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_known_schemas.md)
lists the artifact schemas this version of the package understands:

``` r

head(mm_json_known_schemas())
#>                                       name version
#> 1                                  formula      v0
#> 2      mixedmodels.compiled_model_artifact       1
#> 3           mixedmodels.model_audit_report       2
#> 4             mixedmodels.random_term_card       1
#> 5 mixedmodels.fixed_effect_inference_table   1.1.0
#> 6      mixedmodels.marginal_quantity_table   1.0.0
```

## Lower-level tools

Most users start with
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md),
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md),
[`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md),
and
[`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md).
The lower-level functions are there when you need them:

- [`mm_parse_formula()`](https://bbuchsbaum.github.io/mixeff/reference/mm_parse_formula.md)
  checks formula syntax and returns the canonical spelling (demonstrated
  in
  [`vignette("lmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm.md)).
- [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  builds a pre-fit model specification;
  [`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
  translates it into named forms and plain-language scopes (also in
  [`vignette("lmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm.md)).
- [`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
  and
  [`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
  expose model-state checks.
- [`mm_json_known_schemas()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_known_schemas.md)
  lists the structured artifact schemas understood by this version of
  the package.

The Rust backend does not change how you write R: you fit the model, get
the numbers, and the status of those numbers stays attached to the
object.

## Where to read next

- [`vignette("lmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm.md)
  — linear mixed models end-to-end: reading random-effects formulas,
  fitting, summaries, extractors, weights, and prediction.
- [`vignette("glmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md)
  — the generalized families and estimators.
- [`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md)
  — convergence statuses, boundary and singular fits, and verification.
- [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
  — p-values, contrasts, term tests, bootstrap methods, confidence
  intervals, and marginal means.
- [`vignette("lme4-migration", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md)
  — the verb-for-verb map and the generated envelope tables.
- On the package site: the [benchmarking
  article](https://bbuchsbaum.github.io/mixeff/articles/benchmarking.html),
  the [inference-method
  glossary](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.html),
  [reporting
  workflows](https://bbuchsbaum.github.io/mixeff/articles/reporting-lmms.html),
  and a [full published-study
  reproduction](https://bbuchsbaum.github.io/mixeff/articles/reproducing-aphantasia.html).
