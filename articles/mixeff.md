# Introduction to mixeff

``` r

library(mixeff)
```

`mixeff` fits linear and generalized linear mixed models in R using the
`mixeff-rs` Rust engine. If you use `lme4`, the formula syntax and
functions such as
[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`predict()`](https://rdrr.io/r/stats/predict.html), and
[`summary()`](https://rdrr.io/r/base/summary.html) will be familiar.
Supported models are listed by
[`supported_models()`](https://bbuchsbaum.github.io/mixeff/reference/supported_models.md);
differences in arguments and output are described in
[`vignette("lme4-migration")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md).

The fitted object stores the model specification, convergence
diagnostics, and the method used for each inference result. You can
inspect a singular fit, find out which tests are available, and save the
results for later use. This vignette demonstrates those steps on a small
repeated-measures dataset, then explains model support and the
requirements for restoring a saved fit.

Tests compare `mixeff` with `lme4` on the datasets shipped with the
package, using documented numerical tolerances. The [benchmarking
article](https://bbuchsbaum.github.io/mixeff/articles/benchmarking.html)
reports timing comparisons and simulation results, with the package
versions and model designs used in each study.

## The dataset

We simulate 10 daily reaction-time measurements for each of 18 subjects.
The subject intercepts and slopes are nearly perfectly correlated by
construction. This example produces a singular fit in both `lme4` and
`mixeff`.

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

## Fitting a model

The formula below estimates a population slope for `days` and allows
each subject to have a different intercept and slope.

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
#> Optimizer: trust_bq; iterations: 30; objective: 1647.7
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
[`summary()`](https://rdrr.io/r/base/summary.html) all work on this fit;
[`vignette("lme4-migration")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md)
lists where their output differs from lme4’s.

## Singular fits

In this example, the fitted random-effect covariance has reduced rank.
Use
[`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
and
[`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md)
to inspect the result:

``` r

fit_status(fit)
#> [1] "converged_reduced_rank"
is_singular(fit)
#> [1] TRUE
```

[`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
gives more detail, and
[`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
compares the requested model with its effective specification and fitted
covariance. See
[`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md)
for the meaning of each status and the consequences for inference.

## Available inference methods

[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
lists the methods available for a fitted model, explains any
restrictions, and gives the function call for each method. The status
and reason codes are defined in the [inference-method
glossary](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.html).

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

For this singular fit, `mixeff` withholds Satterthwaite and
Kenward–Roger tests. Parametric bootstrap tests remain available. The
inference vignette explains these restrictions and demonstrates the
alternatives:
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).

An inference result includes its method and status. If a requested
quantity is unavailable, the function returns a reason code with the
result or raises a typed error.

## Supported models

[`supported_models()`](https://bbuchsbaum.github.io/mixeff/reference/supported_models.md)
reads the package’s support registry
(`inst/support/model-support.json`). Tests check that the registry,
fitting code, and documentation agree.

``` r

sm <- supported_models()
knitr::kable(sm$glmm_families[, c("family", "links", "estimators")],
             caption = "GLMM family/link/estimator cells, from the registry.")
```

| family            | links                  | estimators                    |
|:------------------|:-----------------------|:------------------------------|
| binomial          | logit, probit, cloglog | pirls_profiled, joint_laplace |
| poisson           | log, sqrt              | pirls_profiled, joint_laplace |
| Gamma             | inverse, log           | pirls_profiled, joint_laplace |
| inverse.gaussian  | inverse, log           | pirls_profiled, joint_laplace |
| gaussian          | log, inverse, sqrt     | pirls_profiled, joint_laplace |
| negative_binomial | log                    | pirls_profiled                |

GLMM family/link/estimator cells, from the registry. {.table}

``` r

knitr::kable(sm$features[, c("feature", "lmm", "glmm")],
             caption = "Feature support, from the registry.")
```

| feature | lmm | glmm |
|:---|:---|:---|
| case weights | supported | supported |
| offset | supported | supported |
| simulate | supported | supported |
| refit | supported | supported |
| subset / na.action / contrasts | supported | supported |
| random = (GLMM) | not applicable | refused |
| stateful fixed-effect terms (factor(), poly(), scale(), ns()/bs(), cut(), relevel(), offset(), ^, %in%, parenthesised groups) | supported | supported |
| marginal means verbs (mm_means / mm_comparisons / mm_grid / mm_predictions / test_effect) | supported | refused |
| double-bar \|\| with factor terms | supported | supported |

Feature support, from the registry. {.table}

mixeff does not fit nonlinear mixed models. Only the GLM families and
links listed above are supported. The nAGQ \> 1 option lets you check
how quadrature order affects a profiled-PIRLS fit. It requires a single
scalar random effect and is unavailable with the joint-Laplace
estimator.

For the full tables, including notes on individual features, see
[`vignette("lme4-migration", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md).
When comparing results with lme4, keep two qualifications in mind:

- Estimates are not bit-for-bit identical to `lme4`; the criterion is
  agreement within documented tolerances on the parity datasets shipped
  with the package. Where results are *expected* to differ, the
  differences are classified in `inst/extdata/expected-mismatches.json`,
  with tolerances enforced by the test suite.
- Generalized fits have estimator-specific inference methods; see
  [`vignette("glmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md)
  for the estimators and
  [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
  for the inference routes.

## Saving and restoring a fit

Save a fit with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
reload it with [`readRDS()`](https://rdrr.io/r/base/readRDS.html). The
stored coefficients and diagnostics remain available after reloading:

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

### Operations on a restored fit

Stored coefficients, fitted values, variance components, and audit
records are available after
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) without a live Rust
handle. The checks above compare these stored values before and after
saving.

New LMM contrasts and term tests, conditional random-effect variances,
engine-backed prediction, and prediction uncertainty reuse a live native
handle. After serialization, native queries restore the recorded fitted
state without running an optimizer.
[`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
resets the R cache and clears a dead handle; restoration rebuilds matrix
factors at the saved parameters. Derivative calculations required by an
inference method may still take time.

Restoration requires a matching engine source and build configuration. A
fit without a live handle and without a compatible snapshot keeps its
stored extractors but refuses native queries that require restoration.
Keep the package and engine build with the saved fit when reproducing an
analysis. Explicit refitting, bootstrap, profile intervals and
convergence verification still perform the optimization those operations
request.

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
resets the object’s process-local cache. You can call it after
[`readRDS()`](https://rdrr.io/r/base/readRDS.html), but the extractors
and inference functions shown below also work on the plain
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) result. Stored
extractors read durable values; operations requiring a native model
restore the saved state as described above when no live handle remains.

``` r

fixef(restored)
#>      (Intercept)             week treatmentcoached 
#>        7.6828778       -0.2783994       -0.8994747
head(predict(restored))
#>        1        2        3        4        5        6 
#> 7.585932 7.307533 7.029134 6.750734 6.472335 6.193935
reporting_table(restored, "fixed_effects")
#>              term   estimate  std_error        df  statistic statistic_name
#>       (Intercept)  7.6828778 0.19646018 12.565500  39.106539              t
#>              week -0.2783994 0.02595083 59.000150 -10.727955              t
#>  treatmentcoached -0.8994747 0.26225014  9.999671  -3.429835              t
#>       p_value        method    status reliability
#>  1.665335e-14 satterthwaite available    moderate
#>  1.776357e-15 satterthwaite available    moderate
#>  6.440574e-03 satterthwaite available    moderate
```

Conditional variances for random effects can also be requested after
reload. The engine restores the saved fit without optimization. With
`condVar = TRUE`, each grouping table in this example carries a finite
`postVar` array.

``` r

re <- ranef(restored, condVar = TRUE)
attr(re, "mm_unavailable_reason")
#> NULL
dim(attr(re$clinic, "postVar"))
#> [1]  1  1 12
```

For details of [`vcov()`](https://rdrr.io/r/stats/vcov.html), see
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).
Rebuilding design matrices with
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) and
[`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md) is
demonstrated in
[`vignette("lmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm.md).

### Stored records and schemas

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
#> 5 mixedmodels.fixed_effect_inference_table   1.2.0
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
  — corresponding functions, argument differences, and support tables.
- On the package site: the [benchmarking
  article](https://bbuchsbaum.github.io/mixeff/articles/benchmarking.html),
  the [inference-method
  glossary](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.html),
  [reporting
  workflows](https://bbuchsbaum.github.io/mixeff/articles/reporting-lmms.html),
  and [selected analyses from a submitted
  manuscript](https://bbuchsbaum.github.io/mixeff/articles/reproducing-aphantasia.html).
