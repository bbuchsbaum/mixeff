# Introduction to mixeff

``` r

library(mixeff)
```

If you have ever fit a mixed model in R, you have almost certainly
written [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).
[`mixeff::lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
accepts the same formulas and answers to the same extractors. An `lmer`
script becomes an `lmm` script with two edits.

What `mixeff` adds is mostly *around* the fit. It exposes the
random-effects design before optimization, keeps the method behind every
p-value visible in the output, and produces a fitted object that
survives [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) without
losing its audit trail. The goal is not to replace your statistical
judgment; it is to make the fitted object carry enough information that
you can inspect, report, save, reload, and compare models without
guessing which details were used.

This vignette is a guided tour at a deliberate pace. The shorter
elevator pitch is in
[`vignette("mixeff", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/mixeff.md).

## What problem does it solve?

In ordinary mixed-model work, three questions come up again and again:

- What does this random-effects formula actually mean?
- Are the p-values or tests available, and by what method?
- Can I reconstruct the model state later, after the fit has been saved?

`mixeff` makes those questions part of the fitted object.

## One small model

``` r

fit <- lmm(
  score ~ week + treatment + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)
```

The same object gives the fitted coefficients and the inferential status
of those coefficients.

``` r

coef_table <- summary(fit, method = "auto")$coefficients
knitr::kable(coef_table, digits = 4)
```

|  | Estimate | Std. Error | df | t value | Pr(\>\|t\|) | method |
|:---|---:|---:|---:|---:|---:|:---|
| (Intercept) | 7.6829 | 0.1965 | 12.5650 | 39.1065 | 0.0000 | satterthwaite |
| week | -0.2784 | 0.0260 | 58.9997 | -10.7280 | 0.0000 | satterthwaite |
| treatmentcoached | -0.8995 | 0.2623 | 9.9993 | -3.4298 | 0.0064 | satterthwaite |

## Reading the formula before fitting

If you are unsure what a random-effects expression actually models, the
answer is
[`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
followed by
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md).
[`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
builds the pre-fit specification — the same specification the optimizer
will receive — and
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
translates each random term into a named form, a plain-language scope,
and a parameter count.

``` r

spec <- compile_model(score ~ week + treatment + (1 | clinic), clinic_visits)
explain_model(spec)
#> Random effects explanation:
#>   formula: score ~ 1 + week + treatment + (1 | clinic)
#> 
#> Random effects:
#>   r0:
#>     wrote:      (1 | clinic)
#>     canonical:  (1 | clinic)
#>     named form: re(group = clinic, intercept = TRUE, slopes = NULL, cov = "scalar")
#>     scope:      `clinic` units may differ in average outcome.
#>     covariance: scalar; theta parameters: 1
#>     support:    sufficient; group levels: 12; min rows/group: 6; median rows/group: 6
#>     variation:  intercept=not_assessed
#> 
#> Design notes:
#>   scope_note: `week` varies within `clinic`, so a `clinic`-level slope is structurally possible
```

[`vignette("demystifying-formulas", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/demystifying-formulas.md)
works through the random-effects spellings in detail: scalar versus
diagonal versus full covariance, `||` shorthand, nested grouping, and
the difference between “this formula cannot be estimated” and “this
formula can be estimated but the data are sparse”.

A lower-level utility —
[`mm_parse_formula()`](https://bbuchsbaum.github.io/mixeff/reference/mm_parse_formula.md)
— exists for checking that a formula string parses at all and reducing
two equivalent spellings to the same canonical string. It is the
primitive that equivalence-class testing uses, not a reader-facing
explanation, so it is not what you want when the question is “what does
this formula mean?”.

## Reporting tables

Use
[`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md)
when you want a data-frame result instead of printed console output. The
default view is compact; use `view = "audit"` when you want the full
provenance columns.

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
reporting_table(fit, "fixed_effects")
#>                term   estimate  std_error  statistic statistic_name
#>         (Intercept)  7.6828778 0.19646018  39.106539              z
#>                week -0.2783994 0.02595083 -10.727955              z
#>  treatment: coached -0.8994747 0.26225014  -3.429835              z
#>       p_value            method    status reliability
#>  0.0000000000 asymptotic_wald_z available         low
#>  0.0000000000 asymptotic_wald_z available         low
#>  0.0006039485 asymptotic_wald_z available         low
reporting_table(fit, "fixed_effects", view = "audit")$table[, c("term", "source", "status")]
#>                 term                       source    status
#> 1        (Intercept) fixed_effect_inference_table available
#> 2               week fixed_effect_inference_table available
#> 3 treatment: coached fixed_effect_inference_table available
```

## Saving and reloading

The fitted object stores the values needed by the main extractors, so an
RDS round trip preserves the pieces you usually report.

``` r

path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
restored <- readRDS(path)

fixef(restored)
#>      (Intercept)             week treatmentcoached 
#>        7.6828778       -0.2783994       -0.8994747
reporting_table(restored, "fixed_effects")
#>                term   estimate  std_error  statistic statistic_name
#>         (Intercept)  7.6828778 0.19646018  39.106539              z
#>                week -0.2783994 0.02595083 -10.727955              z
#>  treatment: coached -0.8994747 0.26225014  -3.429835              z
#>       p_value            method    status reliability
#>  0.0000000000 asymptotic_wald_z available         low
#>  0.0000000000 asymptotic_wald_z available         low
#>  0.0006039485 asymptotic_wald_z available         low
```

## Lower-level tools

Most users should start with
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md),
[`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md),
and
[`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md).
The lower-level functions are there when you need them:

- [`mm_parse_formula()`](https://bbuchsbaum.github.io/mixeff/reference/mm_parse_formula.md)
  checks formula syntax.
- [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  builds a pre-fit model specification.
- [`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
  and
  [`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
  expose model-state checks.
- [`mm_json_known_schemas()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_known_schemas.md)
  lists the structured artifact schemas understood by this version of
  the package.

The computational backend is intentionally not the opening story for
most R users. It matters because it lets `mixeff` keep a structured
audit trail, but the user-facing reason to use the package is simpler:
fit the model, get the numbers, and keep the status of those numbers
attached to the object.

## What’s next?

Use
[`vignette("mixeff", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/mixeff.md)
for the fastest end-to-end path. Use
[`vignette("lmm-basics", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm-basics.md)
for a slower fitted-model walkthrough. Use
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
for p-values, contrasts, term tests, and model comparisons. Use
[`vignette("demystifying-formulas", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/demystifying-formulas.md)
for random-effects syntax.
