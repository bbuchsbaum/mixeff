# Demystifying Random-Effects Formulas

``` r

library(mixeff)
```

Random-effects formulas are compact, and compactness hides assumptions.
The four characters of `(x | g)` commit the analyst to a specific
covariance structure, a specific number of free parameters, and a
specific identification claim about what the data can support. Two
formulas that differ by a single character — `||` for `|`, say — fit
different models.

`mixeff` shares the formula language with `lme4`, so the notation you
already read transfers directly. What it adds is a way to read the
formula *before* the optimizer runs: the named form of every random
term, the plain-language scope, the parameter count, and the design
facts that support or refuse the formula on this particular data.

The goal of this vignette is not to pick a model for you. It is to show
what the model you wrote can and cannot express, and to make the
difference between nearby spellings legible.

We use a small repeated-measures dataset with one row per subject per
time point.

``` r

head(study)
#>   subject time        dose          x          z    score         y hit
#> 1       1    1 -0.59103110  1.2023289  0.9569158 2.085034 0.6487779   1
#> 2       2    1  0.02659437 -1.0278654  0.6381313 2.469368 0.6684857   1
#> 3       3    1 -1.51655310  0.9382700  1.7919729 2.546363 1.9076263   1
#> 4       4    1 -1.36265335 -0.5431547  1.1717927 1.768118 0.2011020   0
#> 5       5    1  1.17848916  0.5130951  1.4401183 2.062263 0.9445929   1
#> 6       6    1 -0.93415132 -0.3525909 -1.3062800 3.130125 0.7523965   1
```

## What does `(1 | subject)` say?

A random intercept is the simplest random-effects structure: subjects
may differ in average score, but the fixed `time` effect is the same
across subjects. It is often the right choice, and it is always a
reasonable starting point.
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
spells out what the model assumes and reports that no random slopes were
added.

``` r

intercept_only <- compile_model(score ~ time + (1 | subject), study)
explain_model(intercept_only)
#> Random effects explanation:
#>   formula: score ~ 1 + time + (1 | subject)
#> 
#> Random effects:
#>   r0:
#>     wrote:      (1 | subject)
#>     canonical:  (1 | subject)
#>     named form: re(group = subject, intercept = TRUE, slopes = NULL, cov = "scalar")
#>     scope:      `subject` units may differ in average outcome.
#>     covariance: scalar; theta parameters: 1
#>     support:    sufficient; group levels: 30; min rows/group: 5; median rows/group: 5
#>     variation:  intercept=not_assessed
#> 
#> Design notes:
#>   scope_note: `time` varies within `subject`, so a `subject`-level slope is structurally possible
```

Each random-effects card has two header lines worth pausing on. `wrote:`
is the term as you typed it. `canonical:` is what the parser hands
downstream — what the optimizer, the audit, and the inference contract
will actually see. They are identical here because `(1 | subject)` is
already in canonical form, and identity is a positive signal: nothing
was rewritten. The two diverge for more elaborate formulas.
`(1 | clinic/site)` is rewritten as `(1 | clinic) + (1 | clinic:site)`;
implicit intercepts are made explicit; double-bar shorthand expands to
its split-block form. The `canonical:` line is where you read off the
model the engine will actually fit.

The remaining lines decode the term further. `named form:` restates it
as a function call that names the grouping factor, the intercept, the
slopes (if any), and the covariance family. `scope:` is a plain-language
sentence. `covariance` and `support` report parameter cost and design
sufficiency.

The `Design notes:` section at the bottom records facts about the data
that the fitted formula did not use. Here `time` varies within each
subject, so a subject-level `time` slope is structurally possible. That
is reported as a note, not a correction: whether to add
`(time | subject)` is a scientific question about whether subjects
plausibly differ in their `time` slopes, not a mechanical one about what
the package can do.

## What changes when you add a random slope?

`(1 + time | subject)` requests a two-dimensional random effect for each
subject: one baseline coefficient and one `time` coefficient, with a
fitted covariance between them.

``` r

full <- compile_model(score ~ time + (1 + time | subject), study)
explain_model(full)
#> Random effects explanation:
#>   formula: score ~ 1 + time + (1 + time | subject)
#> 
#> Random effects:
#>   r0:
#>     wrote:      (1 + time | subject)
#>     canonical:  (1 + time | subject)
#>     named form: re(group = subject, intercept = TRUE, slopes = time, cov = "full")
#>     scope:      `subject` units differ in baseline and `time` slope; the model estimates whether these are associated.
#>     covariance: full; theta parameters: 3
#>     support:    sufficient; group levels: 30; min rows/group: 5; median rows/group: 5
#>     variation:  intercept=not_assessed; time=present
```

After fitting,
[`parameterization()`](https://bbuchsbaum.github.io/mixeff/reference/parameterization.md)
shows the model’s *internal* coordinates. A brief detour explains what
they are.

The random-effects covariance matrix Σ — here a 2×2 matrix describing
the variance of the intercept, the variance of the `time` slope, and
their covariance — is not optimized directly. Instead, the engine
searches over a vector θ that maps into a lower-triangular factor Λ via
Σ = σ² Λ Λᵀ. The reasons are practical: θ-space has the minimum number
of free parameters, it is better-conditioned than searching over Σ
directly, and the implied Σ is positive-semidefinite by construction.
This is the same parameterization Bates et al. (2015, *JSS* 67(1))
describe for `lme4`; `MixedModels.jl` and the engine behind `mixeff`
inherit it.

For a 2×2 random-effects matrix, Λ has three free entries — the three
rows you see below. `theta_value` is what the optimizer worked with;
`lambda_value` places that same number back into the matrix. They match
here because the mapping is the identity for this structure; a mismatch
would point to a parameterization bug rather than a modelling fact.

``` r

full_fit <- lmm(
  score ~ time + (1 + time | subject),
  study,
  control = mm_control(verbose = -1)
)
parameterization(full_fit)$table[, c("theta_name", "theta_value", "lambda_value")]
#>                     theta_name theta_value lambda_value
#> 1 theta[0:intercept,intercept]   1.0706644    1.0706644
#> 2      theta[0:time,intercept]   0.1222131    0.1222131
#> 3           theta[0:time,time]   0.5466997    0.5466997
```

A practical use of this table is boundary-fit detection: when a θ entry
is pinned at zero, the corresponding variance or covariance has reached
the edge of its parameter space, and the random-effects covariance is
reduced-rank. The headliner vignette walks through that case end-to-end.

## What do split blocks and `||` mean?

The formula `(1 | subject) + (0 + time | subject)` uses two separate
blocks for the same grouping factor. That fixes the intercept-slope
covariance to zero.

``` r

split_blocks <- compile_model(
  score ~ time + (1 | subject) + (0 + time | subject),
  study
)
explain_model(split_blocks)
#> Random effects explanation:
#>   formula: score ~ 1 + time + (1 | subject) + (0 + time | subject)
#> 
#> Random effects:
#>   subject has 2 separate random-effect blocks.
#>   r0:
#>     wrote:      (1 | subject)
#>     canonical:  (1 | subject)
#>     named form: re(group = subject, intercept = TRUE, slopes = NULL, cov = "scalar")
#>     scope:      `subject` units may differ in average outcome.
#>     covariance: scalar; theta parameters: 1
#>     support:    sufficient; group levels: 30; min rows/group: 5; median rows/group: 5
#>     variation:  intercept=not_assessed
#>   r1:
#>     wrote:      (0 + time | subject)
#>     canonical:  (0 + time | subject)
#>     named form: re(group = subject, intercept = FALSE, slopes = time, cov = "scalar")
#>     scope:      `subject` units may differ in their `time` slope.
#>     covariance: scalar; theta parameters: 1
#>     support:    sufficient; group levels: 30; min rows/group: 5; median rows/group: 5
#>     variation:  time=present
#> 
#> Relationship between blocks:
#>   r0 <-> r1 (Intercept <-> time): separate random-effect blocks fix the covariance between `Intercept` and `time` to zero.
#> 
#> Design notes:
#>   random_slope_without_intercept: random slope term omits a random intercept; this leaves baseline grouping dependence unmodeled unless represented elsewhere
#>   covariance_assumption: the covariance between 'Intercept' and 'time' is fixed at zero by separate random-effect blocks
```

`(1 + time || subject)` is a shorthand for the same split-block model.

``` r

double_bar <- compile_model(score ~ time + (1 + time || subject), study)
explain_model(double_bar)
#> Random effects explanation:
#>   formula: score ~ 1 + time + (1 + time || subject)
#> 
#> Random effects:
#>   subject has 2 separate random-effect blocks.
#>   r0:
#>     wrote:      (1 + time || subject)
#>     canonical:  (1 | subject)
#>     named form: re(group = subject, intercept = TRUE, slopes = NULL, cov = "scalar")
#>     scope:      `subject` units may differ in average outcome.
#>     covariance: scalar; theta parameters: 1
#>     support:    sufficient; group levels: 30; min rows/group: 5; median rows/group: 5
#>     variation:  intercept=not_assessed
#>   r1:
#>     wrote:      (1 + time || subject)
#>     canonical:  (0 + time | subject)
#>     named form: re(group = subject, intercept = FALSE, slopes = time, cov = "scalar")
#>     scope:      `subject` units may differ in their `time` slope.
#>     covariance: scalar; theta parameters: 1
#>     support:    sufficient; group levels: 30; min rows/group: 5; median rows/group: 5
#>     variation:  time=present
#> 
#> Relationship between blocks:
#>   r0 <-> r1 (Intercept <-> time): double-bar syntax fixes the covariance between `Intercept` and `time` to zero.
#> 
#> Design notes:
#>   covariance_assumption: the covariance between 'Intercept' and 'time' is fixed at zero by || syntax
```

One caveat when a **factor** appears inside `||`: mixeff decorrelates
*everything* in the block, including the factor’s level contrasts — each
treatment-coded contrast receives an independent variance and no
within-factor covariances are estimated (announced via a
`covariance_assumption` diagnostic with reason
`double_bar_factor_term`). Implementations differ here — `lme4`’s `||`
leaves factor terms intact with a full within-factor covariance block —
so when matching an external fit, write the intended expansion
explicitly: `(1 | g) + (0 + f | g) + (0 + x | g)` keeps `f`’s level
covariances, while `(1 + f + x || g)` does not.

## What are the three kinds of help?

The audit surface separates three situations that are often mixed
together in ordinary model output.

First, a structural impossibility: a requested random slope cannot be
estimated if that variable does not vary within the group. Here `dose`
is constant within each subject.

``` r

impossible <- compile_model(score ~ dose + (1 + dose | subject), between_study)
explain_model(impossible)
#> Random effects explanation:
#>   formula: score ~ 1 + dose + (1 + dose | subject)
#> 
#> Random effects:
#>   r0:
#>     wrote:      (1 + dose | subject)
#>     canonical:  (1 + dose | subject)
#>     named form: re(group = subject, intercept = TRUE, slopes = dose, cov = "full")
#>     scope:      `subject` units differ in baseline and `dose` slope; the model estimates whether these are associated.
#>     covariance: full; theta parameters: 3
#>     support:    sufficient; group levels: 30; min rows/group: 5; median rows/group: 5
#>     variation:  dose=absent; intercept=not_assessed
#> 
#> Possible repairs, not applied automatically:
#>   1. structural_refusal: `dose` does not vary within `subject`, so a `subject`-level `dose` slope cannot be estimated from this design.
```

Second, low information budget: a requested full covariance can be
estimable in principle while still having little grouping-level support.
The package reports the parameter count and observed levels.

``` r

low_info <- compile_model(score ~ time + (1 + time | subject), low_study)
diagnostics(low_info)$table[, c("code", "severity", "message")]
#>                  code severity
#> 1 covariance_too_rich  warning
#> 2 covariance_too_rich  warning
#>                                                                              message
#> 1 4 levels are below the v0 full-covariance threshold 15 for 3 covariance parameters
#> 2 4 levels are below the v0 full-covariance threshold 15 for 3 covariance parameters
```

Third, unmodeled-but-possible: a fixed effect varies within group, but
the random-effects formula does not include the corresponding slope.
This is the quiet design note you saw for `(1 | subject)`.

``` r

diagnostics(intercept_only)$table[, c("code", "severity", "message")]
#>         code severity
#> 1 scope_note     info
#> 2 scope_note     info
#>                                                                               message
#> 1 `time` varies within `subject`, so a `subject`-level slope is structurally possible
#> 2 `time` varies within `subject`, so a `subject`-level slope is structurally possible
```

## How do you inspect nearby formulas?

[`random_options()`](https://bbuchsbaum.github.io/mixeff/reference/random_options.md)
is a map of nearby random-effect spellings for one grouping factor. It
marks the current spelling and reports parameter costs and support
facts. It does not rank the rows.

``` r

random_options(intercept_only, group = subject, slope = time)
#> Random-effect options for group: subject
#> Current model:
#>   (1 | subject) <- this is what you wrote
#>   `subject` units may differ in average outcome.
#> Nearby options:
#>   (1 | subject) <- this is what you wrote
#>     varying coefficients: intercept
#>     covariance family:    scalar
#>     theta parameters:     1
#>     design status:        sufficient
#>     plain meaning:        `subject` units may differ in average outcome.
#>   (0 + time | subject)
#>     varying coefficients: time
#>     covariance family:    scalar
#>     theta parameters:     1
#>     design status:        sufficient
#>     plain meaning:        `subject` units may differ in their `time` slope.
#>   (1 | subject) + (0 + time | subject)
#>     varying coefficients: intercept, time
#>     covariance family:    diagonal via separate blocks
#>     theta parameters:     2
#>     design status:        sufficient
#>     plain meaning:        `subject` units may differ in average outcome. `subject` units may differ in their `time` slope. separate random-effect blocks fix the covariance between `Intercept` and `time` to zero.
#>   (1 + time || subject)
#>     varying coefficients: intercept, time
#>     covariance family:    diagonal via separate blocks
#>     theta parameters:     2
#>     design status:        sufficient
#>     plain meaning:        `subject` units may differ in average outcome. `subject` units may differ in their `time` slope. double-bar syntax fixes the covariance between `Intercept` and `time` to zero.
#>   (1 + time | subject)
#>     varying coefficients: intercept, time
#>     covariance family:    full
#>     theta parameters:     3
#>     design status:        sufficient
#>     plain meaning:        `subject` units differ in baseline and `time` slope; the model estimates whether these are associated.
```

The same information is available as a data frame if you want to build a
custom report.

``` r

opts <- random_options(intercept_only, group = subject, slope = time)
opts$options[, c("formula", "theta_parameters", "design_status", "current")]
#>                                formula theta_parameters design_status current
#> 1                        (1 | subject)                1    sufficient    TRUE
#> 2                 (0 + time | subject)                1    sufficient   FALSE
#> 3 (1 | subject) + (0 + time | subject)                2    sufficient   FALSE
#> 4                (1 + time || subject)                2    sufficient   FALSE
#> 5                 (1 + time | subject)                3    sufficient   FALSE
```

## Why is there no recommendation row?

`mixeff` treats formula explanations as an audit problem, not a
model-selection problem. The package reports the current model, nearby
spellings, assumptions, and data-support facts. Choosing among
scientifically different random-effects structures remains part of the
analysis design.

[`compare_covariance()`](https://bbuchsbaum.github.io/mixeff/reference/compare_covariance.md)
gives another view of the same principle: full, diagonal, and scalar
covariance families are displayed side by side, without a preferred row.

``` r

compare_covariance(full)
#> Covariance comparison:
#>   r0 / subject / full <- current
#>     basis:            intercept, time
#>     theta parameters: 3
#>     assumes zero:     none
#>     design status:    sufficient
#>   r0 / subject / diagonal
#>     basis:            intercept, time
#>     theta parameters: 2
#>     assumes zero:     off-diagonal covariances
#>     design status:    sufficient
#>   r0 / subject / scalar
#>     basis:            intercept, time
#>     theta parameters: 1
#>     assumes zero:     off-diagonal covariances
#>     design status:    sufficient
```

## Where do fitted-model changes appear?

After fitting,
[`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
records the journey from the formula you wrote to the model the engine
actually solved. Each row carries the stage at which something happened,
the term it affected, the status of the change, and three side-by-side
renderings: what was *requested*, what was *effective* (after design
parsing and identifiability checks), and what was *fitted* (after the
optimizer landed). A reduced-rank covariance estimate appears here as a
labelled fact, not an instruction to rewrite the formula.

The full row carries a lot of text. We show it in two passes so neither
pass wraps awkwardly in the rendered vignette. The first pass names the
stage of each change and its status:

``` r

knitr::kable(
  changes(full_fit)$table[, c("stage", "term_id", "status", "detail")]
)
```

| stage            | term_id | status    | detail                          |
|:-----------------|:--------|:----------|:--------------------------------|
| semantic_ir      |         | unchanged | formula display                 |
| certificate_time | r0      | full_rank | requested rank 2; fitted rank 2 |

The second pass shows the requested → effective → fitted arc itself,
where the three columns are the substance of the report:

``` r

knitr::kable(
  changes(full_fit)$table[, c("stage", "requested", "effective", "fitted")]
)
```

| stage | requested | effective | fitted |
|:---|:---|:---|:---|
| semantic_ir | score ~ 1 + time + (1 + time \| subject) | score ~ 1 + time + (1 + time \| subject) | converged_interior |
| certificate_time | intercept, time | intercept, time | full_rank |

For a first end-to-end fit workflow, read
[`vignette("lmm-basics", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm-basics.md).
