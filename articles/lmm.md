# Linear Mixed Models

``` r

library(mixeff)
```

You use a linear mixed model when observations are not independent:
visits within clinics, trials within people, students within classrooms.
The standard R answer is
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).
[`mixeff::lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
accepts the same formula language — `(1 | group)`, `(x | group)`,
`(1 + x || group)`, `(1 | a/b)`, crossed blocks — and the fitted object
answers to the common generics:
[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`anova()`](https://rdrr.io/r/stats/anova.html), and friends. The
supported surface is documented and closed rather than open-ended:
[`vignette("mixeff", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/mixeff.md)
states the contract, and
[`vignette("lme4-migration", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md)
gives the exact verb-and-argument map.

What [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
adds is that the fitted object carries the design, the convergence
status, and the inference labels with it. You can audit, summarize,
save, and reload without recomputing.

This vignette is the LMM walkthrough: reading a random-effects formula
before you fit, fitting, reading the summary, the standard extractors,
prediction, and weights — followed by a longer section on what the
random-effects spellings themselves mean and the compile-time design
help.

## The data

The example simulates reaction times (ms) for 18 subjects over ten days
of sleep deprivation. Each subject has a different baseline and a
different response to continued sleep loss. The data are generated
inside this vignette, so the walkthrough runs without an example-data
package.

``` r

str(sleep)
#> 'data.frame':    180 obs. of  3 variables:
#>  $ Days    : int  0 1 2 3 4 5 6 7 8 9 ...
#>  $ Subject : Factor w/ 18 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Reaction: num  224 248 293 300 292 ...
#>  - attr(*, "out.attrs")=List of 2
#>   ..$ dim     : Named int [1:2] 10 18
#>   .. ..- attr(*, "names")= chr [1:2] "Days" "Subject"
#>   ..$ dimnames:List of 2
#>   .. ..$ Days   : chr [1:10] "Days=0" "Days=1" "Days=2" "Days=3" ...
#>   .. ..$ Subject: chr [1:18] "Subject=1" "Subject=2" "Subject=3" "Subject=4" ...
```

The key structure: `Subject` appears 10 times in the data — once per
day. Those 10 observations are *not* independent. They share the
subject’s baseline speed, their individual sensitivity to sleep loss,
and any other between-person trait that we haven’t measured. If you fit
a plain OLS regression, those correlations inflate your Type I error
rate and make standard errors too small.

A mixed model handles this by giving each subject its own intercept
*and* its own slope — their personal baseline and their personal rate of
slowing — and then estimating a population distribution over those
person-level parameters.

![](lmm_files/figure-html/individual-trajectories-1.png)

The grey lines are individuals; the blue line is the population average.
Subjects vary in their starting point *and* in how fast they slow down —
which is exactly what a random-intercept-and-slope model captures.

## Step 1: Read the formula before you fit

Before fitting, compile the model to see what `mixeff` understands about
your formula. This separates formula interpretation from optimization;
you can catch mis-specified random effects without paying the cost of a
full fit.

``` r

spec <- compile_model(Reaction ~ Days + (Days | Subject), sleep)
explain_model(spec)
#> Random effects explanation:
#>   formula: Reaction ~ 1 + Days + (1 + Days | Subject)
#> 
#> Random effects:
#>   r0:
#>     wrote:      (Days | Subject)
#>     canonical:  (1 + Days | Subject)
#>     named form: re(group = Subject, intercept = TRUE, slopes = Days, cov = "full")
#>     scope:      `Subject` units differ in baseline and `Days` slope; the model estimates whether these are associated.
#>     covariance: full; theta parameters: 3
#>     support:    sufficient; group levels: 18; min rows/group: 10; median rows/group: 10
#>     variation:  Days=present; intercept=not_assessed
```

This pre-fit explanation confirms:

- one fixed effect (`Days`) plus an intercept,
- one random-effects block with a correlated intercept and slope per
  subject,
- the formula has been canonicalized to the explicit
  `(1 + Days | Subject)` form that the optimizer will receive.

Use `audit(spec)` when you want the deeper design audit. Its compact
print starts with the audit summary and requested model;
`print(audit(spec), full = TRUE)` shows the complete upstream report.

If you had written `(1 | Subject)` by mistake (random intercepts only),
the explanation would show only one random-effects column per subject —
useful to verify before a long fit. The second half of this vignette
works through the random-effects spellings, canonical forms, and design
diagnostics in detail.

## Step 2: Fit

``` r

fit <- lmm(Reaction ~ Days + (Days | Subject), sleep)
```

[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) returns
an `mm_lmm` object. The bundled Rust engine performs the optimization;
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) is the R
entry point and result container.

Printing the fit gives the formula, convergence status, likelihood
summary, residual scale, and fixed effects.

``` r

fit
#> Linear mixed model fit by REML
#> Formula: Reaction ~ Days + (Days | Subject)
#> Fit status: converged_interior
#> Optimizer: trust_bq; iterations: 438; objective: 1715.67
#> nobs: 180, sigma: 23.5051, logLik: -857.837
#> Fixed effects:
#> (Intercept)        Days 
#>    256.2570     10.1502 
#> Audit verbs: audit(), diagnostics(), inference_table(), model_report()
```

## Step 3: Read the summary

``` r

summary(fit, tests = "coefficients")
#> Linear mixed model fit by REML
#> Formula: Reaction ~ Days + (Days | Subject)
#> Fit status: converged_interior
#> 
#> Variance components:
#>    group        name variance  std_dev correlation
#>  Subject (Intercept) 652.6850 25.54770            
#>  Subject        Days  29.5779  5.43856       -0.08
#> Residual std. dev.: 23.5051
#> 
#> Fixed effects:
#>              Estimate Std. Error       df   t value  Pr(>|t|)        method
#> (Intercept) 256.25687   6.845701 16.99932 37.433253   < 1e-16 satterthwaite
#> Days         10.15023   1.419601 16.99972  7.150056 1.626e-06 satterthwaite
#> 
#> Inference status:
#>         term        method    status reliability
#>  (Intercept) satterthwaite available    moderate
#>         Days satterthwaite available    moderate
#>                             reliability_reason
#>  satterthwaite_finite_difference_approximation
#>  satterthwaite_finite_difference_approximation
#> 
#> Notes:
#>   Satterthwaite denominator df computed from finite-difference vcov_beta Jacobian and deviance Hessian over varpar
```

The summary has four blocks, in print order.

**Fit status.** Check this line first: `converged_interior` is the
status every other block of the summary assumes. When a fit lands
somewhere else — a variance component at zero, a reduced-rank covariance
— the status says so and the available inference methods change with it.
What the statuses mean, and what a boundary does to the rest of the
output, is the subject of
[`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md).

**Variance components.** `Subject (Intercept)` is the estimated
between-subject spread in baseline reaction time; `Subject Days` is the
estimated between-subject spread in sensitivity to sleep loss;
`correlation` is their estimated association. The residual standard
deviation is the within-subject noise left over.

**Fixed effects.** `Estimate` is the population-level coefficient. The
`Days` estimate is the fitted average change in reaction time per day;
the simulation used a population slope of 10.5 ms per day. The `method`
column names how each p-value was computed — here `satterthwaite`, a
finite-sample *t* test whose `df` column is on the scale of the 18
subjects, not the 180 raw rows. Which methods exist, when each applies,
and how to choose among them is covered in
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).

**Inference status.** One row per coefficient stating how much to trust
the test: `status` says whether it was computed, `reliability` grades
it, and `reliability_reason` names the engine’s reason for that grade —
here `satterthwaite_finite_difference_approximation`, the standard route
for this method. The Rust engine assigns these grades and reasons; the R
side passes them through unchanged. The closed vocabularies behind these
columns are defined in the [inference-method
glossary](https://bbuchsbaum.github.io/mixeff/articles/inference-method-glossary.html)
on the package site.

## Step 4: Extract components

The lme4-compatible extractors work on `mm_lmm` objects.

``` r

fixef(fit)
#> (Intercept)        Days 
#>   256.25687    10.15023
sigma(fit)
#> [1] 23.50513
logLik(fit)
#> 'log Lik.' -857.8369 (df=6)
```

[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
reports fitted variance components, and
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
returns conditional random effects by grouping factor.

``` r

VarCorr(fit)
#> Variance components:
#>    group        name variance  std_dev correlation
#>  Subject (Intercept) 652.6850 25.54770            
#>  Subject        Days  29.5779  5.43856       -0.08
#> Residual std. dev.: 23.5051
head(ranef(fit)$Subject)
#>   (Intercept)      Days
#> 1   -8.047967  2.798294
#> 2    4.457804  4.828437
#> 3  -17.180029  4.799092
#> 4   36.246461  4.003065
#> 5   -4.147918  1.625769
#> 6  -23.075829 -8.536017
```

``` r

confint(fit, method = "asymptotic")
#> Confidence intervals:
#>                 2.5 %    97.5 %
#> (Intercept) 242.83954 269.67420
#> Days          7.36786  12.93259
#> method: wald_asymptotic_from_stored_standard_errors
#> status: Wald (asymptotic) intervals from stored standard errors (engine-certified profile intervals: method = "profile")
```

[`confint()`](https://rdrr.io/r/stats/confint.html) with
`method = "asymptotic"` gives Wald intervals. The `(Intercept)` interval
is the population mean reaction time on Day 0; the `Days` interval is
the range of plausible slopes. The other interval methods, and when each
is available, are covered in
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).

Design matrices can be rebuilt from the stored formula and model frame —
including after a save/reload cycle
([`vignette("mixeff", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/mixeff.md)
covers persistence).
[`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md)
provides a small familiar subset for code that expects lme4-style names.

``` r

X <- model.matrix(fit, type = "fixed")
Z <- model.matrix(fit, type = "random")

dim(X)
#> [1] 180   2
dim(Z)
#> [1] 180  36
class(Z)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"

getME(fit, c("theta", "beta", "cnms"))
#> $theta
#> [1]  1.08689842 -0.01922381  0.23057757
#> 
#> $beta
#> (Intercept)        Days 
#>   256.25687    10.15023 
#> 
#> $cnms
#> $Subject
#> [1] "(Intercept)" "Days"       
#> 
#> attr(,"class")
#> [1] "mm_cnms" "list"
```

For report-ready data-frame tables of the same fit — design summary,
random terms, fixed effects — see the [reporting
article](https://bbuchsbaum.github.io/mixeff/articles/reporting-lmms.html)
on the package site.

## Step 5: Prediction and residuals

For fitted data, [`predict()`](https://rdrr.io/r/stats/predict.html)
returns in-sample fitted values. Use `re.form = NA` for the fixed-effect
part only.

``` r

prediction_check <- data.frame(
  Reaction = sleep$Reaction,
  fitted = predict(fit),
  fixed_only = predict(fit, re.form = NA),
  residual = residuals(fit)
)
head(prediction_check)
#>   Reaction   fitted fixed_only   residual
#> 1 224.4814 248.2089   256.2569 -23.727494
#> 2 248.2831 261.1574   266.4071 -12.874274
#> 3 292.6939 274.1059   276.5573  18.588001
#> 4 299.7000 287.0545   286.7075  12.645567
#> 5 291.9349 300.0030   296.8578  -8.068110
#> 6 305.1412 312.9515   307.0080  -7.810256
```

### Prediction standard errors and intervals

Population-level prediction standard errors and intervals are available
via `re.form = NA` (the Wald SE of the fixed-effect linear predictor):

``` r

pop <- predict(fit, re.form = NA, se.fit = TRUE)
head(pop$se.fit)
#>        1        2        3        4        5        6 
#> 6.845701 6.651902 6.757404 7.148969 7.783542 8.607545
head(predict(fit, re.form = NA, interval = "confidence"))
#>        fit      lwr      upr
#> 1 256.2569 242.8395 269.6742
#> 2 266.4071 253.3696 279.4446
#> 3 276.5573 263.3131 289.8016
#> 4 286.7075 272.6958 300.7193
#> 5 296.8578 281.6023 312.1132
#> 6 307.0080 290.1375 323.8785
```

*Conditional* prediction standard errors (the default, `re.form = NULL`)
come from the engine’s prediction-variance payload, which combines the
fixed-effect Wald variance with the random-effect (BLUP) variance and
the fixed/random covariance. The covariance term can be negative, so the
conditional SE can come out *smaller* than the population SE.
Conditional confidence and prediction intervals come from the same
payload.

``` r

pred <- predict(fit, se.fit = TRUE)
head(data.frame(population_se = pop$se.fit, conditional_se = pred$se.fit))
#>   population_se conditional_se
#> 1      6.845701      11.757630
#> 2      6.651902      10.086320
#> 3      6.757404       8.665374
#> 4      7.148969       7.635864
#> 5      7.783542       7.168472
#> 6      8.607545       7.370917
head(predict(fit, interval = "confidence"))
#>        fit      lwr      upr
#> 1 248.2089 225.1644 271.2534
#> 2 261.1574 241.3886 280.9262
#> 3 274.1059 257.1221 291.0898
#> 4 287.0545 272.0884 302.0205
#> 5 300.0030 285.9530 314.0529
#> 6 312.9515 298.5048 327.3982
```

Rows the engine cannot certify — for example an unseen grouping level
predicted with `allow.new.levels = TRUE`, where there is no posterior
variance for the missing level — return `NA` with the engine’s reason in
the `mm_reason` attribute rather than a fabricated number. That refusal
vocabulary is part of the inference contract
([`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)).

## Step 6: Case weights

[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) accepts
case weights: one finite, positive value per row, evaluated in `data`
like the formula variables. A common use is rows that summarize
different numbers of underlying trials.

``` r

sleep$n_trials <- rep(c(1, 2), length.out = nrow(sleep))

fit_w <- lmm(Reaction ~ Days + (Days | Subject), sleep,
             weights = n_trials,
             control = mm_control(verbose = -1))
fixef(fit_w)
#> (Intercept)        Days 
#>   255.72882    10.30777
head(fit_w$weights)
#> [1] 1 2 1 2 1 2
```

The weights are stored on the fit and reused by
[`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md) and
the inference routines. Invalid weights — wrong length, negative, or
non-finite — are refused with a typed `mm_data_error` rather than
silently recycled.

For contrasts, term tests, marginal means, and model comparisons on a
fit like this one, continue with
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).
The rest of this vignette turns to the random-effects formulas
themselves.

## Random-effects formulas in depth

Random-effects formulas are compact, and compactness hides assumptions.
Seven characters — `(x | g)` — commit the analyst to a specific
covariance structure, a specific number of free parameters, and a
specific identification claim about what the data can support. Two
formulas that differ by a single character — `||` for `|`, say — fit
different models.

`mixeff` shares the formula language with `lme4`, so the notation you
already read transfers directly. What it adds is a way to read the
formula *before* the optimizer runs: the named form of every random
term, the plain-language scope, the parameter count, and the design
facts that support or refuse the formula on this particular data.

The goal of this section is not to pick a model for you. It is to show
what the model you wrote can and cannot express, and to make the
difference between nearby spellings legible.

We use a second small repeated-measures dataset with one row per subject
per time point.

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

### What does `(1 | subject)` say?

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
already in canonical form — nothing was rewritten. The two diverge for
more elaborate formulas. `(1 | clinic/site)` is rewritten as
`(1 | clinic) + (1 | clinic:site)`; implicit intercepts are made
explicit; double-bar shorthand expands to its split-block form. The
`canonical:` line is where you read off the model the engine will
actually fit.

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

### What changes when you add a random slope?

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
shows the optimizer’s *internal* coordinates: the vector θ the engine
actually searched over, and where each entry sits in the random-effects
covariance factor.

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

What θ is, why the engine optimizes it rather than the covariance matrix
directly, and how boundary fits are read off this table are covered in
[`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md).

### What do split blocks and `||` mean?

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

One caveat: when a **factor** appears inside `||`, mixeff’s expansion
differs from lme4’s —
[`vignette("lme4-migration", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md)
states the exact semantics and the explicit spelling that reproduces
lme4’s model family.

### Nested and crossed grouping

Nesting is written with `/` and expands mechanically:
`(1 | clinic/ward)` means one intercept block per clinic plus one per
ward-within-clinic. Crossed factors are simply two blocks —
`(1 | subject) + (1 | item)` — and nothing is rewritten, because that
spelling is already canonical.
[`mm_parse_formula()`](https://bbuchsbaum.github.io/mixeff/reference/mm_parse_formula.md)
shows the canonical string for any formula without touching data:

``` r

mm_parse_formula("rt ~ days + (1 | subject/task)")
#> [1] "rt ~ 1 + days + (1 | subject) + (1 | subject:task)"
mm_parse_formula("rt ~ days + (1 | subject) + (1 | item)")
#> [1] "rt ~ 1 + days + (1 | subject) + (1 | item)"
mm_parse_formula("rt ~ days + (1 + days || subject)")
#> [1] "rt ~ 1 + days + (1 + days || subject)"
```

The
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
cards agree with these canonicalizations — they show the model the
engine actually fits — with one display difference:
[`mm_parse_formula()`](https://bbuchsbaum.github.io/mixeff/reference/mm_parse_formula.md)
keeps a `||` term as written, while the cards additionally show it as
the split blocks the engine fits.

### What are the three kinds of help?

The compile-time diagnostics separate three situations that are often
mixed together in ordinary model output.

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

### How do you inspect nearby formulas?

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

### Why is there no recommendation row?

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

### Where do fitted-model changes appear?

After fitting,
[`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
records the journey from the formula you wrote to the model the engine
actually solved. Each row carries the stage at which something happened,
the term it affected, the status of the change, and three side-by-side
renderings: what was *requested*, what was *effective* (after design
parsing and identifiability checks), and what was *fitted* (after the
optimizer landed). A reduced-rank covariance estimate appears here as a
labeled fact, not an instruction to rewrite the formula (boundary and
reduced-rank outcomes are the subject of
[`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md)).

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

## Where to go next

| Question | Where |
|----|----|
| What is the supported envelope, and how do fits persist across sessions? | [`vignette("mixeff", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/mixeff.md) |
| How are p-values computed, and which method is right for my fit? | [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md) |
| How do I compute marginal means and treatment contrasts? | [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md) |
| What do convergence statuses, boundaries, and singular fits mean? | [`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md) |
| Can I fit a binomial/Poisson/negative-binomial GLMM? | [`vignette("glmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md) |
| I’m coming from lme4 — what’s different? | [`vignette("lme4-migration", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lme4-migration.md) |
| How do I write up the results for a paper? | the [reporting article](https://bbuchsbaum.github.io/mixeff/articles/reporting-lmms.html) on the package site |
