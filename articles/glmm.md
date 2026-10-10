# Generalized Linear Mixed Models

``` r

library(mixeff)
```

[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits
generalized linear mixed models by one of two estimators, and the
difference between them decides what inference you get:

- **`method = "joint_laplace"` (the default)** maximizes the same joint
  Laplace objective as
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) (at
  `nAGQ = 1`). It returns the full inference table — standard errors, z
  tests, and Wald confidence intervals — from a covariance the engine
  certifies.
- **`method = "pirls_profiled"`** is a fast profiled PIRLS fitter, the
  analogue of lme4’s `nAGQ = 0` (`glmm(..., nAGQ = 0)` selects it). It
  returns point estimates, fitted values, and variance components, but
  **refuses Wald standard errors, z statistics, p-values, and confidence
  intervals**: its working-Hessian covariance payload is not certified
  for inference (the full contract, including the engine’s accuracy
  caveats for the profiled approximation, is stated in
  [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)).
  Negative-binomial families, `nAGQ > 1`, and
  `inference = "working_hessian"` need this path; when `method` is not
  given,
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  selects it for them and says so with an `mm_estimator_notice` message.

Three inference routes exist, and every refusal message names them:

1.  **Certified Wald** — the default `method = "joint_laplace"` gives
    standard errors, z tests, and Wald intervals from an
    engine-certified covariance.
2.  **Parametric bootstrap** — `confint(fit, method = "bootstrap")` on
    profiled-estimator fits, including negative binomial (which has no
    joint-Laplace route).
3.  **Working-Hessian approximation (opt-in)** —
    `inference_options(fit)` shows the invocation and its caveat.

The contract behind these routes — what each one targets, the replicate
accounting the bootstrap reports, and the full availability/refusal
matrix — is stated once, in
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).
The practical rule: report with the default (`joint_laplace`); use
`pirls_profiled` when speed matters or for its parametric bootstrap.

The example uses the
[`lme4::cbpp`](https://rdrr.io/pkg/lme4/man/cbpp.html) data. As in
`lme4`, a binomial response with grouped counts can be written either as
a two-column `cbind(successes, failures)` response or as a proportion
with case `weights`; the two forms give identical fits.

**Response forms and coercion.** Logical responses are converted to 0/1;
a two-level factor response is accepted with its second level treated as
success, announced via an `mm_factor_coercion` condition; factor
responses with more than two levels are refused with a typed
`mm_data_error`.

``` r

env <- new.env(parent = emptyenv())
utils::data("cbpp", package = "lme4", envir = env)
cbpp <- get("cbpp", envir = env, inherits = FALSE)
cbpp$prop <- cbpp$incidence / cbpp$size

cbpp_family <- binomial(link = "logit")
cbind_formula  <- cbind(incidence, size - incidence) ~ period + (1 | herd)
weight_formula <- prop ~ period + (1 | herd)
```

## Family and Audit

The statistical intent is a binomial-logit GLMM with a herd random
intercept. Before fitting, the formula can be compiled and explained so
the fixed- and random-effect structure is visible before optimization
starts.
([`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
requires a plain-column response, so the audit uses the proportion form;
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) itself
accepts both.)

``` r

cbpp_spec <- compile_model(weight_formula, cbpp)
explain_model(cbpp_spec)
#> Random effects explanation:
#>   formula: prop ~ 1 + period + (1 | herd)
#> 
#> Random effects:
#>   r0:
#>     wrote:      (1 | herd)
#>     canonical:  (1 | herd)
#>     named form: re(group = herd, intercept = TRUE, slopes = NULL, cov = "scalar")
#>     scope:      `herd` units may differ in average outcome.
#>     covariance: scalar; theta parameters: 1
#>     support:    sufficient; group levels: 15; min rows/group: 1; median rows/group: 4
#>     variation:  intercept=not_assessed
```

## Fit

``` r

glmm_fit <- glmm(
  cbind_formula,
  cbpp,
  family = cbpp_family,
  control = mm_control(verbose = -1)
)

glmm_fit
#> Generalized linear mixed model fit
#> Formula: .mm_binomial_response ~ period + (1 | herd)
#> Family/link: binomial/logit
#> Method: joint_laplace (nAGQ = 1)
#> Fit status: converged_interior
#> Optimizer: trust_bq; iterations: 419; objective: 184.053
#> nobs: 56, dispersion: 1, logLik: -92.0263
#> Fixed effects:
#> (Intercept)     period2     period3     period4 
#>   -1.398530   -0.992332   -1.128670   -1.580310 
#> Audit verbs: audit(), diagnostics(), model_report()
fixef(glmm_fit)
#> (Intercept)     period2     period3     period4 
#>  -1.3985327  -0.9923324  -1.1286718  -1.5803138
VarCorr(glmm_fit)
#>  Groups Name        Std.Dev.
#>  herd   (Intercept) 0.64226
```

The weights form reproduces it exactly:

``` r

glmm_fit_w <- glmm(
  weight_formula,
  cbpp,
  family = cbpp_family,
  weights = size,
  control = mm_control(verbose = -1)
)
all.equal(fixef(glmm_fit), fixef(glmm_fit_w))
#> [1] TRUE
```

On this default (joint-Laplace) fit,
[`summary()`](https://rdrr.io/r/base/summary.html) prints the full Wald
inference table:

``` r

summary(glmm_fit, tests = "coefficients")
#> Generalized linear mixed model fit
#> Formula: .mm_binomial_response ~ period + (1 | herd)
#> Family/link: binomial/logit
#> Method: joint_laplace (nAGQ = 1)
#> Fit status: converged_interior
#> 
#>  Groups Name        Std.Dev.
#>  herd   (Intercept) 0.64226 
#> 
#> Fixed effects:
#>               Estimate Std. Error   z value  Pr(>|z|)            method
#> (Intercept) -1.3985327  0.2324725 -6.015905 1.789e-09 asymptotic_wald_z
#> period2     -0.9923324  0.3066426 -3.236121 0.0012117 asymptotic_wald_z
#> period3     -1.1286718  0.3266380 -3.455421 0.0005494 asymptotic_wald_z
#> period4     -1.5803138  0.4274366 -3.697189 0.0002180 asymptotic_wald_z
#> 
#> Wald-z reliability: moderate (asymptotic_wald_z).
```

Core extractors read from the durable R object:

``` r

head(fitted(glmm_fit))
#>          1          2          3          4          5          6 
#> 0.30820789 0.14174843 0.12595768 0.08402885 0.15480117 0.06358014
head(residuals(glmm_fit))
#>          1          2          3          4          5          6 
#> -1.4380337  0.9887413  2.3569962 -0.9368586 -0.2431828 -0.1424220
head(ranef(glmm_fit)[[1L]])
#>   (Intercept)
#> 1   0.5900217
#> 2  -0.2988977
#> 3   0.4062566
#> 4   0.0392777
#> 5  -0.1900155
#> 6  -0.4002694
c(
  logLik = as.numeric(logLik(glmm_fit)),
  deviance = deviance(glmm_fit),
  AIC = AIC(glmm_fit),
  BIC = BIC(glmm_fit)
)
#>    logLik  deviance       AIC       BIC 
#> -92.02628  73.47167 194.05256 204.17932
```

The joint route prints an up-front runtime notice before iterating — the
engine chooses its evaluation budget and announces it; cap it with
`mm_control(max_feval = )`. (The `verbose = -1` control used here
silences the notice for the vignette.)

``` r

confint(glmm_fit)
#> Confidence intervals:
#>                 2.5 %     97.5 %
#> (Intercept) -1.854171 -0.9428949
#> period2     -1.593341 -0.3913240
#> period3     -1.768871 -0.4884731
#> period4     -2.418074 -0.7425534
#> method: wald_asymptotic_from_rust_inference_table
#> status: available
```

## The fast profiled path

`method = "pirls_profiled"` fits faster and withholds Wald inference; on
this fit [`summary()`](https://rdrr.io/r/base/summary.html) prints the
estimates with a note naming the alternatives:

``` r

glmm_fast <- glmm(
  cbind_formula,
  cbpp,
  family = cbpp_family,
  method = "pirls_profiled",
  control = mm_control(verbose = -1)
)
summary(glmm_fast, tests = "coefficients")
#> Generalized linear mixed model fit
#> Formula: .mm_binomial_response ~ period + (1 | herd)
#> Family/link: binomial/logit
#> Method: pirls_profiled (nAGQ = 1)
#> Fit status: converged_interior
#> 
#>  Groups Name        Std.Dev.
#>  herd   (Intercept) 0.64182 
#> 
#> Fixed effects:
#>               Estimate Std. Error statistic p.value       method
#> (Intercept) -1.3604731         NA        NA      NA not_computed
#> period2     -0.9761761         NA        NA      NA not_computed
#> period3     -1.1110755         NA        NA      NA not_computed
#> period4     -1.5596789         NA        NA      NA not_computed
#> 
#> Wald-z reliability: not_available (not_computed).
#> 
#> Notes:
#>   standard errors, z statistics, and p-values are not available from the fast profiled method (pirls_profiled). Re-fit with method = "joint_laplace" for glmer-equivalent Wald inference, or use confint(fit, method = "bootstrap") for parametric-bootstrap intervals. inference_options(fit) lists every route.
```

The joint-Laplace estimates differ slightly from the profiled ones —
different estimators, same model — and match
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) closely:

``` r

glmer_fit <- lme4::glmer(
  cbind(incidence, size - incidence) ~ period + (1 | herd),
  data = cbpp, family = binomial()
)
round(rbind(
  pirls_profiled = fixef(glmm_fast),
  joint_laplace  = fixef(glmm_fit),
  glmer          = lme4::fixef(glmer_fit)
), 4)
#>                (Intercept) period2 period3 period4
#> pirls_profiled     -1.3605 -0.9762 -1.1111 -1.5597
#> joint_laplace      -1.3985 -0.9923 -1.1287 -1.5803
#> glmer              -1.3983 -0.9919 -1.1282 -1.5797
```

## Prediction

[`predict()`](https://rdrr.io/r/stats/predict.html) works on both
estimators. Population-level prediction (`re.form = NA`, on the link or
response scale) matches `glmer` on joint-Laplace fits. Conditional
predictions (`re.form = NULL`, the default) carry per-row engine
certificates: where the engine certifies the payload, `se.fit = TRUE`
returns finite standard errors; where it does not — an unseen grouping
level under `allow.new.levels = TRUE`, for example — the affected rows
come back `NA` with the engine’s reason in an `mm_reason` attribute, not
a silent zero.

``` r

head(predict(glmm_fit, type = "response"))
#>          1          2          3          4          5          6 
#> 0.30820789 0.14174843 0.12595768 0.08402885 0.15480117 0.06358014
head(predict(glmm_fit, re.form = NA, type = "response"))
#>          1          2          3          4          5          6 
#> 0.19804905 0.08387193 0.07397291 0.04839072 0.19804905 0.08387193
```

## Quadrature Sensitivity

`nAGQ` is part of the fit request and is recorded on the object. Values
above one are a profiled-path sensitivity check: as in lme4 they need a
single scalar random-effect term such as `(1 | herd)`, and the
joint-Laplace path is restricted to `nAGQ = 1`. Without an explicit
`method`, `nAGQ > 1` selects `"pirls_profiled"` with an
`mm_estimator_notice` (silenced here by `verbose = -1`).

``` r

glmm_fit_agq3 <- glmm(
  cbind_formula,
  cbpp,
  family = cbpp_family,
  nAGQ = 3L,
  control = mm_control(verbose = -1)
)

data.frame(
  nAGQ = c(glmm_fast$nAGQ, glmm_fit_agq3$nAGQ),
  logLik = c(as.numeric(logLik(glmm_fast)),
             as.numeric(logLik(glmm_fit_agq3))),
  AIC = c(AIC(glmm_fast), AIC(glmm_fit_agq3)),
  check.names = FALSE
)
#>   nAGQ   logLik      AIC
#> 1    1 -92.0543 194.1086
#> 2    3 -92.0451 194.0902
```

On this dataset the quadrature order shifts the coefficients by less
than 0.005 on the logit scale (the table above shows the corresponding
movement in the objective) — visible, but an order of magnitude below
the joint-Laplace standard errors shown earlier.

## What is refused where

The table below runs each verb on the profiled fit and reports what
actually comes back — a result with numbers, a result with the inference
withheld, or a typed error. A no-error return is not the same thing as
an answer, so for the coefficient table the harness also checks whether
the standard errors are finite.

``` r

glmm_status <- function(expr, has_numbers = NULL) {
  cnd <- tryCatch({ force(expr); NULL }, error = function(cnd) cnd)
  if (!is.null(cnd)) {
    return(data.frame(status = "typed error", class = class(cnd)[[1L]],
                      check.names = FALSE))
  }
  if (isFALSE(has_numbers)) {
    return(data.frame(status = "returns, inference withheld",
                      class = NA_character_, check.names = FALSE))
  }
  data.frame(status = "available", class = NA_character_, check.names = FALSE)
}

coef_tab <- summary(glmm_fast, tests = "coefficients")$coefficients
rbind(
  predict           = glmm_status(predict(glmm_fast)),
  fitted            = glmm_status(fitted(glmm_fast)),
  coefficient_tests = glmm_status(coef_tab,
                                  has_numbers = any(is.finite(coef_tab[, "Std. Error"]))),
  confint           = glmm_status(confint(glmm_fast)),
  simulate          = glmm_status(stats::simulate(glmm_fit, nsim = 1L, seed = 1L)),
  refit             = glmm_status(refit(glmm_fit,
                                        stats::simulate(glmm_fit, seed = 1L)))
)
#>                                        status                    class
#> predict                             available                     <NA>
#> fitted                              available                     <NA>
#> coefficient_tests returns, inference withheld                     <NA>
#> confint                           typed error mm_inference_unavailable
#> simulate                            available                     <NA>
#> refit                               available                     <NA>
```

Wald [`confint()`](https://rdrr.io/r/stats/confint.html) on the profiled
fit raises `mm_inference_unavailable` (the working-Hessian payload is
not certified for Wald intervals); the same request on the default
joint-Laplace fit succeeds, as shown above. The bootstrap route
complements it from the other side:
`confint(glmm_fast, method = "bootstrap")` succeeds on the profiled fit
(route 2 above) and is refused on joint-Laplace fits at this engine pin.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws family-aware
responses with lme4’s semantics (new random effects by default,
`re.form = NULL` to condition on the fitted ones), and
[`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md)
re-fits the same GLMM to a new response such as a simulated one.
