# GLMM Fitting and Model Comparison

``` r

library(mixeff)
```

[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits
generalized linear mixed models by one of two estimators, and the
difference between them decides what inference you get:

- **`method = "pirls_profiled"` (the default)** is a fast profiled PIRLS
  fitter. It returns point estimates, fitted values, and variance
  components, but **refuses Wald standard errors, z statistics,
  p-values, and confidence intervals**: its covariance payload is a
  working Hessian that the engine does not certify for inference. The
  engine’s own documentation adds a stronger warning: the profiled path
  is not the same statistical approximation as `glmer()`’s joint Laplace
  fit, and it can be less accurate for overdispersed models and models
  with observation-level random effects.
- **`method = "joint_laplace"`** maximizes the same joint Laplace
  objective as
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) (at
  `nAGQ = 1`). It is slower, and it returns the full inference table —
  standard errors, z tests, and Wald confidence intervals — from a
  covariance the engine certifies.

The practical rule: explore with the default, report with
`joint_laplace`.

The example uses the
[`lme4::cbpp`](https://rdrr.io/pkg/lme4/man/cbpp.html) data. As in
`lme4`, a binomial response with grouped counts can be written either as
a two-column `cbind(successes, failures)` response or as a proportion
with case `weights`; the two forms give identical fits.

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
#> Method: pirls_profiled (nAGQ = 1)
#> Fit status: converged_interior
#> Optimizer: cobyla; iterations: 23; objective: 100.152
#> nobs: 56, dispersion: 1, logLik: -92.0543
#> Fixed effects:
#> (Intercept)     period2     period3     period4 
#>   -1.360470   -0.976176   -1.111080   -1.559680 
#> Audit verbs: audit(), diagnostics(), model_report()
fixef(glmm_fit)
#> (Intercept)     period2     period3     period4 
#>  -1.3604731  -0.9761761  -1.1110755  -1.5596789
VarCorr(glmm_fit)
#> Variance components:
#>  group        name variance  std_dev correlation
#>   herd (Intercept) 0.411937 0.641823
```

The weights form reproduces it exactly:

``` r

glmm_fit_w <- glmm(
  weight_formula,
  cbpp,
  family = cbpp_family,
  weights = cbpp$size,
  control = mm_control(verbose = -1)
)
all.equal(fixef(glmm_fit), fixef(glmm_fit_w))
#> [1] TRUE
```

On this default fit, [`summary()`](https://rdrr.io/r/base/summary.html)
prints the estimates and withholds the inference columns, with a note
naming the alternative:

``` r

summary(glmm_fit, tests = "coefficients")
#> Generalized linear mixed model fit
#> Formula: .mm_binomial_response ~ period + (1 | herd)
#> Family/link: binomial/logit
#> Method: pirls_profiled (nAGQ = 1)
#> Fit status: converged_interior
#> 
#> Variance components:
#>  group        name variance  std_dev correlation
#>   herd (Intercept) 0.411937 0.641823            
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
#>   standard errors, z statistics, and p-values are not available from the fast default method (pirls_profiled). Re-fit with method = "joint_laplace" for glmer-equivalent Wald inference.
```

Core extractors read from the durable R object:

``` r

head(fitted(glmm_fit))
#>          1          2          3          4          5          6 
#> 0.30926668 0.14433698 0.12846211 0.08602001 0.15578121 0.06500146
head(residuals(glmm_fit))
#>            1            2            3            4            5            6 
#> -0.166409542  0.105663020  0.315982339 -0.086020007 -0.019417571 -0.009445903
head(ranef(glmm_fit)[[1L]])
#>   (Intercept)
#> 1  0.55692322
#> 2 -0.32948612
#> 3  0.36897593
#> 4  0.01190522
#> 5 -0.22338936
#> 6 -0.43117113
c(
  logLik = as.numeric(logLik(glmm_fit)),
  deviance = deviance(glmm_fit),
  AIC = AIC(glmm_fit),
  BIC = BIC(glmm_fit)
)
#>   logLik deviance      AIC      BIC 
#> -92.0543 184.1086 194.1086 204.2353
```

## Inference: refit with joint_laplace

``` r

glmm_joint <- glmm(
  cbind_formula,
  cbpp,
  family = cbpp_family,
  method = "joint_laplace",
  control = mm_control(verbose = -1)
)

summary(glmm_joint, tests = "coefficients")
#> Generalized linear mixed model fit
#> Formula: .mm_binomial_response ~ period + (1 | herd)
#> Family/link: binomial/logit
#> Method: joint_laplace (nAGQ = 1)
#> Fit status: converged_interior
#> 
#> Variance components:
#>  group        name variance  std_dev correlation
#>   herd (Intercept) 0.412513 0.642272            
#> 
#> Fixed effects:
#>               Estimate Std. Error   z value  Pr(>|z|)            method
#> (Intercept) -1.3985345  0.2324741 -6.015872 1.789e-09 asymptotic_wald_z
#> period2     -0.9923309  0.3066426 -3.236115 0.0012117 asymptotic_wald_z
#> period3     -1.1286704  0.3266381 -3.455416 0.0005494 asymptotic_wald_z
#> period4     -1.5803119  0.4274367 -3.697184 0.0002180 asymptotic_wald_z
#> 
#> Wald-z reliability: moderate (asymptotic_wald_z).
confint(glmm_joint)
#> Confidence intervals:
#>                 2.5 %     97.5 %
#> (Intercept) -1.854175 -0.9428936
#> period2     -1.593339 -0.3913224
#> period3     -1.768869 -0.4884715
#> period4     -2.418072 -0.7425514
#> method: wald_asymptotic_from_rust_inference_table
#> status: available
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
  pirls_profiled = fixef(glmm_fit),
  joint_laplace  = fixef(glmm_joint),
  glmer          = lme4::fixef(glmer_fit)
), 4)
#>                (Intercept) period2 period3 period4
#> pirls_profiled     -1.3605 -0.9762 -1.1111 -1.5597
#> joint_laplace      -1.3985 -0.9923 -1.1287 -1.5803
#> glmer              -1.3983 -0.9919 -1.1282 -1.5797
```

## Quadrature Sensitivity

`nAGQ` is part of the fit request and is recorded on the object. Values
above one are a profiled-path sensitivity check; the joint-Laplace path
is restricted to `nAGQ = 1`.

``` r

glmm_fit_agq3 <- glmm(
  cbind_formula,
  cbpp,
  family = cbpp_family,
  nAGQ = 3L,
  control = mm_control(verbose = -1)
)

data.frame(
  nAGQ = c(glmm_fit$nAGQ, glmm_fit_agq3$nAGQ),
  logLik = c(as.numeric(logLik(glmm_fit)),
             as.numeric(logLik(glmm_fit_agq3))),
  AIC = c(AIC(glmm_fit), AIC(glmm_fit_agq3)),
  check.names = FALSE
)
#>   nAGQ   logLik      AIC
#> 1    1 -92.0543 194.1086
#> 2    3 -92.0451 194.0902
```

On this dataset the quadrature order shifts the objective in the third
decimal and the coefficients by less than 0.005 on the logit scale —
visible, but an order of magnitude below the joint-Laplace standard
errors shown earlier.

## What is refused where

The table below runs each verb on the default fit and reports what
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

coef_tab <- summary(glmm_fit, tests = "coefficients")$coefficients
rbind(
  predict           = glmm_status(predict(glmm_fit)),
  fitted            = glmm_status(fitted(glmm_fit)),
  coefficient_tests = glmm_status(coef_tab,
                                  has_numbers = any(is.finite(coef_tab[, "Std. Error"]))),
  confint           = glmm_status(confint(glmm_fit)),
  simulate          = glmm_status(stats::simulate(glmm_fit, nsim = 1L)),
  refit             = glmm_status(refit(glmm_fit, fitted(glmm_fit)))
)
#>                                        status                    class
#> predict                             available                     <NA>
#> fitted                              available                     <NA>
#> coefficient_tests returns, inference withheld                     <NA>
#> confint                           typed error mm_inference_unavailable
#> simulate                          typed error mm_inference_unavailable
#> refit                             typed error mm_inference_unavailable
```

[`confint()`](https://rdrr.io/r/stats/confint.html) on the default fit
raises `mm_inference_unavailable` (the working-Hessian payload is not
certified for Wald intervals); the same request on the joint-Laplace fit
succeeds, as shown above.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) and
[`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md) are
not implemented for GLMM fits.
