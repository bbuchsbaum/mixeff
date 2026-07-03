# GLMM Fitting and Model Comparison

``` r

library(mixeff)
```

This vignette shows the current GLMM contract in `mixeff`. The default
path is `glmm(..., method = "pirls_profiled")`, which delegates to the
upstream profiled PIRLS GLMM fitter. `method = "joint_laplace"` requests
the labelled joint route (`fast = FALSE`, `nAGQ = 1`) and is available
in the dependency-light vendored build through the native optimizer.

The example uses the
[`lme4::cbpp`](https://rdrr.io/pkg/lme4/man/cbpp.html) data.
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) does
not yet expose binomial case weights, so the grouped counts are expanded
to Bernoulli rows before fitting.

``` r

env <- new.env(parent = emptyenv())
utils::data("cbpp", package = "lme4", envir = env)
cbpp <- get("cbpp", envir = env, inherits = FALSE)

expand_cbpp <- function(data) {
  rows <- lapply(seq_len(nrow(data)), function(i) {
    successes <- data$incidence[[i]]
    failures <- data$size[[i]] - successes
    data.frame(
      herd = data$herd[[i]],
      period = data$period[[i]],
      y = c(rep.int(1L, successes), rep.int(0L, failures))
    )
  })
  droplevels(do.call(rbind, rows))
}

cbpp_binary <- expand_cbpp(cbpp)
c(n_grouped_rows = nrow(cbpp), n_binary_rows = nrow(cbpp_binary))
#> n_grouped_rows  n_binary_rows 
#>             56            842
```

## Family and Audit

The statistical intent is a binomial-logit GLMM with a herd random
intercept. Before fitting, the same formula can be compiled and
explained so the fixed and random-effect structure is visible before
numerical optimization starts.

``` r

cbpp_family <- binomial(link = "logit")
cbpp_formula <- y ~ period + (1 | herd)

cbpp_spec <- compile_model(cbpp_formula, cbpp_binary)
explain_model(cbpp_spec)
#> Random effects explanation:
#>   formula: y ~ 1 + period + (1 | herd)
#> 
#> Random effects:
#>   r0:
#>     wrote:      (1 | herd)
#>     canonical:  (1 | herd)
#>     named form: re(group = herd, intercept = TRUE, slopes = NULL, cov = "scalar")
#>     scope:      `herd` units may differ in average outcome.
#>     covariance: scalar; theta parameters: 1
#>     support:    sufficient; group levels: 15; min rows/group: 26; median rows/group: 61
#>     variation:  intercept=not_assessed
```

## Fit

``` r

glmm_fit <- glmm(
  cbpp_formula,
  cbpp_binary,
  family = cbpp_family,
  method = "pirls_profiled",
  nAGQ = 1L,
  control = mm_control(verbose = -1)
)

glmm_fit
#> Generalized linear mixed model fit
#> Formula: y ~ period + (1 | herd)
#> Family/link: binomial/logit
#> Method: pirls_profiled (nAGQ = 1)
#> Fit status: converged_interior
#> Optimizer: cobyla; iterations: 19; objective: 555.06
#> Artifact: mixedmodels.compiled_model_artifact v1; crate: 1.0.0-rc.1
#> nobs: 842, dispersion: 1, logLik: -277.53
#> Fixed effects:
#> (Intercept)   period: 2   period: 3   period: 4 
#>   -1.360480   -0.976172   -1.111070   -1.559670 
#> Audit verbs: audit(), diagnostics(), model_report()
fixef(glmm_fit)
#> (Intercept)   period: 2   period: 3   period: 4 
#>  -1.3604769  -0.9761719  -1.1110715  -1.5596742
VarCorr(glmm_fit)
#> Variance components:
#>  group        name variance  std_dev correlation
#>   herd (Intercept) 0.411966 0.641846
```

The summary prints fitted quantities and deliberately avoids coefficient
tests unless the inference contract certifies them.

``` r

summary(glmm_fit)
#> Generalized linear mixed model fit
#> Formula: y ~ period + (1 | herd)
#> Family/link: binomial/logit
#> Method: pirls_profiled (nAGQ = 1)
#> Fit status: converged_interior
#> 
#> Variance components:
#>  group        name variance  std_dev correlation
#>   herd (Intercept) 0.411966 0.641846            
#> 
#> Fixed effects:
#>               Estimate Std. Error statistic p.value       method
#> (Intercept) -1.3604769         NA        NA      NA not_computed
#> period: 2   -0.9761719         NA        NA      NA not_computed
#> period: 3   -1.1110715         NA        NA      NA not_computed
#> period: 4   -1.5596742         NA        NA      NA not_computed
#> 
#> Wald-z reliability: not_available (not_computed). Reason: certified GLMM fixed-effect Wald inference is not implemented for fast_pirls_profiled; fast-PIRLS/profiled covariance geometry remains a working-Hessian payload, while only joint-laplace fits with a passing certified active-subspace Hessian over active beta plus interior theta parameters can report Wald SE/z/p/confint.
#> 
#> Notes:
#>   test statistics and p-values are withheld: the fit's covariance payload does not certify fixed-effect inference. Engine-certified Wald inference is available from a fit with method = "joint_laplace".
```

Core extractors are available from the durable R object.

``` r

head(fitted(glmm_fit))
#>         1         2         3         4         5         6 
#> 0.3092686 0.3092686 0.3092686 0.3092686 0.3092686 0.3092686
head(residuals(glmm_fit))
#>          1          2          3          4          5          6 
#>  0.6907314  0.6907314 -0.3092686 -0.3092686 -0.3092686 -0.3092686
head(ranef(glmm_fit)[[1L]])
#>   (Intercept)
#> 1  0.55693613
#> 2 -0.32949225
#> 3  0.36898286
#> 4  0.01190657
#> 5 -0.22339296
#> 6 -0.43118078
c(
  logLik = as.numeric(logLik(glmm_fit)),
  deviance = deviance(glmm_fit),
  AIC = AIC(glmm_fit),
  BIC = BIC(glmm_fit)
)
#>    logLik  deviance       AIC       BIC 
#> -277.5300  555.0599  565.0599  588.7388
```

## Quadrature Sensitivity

`nAGQ` is part of the fit request and is recorded on the object. Values
above one are a profiled-path sensitivity check; the labelled
joint-Laplace path is restricted to `nAGQ = 1`.

``` r

glmm_fit_agq3 <- glmm(
  cbpp_formula,
  cbpp_binary,
  family = cbpp_family,
  method = "pirls_profiled",
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
#>   nAGQ    logLik      AIC
#> 1    1 -277.5300 565.0599
#> 2    3 -277.5208 565.0415
```

## Current Boundaries

In-sample [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) are available,
and `summary(fit, tests = "coefficients")` returns a Wald-z fixed-effect
table built from the upstream PIRLS/Laplace working-Hessian covariance
payload (reliability label: `moderate`). GLMM
[`predict()`](https://rdrr.io/r/stats/predict.html) (population and
conditional, on the link or response scale),
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
[`drop1()`](https://rdrr.io/r/stats/add1.html), and multi-model
[`anova()`](https://rdrr.io/r/stats/anova.html) are available as R-side
computations on top of that certified covariance.
[`confint()`](https://rdrr.io/r/stats/confint.html) (Wald) requires a
joint-Laplace fit: the fast-PIRLS/profiled covariance is a
working-Hessian payload and the engine does not certify Wald intervals
for it. [`simulate()`](https://rdrr.io/r/stats/simulate.html) and
[`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md) for
GLMMs are not yet implemented and are exposed as ordinary R conditions
rather than silently filled in. The joint-Laplace fit route is available
and labelled separately from the profiled fast-PIRLS default.

``` r

glmm_boundary <- function(expr) {
  cnd <- tryCatch({
    force(expr)
    NULL
  }, error = function(cnd) cnd)
  if (is.null(cnd)) {
    return(data.frame(status = "available", class = NA_character_,
                      message = "", check.names = FALSE))
  }
  data.frame(
    status = "unavailable",
    class = class(cnd)[[1L]],
    message = conditionMessage(cnd),
    check.names = FALSE
  )
}

rbind(
  predict = glmm_boundary(predict(glmm_fit)),
  confint = glmm_boundary(confint(glmm_fit)),
  coefficient_tests = glmm_boundary(summary(glmm_fit, tests = "coefficients")),
  simulate = glmm_boundary(stats::simulate(glmm_fit, nsim = 1L)),
  refit = glmm_boundary(refit(glmm_fit, fitted(glmm_fit))),
  joint_laplace = glmm_boundary(glmm(
    cbpp_formula,
    cbpp_binary,
    family = cbpp_family,
    method = "joint_laplace",
    control = mm_control(verbose = -1)
  ))
)
#>                        status                    class
#> predict             available                     <NA>
#> confint           unavailable mm_inference_unavailable
#> coefficient_tests   available                     <NA>
#> simulate          unavailable              simpleError
#> refit             unavailable              simpleError
#> joint_laplace       available                     <NA>
#>                                                                                                                                                                                                                                                                                                                                                                                                               message
#> predict                                                                                                                                                                                                                                                                                                                                                                                                              
#> confint           `confint(method = "wald")` is unavailable for this GLMM fit because certified GLMM fixed-effect Wald inference is not implemented for fast_pirls_profiled; fast-PIRLS/profiled covariance geometry remains a working-Hessian payload, while only joint-laplace fits with a passing certified active-subspace Hessian over active beta plus interior theta parameters can report Wald SE/z/p/confint
#> coefficient_tests                                                                                                                                                                                                                                                                                                                                                                                                    
#> simulate                                                                                                                                                                                                                                                                                                    no applicable method for 'simulate' applied to an object of class "c('mm_glmm', 'mm_fit', 'mm_compiled')"
#> refit                                                                                                                                                                                                                                                                                                          no applicable method for 'refit' applied to an object of class "c('mm_glmm', 'mm_fit', 'mm_compiled')"
#> joint_laplace
```
