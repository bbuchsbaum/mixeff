# The supported model envelope

`supported_models()` returns the package's machine-readable support
registry: which GLMM family/link/estimator combinations are supported,
which features are supported, refused, or carry mixeff-specific
semantics, and the package's explicit non-goals. The registry ships as
`inst/support/model-support.json` and is the single source of truth for
the migration vignette's support table; the test suite asserts that the
fitting code enforces exactly this envelope (unsupported combinations
refuse with the declared condition classes and reason codes).

## Usage

``` r
supported_models()
```

## Value

A list of class `mm_model_support` with elements `glmm_families` (data
frame: family, links, estimators, notes), `features` (data frame:
feature, lmm, glmm, notes), `non_goals` (character), and `registry` (the
parsed JSON, verbatim).

## Details

mixeff aims for tested statistical agreement with lme4 on this envelope
— a documented, closed surface — not blanket equivalence with everything
`lmer()`/`glmer()` accept.

## Examples

``` r
supported_models()
#> mixeff supported model envelope (registry mixeff.model_support v1.0.0)
#> 
#> GLMM families:
#>  family            links                  estimators                   
#>  binomial          logit, probit, cloglog pirls_profiled, joint_laplace
#>  poisson           log, sqrt              pirls_profiled, joint_laplace
#>  Gamma             log                    pirls_profiled, joint_laplace
#>  negative_binomial log                    pirls_profiled               
#>  notes                                                                                                 
#>  Bernoulli 0/1, cbind(successes, failures), and proportion + trial weights responses.                  
#>                                                                                                        
#>                                                                                                        
#>  NB2 via mm_negative_binomial(); theta estimated or fixed. Intervals via confint(method = "bootstrap").
#> 
#> Features:
#>  feature                                                                                  
#>  case weights                                                                             
#>  offset                                                                                   
#>  simulate                                                                                 
#>  refit                                                                                    
#>  subset / random / custom na.action / custom contrasts (GLMM)                             
#>  marginal means verbs (mm_means / mm_comparisons / mm_grid / mm_predictions / test_effect)
#>  double-bar || with factor terms                                                          
#>  lmm              glmm            
#>  supported        supported       
#>  unsupported      supported       
#>  supported        refused         
#>  supported        refused         
#>  see notes        refused         
#>  supported        refused         
#>  mixeff semantics mixeff semantics
#> 
#> Non-goals: nonlinear mixed models; arbitrary GLM family objects or custom links; adaptive Gaussian quadrature as a general estimator contract (nAGQ > 1 is a profiled-path sensitivity mode) 
#> 
#> Gaussian LMMs support the common random-intercept, random-slope,
#> nested, and crossed structures; see vignette("lme4-migration") for
#> the argument-by-argument map and the feature notes column for
#> mixeff-specific semantics.
```
