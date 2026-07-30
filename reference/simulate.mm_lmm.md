# Simulate from a mixeff LMM

Draws Gaussian responses from the stored fixed effects, random-effect
covariance summaries, and residual scale. Simulation is implemented for
linear mixed models only; calling
[`simulate()`](https://rdrr.io/r/stats/simulate.html) on a GLMM fit
signals a typed `mm_inference_unavailable` error with reason code
`glmm_simulate_unimplemented`.

## Usage

``` r
# S3 method for class 'mm_lmm'
simulate(object, nsim = 1, seed = NULL, re.form = NULL, ...)

# S3 method for class 'mm_glmm'
simulate(object, nsim = 1, seed = NULL, re.form = NULL, ...)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- nsim:

  Number of simulated responses.

- seed:

  Optional random seed.

- re.form:

  Random-effects conditioning. `NULL` simulates new random effects; `NA`
  simulates from the population-level mean only.

- ...:

  Reserved for future methods.

## Value

A data frame of simulated responses.
