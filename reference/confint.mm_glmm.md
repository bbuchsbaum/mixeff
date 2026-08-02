# Confidence intervals for fixed effects of a mixeff GLMM

Two routes are available. `method = "asymptotic"` builds Wald intervals
(`estimate +/- z * SE`) from the Rust fixed-effect inference table; it
requires a certified fit (`method = "joint_laplace"`) or the explicit
working-Hessian opt-in, and refuses the default profiled estimator with
a typed reason. `method = "bootstrap"` runs a parametric bootstrap —
each replicate simulates a response from the fitted model under fresh
random-effect draws and refits the same (effective) estimator — and
returns percentile intervals; it works on profiled-estimator fits for
every supported family, including negative binomial (fixed-theta fits
condition replicates on that theta; estimated-theta fits re-estimate
theta per replicate, propagating its uncertainty). On `joint_laplace`
fits the bootstrap is refused at this engine pin — replicate refits
cannot re-run the joint estimator — and those fits have certified Wald
intervals instead. The bootstrap result carries replicate accounting as
attributes: `mm_bootstrap` (requested, successful, failed, seed,
per-parameter standard errors, Monte Carlo SE, and a reliability grade)
and `mm_estimator`. Failed replicates are counted, never silently
dropped. Profile intervals are not certified for GLMMs and are refused.

## Usage

``` r
# S3 method for class 'mm_glmm'
confint(
  object,
  parm,
  level = 0.95,
  method = c("asymptotic", "wald", "profile", "bootstrap"),
  ...
)
```

## Arguments

- object:

  A fitted `mm_glmm`.

- parm:

  Optional fixed-effect names or indices; defaults to all.

- level:

  Confidence level.

- method:

  `"asymptotic"` (the default; the package-wide name for the closed-form
  Wald interval), its synonym `"wald"`, or `"bootstrap"`.

- ...:

  For `method = "bootstrap"`: `nsim` (number of replicates, default 999)
  and `seed` (non-negative integer; when omitted one is drawn from R's
  RNG so [`set.seed()`](https://rdrr.io/r/base/Random.html) governs
  reproducibility).

## Value

An `mm_confint` matrix of lower/upper bounds (`"asymptotic"`), or a
plain matrix of percentile bounds with bootstrap-accounting attributes
(`"bootstrap"`).
