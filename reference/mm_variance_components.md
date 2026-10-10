# Variance decomposition, Nakagawa R2 and ICC for mixeff fits

`mm_variance_components()` splits the variance of a fitted mixed model
into the components used by Nakagawa and Schielzeth's R2 and the
intraclass correlation: the variance of the fixed-effect linear
predictor, the random-effect variance, and the residual
(distribution-specific plus observation-level) variance. `mm_r2()`
returns the marginal and conditional R2 and `mm_icc()` the adjusted and
unadjusted ICC. The numbers follow
[`insight::get_variance()`](https://easystats.github.io/insight/reference/get_variance.html),
[`performance::r2_nakagawa()`](https://easystats.github.io/performance/reference/r2_nakagawa.html)
and
[`performance::icc()`](https://easystats.github.io/performance/reference/icc.html)
on the equivalent lme4 fit.

## Usage

``` r
mm_variance_components(fit, tolerance = 1e-05)

mm_r2(fit, tolerance = 1e-05)

mm_icc(fit, tolerance = 1e-05)
```

## Arguments

- fit:

  A fitted `mm_lmm` or `mm_glmm`.

- tolerance:

  Variances below this value mark the fit as singular.

## Value

`mm_variance_components()`: a named list with `var.fixed`, `var.random`,
`var.residual`, `var.distribution`, `var.dispersion`, `var.intercept`,
`var.slope` and `cor.slope_intercept` (the layout of
[`insight::get_variance()`](https://easystats.github.io/insight/reference/get_variance.html)).
`mm_r2()`: a list with `R2_conditional` and `R2_marginal`. `mm_icc()`: a
one-row data frame with `ICC_adjusted`, `ICC_conditional` and
`ICC_unadjusted`.

## Details

When the insight and performance packages are installed, mixeff also
registers methods so that `performance::icc(fit)`,
`performance::r2(fit)` and `insight::get_variance(fit)` work directly on
mixeff fits (not `by_group = TRUE`, which needs lme4 internals).

- Fixed variance: `var(X beta)`.

- Random variance: for each random term whose grouping factor has fewer
  levels than observations, the mean over observations of
  `z_i' Sigma z_i`, where `z_i` holds the term's coefficients' columns
  of the fixed-effect design (a random slope without a matching fixed
  column is dropped, as insight does).

- Distribution-specific variance: `sigma^2` for an LMM; `pi^2/3`, `1`
  and `pi^2/6` for binomial logit, probit and cloglog links; for Gamma,
  `sigma^2` (inverse link) or the log-normal approximation
  `log1p(sigma^2)` (log link), with
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) the lme4 \>= 2.1-0
  dispersion `sqrt(phi)`; and the log-normal approximation
  `log(1 + V(mu)/mu^2)` for Poisson (`V(mu) = mu`) and negative-binomial
  (`V(mu) = mu (1 + mu/theta)`) log-link models, with
  `mu = exp(b0 + v0/2)` from the intercept `b0` and the random-effect
  variance `v0` (observation-level terms excluded) of the null model
  (`y ~ 1` with the same random effects), which is refitted. A
  `cbind(successes, failures)` binomial response divides the binomial
  value by the mean number of trials. A Gaussian GLMM with a
  non-identity link uses `sigma^2`, as insight does for every Gaussian
  model: note that this is the residual variance on the response scale,
  while the fixed and random variances are on the link scale, so the R2
  and ICC of such a fit depend on the units of the response.
  Observation-level terms of a Gaussian GLMM do not add to the residual
  variance (insight treats the model as linear). These follow insight
  1.5.

- Inverse-Gaussian GLMMs are refused (`mm_inference_unavailable`, reason
  code `"r2_distribution_variance_undefined"`). insight 1.5 has no
  inverse-Gaussian branch and falls back to
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) itself (the square
  root of the dispersion, not a variance); since the inverse-Gaussian
  dispersion has units of `1/y`, that value, and the R2 and ICC built on
  it, change when the response is rescaled.
  [`insight::get_variance()`](https://easystats.github.io/insight/reference/get_variance.html)
  on a mixeff inverse-Gaussian fit therefore returns `NA` with a
  warning.

- Observation-level random effects (one level per observation) count
  towards the residual variance of a GLMM.

If a random-effect variance is below `tolerance` (a singular fit), the
random variance is not computed: `mm_r2()` and `mm_icc()` return `NA`
for the quantities that need it and signal a warning, as performance
does.

## Examples

``` r
set.seed(1)
g <- factor(rep(seq_len(10), each = 6))
df <- data.frame(x = rnorm(60), g = g)
df$y <- 1 + 0.5 * df$x + rnorm(10)[g] + rnorm(60)
fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
mm_r2(fit)
#> $R2_conditional
#> Conditional R2 
#>      0.7028025 
#> 
#> $R2_marginal
#> Marginal R2 
#>   0.1042069 
#> 
mm_icc(fit)
#>   ICC_adjusted ICC_conditional ICC_unadjusted
#> 1    0.6682297       0.5985956      0.5985956
```
