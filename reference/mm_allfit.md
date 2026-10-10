# Refit a model with every available optimizer

`mm_allfit()` is mixeff's counterpart of
[`lme4::allFit()`](https://rdrr.io/pkg/lme4/man/allFit.html): it refits
`fit` once per optimizer (through
[`update()`](https://rdrr.io/r/stats/update.html) with
`mm_control(optimizer = )`, keeping the rest of the fit's control) and
summarises the fixed effects, covariance parameters, log-likelihood,
convergence status and run time side by side, so a fit's sensitivity to
the optimizer is visible. An optimizer that is not compiled into this
build, or whose refit fails, is reported as a failed row with its error
message – never dropped.

## Usage

``` r
mm_allfit(fit, optimizers = mm_allfit_optimizers(), verbose = FALSE)
```

## Arguments

- fit:

  A fitted `mm_lmm` or `mm_glmm`.

- optimizers:

  Character vector of
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
  optimizer names. The default tries every name
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
  accepts.

- verbose:

  Print each optimizer as it runs.

## Value

An object of class `mm_allfit`: a named list of fits (or conditions for
failed refits). [`summary()`](https://rdrr.io/r/base/summary.html)
returns a list with `which.OK`, `msgs`, `fixef`, `theta`, `llik`,
`fit_status`, and `times`, mirroring `summary(lme4::allFit(...))`.

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:8, each = 5)), x = rep(0:4, 8))
d$y <- 1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(40)
fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
af <- mm_allfit(fit, optimizers = c("auto", "pattern_search", "cobyla"))
summary(af)$llik
#>           auto pattern_search         cobyla 
#>      -56.02412      -56.02412      -56.02412 
```
