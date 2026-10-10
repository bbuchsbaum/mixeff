# Parametric bootstrap of any statistic of a fitted LMM

`mm_bootmer()` is mixeff's counterpart of
[`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html): it
simulates `nsim` responses from the fitted model, refits the model to
each, and evaluates `FUN` on every refit. The result has class
`c("mm_bootmer", "boot")` with the fields
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) reads
(`t0`, `t`, `R`, `sim`, `statistic`, ...), so
`boot::boot.ci(b, type = "perc", index = 1)` works as it does for a
`bootMer` result.

## Usage

``` r
mm_bootmer(
  fit,
  FUN,
  nsim = 1L,
  seed = NULL,
  type = c("parametric", "semiparametric"),
  use.u = FALSE,
  verbose = FALSE
)
```

## Arguments

- fit:

  A fitted `mm_lmm`. GLMM fits are refused (the GLMM parametric
  bootstrap is `confint(fit, method = "bootstrap")`).

- FUN:

  A function of a fitted `mm_lmm` returning a numeric vector (for
  example `fixef`, or `function(m) c(fixef(m), sigma = sigma(m))`).

- nsim:

  Number of bootstrap replicates.

- seed:

  Optional seed; `NULL` uses the current R RNG state, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) governs
  reproducibility.

- type:

  Only `"parametric"` is available. lme4's `"semiparametric"` residual
  bootstrap is refused.

- use.u:

  `FALSE` (default) draws new random effects for each replicate (lme4's
  default, the unconditional bootstrap); `TRUE` conditions on the fitted
  random effects and resamples only the residual noise.

- verbose:

  Print a progress dot per replicate.

## Value

An object of class `c("mm_bootmer", "boot")`: a list with `t0`
(`FUN(fit)`), `t` (an `nsim` x `length(t0)` matrix), `R`, `data`,
`seed`, `statistic`, `sim = "parametric"`, `call`, and `mle`, plus
attributes `bootFail` (failed replicates) and `boot.fail.msgs`.

## Details

The function is named `mm_bootmer()` rather than `bootMer()` because
[`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html) is not an
S3 generic: exporting a `bootMer()` would mask lme4's for `merMod` fits.

Replicates use the same model, `REML` setting, weights, and
[`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
as `fit`; each refit warm-starts from the fitted covariance parameters
(as lme4's
[`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md)
does). A replicate whose refit fails, or whose `FUN` errors, is recorded
as an `NA` row and counted in `attr(, "bootFail")` (lme4's convention) –
never silently dropped.

## See also

[`confint()`](https://rdrr.io/r/stats/confint.html) with
`method = "bootstrap"` for engine-run fixed-effect bootstrap intervals;
[`simulate.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/simulate.mm_lmm.md);
[`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md).

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:8, each = 5)), x = rep(0:4, 8))
d$y <- 1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(40)
fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
b <- mm_bootmer(fit, fixef, nsim = 20, seed = 2)
apply(b$t, 2, sd)
#> (Intercept)           x 
#>   0.3896998   0.0944960 
if (requireNamespace("boot", quietly = TRUE)) {
  boot::boot.ci(b, type = "perc", index = 2)
}
#> Warning: extreme order statistics used as endpoints
#> BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#> Based on 20 bootstrap replicates
#> 
#> CALL : 
#> boot::boot.ci(boot.out = b, type = "perc", index = 2)
#> 
#> Intervals : 
#> Level     Percentile     
#> 95%   ( 0.3645,  0.6936 )  
#> Calculations and Intervals on Original Scale
#> Warning : Percentile Intervals used Extreme Quantiles
#> Some percentile intervals may be unstable
```
