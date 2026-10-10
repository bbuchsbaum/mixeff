# Backward elimination of random and fixed effects (lmerTest's `step()`)

[`step()`](https://bbuchsbaum.github.io/mixeff/reference/step.md) for
`mm_lmm` fits follows
[`lmerTest::step()`](https://rdrr.io/pkg/lmerTest/man/step.html): it
first removes random-effect terms one at a time while the least
significant has `p > alpha.random`, then removes fixed-effect terms
(respecting marginality: a term is only dropped when no higher-order
term contains it) while the least significant has `p > alpha.fixed`,
refitting after every removal.

## Usage

``` r
# S3 method for class 'mm_lmm'
step(
  object,
  ddf = c("Satterthwaite", "Kenward-Roger"),
  alpha.random = 0.1,
  alpha.fixed = 0.05,
  reduce.fixed = TRUE,
  reduce.random = TRUE,
  keep,
  ...
)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- ddf:

  Denominator degrees of freedom for the fixed-effect F tests.

- alpha.random, alpha.fixed:

  Elimination thresholds.

- reduce.fixed, reduce.random:

  Whether to eliminate fixed / random terms.

- keep:

  Character vector of fixed-effect terms never to remove.

- ...:

  Unused.

## Value

An object of class `mm_step` with the elimination tables `random` and
`fixed` and the final model in `$model` (also returned by
[`lmerTest::get_model()`](https://rdrr.io/pkg/lmerTest/man/get_model.html)
when lmerTest is loaded).

## Details

Differences from lmerTest, by design: random terms are tested with
[`test_random_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_random_effect.md)'s
boundary-corrected likelihood-ratio test (a chi-bar-square mixture, so
p-values are about half of lmerTest's `ranova()` naive chi-square); only
whole random terms are removed (lmerTest's `ranova()` also tries
reducing `(x | g)` to `(1 | g)` – write the reduced term as its own
term, e.g. `(1 | g) + (0 + x | g)`, to make that reduction available);
and the last random term is never removed (a model without random
effects is not an LMM). Fixed terms are tested with
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md)
(`ddf` = `"Satterthwaite"` or `"Kenward-Roger"`).

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:10, each = 6)), x = rep(0:5, 10),
                z = rnorm(60))
d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(60)
fit <- lmm(y ~ x + z + (1 | g) + (0 + x | g), d,
           control = mm_control(verbose = -1))
s <- step(fit)
s
#> Backward reduced random-effect table:
#>         term statistic   p_value                          method eliminated
#>  (0 + x | g) 0.4933364 0.2412214 boundary_lrt_self_liang_mixture          1
#> 
#> Backward reduced fixed-effect table (ddf: Satterthwaite):
#>  term num_df   den_df     F_value      p_value eliminated
#>     z      1 48.60119  0.03294079 8.567331e-01          1
#>     x      1 48.99997 49.89804401 5.282960e-09          0
#> 
#> Model found:
#> y ~ x + (1 | g)
s$model
#> Linear mixed model fit by REML
#> Formula: y ~ x + (1 | g)
#> Fit status: converged_interior
#> Optimizer: pattern_search; iterations: 24; objective: 170.46
#> nobs: 60, sigma: 0.793663, logLik: -85.2301
#> Fixed effects:
#> (Intercept)           x 
#>    1.677290    0.423798 
#> Audit verbs: audit(), diagnostics(), inference_table(), model_report()
```
