# Simulate responses from a mixeff fit

Draws responses from a fitted `mm_lmm` or `mm_glmm` following the
semantics of
[`lme4::simulate.merMod()`](https://rdrr.io/pkg/lme4/man/simulate.merMod.html):

## Usage

``` r
# S3 method for class 'mm_lmm'
simulate(object, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA, ...)

# S3 method for class 'mm_glmm'
simulate(object, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA, ...)
```

## Arguments

- object:

  A fitted `mm_lmm` or `mm_glmm`.

- nsim:

  Number of simulated responses.

- seed:

  Optional random seed.

- use.u:

  Logical; `TRUE` is the same as `re.form = NULL` and `FALSE` the same
  as `re.form = NA`. Specify at most one of `use.u` and `re.form`.

- re.form:

  `NA` or `~0` (default `NA`) to simulate new random effects, `NULL` to
  condition on the fitted random effects, or a formula such as
  `~ (1 | g)` to condition on those terms' fitted random effects and
  draw new ones for the others (as lme4).

- ...:

  Reserved; lme4 arguments that cannot be honoured (`newdata`,
  `newparams`, ...) are refused with a typed error.

## Value

A data frame with one column per simulation (`sim_1`, ...), rows named
like the model frame. For a binomial
[`cbind()`](https://rdrr.io/r/base/cbind.html) fit each column is a
two-column matrix. Attributes `seed`, `mm_method` and `mm_re_form`
record the seed, the simulation route, and the resolved conditioning
(`"unconditional"`, `"conditional"` or `"partial"`). Following
[`stats::simulate`](https://rdrr.io/r/stats/simulate.html), attribute
`"seed"` is the `.Random.seed` in force before simulation when
`seed = NULL`, otherwise `seed` with the RNG kind as attribute `"kind"`.

## Details

- `re.form = NA` (the default) or `~0`, equivalently `use.u = FALSE`,
  draws **new** random effects from the fitted random-effect covariance
  for every simulation (unconditional, parametric-bootstrap simulation);

- `re.form = NULL`, equivalently `use.u = TRUE`, conditions on the
  fitted conditional modes (BLUPs) and simulates only the response noise
  around `fitted(object)`.

Partial random-effect formulas are refused with a typed
`mm_inference_unavailable` error.

LMM responses are Gaussian with standard deviation `sigma / sqrt(w)` for
prior weights `w` (lme4 ignores the prior weights here). GLMM responses
are family-aware: binomial draws use the prior weights as trial counts
and return proportions (`y / w`, 0/1 values for Bernoulli fits), or a
two-column `(successes, failures)` matrix per simulation when the model
was fitted with a [`cbind()`](https://rdrr.io/r/base/cbind.html)
response, as lme4 does; Poisson and negative-binomial (fitted `theta`)
draws return counts. For the families with a free dispersion
`phi = sigma(object)^2` the draws follow the fitted variance function:
Gamma draws use shape `w / phi` (`Var(y) = phi mu^2 / w`),
inverse-Gaussian draws use
[`statmod::rinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
with shape `w / phi` (`Var(y) = phi mu^3 / w`; needs the statmod
package), and Gaussian non-identity-link draws have standard deviation
`sigma / sqrt(w)`. New random effects are drawn on the absolute scale of
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
as lme4 \>= 2.1-0 does (its theta is the absolute random-effect SD for
every GLMM).

Compared with lme4 2.1-0's `simulate.merMod()`: Gamma draws and Gaussian
draws without prior weights agree draw for draw (same RNG use) when the
estimates agree. lme4 ignores prior weights for Gaussian responses (with
a warning); mixeff uses them, as for LMMs. For the inverse Gaussian lme4
passes shape `w / sigma(object)` – `sqrt(phi)`, not `phi` – so its draws
have variance `sigma mu^3 / w` rather than the fitted `phi mu^3 / w`;
mixeff uses the fitted variance, as `stats::inverse.gaussian()$simulate`
does. A binomial fit whose response was a factor is simulated as 0/1
numeric (lme4 returns a factor).

The random-number stream is consumed in lme4's order (all random-effect
draws first, terms in lme4's internal order, then the response draws),
so when the estimates agree a seeded
[`simulate()`](https://rdrr.io/r/stats/simulate.html) reproduces lme4's
draws. Unlike lme4, a supplied `seed` does not leave the global RNG
state changed.

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rnorm(60), x = rnorm(60),
  g = factor(rep(seq_len(10), each = 6))
)
fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
head(simulate(fit, nsim = 2, seed = 1))
#>         sim_1        sim_2
#> 1  0.80406962 -0.477766041
#> 2  0.78783907 -0.002568654
#> 3  0.14779051  1.098644257
#> 4 -1.60302473 -1.201679482
#> 5  0.67731157  0.655012358
#> 6  0.05604579  0.391297220
# condition on the fitted random effects
head(simulate(fit, nsim = 1, seed = 1, re.form = NULL))
#>        sim_1
#> 1 -0.5275568
#> 2  0.2721454
#> 3 -0.6364811
#> 4  1.4856869
#> 5  0.4271580
#> 6 -0.6025502

df$count <- rpois(60, 2)
gfit <- glmm(count ~ x + (1 | g), df, family = poisson(),
             control = mm_control(verbose = -1))
head(simulate(gfit, nsim = 2, seed = 1))
#>   sim_1 sim_2
#> 1     2     2
#> 2     2     1
#> 3     2     1
#> 4     2     5
#> 5     2     2
#> 6     2     1
```
