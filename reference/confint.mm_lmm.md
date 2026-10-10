# Confidence intervals for a mixeff LMM

`method = "asymptotic"` (the default; synonyms `"wald"` and lme4's
`"Wald"`) gives Wald intervals for the fixed effects from the stored
standard errors. `method = "profile"` gives profile-likelihood intervals
with lme4's rows and names: `.sig01`, `.sig02`, ... (random-effect
standard deviations and correlations, numbered as lme4 \>= 2.0 does:
term by term, each term's standard deviations before its correlations),
`.sigma`, then the fixed effects. As in lme4, every fit is profiled on
the ML deviance: a REML fit's intervals are those of its ML refit. The
engine's relative-Cholesky `theta1`, ... intervals (on the fit's own
criterion) are kept in `attr(ci, "mm_profile")$table` and can be
requested by name through `parm`. When one parameter's profile is
irregular (non-monotone, not bracketing the estimate, or a failed
refit), its undetermined bounds are `NA` – lme4 would warn and
interpolate linearly, often reporting `[-1, 1]` for a correlation – and
the reason is in `attr(ci, "mm_profile")$table` (`status`,
`reason_code`, `reason`); printing the result adds a note. See
[`profile.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/profile.mm_lmm.md).
`method = "bootstrap"` (lme4's `"boot"`) gives parametric-bootstrap
intervals for the fixed effects, controlled by
[`bootstrap_control()`](https://bbuchsbaum.github.io/mixeff/reference/bootstrap_control.md).

## Usage

``` r
# S3 method for class 'mm_lmm'
confint(
  object,
  parm,
  level = 0.95,
  method = c("asymptotic", "wald", "bootstrap", "profile"),
  bootstrap = NULL,
  interval = c("percentile", "basic"),
  threads = 1L,
  ...
)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- parm:

  Parameters to report: names (fixed-effect names, and for `"profile"`
  also `".sig01"`, ..., `".sigma"`, `"theta1"`, ...), lme4's `"theta_"`
  / `"beta_"` shortcuts (`"profile"` only), or indices. Defaults to all
  fixed effects (all lme4 rows for `"profile"`).

- level:

  Confidence level.

- method:

  `"asymptotic"`, `"wald"`, `"bootstrap"`, or `"profile"`.

- bootstrap:

  Optional
  [`bootstrap_control()`](https://bbuchsbaum.github.io/mixeff/reference/bootstrap_control.md)
  for `method = "bootstrap"` (its `threads` field sets the bootstrap's
  worker threads).

- interval:

  Bootstrap interval type, `"percentile"` or `"basic"`.

- threads:

  Worker threads for `method = "profile"` (default 1). The per-parameter
  profiles run concurrently; the result is identical for every value.
  The default is 1 because CRAN policy limits packages to at most two
  threads unless the user explicitly asks for more; passing `threads` is
  that request.

- ...:

  Unused.

## Value

An `mm_confint` matrix of lower/upper bounds.

## Details

lme4's `confint.merMod()` defaults to `method = "profile"`; mixeff keeps
the Wald default because the profile refits the model many times (about
the cost of lme4's profile; several seconds for a modest model) and can
refuse on near-singular fits, while the Wald route is free. Request
`method = "profile"` to reproduce lme4's default output.

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(60))
d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(60)
fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
confint(fit)
#> Confidence intervals:
#>                 2.5 %    97.5 %
#> (Intercept) 0.7432779 2.2179314
#> x           0.3087413 0.8061215
#> method: wald_asymptotic_from_stored_standard_errors
#> status: Wald (asymptotic) intervals from stored standard errors (engine-certified profile intervals: method = "profile")
confint(fit, method = "profile")
#> Confidence intervals:
#>                 2.5 %    97.5 %
#> .sig01      0.6951562 1.8545088
#> .sigma      0.6630609 0.9829297
#> (Intercept) 0.7087136 2.2526274
#> x           0.3069545 0.8087305
#> method: profile_likelihood
#> status: available
```
