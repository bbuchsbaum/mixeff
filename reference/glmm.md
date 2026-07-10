# Fit a generalized linear mixed model

`glmm()` validates the R-side family/link request, compiles the model
formula, and delegates the numerical fit to the upstream Rust
`GeneralizedLinearMixedModel`. The default `method = "pirls_profiled"`
is the labelled fast-PIRLS path. `method = "joint_laplace"` uses the
upstream labelled joint Laplace route (`fast = FALSE`, `nAGQ = 1`)
backed by the native dependency-light optimizer in this vendored build.

## Usage

``` r
glmm(
  formula,
  data,
  family,
  random = NULL,
  weights = NULL,
  offset = NULL,
  subset = NULL,
  na.action = na.omit,
  contrasts = NULL,
  method = c("pirls_profiled", "joint_laplace"),
  nAGQ = 1L,
  inference = c("auto", "none", "asymptotic", "bootstrap"),
  control = mm_control(),
  ...
)
```

## Arguments

- formula:

  A two-sided lme4-style formula.

- data:

  A `data.frame`.

- family:

  A supported GLMM family object or family constructor. The certified
  1.0 surface is: [`binomial()`](https://rdrr.io/r/stats/family.html)
  with `"logit"`, `"probit"`, or `"cloglog"` links;
  [`poisson()`](https://rdrr.io/r/stats/family.html) with `"log"` or
  `"sqrt"` links; [`Gamma()`](https://rdrr.io/r/stats/family.html) with
  `"log"` link; and negative binomial (NB2, `"log"` link) via
  [`mm_negative_binomial()`](https://bbuchsbaum.github.io/mixeff/reference/mm_negative_binomial.md)
  (theta estimated, like
  [`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html)) or
  `MASS::negative.binomial(theta)` (fixed theta).

- random:

  Reserved for the native random-effect constructor path.

- weights:

  Optional prior weights. For binomial models these are trial counts for
  proportion responses; weights must be positive and finite.

- offset:

  Optional fixed linear-predictor offset; values must be finite.

- subset, na.action, contrasts:

  Reserved for future parity with
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md).

- method:

  GLMM estimation method. `"pirls_profiled"` is the default fast-PIRLS
  profiled path. `"joint_laplace"` requests the labelled joint Laplace
  route and requires `nAGQ <= 1`. The joint route tracks the lme4
  joint-Laplace reference far more closely than the profiled path on
  high-baseline models, at a higher optimizer cost; cap that cost with
  `mm_control(max_feval = )`. The default profiled path is **not**
  glmer's estimator and its coefficients do not match `glmer()` exactly;
  when `method` is left at its default, `glmm()` emits an informational
  notice to that effect (suppress with `mm_control(verbose = -1)`). Use
  `method = "joint_laplace"` for glmer-equivalent estimates.

- nAGQ:

  Number of adaptive Gauss-Hermite quadrature points. `1` is the Laplace
  setting. Values above `1` are allowed on the profiled path and are
  rejected for `method = "joint_laplace"` in the R wrapper.

- inference:

  Requested inference method.

- control:

  A list from
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md).

- ...:

  Reserved for future use.

## Value

An object of class `mm_glmm`, also inheriting from `mm_fit` and
`mm_compiled`.

## Details

Optimization runs inside a single native call with no progress output:
the pre-fit explanation block (when `verbose >= 0`) is the last thing
printed before the fitted result returns, and the call cannot be
interrupted from R. Every optimizer budget is bounded, so fits always
terminate; runtime on large problems is governed by
`mm_control(max_feval = )`.

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rbinom(120, 1, 0.5),
  x = rnorm(120),
  g = factor(rep(seq_len(12), each = 10))
)
fit <- glmm(y ~ x + (1 | g), df, family = binomial(),
            control = mm_control(verbose = -1))
fixef(fit)
#> (Intercept)           x 
#> -0.10067766  0.01797541 
# glmer-equivalent (joint Laplace) estimates:
fit_joint <- glmm(y ~ x + (1 | g), df, family = binomial(),
                  method = "joint_laplace",
                  control = mm_control(verbose = -1))
fixef(fit_joint)
#> (Intercept)           x 
#> -0.10067766  0.01797541 
```
