# Fit a generalized linear mixed model

`glmm()` validates the R-side family/link request, compiles the model
formula, and delegates the numerical fit to the upstream Rust
`GeneralizedLinearMixedModel`. The default `method = "joint_laplace"` is
the labeled joint Laplace route (`fast = FALSE`, `nAGQ = 1`), the
estimator [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)
uses by default, backed by the native dependency-light optimizer in this
vendored build. `method = "pirls_profiled"` is the labeled fast-PIRLS
profiled path (lme4's `nAGQ = 0`).

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
  na.action = getOption("na.action"),
  contrasts = NULL,
  method = c("joint_laplace", "pirls_profiled"),
  nAGQ = 1L,
  inference = c("auto", "none", "asymptotic", "bootstrap", "working_hessian"),
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

  A supported GLMM family object or family constructor. The supported
  surface (every family/link pair the engine fits) is:
  [`binomial()`](https://rdrr.io/r/stats/family.html) with `"logit"`,
  `"probit"`, or `"cloglog"` links;
  [`poisson()`](https://rdrr.io/r/stats/family.html) with `"log"` or
  `"sqrt"` links; [`Gamma()`](https://rdrr.io/r/stats/family.html) with
  `"inverse"` (R's default) or `"log"` links;
  [`inverse.gaussian()`](https://rdrr.io/r/stats/family.html) with
  `"inverse"` or `"log"` links (R's default `"1/mu^2"` link is not
  available and is refused);
  [`gaussian()`](https://rdrr.io/r/stats/family.html) with non-identity
  `"log"`, `"inverse"`, or `"sqrt"` links (a Gaussian identity-link
  model is an LMM: use
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)); and
  negative binomial (NB2, `"log"` link) via
  [`mm_negative_binomial()`](https://bbuchsbaum.github.io/mixeff/reference/mm_negative_binomial.md)
  (theta estimated, like
  [`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html)) or
  `MASS::negative.binomial(theta)` (fixed theta). Binomial `"cauchit"`,
  `"log"`, and `"identity"` and Poisson `"identity"` links are not
  available in the engine and are refused with a typed
  `mm_inference_unavailable` error naming the supported set.

- random:

  Reserved for the native random-effect constructor path.

- weights:

  Optional prior weights: a column of `data` (as in lme4,
  `weights = trials`) or a numeric vector. For binomial models these are
  trial counts for proportion responses; weights must be positive and
  finite.

- offset:

  Optional fixed linear-predictor offset: a column of `data` or a
  numeric vector; values must be finite.
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the formula
  are also supported and are added to it (as in
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html)).

- subset, na.action, contrasts:

  As in [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md):
  `subset` selects rows (evaluated in `data`), `na.action` controls
  missing values (default `getOption("na.action")`, i.e. `na.omit`, as
  in `glmer()`; only model variables, `weights` and `offset` count;
  dropped rows are announced with a typed `mm_rows_dropped` message and
  recorded in `na.action(fit)`; `na.exclude` pads
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`residuals()`](https://rdrr.io/r/stats/residuals.html)/[`predict()`](https://rdrr.io/r/stats/predict.html)
  back to the original rows; `na.fail` and `na.pass` are refused with a
  typed `mm_data_error`), and `contrasts` is honoured only when it names
  the engine's coding (`contr.treatment` for unordered, `contr.poly` for
  ordered factors).

- method:

  GLMM estimation method. `"joint_laplace"` (the default) is the joint
  Laplace route, matching `glmer()`'s default estimator, and requires
  `nAGQ = 1`; it certifies Wald standard errors, tests, and intervals.
  Its optimizer cost is higher than the profiled path's; cap it with
  `mm_control(max_feval = )`. `"pirls_profiled"` is the fast profiled
  PIRLS path (equivalent to lme4's `nAGQ = 0` fast estimate); its
  coefficients do not match `glmer(nAGQ = 1)` exactly and Wald inference
  is withheld. When `method` is not supplied and the request needs the
  profiled path – a negative-binomial family (the joint route is not
  available for it yet), `nAGQ > 1`, or `inference = "working_hessian"`
  – `glmm()` uses `"pirls_profiled"` and says so with an
  `mm_estimator_notice` message (silence with
  `mm_control(verbose = -1)`). An explicit `method = "joint_laplace"`
  with a negative-binomial family or `nAGQ != 1` is refused, never
  swapped.

- nAGQ:

  Number of adaptive Gauss-Hermite quadrature points. `1` (the default)
  is the Laplace approximation. `0` requests lme4's PIRLS-only fast
  estimate and selects `method = "pirls_profiled"` (refused together
  with an explicit `method = "joint_laplace"`). Values above `1` run on
  the profiled path only and, as in lme4, require a model with a single
  scalar random-effect term (e.g. `(1 | g)`); other models are refused.

- inference:

  Requested inference posture. The default `"auto"` keeps the certified
  contract: Wald standard errors, tests, and intervals are available
  only when the engine certifies them (currently
  `method = "joint_laplace"`, the default); the profiled estimator
  withholds them with a typed refusal. `"working_hessian"` is an
  explicit opt-in that unlocks the UNCERTIFIED profiled working-Hessian
  approximation on every inference route; each resulting row is labelled
  `wald_z_working_hessian` with reliability `moderate`. Its standard
  errors came within 1% of `glmer()`'s on the package's reference
  dataset, but nothing certifies that agreement in general, so treat it
  as an exploration and screening tool, not a reporting route; see
  [`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md).
  `"none"`, `"asymptotic"`, and `"bootstrap"` are accepted and recorded
  but currently equivalent to `"auto"`.

- control:

  A list from
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md).
  Besides the optimizer settings it carries lme4 2.1-0's
  `glmerControl()` dispersion controls (`disp_method`,
  `disp_dof_correction`, `max_phi_iter` for lme4's `maxPhiIter`), which
  affect Gamma, inverse-Gaussian and Gaussian non-identity-link fits.
  The control is stored on the fit and reused by
  [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md),
  [`update()`](https://rdrr.io/r/stats/update.html) and the internal
  refits (bootstrap, profiles, predictions).

- ...:

  Reserved for future use.

## Value

An object of class `mm_glmm`, also inheriting from `mm_fit` and
`mm_compiled`.

## Details

Optimization runs inside a single native call with no progress output:
the pre-fit explanation block (when `verbose >= 0`) is the last thing
printed before the fitted result returns. The fit checks for a user
interrupt (Ctrl-C / Esc) between optimizer evaluations and stops with a
typed `mm_interrupted` error (the check happens at evaluation
boundaries, so a single very expensive evaluation finishes first).
Evaluation budgets are bounded (a bounded budget caps optimizer
iterations; it does not prove every native evaluation terminates);
runtime on large problems is governed by `mm_control(max_feval = )`.

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rbinom(120, 1, 0.5),
  x = rnorm(120),
  g = factor(rep(seq_len(12), each = 10))
)
# Default: glmer-equivalent joint Laplace estimates.
fit <- glmm(y ~ x + (1 | g), df, family = binomial(),
            control = mm_control(verbose = -1))
fixef(fit)
#> (Intercept)           x 
#> -0.10067766  0.01797541 
# Fast profiled PIRLS path (lme4's nAGQ = 0):
fit_fast <- glmm(y ~ x + (1 | g), df, family = binomial(),
                 method = "pirls_profiled",
                 control = mm_control(verbose = -1))
fixef(fit_fast)
#> (Intercept)           x 
#> -0.10067766  0.01797541 
```
