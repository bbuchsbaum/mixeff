# Predict from a fitted mixeff GLMM

GLMM predictions are computed on the R side from the stored fixed
effects (population, `re.form = NA`) or fixed effects plus conditional
modes (`re.form = NULL`), then mapped through the family link. This
mirrors
[`lme4::predict.merMod`](https://rdrr.io/pkg/lme4/man/predict.merMod.html)
for generalized models: `type = "link"` returns the linear predictor and
`type = "response"` the mean.

## Usage

``` r
# S3 method for class 'mm_glmm'
predict(
  object,
  newdata = NULL,
  re.form = NULL,
  allow.new.levels = FALSE,
  type = c("response", "link"),
  se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A fitted `mm_glmm` object.

- newdata:

  Optional new data. Must be a `data.frame` containing every variable
  referenced by the model's formula. Categorical levels must either
  match the training factor levels or trigger the `allow.new.levels`
  policy.

- re.form:

  Random-effects conditioning, following lme4's basic convention. `NULL`
  returns conditional predictions; `NA` (or `~0`) returns
  population-level (fixed-effect) predictions. Conditioning on a subset
  of grouping factors via a one-sided formula is not supported by the
  current Rust contract and raises `mm_inference_unavailable`.

- allow.new.levels:

  When `FALSE` (default), unseen grouping levels in `newdata` raise
  `mm_inference_unavailable` through the Rust `NewReLevels::Error`
  policy. When `TRUE`, unseen levels are replaced by the population mean
  (zero random effect), matching
  `lme4::predict(allow.new.levels = TRUE)`.

- type:

  Prediction scale. Gaussian LMMs use the same values for `"response"`
  and `"link"`.

- se.fit:

  Logical; when `TRUE`, returns a list with `fit` and `se.fit`. For
  population predictions (`re.form = NA`) the standard error is the Wald
  SE of the fixed-effect linear predictor, `sqrt(diag(X V X'))`. For
  conditional predictions (`re.form = NULL`) the SE comes from the
  engine prediction-variance payload, which adds the random-effect
  contribution (BLUP variance and the fixed/random covariance). Rows the
  engine cannot certify — e.g. unseen grouping levels under
  `allow.new.levels = TRUE` — return `NA` with the engine's reason in
  the `mm_reason` attribute.
  ([`lme4::predict.merMod`](https://rdrr.io/pkg/lme4/man/predict.merMod.html)
  offers no conditional SE at all.)

- interval:

  Interval type: `"confidence"` for the fitted mean or `"prediction"`
  for a new observation (adds the residual variance). Population
  (`re.form = NA`) intervals are `fit +/- z*se` computed R-side;
  conditional (`re.form = NULL`) bounds come from the engine
  prediction-variance payload. Returns a matrix with `fit`/`lwr`/`upr`.

- level:

  Confidence level for `interval` / `se.fit` intervals.

- ...:

  Reserved for generic compatibility.

## Value

A numeric vector, or a list with `fit` and `se.fit` when
`se.fit = TRUE`.

## Details

In-sample response predictions reuse the engine's certified fitted
values.

Standard errors and confidence intervals: population (`re.form = NA`)
SEs are the fixed-effect Wald SE mapped through the link by the delta
method; conditional (`re.form = NULL`) SEs and confidence bounds come
from the engine prediction-variance payload. The engine certifies these
rows for `method = "joint_laplace"` fits and for default
`pirls_profiled` fits whose post-fit profiled-optimum certificate is
issued (per-row status `"available"`). Uncertified fits (e.g. singular
fits, or fits whose certificate fails) keep status `"degraded"`, and
their conditional SEs and bounds are withheld as `NA` with the engine's
reason in the `mm_reason` attribute — consistent with the package's "no
fake certainty" contract.

Prediction (future-observation) intervals (`interval = "prediction"`)
are available for conditional, response-scale predictions: the engine
returns quantiles of the plug-in predictive distribution (the family
conditional distribution mixed over link-scale fitted-mean uncertainty),
so bounds are integers for count families and support points for
Bernoulli. They are refused with a typed condition on the link scale
(future observations are response-scale objects), for population-level
requests, and for grouped binomial fits (the future trial count is not
representable in `newdata`).
