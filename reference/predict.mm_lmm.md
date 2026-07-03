# Predict from a fitted mixeff LMM

Predictions follow the lme4 generic shape. In-sample predictions reuse
the cached fitted/fixed values; new-data predictions are dispatched
through the Rust `predict_new` contract.

## Usage

``` r
# S3 method for class 'mm_lmm'
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

# S3 method for class 'mm_lmm'
fitted(object, ...)

# S3 method for class 'mm_lmm'
residuals(
  object,
  type = c("response", "pearson", "deviance", "working"),
  scaled = FALSE,
  ...
)

# S3 method for class 'mm_glmm'
fitted(object, ...)

# S3 method for class 'mm_glmm'
residuals(object, type = c("response"), ...)
```

## Arguments

- object:

  A fitted `mm_lmm` object.

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

- scaled:

  Logical; when `TRUE`, residuals are divided by the residual scale.

## Value

A numeric vector, or a list with `fit` and `se.fit` when
`se.fit = TRUE`.
