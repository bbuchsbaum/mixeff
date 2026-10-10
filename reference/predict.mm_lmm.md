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
  na.action = stats::na.pass,
  offset = NULL,
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
residuals(
  object,
  type = c("deviance", "pearson", "working", "response"),
  scaled = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted `mm_lmm` object.

- newdata:

  Optional new data. Must be a `data.frame` containing every variable
  the requested prediction needs (the fixed-effect variables, plus the
  grouping and slope variables of the random terms used); other columns
  are ignored. Grouping columns may be factor, character, or integer
  (coerced the way they were at fit time); fixed-effect factor and
  character columns are matched to the training levels, and a level not
  seen when fitting is refused. Stateful terms
  ([`poly()`](https://rdrr.io/r/stats/poly.html),
  [`scale()`](https://rdrr.io/r/base/scale.html), `ns()`, ...) are
  evaluated with the training basis, as in
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

- re.form:

  Random-effects conditioning, following lme4. `NULL` returns
  conditional predictions; `NA` (or `~0`) returns population-level
  (fixed-effect) predictions; a one-sided formula naming some of the
  model's random terms (e.g. `~ (1 | subject)`) conditions on those
  terms' conditional modes only. Standard errors and intervals are not
  available for such partial conditioning.

- allow.new.levels:

  When `FALSE` (default), unseen grouping levels in `newdata` raise
  `mm_inference_unavailable` through the Rust `NewReLevels::Error`
  policy. When `TRUE`, unseen levels are replaced by the population mean
  (zero random effect), matching
  `lme4::predict(allow.new.levels = TRUE)`.

- type:

  Prediction scale. Gaussian LMMs use the same values for `"response"`
  and `"link"`. For
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), the residual
  type, as in lme4: LMMs default to `"response"` (`y - mu`, also
  `"working"`), with `"pearson"`/`"deviance"` the weighted residuals
  `sqrt(w) * (y - mu)`; GLMMs default to `"deviance"` and also offer
  `"pearson"` (`(y - mu) * sqrt(w / V(mu))`), `"working"` and
  `"response"`.

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

- na.action:

  Missing-value handling for `newdata`, as in
  [`lme4::predict.merMod()`](https://rdrr.io/pkg/lme4/man/predict.merMod.html):
  the default [stats::na.pass](https://rdrr.io/r/stats/na.fail.html)
  returns `NA` for rows with a missing value in a variable the
  prediction needs;
  [stats::na.omit](https://rdrr.io/r/stats/na.fail.html) drops them and
  [stats::na.exclude](https://rdrr.io/r/stats/na.fail.html) pads them
  back as `NA`. In-sample predictions of a fit made with
  `na.action = na.exclude` are padded to the original rows, like
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html).

- offset:

  Offset for `newdata` rows when the model was fitted with the
  `offset =` argument: a numeric vector (one value per `newdata` row) or
  an expression evaluated in `newdata`. Offsets written in the formula
  (`offset(log(t))`) are evaluated from `newdata` automatically. (lme4
  silently drops an `offset =` argument offset for new data; mixeff
  refuses instead of guessing.)

- ...:

  Reserved for generic compatibility.

- scaled:

  Logical; when `TRUE`, residuals are divided by `sigma(object)` (once),
  as in lme4.

## Value

A numeric vector, or a list with `fit` and `se.fit` when
`se.fit = TRUE`.
