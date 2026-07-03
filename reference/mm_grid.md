# Marginal grids, predictions, means, and comparisons

These helpers provide a small native marginal-quantities surface for
Gaussian LMM fits. They cover the common population-level workflow:
construct a reference grid, evaluate fixed-effect predictions, average
them into marginal means, and compare those means by simple differences.

## Usage

``` r
mm_grid(fit, specs, by = NULL, at = list(), cov.reduce = mean, ...)

# S3 method for class 'mm_lmm'
mm_grid(fit, specs, by = NULL, at = list(), cov.reduce = mean, ...)

mm_predictions(
  fit,
  grid = NULL,
  specs = NULL,
  by = NULL,
  at = list(),
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  level = 0.95,
  target = c("population"),
  scale = c("response", "link"),
  ...
)

# S3 method for class 'mm_lmm'
mm_predictions(
  fit,
  grid = NULL,
  specs = NULL,
  by = NULL,
  at = list(),
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  level = 0.95,
  target = c("population"),
  scale = c("response", "link"),
  ...
)

mm_means(
  fit,
  specs,
  by = NULL,
  at = list(),
  grid = NULL,
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  level = 0.95,
  weights = c("equal", "proportional"),
  target = c("population"),
  scale = c("response", "link"),
  ...
)

# S3 method for class 'mm_lmm'
mm_means(
  fit,
  specs,
  by = NULL,
  at = list(),
  grid = NULL,
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  level = 0.95,
  weights = c("equal", "proportional"),
  target = c("population"),
  scale = c("response", "link"),
  ...
)

mm_comparisons(
  fit,
  specs,
  by = NULL,
  at = list(),
  grid = NULL,
  comparison = c("difference", "ratio", "odds_ratio"),
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  level = 0.95,
  weights = c("equal", "proportional"),
  target = c("population"),
  scale = c("response", "link"),
  ...
)

# S3 method for class 'mm_lmm'
mm_comparisons(
  fit,
  specs,
  by = NULL,
  at = list(),
  grid = NULL,
  comparison = c("difference", "ratio", "odds_ratio"),
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  level = 0.95,
  weights = c("equal", "proportional"),
  target = c("population"),
  scale = c("response", "link"),
  ...
)
```

## Arguments

- fit:

  A fitted `mm_lmm`.

- specs:

  Character vector, or a one-sided formula such as `~ trt` or
  `~ trt | group`, naming the displayed marginal dimensions.

- by:

  Optional character vector of grouping variables for marginal summaries
  or pairwise comparisons.

- at:

  Named list of fixed-predictor values to force in the grid.

- cov.reduce:

  Function used to reduce numeric fixed predictors that are not
  explicitly gridded.

- ...:

  Reserved for future methods.

- grid:

  Optional object returned by `mm_grid()`.

- method:

  Requested inference method, passed to
  [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md).

- level:

  Confidence level for intervals computed from contrast standard errors.

- target:

  Prediction target. Only `"population"` is implemented.

- scale:

  Prediction scale. Gaussian LMMs have identical `"link"` and
  `"response"` scales.

- weights:

  Averaging weights for `mm_means()` and `mm_comparisons()`. `"equal"`
  weights reference-grid cells equally; `"proportional"` weights cells
  by observed fixed-factor frequencies.

- comparison:

  Comparison scale. Only `"difference"` is implemented.

## Value

`mm_grid()` returns an `mm_grid` object. The other helpers return an
`mm_marginal_quantity` object with a contract-shaped `table`.

## Details

The returned tables use the `mixedmodels.marginal_quantity_table` row
contract. Inference is routed through
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
so rows retain the same method, status, reliability, estimability, and
reason fields as fixed-effect contrasts. Ordinary full-rank LMMs use the
versioned `mixedmodels.fixed_effect_covariance_matrix` payload for
fixed-effect uncertainty; rank-deficient or otherwise uncertified fits
surface explicit unavailable status and reasons instead of partial
covariance numbers.
