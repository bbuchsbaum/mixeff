# Compare fitted mixeff models

`compare()` is the namespace-qualified model-comparison front door. For
LMMs it reports likelihood, information criteria, and asymptotic
likelihood-ratio comparisons. REML fits are refit by ML when
`refit_for_comparison = "auto"` or `"ml"`; `"error"` refuses that
comparison.

## Usage

``` r
compare(object, ...)

# S3 method for class 'mm_lmm'
compare(
  object,
  ...,
  target = c("fixed_effects", "random_effects", "prediction"),
  method = c("auto", "lrt", "bootstrap", "aic", "kenward_roger", "satterthwaite"),
  refit_for_comparison = c("auto", "error", "ml"),
  nsim = 0L,
  seed = NULL,
  threads = 1L
)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- ...:

  Additional fitted `mm_lmm` objects.

- target:

  Comparison target label.

- method:

  `"auto"` / `"lrt"` for asymptotic likelihood-ratio rows, `"aic"` for
  information criteria only, `"bootstrap"` for a small
  parametric-bootstrap LRT when `nsim > 0`, or `"kenward_roger"` /
  `"satterthwaite"` for pbkrtest's `KRmodcomp()` / `SATmodcomp()`: an F
  test of the larger model's fixed effects restricted to the smaller
  model's (two models with identical random effects, nested in their
  fixed effects). The restriction matrix is built from the two design
  matrices as pbkrtest does; the test runs on the REML fit of the larger
  model (refitted by REML if it was fitted by ML, as pbkrtest does). The
  Kenward-Roger F is pbkrtest's scaled statistic (the main `Ftest` row
  of `KRmodcomp()`, which lmerTest's `anova(ddf = "Kenward-Roger")` also
  reports): \\F = \lambda F_U\\ referred to \\F(q, \nu)\\ with the KR
  denominator df, where \\\lambda\\ is the KR scaling factor
  (`f_scaling`; 1 for a single restriction). The unscaled statistic and
  its p-value (pbkrtest's `FtestU` row) are kept as `unscaled_statistic`
  / `unscaled_p_value`. The F row replaces the LRT row
  (`statistic_name = "F"`, `df` = numerator df, `den_df` = denominator
  df); the full result is in `$fixed_f`.

- refit_for_comparison:

  How to handle REML fits.

- nsim:

  Number of bootstrap simulations for `method = "bootstrap"`.

- seed:

  Optional bootstrap seed.

- threads:

  Worker threads for the bootstrap refits (default 1; see
  [`bootstrap_control()`](https://bbuchsbaum.github.io/mixeff/reference/bootstrap_control.md)).
  The result is identical for every value.

## Value

An `mm_model_comparison` object with a data-frame `table`.
