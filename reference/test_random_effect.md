# Test a random-effect variance component

`test_random_effect()` exposes the boundary-aware likelihood-ratio route
for random-effect variance components. The v1 certified route is a
nested ML comparison that adds exactly one variance/covariance parameter
and reports the Self-Liang 50:50 mixture reference distribution. It is
intentionally separate from
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md),
which tests fixed effects.

## Usage

``` r
test_random_effect(
  fit,
  term,
  method = c("boundary_lrt"),
  refit_for_comparison = c("auto", "error", "ml"),
  ...
)

# S3 method for class 'mm_lmm'
test_random_effect(
  fit,
  term,
  method = c("boundary_lrt"),
  refit_for_comparison = c("auto", "error", "ml"),
  ...
)
```

## Arguments

- fit:

  A fitted `mm_lmm`.

- term:

  Random-effect term to test. This can be the term id (`"r0"`), the
  original random-effect fragment such as `"(1 | subject)"`, or a unique
  grouping factor name such as `"subject"`.

- method:

  Currently `"boundary_lrt"`.

- refit_for_comparison:

  How to handle REML fits. `"auto"` and `"ml"` refit to ML; `"error"`
  refuses.

- ...:

  Reserved for future methods.

## Value

An `mm_random_effect_test` object with a one-row `table`.
