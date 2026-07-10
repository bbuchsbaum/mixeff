# Profile a fitted linear mixed model

Computes profile-likelihood intervals for the model's parameters via the
engine's certified profile payload and returns them as an `mm_profile`
object: `$table` has one row per profiled parameter (`parameter`,
`estimate`, `lower`, `upper`, `regularity`, `reason_code`). Under REML,
fixed-effect coefficients are not profiled (upstream contract); their
rows carry `reason_code = "profile_beta_unavailable_under_reml"` rather
than being silently dropped. Use
[`confint()`](https://rdrr.io/r/stats/confint.html) with
`method = "profile"` for the matrix form.

## Usage

``` r
# S3 method for class 'mm_lmm'
profile(fitted, which = NULL, level = 0.95, ...)
```

## Arguments

- fitted:

  A fitted `mm_lmm`.

- which:

  Optional character vector of parameter names to keep (coefficient
  names, `"sigma"`, `"theta1"`, ...).

- level:

  Confidence level for the reported interval endpoints.

- ...:

  Unused; for generic consistency.

## Value

An `mm_profile` object with `$table`, `$level`, `$fit_criterion`, and
`$notes`.

## See also

[`confint()`](https://rdrr.io/r/stats/confint.html) with
`method = "profile"` for the matrix form.
