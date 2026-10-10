# Profile a fitted linear mixed model

Computes profile-likelihood intervals for the model's parameters via the
engine's certified profile payload and returns them as an `mm_profile`
object: `$table` has one row per profiled parameter (`parameter`,
`parameter_kind`, `estimate`, `lower`, `upper`, `regularity`, `status`,
`reason_code`, `reason`, `profiled_criterion`). Rows use lme4's names:
`.sig01`, ... (kind `"sd"` or `"cor"`, numbered in lme4's term order),
`.sigma`, and the fixed effects, all profiled on the ML deviance as lme4
does (a REML fit is refitted by ML); the engine's relative-Cholesky
`theta1`, ... rows follow, profiled on the fit's own criterion. Use
[`confint()`](https://rdrr.io/r/stats/confint.html) with
`method = "profile"` for the matrix form.

## Usage

``` r
# S3 method for class 'mm_lmm'
profile(fitted, which = NULL, level = 0.95, threads = 1L, ...)
```

## Arguments

- fitted:

  A fitted `mm_lmm`.

- which:

  Optional character vector of parameter names to keep (coefficient
  names, `".sig01"`, `".sigma"`, `"theta1"`, ...).

- level:

  Confidence level for the reported interval endpoints.

- threads:

  Worker threads (default 1); see
  [`confint.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/confint.mm_lmm.md).

- ...:

  Unused; for generic consistency.

## Value

An `mm_profile` object with `$table`, `$level`, `$fit_criterion`, and
`$notes`.

## Details

A parameter whose profile is irregular keeps its row: `status` is
`"non_monotone"` (the signed-root deviance is not monotone),
`"not_bracketing"` (the profile does not reach the cutoff on one side),
or `"failed"` (a constrained refit failed); `reason_code` holds the
matching `regularity` code (`"non_monotone_profile"`,
`"profile_not_bracketing"`, `"profile_failed"`) and `reason` the
engine's explanation. Bounds the profile cannot determine are `NA`; a
bound inside the regular part of the profile is still reported. lme4
instead warns and falls back to linear interpolation (often giving
`[-1, 1]` for a correlation); mixeff reports `NA` plus the reason rather
than an interpolated bound. The other parameters' intervals are
unaffected.

## See also

[`confint()`](https://rdrr.io/r/stats/confint.html) with
`method = "profile"` for the matrix form.
