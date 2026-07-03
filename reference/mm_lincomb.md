# Wald inference on a linear combination of fixed effects

Convenience helper for the common case of testing \\H_0:\\ c^\top \beta
= 0\\ where `c` is a sparse, named weight vector. The estimate is
\\c^\top \hat\beta\\, the standard error is \\\sqrt{c^\top V c}\\ where
`V` is the model's fixed-effect covariance, the statistic is the Wald
ratio, and the interval is the symmetric Wald CI at `level`.

## Usage

``` r
mm_lincomb(fit, weights, level = 0.95, method = NULL, ...)

# Default S3 method
mm_lincomb(fit, weights, level = 0.95, method = NULL, ...)

# S3 method for class 'mm_glmm'
mm_lincomb(fit, weights, level = 0.95, method = NULL, ...)

# S3 method for class 'mm_lmm'
mm_lincomb(
  fit,
  weights,
  level = 0.95,
  method = c("auto", "satterthwaite", "kenward_roger", "asymptotic"),
  ...
)
```

## Arguments

- fit:

  A fitted `mm_lmm` or `mm_glmm`.

- weights:

  A named numeric vector (or named list / single-row data.frame
  coercible to one). Names must match `names(fixef(fit))` exactly.

- level:

  Confidence level for the Wald interval. Default 0.95.

- method:

  For `mm_lmm`, the degrees-of-freedom method passed to
  [`df_for_contrast()`](https://bbuchsbaum.github.io/mixeff/reference/df_for_contrast.md).
  Defaults to `"auto"` (Satterthwaite when available). For `mm_glmm`,
  only `"asymptotic"` is accepted.

- ...:

  Reserved for future methods.

## Value

A single-row data.frame with columns `estimate`, `std_error`,
`statistic`, `statistic_name` (`"t"` or `"z"`), `df`, `p_value`,
`lower`, `upper`, and `method`. The result carries an `"mm_status"`
attribute reflecting the underlying vcov reliability (`status`,
`method`, `reliability`, `reason`).

## Details

For `mm_glmm`, the statistic is the asymptotic Wald *z* (no df). For
`mm_lmm`, the default is Satterthwaite-approximated *t* via
[`df_for_contrast()`](https://bbuchsbaum.github.io/mixeff/reference/df_for_contrast.md);
pass `method = "asymptotic"` to force Wald *z*.

Weight names must be a subset of `names(fixef(fit))`. Coefficients not
named in `weights` contribute zero. Pass the long-form
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
front door if you need multiple contrasts or non-default rhs.

## See also

[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
for the long-form, Rust-routed contrast surface with full estimability /
reliability reporting.

## Examples

``` r
if (FALSE) { # \dontrun{
# Difference-in-differences contrast at a focal SOA = 25 ms
# (Loo et al. 2026 aphantasia primary estimand, glmm path)
soa_s_25 <- (log(0.025) - mean(fit$data$soa_log)) / sd(fit$data$soa_log)
mm_lincomb(fit, c(
  "group: aphant:mask: masked"        = 1,
  "group: aphant:mask: masked:soa_s"  = soa_s_25
))
} # }
```
