# Audit a fitted mixeff model

`audit()` is the post-fit alias for
[`audit_design()`](https://bbuchsbaum.github.io/mixeff/reference/audit_design.md).
It renders the same upstream-authored audit report, now backed by the
fitted artifact carried by an `mm_fit`.

## Usage

``` r
audit(fit, ...)

# S3 method for class 'mm_fit'
audit(fit, ...)
```

## Arguments

- fit:

  A fitted `mm_fit`, usually from
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md).

- ...:

  Reserved for future methods.

## Value

An `mm_audit` object; see
[`audit_design()`](https://bbuchsbaum.github.io/mixeff/reference/audit_design.md).
