# Inspect the optimizer certificate

Inspect the optimizer certificate

## Usage

``` r
optimizer_certificate(fit, ...)

# S3 method for class 'mm_compiled'
optimizer_certificate(fit, ...)
```

## Arguments

- fit:

  A compiled `mm_spec` or fitted `mm_fit`.

- ...:

  Reserved for future methods.

## Value

An `mm_optimizer_certificate` object containing the raw certificate and
a compact table view.
