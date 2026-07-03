# Inspect mixeff diagnostics and fit status

`diagnostics()` returns the structured diagnostics carried by a compiled
spec or fitted model artifact. `fit_status()` is the compact status
string recorded by the optimizer certificate for fitted models.

## Usage

``` r
diagnostics(fit, severity = NULL, stage = NULL, ...)

# S3 method for class 'mm_compiled'
diagnostics(fit, severity = NULL, stage = NULL, ...)

fit_status(fit, ...)

# S3 method for class 'mm_fit'
fit_status(fit, ...)

# S3 method for class 'mm_compiled'
fit_status(fit, ...)
```

## Arguments

- fit:

  A compiled `mm_spec` or fitted `mm_fit`.

- severity:

  Optional character vector used to filter diagnostics by severity.

- stage:

  Optional character vector used to filter diagnostics by stage.

- ...:

  Reserved for future methods.

## Value

`diagnostics()` returns an `mm_diagnostics` object containing the raw
diagnostic list and a data-frame view. `fit_status()` returns a
length-one character string.
