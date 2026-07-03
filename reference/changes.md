# Show requested, effective, and fitted model-state changes

`changes()` summarizes the transitions recorded in the compiler
artifact: requested formula to effective formula, design-time reductions
or covariance transitions, and fitted covariance rank/status from the
optimizer certificate pass.

## Usage

``` r
changes(fit, ...)

# S3 method for class 'mm_compiled'
changes(fit, ...)
```

## Arguments

- fit:

  A compiled `mm_spec` or fitted `mm_fit`.

- ...:

  Reserved for future methods.

## Value

An `mm_change_log` object with a data-frame `table` and the raw artifact
fragments used to build it.
