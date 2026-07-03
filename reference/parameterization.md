# Inspect covariance parameterization

`parameterization()` exposes the fitted theta/Lambda mapping recorded in
the compiler artifact. It is the R table view of the upstream theta-map
and covariance-parameter trace records.

## Usage

``` r
parameterization(fit, ...)

# S3 method for class 'mm_compiled'
parameterization(fit, ...)
```

## Arguments

- fit:

  A compiled `mm_spec` or fitted `mm_fit`.

- ...:

  Reserved for future methods.

## Value

An `mm_theta_map` object with a data-frame `table` and raw trace
records.
