# Inspect random-effect blocks

`random_blocks()` summarizes the random-effect block structure recorded
in the compiler artifact: grouping factor, basis, covariance family,
theta parameter count, level counts, and design-support status.

## Usage

``` r
random_blocks(fit, ...)

# S3 method for class 'mm_compiled'
random_blocks(fit, ...)
```

## Arguments

- fit:

  A compiled `mm_spec` or fitted `mm_fit`.

- ...:

  Reserved for future methods.

## Value

An `mm_random_blocks` object with a data-frame `table`.
