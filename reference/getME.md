# Extract low-level model components

`getME()` provides a small, honest subset of the familiar lme4
extractor. The fixed-effect design (`"X"`), random-effect design
(`"Z"`), relative covariance factor (`"Lambda"` / `"Lambdat"`), grouping
factors (`"flist"`), random coefficient names (`"cnms"`), response
(`"y"`), fixed coefficients (`"beta"` / `"fixef"`), and theta vector
(`"theta"`) are rebuilt lazily from the serialized R object.

## Usage

``` r
getME(object, name, ...)

# S3 method for class 'mm_lmm'
getME(object, name, ...)
```

## Arguments

- object:

  A fitted `mm_lmm` object.

- name:

  Component name, or a character vector of names.

- ...:

  Reserved for future methods.

## Value

The requested component, or a named list for multiple names.
