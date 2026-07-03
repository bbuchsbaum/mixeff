# Drop one fixed-effect term at a time from a GLMM

Refits reduced fixed-effect GLMMs (random-effect terms preserved
exactly) and compares each to the full fit by asymptotic
likelihood-ratio test, mirroring `drop1(glmerMod, test = "Chisq")`.

## Usage

``` r
# S3 method for class 'mm_glmm'
drop1(object, scope = NULL, test = c("none", "Chisq"), ...)
```

## Arguments

- object:

  A fitted `mm_glmm`.

- scope:

  Optional character vector of fixed-effect terms to drop.

- test:

  `"Chisq"` reports asymptotic LRT rows; `"none"` reports information
  criteria only.

- ...:

  Reserved for future methods.

## Value

An `mm_drop1` object.
