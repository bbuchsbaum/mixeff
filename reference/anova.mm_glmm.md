# Analysis of deviance for GLMMs

With two or more fitted models,
[`anova()`](https://rdrr.io/r/stats/anova.html) performs a sequential
likelihood-ratio comparison (like `anova(glmer1, glmer2)`). For a single
model, fixed-effect tests are routed to
[`drop1()`](https://rdrr.io/r/stats/add1.html) (term LRTs),
[`summary()`](https://rdrr.io/r/base/summary.html) (Wald z), or
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
(custom Wald contrasts), which the GLMM contract supports directly.

## Usage

``` r
# S3 method for class 'mm_glmm'
anova(object, ...)
```

## Arguments

- object:

  A fitted `mm_glmm`.

- ...:

  Additional fitted models to compare.

## Value

An `mm_model_comparison` object (multi-model case).
