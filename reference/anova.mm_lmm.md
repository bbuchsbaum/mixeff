# Analysis of variance for a mixeff LMM

With no extra models, [`anova()`](https://rdrr.io/r/stats/anova.html)
returns a term-level table for `object`. With further fitted models in
`...`, it defers to
[`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)
for a nested model comparison.

## Usage

``` r
# S3 method for class 'mm_lmm'
anova(
  object,
  ...,
  type = c("III", "II", "I", "block"),
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  refit_for_comparison = c("auto", "error", "ml")
)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- ...:

  Optional additional fitted models; triggers
  [`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md).

- type:

  Term-hypothesis type. `"III"` (default) tests marginal Type III
  hypotheses: a term's own contrast columns plus the equally-weighted
  average over the levels of every term that contains it, which makes
  the test invariant to which factor level is the reference and matches
  SAS / `car` / `lmerTest` Type III. `"II"` respects marginality (each
  term adjusted for all terms that do not contain it). `"I"` is
  sequential in [`terms()`](https://rdrr.io/r/stats/terms.html) order
  (main effects, then two-way interactions, and so on). `"block"` tests
  the raw coefficient block for each term — under treatment coding that
  is the simple effect at the other factors' reference levels, which is
  a legitimate quantity but is *not* Type III on unbalanced designs; it
  is the hypothesis mixeff computed for `"III"` before engine `1f3f689`.

- method:

  Degrees-of-freedom / test method.

- refit_for_comparison:

  Passed to
  [`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)
  when `...` is used.

## Value

An `mm_anova` object.
