# Analysis of variance for a mixeff LMM

With no extra models, [`anova()`](https://rdrr.io/r/stats/anova.html)
returns a term-level F-test table for `object` shaped like lmerTest's: a
data frame of class `"anova"` with columns `Sum Sq`, `Mean Sq`, `NumDF`,
`DenDF`, `F value` and `Pr(>F)` (rows named by term;
`Mean Sq = F * sigma^2` and `Sum Sq = NumDF * Mean Sq`, as in lmerTest).
With further fitted models in `...`, it returns lme4's likelihood-ratio
table (`npar`, `AIC`, `BIC`, `logLik`, `-2*log(L)`, `Chisq`, `Df`,
`Pr(>Chisq)`, rows named after the arguments; lme4 before 2.0 called the
`-2*log(L)` column `deviance`), computed by
[`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)
(REML fits are refit by ML, as in lme4).

## Usage

``` r
# S3 method for class 'mm_lmm'
anova(
  object,
  ...,
  type = c("III", "II", "I", "block"),
  method = c("auto", "satterthwaite", "kenward_roger", "bootstrap", "asymptotic", "none"),
  refit_for_comparison = c("auto", "error", "ml"),
  ddf = NULL
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
  lmerTest's spellings `3`, `2`, `1` (numeric or character) are
  accepted.

- method:

  Degrees-of-freedom / test method.

- refit_for_comparison:

  Passed to
  [`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)
  when `...` is used.

- ddf:

  lmerTest's spelling of the method: `"Satterthwaite"`
  (`method = "satterthwaite"`), `"Kenward-Roger"`
  (`method = "kenward_roger"`), or `"lme4"`, which returns lme4's own
  single-model table (sequential `npar`, `Sum Sq`, `Mean Sq`, `F value`,
  no denominator df or p-values). Supply `ddf` or `method`, not both.

## Value

A data frame of class `c("mm_anova", "anova", "data.frame")` for one
model, or
`c("mm_anova_comparison", "mm_model_comparison", "anova", "data.frame")`
for several.

## Details

mixeff's provenance stays on the result: `x$table` is the full
term-level (or
[`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md))
table with method, status, reliability and reason columns, and `x$type`
/ `x$requested_method` record the request. Rows the engine could not
certify keep `NA` statistics and their reason is printed under the
table.
