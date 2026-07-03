# Tidy, glance, and augment methods for mixeff fits

These implement the broom / broom.mixed generics for
[`mm_lmm`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
[`mm_glmm`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits,
so `tidy()`, `glance()`, and `augment()` work on mixeff models the same
way they do on `lme4` fits.

## Arguments

- x, data:

  A fitted `mm_lmm` or `mm_glmm` (and, for `augment()`, optional data to
  augment; defaults to the model frame).

- effects:

  Which terms to return: any of `"fixed"`, `"ran_pars"`, `"ran_vals"`.

- conf.int:

  Logical; add Wald `conf.low`/`conf.high` for fixed effects.

- conf.level:

  Confidence level for `conf.int`.

- ...:

  Unused; for generic compatibility.

## Value

A data frame.

## Details

`tidy()` returns one row per model term. `effects = "fixed"` yields the
fixed-effect coefficients (`estimate`, `std.error`, `statistic`, and,
for GLMMs, a Wald `p.value`); `effects = "ran_pars"` yields the
variance- component standard deviations (`sd__<term>`), correlations
(`cor__<a>.<b>`), and the residual SD (`sd__Observation`);
`effects = "ran_vals"` yields the conditional modes. `glance()` returns
a one-row model-summary frame; `augment()` returns the model frame with
`.fitted` and `.resid` columns.

These methods are registered with generics when the package is loaded;
call them via
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) /
[`broom.mixed::tidy()`](https://generics.r-lib.org/reference/tidy.html)
etc.
