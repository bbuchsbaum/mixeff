# Optional emmeans support for mixeff LMMs

These methods let `emmeans` build reference grids for `mm_lmm` objects
when the optional `emmeans` package is installed. They expose the same
fixed-effect design surface used by
[`mm_grid()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
and
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md).

## Usage

``` r
recover_data.mm_lmm(object, data = NULL, ...)

emm_basis.mm_lmm(
  object,
  trms,
  xlev,
  grid,
  method = c("auto", "satterthwaite", "kenward_roger", "asymptotic", "none"),
  ...
)

recover_data.mm_glmm(object, data = NULL, ...)

emm_basis.mm_glmm(object, trms, xlev, grid, ...)
```

## Arguments

- object:

  A fitted `mm_lmm`.

- data:

  Optional data override supplied by `emmeans`.

- trms, xlev, grid, ...:

  Arguments supplied by `emmeans`.

## Value

Objects expected by
[`emmeans::recover_data()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.html)
and
[`emmeans::emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.html).

## Details

The current bridge is intentionally narrow: Gaussian LMMs only and
population fixed-effect means only. When the fitted artifact carries an
available `mixedmodels.fixed_effect_covariance_matrix` payload,
`emmeans` receives that full fixed-effect covariance matrix. Native
[`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
[`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md),
and
[`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
remain the contract-preserving mixeff surface because they preserve
row-level status and reason fields.
