# The wrapper's formula manifest

A versioned record of what `mixeff` currently supports — formula syntax
surface, schema versions per artifact type, and capability flags. The
manifest is the wrapper's machine-readable answer to *"what does this
build know how to do?"*. Every `mm_fit` object created by future phases
will store a snapshot of `mm_formula_manifest()` at construction time so
the wrapper's
[`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md) verb
(Phase 1+) can answer the same question for old fits without consulting
the Rust handle.

## Usage

``` r
mm_formula_manifest()
```

## Value

A named list with the following elements:

- `mixeff_rust_version`:

  Version of the bundled extendr crate.

- `crate_version`:

  Version of the bundled `mixedmodels` upstream crate.

- `schema_versions`:

  Named list, one entry per artifact schema the wrapper currently emits
  or consumes.

- `formula_features`:

  Named list with `operators`, `intercept_forms`, `random_term_forms`,
  and `transformations` — the lme4-style syntax surface.

- `capabilities`:

  Named list of logical flags (`parse_formula`, `compile_model`,
  `audit_design`, `explain_model`, `random_options`,
  `compare_covariance`, `fit_lmm`, `fit_glmm`, `audit`, `changes`,
  `diagnostics`, `fit_status`, `parameterization`, `roles`, `as_json`,
  `simulate`, `inference`, `model_comparison_table`,
  `fit_summary_payload`, `marginal_quantity_table`,
  `marginal_quantities`, `verify_convergence`). The
  `marginal_quantity_table` schema may be available before the
  corresponding `marginal_quantities` verbs are implemented.

## Details

Capability flags evolve over phases; gate behavior on flags rather than
on package version.

## Examples

``` r
m <- mm_formula_manifest()
m$schema_versions$formula
#> [1] "v0"
m$capabilities$parse_formula
#> [1] TRUE
```
