# The wrapper's formula manifest

A versioned record of what `mixeff` currently supports — formula syntax
surface, schema versions per artifact type, and capability flags. The
manifest is the wrapper's machine-readable answer to *"what does this
build know how to do?"*. Fitted objects do not carry the whole manifest:
each one stamps the narrower `fit$schema` record (artifact schema name
and version, crate version, package version), which the serialization
and reporting paths read back. What lets
[`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md)
describe a fit after [`readRDS()`](https://rdrr.io/r/base/readRDS.html),
with no live Rust handle, is the durable compiled artifact stored beside
it.

## Usage

``` r
mm_formula_manifest()
```

## Value

A named list with the following elements:

- `mixeff_rust_version`:

  Version of the bundled extendr crate.

- `crate_version`:

  Version of the bundled `mixeff-rs` upstream crate.

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

Capability flags evolve across releases; gate behavior on flags rather
than on package version.

## Examples

``` r
m <- mm_formula_manifest()
m$schema_versions$formula
#> [1] "v0"
m$capabilities$parse_formula
#> [1] TRUE
```
