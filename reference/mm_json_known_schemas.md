# Closed list of schema/version pairs the wrapper understands

Returns the known-schema table that backs
[`mm_json_negotiate()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_negotiate.md).
New schemas appear here as later phases add artifacts (compiled-model,
audit, theta_map, certificate, inference, reproducibility, prediction).

## Usage

``` r
mm_json_known_schemas()
```

## Value

A data frame with two character columns: `name` and `version`.

## See also

[`mm_json_negotiate()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_negotiate.md).

## Examples

``` r
mm_json_known_schemas()
#>                                        name version
#> 1                                   formula      v0
#> 2       mixedmodels.compiled_model_artifact       1
#> 3            mixedmodels.model_audit_report       2
#> 4              mixedmodels.random_term_card       1
#> 5  mixedmodels.fixed_effect_inference_table   1.0.0
#> 6       mixedmodels.marginal_quantity_table   1.0.0
#> 7        mixedmodels.model_comparison_table   1.0.0
#> 8                  mixedmodels.boundary_lrt   1.0.0
#> 9                   mixedmodels.fit_summary   1.0.0
#> 10        mixedmodels.profile_likelihood_ci   1.0.0
```
