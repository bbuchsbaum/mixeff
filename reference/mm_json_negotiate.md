# Negotiate a JSON schema header against what `mixeff` supports

Every artifact crossing the Rust-R bridge carries a header that names a
schema (e.g., `formula`, `artifact`, `audit`) and the version of that
schema (`v0`, `v1`, ...). `mm_json_negotiate()` validates a header
against the closed set of `(schema_name, schema_version)` pairs the
current wrapper build understands, and raises a typed `mm_schema_error`
on mismatch.

## Usage

``` r
mm_json_negotiate(header)
```

## Arguments

- header:

  A list with at least `schema_name` and `schema_version` as length-1
  character strings. Additional fields (e.g., `crate_version`,
  `package_version`) are accepted and ignored at this layer; downstream
  code records them on the `mm_fit` provenance.

## Value

Invisibly returns `TRUE` on a clean match.

## Details

This is the *fast-fail* primitive: any code path that consumes a Rust
artifact should call `mm_json_negotiate()` before parsing the body, so a
version skew between the Rust crate and the R wrapper produces a single
clean error rather than a confusing field-by-field decode failure.

## Errors

Any of the following raise an `mm_schema_error` (also inheriting from
`mm_condition` and `error`):

- `header` is not a list, or is missing `schema_name` /
  `schema_version`, or those fields are not length-1 character.

- `schema_name` is not in the wrapper's known set (see
  [`mm_json_known_schemas()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_known_schemas.md)).

- `schema_name` is known but `schema_version` does not match what the
  wrapper expects.

The condition object carries the original `header` in its `input` field.
(The field is *not* called `header` because rlang reserves that name on
conditions for `cnd_header()` formatting.)

## See also

[`mm_json_known_schemas()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_known_schemas.md)
for the closed set,
[`mm_formula_manifest()`](https://bbuchsbaum.github.io/mixeff/reference/mm_formula_manifest.md)
for the broader capability record.

## Examples

``` r
mm_json_negotiate(list(schema_name = "formula", schema_version = "v0"))
if (FALSE) { # \dontrun{
# Raises mm_schema_error:
mm_json_negotiate(list(schema_name = "formula", schema_version = "v99"))
} # }
```
