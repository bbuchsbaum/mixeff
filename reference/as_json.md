# Serialize a mixeff spec or fit to JSON

`as_json()` returns a JSON string that carries the parsed object's
public R-side state and the raw compiler artifact JSON.
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) remains the primary
persistence path for fitted objects;
[`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
restores process-local caches after deserialization.

## Usage

``` r
as_json(x, pretty = FALSE, ...)

# S3 method for class 'mm_compiled'
as_json(x, pretty = FALSE, ...)
```

## Arguments

- x:

  A compiled `mm_spec` or fitted `mm_fit`.

- pretty:

  Logical; pretty-print JSON when `TRUE`.

- ...:

  Reserved for future methods.

## Value

A length-one character string containing JSON.
