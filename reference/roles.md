# Declare or inspect design roles

`roles()` has two Phase 1.F uses. With named string arguments it
constructs a declared-role object, e.g.
`roles(subject = "sampled_unit")`. With a compiled spec or fit as its
only unnamed argument, it returns the observed role fallback inferred
from the formula and model frame.

## Usage

``` r
roles(...)
```

## Arguments

- ...:

  Either named role strings or one unnamed `mm_spec` / `mm_fit`.

## Value

An `mm_roles` object with a data-frame `table`.
