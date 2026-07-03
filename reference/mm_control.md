# Control mixeff fitting behavior

`mm_control()` collects small R-side controls for
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md).
`verbose = -1` suppresses the pre-fit
[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
message; non-negative values emit it once before optimization (it
travels on the message stream, so
[`suppressMessages()`](https://rdrr.io/r/base/message.html) and knitr's
`message = FALSE` also quiet it).

## Usage

``` r
mm_control(
  verbose = 0L,
  max_feval = NULL,
  optimizer = NULL,
  start = NULL,
  ftol_rel = NULL,
  ftol_abs = NULL,
  xtol_rel = NULL
)
```

## Arguments

- verbose:

  Integer verbosity level. Use `-1` to suppress the automatic model
  explanation (and the GLMM estimator notice).

- max_feval:

  Optional positive integer capping the optimizer's objective
  evaluations. Most useful for
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) with
  `method = "joint_laplace"`, whose native joint optimizer otherwise
  runs to an engine-chosen budget. `NULL` (default) leaves the engine
  default in place.

- optimizer:

  Optional optimizer name, overriding the driver's automatic choice. One
  of `"auto"` (default behaviour), `"bobyqa"`, `"newuoa"`, `"cobyla"`,
  `"pattern_search"`, `"trust_bq"`, or the PRIMA variants
  (`"prima_bobyqa"`, `"prima_cobyla"`, `"prima_lincoa"`,
  `"prima_newuoa"`). An unsupported or not-compiled choice raises a
  typed error rather than silently falling back. `NULL`/`"auto"` keep
  automatic selection.

- start:

  Optional numeric warm-start vector for the covariance parameters
  (theta). Its length must match the model's theta dimension (the engine
  validates this). `NULL` (default) cold-starts.

- ftol_rel, ftol_abs:

  Optional positive relative/absolute convergence tolerances on the
  objective. `NULL` keeps the engine default.

- xtol_rel:

  Optional positive relative convergence tolerance on the optimizer
  parameters. `NULL` keeps the engine default.

## Value

A list of class `mm_control`.

## Details

By default the fit driver selects the optimizer and its tolerances
automatically (see
[`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md)
to inspect what ran). The `optimizer`, `start`, and `ftol_*`/`xtol_rel`
arguments are a narrow, opt-in escape hatch — for recourse when the
default fails to converge, for warm starts, and for explicit tolerance
overrides. Any override you supply is recorded in the optimizer
certificate, so the fit stays auditable.

## See also

[`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md)
to inspect which optimizer ran and whether a caller override was
applied.
