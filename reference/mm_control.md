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
  xtol_rel = NULL,
  optCtrl = NULL,
  disp_method = NULL,
  disp_dof_correction = NULL,
  max_phi_iter = NULL
)
```

## Arguments

- verbose:

  Integer verbosity level. Use `-1` to suppress the automatic model
  explanation and the fit notices (the GLMM estimator notice, the
  `mm_rows_dropped` missing-value notice, grouping coercion and scaling
  advisories).

- max_feval:

  Optional positive integer capping the optimizer's objective
  evaluations. Most useful for
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) with
  `method = "joint_laplace"`, whose native joint optimizer otherwise
  runs to an engine-chosen budget. `NULL` (default) leaves the engine
  default in place.

- optimizer:

  Optional optimizer name, overriding the driver's automatic choice. One
  of `"auto"` (default behavior), `"bobyqa"`, `"newuoa"`, `"cobyla"`,
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

- optCtrl:

  Optional named list in the style of lme4's `lmerControl(optCtrl = )`,
  translated to the engine's controls: `maxfun`/`maxeval` -\>
  `max_feval`; `ftol_rel`/`ftol_abs`/`xtol_rel` -\> the arguments of the
  same name; `xtol_abs` (scalar or one per theta) -\> the absolute
  parameter tolerance; `rhobeg` (scalar or one per theta) -\> the
  optimizer's initial step. Unknown names are refused with an
  `mm_arg_error` rather than ignored. lme4's `check.conv.*` options have
  no counterpart: mixeff does not emit post-hoc convergence warnings;
  the fit's convergence status is the engine's typed certificate
  ([`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md),
  [`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md)).

- disp_method, disp_dof_correction, max_phi_iter:

  lme4 2.1-0's `glmerControl()` dispersion controls for
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits
  with a free dispersion parameter (Gamma, inverse Gaussian, Gaussian
  with a non-identity link). `disp_method = "moment"` (the default)
  profiles the dispersion `phi` in a damped fixed-point loop around
  PIRLS, with the working weights divided by `phi`; `"old/buggy"`
  restores lme4 \< 2.1 (working weights with `phi = 1`, theta relative
  to [`sigma()`](https://rdrr.io/r/stats/sigma.html)).
  `disp_dof_correction` (default `TRUE`) divides the deviance by
  `n - rank([X, Z])` rather than `n` in the moment estimator, and
  `max_phi_iter` (lme4's `maxPhiIter`, default 100) caps the `phi`
  iterations per objective evaluation. Both only matter for `"moment"`.
  `NULL` keeps the engine default (the lme4 2.1-0 default). Invalid
  values raise an `mm_arg_error`. For families without a free dispersion
  (binomial, Poisson, negative binomial) and for
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) the
  settings have no effect: they are kept on the control, and the fit
  announces that they were ignored with an `mm_control_ignored_notice`
  message (silenced by `verbose = -1`). lme4 ignores them silently.

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
