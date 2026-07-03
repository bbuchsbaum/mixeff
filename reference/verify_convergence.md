# Verify convergence of a fitted linear mixed model

`verify_convergence()` re-runs the fit under the engine's bounded
verification workflow and reports whether the extra runs agree with the
fitted optimum: a restart from the optimum, one or more jittered
restarts, and (opt-in) an alternate-optimizer consensus pass. Agreement
is judged by the engine against the objective/theta/beta tolerances
below; the verdict (`status`), the per-run deltas, and the wording are
all owned by the Rust contract — R only formats them.

## Usage

``` r
verify_convergence(fit, ...)

# Default S3 method
verify_convergence(fit, ...)

# S3 method for class 'mm_lmm'
verify_convergence(
  fit,
  ...,
  restart = TRUE,
  jitter_starts = 1L,
  jitter_scale = 1e-04,
  consensus = FALSE,
  max_feval = 500L,
  objective_tolerance = 1e-05,
  theta_tolerance = 0.001,
  beta_tolerance = 1e-04
)
```

## Arguments

- fit:

  A fitted `mm_lmm` from
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md).

- ...:

  Reserved for future methods.

- restart:

  Logical; re-optimize starting from the fitted optimum and compare
  against it.

- jitter_starts:

  Number of restarts from jittered copies of the fitted covariance
  parameters.

- jitter_scale:

  Relative scale of the jitter applied to theta.

- consensus:

  Logical; also refit with an engine-chosen alternate optimizer and
  compare. Default `FALSE`: this vendored build compiles without the
  optional `nlopt` backend, and for some models the engine's alternate
  choice is an nlopt optimizer — its absence would then be reported as a
  non-agreeing run (status `fragile`) that reflects the build, not the
  fit. Enable it when you want the consensus pass and will read the
  per-run diagnostics.

- max_feval:

  Positive integer cap on objective evaluations per verification run.

- objective_tolerance, theta_tolerance, beta_tolerance:

  Positive agreement tolerances on the objective value, the covariance
  parameters, and the fixed effects.

## Value

An object of class `mm_convergence_verification` carrying:

- `status`:

  the engine verdict: `not_run`, `restart_agrees`,
  `optimizer_consensus`, `fragile`, or `unstable`

- `message`:

  the engine's one-line summary

- `table`:

  a data frame with one row per verification run (label, optimizer,
  return code, objective/theta/beta deltas, agreement)

- `reference`:

  the reference optimum the runs were compared to

- `tolerances`:

  the agreement tolerances that were applied

- `raw`:

  the parsed engine payload

## Details

The verifier refits the model from the stored specification before it
starts, so a call costs roughly `2 + jitter_starts` fits (plus consensus
runs when enabled).

## See also

[`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md)
for what the original fit ran;
[`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
to refit with a different optimizer or tolerances.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lmm(y ~ t + (1 | s), df)
verify_convergence(fit)
} # }
```
