# Bootstrap Fixed-Effect Payload Surface

Status: design plan; no R-side bootstrap fixed-effect p-value computation.

## Boundary

`mixeff` must treat fixed-effect bootstrap p-values as Rust-owned inference
rows. R may request, format, validate, and preserve metadata, but it must not
derive a fixed-effect p-value from full-model bootstrap distributions or from
R-side refit summaries.

The current Rust contract distinguishes bootstrap distributions from bootstrap
hypothesis tests. A printable p-value requires a certified `fixed_effect_null`
payload and a returned `mixedmodels.fixed_effect_inference_table` row with
`method = "bootstrap"` and `status = "available"`.

Until the Rust bridge exposes that callable payload path, the current behavior
is correct:

```r
contrast(fit, L, method = "bootstrap")
```

returns the Rust row with `status = "not_assessed"`, `p_value = NA`, and the
Rust reason explaining that a certified bootstrap payload is required.

## Recommended Public API

Use `contrast(..., method = "bootstrap", bootstrap = bootstrap_control())` as
the primary surface. Do not introduce a separate `bootstrap_contrast()` front
door in the first pass.

Reasons:

- `contrast()` is already the primitive inferential operation.
- Bootstrap is a method choice, not a different estimand.
- The result can reuse the same `mm_contrast` table and the same
  `method/status/reliability/reason/notes` columns.
- A control object keeps simulation settings explicit without overloading
  `contrast()` with many scalar arguments.

Proposed user call:

```r
contrast(
  fit,
  L,
  rhs = 0,
  method = "bootstrap",
  bootstrap = bootstrap_control(nsim = 999, seed = 1)
)
```

`bootstrap_control()` should initially accept:

- `nsim`: requested replicate count.
- `seed`: reproducibility seed; `NULL` means Rust records generated seed state.
- `failed_refit`: one of `"error"`, `"discard"`, `"record_unavailable"`.
- `min_success`: minimum successful refits required before a p-value can be
  reported.
- `parallel`: reserved configuration passed through only after Rust certifies
  reproducible parallel execution.

`contrast(..., method = "bootstrap")` without `bootstrap =` should keep the
current Rust-owned unavailable row rather than silently choosing defaults.

## Bridge Contract

Add a Rust bridge only after upstream exposes a single call that accepts the
model, fixed-effect hypothesis, and bootstrap run control and returns a
versioned inference table.

Candidate bridge:

```text
mm_fixed_effect_bootstrap_contrast_json(
  formula,
  reml,
  data_payload,
  control_json,
  contrast_payload,
  bootstrap_control_json
) -> mixedmodels.fixed_effect_inference_table
```

R responsibilities:

- validate `L` dimensions, row labels, `rhs` length, and control object shape;
- serialize the fit/data payload already used by `mm_fixed_effect_contrast_json`;
- pass bootstrap settings through unchanged;
- parse Rust's returned inference table;
- preserve run metadata in `mm_contrast$raw` or a dedicated
  `mm_contrast$bootstrap` field.

Rust responsibilities:

- construct the fixed-effect null target;
- simulate under that null target;
- refit and account for failed refits;
- compute the continuity-corrected p-value;
- compute/report Monte Carlo standard error;
- decide row status and reliability;
- emit all unavailable reasons and reproducibility metadata.

## Result Shape

Keep `mm_contrast$table` stable. Bootstrap-specific details belong in list
columns or raw metadata, not in extra top-level scalar columns that only some
methods use.

Expected available row fields:

- `estimate`: tested difference, `L beta_hat - rhs`.
- `std_error`: Rust-provided row standard error if part of the bootstrap row.
- `df`: `NA` unless Rust supplies a meaningful denominator df.
- `statistic`: Rust-owned test statistic.
- `statistic_name`: likely `"t"` or `"F"` depending on scalar/multi-df row.
- `p_value`: continuity-corrected bootstrap p-value.
- `method`: `"bootstrap"`.
- `status`: `"available"`.
- `reliability`: Rust grade.
- `reason`: `NA`.
- `notes`: includes MCSE and replicate-accounting notes.

The object should also expose run metadata:

```r
attr(ct$table, "bootstrap")
# or
ct$bootstrap
```

with at least:

- requested and successful replicate counts;
- failed-refit count and policy;
- seed/reproducibility state;
- MCSE;
- payload schema name/version;
- null-target summary.

## Tests And Fixtures

Before enabling available bootstrap p-values in R, add fixtures or bridge smoke
tests that pin:

- explicit `method = "bootstrap"` with no payload returns Rust unavailable
  reason, as it does now;
- available scalar contrast row from a Rust-certified `fixed_effect_null`
  payload;
- nonzero `rhs` uses `L beta_hat - rhs` consistently;
- failed-refit handling records counts and either withholds or reports p-values
  according to Rust status;
- too few successful replicates yields missing p-value with Rust reason;
- seed/reproducibility metadata survives `saveRDS()` / `readRDS()`;
- `summary()`, `contrast()`, and future `test_effect()` never reconstruct
  bootstrap p-values from full-model bootstrap output;
- existing `parametric_bootstrap()` LRT output remains separate from
  fixed-effect bootstrap contrast output.

## Non-Goals For First Pass

- No R implementation of fixed-effect bootstrap p-value formulas.
- No automatic bootstrap selection from `method = "auto"`.
- No reuse of model-comparison `parametric_bootstrap()` results for fixed
  effects.
- No parallel bootstrap API until Rust records reproducible execution state.
- No `bootstrap_contrast()` alias unless user testing shows the control-object
  form is too hard to discover.

## Open Questions For Upstream

- What is the final schema name/version for the bootstrap run metadata payload?
- Does the bootstrap row carry MCSE only in `notes`, or should it become a
  structured field in the fixed-effect inference schema?
- Should `min_success` be an R control field or a Rust policy default recorded
  in metadata?
- Should multi-df bootstrap tests land in the same first bridge, or should the
  first R surface only enable scalar contrasts?
