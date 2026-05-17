# Bootstrap uptake plan — wiring the new Rust engine into `mixeff`

Status: implementation plan; informs the R-side bridge, parser, and
API work needed to consume the upstream bootstrap and inference-vocabulary
changes.

## What just shipped upstream

The Rust engineer landed (per their summary message):

**Bootstrap surface**

- **Multi-df fixed-effect bootstrap.** Effective-rank Wald/F statistic
  from `L V L'`, generalized inverse when needed, numerator df =
  effective restriction rank. Scalar contrasts keep absolute studentized
  t. Symbol: same row machinery; multi-row `L` now accepted.
- **Bootstrap likelihood-ratio test engine.** Validates ML-only
  reduced/alt comparison, simulates from reduced model, refits both per
  replicate. Returns `BootstrapLikelihoodRatioTest` (observed statistic,
  p-value, MCSE, replicate stats, metadata, serde). REML refused.
- **Cluster bootstrap.** `DataFrame::cluster_resample` (row-level
  resampling with stable replicate-local relabeling for duplicate
  clusters); `LinearMixedModel::cluster_resample_full_model_contrast_payload`
  (rebuilds/refits and returns full-model interval payloads for scalar
  contrasts).
- **`BootstrapTargetKind` extended** with `LikelihoodRatio` and
  `ClusterResample`; `BootstrapLikelihoodRatioTest` and
  `ClusterResampleDraw` exported.
- **Replicate stats** are now in payloads (tests reconstruct p-values
  and MCSE from them) — request B4 from
  `upstream-bootstrap-coverage-issue.md` is satisfied.
- **Bench entrypoint** at `examples/bench_bootstrap.rs` emitting JSON
  rows for the three bootstrap modes.

**Diagnostic vocabulary surface** (per the engineer's summary of
`upstream-diagnostic-vocabulary-issue.md`):

- **R2 done** — typed inference reason codes; boundary-specific
  Satterthwaite unavailable case; stable `reason_code` plus headline
  reason; implementation detail retained separately.
- **R3 done** — closed `FixedEffectReliabilityReasonCode` enum on every
  inference row, with variants:
  `interior_converged_well_specified`,
  `asymptotic_wald_z_at_boundary`,
  `degrees_of_freedom_unavailable_so_z_substituted`,
  `satterthwaite_finite_difference_approximation`,
  `kenward_roger_approximation`,
  `bootstrap_monte_carlo_replicates`,
  `inference_unavailable_by_policy`,
  `contrast_not_estimable`,
  `standard_error_unavailable`.
- **R4 done** — duplicate `covariance_reduced` emissions removed; display
  aggregation dedupes.
- **R1 partially open** — bare `r0`, `theta[i]`, `varpar[i]` still leak
  through prose; tracked upstream as `bd-01KQF5NBJYXGNN09HFDNMVSPX2`.

## Ground state in mixeff

`src/rust/Cargo.toml` declares
`mixedmodels = { path = '/Users/bbuchsbaum/code/rust/mixedmodels' }` —
a path dependency, not a vendor. mixeff already builds against live
upstream. **No vendor refresh is required**; a clean rebuild
(`devtools::install(".")` or equivalent) is enough to pull the new
symbols into mixeff's compilation.

What is **not** yet in mixeff:

- Bridge entrypoints (`src/rust/src/lib.rs` + `extendr` wrappers) to
  call the new bootstrap LRT and cluster-resample functions.
- R-side JSON parsing for the new fields (`reason_code`,
  `reliability_reason_code`, `replicate_statistics`,
  multi-df numerator df, cluster-resample metadata).
- The user-facing verbs the API plan in `bootstrap-positioning.md`
  promised: `inference_options(fit, term)`, `test_effect(method = "bootstrap")`
  for multi-df, `test_effect(method = "bootstrap_lrt")`,
  `test_effect(method = "cluster_bootstrap")`, `confint(method = "bootstrap")`.

This plan converts the upstream changes into those user-facing surfaces.

## Phases

### Phase 0 — Verify

**P0.1.** Rebuild mixeff against the new upstream. Run `devtools::test()`
end-to-end. Expected: existing tests still pass — none of the upstream
changes are advertised as breaking. Any failure here is a compatibility
issue with current bridge code that needs investigation before further
work.

**P0.2.** Inspect a fresh `inference_table(fit)` JSON to confirm
`reason_code` and `reliability_reason_code` fields are present in the
wire payload. (The existing R parser silently ignores unknown fields,
so no breakage is expected; we just verify they arrived.)

**P0.3.** Inspect a fresh contrast bootstrap JSON to confirm
`replicate_statistics` is now populated. Earlier audit (saved in
`upstream-bootstrap-coverage-issue.md` §B4 evidence) showed the field
absent; the engineer's summary says tests now reconstruct from it, so
expect it present.

**P0.4.** Quick smoke test from R that the new Rust functions are
*reachable* via `extendr_api::call_r` or by writing a one-shot rextendr
test wrapper, before committing to bridge entrypoints.

Estimated effort: half-day. Output: a one-page "verification report"
appended to this plan.

### Phase 1 — Bridge entrypoints

Add to `src/rust/src/lib.rs` and the auto-generated
`R/extendr-wrappers.R`:

**P1.1. `mm_bootstrap_lrt_json`.** Inputs: reduced formula + alt
formula + the usual data-payload columns + `BootstrapOptions` JSON.
Output: serialised `BootstrapLikelihoodRatioTest`. Must reject REML
fits with the same typed reason the Rust engine uses.

**P1.2. `mm_cluster_resample_contrast_json`.** Inputs: formula +
columns + grouping-factor name + `L` + `rhs` + `BootstrapOptions`
JSON. Output: serialised cluster-resample contrast payload. The
`ClusterResampleDraw` per replicate is recorded in the run metadata
for reproducibility.

**P1.3. Verify `mm_fixed_effect_bootstrap_contrast_json` accepts
multi-row `L`.** The signature already takes `nrow(L)` and `ncol(L)`;
likely the only change is that the Rust side no longer rejects
`nrow(L) > 1`. Add a smoke test with a 2-df hypothesis.

**P1.4. Decide on `full_model_distribution` exposure for parametric
fixed-effect CIs.** The engineer's summary names
`cluster_resample_full_model_contrast_payload` (cluster-target
intervals) but does not explicitly mention a parametric
`full_model_distribution` bridge entrypoint. Three possibilities:
- The parametric target is exposed; just confirm.
- It is reachable but not yet bridged; add
  `mm_full_model_distribution_contrast_json`.
- It remains "Later scope" upstream and `confint(method = "bootstrap")`
  is supplied by the cluster path only in v0.

Resolve in P0 verification, plan branch accordingly.

### Phase 2 — R parser layer

Update `R/inference.R`, `R/json.R`, and any related parsing helpers:

**P2.1. Surface `reliability_reason_code` on every inference row.**
Add a column of the same name to `inference_table(fit)$table`. The
existing `reliability` column ("low"/"moderate"/"high") stays; the new
column carries the closed-enum machine-readable warrant. Print methods
show both.

**P2.2. Surface `reason_code` on unavailable / not-assessed rows.**
Currently the `reason` column carries the prose; add a parallel
`reason_code` column with the stable typed code. Brittle string-matching
in `mm_inference_unavailable_reason()` and friends is replaced by joins
on the typed code.

**P2.3. Parse `replicate_statistics` into a list-column or attribute on
the bootstrap row.** Decision: list-column on `details$bootstrap`, named
`replicate_statistics`, length = `successful_replicates`. Avoids
inflating the headline data frame; available for plotting / BCa /
custom MC-SE.

**P2.4. Parse multi-df bootstrap rows.** Numerator df = effective
restriction rank, statistic name = `f` (or whatever the engine emits
for the joint case), null-target metadata reflects the joint
constrained fit. The `mm_json_parse_fixed_effect_inference_table()`
parser likely already accepts these but the row-shape assertions may
need adjustment.

**P2.5. New parser for `BootstrapLikelihoodRatioTest`.** Distinct from
the contrast-row parser; produces an `mm_bootstrap_lrt` object with
`statistic`, `p_value`, `mcse`, `replicate_statistics`, plus the
reduced/alt formulas and ML-only flag. Print method labels it
"Bootstrap likelihood-ratio test (parametric, ML)".

**P2.6. New parser for cluster-resample payload.** Adds an
`mm_cluster_resample` object with `intervals` block (percentile / basic),
`group`, `cluster_resample_draw_summary`, run metadata. Print method
makes the resampled grouping factor explicit.

### Phase 3 — User-facing verbs

The API surfaces promised by `bootstrap-positioning.md`:

**P3.1. `inference_options(fit, term)`.** *The audit verb.* Returns a
data frame:

| method | expected_status | expected_reliability_reason | approx_cost | current |
|---|---|---|---|---|
| `asymptotic_wald_z` | `available` | `asymptotic_wald_z_at_boundary` | immediate | TRUE |
| `satterthwaite` | `not_assessed` | (unavailable; see reason_code) | immediate | |
| `kenward_roger` | `not_assessed` | (unavailable) | immediate | |
| `bootstrap` | `available` | `bootstrap_monte_carlo_replicates` | ~30s @ nsim=1k | |
| `bootstrap_lrt` | `available` (multi-df only) | `bootstrap_monte_carlo_replicates` | ~60s @ nsim=1k | |
| `cluster_bootstrap` | `available` (≥10 groups) | `bootstrap_monte_carlo_replicates` | ~60s @ nsim=1k | |

The expected-status column reads from the new `reliability_reason_code`
machinery — no trial-and-error. No "recommended" row, consistent with
`random_options()`. This is the verb that *operationalises* the
"audit-then-bootstrap-when-defensible" thesis.

**P3.2. `test_effect(method = "bootstrap")` — single-df.** Build the
basis vector for the term, dispatch to the existing
`mm_rust_contrast_table(method = "bootstrap")` path. Returns a
term-level `mm_effect_test` row.

**P3.3. `test_effect(method = "bootstrap")` — multi-df.** Detect
multi-df term, build multi-row `L`, dispatch through the same
contrast-bootstrap entrypoint (now multi-df-capable post-P1.3). Returns
the joint Wald/F bootstrap row.

**P3.4. `test_effect(method = "bootstrap_lrt")`.** Construct the
reduced model (drop the term), dispatch through `mm_bootstrap_lrt_json`.
Returns an `mm_effect_test` row with `method = "bootstrap_lrt"` and
the LRT-specific fields (statistic, p, mcse). Refuses on REML fits
with the same typed reason the engine uses.

**P3.5. `test_effect(method = "cluster_bootstrap")`.** Current upstream
schema `1.0.0` makes `cluster_resample` an estimator-distribution target
for intervals/diagnostics, not a certified fixed-effect p-value target.
The R surface therefore accepts the method but returns `status =
"not_assessed"`, `p_value = NA`, and `reason_code =
"bootstrap_cluster_resample_p_value_unavailable"` for single-grouping
models. Multi-grouping-factor case: refuse with `reason_code =
"cluster_bootstrap_multifactor_ambiguous"` unless the user declares
`group = "subj"` (or whichever) explicitly. Do not derive a p-value in R
from cluster-resample full-model distributions.

**P3.6. `confint(fit, parm, method = "bootstrap")`.** Routes to whichever
full-model-CI entrypoint exists per P1.4. Returns the
percentile/basic intervals from the run payload, with the run metadata
attached as attributes.

**P3.7. Update `summary()` to surface `reliability_reason_code`.** The
"Inference status:" block now shows `(term, method, status, reliability,
reliability_reason_code, reason_code)` instead of the current
`(term, method, status, reason)`.

### Phase 4 — Vignettes and benchmark

**P4.1. Rewrite Panel C of `vignettes/mixeff.Rmd`.** Currently the
panel shows `asymptotic_wald_z` / `reliability = low` / Satterthwaite
refusing. With the new machinery the panel becomes:

  1. `inference_options(fit, "days")` — show the audit map.
  2. `summary(fit)$coefficients` — show the asymptotic Wald row with
     the new `reliability_reason_code = asymptotic_wald_z_at_boundary`.
  3. `test_effect(fit, "days", method = "bootstrap", nsim = 200)` —
     show the bootstrap row, status=available, reliability_reason_code
     = `bootstrap_monte_carlo_replicates`, run payload visible
     (boundary_rate, mcse, successful_replicates).

This *enacts* the audit-first thesis: the user sees the audit (P4.1.1),
then the asymptotic answer with its labelled weakness (P4.1.2), then
the bootstrap answer with its provenance (P4.1.3).

**P4.2. New vignette: `vignettes/inference-where-lme4-says-no.Rmd`.**
Three worked examples on different data:

  1. Variance-component LRT (Stram & Lee 1994 territory) on a
     `(1 + x | g)` vs `(1 | g)` comparison via
     `compare(method = "bootstrap_lrt")`.
  2. Multi-df bootstrap on a 3-level factor or a 3-way interaction
     simple-slopes hypothesis.
  3. Cluster bootstrap on a singular fit, contrasted with the
     parametric bootstrap (which inherits the boundary covariance) on
     the same data.

Each example shows lme4's behavior, mixeff's audit, and the bootstrap
answer with its labelled reliability reason.

**P4.3. Benchmark table.** Run `examples/bench_bootstrap.rs` (Rust
side) and a parallel R script that times `lme4::bootMer`,
`pbkrtest::PBmodcomp`, and `RLRsim::exactRLRT` on the same fixtures.
Produce a small table for the new vignette and for `README.md`. If
the speedup is meaningful (≥ 5x), promote it; if marginal, keep the
table but lead with "honest" rather than "fast."

### Phase 5 — Holding pattern items

These are *not* in scope for the first wave; tracked but parked:

- **R1 user-facing diagnostic wording.** Holds at upstream
  `bd-01KQF5NBJYXGNN09HFDNMVSPX2`. Until that lands, the R wrapper can
  optionally add a thin paraphrase layer mapping `r0` → original
  formula fragment; defer until upstream wording is finalised so we
  don't fork strings.
- **Variance-component bootstrap CIs.** Out of scope per
  `upstream-bootstrap-coverage-issue.md` §Scope.
- **BCa / studentized / profile-bootstrap intervals.** Layered on top
  of the run-payload surface; deliberate follow-up.
- **GLMM bootstrap.** Depends on upstream GLMM fitting.

## Sequencing and dependencies

```
P0 (verify)            ─┐
P1 (bridge)            ─┼─►  P2 (parsers)  ─►  P3 (verbs)  ─►  P4 (vignettes)
                         │
                         └─►  P3.1 inference_options can start in parallel
                              with P1/P2 if we read existing JSON only
```

Risk-free first commits: P0.1 (rebuild + tests), P0.2/P0.3 (inspect
JSON shape), P2.1 (surface `reliability_reason_code` in the existing
parser — purely additive). These three together deliver visible value
in the smallest atomic step and de-risk the bigger bridge work.

After that the natural slice is P1.1 + P3.4 (bootstrap LRT bridge +
verb) — small, well-scoped, and unblocks the variance-component LRT
demo for the new vignette.

P3.1 (`inference_options`) is the largest user-facing payoff; can be
written immediately using *existing* JSON since the `reason_code`
fields are already populated upstream.

## Cross-references

- `planning/bootstrap-positioning.md` — strategic thesis; this plan
  is the implementation companion.
- `planning/upstream-bootstrap-coverage-issue.md` — the upstream
  request; the engineer's response largely closes B1 / B3 / B4 and
  partially B2 (cluster-CI exposed; parametric exposure to be
  confirmed in P0/P1.4).
- `planning/upstream-diagnostic-vocabulary-issue.md` — the upstream
  request; R2 / R3 / R4 done, R1 still tracked as
  `bd-01KQF5NBJYXGNN09HFDNMVSPX2` upstream.
- `planning/bootstrap_fixed_effect_payload_surface.md` — the
  long-standing downstream policy that every fixed-effect bootstrap
  p-value mixeff prints routes through the Rust contract; this plan
  follows that policy for every verb in Phase 3.
- `vignettes/mixeff.Rmd` — Panel C target for P4.1.
