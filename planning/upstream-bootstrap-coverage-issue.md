# Upstream issue draft — extend `fixed_effect_null` bootstrap coverage

Filed against `/Users/bbuchsbaum/code/rust/mixedmodels`. The upstream repo
tracks work in `mote` (per its `AGENTS.md`); the canonical form is a mote
bead. A `mote new` invocation is given at the end; the body is the issue
text.

---

## Title

`Complete fixed_effect_null bootstrap coverage: term-level joint hypotheses, full-model-distribution intervals, cluster resampling, and replicate statistics in JSON`

## Tags

`inference, bootstrap, contract, downstream-r-layer`

## Priority

`p2` — none of these block existing functionality. They are the build-out
required for `mixeff` to ship an "elegant API for p-values on degenerate
fits" — which is one of the package's headline value propositions
(`planning/vision.md`: *"refusal is a feature ... when inference is not
defensible, the software says so by name"*). On a singular fit,
asymptotic Wald and Satterthwaite are degraded or refuse; bootstrap is
the path that lets the package say "yes, here is a labelled, defensible
p-value" instead of "no certified method." Sequencing notes at the end.

## Body

### Origin

Surfaced while building the bootstrap surface for `mixeff`. The contrast-
level bootstrap path (`mm_fixed_effect_bootstrap_contrast_json` →
`LinearMixedModel::fixed_effect_null_bootstrap_inference_table`) is in
place and produces clean, labelled p-values on a deliberately-singular
fit (`rt ~ days + (1 + days | subj)`, n = 18, near-perfectly-correlated
intercept and slope). Confirmed end-to-end in R:

```r
contrast(fit, c(0, 1), method = "bootstrap",
         bootstrap = bootstrap_control(nsim = 200, seed = 1))$table[
  , c("contrast", "estimate", "p_value", "method", "status", "reliability")]
#>   contrast estimate     p_value    method    status reliability
#> 1       c1 9.250949 0.004975124 bootstrap available         low
```

The contract document for this surface
(`docs/bootstrap_fixed_effect_contract.md`) explicitly lists four items
under "Later scope":

> - multi-df fixed-effect hypotheses
> - model-comparison/bootstrap LRT rows
> - adaptive replicate escalation
> - parallel execution
> - GLMM bootstrap calibration

Three of those (multi-df, LRT, parallel) plus two further items
(full-model-distribution intervals; cluster resampling; per-replicate
statistics in JSON) are what the R wrapper needs to round out
`test_effect(method = "bootstrap")`, `confint(method = "bootstrap")`,
and a future cluster-bootstrap method label. This issue groups the
ones the R wrapper would benefit from most, in priority order.

### Scope

In scope: bootstrap inference for *fixed-effect* coefficients,
contrasts, and term-level joint hypotheses on Gaussian LMMs.

Out of scope (each deserves a separate request):

- **Variance-component bootstrap CIs.** Percentile / basic intervals on
  `theta`, `sigma`, or random-effect variances near zero behave badly;
  honest reporting needs a labelled interval kind plus boundary-mass
  / coverage commentary. We will request these explicitly when the
  fixed-effect surfaces below are stable.
- **BCa / studentized / profile-bootstrap intervals.** These are
  legitimate next-generation interval kinds but are layered on top of
  the run-payload surface this issue requests. Out of scope here.
- **GLMM bootstrap calibration.** Listed as "Later scope" in the
  contract; depends on upstream GLMM fitting. Out of scope here.
- **R-side reconstruction of any p-value or interval from full-model
  refit summaries.** Per
  `mixeff/planning/bootstrap_fixed_effect_payload_surface.md`, every
  fixed-effect bootstrap p-value and interval mixeff prints will route
  through this Rust contract, not be reassembled R-side. The shape of
  this issue follows from that policy.

### Reproducer

```r
set.seed(3)
n_subj <- 18L
days   <- 0:9
b0     <- rnorm(n_subj, sd = 30)
b1     <- 0.5 * b0 / 30 * 10 + rnorm(n_subj, sd = 0.5)

sleep_like <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
  data.frame(subj = factor(i), days = days,
    rt = 250 + b0[i] + (10 + b1[i]) * days +
         rnorm(length(days), sd = 20))
}))

fit <- lmm(rt ~ days + (1 + days | subj), sleep_like,
           control = mm_control(verbose = -1))
# fit_status(fit)  -> "converged_reduced_rank"
# is_singular(fit) -> TRUE
```

`lme4::lmer` returns `isSingular(m) == TRUE` on the same data, so the
singularity is intrinsic to the design.

---

### B1. Term-level joint null bootstrap (multi-df)

#### Observed

```r
test_effect(fit, "days", method = "bootstrap")
#> term: days, method: bootstrap, status: not_assessed
#> reason: parametric bootstrap fixed-effect inference requires a certified
#>   fixed_effect_null bootstrap payload; call test_contrast_with_bootstrap_payload
#>   with replicate accounting, failed-refit policy, Monte Carlo uncertainty,
#>   and reproducibility state
```

The single-df case (`days` is a 1-df numeric term) can be handled R-side
by routing through the existing contrast entrypoint. The multi-df case
(factor with k > 2 levels, an interaction, a polynomial term, etc.) needs
a Rust-owned joint null payload because the joint statistic depends on
the joint covariance under the constrained null model.

Note this is **not** the same as the bootstrap likelihood-ratio test that
the R wrapper can already produce R-side via `compare(reduced, full,
method = "bootstrap")`. That LRT path drops the term, refits, and
compares deviances — a different statistic with a different reference
distribution. B1 asks for the *Wald-joint* alternative: simulate under
the joint constrained null, evaluate the same Wald-form joint statistic
on each replicate, return one term-level inference row of the same
shape `test_effect()` returns today. The two paths complement each
other and the R wrapper would expose both, labelled distinctly
(`method = "bootstrap"` vs `method = "bootstrap_lrt"`).

#### Requested

Extend `fixed_effect_null_bootstrap_inference_table` (or add a sibling
`fixed_effect_null_bootstrap_term_inference_table`) so it accepts a
multi-row `L` matrix and produces:

- a Wald-type joint statistic under each replicate (e.g. quadratic form
  `(L beta_b - rhs)' (L vcov_beta_b L')^+ (L beta_b - rhs)` against
  the observed value);
- an exact-rank-aware reference distribution under the joint null
  (constrained refit on each simulated response);
- a single `bootstrap` inference row of shape compatible with the
  existing scalar `mixedmodels.fixed_effect_inference_table` row,
  reporting `numerator_df`, joint statistic, p-value, MC-SE, and the
  same run metadata fields that exist today (boundary_rate,
  failed_refit_policy, etc.).

The R wrapper then dispatches `test_effect(fit, term, method = "bootstrap")`
to scalar-bootstrap for 1-df terms and joint-bootstrap for multi-df terms
under the same surface; the user writes one verb regardless of arity.

#### Acceptance criteria

- A new bridge entrypoint accepts `L` with `nrow(L) >= 1` and returns a
  certified bootstrap inference row for the joint hypothesis
  `L beta = rhs`.
- For `nrow(L) == 1` the result is bit-equal (within seed & tolerance)
  to the existing scalar contrast bootstrap.
- The `null_target` payload metadata includes the constrained null fit's
  `theta`, `sigma`, and `coefficient_count` consistent with the scalar
  case.
- Snapshot fixture covers a 2-df numeric term (e.g. polynomial) and a
  k-level factor.

---

### B2. `full_model_distribution` target exposed for bootstrap CIs

#### Observed

The contract (`docs/bootstrap_fixed_effect_contract.md` "Bootstrap
Targets") declares two simulation targets:

| Target | Meaning | May produce p-value? |
|---|---|---|
| `full_model_distribution` | Simulate from the fitted model. | No |
| `fixed_effect_null` | Simulate from the constrained null. | Yes |

Only `fixed_effect_null` is currently exposed at the bridge level. The
contract calls out that `full_model_distribution` is "useful for
replicate distributions, percentile intervals, diagnostics, and smoke
tests" — exactly what `confint(fit, parm, method = "bootstrap", ...)`
needs in the R wrapper. R cannot synthesize percentile intervals from
the null-bootstrap output without re-deriving them under the alternative
model, which would duplicate engine-owned simulation logic.

Doing this honestly R-side is not "take quantiles of refit estimates."
The contract for the null-target run payload already mandates
`failed_refits`, `failed_refit_policy`, `boundary_count`,
`successful_replicates`, `boundary_rate`, `seed_record`, `mcse`, and a
non-trivial validity envelope (`successful_replicates` thresholds, MCSE
on tail areas, etc.). A defensible CI from a `full_model_distribution`
run inherits the same machinery — there is no honesty contract for an
R-side `take quantiles(coef$days for sim in 1:nsim)` that wouldn't
duplicate this. Exposing the target at the bridge keeps that machinery
single-source.

#### Requested

Add a bridge entrypoint that runs simulations from the `full_model_distribution`
target and returns the resulting per-coefficient (or per-contrast) replicate
distributions, plus percentile, basic, and (where applicable) BCa interval
endpoints. The shape can mirror the existing run payload, with `target_kind
= "full_model_distribution"` and a clearly-labelled `intervals` block:

```
intervals: [
  { contrast: "c1", level: 0.95, kind: "percentile", lower: 7.42, upper: 11.10 },
  { contrast: "c1", level: 0.95, kind: "basic",      lower: 7.40, upper: 11.10 }
]
```

#### Acceptance criteria

- New bridge entrypoint produces a `BootstrapRunPayload` with
  `target_kind = "full_model_distribution"` and an `intervals` block.
- The payload also surfaces `replicate_statistics` (see B4 below) so the
  R wrapper or other clients can recompute or audit the intervals.
- Snapshot fixture covers a single-coefficient and a two-coefficient
  contrast.
- The payload's `notes` repeats the contract caveat: full-model
  distribution does **not** certify hypothesis-test p-values; it is
  for intervals and diagnostics only.

---

### B3. Group / cluster resampling as a third simulation target

#### Observed

Both existing targets (`full_model_distribution`, `fixed_effect_null`)
are *parametric*: they simulate a new response vector from a fitted
generative model. Neither is robust to misspecification of the
random-effects covariance — which is exactly the case on a singular fit,
where the parametric bootstrap inherits the boundary covariance into
every replicate (the reproducer above produces `boundary_rate ≈ 0.36`
even at `nsim = 50`).

A cluster-resampling target — sample observed *clusters* (grouping
levels) with replacement, refit on the resampled dataset — yields a
non-parametric bootstrap for the random-effects part, while leaving the
parametric fixed-effect estimator intact. This is one of the standard
distribution-free alternatives in the mixed-model bootstrap literature
(e.g. Field & Welsh 2007 *Bootstrapping clustered data*).

Validity is conditional, not universal. Cluster bootstrap is credible
only when clusters are exchangeable and numerous enough that the
empirical distribution of clusters is a reasonable proxy for the
population. With few clusters (n_groups in the single digits) the
bootstrap distribution itself is coarse; with crossed grouping factors
("subjects" and "items" both random) the resampling unit is ambiguous
and a single-factor cluster bootstrap is not automatically the right
target. The contract should make these scope conditions explicit so the
R wrapper can label them in `reliability_reason` when fired.

Cluster resampling is **not** "simulate-and-refit on resampled
subjects" — it is mechanically distinct in three places that matter:

1. **Cluster relabeling.** When a level is sampled twice, its rows
   need fresh, distinct group keys in the resampled dataset; otherwise
   the model treats two duplicate clusters as one large cluster and
   the bootstrap target is silently the wrong target. The bridge has
   to mint and record stable replicate-local labels.
2. **Multi-grouping-factor disambiguation.** For `(... | a) + (... | b)`
   or crossed grouping factors, "the cluster" is ambiguous. The user
   has to declare which factor is being resampled, and the run payload
   has to record it; resampling one factor while leaving the other
   intact has different inference properties from resampling both.
3. **Inference contract.** Cluster resampling produces a pivotal
   distribution under different conditions than parametric bootstrap;
   for hypothesis testing it requires either centering the statistic
   at the observed value (basic bootstrap) or simulating under a
   constrained null built from the resampled clusters. The contract
   should specify which of these the engine implements and how p-values
   are computed.

#### Requested

Add a third bootstrap target — `cluster_resample` (name negotiable) —
that resamples grouping levels with replacement for one named grouping
factor, refits the LMM on the resampled dataset, and emits the same
run-payload shape as the existing targets. The user's input is the
grouping factor name and the existing `requested_replicates`,
`failed_refit_policy`, `seed`.

For a multi-grouping-factor model, the user picks which grouping factor
to resample over (typically the highest-level factor — subjects in the
reproducer above). The contract should be explicit that
`cluster_resample` is valid both for hypothesis testing
(distribution-free under the null) and interval estimation, and is
particularly recommended when the parametric fit is at a covariance
boundary.

#### Acceptance criteria

- New `BootstrapTarget` variant `ClusterResample { group: String }`
  (or equivalent) added to the closed enum.
- Bridge entrypoint accepts a grouping-factor name; payload reports
  the resampled group, the number of distinct levels per replicate,
  and (optionally) a `level_resample_distribution` summary for
  reproducibility.
- Resampled-cluster relabeling is engine-owned and reproducible: when
  level `L` is sampled `k > 1` times, the engine mints `k` stable
  replicate-local relabels; the run payload records the relabeling
  scheme so two replicates run with the same `seed_record` produce
  identical relabelings.
- The contract makes the inference choice explicit: the run payload
  reports the pivotal-distribution policy (basic-bootstrap centering
  vs. constrained-null cluster bootstrap) and the p-value computation
  method matches.
- Validity is documented in `bootstrap_fixed_effect_contract.md` —
  this is a hypothesis-test-capable target alongside `fixed_effect_null`.
- Snapshot fixtures: (a) single-grouping-factor model, cluster bootstrap
  on the reproducer above produces a finite, reproducibility-stable
  p-value; (b) two-grouping-factor model, user resamples the
  higher-level factor, with the unresampled factor preserved.

---

### B4. Per-replicate statistics in the bootstrap JSON

#### Observed

The bootstrap row currently exposes (in `details$bootstrap` for the
scalar-contrast case):

```r
boot <- ct$table$details[[1]]$bootstrap
names(boot)
#>  [1] "target_kind"          "target_label"         "contrast_label"
#>  [4] "requested_replicates" "completed_replicates" "successful_replicates"
#>  [7] "failed_refits"        "failed_refit_policy"  "boundary_count"
#> [10] "boundary_rate"        "seed_rng"             "seed"
#> [13] "finite_statistic_count" "mcse"               "null_target"
```

Notably absent: `replicate_statistics` — the per-replicate vector of
the bootstrap test statistic. The contract
(`docs/bootstrap_fixed_effect_contract.md` "Run Payload" table) lists
this field as part of the run payload schema:

> `replicate_statistics` | finite/non-finite bootstrap statistic values
> or a durable reference to them

The current bridge JSON does not include it. Without per-replicate
statistics, R cannot:

- recompute MC-SE under a different policy (e.g. for tail-area p-values
  near 0 or 1, where the standard MC-SE formula needs care);
- compute BCa intervals or other resample-summary intervals;
- plot the bootstrap distribution for the user;
- run downstream sanity checks (e.g. distribution shape, outlier
  detection on the simulated null).

#### Requested

Populate the `replicate_statistics` field in the run payload with the
finite (and, where useful, non-finite) per-replicate statistic values.
For large `nsim`, a "durable reference" — e.g. an embedded sidecar
identifier — is acceptable per the contract; the simplest realization
is to inline the vector for `nsim <= some_threshold` (say 10 000) and
fall back to a referenced sidecar above that threshold.

#### Acceptance criteria

- `replicate_statistics` is populated in the JSON returned by
  `mm_fixed_effect_bootstrap_contrast_json` (and any joint-test
  successor from B1; and the full-model-distribution entrypoint from
  B2).
- Vector length equals `successful_replicates` for the inline case;
  non-finite values use a stable JSON sentinel (or a parallel
  `non_finite_indices` list).
- Snapshot fixture: scalar contrast bootstrap on the reproducer with
  `nsim = 50`; the parsed JSON contains a 50-length `replicate_statistics`
  vector that recovers `mcse` to within numerical tolerance when R
  computes it independently.

---

### Sequencing recommendation

**B4** is the smallest, lowest-risk, and unblocks the most R-side
downstream work (MC-SE recomputation, plotting, BCa). Land it first.

**B2** (`full_model_distribution` intervals) is well-bounded — the
target is already named in the contract; the work is exposing it at the
bridge level. Land it second; it directly enables
`confint(fit, parm, method = "bootstrap", nsim, seed)` in the R wrapper.

**B1** (multi-df joint null) is the largest of the four because it
requires constrained refit machinery for joint hypotheses. Schedule
after B2 and B4. The R wrapper can ship `test_effect(method = "bootstrap")`
for the single-df case before B1 lands; B1 then upgrades it to all-cases.

**B3** (cluster resample) is independent of the other three and can be
scheduled in parallel. It is the most semantically novel — adds a third
simulation target — but the most-recommended in the literature for
singular or boundary fits.

### References

- `docs/bootstrap_fixed_effect_contract.md` — base contract; the
  "Bootstrap Targets" and "Run Payload" tables are the anchor for B2
  and B4 respectively.
- `docs/bootstrap_fixed_effect_contract.md` "Later scope" list (lines
  31–37) — already names B1 (multi-df) and adjacent items.
- `src/model/linear.rs` — current
  `fixed_effect_null_bootstrap_inference_table` entrypoint
  (around line 5253 as of writing).
- `src/model/linear.rs` — `FixedEffectBootstrapOptions` struct
  (around line 10189 as of writing); will need a `target` and possibly a
  `cluster_group` field for B2/B3.
- `mixeff/planning/bootstrap_fixed_effect_payload_surface.md` — the
  downstream-side design plan, which already commits to
  `contrast(method = "bootstrap")` as the surface and waits on B1/B2/B3
  for full coverage.
- `mixeff/vignettes/mixeff.Rmd` — the four-panel "Why mixeff?" vignette
  on the singular fit; the headline panel for inference is currently
  asymptotic Wald + reliability=low because that is the only
  *immediately-available* labelled p-value. The bootstrap path is the
  defensible one we want to advertise once the API gaps close.

---

## Mote invocation

```sh
# Run from /Users/bbuchsbaum/code/rust/mixedmodels
mote new "Complete fixed_effect_null bootstrap coverage: term-level joints, full-model-distribution intervals, cluster resampling, replicate statistics" \
  --priority 2 \
  --body "$(cat /Users/bbuchsbaum/code/mixeff/planning/upstream-bootstrap-coverage-issue.md)"
```

Once filed, link the resulting bead in
`mixeff/planning/bootstrap_fixed_effect_payload_surface.md` so the
downstream design plan tracks the upstream dependency.
