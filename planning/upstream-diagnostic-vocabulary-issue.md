# Upstream issue draft — user-facing diagnostic & inference vocabulary

Filed against `/Users/bbuchsbaum/code/rust/mixedmodels`. The upstream repo
tracks work in `mote` (per its `AGENTS.md`); the canonical form is a mote
bead. A `mote new` invocation is given at the end; the body is the issue
text.

---

## Title

`Surface user-facing wording and stable codes in audit, certificate, and inference responses`

## Tags

`audit, inference, contract, diagnostics, downstream-r-layer`

## Priority

`p2` — none of these block downstream functionality, but together they
determine whether the `mixeff` "Why mixeff?" vignette can credibly claim
"every printed number carries provenance" to a working scientist (rather
than a contract author). Sequencing notes at the end.

## Body

### Origin

These four asks surfaced while rewriting the `mixeff::vignette("mixeff")`
"Why mixeff?" page as a four-panel poster on a deliberately-singular
fit (`rt ~ days + (1 + days | subj)`, n = 18, near-perfectly-correlated
intercept and slope). The vignette is meant to advertise that the package
"names the convergence outcome" and that "every reported number carries
its method." When the four panels are read by a working scientist
rather than by a contract author, three kinds of vocabulary leak through
to the user that should not:

1. Term and parameter identifiers that only the engine knows
   (`r0`, `theta[2]`, `varpar[2]`).
2. An implementation-detail error message in place of an existing stable
   reason code (Satterthwaite at boundary; the contract already names the
   code, but the R-facing JSON exposes the finite-difference prose).
3. A `reliability` label without the principled reason that justifies the
   label.

A fourth ask is investigative: a single semantic event (a reduced-rank
covariance) is reported three times in the artifact, once per pipeline
stage that observes it.

The R layer can paraphrase, but the contract is clear that these strings
should be authored once, in Rust, and mirrored downstream
(`compiler_contract_v0_prd.md`; `random_term_card_prd.md` FR4
"Per-block English authorship"). We would rather hold the line than
fork wording in the R wrapper.

### Reproducer

In R, with `mixeff` loaded:

```r
set.seed(3)
n_subj <- 18L
days   <- 0:9
b0     <- rnorm(n_subj, sd = 30)
b1     <- 0.5 * b0 / 30 * 10 + rnorm(n_subj, sd = 0.5)

d <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
  data.frame(
    subj = factor(i),
    days = days,
    rt   = 250 + b0[i] + (10 + b1[i]) * days +
           rnorm(length(days), sd = 20)
  )
}))

fit <- lmm(rt ~ days + (1 + days | subj), d,
           control = mm_control(verbose = -1))
# fit_status(fit)  -> "converged_reduced_rank"
# is_singular(fit) -> TRUE
```

The same data fitted by `lme4::lmer` returns `lme4::isSingular(m) == TRUE`,
so the singularity is intrinsic to the design, not an artefact of the
mixedmodels engine. The four observations below are taken verbatim from the
R-facing surfaces (`diagnostics()`, `inference_table()`, `test_effect()`,
`as_json(fit)$artifact_json`) on this fit.

---

### R1. Human-readable wording in `DiagnosticReport.message` and certificate boundary text

#### Observed

```r
diagnostics(fit)$table[, c("code", "message")]
#>                 code                                                            message
#> 1 covariance_reduced fitted covariance for r0 has effective rank 1 of requested rank 2
#> 2 boundary_parameter                                    theta[2] is on its lower bound
#> 3 covariance_reduced fitted covariance for r0 has effective rank 1 of requested rank 2
```

The user wrote `(1 + days | subj)`. They have no map from `r0` to that
fragment, no map from `theta[2]` to "the variance parameter for the random
`days` slope", and no narrative for what "rank 1 of requested rank 2"
means for *their* model. The structured fields (`term_id`, parameter
indices) are correct and useful for tooling; the issue is that the
`message` strings are written in engine vocabulary.

#### Requested

For every diagnostic emitted by `audit_design()` and the optimizer
certificate, the `message` field should cite (a) the formula fragment the
diagnostic refers to and (b) the parameter's role in user vocabulary, in
addition to whatever structured fields already exist. Suggested rewordings
on the reproducer above:

```
covariance_reduced  | The random-effects covariance for `(1 + days | subj)`
                      collapsed from rank 2 to rank 1: the intercept and
                      `days` slope are estimated as perfectly correlated,
                      which is the boundary of the parameter space.

boundary_parameter  | The variance parameter for the random `days` slope
                      is at its lower bound 0; the slope is being estimated
                      as having no variation across subjects.
```

The mapping from `term_id` to `original_fragment` already exists in the
`RandomTermCard` schema (`random_term_card_prd.md` §FR2). The mapping
from `theta_index` to a user-facing role is one block away — a thin layer
that, given a `theta[i]`, names the term it parameterises and whether it
is a variance, a Cholesky off-diagonal, or a residual.

#### Acceptance criteria

- Every diagnostic emitted by `audit_design()` and `optimizer_certificate`
  has a `message` field whose text references either the user's formula
  fragment, the user-facing parameter role, or both. No bare `r0` /
  `theta[i]` appears in `message` text.
- Internal identifiers (`term_id`, `theta_index`, parameter slot) remain
  in their existing structured fields. Tooling that joins on those fields
  is unaffected.
- Wording is authored once in Rust, per `random_term_card_prd.md` §FR4,
  and round-trips through JSON unchanged.
- A snapshot fixture at the level of `audit_design()` for a singular fit
  asserts the new wording; the R wrapper consumes it without
  post-processing.

---

### R2. Surface existing stable Satterthwaite reason codes; demote finite-difference prose

#### Observed

```r
te <- test_effect(fit, "days", method = "satterthwaite")
te$table[, c("term", "method", "status", "reason")]
#>   term        method       status
#> 2 days satterthwaite not_assessed
#>                                                                                                                                                                                                     reason
#> 2 Satterthwaite fixed-effect inference could not compute vcov_beta derivatives:
#>     Invalid argument: cannot compute central finite-difference derivative for varpar[2]:
#>     value is at or too near lower bound 0
```

This is a *real* limitation, not a bug: Satterthwaite df require
derivatives of the variance parameters, and at a boundary fit the
derivative for at least one variance on the boundary is not defined.
However, `docs/satterthwaite_scalar_contract.md` already names two
relevant stable codes:

- `satterthwaite_varpar_deviance_unavailable`
- `satterthwaite_varpar_covariance_unavailable`

and the same document calls out (line 290) that "a boundary or singular
case [is] expected to produce low reliability or a structured unavailable
reason depending on derivative stability." On the reproducer, the
unavailable response that reaches R has neither a stable code nor a
structured reason — it has the finite-difference prose. Tooling cannot
dispatch on it; the R wrapper cannot translate it into a user-facing
explanation; future Kenward-Roger work will inherit the same shape.

#### Requested

For Satterthwaite (and any other inference method that can hit the same
class of failure), the response that reaches the R wrapper should carry:

- a `code` field drawn from the closed inference-unavailability enum.
  If a more specific variant is warranted than the two existing codes —
  e.g. `satterthwaite_unavailable_at_boundary` — please add it; the
  current case is principled and recurring enough to deserve its own
  variant.
- a `description` field with the *principled* one-sentence reason
  (e.g. "Satterthwaite degrees of freedom require derivatives of the
  variance parameters with respect to the data; at a boundary fit, the
  derivative for at least one variance on the boundary (here, the
  random `days` slope) is not defined").
- the existing finite-difference text demoted to a `detail` field on the
  same response, retained for debugging but no longer the headline.

The same shape generalizes naturally to upcoming methods —
`kenward_roger_unavailable_at_boundary`,
`kenward_roger_unavailable_for_glmm`, etc. — and matches the
"refusal-as-feature" decision tree in
`mixed_model_compiler_inference_contract.md` §"Refusal vs
ConvergedPenalised".

#### Acceptance criteria

- The existing Satterthwaite stable codes appear in the inference response
  for the reproducer above (or a new variant, if warranted).
- A `description` field is populated with the principled reason; it does
  not contain the words "central finite-difference" or "varpar[i]".
- The current finite-difference text remains accessible under a `detail`
  (or `notes`) field on the same response.
- A test fixture asserts the new shape on the boundary case.

---

### R3. `reliability_reason` on every inference row

#### Observed

```r
inference_table(fit)$table[, c("term", "method", "status", "reliability", "reason")]
#>          term            method    status reliability reason
#> 1 (Intercept) asymptotic_wald_z available         low   <NA>
#> 2        days asymptotic_wald_z available         low   <NA>
```

`reliability` is calibrated information, but the `reason` column is
populated only on `unavailable`/`not_assessed` rows. A reader who sees
`low` cannot tell *why* the row is low: was it asymptotic Wald applied to
a boundary fit (so SEs are likely understated)? Was t the requested
distribution but df were unavailable (so z was substituted)? Was it the
`covariance_too_rich` low-grouping-levels case? The verdict travels; the
warrant does not.

The R wrapper can paraphrase, but only by guessing — the upstream knows
which rule fired.

#### Requested

Extend the inference contract so every inference row has access to a
`reliability_reason` field — drawn from a closed enum of stable codes —
populated regardless of `status`. Suggested initial variants:

- `asymptotic_wald_z_at_boundary` — asymptotic Wald applied where the
  fit's covariance is reduced-rank or one or more variance parameters are
  on the boundary; SEs typically understate uncertainty.
- `degrees_of_freedom_unavailable_so_z_substituted` — t was the requested
  distribution but df were unavailable; standard normal substituted.
- `low_grouping_levels` — fewer grouping levels than the contract's
  threshold for full-covariance reliability (already cited in
  `covariance_too_rich`).
- `interior_converged_well_specified` — the "this is fine" reason, kept
  explicit so that a row labelled `high` carries the warrant alongside
  the verdict (no silent assumption that an absent reason means good).

This mirrors the existing `unavailable_reason` field for unavailable
rows; it adds the same shape for available rows.

#### Acceptance criteria

- Every row produced by the inference surface has a `reliability_reason`
  populated, drawn from a closed enum.
- The mapping from (convergence status, requested method, fit features)
  to `reliability_reason` is documented in
  `mixed_model_compiler_inference_contract.md` (or a new section).
- Snapshot fixtures at the inference level cover at least three rows:
  `interior_converged_well_specified`, `asymptotic_wald_z_at_boundary`,
  `degrees_of_freedom_unavailable_so_z_substituted`.

---

### R4. Disambiguate or dedupe `covariance_reduced` reported by multiple pipeline stages

#### Observed

```r
diagnostics(fit)$table[, c("code", "stage", "message")]
#>                 code         stage                                              message
#> 1 covariance_reduced certification fitted covariance for r0 has effective rank 1 ...
#> 2 boundary_parameter certification theta[2] is on its lower bound
#> 3 covariance_reduced certification fitted covariance for r0 has effective rank 1 ...
```

A walk over the artifact JSON shows the `covariance_reduced` code
appearing in three locations: `reductions[…].diagnostics`,
`optimizer_certificate.diagnostics`, and the top-level `diagnostics`
array. The R wrapper aggregates faithfully and shows two of the three.
On a single semantic event ("the random-effects covariance collapsed
from rank 2 to rank 1"), this looks like noise to a reader.

#### Requested

Either:

- **Dedupe at the artifact level.** A single semantic event produces a
  single diagnostic, attached at the most informative stage (likely
  `optimizer_certificate`, since that is when the rank reduction is
  *observed*; design-time can cite it but does not need to re-emit).
- **Or distinguish.** If the three emissions are *intended* to carry
  different information (e.g. the design-time variant is "this is
  flagged as at risk pre-fit", the certificate variant is "this actually
  happened"), they should carry different codes (`covariance_reduced`
  vs `covariance_reduction_observed`) or a discriminating field
  (`source: "audit_design" | "optimizer_certificate"`) so the
  duplication is meaningful.

#### Acceptance criteria

- Either: each (`code`, `term_id`, `message`) triple appears at most
  once across the artifact's diagnostic surfaces;
- Or: surfaces that report related events use distinct codes or carry a
  `source` field, and the reproducer above produces non-redundant rows
  in `diagnostics()`.

---

### Sequencing recommendation

R4 (investigate) is small and likely either trivial or reveals a real
design clarification. Land it first.

R2 (Satterthwaite codes) is well-scoped — the codes already exist in
the contract; the work is propagating them to the R-facing JSON and
demoting the prose. Land it second.

R3 (reliability_reason) is an additive contract field; lands cleanly
once the enum is decided.

R1 (human-readable wording) is the broadest user-experience improvement
and likely the largest implementation diff, since it touches every
diagnostic site that currently emits engine-vocabulary `message` strings.
Land it last but treat it as the headline win — the four-panel
"Why mixeff?" vignette is unconvincing without it.

### References

- `compiler_contract_v0_prd.md`
- `mixed_model_compiler_inference_contract.md` §"Inference Contract",
  §"Refusal vs ConvergedPenalised decision tree"
- `random_term_card_prd.md` §FR2 `term_id`/`original_fragment`,
  §FR4 "Per-block English authorship"
- `satterthwaite_scalar_contract.md` lines ~250–290 (existing stable
  codes; boundary-case acceptance language)
- `mixeff/vignettes/mixeff.Rmd` (the four-panel poster that surfaced
  these gaps)

---

## Mote invocation

```sh
# Run from /Users/bbuchsbaum/code/rust/mixedmodels
mote new "Surface user-facing wording and stable codes in audit, certificate, and inference responses" \
  --priority 2 \
  --body "$(cat /Users/bbuchsbaum/code/mixeff/planning/upstream-diagnostic-vocabulary-issue.md)"
```

Once filed, link the resulting bead in `mixeff/planning/PRD.md` §5.2.4
(stable diagnostic codes) and §5.2.5 (random term card) so the
downstream contract notes the dependency.
