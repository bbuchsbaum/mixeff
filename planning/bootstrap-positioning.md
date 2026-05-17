# Bootstrap positioning — strategic thesis for `mixeff`

Status: strategic note; informs vignette priorities, the upstream
bootstrap-coverage issue, and the R-side bootstrap API plan.

## Thesis

> mixeff should expose when asymptotics are weak, and then offer a fast
> audited bootstrap path when that is the defensible route.

The order of operations is load-bearing. *First* the package tells the
user when asymptotic methods are unreliable on their fit; *then* the
package offers bootstrap as the defensible alternative — when it is in
fact defensible. The framing leaves room for a third outcome: cases
where neither asymptotic Wald nor bootstrap is the right answer (the
fit is so degenerate that no inference is honest), where the package's
job is to refuse cleanly. Bootstrap is the answer "when that is the
defensible route," not always.

Bootstrap converts a large class of mid-90s-vintage lme4 pain points
from "unavailable or awkward" into "available, labeled, auditable."
Bootstrap does not make a model true, rescue weak identification, or
remove small-sample limits. It makes the defensible path ergonomic
when the defensible path runs through a bootstrap.

Bootstrap is positioned as a first-class inference route alongside
Wald and Satterthwaite/KR, not a fallback trick — but it is offered
*conditionally on the audit*, not unconditionally.

Speed is what enables the pitch. A bootstrap that takes ten minutes for
n = 200 is an "ok if you really want it" feature; a bootstrap that takes
five seconds is a default. The Rust engine is the only thing that can
make this the default story — pure-R bootstrap (via lme4 + simulate +
refit) is too slow to be the headline. If the upstream Rust bootstrap
is meaningfully faster than `lme4::bootMer` / `pbkrtest::PBmodcomp` /
`RLRsim` (the working-scientist baselines) then mixeff has a positioning
that lme4 cannot answer.

## Why this matters for mixeff specifically

lme4's most-cited stance is "no p-values" — rooted in the genuine
unreliability of asymptotic methods on the kinds of fits applied users
write down. mixeff has already committed to a different stance via its
"no fake certainty" principle: when a method is unavailable, return `NA`
with a stable reason. Bootstrap closes the loop on this principle by
giving the package a *labelled, defensible inference route* that works
where Wald-asymptotic and Satterthwaite/KR refuse — without weakening
the honesty contract.

The headline lme4 user complaint becomes:

> "lme4 says no p-values because asymptotics are weak. mixeff exposes
> when asymptotics are weak, and offers a fast audited bootstrap path
> when that is the defensible route — and refuses cleanly when neither
> is."

That is the inversion worth pursuing.

## Tier A — bootstrap converts unavailable/awkward into available + labeled

**1. Fixed-effect p-values, including in singular fits.**
Asymptotic Wald z is biased downward at boundary fits; Satterthwaite
refuses (the derivative does not exist on the boundary); Kenward-Roger
inherits the same problem. A model-conditioned bootstrap p-value (under
either the constrained null or the fitted model) does not need
derivatives at the variance parameter. The user gets a labelled p-value
where lme4 prints t and refuses to call it more than that. *Caveat: a
parametric bootstrap p-value is conditional on the fitted model;
"calibrated tail area" should be qualified, not asserted.*

**2. Variance-component tests** (*"is the random slope ≠ 0?"*).
The asymptotic chi-square LRT is wrong at the boundary of the parameter
space (Self & Liang 1987, Stram & Lee 1994 — 50:50 mixture, etc.).
Currently solved by `RLRsim`, `pbkrtest::PBmodcomp`, or hand-rolled
loops. Fast Rust bootstrap LRT replaces all three at the price of a
single labelled method. *Caveat: variance-component nulls are
boundary/nonstandard; even bootstrap LRT often has low power on small
group counts. The package should label the test as
boundary-of-parameter-space and report the run-payload boundary rate.*

**3. Confidence intervals on fixed-effect coefficients and contrasts.**
`confint(merMod, method = "profile")` is famously slow and sometimes
refuses; `confint(method = "boot")` is even slower. Fast Rust bootstrap
CIs (percentile, basic) on the fixed-effect side resolve the practical
problem. *Caveat: variance-component CIs are a separate territory —
percentile/basic intervals near zero variance behave badly; honest
reporting needs labelled interval kind and boundary-mass commentary.
Out of scope for the first wave; tracked as a follow-up.*

**4. Complex multi-df fixed-effect contrasts** — simple slopes in a
3-way interaction, polynomial-component tests, joint hypotheses on
factor terms with k > 2 levels. Bootstrap covers all of these once the
upstream multi-df null target lands (request B1 in
`upstream-bootstrap-coverage-issue.md`). *Caveat: the bootstrap
statistic and null target must be certified, not assembled casually
from full-model refits — which is why this is an upstream contract
ask, not an R-layer build.*

**5. Power / sample-size analysis** (`simr`-style).
This may be the largest product-level upside. simr's whole architecture
is *literally* simulate → refit → record-p, repeated nsim times. A
fast Rust engine turns 8-hour power analyses into 5-minute ones — and
that is genuinely infrastructure-level leverage, not just a faster
inference verb. lme4 has no answer beyond `simr`. mixeff could host
the next-generation answer if the engine is fast enough.

## Tier B — strong leverage, not a silver bullet

**6. Convergence / optimizer-instability diagnostics.**
Bootstrap turns optimizer fragility into a *distribution* of fragility:
the run payload's `boundary_rate`, `failed_refits`, `boundary_count`,
optimizer codes. "The optimizer warned once on the original fit"
becomes "8 of 200 bootstrap refits also warned, here's the spread of
estimates." That is real, defensible evidence — but it is evidence
*about robustness*, not a fix for a misspecified model.

**7. Robustness to non-normality / mis-specified residuals.**
Cluster (subject) bootstrap is non-parametric for the random-effects
part: it resamples observed clusters rather than simulating from a
fitted model. Parametric bootstrap from a Gaussian LMM still inherits
the Gaussian assumption. *Caveat (important): cluster bootstrap is
credible only when clusters are exchangeable and numerous enough that
the empirical distribution of clusters proxies the population. With
few clusters the bootstrap distribution itself is coarse; with crossed
grouping factors ("subjects" and "items" both random) the resampling
unit is ambiguous and a single-factor cluster bootstrap is not
automatically the right target.*

**8. Random-effects structure decisions** ("keep it maximal" vs.
parsimonious). Bootstrap LRT gives a valid test for "should the random
slope variance be in the model" where asymptotic LRT fails. But the
test is one piece of evidence, not a model-selection oracle. The
package's "no model selection or random-effects recommendation engine"
non-goal (`PRD.md` §3) still holds; bootstrap LRT exposes the question
honestly, doesn't answer it for the user.

**9. GLMM inference.** Bootstrap is often the right *direction* for
GLMM p-values — quasi-likelihood / PQL has known biases, and asymptotic
GLMM inference is more fragile than LMM inference. Calling it "the gold
standard" overclaims: separation, rare events, overdispersion, zero
inflation, and misspecified random-effects structure can dominate
whatever inference method is chosen. Downstream of upstream GLMM
fitting; the architecture should anticipate it without claiming a
pre-existing answer.

## Tier C — bootstrap is enabling, not the protagonist

These are mixeff's *other* pillars; bootstrap doesn't help with them
but doesn't compete either:

**10. "What does my formula assume?"** — `compile_model()`,
`explain_model()`, `random_options()`, the *audit-before-fit* pillar.
**11. "Is my fit reproducible?"** — `saveRDS` round-trip plus
`revive()`, the *fit-as-record* pillar.
**12. "What did the optimizer actually do?"** — `changes()`,
`diagnostics()`, `fit_status()`, the *no-silent-surgery* pillar.

These are layered, not displaced, by a strong bootstrap story.

## What this implies for the package

**A. The "Why mixeff?" four-panel poster should make bootstrap the
headline inferential surface, not the fallback.** Currently Panel C
shows `asymptotic_wald_z` / `reliability = low` on the singular fit
because that is the only *immediately-available* labelled p-value. The
strong play is a bootstrap row with finite p-value, status = available,
and the run payload (`boundary_rate`, `mcse`, `successful_replicates`)
visible as the audit trail. Wald z stays as the fast default; bootstrap
becomes the certifiable answer.

**B. A new vignette: "Inference where lme4 says no."** Exhaustive walk
through variance-component LRT, simple-slopes p-values on a 3-way
interaction, fixed-effect CI on a singular fit, all via bootstrap, all
fast. This is the page that converts lme4 users.

**C. Speed becomes a benchmark we should track.** "How long does our
nsim = 1000 bootstrap take vs. `lme4::bootMer`, `pbkrtest::PBmodcomp`,
`RLRsim` on the same model?" A 10× table is the strong pitch. A 2×
table is "honest" not "fast" — fine, but materially different
positioning. The benchmark should live in `benchmarks/` and be run on
CI on a small fixture set so we can detect regressions.

## Prerequisites

The strategic argument above does not exist for users until the
ergonomic verbs exist. Specifically:

1. **`test_effect(fit, term, method = "bootstrap", nsim = ...)`** —
   single-df via `contrast()` dispatch, multi-df via `compare()`-routed
   bootstrap LRT (labelled `method = "bootstrap_lrt"` to distinguish
   from a future Wald-joint bootstrap when upstream B1 lands).
2. **`inference_options(fit, term)`** — capability/audit map, not a
   recommender. Tells the user *which* methods will succeed on this
   fit, with expected status/reliability/cost. Same audit-not-selection
   philosophy as `random_options()`.
3. **`confint(fit, parm, method = "bootstrap")`** —
   percentile/basic CIs on fixed effects. Larger than "take quantiles":
   needs failed-refit handling, reproducibility metadata, and
   successful-refit thresholds. Treat as a separate slice from items
   1 and 2.
4. **MC-SE display on bootstrap rows** — small change, but should be
   threaded through whichever bootstrap surfaces actually report
   p-values. Better deliberate than bolted on.

Items 1 and 2 are the headline ergonomic fix and are doable R-side
without upstream changes (the contrast bridge already certifies scalar
bootstrap rows; `compare(method = "bootstrap")` already returns
bootstrap LRT). Items 3 and 4 are separable and can be scheduled
deliberately.

## Out-of-scope claims (so we don't overpromise)

These are things bootstrap **does not** do, and the package's marketing
language should not imply they do:

- Make an unidentifiable variance component identifiable.
- Rescue an underpowered design.
- Substitute for scientific judgement about random-effects structure.
- Calibrate inference for a misspecified mean structure.
- Produce CIs on variance components near zero without a labelled
  interval-kind / boundary-mass story.
- Make GLMM bootstrap automatically defensible in the presence of
  separation, rare events, or overdispersion.

The pitch is **"bootstrap makes the defensible path ergonomic"** —
not "bootstrap solves mixed-model inference." The latter is overclaim;
the former is true and unique.

## Cross-references

- `planning/upstream-bootstrap-coverage-issue.md` — upstream contract
  asks (multi-df joint null, full-model-distribution intervals, cluster
  resampling, replicate statistics).
- `planning/bootstrap_fixed_effect_payload_surface.md` — downstream
  policy: every fixed-effect bootstrap p-value mixeff prints routes
  through the Rust contract, not R-side reconstruction.
- `planning/PRD.md` §3 — the "no model selection or random-effects
  recommendation engine" non-goal that bootstrap LRT must respect.
- `planning/vision.md` — the "refusal is a feature" / "every printed
  claim has a provenance" framings that bootstrap fits inside, not
  on top of.
- `planning/mission.md` — the "no fake certainty" principle that the
  bootstrap pitch must not violate by overstating "calibrated tail
  area."
