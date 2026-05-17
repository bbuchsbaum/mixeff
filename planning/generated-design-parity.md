# Generated-Design Parity Plan

Status: planning note  
Tracking bead: `bd-01KQDEZQXRQD8PDV5FYVREVEMX`  
Pilot generator: `tools/generated-design-pilot.R`  
Parent context: classic `lme4` fixture parity in `tests/fixtures/lme4_baseline_cases.json`

## Purpose

Classic `lme4` datasets are the first parity baseline because they are familiar,
small, and easy to inspect. Generated designs should come later, after the
classic baseline is stable, to cover formula and data regimes that the classic
fixtures do not exercise systematically.

The generated suite should be a contract test layer, not a fuzzing dump. Every
case must state which design axes it covers, which fields are compared, which
tolerances apply, and whether the expected outcome is `match`,
`known_fragile`, or `expected_unavailable`.

## Scope

The first generated-design suite should cover LMM formulas only. It should
compare `mixeff::lmm()` with `lme4::lmer()` for extractor and reporting fields
already covered by the classic parity helpers:

- fixed effects;
- residual scale and likelihood summaries;
- fitted values and residuals;
- `VarCorr()` summaries;
- random effects when labels align;
- supported fitted-data prediction semantics;
- fit status and singular/boundary classification where appropriate.

It should not become a random search over formulas in PR CI. Larger randomized
or stress sweeps belong in an opt-in or scheduled job after the pilot cases are
stable.

## Design Axes

Generated cases should be selected from an explicit matrix:

| Axis | Values to cover | Primary risk |
| --- | --- | --- |
| Fixed-effect rank | full rank, aliased column, empty factor cell | column dropping and coefficient naming |
| Predictors | numeric, categorical, interaction, transformed predictor | formula/model-matrix translation |
| Group balance | balanced, unbalanced, sparse groups | row grouping and support summaries |
| Random structure | intercept, slope, `||`, split block | theta/VarCorr mapping |
| Grouping topology | single, nested, crossed | grouping labels and row order |
| Fit mode | REML, ML | likelihood and comparison semantics |
| Boundary regime | interior, near-boundary, singular | optimizer/status classification |
| Scale | centered/scaled, large magnitude predictor | numeric conditioning |
| Hypothesis rhs | zero, nonzero scalar contrast | inference row bookkeeping |

Each generated case should cover a small number of axes. A case that tries to
cover everything is hard to diagnose when it fails.

## Classification Contract

Every generated case must be classified before it joins CI:

- `match`: expected to match `lme4` within documented tolerances.
- `known_fragile`: valid case where optimizer, boundary, or rank behavior may
  differ; assertions should use a ledger bound or status-specific invariant.
- `expected_unavailable`: well-formed request that `mixeff` intentionally
  refuses or cannot certify yet.

Unclassified failures should fail loudly. Do not silently relax tolerances for
generated cases without a ledger entry explaining why.

## Pilot Generator

`tools/generated-design-pilot.R` defines `mm_generated_design_pilot()`, a small
sourceable generator that returns four deterministic cases:

- `gen_balanced_random_intercept`: balanced groups, numeric predictor,
  categorical predictor, random intercept.
- `gen_unbalanced_random_slope`: unbalanced groups with within-group slope
  variation.
- `gen_crossed_random_intercepts`: crossed subject/item intercepts.
- `gen_boundary_random_slope`: random-slope variance generated as zero and
  classified as `known_fragile`.

The pilot intentionally returns data frames in memory rather than writing large
fixtures. A future test can materialize a manifest-like structure from this
function and feed it through parity helpers adapted from
`tests/testthat/helper-lme4-parity.R`.

## Candidate Test Families

### Differential Checks

For `match` cases, compare against `lme4::lmer()` with the existing tolerance
model. Start with fixed effects, sigma, fitted values, residuals, and
likelihood. Add random-effect and `VarCorr()` checks once labels and scale
normalization are stable.

### Metamorphic Checks

Use generated cases to test relationships that do not require an exact oracle:

- row permutation should preserve estimates after reordering fitted outputs
  back to original row ids;
- equivalent random-effect spellings such as split block versus `||` should
  preserve intended zero-covariance structure;
- centering a numeric predictor should preserve fitted values after translating
  the intercept;
- balanced and unbalanced versions of the same data-generating process should
  remain finite and converged unless intentionally boundary-classified.

### Edge Checks

Keep edge cases small and labelled:

- aliased fixed-effect column: classify as `known_fragile` until column-dropping
  parity is explicit;
- categorical empty cell: classify as `known_fragile` or
  `expected_unavailable`, depending on the current model-matrix contract;
- zero random-slope variance: classify as `known_fragile` and check status
  rather than exact component equality;
- too-few levels for a rich random slope: classify as `expected_unavailable`.

## Proposed Integration Steps

1. Keep `tools/generated-design-pilot.R` as an untested planning/pilot helper
   until classic fixture parity and the current mismatch ledger are stable.
2. Add `tests/testthat/helper-generated-design-parity.R` that converts pilot
   cases into the same shape as `mm_lme4_parity_cases()`.
3. Add one smoke test that validates generator structure without fitting.
4. Add one opt-in or skipped-by-default differential test for the three
   `match` pilot cases.
5. Add one `known_fragile` boundary test that asserts finite outputs plus an
   explicit singular/boundary status classification, not exact random-slope
   variance equality.
6. Feed generated-case rows into the parity scoreboard so case id, seed, axes,
   field, tolerance, observed difference, status, and reason are recorded.
7. Only after those pass reliably, add rank-deficient and empty-cell generated
   cases.

## CI Policy

Generated parity should follow a test pyramid:

- PR CI: structure checks, one or two fast `match` cases, and no randomized
  search.
- Scheduled/opt-in: larger axis matrix and stress cases.
- Manual upstream filing: any case classified as upstream bug must include the
  generator seed, formula, data dimensions, reference versions, observed diffs,
  and minimal reproducer.

## Acceptance Criteria

Generated-design parity is ready to graduate from planning when:

- every generated case has `id`, `formula`, `reml`, `seed`, `axes`,
  `expected_status`, and `notes`;
- every compared field has a tolerance or ledger classification;
- `known_fragile` and `expected_unavailable` cases assert status/reason
  contracts instead of passing by omission;
- the scoreboard records generated-case seed and axes;
- no generated test adds substantial runtime to ordinary PR checks.
