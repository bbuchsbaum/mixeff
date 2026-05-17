# Marginal Quantities Deferrals

Status: planning note  
Tracking bead: `bd-01KQFZ3Q4615SS1HYVJBQYGDFJ`  
Related schema: `inst/schemas/mixedmodels.marginal_quantity_table.schema.json`  

## Purpose

The marginal-quantity table schema is intentionally broader than the first
native `mixeff` verbs. Phase 3.x should define contract-shaped rows for
marginal quantities without implying that every quantity in the schema is
implemented, certified, or ready for user-facing helpers.

This note records the deferred surface so Phase 3.x does not grow into an
uncertified replacement for `emmeans` or `marginaleffects`.

## Phase 3.x Boundary

Phase 3.x may expose contract-shaped rows for LMM marginal predictions, means,
and simple comparisons only when the required uncertainty, target, scale,
weights, and status fields can be populated honestly.

Unsupported but well-formed requests should return rows shaped like
`mixedmodels.marginal_quantity_table` with:

- `status = "not_assessed"` or another stable unavailable status;
- explicit `reason`;
- explicit `target`, `scale`, `weights`, and `comparison`;
- `NULL`/`NA` estimates or inference fields where Rust has not certified the
  quantity.

Verb names are not the contract. The contract is the table shape, row status,
and reason vocabulary.

## Deferred From Phase 3.x

Do not ship native `mixeff` implementations of these in Phase 3.x:

- `mm_slopes()` or other numeric-derivative slope helpers;
- GLMM response-scale marginal means;
- ratios, odds ratios, or risk ratios as user-facing native quantities;
- delta-method or simulation-backed response-scale uncertainty unless Rust
  certifies the method and records the approximation status;
- transformed-predictor marginal effects that require derivative bookkeeping;
- interaction simple slopes that require method-specific uncertainty beyond the
  fixed-effect contrast table;
- link-to-response transformations whose interval or p-value interpretation is
  not represented in the row status/reason fields.

The schema may still include `comparison = "ratio"` and
`comparison = "odds_ratio"` so unsupported requests can be represented as
contract-shaped rows. Their presence in the schema is not implementation
approval.

## Rationale

Numeric derivatives interact with transformed predictors, interactions,
contrasts, link functions, and conditioning targets. A slope value without a
certified uncertainty method is easy to print and hard to defend.

Ratios and odds ratios require scale-specific uncertainty and interpretation.
For GLMMs, response-scale quantities add at least three status axes: link versus
response scale, conditional versus population target, and approximation method.

`mixeff` should keep its audit-first rule: available numbers are labelled with
method and reliability; unavailable numbers are represented with stable
reasons.

## Expected Destination

Move this surface to Phase 4/5 or later:

- after the GLMM fit bridge has stable link/response semantics;
- alongside or after `recover_data.mm_fit()` and `emm_basis.mm_fit()` support;
- after Rust can certify the uncertainty method for response-scale and ratio
  quantities;
- after model-report rendering has examples for `not_assessed` marginal rows.

For breadth, users should continue to use `emmeans` or `marginaleffects` when
they need a wide marginal-effects interface. Native `mixeff` rows remain
contract-first and should prefer explicit unavailable rows over approximate
numbers without provenance.

## Acceptance Checks

Before any deferred quantity becomes user-facing, require:

- a schema row with explicit `target`, `scale`, `weights`, `comparison`,
  `method`, `requested_method`, `status`, `reliability`, `reason`, and
  `details`;
- tests for a supported row and a well-formed unsupported row;
- model-report rendering that preserves the unavailable reason;
- documentation that distinguishes schema vocabulary from implemented verbs;
- no advice language suggesting a preferred marginal-effects workflow.
