# Reporting Artifact Requirements for `mixedmodels`

Status: draft for cross-repo discussion  
Tracking bead: `bd-01KQDT19VGXK56XRQ3M30AGXR8`  
Audience: `mixeff` and `mixedmodels` engineers  
Motivation: Davies & Meteyard-style best-practice reporting for LMMs

## Purpose

`mixeff` should make transparent LMM reporting easy enough that users do not
need to reconstruct the analysis history from memory, console output, or warning
strings. The R package can assemble tables and prose, but the statistical facts
must come from `mixedmodels` artifacts.

This document defines the report-required fields that `mixeff` needs from Rust
to implement a future `model_report()` / `reporting_table()` surface. It is a
discussion checklist for `mixedmodels` engineers, not a request for R-specific
formatting upstream.

The guiding rule is:

> `mixeff` formats. `mixedmodels` authors model semantics, diagnostics,
> inference availability, fit status, and reasons.

## Non-Goals

- Do not ask Rust to generate publication prose or R-specific tables.
- Do not add a model-selection recommendation engine.
- Do not require Rust to know R-side data-cleaning details that happen before
  the model frame is frozen.
- Do not require p-values for every model. Report availability, method, and
  reason row by row.
- Do not block Phase 1 on later diagnostic plots or residual-normality helpers.

## Report Sections `mixeff` Wants to Produce

The user-facing report should have these sections:

1. Analysis and software provenance
2. Model specification
3. Data and design summary
4. Random-effect specification and design support
5. Fit status and optimizer certificate summary
6. Fixed-effect estimates and inference rows
7. Random-effect variance/covariance table
8. Model fit statistics
9. Model-building or comparison ledger, when applicable
10. Reproducibility record
11. Unavailable or caveated outputs with reasons

The rest of this document maps those sections to required upstream artifacts.

## Field Matrix

| Report need | Required fields | Preferred Rust source | Current status / ask |
| --- | --- | --- | --- |
| Software provenance | crate name/version, schema names/versions, compiler contract version, optimizer backend, feature flags relevant to fit/inference | `CompiledModelArtifact`, `ModelAuditReport`, reproducibility record | Mostly present. Confirm every top-level report source carries stable schema metadata and crate version. Add optimizer backend/feature info if not already serialized. |
| R/package provenance | R version, `mixeff` version, package versions, locale, contrasts, call, row names | R-owned fit object | No Rust ask. `mixeff` records this and includes it in report assembly. |
| Formula and fit mode | original formula, canonical/effective formula, public mode, internal fit intent, REML/ML, family/link for GLMMs | `CompiledModelArtifact`, model-state summary | Present in principle. Confirm requested, semantic, supported, and fitted/effective formulas are all queryable after fit. |
| Model changes | requested -> semantic -> supported -> fitted transitions, reduction category, trigger, affected term, inference consequence, reason code | model-state summary, `changes()` payload, diagnostics | Present in design. Reporting needs this as structured rows, not display text. |
| Model frame size | rows used after subset/NA handling, response count, weight count if applicable | artifact/data manifest plus R row map | Split ownership. Rust should report rows actually received/used. R reports raw rows, excluded rows, and exclusion policy. |
| Grouping-unit counts | grouping factor name, role, number of levels, min/median/max rows per level, repeated-unit status | `ModelAuditReport`, `random_term_cards[].design_support`, grouping audit | Cards currently expose levels/min/median. Add or confirm max rows per group somewhere reportable. Report needs all three for unbalanced designs. |
| Within-group variation | for each fixed-effect basis candidate by group: present/absent/not assessed, threshold, role origin | `random_term_cards[].design_support`, diagnostics | Present for requested terms. Discuss whether report should also expose unrequested-but-possible slopes via `ScopeNote` payloads. |
| Random-effect meaning | original fragment, canonical fragment, group, blocks, basis, intercept/slopes, covariance family, theta parameter count, Rust-authored English | `random_term_cards` | Present. `mixeff` should render block `english` verbatim. |
| Fixed zero covariance assumptions | basis pair, reason, syntax source (`||`, split blocks), affected cards/terms | `random_term_cards[].implied_constraints`, `cross_card_constraints` | Present. Confirm constraints are emitted for both `||` and split-block forms. |
| Design-support status | parameter counts, levels per parameter, information-budget status, thresholds used | audit report, random-effect information budget, compiler policy | Mostly present. Reporting needs the threshold values alongside status so "low support" is auditable. |
| Fixed-effect table | label, estimate, standard error, df, statistic, p-value, method, status, reliability, estimability, reason, notes | `fixed_effect_inference_table` nested in fitted artifact | Present upstream. `mixeff` should consume exactly and not compute inference in R. |
| Term or contrast tests | same row schema as fixed-effect table, plus contrast matrix/term label, rhs, method request | future contrast/test endpoint returning `FixedEffectInferenceTable` | Needed after coefficient-row reporting. Rust should evaluate estimability, method prerequisites, df, statistic, p-value, status, and reason. |
| Random-effect variance table | final effective model random-effect SDs, variances, covariances/correlations if modeled, group labels, term/block ids, boundary/reduced-rank tags | preferred: stable Rust VarCorr/report table; fallback: artifact + theta map + fitted covariance summaries | Key discussion item. `mixeff` can format `VarCorr`, but report should not require R to reverse-engineer covariance semantics from theta. Add or confirm a stable random-effect variance/covariance payload. |
| Residual scale | sigma/residual variance, family-specific dispersion where applicable, residual df if meaningful | fitted artifact/model state | Needed for report and power-relevant variance decomposition. Confirm stable field names. |
| Fit statistics | logLik, deviance/REML criterion, AIC, BIC, nobs, df parameters, objective approximation for GLMMs | fitted artifact/model state | Present or computable. Prefer Rust-owned values with method labels. |
| Optimizer certificate | fit status, optimizer name, objective, iteration count, stop reason, convergence checks, boundary active set, reduced-rank evidence, verification record | `OptimizerCertificate`, audit report | Present. Reporting needs compact summary fields plus links/ids to full certificate details. |
| Inference caveats | unavailable method/status/reason by row, boundary-aware availability, derivative prerequisite status, GLMM unsupported finite-sample reason | inference rows, certificate, diagnostics | Present in design. Avoid global bans in R; Rust row status is authoritative. |
| Power-relevant inputs | group counts, random-effect variances, residual variance, ICC-like components when well-defined | random-effect variance table, residual scale, design summary | Derived report can be R-side if Rust exposes stable variance components. Do not require Rust to perform power analysis. |
| Model comparison ledger | model ids, formulas, fit method, ML/REML refit status, comparison method, statistic, df, p-value/AIC/BIC, validity/refusal reason | `mixeff` ledger plus Rust comparison validity/result API | Split ownership. R records sequence and user calls; Rust should validate comparison semantics and return method/status/reason rows. |
| Reproducibility | seed/random-state status, compiler policy, thresholds, optimizer controls, schema versions, crate version, data hash optional | reproducibility record | Present in design. Discuss whether a model-frame checksum/hash is useful and cheap enough. |
| Assumption diagnostics | residual vectors/fitted values, random effects, leverage-like diagnostics later | extractor APIs / future diagnostics | Later phase. Do not block reporting contract, but ensure extractors preserve enough identifiers for diagnostic plots. |

## Immediate Upstream Discussion Items

### 1. Stable Random-Effect Variance/Covariance Payload

The Davies-Meteyard guidance repeatedly emphasizes reporting random-effect
variances in full. `mixeff` should not have to reconstruct these from optimizer
theta slots.

Ask for a stable payload with one row per reportable variance/covariance item:

- `group`
- `term_id`
- `block_id` or block index
- `basis_lhs`
- `basis_rhs`
- `kind`: `variance`, `std_dev`, `covariance`, `correlation`, or equivalent
- `estimate`
- `scale`: user scale vs optimizer scale
- `covariance_family`
- `status`: interior, boundary, reduced_rank, fixed_zero, not_estimated
- `reason`
- links to `theta_map` slots where relevant

This can be a new table or a clearly documented part of an existing artifact.
The important property is that R can format the final effective `VarCorr`
without re-deriving covariance semantics.

### 2. Complete Design Summary for Reporting

Random term cards already carry much of this. Confirm or add:

- `n_obs_used`
- grouping levels per group
- min/median/max rows per group
- rows excluded by Rust, if any
- thresholds used for support decisions
- whether role facts were user-declared or observed from data

R will add raw input row count, subset/NA policy, and row-exclusion provenance.

### 3. Comparison Validity and Result Rows

The report needs to explain model-building when users compare models. `mixeff`
can store the ledger, but Rust should own the statistical validity of a
comparison:

- Are models nested when LRT is requested?
- Did fixed effects differ under REML?
- Was an ML refit required/performed?
- Did random effects differ while fixed effects also changed?
- Is the requested p-value/status available?

This should return structured status/reason fields, not warning text.

### 4. Reproducibility Record Completeness

Confirm the reproducibility record includes enough information for a report
footer:

- crate version and schema versions
- compiler policy and threshold values
- optimizer controls and backend
- seed/random-state status for simulation/bootstrap paths
- fit mode/intent
- relevant Cargo feature flags or backend capability flags

`mixeff` will add R session details.

### 5. Artifact Location and Schema Negotiation

For each report-required payload, agree on:

- schema name
- schema version
- whether it is nested in `CompiledModelArtifact`, `ModelAuditReport`, or
  returned from a table endpoint
- whether it is available on compile-only specs, fitted models, or both
- what unavailable payloads look like

This is more important than the exact nesting choice. The R wrapper can adapt
to nesting, but it needs stable schema negotiation.

## Proposed `mixeff` Consumption Rules

1. If Rust supplies a report field, `mixeff` formats it and preserves the
   Rust method/status/reason.
2. If a Rust field is absent because the schema is older, `mixeff` prints an
   explicit unavailable reason rather than silently computing a substitute.
3. R may compute display-only reshaping, percentages, and column ordering.
4. R may compute values that are inherently R-owned, such as raw input row
   counts, excluded rows from `na.action`, and package/session provenance.
5. R must not compute Satterthwaite, Kenward-Roger, bootstrap p-values,
   covariance reductions, singularity status, or random-term meanings.

## Acceptance Criteria for the Cross-Repo Contract

- A fitted LMM artifact plus audit report contains enough structured data for
  `mixeff` to render a complete reporting table without parsing warning text.
- Fixed-effect inference rows expose method/status/reason/reliability per row.
- Random-effect variance/covariance output can be rendered from a stable
  payload without decoding theta layout in R.
- Random-term explanations and covariance assumptions come from
  `random_term_cards` and `cross_card_constraints`.
- Optimizer and fit status are backed by certificate fields, not raw optimizer
  messages.
- Reproducibility fields include schema/crate/policy/threshold/optimizer
  metadata.
- Model comparison output carries validity and unavailable reasons explicitly.

## Open Questions for `mixedmodels` Engineers

1. Should the random-effect variance/covariance report be a new dedicated
   table endpoint, part of `CompiledModelArtifact`, or part of
   `ModelAuditReport`?
2. Should `max_rows_per_group` be added to `RandomTermCard.design_support`, or
   should `mixeff` read it from a separate grouping-audit table?
3. Should model-frame checksums be included in reproducibility records, or is
   that too expensive / too R-specific for v1?
4. Should comparison validity return a `FixedEffectInferenceTable`-like row
   schema, or a separate `ModelComparisonTable` schema?
5. Which report-required payloads are guaranteed after `compile_model()` versus
   only after `lmm()` / `glmm()` fit?

## Suggested Implementation Order

1. Confirm the field matrix against current Rust artifacts.
2. Add the random-effect variance/covariance payload if no stable table exists.
3. Fill design-summary gaps: especially max rows per group and thresholds.
4. Lock schema names/versions for all report-required payloads.
5. Update `mixeff` parsers and snapshots to consume the fields.
6. Implement `model_report()` / `reporting_table()` in `mixeff`.
7. Add the Davies-Meteyard reporting vignette.
