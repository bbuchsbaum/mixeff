# `mixeff` LMM Reporting Contract

Status: draft  
Tracking bead: `bd-01KQDT0Y71B13TVQ9ARSTRG2CY`  
Related upstream field checklist: `planning/reporting_artifact_requirements.md`  
Motivation: Davies & Meteyard-style best-practice reporting for LMMs

## Purpose

`mixeff` should make complete mixed-model reporting a normal part of the fitted
object workflow. A user who has fit an LMM should be able to produce the core
materials needed for a methods/results report without reconstructing formulas,
random-effect choices, inference methods, convergence state, or software
versions by hand.

The initial user-facing surface is:

```r
model_report(fit)
reporting_table(fit, section = "all")
```

`model_report()` returns a structured S3 report object suitable for printing,
inspection, and future export. `reporting_table()` returns data-frame-compatible
tables for individual sections.

The contract is deliberately reporting-oriented, not recommendation-oriented:
the report says what was requested, what was fit, what changed, what was
available, and what was unavailable with reasons.

## Design Rules

1. Rust-owned claims stay Rust-owned. `mixeff` formats model semantics,
   diagnostics, inference rows, fit status, covariance state, and reason codes
   from versioned Rust artifacts.
2. R-owned provenance stays R-owned. `mixeff` records the call, raw data row
   count, row-exclusion policy, R version, package versions, locale, contrasts,
   and row names.
3. Every report row has provenance. Each table carries `source`, `schema_name`,
   `schema_version`, or an R-owned provenance marker.
4. Unavailable is a reportable status. Missing p-values, unsupported intervals,
   invalid comparisons, and unavailable diagnostics appear as rows with reasons,
   not as silent omissions.
5. No advice creep. The report must not tell users to drop slopes, try a simpler
   model, or prefer a recommended random-effects structure.

## User-Facing Objects

### `model_report(fit, sections = "all")`

Returns class:

```r
c("mm_model_report", "mm_report")
```

Required fields:

- `metadata`: report creation time, package version, crate version, schema
  versions, fit class, and section availability.
- `sections`: named list of section tables.
- `unavailable`: table of sections or fields that could not be produced, with
  reasons.
- `provenance`: compact source map from report sections to artifact paths.

The default print method is compact and points to section tables rather than
dumping all output at once.

### `reporting_table(fit, section = "fixed_effects")`

Returns one data-frame-compatible object for a named section. Initial section
names:

- `overview`
- `model_specification`
- `data_design`
- `random_terms`
- `random_effects`
- `fixed_effects`
- `fit_statistics`
- `optimizer`
- `comparison_ledger`
- `reproducibility`
- `unavailable`

`section = "all"` returns a named list.

## Required Report Sections

### 1. Overview

Purpose: one compact table users can paste into an analysis supplement.

Columns:

- `field`
- `value`
- `source`
- `status`
- `reason`

Required rows:

- model class: LMM or GLMM
- formula as written
- effective formula, if different
- fit method: REML, ML, Laplace, AGQ, etc.
- public mode and internal fit intent
- number of observations used
- fit status
- inference availability summary
- schema/crate/package versions

### 2. Model Specification

Purpose: make the exact fitted model auditable.

Rows:

- original formula
- canonical formula
- requested random structure
- effective random structure
- fixed effects and contrasts policy
- family/link for GLMMs
- weights/offset status where supported
- model changes, each with category and reason

Source:

- Rust artifact and model-state summary for model semantics.
- R fit object for call capture, contrasts policy, and R-side formula context.

### 3. Data and Design

Purpose: report sample size in the units that matter for mixed models.

Rows or columns:

- raw rows supplied to `mixeff`
- rows used after R filtering / `na.action`
- rows received by Rust
- grouping factor
- grouping role
- number of levels
- min rows per group
- median rows per group
- max rows per group
- within-group variation status for relevant slopes
- support status and thresholds

Source:

- R row map for raw/excluded rows.
- Rust audit report and random term cards for grouping facts and support.

### 4. Random Terms

Purpose: explain the random-effects syntax without model advice.

Rows:

- original formula fragment
- canonical fragment
- group
- block id
- basis
- intercept included
- slopes
- covariance family
- theta parameter count
- Rust-authored English
- implied zero-covariance constraints
- design-support status

Source:

- `random_term_cards`
- `cross_card_constraints`

Rendering rules:

- Render Rust-authored `english` fields verbatim.
- Render constraints as assumptions, not warnings.
- For intercept-only terms with possible unmodeled slopes, report the
  `ScopeNote` fact if Rust supplies it; do not escalate it.

### 5. Random Effects

Purpose: satisfy full random-effect variance/covariance reporting.

Rows:

- group
- term/block id
- basis lhs
- basis rhs
- variance or covariance kind
- estimate
- standard deviation where applicable
- correlation where modeled
- covariance family
- boundary/reduced-rank/fixed-zero status
- reason

Source:

- Preferred: stable Rust random-effect variance/covariance payload described
  in `planning/reporting_artifact_requirements.md`.
- Temporary fallback: `VarCorr()` output only where it can be generated from
  Rust-owned covariance summaries without reinterpreting theta semantics in R.

### 6. Fixed Effects

Purpose: report estimates and inference with method labels.

Columns:

- `label`
- `estimate`
- `std_error`
- `df`
- `statistic`
- `statistic_name`
- `p_value`
- `method`
- `status`
- `reliability`
- `estimability`
- `reason`
- `notes`

Source:

- Rust `fixed_effect_inference_table`.

Rules:

- Do not compute p-values in R.
- Do not substitute Wald, Satterthwaite, Kenward-Roger, or bootstrap labels in R.
- Print `NA` with Rust reason when a row is unavailable.

### 7. Fit Statistics

Rows:

- log likelihood
- deviance or REML criterion
- AIC
- BIC
- number of observations
- residual scale / dispersion
- number of estimated parameters, if Rust supplies it
- GLMM approximation metadata where applicable

Source:

- Rust fitted artifact/model state.

### 8. Optimizer and Fit Status

Rows:

- fit status
- optimizer name
- objective value
- iteration count
- stop reason
- convergence verification summary
- active boundary parameters
- reduced-rank summary
- certificate quality/status

Source:

- Rust optimizer certificate and audit report.

Rules:

- Boundary and singular states are reported as fitted model state.
- The report must not print repair advice such as "drop the random slope."

### 9. Comparison Ledger

Purpose: make model-building history reportable when users compare models.

Rows:

- comparison id
- model id
- formula
- fit method
- whether ML refit occurred
- comparison method
- statistic
- df
- p-value or information criterion
- validity status
- refusal/unavailable reason

Source:

- R ledger for user-call sequence and object ids.
- Rust comparison validity/result payload for statistical status and reasons.

This section is empty with `status = "not_applicable"` when no comparison has
been requested.

### 10. Reproducibility

Rows:

- `mixeff` version
- Rust crate version
- schema versions
- compiler policy and thresholds
- optimizer controls
- seed/random-state status
- R version
- platform
- relevant package versions
- contrast options
- locale

Source:

- Rust reproducibility record.
- R session provenance.

### 11. Unavailable and Caveated Outputs

Every unavailable section or row is collected here.

Columns:

- `section`
- `field`
- `status`
- `reason`
- `source`
- `action_taken`

Examples:

- fixed-effect row p-value unavailable because derivative prerequisites failed
- Kenward-Roger unavailable because method was not requested or not supported
- random-effect covariance table unavailable because upstream schema is too old
- comparison refused because REML fixed-effect comparison is invalid

## Snapshot Test Requirements

Reporting tests should assert:

- `model_report()` includes overview, model specification, data/design,
  random terms, fixed effects, fit statistics, optimizer, reproducibility, and
  unavailable sections.
- Fixed-effect rows preserve Rust `method`, `status`, `reliability`, and
  `reason`.
- Random-term rows preserve Rust-authored English and covariance constraints.
- Random-effect variance/covariance rows include all modeled variances and
  correlations/covariances when available.
- Grouping summaries include levels and min/median/max rows per group.
- Boundary/reduced-rank states appear in optimizer/random-effect sections.
- Report output contains no forbidden recommendation language:
  "we recommend", "you should", "try ... instead", "drop the random slope",
  "suggested starting model".
- Saved and revived fits produce the same report sections from stored artifacts
  whenever no live Rust handle is required.

## Phase Placement

The reporting contract is a planning dependency for Phase 2/3 implementation,
not a blocker for basic Phase 1 fitting. A minimal `model_report()` can ship
after these prerequisites are true:

1. fitted artifacts expose fixed-effect inference rows;
2. audit reports expose random term cards and cross-card constraints;
3. optimizer certificates and reproducibility records parse in R;
4. `mixeff` stores R-side row-exclusion and session provenance;
5. either Rust exposes a stable random-effect variance/covariance payload or
   `mixeff` explicitly marks that section unavailable.

The random-effect variance/covariance payload is the only major upstream gap
that should be resolved before treating `model_report()` as best-practice
complete.
