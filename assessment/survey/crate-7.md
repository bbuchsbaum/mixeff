# Crate Survey 7 — Engine Contract Docs: random-effects formulas, random-term card PRD, formula-transform seam

Survey date: 2026-05-31
Source docs:
- `/Users/bbuchsbaum/code/rust/mixeff-rs/docs/random_effects_formulas.md`
- `/Users/bbuchsbaum/code/rust/mixeff-rs/docs/random_term_card_prd.md`
- `/Users/bbuchsbaum/code/rust/mixeff-rs/docs/formula_transform_seam.md`

---

## 1. What the engine supports

### 1.1 Random-effects formula syntax (v0 contract)

| Construct | Engine status |
|---|---|
| `(1 + x | g)` — correlated random intercept + slope | Supported, conforms |
| `(1 + x || g)` — zero-correlation diagonal | Supported; `||` is a single token; zerocorr flag kept on `RandomTerm` IR |
| `(x | g)` — implicit intercept | Supported; intercept added at materialization |
| `(0 + x | g)` — suppressed intercept / cell-means | Supported; `InterceptPolicy::Omitted` first-class |
| `(1 | g1:g2)` — cell grouping | Supported |
| `(1 | g1 & g2)` — legacy interaction | Supported (verbatim `&` syntax preserved) |
| `(1 | a/b)` — nested grouping | Supported; expanded at parse time → `(1|a) + (1|a:b)` (R1) |
| `(1 | a*b)` — crossed grouping | Supported; expanded at parse time → `(1|a)+(1|b)+(1|a:b)` (R2); `CrossingLikelyUnintended` Info emitted |
| Numeric basis column | Conforms |
| Categorical basis column — treatment coding | Conforms; explicit contrast bases supported |
| Categorical basis column — cell-means `0 + factor` | Conforms |
| Interaction basis `x:y` in random term | Partial — treatment-coded and cell-means numeric/categorical products materialize; random-side empty-cell diagnostics still open |
| Duplicate term detection `(b|g) + (b|g)` | Diagnostic emitted (`DuplicateRandomTerm` Warning); canonical/effective merge pending |
| Conflicting covariance `(b|g) + (b||g)` | Error diagnostic emitted; fit refusal still pending |
| Fixed/random redundancy `g + (1|g)` | Diagnosed (`FixedRandomRedundant`); `design_compiled` drop/refit not yet applied |
| Multiple random terms, same grouping, different basis | Correctly preserved as independent blocks; no false merge |

### 1.2 Canonicalization rules (R1–R9, all v0-bound)

All nine canonicalization rules (R1 nesting, R2 crossing, R3 zerocorr-as-flag, R4 InterceptPolicy, R5 duplicate detection, R6 conflicting covariance, R7 same-grouping/different-basis preserved, R8 fixed/random redundancy, R9 source-syntax preservation) are declared and partially implemented. Full conformance status in Appendix A of the contract doc; main outstanding items:

- R4: `InterceptPolicy` enum exists; audit/print/source distinction needs hardening
- R5/R6: diagnostics emitted; canonical/effective rewriting pending
- R8: diagnosed; `design_compiled` automatic drop pending

### 1.3 Diagnostics inventory (complete engine registry)

All codes from the engine's `diagnostics.rs` are registered in the R wrapper's `mm_diagnostic_code_registry` (R/diagnostics.R lines 157–215). The five pedagogical taxonomy codes introduced by `random_term_card_prd.md` are implemented in the engine and registered in the R wrapper:

| Code | Bucket in R wrapper |
|---|---|
| `scope_note` | `design_note` |
| `support_note` | `design_note` |
| `syntax_expansion` | `design_note` |
| `covariance_assumption` | `design_note` |
| `structural_refusal` | `repair` |

Additional codes (all registered): `formula_canonicalized`, `formula_canonicalization_unsupported`, `duplicate_random_term`, `conflicting_covariance`, `crossing_likely_unintended`, `random_slope_without_intercept`, `fixed_random_redundant`, `repeated_unit_unmodeled`, `random_slope_unsupported`, `random_effect_few_levels`, `covariance_too_rich`, `covariance_reduced`, `fixed_effect_column_missing`, `fixed_effect_rank_deficient`, `fixed_effect_empty_cell`, `boundary_parameter`, `near_unit_random_effect_correlation`, `binomial_separation`, `not_identifiable`, `invalid_agq_request`, `optimizer_not_assessed`, `optimizer_nonconvergence`, `optimizer_recovery`, `pirls_failure`, `inference_unavailable`, `serialization_not_assessed`, `unsupported`.

### 1.4 RandomTermCard schema (random_term_card_prd.md)

The engine now exposes `RandomTermCard` structs on `ModelAuditReport.random_term_cards` (schema `mixedmodels.random_term_card` v1). Each card carries: `term_id`, `original_fragment`, `canonical_fragment`, `group`, `blocks` (with per-block `english` sentence authored upstream), `implied_constraints`, `design_support` (group_levels, min/median rows per group, within_group_variation map, InformationBudgetStatus), `role_origin` (declared_by_user, observed_from_data, resolved GroupingRole).

Block decomposition decision (Option B locked): `(1 + x || g)` produces two `RandomTermIr` entries with shared `block_group`; cross-card zero-covariance constraints live on `ModelAuditReport.cross_card_constraints`. `||` and split-block forms produce structurally identical card lists differing only in `original_fragment` and constraint `reason` string.

`MODEL_AUDIT_REPORT_SCHEMA_VERSION` bumped to v2 to accommodate the new field.

The R wrapper surfaces cards through `audit_design()`, `explain_model()`, `random_options()`, `compare_covariance()`, and the reporting layer — all consuming `audit$random_term_cards`.

### 1.5 Formula-transform seam (formula_transform_seam.md)

The seam is decided and frozen. Two classes:

**Engine-owned (stateless/pointwise) — evaluated below the seam:**
- `I(<arith>)` where `<arith>` covers `+`, `-`, `*`, `/`, `^`, unary `-`, parentheses, numeric literals, column references
- Bare pointwise calls: `log`, `log2`, `log10`, `exp`, `sqrt`, `abs` (composable inside `I()`)
- Applied to both sides: `log(y) ~ x + I(x^2) + (1 | g)` is valid
- `predict_new` re-evaluates the same stateless expression on newdata; no stored recipe

**Host-wrapper-owned (stateful/basis) — forbidden in engine, must be pre-evaluated by R:**
- `poly(x, d)`, `scale(x)`, `ns(x, df=)`, `bs(x)`, `cut(x, ...)`, `factor(x)`, `center(x)`, any function with a second argument (e.g. `log(x, base)`)
- Engine emits an actionable refusal naming the construct if these appear in the formula string

**R wrapper obligation** (per the seam contract, currently NOT implemented in R/compile.R or R/data-translate.R):
- Must own the model frame for stateful transforms using R's `model.matrix`/`terms`/`predvars`
- Must reuse training-time bases on `newdata` (predict path)
- Currently `compile_model()` uses `all.vars(formula)` to narrow data and passes raw formula string directly to Rust — no R-side model frame expansion

### 1.6 Feature-gated capabilities

| Feature | Gate | Status |
|---|---|---|
| `joint_laplace` GLMM estimation | `nlopt` Rust feature | Disabled in CRAN/vendored build; raises `mm_fit_error: estimation_method_unavailable` with actionable message |
| PIRLS-profiled GLMM | default (no gate) | Available |
| AGQ metadata (nAGQ > 1) | default | Available (profiled path) |

### 1.7 Model modes (inference contract taxonomy)

Four modes are defined: `as_specified`, `design_compiled` (default), `exploratory`, `predictive`. Mode affects how redundancy (R8) and similar rules behave at design-audit time. The R wrapper does not yet expose a user-facing mode selector; `design_compiled` is the effective mode.

### 1.8 Grouping factor materialization — non-conforming items

Two items remain non-conforming in the engine implementation (tracked in Appendix A of the formula contract):

1. **Composite-level separator**: implementation still joins with `_` instead of the collision-free `\x1E`-escaped key specified in §5.2. Display is correct but internal key collisions are possible on adversarial level labels.
2. **Lexicographic level ordering**: v0 contract requires lexicographic ordering for composite grouping keys; implementation still uses first-appearance order. This leaks into `parmap`, parameter ordering, and output tables. R wrapper is downstream of this non-conformance.

---

## 2. Documented refusals

| Construct | Engine behavior |
|---|---|
| Empty grouping `(x | )` | `FormulaError::EmptyGrouping` |
| Missing `|` in RE block | `FormulaError::MissingBar` |
| Unmatched parens / unknown tokens | `FormulaError::UnmatchedParen` / `FormulaError::UnexpectedToken` |
| `(b|g) + (b||g)` same basis/grouping | `ConflictingCovariance` Error (fit refusal pending) |
| `g + (1|g)` in `as_specified` mode | `NotIdentifiable` Error (refused) |
| Grouping variable missing from data | `MixedModelError::InvalidArgument` |
| Stateful transforms (`poly`, `scale`, `ns`, etc.) in formula | Actionable refusal naming the construct, pointing to host-wrapper precompute |
| `joint_laplace` without nlopt feature | `mm_fit_error: estimation_method_unavailable` |
| Smooth terms `s()`, `te()`, `bs()`, `ns()` | Not parsed in v0; rejected as unknown token |
| Offsets, `I(...)` beyond stateless subset | Refusal with actionable message |
| GAM smooths | Out of scope, not v0 |
| Multivariate `cbind(y1,y2)~` | Out of scope; tracked in `multivariate_shared_theta.md` |
| AR(1)/spatial residual structure syntax | Out of scope (inference contract, vNext) |

---

## 3. Stability labels and engine-side outstanding items

Items marked non-conforming in the engine contract that affect R wrapper outputs:

| Area | Non-conformance | Downstream R impact |
|---|---|---|
| `||` centering for numeric slopes | Not implemented (§4.5); declared/weighted-mean reference rule decided but not coded | `(1 + x || g)` with numeric x is misspecified without centering; independence assumption is at implicit zero, not mean(x). No R-side mitigation either. |
| Composite-level separator | Uses `_` not `\x1E` | Potential label collisions if grouping levels contain `_` |
| Lexicographic level ordering | Uses first-appearance | Parameter ordering is run-order-dependent; fixture drift on row shuffle |
| Random-side empty-cell diagnostics | `FixedEffectEmptyCell` reuse for random side missing | Unbalanced interaction groupings not explicitly diagnosed |
| Canonical/effective formula rewriting | Design-audit reductions diagnosed but not applied at fit | `design_compiled` mode drop of `FixedRandomRedundant` and `RandomSlopeUnsupported` terms not automatic |
| Interaction basis empty-cell diagnostics | Open | Unbalanced factor×factor random slopes may fail silently |

---

## 4. Engine capabilities not yet surfaced in the R wrapper

| Engine capability | R wrapper gap | Classification |
|---|---|---|
| Stateless formula transforms (`I(x^2)`, `log(y)`, etc.) | `compile_model()` passes raw formula string; R-side model frame expansion for stateless transforms is not separately documented. The engine handles these natively, but the R wrapper uses `all.vars()` which extracts `x` and `y` from `I(x^2)` — the R wrapper would strip the transform before Rust sees it, or Rust would see a column named "I(x^2)" that doesn't exist in `data`. This is an untested gap. | partial |
| Stateful transform pre-evaluation (`poly`, `scale`, `ns`) | `compile_model()` calls `all.vars()` to narrow columns but does no `model.matrix()` expansion; stateful transforms in formula strings will hit Rust's refusal. The R wrapper does not catch these upstream with a friendly error, nor does it pre-evaluate them. `predvars` are not stored. | partial |
| Host-wrapper `predvars` contract for `predict(newdata)` | `mm_predict_conditional_newdata()` passes raw column data; no stateful transform basis is stored or replayed | partial |
| `RoleOrigin.declared_by_user = true` path | `roles()` string-form not wired through FFI in Phase 1.F; `declared_by_user` is always `false` in v1 cards | partial |
| Mode selector (`as_specified`, `exploratory`, `predictive`) | No user-facing API; always `design_compiled` | partial |
| `||` centering reference (declared or weighted-mean) | Not implemented upstream or in wrapper | upstream-blocked |
| `cross_card_constraints` on `ModelAuditReport` | R wrapper reads `audit$random_term_cards`; does not yet read a `cross_card_constraints` field if present | partial |
| `median_obs_per_level` on `GroupingAudit` | PRD says to add it; R card consumer uses `design_support.median_rows_per_group` | upstream-blocked |
