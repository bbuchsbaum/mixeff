# Crate Survey: Engine Ceiling — What mixeff-rs Supports and Refuses

**Source docs:** `model_comparison_policy.md`, `multivariate_shared_theta.md`,
`optimizer_profiles.md`, `mixeff_upstream_support_report.md`,
`v1_0_release_roadmap.md` (plus supporting contracts consulted:
`glmm_support_contract.md`, `random_effects_formulas.md`)

**Survey date:** 2026-05-31

---

## 1. Supported Model Classes

### 1.1 Linear Mixed Models (LMM)

- Gaussian response with identity link.
- REML and ML criteria both supported; REML is the default fit mode.
- Single-level and multi-level random effects: scalar random intercepts,
  vector random intercepts/slopes, crossed, nested, composite grouping
  (`g1:g2`, `g1/g2`, `g1*g2`, `g1 & g2`).
- Zero-correlation (`||`) and full-covariance (`|`) blocks.
- Profile likelihood for sigma, scalar/vector theta, and ML beta; REML beta
  profile is explicitly refused.
- Parametric bootstrap LRT (`stats::parametric_bootstrap_lrt`, schema
  `mixedmodels.parametric_bootstrap_lrt` v1.0.0) — LMM only.
- Boundary-sensitive LRT (`BoundaryLikelihoodRatioTest`, 50:50 chi-square
  mixture) — LMM variance-component comparisons.

**Stability:** production / release-blocking parity class.

### 1.2 Generalized Linear Mixed Models (GLMM)

Certified families and links (v1.0 contract):

| Family | Links |
|---|---|
| Bernoulli | logit, probit, cloglog |
| Binomial (incl. grouped proportions with weights) | logit, probit, cloglog |
| Poisson | log, sqrt |
| Gamma | log |

**InverseGaussian and Gaussian-as-GLMM:** implemented but NOT certified;
explicitly excluded from the v1.0 SemVer-covered GLMM contract. Parametric
bootstrap explicitly refuses these families.

**Fitting modes:**
- `fast = true` (default): profiled fast-PIRLS approximation; mirrors
  MixedModels.jl `fast=true`. Faster but less accurate for inference on
  overdispersed or OLRE models.
- `fast = false` (NLopt-gated): joint Laplace (`n_agq <= 1`) or joint AGQ
  (`n_agq > 1`, single scalar RE term only). Without NLopt the engine refuses
  this mode explicitly rather than silently falling back.

**AGQ restriction:** `n_agq > 1` is accepted only for exactly one scalar
random-effects term. Vector RE and multiple RE terms trigger a hard
`invalid_agq_request` diagnostic before optimizer evaluations.

**REML, Satterthwaite, Kenward-Roger:** LMM-only in the current contract.
GLMM artifacts must report finite-sample LMM inference as unsupported.

**Stability:** certified families above are production; fast-PIRLS rows that
match MixedModels.jl fast=true are `documented_divergence` (not lme4 parity).

---

## 2. Inference Methods

### 2.1 Fixed-Effect Inference

The engine owns inference rows via `FixedEffectInferenceTable`
(schema `mixedmodels.fixed_effect_inference_table` v1.0.0). Each row carries:
`label`, `kind` (coefficient/contrast/term), `estimate`, `std_error`,
`numerator_df`, `denominator_df`, `statistic`, `statistic_name` (z/t/f/chi_sq),
`p_value`, `method`, `status`, `reliability`, `estimability`, `reason`, `notes`.

**Auto ladder for eligible Gaussian LMM rows:**
`auto -> satterthwaite -> asymptotic_wald_z -> not_computed`

**Satterthwaite:** implemented for scalar LMM rows. Validated against
`lmerTestR` parity fixtures covering sleepstudy, Penicillin (crossed), and
unbalanced variants. Internal derivatives: `deviance_varpar`, `vcov_beta_varpar`,
`jac_vcov_beta_varpar`, `vcov_varpar`.

**Kenward-Roger:** explicit/opt-in only; not part of the auto ladder in
schema v1.0.0. Supports dense Sigma/G decomposition, adjusted fixed-effect
covariance, `Lb_ddf` denominator-df, scalar KR t rows, multi-df KR F rows.
Current multi-df row parity tracks the unscaled pbkrtest F statistic; rows
with non-unit F.scaling are not yet using pbkrtest's scaled F output.
KR is LMM REML only; explicit KR on ML fits must not silently degrade.

**Asymptotic Wald z:** labeled fallback; always available when estimates and
SE are finite.

**Bootstrap (fixed-effect null):** requires a certified
`fixed_effect_null` target. Produces continuity-corrected p-values with MC SE
notes. Schema `mixedmodels.bootstrap_run` v1.0.0. Explicit only; not selected
by `auto`.

**Wald confidence intervals (`wald_confint`):** listed as a Phase C v1.0
blocker — NOT YET IMPLEMENTED in the engine. Expected on `CoefTable` and
`MixedModelFit`. This is the engine ceiling: R cannot expose it until it lands.

**Row-level refusal categories (not global bans):** rank-deficient fixed
effects, non-estimable contrasts, missing/non-positive SE, predictive or
regularized fit intent, derivative prerequisites unavailable, boundary or
reduced-rank cases where the requested method is indefensible.

### 2.2 Model Comparison

**LRT policy** (`model_comparison_policy.md`):
- Adjacent pair must share response values, family/link, likelihood criterion,
  and nested fixed/random column spaces; df must increase in supplied order.
- REML fixed-effect comparisons: engine reports `requires_ml_refit` and does
  NOT perform the refit.
- Boundary-sensitive RE comparisons: use parametric-bootstrap LRT route.

**Information criteria:** `ModelComparisonTable` reports label, n, df, logLik,
deviance, AIC, BIC, delta AIC, delta BIC. LRT fields are null for non-nested
models with reason code `non_nested_models_lrt_invalid`.

**Not comparable:** different responses, families, links, or mixed ML/REML.
Such rows carry `information_criteria_available = false` with a stable reason
code.

**R wrapper exposes:** `compare()`, `anova.mm_lmm()`, `drop1.mm_lmm()`,
parametric bootstrap LRT between two nested LMMs. These consume the Rust
`ModelComparisonTable` payload.

---

## 3. Formula / Random-Effects Layer

**Accepted surface syntax (v0 contract):**
`(re | g)`, `(re || g)`, `(1 + x | g)`, `(0 + x | g)`,
grouping forms: `Single`, `g1 & g2`, `g1:g2`, `g1/g2` (expanded),
`g1*g2` (expanded).

**Canonicalization rules R1–R9:** nesting expansion, crossing expansion,
zerocorr as flag, intercept policy first-class, duplicate term detection,
conflicting covariance refusal, same-grouping different-basis preserved,
fixed/random redundancy detection, source-syntax preservation.

**Categorical random slopes:** treatment-coded and cell-means (`0 + factor`)
parameterizations are both implemented and audited.

**Deferred (vNext):**
- Smooth terms (`s()`, `te()`, `bs()`, `ns()`).
- Offsets and `I(...)` literal-protection (tokenizer rejects `^` and `,`).
- Backtick identifiers.
- Formula-level contrast specification.
- Roles declaration syntax.
- Residual-structure formula (`ar1(time, subject)` etc.).

**Known silent-acceptance bugs in the parser (v1.0 Phase B blockers):**
trailing `+`/`-`, `y ~ -x1` parses as `1 + x1`, subtraction treated as
addition, adjacent RE blocks without `+`, numeric-literal terms (`2*x1`).
These are engine-level bugs that cap correctness until Phase B lands.

**`FixedTerm::Nested` variant:** dead code — silently dropped in
`fixed_design.rs:853`. Engine-level gap.

---

## 4. Optimizer Profiles

### 4.1 Default (NLopt) Profile

`cargo build --release` — enables NLopt BOBYQA/NEWUOA for LMM. Rust still
owns the profiled objective, PLS factorization, diagnostics, and result
surface. Fastest iteration efficiency.

### 4.2 Native TrustBQ Profile (no-default-features)

`cargo build --release --no-default-features` — dependency-light, pure-Rust.

| Model family | Theta dim | Cross terms | Max evals | Stop policy |
|---|---:|---|---:|---|
| Small theta / vector RE | d <= 3 | full quadratic | 1000 | numeric ftol + stable theta |
| Moderate theta | 4 <= d < 7 | diagonal only | 1000 | numeric ftol + stable theta |
| Crossed / large theta | d >= 7 | diagonal only | 475 | statistical stall band |

TrustBQ stop classes: `smooth_convergence`, `statistical_stall`,
`certificate_accepted`, `budget_exhaustion`.

KKT-guided invalid-boundary restarts are orchestrated above TrustBQ (not
inside it) and surface as recovered convergence in the optimizer certificate.

**R wrapper:** `mixeff` uses `default-features = false` for CRAN (no NLopt
dependency). This means CRAN builds always use TrustBQ for multi-theta LMMs
and the native fallback for GLMMs.

### 4.3 PRIMA

Optional development backend; not a required runtime dependency. Feature-gated.

---

## 5. Documented Refusals and Stability Labels

| Item | Refusal / Stability |
|---|---|
| Multivariate response `cbind(y1,y2)~` | Explicit vNext/post-1.0 deferral. Requires API-breaking surgery on `FeMat` and `Formula.response`. |
| Profile likelihood for GLMM | Explicitly out of scope for v1.0. |
| GLMM parametric bootstrap | Refusal stub for ALL families in current engine. Phase C v1.0 blocker. |
| REML beta profile likelihood | Engine explicitly refuses. |
| REML fixed-effect LRT | Engine reports `requires_ml_refit`; does not perform the refit. |
| KR for crossed/nested beyond scalar-test scope | Out of scope for v1.0. |
| `n_agq > 1` with multiple or vector RE terms | Hard `invalid_agq_request` diagnostic before any optimizer work. |
| InverseGaussian / Normal-as-GLMM bootstrap | Explicit refusal in parametric bootstrap. |
| `fast = false` without NLopt | Engine refuses explicitly rather than silently falling back. |
| AR(1) / spatial residual covariance | vNext; formula syntax not parsed. |
| `I()` / `poly()` / `splines` | Not parsed; tokenizer rejects `^` and `,`. |
| GAM smooths (`s()`, `te()`) | Not parsed in v0. |
| `wald_confint` | Not yet implemented; Phase C blocker. |
| GLMM Satterthwaite / KR inference | Explicitly unsupported by contract. |

---

## 6. Engine Surfaces Not Yet Surfaced by the R Wrapper

These are capabilities the engine exposes but the R wrapper does not yet
fully consume or expose, constituting wrapper-side gaps.

### 6.1 `fixed_effect_inference_table` from Fitted Artifacts (partial gap)

**Engine:** `CompiledModelArtifact.fixed_effect_inference_table` is populated
on fitted LMM artifacts with full Satterthwaite/Wald rows and row-level
status/reason/reliability.

**R wrapper:** `inference_table.mm_lmm()` in `R/revive.R` now reads
`fit$artifact$fixed_effect_inference_table` and parses it via
`mm_json_parse_fixed_effect_inference_table()`. This is partially consumed
but `summary.mm_lmm()` building of p-value columns still relies on this
pathway and needs to be verified. The upstream support report (§ "Current
mixeff Mismatch") identified the old `method = "unavailable"` fallback in
`R/revive.R` as stale — the code was updated but full test coverage of
available Satterthwaite rows in summary output needs confirmation.

### 6.2 Explicit Kenward-Roger via Contrast Bridge Endpoint

**Engine:** KR is implemented for eligible Gaussian REML LMM hypotheses.

**R wrapper:** `contrast()` and `test_effect()` accept `method = "kenward_roger"`
and route to `mm_rust_contrast_table()`. The bridge endpoint `mm_fixed_effect_contrast_json()`
is described in the upstream support report as the recommended surface for new
contrast rows not in the fitted coefficient table. Whether this is fully wired
for KR versus Satterthwaite needs runtime verification.

### 6.3 Bootstrap Fixed-Effect Null Rows

**Engine:** schema `mixedmodels.bootstrap_run` v1.0.0 exists with certified
fixed-effect-null target construction, continuity-corrected p-values, and
MC SE notes.

**R wrapper:** `contrast(..., method = "bootstrap")` and `test_effect(...,
method = "bootstrap")` route to `mm_rust_term_bootstrap_row()`. `bootstrap_control()`
is exported. The bridge exists but the upstream report warns R must not compute
bootstrap p-values from full-model bootstrap replicates — it must only print
when Rust supplies a `bootstrap` row with `status = "available"`.

### 6.4 `CoefTable` / `ModelSummary` Satterthwaite Row Exposure

**Engine:** `CoefTable` and `ModelSummary` never directly surface
Satterthwaite/KR rows — clients must detour through
`FixedEffectInferenceTable`. This is documented as a Phase C gap in the
v1.0 roadmap ("downstream-misleading").

**R wrapper impact:** `summary.mm_lmm()` must join inference rows from
`fixed_effect_inference_table` rather than from `CoefTable` directly.

### 6.5 Schema Negotiation for `fixed_effect_inference_table`

**Engine:** schema `mixedmodels.fixed_effect_inference_table` v1.0.0.

**R wrapper:** the upstream support report (§ "Schema Negotiation") lists
this as needing explicit addition to the R manifest's known schema list. The
wrapper currently knows `formula/v0`, `mixedmodels.compiled_model_artifact/1`,
`mixedmodels.model_audit_report/2`, `mixedmodels.random_term_card/1`.

### 6.6 Random Term Cards and Cross-Card Constraints

**Engine:** `ModelAuditReport` schema v2 exposes `random_term_cards` and
`cross_card_constraints` with `english` (upstream-authored prose),
`design_support` (levels, min/median rows, within-group variation, status),
and `role_origin`.

**R wrapper:** `R/reporting.R` consumes `audit$random_term_cards` and
`audit$cross_card_constraints`, renders `block$english` verbatim, surfaces
`design_support$status`, and displays cross-card constraints. This is
substantially implemented.

### 6.7 Wald Confidence Intervals (`wald_confint`)

**Engine:** Phase C v1.0 blocker — NOT YET IMPLEMENTED. `CoefTable` does not
yet have `lower`/`upper` CI columns.

**R wrapper:** `confint.mm_lmm()` exposes profile likelihood (`method = "profile"`),
bootstrap (`method = "bootstrap"`), and Wald (`method = "wald"`). The Wald path
is upstream-blocked until the engine lands `wald_confint`.

### 6.8 GLMM Parametric Bootstrap

**Engine:** explicit refusal stub for all families. Phase C v1.0 blocker.

**R wrapper:** no GLMM bootstrap route is exposed (correct given the engine
ceiling).

### 6.9 Lexicographic Level Ordering

**Engine:** v0 contract specifies lexicographic ordering for composite
grouping keys; implementation still uses first-appearance order
(`non-conforming` per `random_effects_formulas.md` Appendix A §5.4).

**R wrapper impact:** level ordering affects parmap, parameter ordering, and
output tables. The mismatch is upstream-driven; the R wrapper cannot work
around it.

### 6.10 `||` Centering Reference

**Engine:** v0 contract requires recording a declared or weighted-mean
reference for zero-correlation terms and back-transforming user-facing
quantities. Implementation is listed as missing in the contract Appendix A.

**R wrapper impact:** `||` syntax is accepted and fits, but the centering
reference rule is not enforced by the engine, so coefficients and predictions
may be on an undocumented implicit zero scale.

---

## 7. Feature-Gated Capabilities

| Feature flag | What it enables |
|---|---|
| `nlopt` (default) | NLopt-backed BOBYQA/NEWUOA for LMM; joint Laplace/AGQ for GLMM `fast=false` |
| no default features | TrustBQ for multi-theta LMM, native COBYLA/PatternSearch for GLMM |
| `prima` | Optional PRIMA development backend; not a runtime dependency |
| `unstable-internals` | Enables internal benchmarks and probes (not a public surface) |

**CRAN build:** `mixeff` uses `default-features = false`. All CRAN users get
TrustBQ for LMM and native COBYLA/PatternSearch for GLMM. NLopt is not
available on CRAN builds. `fast = false` GLMM (joint Laplace/AGQ) is
**not accessible on CRAN**.

---

## 8. Numerical Risks That Cap Parity Claims

From `v1_0_release_roadmap.md` Phase B (not yet closed):

- **Pivoted QR uses Modified Gram-Schmidt** (not LAPACK Householder). On
  near-rank-deficient designs the kept columns differ from lme4/Julia.
  Parity claims for rank-deficient designs remain provisional.
- **Hard-coded `1e-30` zero-clamps** in `solve_scaled_vsize2_row` and
  `rdiv_lower_transpose` diverge from the policy-controlled tolerance.
  Silent behavior at sub-microvariance sigma scales.
- **User-reachable `.expect()`** in `recompute_a_blocks` — can panic from
  PIRLS via `update_irls_weights`.
- **Train/predict factor consistency:** predict-time rebuild does not carry
  training `levels` forward; can silently reorder dummy columns.
- **NaN/Inf not rejected** in `DataFrame::add_numeric`; propagates into
  Cholesky.
- **Reference level is first-appearance order** (lme4 uses alphabetical) —
  silent parity hazard for categorical coefficients.

---

## 9. Summary Table: Engine vs R-Wrapper Coverage

| Capability | Engine status | R wrapper status | Classification |
|---|---|---|---|
| Gaussian LMM fit (REML/ML) | production | exposed | works |
| GLMM Bernoulli/Binomial logit/probit/cloglog | production | exposed | works |
| GLMM Poisson log/sqrt | production | exposed | works |
| GLMM Gamma log | production | exposed | works |
| Satterthwaite inference (LMM) | production (schema 1.0.0) | partially surfaced via `inference_table()` and `contrast()` | partial |
| Kenward-Roger inference (explicit, LMM REML) | production (opt-in) | routed via `contrast()` / `test_effect()` | partial |
| Asymptotic Wald z inference | production | labeled fallback in inference rows | works |
| Bootstrap fixed-effect null (LMM) | production (explicit) | `contrast(method="bootstrap")`, `bootstrap_control()` | partial |
| Boundary-sensitive LRT (variance component) | production | `test_random_effect()`, `compare(method="boundary_lrt")` | works |
| Parametric bootstrap LRT (LMM) | production | `compare(method="bootstrap")` | works |
| Profile likelihood CI (LMM) | production (ML beta only) | `confint(method="profile")` | works |
| Wald CI | NOT YET IMPLEMENTED | upstream-blocked | upstream-blocked |
| AIC/BIC model comparison | production | `compare()`, `anova.mm_lmm()` | works |
| GLMM parametric bootstrap | refusal stub | not exposed | out-of-scope-by-design |
| GLMM fast=false (joint Laplace/AGQ) | NLopt-gated | not on CRAN builds | upstream-blocked |
| Multivariate cbind(y1,y2)~ | deferred vNext | not exposed | out-of-scope-by-design |
| Profile likelihood CI (GLMM) | out of scope v1.0 | not exposed | out-of-scope-by-design |
| KR for crossed/nested beyond scalar | out of scope v1.0 | not exposed | out-of-scope-by-design |
| AR(1)/spatial residual covariance | vNext | not exposed | out-of-scope-by-design |
| I() / poly() / splines | not parsed | not exposed | out-of-scope-by-design |
| GAM smooths | not parsed | not exposed | out-of-scope-by-design |
| Random term cards (audit) | production (schema v2) | substantially consumed | works |
| Cross-card constraints | production | consumed in reporting | works |
| Lexicographic level ordering | non-conforming (still first-appearance) | cannot work around | upstream-blocked |
| `||` centering reference | missing implementation | cannot work around | upstream-blocked |
| Formula parser silent-accept bugs | Phase B blockers | cannot work around | upstream-blocked |
