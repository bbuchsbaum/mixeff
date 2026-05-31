# Crate Survey: mixeff-rs Engine Contract

Survey date: 2026-05-31  
Source docs read:
- `compiler_contract_v0_prd.md`
- `mixed_model_compiler_inference_contract.md`
- `glmm_support_contract.md`
- `kenward_roger_contract.md`
- `satterthwaite_scalar_contract.md`
- `fixed_effect_p_values_plan.md`
- `bootstrap_fixed_effect_contract.md`
- `random_effects_formulas.md`
- `optimizer_profiles.md`
- `semver_policy.md`
- `certified_joint_glmm_optimizer_contract.md`
- `mixeff_upstream_support_report.md`

---

## 1. Supported Model Classes

### LMM (Gaussian REML/ML)
- Random intercepts, slopes, and full/diagonal covariance: **supported and parity-certified**
- Fit intents: `confirmatory/as_specified`, `confirmatory/design_compiled`, `exploratory`, `predictive`
- REML and ML objectives: both supported; intent-aware model comparison described
- Profiled ML/REML objective with sparse PLS factorization (lme4-style)

### GLMM
Certified families (release_blocking_parity or row-scoped joint-promotion):
- Bernoulli / Binomial with **logit**, **probit**, **cloglog** links
- Poisson with **log** and **sqrt** links
- Gamma with **log** link

Implemented but NOT certified for 1.0 (experimental, not in SemVer contract):
- `InverseGaussian`
- Gaussian-GLMM non-identity link paths
- Parametric bootstrap explicitly refuses InverseGaussian / Normal

Offsets and observation weights: supported where family semantics define them.

---

## 2. GLMM Fitting Modes

| Mode | Description | Parity class |
|------|-------------|-------------|
| `fast = true` | Profiled fast-PIRLS; inner PIRLS solves β, outer optimizer solves θ | `documented_divergence` for cbpp/contraception/verbagg vs lme4 joint-estimation |
| `fast = false` + NLopt | Joint Laplace (n_agq ≤ 1) or Joint AGQ (single scalar RE, n_agq > 1) | Row-scoped: culcitalogreg Laplace/AGQ rows promoted; cbpp/contraception still below line |
| `fast = false` without NLopt | Explicit unsupported request (does not silently fall back) | N/A |

AGQ with n_agq > 1 is accepted **only** for exactly one scalar random-effects term; vector RE or multiple RE terms reject before any optimizer evaluations with a stable `invalid_agq_request` diagnostic.

---

## 3. Inference Methods

### Fixed-Effect P-Values (LMM)

| Method | Status |
|--------|--------|
| `asymptotic_wald_z` | Implemented; low-reliability label; default for confirmatory fits where Satterthwaite fails |
| `satterthwaite` | **Implemented and certified** for scalar Gaussian LMM contrasts; parity fixtures against lmerTestR; auto ladder: `satterthwaite → asymptotic_wald_z → not_computed` |
| `kenward_roger` | **Implemented** for explicit scalar and multi-df Gaussian REML LMM requests; parity fixtures vs pbkrtest for sleepstudy + Penicillin + Pastes; opt-in in schema 1.0.0 (not auto) |
| `bootstrap` | **Implemented** for explicit `fixed_effect_null` target; scalar studentized-t and multi-df F; certified run-metadata payload; `auto` does not select bootstrap |
| GLMM Satterthwaite/KR | **Explicitly unsupported**; LMM-only finite-sample methods are marked unavailable for GLMMs |

Row-level inference table schema: `mixedmodels.fixed_effect_inference_table` version `1.0.0`.  
Auto ladder: `auto → satterthwaite → asymptotic_wald_z → not_computed` for eligible Gaussian LMM.

### Random-Effect Variance Components
- Boundary-aware LRT (Self-Liang 50:50 mixture) for nested ML comparison: implemented
- Naive ordinary p-values for random effects: refused by design

### Model Comparison
- Intent-aware REML/ML comparison: documented but described as future work in broader sense
- Boundary LRT with stable schema `mixedmodels.boundary_lrt` version `1.0.0`: implemented

---

## 4. Optimizer Profiles

### Default (NLopt enabled)
- LMM: NLopt-backed BOBYQA/NEWUOA
- GLMM: NLopt joint Laplace/AGQ when `fast = false`; native COBYLA/PatternSearch for `fast = true`

### Native TrustBQ (`--no-default-features`)
- Dependency-light; pure Rust
- LMM: scalar native optimizer for 1-theta; TrustBQ for multi-theta
- Policy table by theta dimension: full quadratic cross-terms for d ≤ 3, diagonal only for d ≥ 4, reduced budget for d ≥ 7
- PRIMA: optional dev backend, not a required runtime dependency
- KR parity: certified on default (NLopt) build only; native-optimizer build certified "finite and plausible" but not "certified to pbkrtest" for crossed/nested models

---

## 5. Formula Surface (Random-Effects)

Supported and canonicalized:
- `(x | g)`, `(1 + x | g)`, `(0 + x | g)` — intercept policy first-class
- `(x || g)` — zero-correlation diagonal; NOT invariant to additive shifts (centering reference decided but not yet implemented)
- `(1 | a/b)` → expands to `(1 | a) + (1 | a:b)` (R1, parse-time)
- `(1 | a*b)` → expands to `(1 | a) + (1 | b) + (1 | a:b)` with `CrossingLikelyUnintended` Info (R2)
- `(1 | a:b)` — cell-only, preserved
- `(1 | a) + (1 | b)` — crossed, preserved as independent blocks (R7, no default warning)
- Composite grouping: `g1:g2`, `g1 & g2`
- Categorical and interaction basis columns: treatment-coded and cell-means expansion

Not yet implemented / non-conforming:
- `||` numeric centering to weighted mean reference: policy decided, code missing
- Composite-level keys still use `_` separator (non-conforming vs decided `\x1E`)
- Lexicographic level ordering default: non-conforming (first-appearance used)

Out of scope (v0 non-goals / vNext):
- GAM smooths: `s()`, `te()`, `bs()`, `ns()`
- `I()`, `poly()`, `offset()`
- Multivariate response `cbind(y1,y2) ~`
- Residual-structure formula syntax `residual = ar1(time, group)`

---

## 6. Compiler / Audit Capabilities

### Implemented
- Semantic random-effects IR (grouping factor, basis, intercept policy, covariance form, source syntax)
- Design audit: fixed-effect rank, aliased columns, empty cells, grouping-level counts, within-group variation, information-budget reporting
- Fixed/random redundancy detection (`FixedRandomRedundant`)
- Repeated-unit under-modeling detection (`RepeatedUnitUnmodeled`)
- `maximal_feasible` v0 rule with named thresholds (configurable, deterministic)
- Covariance family sum type: `Scalar | Diagonal | FullCholesky | Structured | ReducedRank`
- Effective covariance rank summaries (rePCA-style) with user-scale loadings
- `explain_model()` prefit explanation
- `audit()` structured health report with random term cards
- `parameterization()` drilldown: source syntax → semantic → basis → θ → Λ → parmap → VarCorr
- `changes()`: requested → semantic → supported → fitted transition records
- `verify_convergence()`: bounded restart, jittered-start, alternate-optimizer checks
- Optimizer certificate: stop evidence, parameter-space context, gradient/Hessian evidence
- GLMM fit metadata: family/link, objective approximation, response-constant convention, fallback_status
- Five worked-example JSON fixtures (sleepstudy, crossed, mixture, redundant, logistic)
- Versioned JSON wire format for all contract objects (serde)
- `FitStatus` enum: `ConvergedInterior | ConvergedBoundary | ConvergedReducedRank | ConvergedPenalised | NotIdentifiable | NotOptimized`

### Partially implemented / tracked open
- Broad `design_compiled` basis-dropping/effective-rewriting: covariance reduction exists; basis-direction dropping open
- `||` centering reference recording and back-transform
- Written-vs-canonical syntax spans on source_syntax
- Exact-duplicate term merge (diagnostic emitted, merge pending)
- Real derivative-backed KKT/Hessian checks (open: bd-01KQ7X05J0PXWDCAF808479XAP)
- Random-side empty-cell diagnostics (interaction basis partial)
- Lexicographic composite-level ordering

### Documented refusals (by design)
- AGQ for vector-valued or multiple RE terms
- Explicit `satterthwaite`/`kenward_roger` for GLMM
- Explicit `kenward_roger` on ML fits (requires REML or certified REML refit path)
- `fast = false` without NLopt (returns explicit unsupported, not silent fallback)
- Confirmatory p-values for exploratory/predictive/regularized/post-selection fits
- GLMM recovery is not default/silent; any recovery must be opt-in and labeled

---

## 7. Stability Labels

| Surface | Status |
|---------|--------|
| `mixeff_rs::model` (LMM/GLMM fitting, FitOptions, inference/bootstrap/prediction entry points) | Stable post-1.0 (SemVer) |
| `mixeff_rs::stats` (varcorr, coeftable, model_summary, lrt, bootstrap, profile) | Stable post-1.0 |
| `mixeff_rs::formula` (parse_formula, Formula AST, FormulaError) | Stable post-1.0 |
| `mixeff_rs::error` (MixedModelError, Result) | Stable post-1.0 |
| `mixeff_rs::compiler` (~40+ IR types, compiled artifact, audit schemas) | **Unstable** (`unstable-internals` feature gate); JSON schemas also excluded from SemVer |
| `mixeff_rs::pathology` | Unstable internals |
| `mixeff_rs::linalg` | Private (pub(crate)) |
| All enums are `#[non_exhaustive]` | Minor version → add variant |

Deferred to v2.0 (explicitly out of 1.0 scope):
- KR beyond current scalar/multi-df tested scope
- Full `I()` / formula-level transforms
- First-class polars/arrow ingestion
- Profile likelihood for GLMM
- Multivariate response

---

## 8. Engine Capabilities Not Yet Surfaced by the R Wrapper

The following engine capabilities are implemented in Rust but the R wrapper
does not yet fully consume them (identified from `mixeff_upstream_support_report.md`
and direct wrapper inspection):

1. **`inference_table()` reads from artifact correctly** — implemented in
   `R/revive.R` (`inference_table.mm_lmm` parses
   `fit$artifact$fixed_effect_inference_table`). The legacy fallback remains for
   old objects. The old "unavailable/not_certified_by_rust" table described in
   the upstream report has been replaced.

2. **`summary.mm_lmm` uses inference rows** — `mm_summary_coefficients()` joins
   from the inference table via `inference_table()`. This is implemented.

3. **`||` centering back-transform** — engine policy decided but not yet
   implemented in Rust (`random_effects_formulas.md` §4.5 status: missing).
   R wrapper cannot surface what the engine has not yet computed.

4. **Lexicographic composite-level ordering** — engine non-conforming;
   R wrapper cannot rely on stable level ordering for composite grouping keys.

5. **Explicit KR via R `contrast()` and `test_effect()`** — wired via
   `mm_fixed_effect_contrast_json` bridge; available but KR remains opt-in
   and not part of auto ladder (schema 1.0.0 policy).

6. **Bootstrap confidence intervals (`confint(..., method = "bootstrap")`)** —
   R wrapper uses `full_model_distribution` target (not null-constrained), which
   produces intervals but NOT certified hypothesis-test p-values. The engine
   distinguishes these; R correctly labels the interval method.

7. **`vcov()` theta covariance** — engine has `vcov_varpar` but R exposes only
   NA matrix with reason `theta_covariance_unavailable`. Not yet bridged.

8. **Random term card `english` prose** — `explain.R` and `audit.R` already
   consume `random_term_cards` and `cross_card_constraints` from the audit
   report. Cards are passed through; engine-authored `english` fields are
   rendered.

9. **`mm_formula_manifest()` capability flags** — the manifest is a `.Call` to
   Rust; capability flags such as `satterthwaite`, `kenward_roger_explicit`,
   `bootstrap_fixed_effect_payload` are not yet enumerated in the R manifest
   docstring (the manifest is Rust-owned). The upstream report recommends
   exposing these as named flags; they are not yet in the wrapper docstring.

10. **`test_random_effect()` / boundary LRT** — implemented in R; calls the
    Rust boundary-LRT bridge. Wired and documented.

11. **Profile likelihood CIs (`confint(..., method = "profile")`)** — present in
    `R/inference.R`; delegates to `mm_profile_confint()`. GLMM profile CI is
    out of scope for 1.0 per semver_policy.md.

12. **Fit-intent `design_compiled` / `as_specified` / `exploratory` / `predictive`
    modes exposed to R** — accessible via `mm_control()` options but the
    user-facing documentation of these modes is thin.

---

## 9. Feature-Gated Capabilities

| Feature | Gate | Behavior without gate |
|---------|------|-----------------------|
| NLopt (default enabled) | `default-features = true` | TrustBQ pure-Rust LMM; native COBYLA/PatternSearch GLMM |
| Joint Laplace/AGQ for GLMM (`fast = false`) | Requires NLopt | Explicit unsupported, labeled |
| PRIMA backend | Optional dev feature | Not exposed in release builds |
| `unstable-internals` | Opt-in Cargo feature | Compiler IR, pathology internals hidden |

For CRAN packaging (`--no-default-features`): all LMM inference, Satterthwaite,
KR, and bootstrap remain available; joint GLMM objective not available; GLMM
fits use native fast-PIRLS path only.

---

## 10. Documented Non-Goals / Out-of-Scope by Design

From `compiler_contract_v0_prd.md` and `semver_policy.md`:
- AR(1), spatial exponential, Matern residual covariance structures
- Kenward-Roger beyond tested scalar/multi-df scope (v2 candidate)
- Automatic regularized covariance selection
- Full model-lattice search (only deterministic `maximal_feasible` v0 rule)
- Influence diagnostics requiring repeated refits
- Dense `n×n` covariance construction in any default diagnostic path
- R package interface (owned by mixeff R wrapper; Rust owns semantics)
- Multivariate response `cbind(y1,y2) ~` (vNext)
- GLMM Satterthwaite/KR (not promised in this contract)
- Naive random-effect p-values
- Profile likelihood for GLMM (1.0 deferred)
- `nlmer` nonlinear models
- Post-selection ordinary p-values without explicit unpenalized-refit contract
