# Crate Survey 8 — Engine Contract vs R Wrapper

Sources: `cross_engine_parity_scoreboard.md`, `julia_parity_fixture_drift.md`,
`difficult_model_certification.md`, plus supporting docs:
`glmm_support_contract.md`, `certified_joint_glmm_optimizer_contract.md`,
`compiler_contract_v0_prd.md`, `satterthwaite_scalar_contract.md`,
`kenward_roger_contract.md`, `bootstrap_fixed_effect_contract.md`,
`profile_likelihood_json_contract.md`, `boundary_lrt_variance_component_contract.md`,
`compiler_verdicts.md`, `random_effects_formulas.md`.

---

## 1. Supported Model Classes

### LMM (Gaussian)

| Feature | Engine status | R wrapper status |
|---|---|---|
| Fit (REML/ML) | Certified, parity fixtures for sleepstudy / penicillin / pastes scales | Surfaced via `lmm()` |
| Weights (case) | Supported | Surfaced via `lmm(weights=)` |
| Interior / boundary / reduced-rank fits | All certified; KKT certificates for scalar and 2x2 blocks | Surfaced: `fit_status`, `is_singular()`, `changes()` |
| `verify_convergence()` restart/jitter/alternate-optimizer checks | Implemented in Rust, structured result attached to optimizer certificate | **NOT surfaced** in R wrapper; no `verify_convergence()` R API |
| Profile-likelihood CIs | Schema `mixedmodels.profile_likelihood_ci 1.0.0`; ML: sigma/theta/beta; REML: sigma/theta only (beta deferred) | Surfaced via `confint(method="profile")` using `mm_lmm_profile_confint_json` |
| Satterthwaite scalar contrast/term | Implemented; auto-fallback ladder: `satterthwaite -> asymptotic_wald_z -> not_computed`; parity fixture vs lmerTest | Surfaced via `contrast()`, `test_effect()`, `df_for_contrast()` |
| Kenward-Roger scalar + multi-df | Implemented; requires REML; certified for sleepstudy + penicillin + pastes; opt-in (not auto default); pbkrtest parity fixture | Surfaced via `contrast(method="kenward_roger")`, `test_effect(method="kenward_roger")` |
| Bootstrap fixed-effect null (parametric) | Implemented; certified `fixed_effect_null` target; `method=bootstrap` scalar and multi-df rows | Surfaced via `contrast(method="bootstrap", bootstrap=bootstrap_control(...))`, `test_effect(method="bootstrap")` |
| Bootstrap full-model distribution / confint | Implemented | Surfaced via `confint(method="bootstrap")` |
| Bootstrap LRT (parametric) | Implemented via `mm_bootstrap_lrt_json` FFI | Surfaced via `test_effect(method="bootstrap_lrt")` |
| Boundary LRT (variance component) | Self-Liang 50:50 mixture; one-parameter only; refuses fixed-effect requests | Surfaced via `test_random_effect(method="boundary_lrt")` |
| Cluster bootstrap (estimator distribution) | Implemented as distribution target; p-values deliberately `not_assessed` in schema 1.0.0 | Surfaced via `test_effect(method="cluster_bootstrap")` but p-values are `not_assessed` by design |
| Model comparison (AIC/BIC/LRT) | Implemented via `mm_compare_models_json` | Surfaced via `compare()` |
| Audit report | `mixedmodels.model_audit_report` schema | Surfaced via `audit_design()` |
| explain_model | Pre-fit and post-fit | Surfaced via `explain_model()` |
| parameterization() drilldown | theta/Lambda/parmap/VarCorr trace | Surfaced via `parameterization()` |
| changes() | requested->semantic->supported->fitted transitions | Surfaced via `changes()` |
| emmeans integration | `emm_basis.mm_lmm` / `recover_data.mm_lmm` | Surfaced |

### GLMM

| Feature | Engine status | R wrapper status |
|---|---|---|
| Binomial (logit, probit, cloglog) | Certified | Surfaced; R maps `binomial()` family |
| Poisson (log, sqrt) | Certified | Surfaced; R maps `poisson()` family |
| Gamma (log) | Certified | Surfaced; R maps `Gamma()` family; R normalizes name to `"gamma"` |
| `InverseGaussian` and Gaussian non-identity links | Implemented but NOT certified for 1.0; parametric bootstrap explicitly refuses these families | **NOT exposed** in R wrapper; `mm_glmm_supported_family_links()` explicitly excludes them |
| `fast=true` (profiled fast-PIRLS) | Certified; `estimation_method = fast_pirls_profiled` | Surfaced as `method = "pirls_profiled"` |
| `fast=false` / joint Laplace (nlopt backend) | Row-scoped; `culcitalogreg` Laplace and AGQ promoted; `cbpp`/`contraception` still below parity line; requires nlopt feature | Exposed as `method = "joint_laplace"` but docstring says "refused in this build because vendored Rust engine compiled without nlopt backend" — **wrapper refuses joint_laplace at runtime** |
| Joint AGQ (`nAGQ > 1`) | Requires single scalar RE; refuses vector RE; requires nlopt | Accepted by `glmm(nAGQ=)` input validation but blocked in practice by the joint_laplace nlopt refusal; AGQ metadata accepted on profiled path only |
| GLMM Wald-z inference | Asymptotic; `summary(mm_glmm, tests='coefficients')` returns Wald-z table | Surfaced via `summary()` and `mm_lincomb.mm_glmm` (asymptotic only) |
| GLMM Satterthwaite / KR | **Engine contract explicitly refuses** (LMM-only) | R wrapper has no GLMM `contrast()` / `test_effect()` dispatch |
| GLMM bootstrap calibration | Deferred in engine bootstrap contract | Not surfaced |
| GLMM profile-likelihood CI | Not certified (no `mm_glmm_profile_confint` FFI) | Not surfaced |
| GLMM `confint()` | Not in NAMESPACE | **Missing**; no `confint.mm_glmm` method |
| GLMM `contrast()` / `test_effect()` | Not defined for `mm_glmm` class | **Missing** |
| GLMM emmeans | `emm_basis.mm_glmm` / `recover_data.mm_glmm` added | Surfaced |
| GLMM `predict(newdata=)` | Surfaced | Surfaced |
| GLMM diagnostics (5 distinguishable failure modes) | Engine exposes structured codes: `optimizer_nonconvergence`, `invalid_agq_request`, `pirls_failure`, `boundary_parameter`, `near_unit_random_effect_correlation`, `binomial_separation` | R wrapper exposes `fit_status` field and `artifact` diagnostics; no dedicated helper to surface the 5-way distinction; partial (test-gap) |
| GLMM `estimation_method` field | Engine emits `fast_pirls_profiled`, `joint_laplace`, `joint_agq`, `fallback_fast_pirls` | Field present in `fit_result` but **not stored as a named slot** on the `mm_glmm` object; buried in `fit` list |
| GLMM `response_constants` / `objective_definition` fields | Engine emits these as structured metadata | Not surfaced as named slots or accessor |
| GLMM `fallback_status` | Engine emits when joint attempt returns fast-PIRLS fallback | Not surfaced |

---

## 2. Formula / Random-Effects Syntax

| Feature | Engine status | R wrapper status |
|---|---|---|
| `(x \| g)`, `(1+x \| g)`, `(0+x \| g)` | Certified | Surfaced |
| `(x \|\| g)` zero-correlation | Implemented; `zerocorr=true` flag | Surfaced via lme4-style `\|\|` |
| `\|\|` centering rule (weighted mean reference) | Specified in `random_effects_formulas.md §4.5`; **non-conforming — not yet implemented** in Rust | **Not surfaced**; no centering applied or reported |
| `(1 \| a:b)`, `(1 \| a/b)`, `(1 \| a*b)` | Canonicalization rules R1/R2 implemented | Surfaced |
| Composite grouping key (collision-free `\x1E` separator) | Specified; **non-conforming** — Rust still uses `_` join | Carries through |
| Lexicographic level ordering | Specified as v0 default; **non-conforming** — first-appearance still used | Carries through; ordering unstable between row-permuted data |
| Categorical basis columns in random slope | Implemented (treatment-coded + cell-means) | Surfaced |
| Interaction basis in random slope | Partially implemented | Partially surfaced |
| Design audit diagnostics | `FixedRandomRedundant`, `RandomSlopeUnsupported`, `RepeatedUnitUnmodeled`, `CovarianceTooRich`, etc. | Surfaced via `audit_design()` |

---

## 3. Inference Methods Summary

| Method | LMM | GLMM |
|---|---|---|
| `auto` (Satterthwaite then Wald fallback) | Certified | Not applicable |
| `satterthwaite` | Certified; parity fixtures vs lmerTest | Not applicable (engine refuses) |
| `kenward_roger` | Certified (opt-in); default build (nlopt) only for crossed/nested; native optimizer "finite and plausible" not certified | Not applicable |
| `asymptotic` (Wald-z) | Available | Available (only option) |
| `bootstrap` (parametric fixed-effect null) | Certified | Not yet |
| `bootstrap` full-model distribution / confint | Available | Not yet |
| `bootstrap_lrt` | Available | Not yet |
| `boundary_lrt` (variance component) | Certified, one-parameter only | Not applicable |
| `cluster_bootstrap` | Estimator distribution only; p-values `not_assessed` | Not applicable |
| `profile` confint | Certified for LMM; REML omits beta intervals | Not implemented for GLMM |

---

## 4. Parity Classes and Cross-Engine Status

From `cross_engine_parity_scoreboard.md` and `glmm_support_contract.md`:

| Fixture / row class | Class | Notes |
|---|---|---|
| `easy_full_rank` LMM | `parity` | lme4 / MixedModels.jl / Rust all agree |
| `reduced_rank_unit_correlation` | `documented_divergence` | Rust: `ConvergedReducedRank`; lme4 + Julia: `ok`; intentional |
| `cbpp`, `contraception`, `verbagg` GLMM (fast-PIRLS) | `documented_divergence` | Profiled objective; MixedModels.jl fast=true match on large rows; not lme4 joint-estimation parity |
| `culcitalogreg` Laplace + AGQ (fast=false) | `release_blocking_parity` | Promoted through joint-optimizer gates |
| `cbpp`, `contraception` fast=false | Still `documented_divergence` | Below promotion tolerance line |
| `gopherdat2` | `documented_divergence` | Near-zero covariance without lme4 singular flag; threshold/convention gap |
| `grouseticks` | `performance_known_slow` | MixedModels.jl fast=true objective contract; known lme4 beta gap; performance issue |

Julia parity fixture drift gate (`julia_parity_fixture_drift.md`): separate script
`scripts/check_julia_parity_fixtures.sh`; not part of default Cargo test suite;
requires working Julia + MixedModels.jl environment.

---

## 5. Documented Refusals (Engine-Level)

- `InverseGaussian` and Gaussian non-identity link GLMMs: implemented but not certified; parametric bootstrap refuses them explicitly.
- `nAGQ > 1` with vector-valued or multiple RE terms: refused before optimizer; stable `invalid_agq_request` diagnostic.
- GLMM Satterthwaite / KR: explicitly unsupported; must report as unavailable.
- GLMM profile-likelihood CI: no certified contract.
- `fast=false` (joint Laplace/AGQ) without nlopt feature: explicit refusal, not silent fallback.
- Fixed-effect comparison via `boundary_lrt`: refused with `boundary_lrt_not_fixed_effect_method`.
- Multi-parameter boundary LRT mixture: refused; routes to parametric bootstrap LRT.
- Kenward-Roger on ML fits: refused with `kenward_roger_requires_reml`.
- Bootstrap p-values from full-model distribution or cluster-resample targets: refused (not a certified null simulation).

---

## 6. Feature-Gated Capabilities (nlopt / prima)

- `nlopt` feature: enables `fast=false` joint Laplace and joint AGQ for GLMMs; required for certified KR parity on crossed/nested models (native optimizer produces "finite and plausible" KR output only, not `pbkrtest`-certified). Default release build enables nlopt; CRAN build and `--no-default-features` do not.
- `prima`: optional development backend; not a required runtime dependency; not part of the public contract.
- The vendored Rust snapshot bundled in the R package is compiled **without nlopt** (per `glmm.R` docstring: "refused in this build because vendored Rust engine compiled without nlopt backend"), meaning `joint_laplace` GLMM and certified crossed/nested KR are currently unavailable to R users.

---

## 7. Engine-Exposed Features Not Yet Surfaced in R Wrapper

| Engine feature | Evidence | Gap class |
|---|---|---|
| `verify_convergence()` — bounded restart + jitter + alternate-optimizer check | `satterthwaite_scalar_contract.md`, `compiler_contract_v0_prd.md`; full `verify_convergence()` API in Rust | No R-side API; no accessor; wrapper-side gap |
| GLMM `estimation_method` named slot | Engine emits `fast_pirls_profiled` / `joint_laplace` / `joint_agq` / `fallback_fast_pirls` | Buried in `fit$fit`; no top-level named slot; partial |
| GLMM `response_constants`, `objective_definition`, `fallback_status` structured fields | `glmm_support_contract.md` §stable summary fields | Not stored as named slots on `mm_glmm`; wrapper-side gap |
| 5-way GLMM failure mode distinction | `glmm_support_contract.md` §Distinguishable Failure Modes | `fit_status` surfaces top-level status; no dedicated R helper for the optimizer/approximation-gap/weak-identification/response-constant/separation distinction |
| GLMM `confint()` | No `confint.mm_glmm` in NAMESPACE | Missing method entirely |
| GLMM `contrast()` / `test_effect()` | No dispatch for `mm_glmm` class | Missing (asymptotic Wald through `mm_lincomb` only) |
| `||` centering reference value and back-transform | `random_effects_formulas.md §4.5` specifies declared/weighted-mean centering; engine not yet conforming | Upstream-blocked (engine non-conforming) |
| Lexicographic composite grouping level ordering | `random_effects_formulas.md §5.4`; engine non-conforming (still first-appearance) | Upstream-blocked |
| Composite grouping key collision-free `\x1E` separator | `random_effects_formulas.md §5.2`; engine still uses `_` | Upstream-blocked |
| `design_compiled` automatic full-to-diagonal reductions applied at fit time | `compiler_contract_v0_prd.md`; partial — diagnostics emitted but effective-basis rewriting not yet applied | Upstream-blocked (partial engine implementation) |
| REML profile-likelihood CI for fixed-effect beta | `profile_likelihood_json_contract.md`; explicitly deferred until REML beta-profile contract certified | Upstream-blocked |

---

## 8. Stability Labels

- **Certified / release-blocking**: LMM fit (all standard covariance structures), Satterthwaite, KR (default build), parametric bootstrap fixed-effect null, boundary LRT (one-parameter), profile-likelihood CI (ML), boundary/singular/reduced-rank certificate vocabulary, GLMM fast-PIRLS for certified families.
- **Row-scoped certified**: `culcitalogreg` Laplace + AGQ fast=false; other GLMM fast=false rows remain `documented_divergence`.
- **Experimental / not certified for 1.0**: `InverseGaussian`, Gaussian non-identity link GLMMs.
- **Deferred (explicitly out of v0/v1 scope)**: AR(1)/spatial residual covariance, multivariate cbind(), GAM smooths, full Kenward-Roger beyond scalar on non-default builds.
