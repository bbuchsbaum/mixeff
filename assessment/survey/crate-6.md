# Crate Survey 6 — Engine Contract vs R Wrapper Parity

Sources read:
- `docs/profile_likelihood_json_contract.md`
- `docs/fixed_effect_p_value_validation.md`
- `docs/glmm_support_contract.md`
- `docs/mixed_model_compiler_inference_contract.md`
- `docs/kenward_roger_contract.md`
- `docs/bootstrap_fixed_effect_contract.md`
- `docs/satterthwaite_scalar_contract.md`
- `docs/optimizer_profiles.md`
- R wrapper: `R/inference.R`, `R/glmm.R`, `R/methods-summary.R`, `R/predict.R`, `R/inference-options.R`

---

## 1. Engine Model Classes and Families

### LMM (Gaussian)

- Profiled ML/REML objective, KKT-certified optimizer.
- Fixed-effect inference methods supported: `asymptotic_wald_z`, `satterthwaite`, `kenward_roger`, `bootstrap` (fixed_effect_null target), `boundary_lrt` (for variance components).
- Auto ladder: `satterthwaite -> asymptotic_wald_z -> not_computed` (KR is opt-in in schema 1.0.0).
- Profile-likelihood CIs: ML fits include sigma, theta, and beta intervals; REML fits include sigma and theta only — beta profile intervals explicitly omitted by engine contract.
- Optimizer backends: NLopt (default, performance) and native TrustBQ (`--no-default-features`, CRAN-safe).
- FitStatus enum: `ConvergedInterior`, `ConvergedBoundary`, `ConvergedReducedRank`, `ConvergedPenalised`, `NotIdentifiable`, `NotOptimized`.

### GLMM

Certified for 1.0:
- Bernoulli / Binomial: logit, probit, cloglog links
- Poisson: log, sqrt links
- Gamma: log link

**Implemented but NOT certified (experimental):**
- `InverseGaussian` — exists in engine, not validated to cross-language parity standard; parametric bootstrap explicitly refuses it; not part of SemVer GLMM contract.
- Gaussian-GLMM non-identity link paths — same status.

**Approximation modes:**
- `fast = true` (default, profiled fast-PIRLS): available without NLopt.
- `fast = false` joint Laplace (`n_agq <= 1`) and joint AGQ (`n_agq > 1`, single-scalar RE only): require NLopt feature. Without NLopt, explicit `unsupported` response (no silent fallback).
- AGQ `n_agq > 1` is refused for vector-valued or multiple RE terms with a stable `invalid_agq_request` diagnostic.

**GLMM inference:** REML, Satterthwaite, and Kenward-Roger are LMM-only by contract. GLMM artifacts must report finite-sample LMM inference as unsupported. Asymptotic Wald-z is the only certified GLMM inference method.

**GLMM fit metadata the engine exposes (contract-stable fields):**
- `estimation_method`: `fast_pirls_profiled`, `joint_laplace`, `joint_agq`, `fallback_fast_pirls`
- `objective_definition`: `profiled_glmm_deviance`, `joint_glmm_laplace_deviance`, `joint_glmm_agq_deviance`
- `response_constants`: `dropped` or `included`
- `n_agq`: requested/effective
- `fallback_status`: set only on uncertified joint attempts

**GLMM diagnostics (stable codes):**
- `optimizer_nonconvergence`, `invalid_agq_request`, `pirls_failure`, `boundary_parameter`, `near_unit_random_effect_correlation`, `binomial_separation`

**GLMM parity class distinctions (not simple lme4 parity):**
- `release_blocking_parity`, `documented_divergence`, `performance_known_slow`, `stress_opt_in`, `unsupported_with_contract`
- Current documented-divergence rows: `cbpp`, `contraception`, `verbagg` (fast-PIRLS, not lme4 joint-estimation parity); `gopherdat2` (diagnostic threshold gap); `grouseticks` (performance_known_slow).

---

## 2. Inference Methods — Engine Status

| Method | Scope | Validation status | Auto? |
|---|---|---|---|
| `asymptotic_wald_z` | LMM + GLMM | Fixture-backed, coefficient-table consistency, bounded H0 simulation smoke test | LMM fallback; GLMM default |
| `satterthwaite` | LMM only (REML or ML) | lmerTestR parity fixtures (sleepstudy, Penicillin, unbalanced, boundary/rank-deficient); auto-enabled after parity certification | Yes (first choice in auto ladder) |
| `kenward_roger` | Gaussian REML LMM only | pbkrtest parity fixtures; scalar + multi-df rows; crossed/nested certified on nlopt build only; native build: finite/plausible but not certified crossed/nested | Opt-in (not in auto ladder, schema 1.0.0) |
| `bootstrap` (fixed_effect_null) | Gaussian LMM scalar/multi-df contrasts | Null target certification contract; continuity-corrected p-value; MCSE; min 30 finite replicates | No (explicit only) |
| `boundary_lrt` | LMM variance components (random effects) | Self-Liang 50:50 mixture; ML-nested comparison; exactly one theta parameter | No (explicit only) |
| Profile-likelihood CI | LMM (ML: sigma+theta+beta; REML: sigma+theta only) | Schema-versioned JSON contract; spline internals are Rust-internal | Via `confint(method="profile")` |

---

## 3. Profile-Likelihood CI Contract Details

Schema: `mixedmodels.profile_likelihood_ci` version `1.0.0`.

- ML fits: intervals for `sigma`, `theta`, and fixed-effect `beta`.
- REML fits: intervals for `sigma` and `theta` only; beta profile explicitly omitted by contract ("until a REML beta-profile contract is certified").
- Payload fields: `schema_name`, `schema_version`, `level`, `fit_criterion`, `intervals` (with estimate, lower, upper, method, regularity, boundary_clamped flag), `profile_rows`, `notes`.
- R must call the Rust payload and deserialize JSON; must not recompute profile intervals or reinterpret spline diagnostics.

---

## 4. Fixed-Effect P-Value Validation Summary

All four analytic methods have explicit validation coverage (per `fixed_effect_p_value_validation.md`):

- Wald-z: explicit scalar tests + bounded H0 smoke test. Residual gap: add row-level table test asserting Wald rows match `coeftable()` explicitly.
- Satterthwaite: lmerTestR fixtures + boundary/rank-deficient unavailable-reason tests + H0 simulation smoke test. Residual gap: broaden simulation cases.
- Kenward-Roger: pbkrtest fixtures (scalar + multi-df) + H0 simulation smoke test. Residual gap: multi-df rows document unscaled-F parity (scaled-F support not yet in row payload when `F.scaling != 1`).
- Bootstrap: null target shape + null simulate/refit/payload row construction + continuity-corrected p-value. Residual gap: larger/adaptive bootstrap calibration deferred.
- Unsupported cases: labeled `unavailable` rows with stable reason codes for rank-deficient, predictive, regularized, post-selection, missing-SE, boundary, method-prerequisite failures.

Simulation follow-ups: bounded H0 smoke tests are done. Larger calibration studies deferred with documented rationale.

---

## 5. Optimizer Backends

| Profile | Command | LMM optimizer | GLMM optimizer | Notes |
|---|---|---|---|---|
| Default | `cargo build --release` | NLopt BOBYQA/NEWUOA | NLopt + native fallback | Performance profile; enables joint Laplace/AGQ |
| Native TrustBQ | `--no-default-features` | TrustBQ (pure Rust) | Native COBYLA/PatternSearch | CRAN-safe; no CMake/system deps; KR crossed/nested: finite but not certified against pbkrtest |
| PRIMA | optional dev backend | — | — | Not a required runtime dependency |

TrustBQ policy: small theta (d<=3) full quadratic; moderate (4<=d<7) diagonal only; crossed/large (d>=7) diagonal + 475 eval budget + statistical stall band. KKT-guided boundary restarts above TrustBQ. Certificate-aware stopping for scalar and 2x2 covariance certificates.

---

## 6. Engine Capabilities Not Yet Surfaced by R Wrapper

### 6a. GLMM `estimation_method` metadata

The engine contract specifies stable wire fields: `estimation_method`, `objective_definition`, `response_constants`, `n_agq`, `fallback_status`. The R wrapper stores `fit_result$method` and `fit_result$n_agq` but does NOT extract `objective_definition`, `response_constants`, or `fallback_status` into the `mm_glmm` object. These are available in `artifact` / `fit_result` but have no named R-side fields and no print exposure.

### 6b. GLMM `test_effect()` / `contrast()` are not implemented for `mm_glmm`

The engine contract states GLMM inference is Wald-z only (no Satterthwaite/KR), but there is no `test_effect.mm_glmm` or `contrast.mm_glmm` dispatch. The `summary.mm_glmm(tests="coefficients")` provides Wald-z, but targeted contrast testing for GLMMs (e.g., comparing levels of a factor) has no R surface. The engine's Wald-z machinery (`test_contrast_with_method(..., AsymptoticWaldZ)`) could be called from the R layer.

### 6c. GLMM `predict()` is explicitly refused

`predict.mm_glmm` raises `mm_inference_unavailable`. The engine contract does not prohibit GLMM prediction — the Rust GLMM fit stores `fitted` values and residuals. This is a deliberate R wrapper deferral, not an engine limitation.

### 6d. Profile CI restricted to LMM

`confint(method="profile")` is implemented only for `mm_lmm` objects. The engine GLMM contract does not specify a profile CI surface for GLMMs (and the PRD §3 non-goals explicitly defers GLMM profile-LL CIs to v2/out-of-scope). This is by design.

### 6e. KR parity qualification for TrustBQ build

The KR contract notes: on `--no-default-features` builds, crossed/nested KR parity against pbkrtest is NOT certified — only scalar random-slope cases are, with finiteness tolerances only. The R wrapper has no user-visible indication of which optimizer build is active or that KR may be less reliable on non-NLopt builds.

### 6f. GLMM `InverseGaussian` / non-identity Gaussian paths are silently out of scope

The engine marks these experimental. The R wrapper's `mm_glmm_supported_family_links()` simply does not include them, so any attempt raises `mm_inference_unavailable` with `unsupported_glmm_family_link`. This is correct behavior — the refusal is clean.

### 6g. GLMM five-way failure-mode distinction not exposed in R

The engine contract mandates five distinguishable GLMM failure signals: optimizer failure, approximation gap (fast-PIRLS vs joint), weak identification, response-constant convention difference, and separation-like behavior. The R `mm_glmm` object surfaces `fit_status` and `artifact` but has no structured R-side fields or print-method sections that distinguish these five from each other. A user inspecting `fit_status` sees a string; the distinguishable sub-reasons live in `artifact` without named accessors.

### 6h. Satterthwaite/KR reliability grade not displayed by default in `summary()`

The engine emits `reliability` and `reliability_reason` per inference row. The `summary.mm_lmm` print method shows a coefficient table but the reliability grade is an additional column that is suppressed unless the user inspects `inference_table(fit)$table` directly.

---

## 7. Documented Refusals (Engine-Side)

- GLMM `n_agq > 1` with vector or multiple RE terms: `invalid_agq_request` stable diagnostic.
- GLMM joint Laplace/AGQ without NLopt: explicit `unsupported` (no silent fallback to fast-PIRLS).
- REML profile beta intervals: explicitly omitted; reason code `profile_beta_unavailable_under_reml`.
- Bootstrap fixed-effect p-values from non-null targets (`full_model_distribution`, `cluster_resample`): contract prohibits; `cluster_resample` returns `not_assessed` with `bootstrap_cluster_resample_p_value_unavailable`.
- KR on ML fits: returns `method = kenward_roger` unavailable row with `kenward_roger_requires_reml`.
- KR on GLMM: GLMM artifacts must report finite-sample LMM inference as unsupported.
- Rank-deficient / aliased contrasts: `not_estimable` status.
- `ConvergedPenalised` fits: Wald/profile statistics must not be promoted as MLE-valid.

---

## 8. Stability Labels

- Schema `1.0.0` methods (covered by SemVer contract): `asymptotic_wald_z`, `satterthwaite`, `boundary_lrt`, profile-likelihood CI, bootstrap with certified `fixed_effect_null` target.
- Opt-in in schema 1.0.0 (not in auto): `kenward_roger`, `bootstrap` (explicit only).
- Experimental / not in 1.0 SemVer contract: `InverseGaussian`, Gaussian-GLMM non-identity link, joint AGQ for multi-scalar RE, PRIMA optimizer backend.
- Deferred to v2 / out-of-scope (per PRD §3): GLMM profile-LL CIs, multivariate cbind(), KR beyond scalar, AR(1)/spatial residual covariance, nlmer, GAM smooths, model-selection engine.
