# Crate Survey: GLMM Engine Contract (crate-4)

Sources: `docs/glmm_support_contract.md`, `docs/certified_joint_glmm_optimizer_contract.md`,
`docs/difficult_model_certification.md`, cross-checked against `R/glmm.R`, `R/inference.R`,
`R/inference-options.R`, `R/predict.R`, `R/methods-print.R`, `R/methods-summary.R`, `NAMESPACE`.

Survey date: 2026-05-31.

---

## 1. Supported Model Classes

The engine (as of this contract snapshot) certifies **one model class** for GLMM:

| Class | Status |
|-------|--------|
| Generalized linear mixed model (GLMM) | Certified for listed families/links (see §2) |
| LMM | Separate contract; REML/Satterthwaite/KR are LMM-only |
| Multivariate / cbind(y1,y2) | Not supported (PRD non-goal) |
| Nonlinear (nlmer) | Not supported |

---

## 2. Certified Families and Links

### Fully certified (1.0 SemVer surface):

| Family | Links |
|--------|-------|
| Bernoulli | logit, probit, cloglog |
| Binomial | logit, probit, cloglog (includes grouped-binomial trial weights) |
| Poisson | log, sqrt |
| Gamma | log |

### Implemented but NOT certified (experimental, outside SemVer):

| Family | Reason |
|--------|--------|
| InverseGaussian | Not validated to cross-language parity standard; finite-sample inference surface intentionally incomplete; parametric bootstrap explicitly refuses it |
| Gaussian-GLMM (non-identity link) | Same as above |

The R wrapper (`R/glmm.R: mm_glmm_supported_family_links()`) only exposes the certified surface:
`binomial` {logit, probit, cloglog}, `poisson` {log, sqrt}, `Gamma` {log}. The experimental
families are correctly excluded. Note: the engine contract also names "Bernoulli" as a distinct
family, but the R wrapper folds Bernoulli under `binomial` (standard R convention; a 0/1
response with binomial family is functionally Bernoulli).

---

## 3. Offsets and Weights

The engine contract states:
- Offsets are supported as fixed linear-predictor offsets.
- Observation weights are supported where family semantics define them, including binomial
  trial weights.

**Wrapper gap**: `R/glmm.R` currently treats `weights` as reserved and immediately errors with
`mm_fit_error` if non-NULL. Offsets are not plumbed through the GLMM call path at all. The FFI
call (`wrap__mm_fit_glmm_json`) does not accept an offsets argument. Both capabilities exist in
the engine but are not yet surfaced in the R wrapper.

---

## 4. Fitting Methods and Approximation Semantics

### Fast-PIRLS (profiled) path — `method = "pirls_profiled"` / `fast = true`:

- Default and only available method in the CRAN build (no nlopt feature).
- Profiles fixed effects through PIRLS while optimizing covariance parameters on the profiled
  GLMM objective.
- Matches MixedModels.jl `fast=true` family; diverges from `lme4` joint-estimation coefficients.
- `response_constants = dropped` convention (profiled objective).
- Labelled `documented_divergence` in the parity scorecard — deliberate, not a failure.

### Joint Laplace / Joint AGQ path — `method = "joint_laplace"` / `fast = false`:

- Requires the optional `nlopt` feature.
- Estimates `[β; θ]` jointly on the full Laplace or AGQ marginal deviance with response
  constants retained (`response_constants = included`).
- Falls back deterministically to fast-PIRLS if the joint step fails to certify; fallback is
  labelled, never silently promoted.
- Row-scoped certification: `culcitalogreg` (Laplace + AGQ) are promoted; `cbpp`, `contraception`,
  `cbpp`, `verbagg` remain `documented_divergence` as of this contract.
- AGQ (`n_agq > 1`) accepted only for exactly one scalar random-effects term; multi-term or
  vector-valued terms must reject before any optimizer evaluations.

**Wrapper behavior**: `method = "joint_laplace"` is accepted as an argument but the CRAN build
refuses it with a typed `mm_fit_error` because nlopt is not compiled in. `nAGQ > 1` with
`method = "joint_laplace"` is also explicitly rejected before the FFI call. The `pirls_profiled`
nAGQ handling passes the value through but the profiled path does not genuinely honour AGQ;
the engine records the effective n_agq in the fit summary.

---

## 5. Inference Methods

### Certified for GLMM:

| Method | Engine support | R wrapper | Notes |
|--------|---------------|-----------|-------|
| Asymptotic Wald-z | Yes (from stored vcov) | `summary(fit, tests="coefficients")` | Standard for GLMM |
| `mm_lincomb` Wald-z | Yes | `mm_lincomb.mm_glmm` (asymptotic only) | Documented refusal for non-asymptotic |

### Explicitly NOT certified for GLMM (engine contract, §Approximation Semantics):

| Method | Reason |
|--------|--------|
| REML | LMM-only by contract |
| Satterthwaite df | LMM-only by contract |
| Kenward-Roger | LMM-only by contract |
| Bootstrap (parametric, fixed-effect-null) | LMM-only Rust paths; InverseGaussian/Normal explicitly refused |
| Bootstrap LRT | LMM-only Rust paths |
| Profile-likelihood CIs | Separate `lmm_profile_confint` FFI exists but has no GLMM equivalent |

**Wrapper alignment**: `test_effect` and `contrast` dispatch only for `mm_lmm`; there are no
`test_effect.mm_glmm` or `contrast.mm_glmm` methods. `confint.mm_glmm` is not defined;
`confint` would fall through to a non-method error. `inference_options` is `mm_lmm`-only.
These refusals are correct per the engine contract.

---

## 6. Diagnostic Codes (Engine Contract)

The engine contract defines six stable diagnostic codes for GLMM:

| Code | Meaning |
|------|---------|
| `optimizer_nonconvergence` | Optimizer stopped without acceptable convergence criterion |
| `invalid_agq_request` | Rejected AGQ shape (multi-term or vector RE with n_agq > 1) |
| `pirls_failure` | Final PIRLS update failed after optimizer selection |
| `boundary_parameter` | Theta on lower bounds |
| `near_unit_random_effect_correlation` | Absolute correlation near 1 |
| `binomial_separation` | Conservative fixed-effect separation diagnostic |

Additionally, five distinguishable failure modes must be separately readable from artifact alone:
optimizer failure, approximation gap, weak identification, response-constant convention, and
separation-like behavior.

**Wrapper alignment**: The R wrapper stores `fit_status` from the engine artifact and prints it
via `print.mm_glmm` and `summary.mm_glmm`. The optimizer certificate is accessible via
`optimizer_certificate(fit)` and `diagnostics(fit)`. The artifact JSON is stored on `fit$artifact`.
However, the specific diagnostic codes (`binomial_separation`, `pirls_failure`, etc.) are not
surfaced as typed R-level warnings or conditions — they live in the artifact JSON and require
manual inspection. There is no R-level helper that extracts and presents the five failure-mode
signals as structured R objects.

---

## 7. Optimizer Policy

| Condition | Default build (CRAN; no nlopt) | Full build (nlopt enabled) |
|-----------|-------------------------------|---------------------------|
| LMM optimizer | TrustBQ | TrustBQ |
| GLMM optimizer | COBYLA / PatternSearch (native) | NLopt (enables joint Laplace/AGQ) |
| PRIMA | Not required at runtime | Optional development backend |

The artifact must record: optimizer name and backend, objective approximation boundary,
n_agq, optimizer certificate status, return code, objective value, function evaluations,
boundary evidence, and diagnostics.

---

## 8. Prediction

| Route | Engine support | R wrapper |
|-------|---------------|-----------|
| LMM conditional predictions (`predict_new`) | Yes (Rust `predict_new` FFI) | `predict.mm_lmm` with re.form |
| LMM population predictions (fixed-only) | Yes (R-side fixed design) | `predict.mm_lmm` with `re.form = NA` |
| GLMM predictions | Not yet certified | `predict.mm_glmm` → `mm_inference_unavailable` |

The engine's GLMM prediction capability is not documented as certified in these contracts.
The wrapper explicitly refuses with a stable typed error.

---

## 9. LMM-Only Features (Not Applicable to GLMM)

The following are LMM-only per engine contract and correctly not wired for GLMM:

- REML estimation
- Satterthwaite / Kenward-Roger degrees of freedom
- Parametric bootstrap (fixed-effect-null and full-model)
- Bootstrap LRT
- Cluster bootstrap
- Profile-likelihood confidence intervals (`mm_lmm_profile_confint_json`)
- Conditional variance of random effects (`mm_lmm_cond_var_json`) — `ranef(condVar=TRUE)` on
  `mm_glmm` returns typed `random_effect_conditional_variance_unavailable_for_glmm`
- KKT-guided boundary restart recovery (LMM covariance-space mechanism only)
- `test_random_effect` (boundary LRT) — only dispatches for `mm_lmm`
- `simulate.mm_lmm` — only wired for LMM

---

## 10. Wrapper-Side Gaps (Engine Capability Not Yet Surfaced in R)

These are capabilities the engine exposes or the contract mandates that the R wrapper does
not yet surface:

1. **Offsets** — engine supports fixed linear-predictor offsets; R wrapper errors if any
   non-NULL offset-like argument is passed; no FFI slot for offsets in GLMM call.

2. **Weights for Binomial grouped-binomial** — engine supports observation weights including
   binomial trial weights; R wrapper treats `weights` as reserved and errors if non-NULL.

3. **`estimation_method` metadata field** — engine mandates this stable field
   (`fast_pirls_profiled`, `joint_laplace`, `joint_agq`, `fallback_fast_pirls`) in fit
   summary payloads. The R wrapper stores `fit$method` (the requested method string), not the
   engine-resolved `estimation_method`. Under fast-PIRLS the two coincide; under a joint
   fallback they would diverge. No R-level accessor exposes this field.

4. **`fallback_status` field** — engine mandates recording when joint attempt returned
   deterministic fast-PIRLS fallback. Not exposed as a named R-level field.

5. **`response_constants` field** (`dropped` vs `included`) — engine mandates this for
   objective comparability. Not accessible as a named R-level accessor; lives in raw artifact.

6. **`objective_definition` field** — engine mandates (`profiled_glmm_deviance`,
   `joint_glmm_laplace_deviance`, `joint_glmm_agq_deviance`). Not exposed as named field.

7. **Structured diagnostic code surfacing** — the six stable diagnostic codes are emitted into
   the artifact JSON but are not promoted to typed R conditions, warnings, or a structured
   diagnostic table. A user must inspect `fit$artifact` manually.

8. **`InverseGaussian` / Gaussian-GLMM experimental path** — exists in engine but is not
   exposed in the R wrapper (correct per the 1.0 certified surface; noting it here as a
   future-surfacing candidate).

9. **`inference_options` for `mm_glmm`** — no method exists; the audit verb is LMM-only.
   For GLMM the available inference is limited to asymptotic Wald-z, but there is no
   `inference_options.mm_glmm` to document this.

---

## 11. Feature-Gated Capabilities

| Feature gate | Capability | Build status |
|-------------|------------|-------------|
| `nlopt` (default in upstream; disabled in CRAN build) | Joint Laplace estimation, Joint AGQ estimation | Off in CRAN wrapper; wrapper correctly refuses with typed error |
| `prima` (optional dev backend) | PRIMA optimizer backend | Not a runtime requirement; not exposed in R wrapper |

---

## 12. Parity Status Summary

| Dataset/row | Classification | Notes |
|-------------|---------------|-------|
| `culcitalogreg` (Laplace, AGQ) | `release_blocking_parity` (joint path, nlopt-gated) | Promoted through joint gate |
| `cbpp` | `documented_divergence` | Below promotion line on fitted estimates |
| `contraception` | `documented_divergence` | fast-PIRLS / profiled-objective row |
| `verbagg` | `documented_divergence` | fast-PIRLS / profiled-objective row |
| `gopherdat2` | `documented_divergence` | Near-zero covariance; diagnostic threshold gap |
| `grouseticks` | `performance_known_slow` | Numerical claim separate; performance issue tracked |
