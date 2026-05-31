# Rust Engine Contract Survey — crate-2

Sources read: `mixed_model_compiler_inference_contract.md`,
`compiler_contract_v0_prd.md`, `glmm_support_contract.md`,
`satterthwaite_scalar_contract.md`, `kenward_roger_contract.md`,
`bootstrap_fixed_effect_contract.md`, `boundary_lrt_variance_component_contract.md`,
`profile_likelihood_json_contract.md`, `random_effects_formulas.md` (§1–4),
`optimizer_profiles.md`, `v1_0_release_roadmap.md`.

Survey date: 2026-05-31. Crate status per roadmap: "Phase A publishability landed;
Phase B–D hardening in progress."

---

## 1. Model Classes Supported

### 1.1 LMM (Gaussian)

The engine's primary certified model class. Profiled ML/REML over the Cholesky
covariance parameterization, sparse penalized least-squares PLS solves, and the
full compiler-contract pipeline:

- Arbitrary crossed and nested random-effect grouping structures via lme4-style
  formula syntax.
- Random intercepts, random slopes (numeric and categorical basis columns),
  random interactions, zero-correlation (`||`) blocks.
- Grouping expressions: single factor, `a:b` (cell), `a/b` (nested expansion),
  `a*b` (crossed expansion), `a&b` (legacy interaction alias).
- Covariance families: Full Cholesky, Diagonal, Scalar — as a sum type
  (`ThetaMap`) with explicit family transitions recorded.
- Observation weights are supported.
- Offsets: not parsed in v0 (formula-layer non-goal).
- Response transforms `I(...)`, `poly(...)`, `log(y)`: not parsed in v0.

### 1.2 GLMM

Supported families and links (certified for schema 1.0 with `release_blocking_parity`):

| Family | Links |
|---|---|
| Bernoulli | logit, probit, cloglog |
| Binomial (grouped) | logit, probit, cloglog |
| Poisson | log, sqrt |
| Gamma | log |

`InverseGaussian` and Gaussian-GLMM non-identity paths are **implemented but
NOT certified for 1.0** — treated as experimental, not under SemVer guarantee.

Offsets are supported as fixed linear-predictor addends. Trial weights are
supported for binomial grouped proportions.

---

## 2. Estimation Methods

### 2.1 LMM

- **ML** and **REML** profiled objectives. Selection is caller-specified (`REML =
  TRUE/FALSE`).
- **Optimizer backends** (two profiles):
  - Default build (NLopt enabled): BOBYQA/NEWUOA via NLopt for multi-theta
    LMMs; scalar models use native scalar optimizer.
  - `--no-default-features` (TrustBQ): native pure-Rust trust-region optimizer.
    KR parity against `pbkrtest` is certified on default build only; TrustBQ
    crossed/nested KR output is "finite and plausible", not certified.
  - PRIMA is an optional development backend, not a required runtime dependency.

### 2.2 GLMM

- **`fast = true` / `pirls_profiled`**: profiled fast-PIRLS approximation.
  Certified default path. Matches MixedModels.jl `fast=true` family, not lme4
  joint-estimation. Objective convention: `response_constants = dropped`.
- **`fast = false` / `joint_laplace`** (NLopt feature-gated): joint Laplace
  (`n_agq <= 1`) or joint AGQ (`n_agq > 1`, single scalar RE only) when NLopt
  is enabled. Without NLopt this is an explicit unsupported request, not a
  silent fallback.
- AGQ (`n_agq > 1`) is accepted only for models with exactly one scalar
  random-effects term. Vector-valued or multi-term models reject AGQ before
  any optimizer evaluations.

---

## 3. Inference Methods

### 3.1 Fixed-Effect Inference (LMM only)

All certified inference is LMM-only. GLMM artifacts must report Satterthwaite
and Kenward-Roger as unsupported.

**Auto ladder** (schema 1.0):
```
auto -> satterthwaite -> asymptotic_wald_z -> not_computed
```
Kenward-Roger is opt-in in schema 1.0; the auto ladder may promote it in a later
major schema version.

**Satterthwaite** (`method = satterthwaite`): implemented and certified.
- Scalar fixed-effect contrasts and term-level multi-df F rows for Gaussian
  REML or ML LMMs.
- Artifacts: `deviance_varpar`, `vcov_beta_varpar`, `jac_vcov_beta_varpar`
  (central finite differences, lower-bound stencil rejection), `vcov_varpar`
  (Moore-Penrose pseudoinverse of deviance Hessian).
- Parity fixtures: `sleepstudy` (intercept, slope), `Penicillin` (crossed
  intercept), unbalanced variant, boundary/rank-deficient stubs.
- Reliability grades: `moderate`, `low`, `not_available`.
- No silent fallback: explicit Satterthwaite requests that fail prerequisites
  return labeled unavailable rows, not Wald rows.

**Kenward-Roger** (`method = kenward_roger`): implemented and certified for
supported model classes; opt-in for schema 1.0.
- Gaussian REML LMMs only. ML fits return a KR-labeled unavailable row.
- Artifacts: `KenwardRogerSigmaG`, `KenwardRogerAdjustedVcov` (PhiA = Phi + 2
  Gamma, P, Q, W, IE2 eigenvalue diagnostics), `KenwardRogerLbDdf`.
- Scalar t-rows with adjusted SE and KR denominator df.
- Multi-df F-rows with numerator df = numerical restriction rank. F-scaling
  from pbkrtest's non-unit scaling is documented but not yet applied to row
  payloads.
- Parity fixtures: `sleepstudy` scalar rows, `Penicillin` crossed intercept,
  `Pastes` nested intercept (default NLopt build only for crossed/nested).
- KR is unavailable for GLMMs, residual covariance structures, and non-REML
  fits.

**Bootstrap — fixed-effect null** (`method = bootstrap`): implemented and
certified for scalar LMM contrasts.
- Produces a certified `fixed_effect_null` target (projected beta under `L
  beta = rhs`, reuse-fitted-covariance policy). Full-model and cluster-resample
  targets do not certify hypothesis-test p-values.
- Parametric simulation from Rust-owned fitted model state; refit through Rust
  LMM engine.
- Continuity-corrected p-value `(r + 1) / (B + 1)` with MCSE.
- Minimum 30 finite replicate statistics required for an available row.
- Multi-df bootstrap F-rows supported.
- Durable `BootstrapRunPayload` JSON with target, replicate accounting, failed-
  refit policy, boundary rate, seed record.
- `auto` does not select bootstrap in schema 1.0.

**Asymptotic Wald z** (`method = asymptotic_wald_z`): labeled fallback in the
auto ladder. Low reliability, no finite-sample correction.

**None** (`method = none`): explicit no-inference path.

### 3.2 Variance-Component Inference

**Boundary LRT** (`method = boundary_lrt`): certified v1 route.
- Nested ML model comparison adding exactly one variance/covariance parameter.
- Reference distribution: Self-Liang 50:50 mixture (0.5 * point mass at 0 +
  0.5 * chi-square(1)).
- Refusals: multi-parameter boundary comparisons → `boundary_lrt_mixture_
  weights_not_certified`; callers are directed to `parametric_bootstrap_lrt`.
- Fixed-effect contrasts routed to `boundary_lrt` return `boundary_lrt_not_
  applicable_to_fixed_effects`.

**Parametric bootstrap LRT** (`stats::parametric_bootstrap_lrt`, schema
`mixedmodels.parametric_bootstrap_lrt` 1.0.0): available for multi-parameter
boundary variance comparisons; calibrates reference by simulation from the null
fit. Separate surface from fixed-effect bootstrap.

### 3.3 Profile Likelihood (LMM)

- Schema `mixedmodels.profile_likelihood_ci` version 1.0.0.
- ML fits: profiles `sigma`, `theta`, and active fixed-effect `beta` parameters.
- REML fits: profiles `sigma` and `theta` only. Beta profile intervals are
  explicitly refused until a REML-beta-profile contract is certified.
- Spline internals are Rust-owned; only the JSON payload crosses the boundary.

### 3.4 Confidence Intervals

- Profile CI: as above.
- Wald CI: trivially derivable from `std_errors`; not a Rust-certified surface.
  The R wrapper documents these as `not_certified_by_rust_inference_contract`.
- Bootstrap CI (full-model distribution): implemented as estimator-distribution
  target with percentile and basic intervals. Does not certify hypothesis-test
  p-values.

### 3.5 GLMM Inference

GLMM artifacts must expose:
- `estimation_method`: `fast_pirls_profiled`, `joint_laplace`, `joint_agq`,
  `fallback_fast_pirls`.
- `objective_definition`: `profiled_glmm_deviance`, `joint_glmm_laplace_
  deviance`, `joint_glmm_agq_deviance`.
- `response_constants`: `dropped` (fast path) or `included` (joint).
- `n_agq`: requested/effective quadrature count.

Satterthwaite, Kenward-Roger, profile-likelihood beta intervals, and REML are
LMM-only in this contract. GLMM bootstrap is a **refusal stub for all families**
in the current implementation (v1 roadmap gap).

---

## 4. FitStatus Taxonomy

```rust
pub enum FitStatus {
    ConvergedInterior,
    ConvergedBoundary,
    ConvergedReducedRank,
    ConvergedPenalised,
    NotIdentifiable,
    NotOptimized,
}
```

`ConvergedPenalised` requires explicit opt-in to a penalty method; the engine
does not silently substitute a penalized path when MLE does not exist. The
refusal/penalise decision tree is specified in the contract.

---

## 5. Documented Refusals

| Situation | Engine response |
|---|---|
| AGQ with vector or multi-term RE | `invalid_agq_request` before optimizer evaluations |
| Binomial/complete separation | `binomial_separation` diagnostic; status `NotIdentifiable` or `ConvergedPenalised` if opt-in penalty |
| KR on GLMM or ML fit | `kenward_roger_requires_reml` labeled unavailable row |
| Satterthwaite on boundary/reduced-rank fit | `satterthwaite_unavailable_at_boundary` |
| REML beta profile CI | Refused with documented reason |
| Boundary LRT on fixed effects | `boundary_lrt_not_applicable_to_fixed_effects` |
| Boundary LRT on multi-parameter comparison | `boundary_lrt_mixture_weights_not_certified` |
| `joint_laplace` without NLopt | Typed `mm_fit_error` (not silent fallback) |
| Offsets, `I(...)`, `poly(...)`, GAM smooths | Not parsed; `FormulaError` |
| `I(x^2)`, backtick identifiers, formula contrasts | Not parsed in v0 |
| GLMM bootstrap | Refusal stub for all families in current implementation |
| Cluster bootstrap p-values | `bootstrap_cluster_resample_p_value_unavailable` (estimator distribution only) |
| Conditional variance for GLMM `ranef(condVar=TRUE)` | Returns NA with `random_effect_conditional_variance_unavailable_for_glmm` |

---

## 6. Stability Labels

- **Certified / release_blocking_parity**: Gaussian LMM (ML/REML), Satterthwaite,
  Kenward-Roger (default build), Boundary LRT, profile likelihood (ML; REML
  theta/sigma), parametric bootstrap LMM, GLMM families listed in §1.2.
- **Implemented, not certified for 1.0**: `InverseGaussian`, Gaussian-GLMM
  non-identity link, GLMM bootstrap, full-model cluster-resample p-values,
  REML beta profile CI, F-scaling from pbkrtest non-unit scaling in KR multi-df
  rows.
- **Documented divergence** (not lme4 parity): fast-PIRLS/profiled-objective
  GLMM rows (`cbpp`, `contraception`, `verbagg`). These match MixedModels.jl
  `fast=true`, not lme4 joint estimation.
- **Feature-gated (nlopt)**: `joint_laplace`, `joint_agq`, NLopt-backed BOBYQA/
  NEWUOA. KR crossed/nested parity certified on default (NLopt) build only.
- **Explicitly deferred to vNext / out of v0 scope**: residual covariance
  structures (AR(1), spatial Exp, Matern), full model-lattice regularized search,
  Kenward-Roger for GLMM, full multivariate `cbind(y1,y2)~` response, GAM
  smooths, `nlmer`, `poly()`/`I()`/`splines` in formulas.

---

## 7. Engine Surfaces Not Yet Surfaced in the R Wrapper

These are capabilities the Rust engine documents or implements that the R wrapper
does not yet expose or exposes only partially:

| Engine capability | R wrapper status | Gap classification |
|---|---|---|
| `explain_model()` / `ModelExplanation` prefit semantic explanation | `explain_model()` is called internally on every `lmm()` fit (printed when `verbose >= 0`) and returns an `mm_explanation` object; no standalone user-facing pre-fit path documented | partial |
| `audit_report()` / `ModelAuditReport` structured health report | Not publicly exposed as a user-facing function; artifact is stored on `fit$artifact` | partial |
| `parameterization()` / `ParameterizationDrilldown` (theta/Lambda/parmap/VarCorr trace) | Not exposed as a public R function | in-scope-missing |
| `changes()` / `Vec<ModelStateChange>` (requested→effective model changes) | Not exposed as a public R function | in-scope-missing |
| `verify_convergence()` bounded restart/jitter/alternate-optimizer checks | Not exposed to R; result is inside `fit$artifact` | partial |
| `inference_table()` structured fixed-effect inference table (direct contract surface) | Exposed via `inference_table(fit)` in R | works |
| `random_options()` | Exposed as user-facing function | works |
| `inference_options()` | Exposed; enumerates available methods per fit | works |
| `test_random_effect()` boundary LRT | Exposed | works |
| KR opt-in via `method = "kenward_roger"` | Exposed in `contrast()`, `test_effect()`, `df_for_contrast()` | works |
| Bootstrap fixed-effect null via `method = "bootstrap"` with `bootstrap_control()` | Exposed | works |
| Profile CI via `confint(method = "profile")` | Exposed; REML beta refused as per contract | works |
| `ConvergedPenalised` / `NotIdentifiable` distinction in fit status | `fit$fit_status` stores the string; no R-level helper to inspect penalty method details | partial |
| Deterministic `maximal_feasible` / `design_compiled` compiler policy | No R-level argument to `lmm()` to select fit intent (`as_specified` vs `design_compiled` vs `exploratory`) | in-scope-missing |
| Cluster bootstrap `full_model_distribution` intervals (`confint(method="bootstrap")`) | Exposed with documented limitation that p-values are not certified | works |
| Conditional variance `ranef(condVar=TRUE)` for LMM | Exposed; routes through Rust `cond_var()` bridge | works |
| GLMM `summary(tests="coefficients")` Wald-z table | Exposed (recent commit) | works |
| `emmeans` integration (`emm_basis.mm_glmm`, `recover_data.mm_glmm`) | Exposed (recent commit) | works |
| Formula parser silent-accept bugs (`y ~ x-z` treated as addition, trailing operators, `2*x` → column "2") | Not guarded at R layer | upstream-blocked |
| Non-treatment contrast coding (Sum, Helmert, Polynomial) | Decorative labels only in Rust; only treatment coding built | upstream-blocked |
| Reference level ordering (first-appearance vs lme4 alphabetical) | Silent parity hazard in Rust; not guarded at R layer | upstream-blocked |
| Prediction intervals (`interval = "confidence"/"prediction"`) | Raises `mm_inference_unavailable`; not certified | in-scope-missing |
| Prediction SE (`se.fit = TRUE`) | Returns NA with unavailable-reason attribute | in-scope-missing |
| Partial conditioning (`re.form = ~(1|subject)`) | Raises `mm_inference_unavailable` | partial |
| GLMM bootstrap inference | Rust refusal stub; R propagates the refusal | upstream-blocked |
| `compare()` bootstrap LRT (`method = "bootstrap"`) | Exposed via `compare(method="bootstrap", nsim=...)` | works |

---

## 8. Feature-Gated Capabilities Summary

| Feature gate | What it enables |
|---|---|
| `nlopt` (default) | NLopt-backed BOBYQA/NEWUOA for LMMs; `joint_laplace` and `joint_agq` GLMM paths; KR crossed/nested parity certification |
| `--no-default-features` (TrustBQ) | Native dependency-light LMM (TrustBQ); GLMM uses native COBYLA/PatternSearch fallback; no NLopt or CMake required |
| `prima` | Optional development backend; not a required runtime dependency |

The R package bundles the vendored Rust snapshot; the feature set is determined
at package build time. The `glmm.R` wrapper explicitly refuses `method =
"joint_laplace"` when the vendored build lacks NLopt, raising a typed error
rather than silently falling back.

---

## 9. Pathology Corpus / Contract Version

Pathology-corpus TOML fixtures carry `contract_version = "v0.3"` (2026-04-29).
Separation-stratum fixtures accept `{NotIdentifiable, NotOptimized,
ConvergedPenalised}` as valid status sets; landing outside that set is a
contract regression.

---

## 10. Key Engine Ceilings for the R Wrapper

1. GLMM inference ceiling: Wald-z asymptotic rows only for GLMMs. No
   Satterthwaite, no KR, no profile-beta CI. This is a hard engine ceiling, not
   an R-wrapper gap.
2. REML beta profile CI: engine refuses. R wrapper correctly propagates refusal.
3. Multi-parameter boundary LRT: engine refuses and directs to parametric
   bootstrap LRT. R wrapper routes correctly.
4. AGQ with multiple random-effect terms: engine refuses before optimizer.
5. Non-treatment contrast coding: only treatment coding is actually built in the
   engine; other coding labels are decorative. This limits the R wrapper's
   ability to implement sum/Helmert/polynomial contrast coding even if it wanted
   to.
6. Prediction SEs and intervals: not implemented in the engine for either LMM or
   GLMM.
7. Formula terms `I(...)`, `poly(...)`, `offset(...)`, GAM smooths: not parsed
   by the engine.
8. KR on TrustBQ build: crossed/nested KR is "finite and plausible" but not
   certified against pbkrtest on the no-default-features build.
