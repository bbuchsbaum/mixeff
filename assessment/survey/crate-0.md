# Crate Survey: mixeff-rs Engine Capability Inventory

Survey date: 2026-05-31
Sources: docs/guide/05_what_is_supported.md, docs/guide/04_when_the_crate_refuses.md,
         docs/guide/01_getting_started.md, docs/glmm_support_contract.md,
         docs/mixed_model_compiler_inference_contract.md, docs/profile_likelihood_json_contract.md
R wrapper sources: R/glmm.R, R/inference.R, R/inference-options.R, R/predict.R, R/compare.R, NAMESPACE

---

## Model Classes (Engine)

| Class | Rust type | Engine status |
|---|---|---|
| Linear mixed model | `LinearMixedModel` / `LinearMixedModelBuilder` | Stable |
| Generalized linear mixed model | `GeneralizedLinearMixedModel` | Stable for certified family/link matrix |

Both are surfaced in the R wrapper (`lmm()`, `glmm()`).

---

## GLMM Families and Links

### Engine-certified 1.0 surface

| Family | Links | Default |
|---|---|---|
| Bernoulli | Logit, Probit, Cloglog | Logit |
| Binomial | Logit, Probit, Cloglog | Logit |
| Poisson | Log, Sqrt | Log |
| Gamma | Log, Inverse | Inverse (engine default) |

### Engine-implemented but NOT certified for 1.0

| Family | Status |
|---|---|
| InverseGaussian | Implemented; not validated to cross-language parity standard; bootstrap explicitly refused |
| Normal (non-identity link, as GLMM) | Log, Inverse, Sqrt links implemented; not certified; bootstrap explicitly refused |

### Normal + Identity

Explicitly refused by engine at construction time with actionable message: "use LinearMixedModel".

### R wrapper family/link gate (`mm_glmm_supported_family_links`)

The R wrapper imposes a stricter gate than the engine exposes:

| Family | Wrapper-allowed links |
|---|---|
| binomial | logit, probit, cloglog |
| poisson | log, sqrt |
| Gamma | **log only** |

**Wrapper gap — Gamma/Inverse link**: The engine certifies `Gamma` with the `Inverse` link as a
1.0 surface (it is the canonical default). The R wrapper blocks it (`mm_glmm_supported_family_links`
allows only `"log"` for `Gamma`). This is a wrapper-side restriction, not an engine limitation.

**Wrapper gap — InverseGaussian/Normal-GLMM**: Not exposed in the R wrapper at all. This matches
the engine's own "not certified" designation — an intentional exclusion from the R certified surface.

---

## Formula DSL

| Construct | Engine status | R wrapper status |
|---|---|---|
| Additive fixed effects `y ~ x1 + x2` | Stable | Surfaced |
| Main effects + interaction `x1 * x2` | Stable | Surfaced |
| Interaction only `x1 : x2` | Stable | Surfaced |
| Nesting `x1 / x2` | Stable | Surfaced |
| Explicit intercept `0+`, `-1`, `1+` | Stable | Surfaced |
| Correlated RE `(re \| g)` | Stable | Surfaced |
| Zero-correlation RE `(re \|\| g)` | Stable | Surfaced |
| Interaction grouping `(re \| g1 & g2)` | Stable | Surfaced |
| Cell-level grouping `(re \| g1:g2)` | Stable | Surfaced |
| Nested grouping `(re \| g1/g2)` | Stable | Surfaced |
| Main+cell grouping `(re \| g1*g2)` | Stable | Surfaced |
| Minimal stateless `I()` subset | Stable | Surfaced |
| Full `I()` / model.matrix transforms | Out of scope | Not surfaced (correct) |

---

## Estimation Paths

| Path | Engine status | Feature gate | R wrapper status |
|---|---|---|---|
| LMM profiled (RE)ML via blocked-Cholesky PLS; auto-dispatched optimizer | Stable | None (default build) | Surfaced via `lmm()` |
| GLMM profiled fast-PIRLS (`fast=true`) | Stable, labelled | None | Surfaced via `glmm(method="pirls_profiled")` |
| GLMM joint Laplace (`fast=false`, `n_agq <= 1`) | Stable, labelled | **`nlopt` feature required** | `glmm(method="joint_laplace")` is wired but raises `mm_fit_error` because the vendored Rust is compiled without nlopt |
| GLMM AGQ (`fast=false`, `n_agq > 1`) | Stable, labelled | **`nlopt` feature required** | Same refusal path |

### Optimizer auto-dispatch (LMM)

The engine selects among: `PatternSearch`, `TrustBq`, `NloptBobyqa`, `NloptNewuoa`, `Cobyla`,
or `PrimaBobyqa` (if `prima` feature). Callers do not choose. The chosen optimizer and convergence
outcome are recoverable from `MixedModelFit::opt_summary`. The R wrapper exposes the `artifact`
field which carries this metadata.

---

## Inference and Post-Fit Summaries

### LMM inference (engine)

| Surface | Engine status |
|---|---|
| Point estimates (`coef`, `vcov`, `stderror`) | Stable |
| Random effects, fitted, logLik, AIC, BIC | Stable |
| Variance components (`VarCorr`) | Stable |
| Model summary (markdown/HTML/LaTeX) | Stable |
| Wald CIs (`CoefTable::wald_confint`) | Stable |
| Satterthwaite df (scalar contracts, crossing/nested fixture-expanding) | Stable for Gaussian REML LMMs |
| Kenward-Roger df (scalar-test scope) | Stable for scalar tests; beyond scalar is out of scope |
| Profile-likelihood CIs — sigma, theta, ML beta | Stable for LMM; REML omits beta intervals by contract |
| Profile-likelihood CIs — GLMM | Out of scope |
| Parametric bootstrap — LMM | Stable |
| Parametric bootstrap — GLMM (Bernoulli, Binomial, Poisson, Gamma) | Stable |
| Parametric bootstrap — InverseGaussian and Normal-as-GLMM | **Explicitly refused** |
| LRT (`LikelihoodRatioTest`, `ModelComparisonTable`) | Stable, typed taxonomy with reason codes |
| Boundary LRT (`BoundaryLikelihoodRatioTest`, Self-Liang mixture) | Stable |

### LMM inference — R wrapper coverage

| Surface | R wrapper status |
|---|---|
| `coef`, `fixef`, `vcov`, `stderror` | Surfaced |
| `ranef`, `fitted`, `logLik`, `AIC`, `BIC` | Surfaced |
| `VarCorr` | Surfaced |
| `summary()` / Wald-z table | Surfaced |
| `contrast()` / Satterthwaite / Kenward-Roger | Surfaced |
| `test_effect()` — Satterthwaite, KR, bootstrap, bootstrap_lrt | Surfaced |
| `test_random_effect()` — boundary LRT | Surfaced |
| `estimability()` | Surfaced |
| `df_for_contrast()` | Surfaced |
| `confint()` — Wald, bootstrap, profile | Surfaced |
| Profile CI (REML omits beta by contract) | Surfaced with labeled refusal |
| `parametric_bootstrap()` | Surfaced (LMM model-comparison bootstrap) |
| `compare()` | Surfaced |
| `inference_options()` | Surfaced |
| `inference_table()` | Surfaced |
| `anova()`, `drop1()` | Surfaced |

### GLMM inference — R wrapper coverage

| Surface | R wrapper status | Gap? |
|---|---|---|
| `summary(object, tests="coefficients")` Wald-z table | Surfaced | None |
| `coef`, `fixef`, `fitted`, `residuals`, `logLik`, `AIC`, `BIC` | Surfaced | None |
| `ranef`, `VarCorr` | Surfaced | None |
| `confint()` for GLMM | **Not registered** — only `confint.mm_lmm` is in NAMESPACE | Wrapper gap |
| `test_effect()` for GLMM | **Not registered** — only `test_effect.mm_lmm` is in NAMESPACE | Wrapper gap |
| `test_random_effect()` for GLMM | Not registered | Wrapper gap |
| `inference_options()` for GLMM | Not registered | Wrapper gap |
| `inference_table()` for GLMM | Not registered | Wrapper gap |
| `predict.mm_glmm` | Registered but **raises `mm_inference_unavailable` unconditionally** | Wrapper gap |
| GLMM parametric bootstrap (LRT comparison) | `parametric_bootstrap()` accepts only `mm_lmm` | Wrapper gap |
| Profile CI for GLMM | Out of scope by engine design | Out of scope (correct) |

---

## Refusal Taxonomy (Engine)

The engine returns typed errors or typed inference refusals — never fabricated numbers.

### Construction-time refusals

| Code | Trigger |
|---|---|
| `UnsupportedFamilyLink` | Family/link pair outside the supported matrix |
| `InvalidArgument` | Normal + Identity passed to GLMM (actionable: use LMM) |
| `ConstantResponse` | Response column is constant |
| `NoRandomEffects` | Formula has no random-effects terms |
| `RankSaturatedFixedEffects` | Fixed-effect design is rank-saturated |
| `RankDeficient` | Fixed-effect design is rank-deficient |

Each code is stable across releases (SemVer-covered after 1.0) via `MixedModelError::code()`.

### Inference-time refusals (typed, not thrown)

| Situation | Mechanism |
|---|---|
| Variance-component boundary test | `BoundaryLikelihoodRatioTest` + `BoundaryLrtStatus`; chi-sq mixture is named |
| Model comparison | `ModelComparisonClass` + `ModelComparisonReasonCode` |
| Profile CI — unidentifiable side | Typed refusal; no spline extrapolation |
| AGQ with vector-valued RE | `invalid_agq_request` diagnostic + refusal before optimizer evaluations |
| GLMM separation-like behavior | `binomial_separation` diagnostic; `ConvergedPenalised` or `NotIdentifiable` |

### GLMM diagnostic codes (stable)

`optimizer_nonconvergence`, `invalid_agq_request`, `pirls_failure`, `boundary_parameter`,
`near_unit_random_effect_correlation`, `binomial_separation`.

### Five distinguishable GLMM failure modes

The engine contract guarantees that optimizer failure, approximation gap (fast-PIRLS vs. joint),
weak identification, response-constant convention difference, and separation-like behavior are
each separately readable from the artifact — they must not collapse into a single "did not converge".

---

## Fit Status Taxonomy

```
FitStatus: ConvergedInterior | ConvergedBoundary | ConvergedReducedRank |
           ConvergedPenalised | NotIdentifiable | NotOptimized
```

The R wrapper surfaces `fit_status` as a character field on the fit object and checks it in
`is_singular()`, `mm_boundary_df_method_unavailable()`, and `inference_options()`.

---

## Stability Labels

| Label | Meaning |
|---|---|
| Stable | Stable surface, parity-tested, SemVer-covered after 1.0 |
| Stable, labelled | Stable but approximation/path is explicit in fit metadata |
| Refused | Typed error or typed inference refusal |
| Out of scope | Deferred to 2.0 |

---

## Explicitly Out of Scope (Engine — 2.0 candidates)

- Multivariate response `cbind(y1, y2) ~`
- Profile-likelihood CIs for GLMMs
- Parametric bootstrap for InverseGaussian and Normal-as-GLMM
- Full `I()` / arbitrary formula transforms beyond the minimal stateless subset
- First-class polars/arrow ingestion
- Kenward-Roger beyond the current scalar-test scope
- AR(1)/spatial residual covariance structures (designed but not in stable surface)
- Confirmatory vs. regularized mode switch (`random_strategy` argument) — designed in contract but not yet in R API

---

## Engine Capabilities Not Yet Surfaced in R Wrapper (Wrapper-Side Gaps)

1. **Gamma/Inverse link for GLMM** — engine certifies it; R wrapper blocks it.
2. **`predict.mm_glmm`** — engine has no certified `predict_new` contract for GLMM yet; R wrapper hard-refuses. This is likely correct but should be labeled as `upstream-blocked` pending a GLMM predict contract.
3. **GLMM fixed-effect inference** (`test_effect`, `contrast`, `confint`, `inference_options`, `inference_table`) — only wired for `mm_lmm`; no `mm_glmm` methods registered. The engine GLMM artifact carries Wald-z via `summary()` but the full inference dispatch (Satterthwaite, KR, bootstrap) is LMM-only by engine contract. Asymptotic Wald-z for GLMM contrasts is engine-available but not wired on the R side.
4. **GLMM parametric bootstrap** — `parametricbootstrap_glmm` exists in the engine (stable for Bernoulli, Binomial, Poisson, Gamma); `parametric_bootstrap()` in R only accepts `mm_lmm`.
5. **`joint_laplace` / AGQ paths** — wired in glmm.R but always refused because the vendored build lacks the `nlopt` feature. Users see a typed `mm_fit_error` rather than a silent fallback.
6. **`n_agq > 1` on `pirls_profiled`** — accepted by the R validator but the engine only uses AGQ metadata for the joint path; behavior on the profiled path with `n_agq > 1` is not documented in the R layer.
7. **Conditional variance for random effects** — `mm_lmm_cond_var_json` FFI exists in lib.rs; R-side wiring is present in the bridge but the public `ranef(..., condVar = TRUE)` surface may be partial (not confirmed as complete).
8. **REML vs. ML profile CI beta** — engine omits beta intervals under REML by explicit contract; R wrapper handles this with a labeled refusal row. Correct behavior.

---

## GLMM Approximation Semantics Note

`fast=true` (default, `pirls_profiled`) is a profiled fast-PIRLS approximation —
**not** the same as `lme4::glmer`'s joint Laplace fit. Objective values differ
by a response-constant convention (`dropped` vs `included`). The engine records
this in `estimation_method`, `objective_definition`, and `response_constants`
fields so the R layer can explain the divergence without guessing.

Documented-divergence rows (cbpp, contraception, verbagg) are deliberate 1.0
release exclusions from `release_blocking_parity` — they track MixedModels.jl
`fast=true` behavior, not lme4 joint-estimation parity.
