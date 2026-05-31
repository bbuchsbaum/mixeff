# Crate Survey 3 — Engine Contract Summary

Survey date: 2026-05-31
Sources consulted:
- `docs/satterthwaite_scalar_contract.md`
- `docs/kenward_roger_contract.md`
- `docs/glmm_support_contract.md`
- `docs/mixed_model_compiler_inference_contract.md`
- `docs/random_effects_formulas.md`
- `docs/optimizer_profiles.md`
- `docs/bootstrap_fixed_effect_contract.md`
- `docs/profile_likelihood_json_contract.md`
- `docs/boundary_lrt_variance_component_contract.md`
- `docs/certified_joint_glmm_optimizer_contract.md`
- `docs/mixeff_upstream_support_report.md`
- `docs/v1_0_release_roadmap.md`

---

## 1. Supported Model Classes and Families

### LMM (Gaussian)
The engine fully supports Gaussian linear mixed models via profiled ML/REML. The
full inference stack is certified for this class:

- Fitted artifacts carry `fixed_effect_inference_table` (schema
  `mixedmodels.fixed_effect_inference_table` v1.0.0) after every LMM fit.
- The auto ladder for eligible scalar coefficient rows is:
  `auto → satterthwaite → asymptotic_wald_z → not_computed`.
- Satterthwaite is implemented for scalar fixed-effect contrasts (Gaussian LMM,
  any criterion). parity pinned in
  `tests/fixtures/compiler_contract/satterthwaite_lmer_test_parity_v1.json`.
- Kenward-Roger is implemented for explicit REML-fitted Gaussian LMM requests
  (scalar and multi-df). KR is opt-in only; it is not part of the `auto` ladder
  in schema 1.0.0. parity pinned in
  `tests/fixtures/compiler_contract/kenward_roger_pbkrtest_parity_v1.json`.
- Parametric bootstrap fixed-effect null tests are implemented for Gaussian LMM
  scalar contrasts (explicit only, not auto).
- Profile-likelihood CIs are implemented for ML fits (sigma, theta, and active
  beta). REML fits omit beta profile intervals by contract.
- Boundary LRT variance-component tests are implemented (Self-Liang 50:50
  mixture; single added parameter only).
- Random-effect conditional modes (ranef), VarCorr, and model audit with
  random-term cards are implemented.
- Fit status enum: `ConvergedInterior`, `ConvergedBoundary`,
  `ConvergedReducedRank`, `ConvergedPenalised`, `NotIdentifiable`,
  `NotOptimized`.

### GLMM
Supported families and links (certified for release blocking parity):

| Family | Links |
|---|---|
| Bernoulli | logit, probit, complementary log-log |
| Binomial (grouped) | logit, probit, complementary log-log; trial weights supported |
| Poisson | log, square-root |
| Gamma | log |

`InverseGaussian` and Gaussian-GLMM non-identity links exist in the engine but
are NOT certified for 1.0. They are experimental; parametric bootstrap
explicitly refuses them.

GLMM fitting modes (when NLopt feature is enabled):
- `fast = true` (default): profiled fast-PIRLS approximation, MixedModels.jl
  `fast=true` family. Faster but can diverge from `lme4` joint coefficients.
- `fast = false` (NLopt required): labelled joint Laplace (`n_agq <= 1`) or
  joint AGQ (`n_agq > 1`). Joint path is row-scoped — only rows that pass the
  joint objective/certificate/scorecard gate are `release_blocking_parity`.
  Current certified joint rows: `culcitalogreg` Laplace + AGQ.
  `cbpp` and `contraception` remain `documented_divergence`.
- AGQ (`n_agq > 1`) restricted to exactly one scalar random-effects term.
  Vector-valued or multiple RE terms refuse AGQ.

GLMM inference limitations (by contract):
- REML, Satterthwaite, and Kenward-Roger are LMM-only; GLMM artifacts report
  these as explicitly unsupported.
- GLMM parametric bootstrap is a refusal stub for all families (v1.0 gap,
  Phase C blocker in the roadmap).
- Profile likelihood is LMM-only in certified scope.

---

## 2. Supported Inference Methods

| Method | Scope | Status |
|---|---|---|
| `asymptotic_wald_z` | LMM and GLMM fixed effects | Implemented; auto fallback for LMM, default for GLMM |
| `satterthwaite` | Gaussian LMM scalar contrasts | Implemented; auto primary for eligible LMM rows |
| `kenward_roger` | Gaussian REML LMM (scalar + multi-df) | Implemented; opt-in only |
| `bootstrap` | Gaussian LMM scalar/multi-df (explicit) | Implemented; requires certified `fixed_effect_null` target |
| `boundary_lrt` (Self-Liang) | Single added variance component | Implemented; refuses multi-parameter and fixed-effect routes |
| Profile-likelihood CI | LMM ML (sigma, theta, beta) | Implemented; REML beta omitted |
| `not_computed` | All refusal paths | Implemented with stable reason codes |

Multi-df KR F test: implemented. Current parity tracks unscaled `KRmodcomp()`
F statistic. Rows where `F.scaling != 1` are explicitly documented as not yet
using pbkrtest's scaled F output — stored in fixture but not in row payload.

KR parity on native (`--no-default-features`) build: for crossed/nested designs
only finiteness with realistic tolerances is asserted, not certified pbkrtest
parity. The `nlopt` feature (default build) is required for certified crossed/
nested KR parity.

---

## 3. Formula / Random Effects Syntax

Supported (v0 contract):
- `(1 + x | g)` — correlated random intercept/slope
- `(1 + x || g)` — zero-correlation (diagonal) block
- `(x | g)` — implicit intercept (intercept added at materialization)
- `(0 + x | g)` — explicit intercept omission
- `(b | g1:g2)` — cell grouping
- `(b | g1/g2)` — nesting expansion → `(b|g1) + (b|g1:g2)`
- `(b | g1*g2)` — crossing expansion → `(b|g1) + (b|g2) + (b|g1:g2)` with
  `CrossingLikelyUnintended` Info diagnostic
- `(b | g1 & g2)` — legacy interaction syntax
- Categorical basis columns (treatment-coded and cell-means)
- Interaction basis columns

Explicitly not parsed in v0:
- `I(...)`, `poly(...)`, `offset(...)`, `s()`, `te()`, `bs()`, `ns()`
- Backtick identifiers
- Residual structure formula syntax (`residual = ar1(...)`)
- Multivariate response (`cbind(y1,y2) ~ ...`)

Known formula parser silent-acceptance bugs (v1.0 Phase B blockers):
- Trailing `+`/`-` parses cleanly
- `y ~ -x1` parsed as `1 + x1` (should error)
- `y ~ x1 - x2` parsed as `1 + x1 + x2` (subtraction treated as addition;
  lme4 uses `-` for term removal)
- Adjacent RE blocks without `+` parse cleanly
- `2*x1` becomes column `"2"`, error deferred to design-build time

`||` centering for numeric bases: **not yet implemented**. The v0 contract
specifies using a declared or weighted-mean reference and recording it — this is
listed as non-conforming in the formula doc Appendix A.

---

## 4. Optimizer Profiles

| Profile | Feature flag | LMM optimizer | GLMM optimizer |
|---|---|---|---|
| Default (NLopt) | `default-features = true` | NLopt BOBYQA/NEWUOA | NLopt joint Laplace/AGQ or fast-PIRLS |
| Native TrustBQ | `--no-default-features` | TrustBQ (pure Rust) | native COBYLA/PatternSearch |

PRIMA is an optional development backend, not a required runtime dependency.

TrustBQ is not a fallback — it is a first-class profile for dependency-light
downstream builds. Same model object, diagnostics, covariance summaries, and
inference surface as the NLopt path.

---

## 5. Documented Refusals (Stable Reason Codes)

Satterthwaite refusals:
`satterthwaite_varpar_deviance_unavailable`,
`satterthwaite_vcov_beta_derivative_unavailable`,
`satterthwaite_varpar_covariance_unavailable`,
`satterthwaite_boundary_derivative_unstable`,
`satterthwaite_nonpositive_contrast_variance`,
`satterthwaite_nonpositive_df_denominator`,
`satterthwaite_nonfinite_df`,
`satterthwaite_rank_deficient_contrast`,
`satterthwaite_unvalidated_against_reference`

Kenward-Roger refusals:
`kenward_roger_requires_reml`,
`kenward_roger_sigma_g_unavailable`,
`kenward_roger_sigma_not_positive_definite`,
`kenward_roger_adjusted_vcov_unavailable`,
`kenward_roger_adjusted_vcov_non_positive`,
`kenward_roger_information_singular`,
`kenward_roger_lbddf_unavailable`,
`kenward_roger_nonfinite_df`,
`kenward_roger_f_scaling_unavailable`,
`kenward_roger_unvalidated_against_pbkrtest`

Bootstrap refusals:
`bootstrap_null_target_unavailable`,
`bootstrap_null_fit_failed`,
`bootstrap_replicate_accounting_unavailable`,
`bootstrap_successful_replicates_too_few` (minimum 30),
`bootstrap_observed_statistic_nonfinite`,
`bootstrap_replicate_statistic_nonfinite`,
`bootstrap_failed_refit_policy_unavailable`,
`bootstrap_mcse_unavailable`,
`bootstrap_boundary_rate_too_high`

Boundary LRT refusals:
`boundary_lrt_requires_variance_component_comparison`,
`boundary_lrt_not_fixed_effect_method`,
`boundary_lrt_mixture_weights_not_certified` (multi-parameter → route to
parametric bootstrap LRT instead)

GLMM refusals:
`optimizer_nonconvergence`, `invalid_agq_request`,
`pirls_failure`, `boundary_parameter`,
`near_unit_random_effect_correlation`, `binomial_separation`

---

## 6. Stability Labels (Engine View)

| Capability | Engine status |
|---|---|
| Gaussian LMM Satterthwaite scalar | Certified; parity-fixture backed |
| KR scalar + multi-df, nlopt build | Certified; parity-fixture backed |
| KR crossed/nested, nlopt build | Certified (Penicillin, Pastes) |
| KR crossed/nested, native build | Finite + plausible only; NOT certified |
| KR F scaling (non-unit `F.scaling`) | Not yet in row payload; stored in fixture only |
| GLMM fast-PIRLS (cbpp, contraception) | `documented_divergence`; NOT `lme4` parity |
| GLMM joint Laplace (culcitalogreg) | `release_blocking_parity` |
| GLMM joint AGQ (culcitalogreg) | `release_blocking_parity` |
| GLMM bootstrap | Refusal stub; all families |
| InverseGaussian / Gaussian-GLMM | Experimental; not certified |
| Profile-likelihood CI LMM ML | Certified |
| Profile-likelihood CI LMM REML beta | Explicitly omitted by contract |
| Profile-likelihood CI GLMM | Not in scope |
| Multivariate `cbind(y1,y2)~` | Post-1.0; API-breaking |
| AR(1) / spatial residual covariance | Designed but not implemented |
| `I()`, `poly()`, formula transforms | Not parsed in v0 |
| `||` numeric centering | Contract specified; implementation missing |
| GLMM parametric bootstrap | Phase C blocker; not yet delivered |

---

## 7. Engine Surfaces the R Wrapper Has Not Yet Consumed (Wrapper-Side Gaps)

Based on comparing engine contracts against `R/inference.R`, `R/revive.R`,
`R/glmm.R`, and `R/reporting.R`:

1. **KR multi-df F scaling**: The engine stores scaled F statistic and
   `F.scaling` in the fixture but the row payload currently carries the unscaled
   statistic. The R wrapper does not distinguish or warn about this. When
   `F.scaling != 1`, the row p-value is computed from the unscaled F — a
   known discrepancy from pbkrtest. The wrapper does not surface this caveat.

2. **`||` centering note**: The engine contract requires recording and
   emitting the centering reference value for numeric `||` terms as a
   `FormulaCanonicalized` diagnostic. The implementation is missing in the
   engine; the R wrapper also has no path to surface this when it lands.

3. **GLMM `fast=false` fallback labeling**: The engine emits
   `fallback_status = fallback_fast_pirls` when a joint optimizer attempt
   returns to the fast path. The R wrapper's `glmm()` function surfaces
   `estimation_method` and `objective_definition` fields but there is no R-side
   test or print method that distinctly warns users when `fallback_status` is
   present — the approximation-gap distinction between a certified joint path
   and a fallback-PIRLS path may not be visible.

4. **GLMM inference refusal not forwarded**: When a user calls
   `inference_table()` or `test_effect()` on a GLMM fit, the engine contract
   says GLMM artifacts must report Satterthwaite/KR/REML as explicitly
   unsupported. It is unclear from the current R code whether the `mm_glmm`
   inference dispatch routes correctly to an explicit unsupported reason vs.
   silently returning empty/NA rows.

5. **Bootstrap `cluster_resample` vs. `fixed_effect_null` distinction**: The
   engine contract distinguishes three bootstrap targets and only
   `fixed_effect_null` can produce a p-value. The R wrapper has a
   `cluster_bootstrap` method stub that returns unavailable tables — correct —
   but documentation and inference-options display may not clearly tell users
   why `cluster_bootstrap` produces no p-value while `bootstrap` can.

6. **`CoefTable` / `ModelSummary` Satterthwaite surfacing**: The v1.0 roadmap
   (Phase C item 14) notes that `CoefTable` and `ModelSummary` never surface
   Satterthwaite/KR rows — clients must use `FixedEffectInferenceTable` via
   `inference_table()`. The R wrapper's `summary.mm_lmm` path passes through
   the inference table correctly for the default `method = "auto"` case, but
   the engine roadmap notes this as a downstream-misleading gap at the Rust
   level.

7. **Wald CI (`wald_confint`)**: The engine roadmap lists `wald_confint(model,
   level)` as a Phase C item not yet delivered. The R wrapper's
   `confint.mm_lmm` has `method = "wald"` which falls through to Wald CI
   construction using `std_error` from the fit object — but this path is
   R-side math rather than a Rust-certified Wald CI endpoint with a stable
   schema.

8. **Parametric bootstrap LRT (`parametric_bootstrap_lrt`)**: The engine exposes
   `stats::parametric_bootstrap_lrt` for variance-component comparisons that
   have more than one added parameter (where boundary LRT refuses). The R
   wrapper's `test_random_effect()` routes multi-parameter refusals to a
   reason code, but there is no `method = "parametric_bootstrap_lrt"` dispatch
   path in the R layer that calls this engine endpoint.

---

## 8. Feature-Gated Capabilities

| Capability | Feature flag | Effect |
|---|---|---|
| NLopt BOBYQA/NEWUOA (LMM) | `nlopt` (default) | Faster/more accurate LMM optimization |
| Joint Laplace / joint AGQ (GLMM) | `nlopt` (default) | Only certified joint GLMM paths require this |
| Native TrustBQ (LMM) | `--no-default-features` | Dependency-light; same inference surface |
| PRIMA optimizer | `prima` feature (optional) | Development backend; not required at runtime |
| Unstable internals | `unstable-internals` | For benchmarks/probes; not public API |

The CRAN-compatible build of the `mixeff` R package must use
`--no-default-features` (NLopt excluded). This means: on a CRAN build,
`fast = false` / joint Laplace / joint AGQ paths are unavailable and must be
refused explicitly. The R wrapper's `glmm()` function already notes this: the
`joint_laplace` method is guarded and refused if the NLopt backend is absent.
However, certified crossed/nested KR parity also degrades on the native build
(finiteness only, not pbkrtest-certified) — the R wrapper does not currently
surface this distinction.
