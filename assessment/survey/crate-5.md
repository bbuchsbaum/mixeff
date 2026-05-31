# Crate Survey 5: Bootstrap Fixed-Effect Contract + Boundary LRT Contract

Sources read:
- `docs/bootstrap_fixed_effect_contract.md`
- `docs/boundary_lrt_variance_component_contract.md`
- `docs/mixed_model_compiler_inference_contract.md` (architecture context)
- `docs/glmm_support_contract.md`
- `docs/satterthwaite_scalar_contract.md`
- `docs/kenward_roger_contract.md`
- `docs/profile_likelihood_json_contract.md`
- `R/inference.R`, `R/glmm.R`, `R/predict.R`, `NAMESPACE`

---

## 1. Engine Capabilities Summary

### Model Classes

| Class | Status |
|---|---|
| Gaussian LMM (REML/ML) | Fully supported; primary target |
| GLMM Bernoulli / logit | Certified for 1.0 |
| GLMM Binomial / logit (incl. grouped) | Certified for 1.0 |
| GLMM Bernoulli+Binomial / probit, cloglog | Certified for 1.0 |
| GLMM Poisson / log, sqrt | Certified for 1.0 |
| GLMM Gamma / log | Certified for 1.0 |
| GLMM InverseGaussian / Gaussian-GLMM non-identity | Implemented but NOT certified for 1.0; experimental |

### Inference Methods (LMM)

| Method | Engine Status |
|---|---|
| Satterthwaite (scalar + multi-df) | Implemented; auto-ladder: `auto -> satterthwaite -> asymptotic_wald_z` |
| Kenward-Roger (scalar + multi-df) | Implemented; opt-in only for schema 1.0.0 (`method = "kenward_roger"`); certified on `nlopt` build for crossed/nested; native-optimizer build only asserts finiteness (not pbkrtest parity) |
| Bootstrap fixed-effect null (certified `fixed_effect_null` target) | Implemented via `fixed_effect_null_bootstrap_inference_table()`; explicit `method = bootstrap` only; never auto-selected in schema 1.0.0 |
| Bootstrap full-model distribution (estimator intervals, NOT p-values) | Implemented; produces percentile/basic CIs but does not certify hypothesis-test p-values |
| Bootstrap cluster resample (estimator distribution, NOT p-values) | Implemented; no p-value certification in schema 1.0.0 |
| Bootstrap parametric LRT (`stats::parametric_bootstrap_lrt`) | Implemented; separate surface from fixed-effect inference rows |
| Boundary LRT (variance component, Self-Liang mixture) | Implemented; v1 certified route: exactly 1 variance/covariance parameter added, nested ML comparison |
| Asymptotic Wald z (fallback) | Implemented; low reliability; auto fallback when Satterthwaite prerequisites fail |
| Profile likelihood CIs | Implemented; ML fits expose sigma/theta/beta; REML fits expose sigma/theta only (beta profile not certified for REML) |

### Inference Methods (GLMM)

| Method | Engine Status |
|---|---|
| Asymptotic Wald z | Supported |
| Satterthwaite | Explicitly unsupported for GLMM; contract states LMM-only |
| Kenward-Roger | Explicitly unsupported for GLMM; contract states LMM-only |
| REML | Explicitly unsupported for GLMM |
| Profile CI (GLMM) | Explicitly refused for InverseGaussian/Normal GLMM paths |

### Documented Refusals (stable)

- `boundary_lrt_requires_variance_component_comparison` / `boundary_lrt_not_fixed_effect_method`: Boundary LRT is not a fixed-effect p-value method.
- `boundary_lrt_mixture_weights_not_certified`: Multi-parameter boundary comparisons are refused; callers are pointed at `stats::parametric_bootstrap_lrt`.
- `invalid_agq_request`: AGQ with `n_agq > 1` is refused for vector-valued or multi-term random effects.
- `bootstrap_null_target_unavailable` and 8 other `bootstrap_*` codes: Bootstrap rows are unavailable when certified null target cannot be constructed, too few replicates, etc.
- `kenward_roger_requires_reml` and 9 other `kenward_roger_*` codes: KR refuses ML fits or missing adjusted-covariance artifacts.
- `satterthwaite_*` family (8 codes): Satterthwaite refuses when derivatives, vcov_varpar, or df prerequisites fail.
- GLMM Satterthwaite/KR: artifacts must report LMM-only methods as unsupported.
- REML beta profile CI: explicitly omitted until a REML beta-profile contract is certified.

### Feature-Gated Capabilities

- `nlopt` feature: Required for `method = "joint_laplace"` / `method = "joint_agq"` GLMM path. Default build enables NLopt. Dependency-light `--no-default-features` uses native COBYLA/PatternSearch for GLMM.
- KR crossed/nested parity: Only certified on the `nlopt` build; `--no-default-features` build gives finite rows but not pbkrtest-certified parity for crossed/nested models.
- `prima` backend: Optional development backend; not a required runtime dependency.

### Reliability Labels

Bootstrap: `low` or `moderate` only (never `high` in schema 1.0.0); `not_available` for prerequisites failure.
Satterthwaite: `low`, `moderate`, `not_available`.
KR: same pattern.

### Bootstrap Targets (stable wire labels)

| `target_kind` wire label | Certifies p-value? |
|---|---|
| `fixed_effect_null` | Yes (explicit `method = bootstrap` only) |
| `full_model_distribution` | No |
| `cluster_resample` | No |

### Structured Bootstrap Row Fields (v1 contract)

Every available bootstrap row carries `details.bootstrap` with:
`mcse`, `requested_replicates`, `completed_replicates`, `successful_replicates`, `failed_refits`, `failed_refit_policy`, `boundary_count`, `boundary_rate`, `seed_record` (`seed`, `seed_rng`), null-target summary. Human-readable `notes` remain, but R should prefer structured fields for programmatic decisions.

---

## 2. Wrapper-Side Gaps (engine exposes; R may not surface fully)

### GAP-1: Structured `details.bootstrap` fields partially surfaced

The contract requires that R prefer structured `details.bootstrap` fields (mcse, boundary_rate, seed_record, null-target summary, replicate accounting) over prose notes for programmatic decisions. The R wrapper does read `seed_record$seed`, `boundary_rate`, and `mcse` for the `mm_bootstrap_confint_summary()` print helper, but there is no dedicated accessor that returns the full `details.bootstrap` structured payload from a `mm_contrast` or `mm_effect_test` row to the user. The `details` list column is accessible via `x$table$details`, but there is no documented helper and no test asserting the structured fields are present. Severity: minor — the data is reachable, just not ergonomically surfaced.

### GAP-2: `details.bootstrap.target_kind` / null-target summary not exposed as a documented slot

The contract pins `target_kind` (`fixed_effect_null` / `full_model_distribution`) as a stable wire label that downstream code must use to distinguish certified p-value runs from distribution runs. The R wrapper stores this in `details` (list column), but it does not document it and has no accessor. A consumer could accidentally treat a `full_model_distribution` run as a p-value source. Severity: minor.

### GAP-3: Multi-parameter boundary LRT path points to `stats::parametric_bootstrap_lrt` but that function is not exported by the R wrapper

The engine contract says multi-parameter boundary comparisons should be routed to `stats::parametric_bootstrap_lrt` (schema `mixedmodels.parametric_bootstrap_lrt`, version 1.0.0). The R wrapper has `bootstrap_lrt` as a `method` in `test_effect()` and an internal `compare()` path using `parametric_bootstrap_lrt` as a label, but there is no exported top-level function named `parametric_bootstrap_lrt` and no documented user-facing route for the multi-parameter case. The engine's refusal message points callers at a function that doesn't exist in the public R API. Severity: major — the gap between refusal message and available user surface is confusing.

### GAP-4: GLMM inference methods `test_effect` / `contrast` / `confint` not exported for `mm_glmm`

The NAMESPACE exports `test_effect.mm_lmm`, `contrast.mm_lmm`, `confint.mm_lmm`, `estimability.mm_lmm`, and `df_for_contrast.mm_lmm`, but none of these have `mm_glmm` dispatch. The engine contract for GLMM explicitly limits inference to asymptotic Wald z, but that method is valid and the wrapper only exposes it through `inference_table` and `mm_lincomb.mm_glmm`. A user who calls `contrast(glmm_fit, ...)` will get a method-not-found error. Severity: major — Wald-z contrasts are certified by the engine for GLMM but not surfaced through the primary contrast/test_effect interface.

### GAP-5: KR `F.scaling` not yet applied in multi-df rows

The KR contract explicitly notes that multi-df rows with `F.scaling != 1` from pbkrtest remain documented as using the unscaled F statistic (`f_scaling = 1.0`) until scaled-F support is promoted. The R wrapper faithfully passes through whatever Rust returns, but users comparing against `pbkrtest::KRmodcomp()` on models with non-unit F scaling will see different statistics. This is a documented engine-side limitation, not a wrapper bug. Severity: minor — documented divergence.

### GAP-6: Profile CI for REML beta not surfaced

The profile CI contract says REML fits omit fixed-effect beta profile intervals. The R wrapper correctly routes through `mm_profile_confint` and defers to Rust, but there is no user-facing message explaining why beta intervals are absent for REML fits. A user calling `confint(reml_fit, method = "profile")` will get sigma/theta CIs only with no diagnostic explaining the omission. Severity: minor.

### GAP-7: `ranef(condVar = TRUE)` for LMM returns unavailable-reason stub

The R wrapper returns `mm_attach_ranef_postvar_unavailable` for both LMM and GLMM when `condVar = TRUE`. The engine contract does not explicitly certify conditional variance output in these survey docs, so this is consistent with the engine's current certification state. However, the test file `tests/testthat/test-ranef-condvar.R` exists in the untracked files list, suggesting this is being actively worked. Severity: test-gap.

### GAP-8: `joint_laplace` / AGQ GLMM path blocked in current R build

The glmm.R wrapper documents that `method = "joint_laplace"` requires the optional upstream `nlopt` backend and raises `mm_fit_error` when unavailable. This is correct and intentional. However the `nAGQ > 1` AGQ path is checked only after method validation; the R wrapper gates `method = "joint_laplace"` but the vendored Rust build may already gate this at the C boundary. Severity: works — the refusal is correct and tested.
