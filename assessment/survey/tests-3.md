# Test Survey: tests-3 — confint/profile, boundary-LRT, boundary-singular

**Survey date:** 2026-05-31  
**Files covered:**
- `tests/testthat/test-confint-profile.R`
- `tests/testthat/test-boundary-lrt.R`
- `tests/testthat/test-boundary-singular.R`

**Source files consulted:**
- `R/inference.R` (confint, boundary_lrt, test_random_effect implementations)
- `R/methods-extract.R` (VarCorr, is_singular, mm_varcorr_boundary_flag)
- `R/methods-print.R` (print.mm_lmm, print.mm_varcorr)
- `NAMESPACE` (exported S3 methods)
- `planning/PRD.md` §9.4, §9.5.6
- `tests/testthat/_snaps/boundary-singular.md` (existing snapshots)

---

## File 1: test-confint-profile.R

### What is tested

| # | Test name | What it asserts |
|---|-----------|-----------------|
| 1 | `confint(method='profile') under ML returns beta/sigma/theta rows` | `mm_confint` class; `attr(ci, "method") == "profile_likelihood"`; `attr(ci, "fit_criterion") == "ML"`; `(Intercept)` and `Days` in rownames; `sigma` in rownames; payload `$table` contains `parameter_kind` values `beta`, `sigma`, `theta`; beta rows have finite, monotonic (lower ≤ estimate ≤ upper) bounds; VC rows with finite bounds bracket the estimate. |
| 2 | `confint(method='profile') under REML omits beta with reason_code` | `fit_criterion == "REML"`; beta rows in payload carry `reason_code == "profile_beta_unavailable_under_reml"`; beta lower/upper are `NA`; non-beta rows carry no reason_code. |
| 3 | `profile CI for beta agrees with Wald CI on well-behaved ML fit` | Relative difference `< 5 %` on `(Intercept)` and `Days` for both lower and upper bounds. |
| 4 | `profile CI parm subsetting filters returned rows` | `confint(fit, parm = "(Intercept)")` returns exactly `rownames == "(Intercept)"` and that row is a subset of the full-CI matrix. |
| 5 | `profile CI surfaces a typed refusal on a boundary singular fit` | On Dyestuff2 boundary fit: either a structured error with class in `{mm_inference_unavailable, mm_schema_error, mm_bridge_error, mm_fit_error}`, or a payload where every NA-bound row carries a non-empty `regularity` note. |

### Dataset / fixtures

- `sleepstudy` (lme4). Loaded with `skip_if_not_installed("lme4")` + `skip("sleepstudy dataset unavailable")` fallback.
- `Dyestuff2` (lme4). Same guards.

### Skip conditions

- All tests: `skip_if_not_installed("lme4")` via `mm_skip_if_no_lme4_local()` (called inside `mm_sleepstudy_data()`).
- Test 5 explicitly calls `testthat::skip_if_not_installed("lme4")` a second time.
- No `skip_on_cran()` guards.

### Tolerances / assertion style

- Monotonicity: `lower ≤ estimate + 1e-8`, `estimate ≤ upper + 1e-8`.
- Wald/profile agreement: relative `< 0.05` (5 %).
- No numeric parity against lme4 profile CI values.

### What is NOT tested (gaps)

1. **Non-default confidence level.** Only `level = 0.95` is exercised. A level of `0.90` or `0.99` changes the cutoff passed to the FFI and the column-name formatting (`"5.0 %"` / `"95.0 %"`) — both untested.
2. **`boundary_clamped_lower` flag propagation.** The table column `boundary_clamped_lower` is set by `mm_translate_profile_row` but no test checks that it is `TRUE` on a theta/sigma row that the upstream engine clamped to zero.
3. **Theta row naming.** Theta rows are named `theta1`, `theta2`, … by `mm_map_profile_parameter`. No test verifies these row names appear in `rownames(ci)` or that `parm = "theta1"` subsets correctly.
4. **`parm` subsetting for sigma / theta by name.** Test 4 only checks `(Intercept)` (a beta). Subsetting by `"sigma"` or `"theta1"` is untested.
5. **`parm` with names not in the fit.** The filter `table$parameter %in% parm` would silently return zero rows for a misspelled parameter name. No test for that error path or empty-result handling.
6. **`print.mm_confint` output.** The `print` S3 method is exported but no test checks that the printed matrix or its summary line is well-formed.
7. **Empty intervals list from payload.** When `payload$intervals` is empty, `mm_empty_profile_table()` is returned. No test exercises this path (e.g., a near-trivial model with no profile support).
8. **Schema negotiation failure.** No test exercises `mm_json_negotiate` raising `mm_schema_error` when the crate returns an unrecognised schema version.
9. **`confint` is not registered for `mm_glmm`.** `NAMESPACE` only exports `confint.mm_lmm`. There is no test asserting that `confint(glmm_fit, method="profile")` raises an informative error or routes correctly — the gap exists at the registration level as well as the test level.

---

## File 2: test-boundary-lrt.R

### What is tested

| # | Test name | What it asserts |
|---|-----------|-----------------|
| 1 | `test_random_effect() exposes certified one-component boundary LRT` | `mm_random_effect_test` class; `status == "available"`; `method == "boundary_lrt_self_liang_mixture"`; `statistic_name == "chi_bar_square"`; finite statistic and p-value in `[0,1]`; `ordinary_chisq_dof == 1L`; `reference_distribution == "0.5 * chi-square(0) + 0.5 * chi-square(1)"`; `details$mixture` has two components each with `weight == 0.5`; Self and Liang citation present; `print(out)` contains `"chi-square(0)"`; `reporting_table(out)` has correct `reference_distribution`. |
| 2 | `boundary_lrt is refused on fixed-effect inference surfaces` | `test_effect(fit, "x", method = "boundary_lrt")` returns `mm_effect_test` with `status == "unsupported"` and `reason_code == "boundary_lrt_not_applicable_to_fixed_effects"`; same for `contrast(fit, c(0,1), method = "boundary_lrt")`. |
| 3 | `boundary_lrt refuses multi-parameter random-effect geometry` | On a `(1 + x | subject)` fit (3 theta parameters): `status == "not_assessed"`, `reason_code == "boundary_lrt_mixture_weights_not_certified"`, `p_value` is NA, reason message contains `"certifies only one added boundary"`, `theta_parameters == 3L`. |

### Dataset / fixtures

All three tests use `mk_boundary_lrt_data()` — a fully synthetic dataset generated with `set.seed(617L)`, no external dependency. No `skip_if_not_installed` / `skip_on_cran` guards.

### Tolerances / assertion style

No numeric tolerances — the test is structural: correct classes, status codes, reason_code strings, mixture weights exactly `c(0.5, 0.5)`.

### What is NOT tested (gaps)

1. **REML input to `test_random_effect`.** The implementation re-fits under ML if needed (via `mm_boundary_lrt_ml_fit`), but no test verifies that a REML-fitted input is accepted (and silently re-fitted) rather than erroring, nor that the re-fit path produces valid output.
2. **`refit_for_comparison = FALSE`.** The `test_random_effect.mm_lmm` signature accepts a `refit_for_comparison` argument; the existing test only uses the default. Edge case: when the full fit is already ML but passing `FALSE` skips the refit.
3. **Large p-value / zero-statistic case.** No test checks the degenerate case where the reduced model has a *larger* log-likelihood than the full model (chi-bar-square statistic = 0, p = 1.0). The p-value clamp/floor behavior is untested.
4. **Two-component model with certified mixture weights.** The multi-component case is tested only for refusal. No positive test exists for a geometry where the crate does certify mixture weights for 2+ theta parameters (when/if that capability is added).
5. **`test_random_effect` on a model with multiple random grouping factors.** Crossed or nested multi-group models (e.g., `(1|subject) + (1|item)`) are not covered; the tested term lookup and reduced-formula construction are untested in that setting.
6. **`test_random_effect` on a non-existent term name.** No test checks that requesting an unknown group name raises an informative error rather than silently returning garbage.
7. **`reporting_table(out)` beyond `reference_distribution` column.** Only `reference_distribution` is spot-checked; other table columns (statistic, p_value, method, reason) are unchecked in `reporting_table` output.
8. **`print(out)` full snapshot.** `print(out)` is only checked for the substring `"chi-square(0)"` via `expect_match`. No snapshot test locks down the full print layout.

---

## File 3: test-boundary-singular.R

### What is tested

| # | Test name | What it asserts |
|---|-----------|-----------------|
| 1 | `VarCorr table carries a boundary flag on near-zero std_devs` | `mm_varcorr` class; `"boundary"` column present in `vc$table`; `any(vc$table$boundary) == TRUE` on Dyestuff2. |
| 2 | `snapshot: print(VarCorr(fit)) tags boundary components with [boundary]` | Snapshot of `print(VarCorr(fit))` output. Locked to existing snap in `_snaps/boundary-singular.md`: shows `[boundary]` tag in table row and explanation footer. |
| 3 | `snapshot: print(fit) on singular fit names rank and points to audit verbs` | On `is_singular(fit) == TRUE`: matches `"covariance matrix is rank-deficient"`, `"Use changes(fit)"`, `"Use random_options(spec, group"` in printed output; asserts absence of forbidden phrases (`"Try (1 | "`, `"Drop the random slope"`, `"suggested starting model"`, `"we recommend"`, `"you should"`, `"try .* instead"`); snapshot of `cat(printed)`. If `is_singular(fit)` is FALSE, test is skipped. |

### Dataset / fixtures

All tests use `mk_boundary_fit()`: Dyestuff2 from lme4, REML fit of `Yield ~ 1 + (1|Batch)`. Guards: `skip_if_not_installed("lme4")`, `skip("Dyestuff2 dataset is unavailable")`. Test 3 has an additional conditional skip if `is_singular(fit)` is FALSE.

### Skip conditions

- `skip_if_not_installed("lme4")` in `mk_boundary_fit()` (affects all three tests).
- No `skip_on_cran()`.
- Test 3 conditionally skips if boundary convergence did not occur on this build.

### Tolerances / assertion style

- Structural: class checks, column presence, boolean predicates.
- Two `expect_snapshot` calls (tests 2 and 3) — sensitive to exact print output.

### What is NOT tested (gaps)

1. **`is_singular()` on a non-singular (full-rank) fit.** No test verifies that `is_singular()` returns `FALSE` for a well-conditioned model. The only exercised path is the boundary case.
2. **`is_singular()` return type and class.** No test checks the type (`logical(1)`) or that it does not return `NA` on a successfully converged model.
3. **`VarCorr` on a non-boundary fit.** `vc$table$boundary` should be all-FALSE for a full-rank fit. This is untested; a regression could incorrectly flag non-zero components.
4. **`VarCorr` with multiple random groups.** Only a single-group intercept model is tested. The boundary-flag logic for multi-group models (e.g., `(1|subject) + (1|item)`) — including that only truly near-zero components are flagged — is untested.
5. **`VarCorr` with a correlated intercept+slope boundary.** When a correlation parameter hits ±1 (a different boundary condition than std_dev = 0), `mm_varcorr_boundary_flag` behaviour is untested.
6. **`VarCorr.mm_glmm`.** `NAMESPACE` registers `VarCorr.mm_glmm` (aliased to `VarCorr.mm_lmm`). No test verifies this method on a GLMM boundary fit.
7. **`print(VarCorr(fit))` on a non-singular fit.** No snapshot or structural test verifies that the `[boundary]` annotation is *absent* when none of the std_devs are near zero — regression protection for spurious tagging is missing.
8. **The `regularity` note in `print(fit)` for non-boundary rank-deficient fits with multiple reduced dimensions.** The snapshot covers rank 1 reduced to rank 0 (`r0: requested rank 1; fitted effective rank 0`). A model reduced from rank 3 to rank 1 (two unsupported dimensions) is untested.
9. **Forbidden-phrase coverage for `print(VarCorr)`.** The advice-creep check (R9 contract) is applied to `print(fit)` but not to `print(VarCorr(fit))`. The VarCorr output could independently introduce forbidden recommendation language without being caught.
10. **`changes(fit)` output for a boundary/reduced-rank fit.** The singular print test checks that `"Use changes(fit)"` appears, but no test calls `changes(fit)` and verifies it returns structured diagnostic content (not just the instruction to call it).

---

## Cross-cutting observations

- **No `skip_on_cran()` in any of the three files.** All tests depend on lme4/Dyestuff2 or on FFI calls to the compiled crate. If these files are included in CRAN submission without `skip_on_cran()`, they will run on CRAN infrastructure and may fail or time out on profile CI.
- **Confidence level hardcoded at 0.95 throughout test-confint-profile.R.** The column-name formatter `sprintf("%.1f %%", ...)` and the FFI `as.numeric(level)` argument are both untested at any other level.
- **No numeric parity against lme4 for profile CI.** The `lme4-inference-confint-anova.md` survey documents that lme4's profile CIs are the gold standard for VC parameters. mixeff has a Wald/profile agreement test but no lme4 parity check; this is a test-gap (not an out-of-scope issue) since the bead's done condition (item c) specifies agreement within stated tolerance.
- **`boundary_clamped_lower` is structurally present but never asserted.** This flag is populated by `mm_translate_profile_row` and could silently regress to always-FALSE without any test catching it.
