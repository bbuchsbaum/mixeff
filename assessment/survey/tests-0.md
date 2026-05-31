# Test Coverage Survey: tests-0
**Files surveyed:** `tests/testthat/test-lmm.R`, `tests/testthat/test-lme4-parity.R`,
`tests/testthat/test-lmerTest-parity.R`, `tests/testthat/helper-lme4-parity.R`

---

## test-lmm.R

### What is tested

| # | Test name | Assertions |
|---|-----------|------------|
| 1 | `lmm() fits an LMM and stores flat extractor fields` | S3 classes, model_frame names, schema_name, `fixef` names, `theta` length, `sigma` finiteness, `logLik` finiteness, `nobs`, `df.residual` |
| 2 | `print.mm_lmm exposes artifact provenance and audit entry points` | Output contains schema string, crate version, audit verb list |
| 3 | `lmm() auto-prints explain_model unless verbose is -1` | Presence/absence of "Random effects" in output at default vs. `verbose=-1` |
| 4 | `standard extractors return stored fit quantities` | `logLik` class+attrs, `AIC`/`BIC`/`deviance` identity with stored fields, `formula`, `model.frame`, `fitted == predict`, `residuals`, `predict(re.form=NA)`, `predict(se.fit=TRUE)` returns named list with all-NA SEs and `mm_unavailable_reason` attr |
| 5 | `lmm() accepts positive case weights and preserves them for inference` | `fit$weights` preserved, `refit` inherits weights, `vcov` status="available", `diag(vcov) == std_errors^2`, contrast `std_error` matches `sqrt(vcov["x","x"])` |
| 6 | `lmm() rejects invalid case weights` | Wrong-length weights error class `mm_data_error`, non-positive weights error class `mm_data_error` |
| 7 | `random-effect and variance-component extractors are shaped like lme4 basics` | `ranef` class, names, column names, rownames; `coef` class, names, intercept identity (ranef+fixef); `VarCorr` class, table rows, `residual_sd == sigma` |
| 8 | `revived extractor paths return typed values` | `ranef(condVar=TRUE)` postVar array dimensions/finiteness; `predict(newdata=df)` reproduces `fitted`; `predict(re.form=~(1|subject))` raises `mm_inference_unavailable` |

### Tolerances / assertions used
- `expect_equal` with `tolerance = 1e-8` for residual/fitted arithmetic identity
- `expect_equal` with `tolerance = 1e-12` for `diag(vcov) == std_errors^2`
- `expect_identical` / `expect_s3_class` / `expect_named` for structural checks
- Error class checks (`mm_data_error`, `mm_inference_unavailable`) via `expect_error(class=...)`

### Skips
- None (`skip_on_cran`, `skip_if_not_installed` absent from this file)

### What is NOT tested in this file
- ML (`REML = FALSE`) fits — all fixtures use default `REML = TRUE`
- Singular/boundary fits (`is_singular`, `theta` clamped to zero)
- `vcov` parity against lme4 (only self-consistency checked)
- `getME` beyond `"theta"` (e.g. `"Zt"`, `"Lambda"`, `"X"`)
- `model.matrix` as a standalone extractor (only exercised inside parity helper)
- `contrast()` with non-asymptotic methods (`satterthwaite`, `kenward_roger`)
- `test_effect()` / `anova()` (not called at all)
- `confint()` / `profile()` (deferred; not called)
- `simulate()` / `refit()` beyond trivial weight-preservation smoke test
- `revive()` round-trip (only condVar/predict tested; no JSON save-reload)
- `audit()` / `diagnostics()` / `inference_table()` / `model_report()` — only their names appear in print output
- Multiple random-effects groups (all fixtures use a single grouping factor)
- Random slope (`(1 + x | g)`) structure
- Crossed random effects (`(1|g1) + (1|g2)`)
- `allow_new_levels` prediction path
- `parametric_bootstrap` / `bootstrap_control`
- `test_random_effect`
- Numeric response edge cases (all-zero response, single observation per group)

---

## test-lme4-parity.R

### What is tested

| # | Test name | What it asserts |
|---|-----------|-----------------|
| 1 | Mismatch ledger is shipped and well-formed | `expected-mismatches.json` exists, schema_version=1, required keys in every entry, valid status values, abs/rel bounds present for mismatch entries, no duplicate (case_id, field) pairs |
| 2 | Upstream mismatch report template is shipped | `upstream-mismatch-report-template.md` exists and contains 14 required section strings |
| 3 | Classic lme4 parity fixture manifest is valid | manifest schema_version, ≥6 cases, unique IDs, required keys in each case |
| 4 | Classic lme4 LMM cases match core extractors within documented tolerances | Via `mm_expect_core_lme4_parity`: `fixef`, `sigma`, `logLik`, `deviance`, `AIC`, `BIC`, `fitted`, `residuals`, `nobs`, `df.residual`, `model.matrix`, `VarCorr` (std_dev + residual_sd), `theta` (when labels align) — for all LMM cases in the fixture |
| 5 | Classic lme4 LMM cases match random-effect modes | Via `mm_expect_ranef_lme4_parity`: `ranef` values at common group/level keys |
| 6 | Classic lme4 LMM cases match supported prediction semantics | Via `mm_expect_prediction_lme4_parity`: conditional `predict`, `predict(re.form=NA)`, `fitted==predict`, `se.fit=TRUE` structure, `newdata` conditional/population, `re.form=~0`, error on `interval="confidence"` |

**Fixture cases exercised (LMM filter):**
- `sleepstudy_random_intercept` (RI, REML)
- `sleepstudy_random_intercept_slope` (RS, REML)
- `dyestuff_random_intercept` (balanced RI, REML)
- `dyestuff2_singular_random_intercept` (`known_boundary`, REML)
- `penicillin_crossed_intercepts` (crossed RI, REML)
- `pastes_two_intercepts` (two RI groups, REML)
- `cake_recipe_temp` (categorical×continuous interaction + RI, REML)

**Default tolerances from fixture:**
`fixef=1e-5`, `sigma=1e-4`, `logLik=1e-4`, `deviance=1e-4`, `AIC=1e-4`, `BIC=1e-4`, `fitted=1e-4`, `residuals=1e-4`, `theta=1e-3`, `varcorr=1e-3`, `ranef=1e-4`, `model_matrix=0` (exact)

GLMM cases (`cbpp_binomial_logit_profiled_pirls`, `grouseticks_poisson_log_profiled_pirls`) are present in the fixture but filtered out by the `model="lmm"` default in `mm_lme4_parity_cases()`; they are not exercised by these three tests.

### Skips
- Each parity test block calls `mm_skip_if_no_lme4()` → `skip_if_not_installed("lme4")` — all three parity test blocks skip on CRAN or any environment without lme4
- Individual cases skip if dataset is unavailable (`mm_lme4_case_data` emits a skip when the dataset is not found)
- `ranef` parity silently skips cases with no common group label keys

### What is NOT tested in this file
- `vcov` parity against lme4 (helper computes lme4 vcov internally for `asymptotic` reference but never asserts `vcov(fit)` ≈ `vcov(ref)` as a named check)
- `coef` parity (not compared to lme4)
- `getME` beyond `"theta"` (no `"Zt"`, `"Lambda"`, `"Z"`, `"b"` comparison)
- ML (`REML=FALSE`) LMM cases — all seven LMM fixture cases have `"reml": true`
- GLMM parity (two GLMM fixture cases are filtered out; no test exercises them here)
- `is_singular` return value verified against lme4 `isSingular` for the `dyestuff2` known-boundary case
- `ranef(condVar=TRUE)` postVar parity vs lme4 (condVar is tested in `test-lmm.R` but not compared to lme4 values)
- `predict` on genuinely new groups (`allow_new_levels=TRUE`)
- `predict` on unseen continuous covariate values outside training range (extrapolation smoke test)
- `contrast()` / `test_effect()` / `anova()` parity against lme4 extractors (those are in `test-lmerTest-parity.R`)
- `confint` parity (deferred by design per PRD §3; not present)
- Ledger entries for `dyestuff2` boundary case — no assertion that `is_singular` or boundary warning is surfaced correctly

---

## test-lmerTest-parity.R

### What is tested

| # | Test name | Scope |
|---|-----------|-------|
| 1 | Scalar fixed-effect contrasts match lmerTest references | `contrast()` estimate, std_error, df, statistic, p_value, statistic_name for methods `satterthwaite`/`kenward_roger`/`asymptotic` on 5 cases (sleepstudy RI, sleepstudy RS, dyestuff, penicillin, cake) |
| 2 | Single-df term rows match lmerTest ANOVA F equivalents | `test_effect()` and `anova()` term rows for `"Days"`, methods `satterthwaite`/`kenward_roger`, on sleepstudy RI + RS |
| 3 | Multi-df Kenward-Roger term rows match lmerTest ANOVA | `test_effect()` and `anova()` for terms `"recipe"` and `"recipe:temp"` with `kenward_roger` on `cake_recipe_temp` |

**Tolerances used (in `mm_expect_scalar_lmerTest_parity`):**
`estimate=1e-6`, `std_error=1e-3`, `df=1e-2`, `statistic=1e-3`, `p_value=1e-3`

**Tolerances used (in `mm_expect_term_lmerTest_parity`):**
`DenDF` tolerance `1e-2`, F-equivalent `1e-3`, p-value `1e-3`; effective rank exact

### Skips
- All three tests call `mm_skip_if_no_lmerTest()` → `skip_if_not_installed("lmerTest")`
- KR tests additionally call `mm_skip_if_no_pbkrtest()` → `skip_if_not_installed("pbkrtest")`
- No `skip_on_cran` — these tests run on CRAN if lmerTest/pbkrtest are available (they are Suggests)

### What is NOT tested in this file
- `anova()` type I and type II (only type III exercised via lmerTest)
- Multi-df Satterthwaite ANOVA (only KR is tested for multi-df; `dyestuff`/`penicillin`/`pastes` intercept-only models are not exercised for ANOVA term rows)
- `pastes_two_intercepts` inference (not in the lmerTest case list)
- `dyestuff2_singular_random_intercept` inference (not in case list; boundary behavior under test_effect is untested)
- `asymptotic` method for multi-df (F-equivalent) terms (only scalar contrast exercises asymptotic)
- Non-zero `rhs` values other than the artificially shifted RHS used by the helper
- `contrast()` with a user-supplied contrast matrix (multiple rows / multi-df direct contrast)
- `mm_lincomb` parity against lmerTest
- `df_for_contrast` directly
- KR detail metadata for non-scalar (multi-row) contrasts
- Intercept-only models under lmerTest (dyestuff, penicillin — only in scalar contrast block, not in ANOVA term block)

---

## helper-lme4-parity.R

### What this file provides (not a test file itself)

- Path resolution for fixture JSON (`mm_lme4_parity_manifest_path`, `mm_lme4_parity_manifest`)
- Case filtering (`mm_lme4_parity_cases`)
- Skip helpers: `mm_skip_if_no_lme4`, `mm_skip_if_no_lmerTest`, `mm_skip_if_no_pbkrtest`
- Reference version collector (`mm_reference_versions`)
- Data/formula/tolerance accessors for fixture cases
- Fit-pair constructors: `mm_lme4_fit_pair`, `mm_lmerTest_fit_pair`
- Numeric payload extraction and close-comparison utilities (`mm_numeric_payload`, `mm_expect_numeric_close`, `mm_expect_p_value_close`)
- Reference computation helpers: `mm_lmerTest_reference_contrast`, `mm_lmerTest_anova_table`, `mm_lmerTest_anova_row`, `mm_mixeff_term_row`
- Full parity expectation functions: `mm_expect_core_lme4_parity`, `mm_expect_varcorr_lme4_parity`, `mm_expect_ranef_lme4_parity`, `mm_expect_prediction_lme4_parity`, `mm_expect_scalar_lmerTest_parity`, `mm_expect_term_lmerTest_parity`

**Weaknesses in the helper itself:**
- `mm_expect_prediction_lme4_parity` tests `predict(re.form=~(1|subject))` raises an error in `test-lmm.R` but does NOT verify the same error is raised for each parity case (formula-form `re.form` is silently skipped)
- `mm_expect_ranef_lme4_parity` silently skips (not fails) when no common group keys are found — a naming mismatch would be invisible
- `mm_lmerTest_scalar_rhs` shifts the RHS by one SE to ensure the test is non-trivial, but the exact shift is arbitrary and not documented in assertions
- `mm_assert_parity` delegates to `mm_scoreboard_record` which accumulates global state; test isolation is not enforced (reset is manual via `mm_scoreboard_reset`)

---

## Cross-cutting gaps

1. **No ML (REML=FALSE) LMM fixture case** — all seven LMM fixture cases use `"reml": true`; logLik/AIC comparison semantics differ under ML
2. **vcov parity not asserted** — `stats::vcov(fit)` is only checked for self-consistency; no comparison against `stats::vcov(ref)` from lme4
3. **GLMM parity tests absent from these four files** — two GLMM fixture cases exist in the manifest but are not exercised by any test in this group
4. **Singular/boundary case (`dyestuff2`) has no behavioral assertions** — included in the LMM loop but `is_singular`, warnings, and structured diagnostics for boundary fits are not verified
5. **`coef()` parity not checked** — no comparison of `coef(fit)` vs `coef(ref)` from lme4
6. **`ranef(condVar=TRUE)` postVar values not compared to lme4** — only shape/finiteness verified
7. **New-level prediction path (`allow_new_levels`) untested**
8. **`confint` smoke test absent** — not even an error check that it raises `mm_inference_unavailable`
9. **`test_effect` / `anova` not exercised for intercept-only models** (dyestuff, penicillin, pastes)
10. **Asymptotic method not exercised for multi-df ANOVA terms**
11. **`parametric_bootstrap` / `test_random_effect` have zero coverage in these files**
