# Survey: tests-7 — Tutorial-Derived Integration & Structural Test Files

Survey date: 2026-05-31
Files covered:
- `tests/testthat/test-brown-2021-lme-tutorial.R`
- `tests/testthat/test-bw-lme-tutorial.R`
- `tests/testthat/test-iamciera-lme4tutorial.R`
- `tests/testthat/test-sdamr-lmm-companion.R`
- `tests/testthat/test-codingclub-mixed-models.R`
- `tests/testthat/test-aphantasia-reproduction.R`
- `tests/testthat/test-pw2-crossed-nested.R`
- `tests/testthat/test-factor-contrast-semantics.R`

---

## 1. test-brown-2021-lme-tutorial.R

### What is tested
Five LMM models from Brown (2021) "An Introduction to Linear Mixed-Effects Modeling in R" (OSF project v6qag):

| Case | Formula pattern | REML | Slow? |
|------|----------------|------|-------|
| `figure_random_intercepts` | `yvar ~ xvar + (1\|PID)` | TRUE | no |
| `figure_random_slopes` | `yvar ~ xvar + (1+xvar\|PID)` | TRUE | no |
| `rt_modality_full` | `RT ~ modality + (1+modality\|PID) + (1+modality\|stim)` | TRUE | yes |
| `rt_modality_reduced` | intercept-only fixed, two correlated-slope RE groups | TRUE | yes |
| `rt_interaction` | two-way interaction, correlated and uncorrelated slopes, crossed grouping | TRUE | yes |

For fast (non-slow) cases, parity is asserted on: `fixef`, `sigma`, `logLik`, `AIC`, `BIC`, `fitted`, `residuals`, `VarCorr` (diagonal SDs + residual SD).

Slow cases run under `MIXEFF_RUN_SLOW_PARITY=true` and assert the same quantities on the large real-data RT dataset (21,679 rows, 53 subjects, 543 stimuli). A separate test asserts LRT statistic and p-value parity for the full vs. reduced comparison (also slow-gated).

### Tolerances
- `figure_*` cases: fixef 1e-3, scalar 1e-2, fitted 1e-2, varcorr 5e-2
- `rt_modality_full`: fixef 2e-4, scalar 2e-3, fitted 2e-3, varcorr 3e-1
- `rt_modality_reduced`: fixef 1, scalar 2e-3, fitted 1, varcorr 3e-1 (very loose — near-degenerate RE)
- `rt_interaction`: fixef 5e-4, scalar 5e-3, fitted 5e-3, varcorr 2e-2

### Skip conditions
- `mm_skip_if_no_lme4()` (`skip_if_not_installed("lme4")`) on every model-fitting test
- Skip on missing fixture file via `brown_fixture_path()` (calls `skip()`)
- Slow cases: `skip_if_not(MIXEFF_RUN_SLOW_PARITY == "true")`
- No `skip_on_cran()` present

### What is NOT tested
- `ranef()` BLUP parity against lme4 for any Brown case
- `nobs()`, `df.residual()` extractors
- `deviance()` parity
- `model.matrix()` / fixed design matrix parity
- `theta` (Cholesky factor) parity
- `predict()` on new data (out-of-sample rows, re.form=NA)
- LRT cases with non-slow models (only slow LRT test exists)
- `compare()` table structure assertions beyond LRT stat and p-value (e.g. AIC column in comparison table)
- Boundary/singular fit handling for `rt_modality_reduced` (loose tolerances mask near-singular variance — no diagnostic assertion that the singular flag is surfaced)
- Error path: formula with grouping factor absent from data
- `print()` / `summary()` smoke test for these tutorial fits

---

## 2. test-bw-lme-tutorial.R

### What is tested
Seven LMM models from Bodo Winter's "A Very Basic Tutorial for Performing Linear Mixed Effects Analyses" using the `bw_politeness_data.csv` (83 complete rows, 6 subjects, 7 scenarios):

| Case | Formula | REML |
|------|---------|------|
| `random_intercepts_attitude` | `frequency ~ attitude + (1\|subject) + (1\|scenario)` | TRUE |
| `random_intercepts_attitude_gender` | adds `gender` | TRUE |
| `ml_null_gender` | ML, gender only | FALSE |
| `ml_full_attitude_gender` | ML, attitude + gender | FALSE |
| `ml_interaction` | ML, `attitude * gender` | FALSE |
| `random_slopes_attitude` | ML, correlated slopes for attitude | FALSE |
| `random_slopes_null_gender` | ML, null model with random slopes | FALSE |

Parity fields: `fixef`, `sigma`, `logLik`, `AIC`, `BIC`, `fitted`, `residuals`, `VarCorr`.

LRT test: three model comparisons (`attitude`, `interaction`, `random_slope_attitude`) asserting LRT statistic (tol 1e-3) and p-value (tol 1e-4) against `lme4::anova()`.

### Tolerances
fixef 2e-5 to 2e-2; fitted 1e-3 to 1e-2; varcorr 1e-2 to 2e-2 (tightest suite in this group).

### Skip conditions
- `mm_skip_if_no_lme4()` on fitting tests
- Missing fixture via `bw_politeness_path()` → `skip()`
- No `skip_on_cran()`, no slow gate

### What is NOT tested
- `ranef()` BLUP parity
- `nobs()`, `df.residual()`, `deviance()`, `model.matrix()`, `theta` parity
- `predict()` on new data
- Missing data path: raw 84-row dataset includes 1 NA frequency; only `complete = TRUE` path is used in fitting; incomplete-data behavior of `lmm()` itself is not exercised
- Contrast coding: attitude is an unordered factor — the direction of the treatment contrast is never asserted (only the magnitude)
- `print()` / `summary()` smoke tests
- `compare()` with REML=TRUE models (LRT section only compares ML models, as correct; no test verifying that comparing REML models raises a diagnostic warning/error)

---

## 3. test-iamciera-lme4tutorial.R

### What is tested
Five LMM models from the iamciera/lme4tutorial GitHub repository using stomata density data (727 rows, 75 `il` lines, 16 trays, 10 rows, 5 columns):

| Case | Formula | Slow? |
|------|---------|-------|
| `drop_il` | `trans_abs_stom ~ 1 + (1\|tray) + (1\|row) + (1\|col)` | no |
| `max_model` | `~ il + (1\|tray) + (1\|row) + (1\|col)` | yes |
| `drop_col` | `~ il + (1\|tray) + (1\|row)` | yes |
| `drop_row` | `~ il + (1\|tray) + (1\|col)` | yes |
| `drop_tray` | `~ il + (1\|row) + (1\|col)` | yes |

All 4 slow cases run under `MIXEFF_RUN_SLOW_PARITY=true`. A separate slow-gated test runs the backwards-selection LRT comparisons (drop_col, drop_row, drop_tray, drop_il vs. max_model), asserting the ML-refit flag (`mm_cmp$refit`), LRT statistic and p-value.

Parity fields: `fixef`, `sigma`, `logLik`, `AIC`, `BIC`, `fitted`, `residuals`, `VarCorr`.

### Tolerances
Uniform 1e-4 fixef, 1e-4 scalar, 1e-4 to 2e-4 fitted, 1e-3 to 3e-3 varcorr.

### Skip conditions
- `mm_skip_if_no_lme4()` on all model tests
- Missing fixture → `skip()`
- Slow gate: `MIXEFF_RUN_SLOW_PARITY=true`

### What is NOT tested
- `ranef()` BLUP parity (75-level `il` factor — meaningful check)
- `nobs()`, `df.residual()`, `deviance()`, `model.matrix()`, `theta`
- The `drop_il` intercept-only model: only 1 fixed effect, no slope — `compare()` when one model has no continuous slope is not exercised in the non-slow path
- `trans_epi_count` column is loaded but never used in any model — the epi_count response is an undocumented dead weight in the fixture
- `predict()` with new data
- Error path: high-dimensional `il` factor (75 levels) passed to `lmm()` with a fixed-slope formula — potential aliasing/rank deficiency is not probed

---

## 4. test-sdamr-lmm-companion.R

### What is tested
Five LMM models from Spekenbrink's "Statistics: Data Analysis and Modelling R companion":

| Case | Data | Formula pattern | Special |
|------|------|----------------|---------|
| `anchoring_random_intercept` | anchoring (4632 rows, 31 referrers) | `everest_feet ~ anchor + (1\|referrer)` | custom sum contrast on anchor |
| `anchoring_random_slope` | same | `~ anchor + (1+anchor\|referrer)` | correlated slope; relative varcorr comparison |
| `anchoring_uncorrelated_numeric_slope` | same (numeric `anchor_contrast`) | `~ anchor_contrast + (1\|referrer) + (0+anchor_contrast\|referrer)` | uncorrelated; loose tol |
| `speeddate_maximal_crossed` | speeddate (1509 rows, 102 iid/pid) | maximal crossed, 3-way interaction | slow |
| `speeddate_uncorrelated_crossed` | same | double-bar `\|\|` uncorrelated RE | slow |

The suite also tests a manual LRT for the anchoring random-slope vs. uncorrelated model.

Contrast semantics: `anchoring_random_intercept` and `anchoring_random_slope` use `contrasts(dat$anchor) <- c(1/2, -1/2)` (sum-like half-coding), and `sdamr_expected_fixef()` / `sdamr_expected_varcorr()` re-derive the expected values in lme4's parameterisation to compare against mixeff's cell-means or treatment output.

Parity fields: `fixef`, `sigma`, `logLik`, `AIC`, `BIC`, `fitted`, `residuals`, `VarCorr` (with mode="relative" for near-zero variance components).

### Tolerances
Very mixed: `anchoring_uncorrelated_numeric_slope` fixef tol=0.5, fitted tol=10 (near-zero variance component makes numeric comparison fragile). Slow speeddate cases: 1e-4 tight.

### Skip conditions
- `mm_skip_if_no_lme4()`
- Missing fixture → `skip()`
- `MIXEFF_RUN_SLOW_PARITY=true` for speeddate models

### What is NOT tested
- `ranef()` BLUP parity for any case
- `nobs()`, `df.residual()`, `deviance()`, `model.matrix()`, `theta`
- `predict()` on new data
- The manual LRT tolerance for the anchoring case is only 1e-1 for statistic and 1e-3 for p-value — the assertion is very weak for the statistic; a tighter check is absent
- `compare()` for the slow crossed models (only anchoring LRT uses `compare()`)
- Double-bar `||` uncorrelated RE in speeddate: correctness of VarCorr off-diagonal being zero is not explicitly asserted (only diagonal SDs are checked)
- The `anchor_contrast` numeric column is used only in the uncorrelated case but is never verified to equal the expected linear contrast values; its derivation in the data-loading helper is untested
- Incomplete speeddate data: 53 NA rows dropped; `lmm()` behavior on a data frame with NAs (before calling `na.omit()`) is not probed

---

## 5. test-codingclub-mixed-models.R

### What is tested
Seven LMM models from the Coding Club "Introduction to Linear Mixed Models" tutorial using the dragons dataset (480 rows, 8 mountain ranges, 3 sites, 24 mountain:site samples):

| Case | Formula | REML |
|------|---------|------|
| `first_random_intercept` | `testScore ~ bodyLength2 + (1\|mountainRange)` | TRUE |
| `crossed_site_warning_example` | adds `(1\|site)` (not truly nested) | TRUE |
| `nested_sample_explicit` | `(1\|mountainRange) + (1\|sample)` where sample = mountainRange:site | TRUE |
| `nested_slash_equivalent` | `(1\|mountainRange/site)` | TRUE |
| `nested_random_slope` | `(1+bodyLength2\|mountainRange/site)` | TRUE |
| `ml_full_body_length` | ML version of nested_sample_explicit | FALSE |
| `ml_reduced_intercept` | ML intercept-only | FALSE |

LRT test: full vs. reduced ML models asserting LRT stat (tol 1e-4) and p-value (tol 1e-4) via `compare()`.

Parity fields: `fixef`, `sigma`, `logLik`, `AIC`, `BIC`, `fitted`, `residuals`, `VarCorr`.

None of the cases are slow-gated.

### Tolerances
`nested_sample_explicit` and `nested_slash_equivalent` varcorr tol 1e-2 to 3e-2, others tighter.

### Skip conditions
- `mm_skip_if_no_lme4()`
- Missing fixture → `skip()`

### What is NOT tested
- `ranef()` parity — especially relevant for the nested models where BLUP structure differs
- `nobs()`, `df.residual()`, `deviance()`, `model.matrix()`, `theta`
- `predict()` on new data
- Assertion that `(1|mountainRange/site)` and `(1|mountainRange) + (1|sample)` produce the same `fixef` and `logLik` — both formulas are run separately against lme4 but never compared to each other directly
- The `crossed_site_warning_example` case is labelled "warning_example" because lme4 emits an ambiguity warning when `site` is used without nesting notation; whether mixeff surfaces an equivalent diagnostic is not tested
- VarCorr off-diagonal correlations are never asserted (only diagonal SDs)

---

## 6. test-aphantasia-reproduction.R

### What is tested
Full multi-model reproduction from a real aphantasia masked-image recognition experiment (anonymised). The fixture lives in `inst/extdata/aphantasia/` and contains `reference.json`, `trials.rds`, `metadata.rds`.

**Fixture shape test** (always runs, no skip):
- Schema name, trial count, metadata row count, participant ID anonymisation (hex format), presence of core model IDs, S1 count = 5, primary/combined trial counts.

**Core fit-side reproduction** (`MIXEFF_RUN_APHANTASIA=true`):
Models: `primary`, `sensitivity`, `intact`, `combined` (GLMM, binomial), `rt` (LMM), `S7_age_covariate`, `S9_age_matched_subset`, `S9_age_matched_subset_age_covariate`.
Asserts: `fixef` absolute tolerance (from `reference.json`), `logLik` relative tolerance, `AIC` relative tolerance, all vs. cached lme4 references.

**S1 random-effects stability** (`MIXEFF_RUN_APHANTASIA_STRESS=true`):
Five supplemental models (`S1_intercept_only`, `S1_current_uncorrelated_slopes`, `S1_correlated_slopes`, `S1_item_mask_slope`, `S1_maximal`). Same fixef/logLik/AIC assertions.

**GLMM inference vcov gate** (`MIXEFF_RUN_APHANTASIA=true` + runtime check that `vcov` is available):
Asserts `mm_lincomb()` DiD estimates at centered SOA and 25ms SOA match the reference within 2e-2 absolute. The gate (`aphantasia_has_glmm_full_vcov`) dynamically skips if the vcov payload is unavailable.

### Tolerances
From `reference.json$tolerances`: separate `fixef_abs`, `logLik_rel`, `AIC_rel` for LMM and GLMM models.

### Skip conditions
- `skip_on_cran()` on all three model-fitting tests
- `mm_skip_if_no_lme4()`
- `MIXEFF_RUN_APHANTASIA=true` for core and inference tests
- `MIXEFF_RUN_APHANTASIA_STRESS=true` for S1 stress tests
- Runtime vcov gate for inference test

### What is NOT tested
- `VarCorr` parity (no assertion on random-effect variance components against lme4 reference)
- `sigma` parity for the RT LMM model
- `fitted` / `residuals` parity for any aphantasia model
- `predict()` on new data for either the LMM or GLMM aphantasia models
- The inference test only checks two linear combinations (centered-SOA DiD, 25ms DiD); the full `primary_dd` reference table may have additional rows not being compared
- No test that the vcov-unavailable path returns a well-formed `mm_inference_unavailable` condition (only the available-vcov path is exercised when the gate passes)
- `BIC` is not asserted (only logLik and AIC)
- The sensitivity, intact, combined, and S7/S9 models are asserted on fixef/logLik/AIC but not on inference quantities

---

## 7. test-pw2-crossed-nested.R

### What is tested
Two LMM models from Patrick Ward's "Crossed vs Nested Random Effects" tutorial using a 60-row simulated sports dataset (6 teams, 37 players, some players on multiple teams):

| Case | Formula | REML |
|------|---------|------|
| `crossed` | `player_value ~ 1 + (1\|team) + (1\|player_id)` | TRUE |
| `nested` | `player_value ~ 1 + (1\|team) + (1\|team:player_id)` | TRUE |

Parity fields: `fixef`, `sigma`, `logLik`, `AIC`, `fitted`, `VarCorr` (note: `BIC` and `residuals` are absent from `pw2` parity assertions, unlike other tutorial suites).

Structural test: asserts that `ranef(fit)` exposes groups named `team` and `team:player_id` for the nested model.

### Tolerances
fixef 1e-4, scalar 1e-3, fitted 5e-3, varcorr 5e-2.

### Skip conditions
- `mm_skip_if_no_lme4()`
- Missing fixture → `skip()`

### What is NOT tested
- `BIC` parity (explicitly absent — the only tutorial file missing it)
- `residuals` parity (also absent from the parity loop)
- `ranef()` BLUP value parity (only group-name structure is checked, not BLUP magnitudes)
- `nobs()`, `df.residual()`, `deviance()`, `model.matrix()`, `theta`
- `predict()` on new data
- LRT / `compare()` between the crossed and nested models (the tutorial's primary pedagogical point — distinguishing crossed from nested — is not tested numerically)
- Whether the fixture satisfies the partial-crossing property is tested structurally (players on >1 team), but not quantitatively (e.g. the exact number of multi-team players)
- Error path: passing `(1|player_id:team)` (reversed interaction) vs. `(1|team:player_id)` — same grouping but potentially different label handling

---

## 8. test-factor-contrast-semantics.R

### What is tested
One test: a purely synthetic dataset (`factor_contrast_design()`, 18 groups × 2 factor levels × 4 reps = 144 rows, seed 8) with a Helmert-like half-unit contrast (`contrasts(dat$f) <- matrix(c(-0.5, 0.5), ...)`).

Tests the no-intercept random-slope formula `y ~ f + (0 + f | g)`:
- `VarCorr(fit)$table$name` contains `"f: a"` and `"f: b"` (cell-means coding)
- `ranef(fit)$g` column names contain `"f: a"` and `"f: b"`
- Neither contains `"half"` (the contrast column name must not leak into the RE labels)
- When lme4 is available: lme4 uses `fa`/`fb` labels; mixeff uses `f: a`/`f: b` — the test only asserts set-equality of the lme4 names, not parity of SD values

### Tolerances
No numeric tolerances; purely structural label assertions.

### Skip conditions
- lme4 comparison block uses `requireNamespace("lme4", quietly = TRUE)` (soft; no skip — the comparison is conditional but the primary assertions always run)
- No `skip_on_cran()`, no slow gate

### What is NOT tested
- `fixef` parity vs. lme4 for this contrast design
- `VarCorr` SD parity vs. lme4 (the contrast affects the reparameterisation; only names are checked)
- `fitted` / `residuals` parity
- The symmetric case: `(0 + f | g)` with treatment coding (no custom contrasts) — ensures the default-contrast path is exercised
- Multiple factor levels (>2) with no-intercept slopes
- Interaction `(0 + f:x | g)` where f is a factor and x is numeric
- Explicit `(1 + f | g)` alongside `(0 + f | g)` to verify the intercept is not double-counted

---

## Cross-cutting gaps (present across multiple files)

1. **`ranef()` BLUP value parity** — no tutorial file asserts that BLUPs numerically match lme4. Only `test-pw2` checks group names structurally; all others skip ranef entirely.
2. **`deviance()` parity** — absent from all 8 files; present only in `helper-lme4-parity.R` via `mm_expect_core_lme4_parity` (not called from these tutorial tests).
3. **`model.matrix()` / design matrix parity** — absent from all 8 files.
4. **`theta` (Cholesky factor) parity** — absent from all 8 files.
5. **`predict()` on new data** — absent from all 8 files; only in-sample `fitted()` is compared.
6. **`compare()` table columns beyond LRT stat/p-value** — the AIC/BIC/logLik columns of the comparison table are never asserted.
7. **Singular/boundary fit diagnostics** — no test checks that a near-zero variance component triggers the expected structured diagnostic (rather than silently producing zero).
8. **`print()` / `summary()` smoke tests** — no tutorial file calls `print()` or `summary()` on a fit object.
9. **`nobs()` / `df.residual()` parity** — absent from all 8 files.
10. **REML vs. ML `compare()` guard** — no test verifies that comparing two REML-fitted models raises a diagnostic or auto-refits.
