# Test Coverage Survey: GLMM Test Files
**Family:** tests-1
**Files surveyed:** `tests/testthat/test-glmm.R`, `tests/testthat/test-glmm-summary-tests.R`, `tests/testthat/test-emmeans-glmm.R`
**Survey date:** 2026-05-31

---

## 1. `test-glmm.R`

### What is tested

This file provides the lme4 parity harness for `mm_glmm` fits. It defines helper
functions that are local to this file (not in the `helper-*.R` infrastructure) and
two `test_that` blocks:

1. **Manifest invariants** — verifies that exactly two parity cases are declared
   (`cbpp_binomial_logit_profiled_pirls`, `grouseticks_poisson_log_profiled_pirls`),
   that both carry `model = "glmm"`, `expected_status = "expected_mismatch"`, and
   that ledger entries exist for the four fields `fixef`, `theta`, `logLik`,
   `deviance`.

2. **lme4 parity within classified bounds** — for each case, fits both `glmm()`
   and `lme4::glmer()` and asserts:
   - `expect_s3_class(fit, "mm_glmm")`
   - `nobs` equality (exact integer match)
   - `fit$family$family` and `fit$family$link` match the case spec
   - `fit$method` matches the case spec
   - `fixef` within `tol$fixef` (1e-4)
   - `theta` within `tol$theta` (1e-3)
   - `logLik` within `tol$logLik` (1e-4)
   - `deviance` vs `−2 * logLik(lme4)` within `tol$deviance` (1e-4)

### Tolerances and assertions

| Quantity | Tolerance | Assertion style |
|---|---|---|
| `fixef` | 1e-4 | `mm_assert_parity` (scoreboard-aware) |
| `theta` | 1e-3 | `mm_assert_parity` |
| `logLik` | 1e-4 | `mm_assert_parity` |
| `deviance` | 1e-4 | `mm_assert_parity` vs `−2 * logLik(lme4)` |
| `nobs` | exact | `expect_identical` |
| class | exact | `expect_s3_class` |
| family/link/method | exact | `expect_identical` |

### Skip conditions

- Both test blocks call `mm_skip_if_no_lme4()` (`skip_if_not_installed("lme4")`).
- No `skip_on_cran`.

### What is NOT tested (gaps)

- **AIC / BIC parity**: not asserted against `lme4::AIC`/`BIC` for GLMM cases.
- **sigma/dispersion**: not asserted (binomial/Poisson sigma is 1.0 by contract but
  the assertion is absent; for Gamma it is non-trivial).
- **fitted values parity**: no `fitted(fit)` vs `fitted(lme4)` comparison.
- **residuals parity**: no `residuals()` comparison.
- **VarCorr parity**: no `VarCorr(fit)` vs `lme4::VarCorr(ref)` comparison.
- **ranef parity**: no `ranef(fit)` vs `lme4::ranef(ref)` comparison.
- **model.matrix parity**: not checked for GLMM cases.
- **Gamma family case**: only binomial (cbpp) and Poisson (grouseticks) cases exist
  in the fixture; no Gamma/log parity case.
- **Non-default links**: probit and cloglog (binomial) and sqrt (Poisson) have no
  parity fixture cases.
- **Prediction on GLMM**: `predict.mm_glmm` intentionally raises
  `mm_inference_unavailable` — this deliberate refusal is not tested here (covered
  partially in `test-phase4.R`).
- **Multi-random-effects theta vector ordering**: the grouseticks case has three
  RE terms; no test verifies that the theta vector ordering matches lme4's convention.

---

## 2. `test-glmm-summary-tests.R`

### What is tested

Exercises `summary.mm_glmm` with `tests = "coefficients"` and `tests = "none"` on
a single synthetic binomial GLMM (`y ~ x * g + (1 | subject)`, seed 41). Five
`test_that` blocks:

1. **No error on `tests='coefficients'`** — `expect_no_error` + `expect_s3_class`.

2. **`tests='none'` path** — verifies `summary.mm_glmm` class, `inference` is NULL,
   coefficient table has `"Estimate"` column and `df` is NA/numeric.

3. **Wald z arithmetic** — for each named fixed effect, calls `mm_lincomb(fit, …)` and
   asserts `Estimate` (tol 1e-12), `Std. Error` (tol 1e-12), `z value` (tol 1e-10),
   and `Pr(>|z|)` (tol 1e-10) match the `mm_lincomb` result.

4. **`vcov_status` structure** — checks `is.list`, required keys present, and that
   for the PIRLS-Laplace path `status = "available"` and
   `method = "pirls_laplace_working_hessian"`.

5. **Print reliability notice** — `capture.output(print(s))` checks for
   `"Fixed effects:"` and `"Wald-z reliability: moderate"`.

6. **Column labels** — verifies `"Estimate"`, `"Std. Error"`, `"z value"`,
   `"Pr(>|z|)"` in `colnames(s$coefficients)`.

### Tolerances and assertions

All numeric checks use `expect_equal` with either 1e-12 or 1e-10. Structural checks
use `expect_true`, `expect_null`, `expect_identical`, `expect_named`.

### Skip conditions

None — no `skip_*` calls. Runs on CRAN.

### What is NOT tested (gaps)

- **Poisson and Gamma families**: fixture uses only binomial/logit. The Wald-z path
  for non-binomial families (where `dispersion != 1`) is not exercised.
- **`tests = "coefficients"` after `revive()`**: the `vcov` matrix may be degraded
  after serialization/revival; no test checks that `summary(revived_fit, tests="coefficients")`
  produces NA statistics with an appropriate `vcov_status`.
- **`summary()` with zero fixed effects**: edge case not tested.
- **Interaction term labels**: the synthetic model has `x:gtreat` but the column label
  format for interaction terms is not snapshot-tested.
- **`print.summary.mm_glmm` full output structure**: only two strings are checked in
  the captured output; no snapshot test anchors the full print layout.
- **`dispersion != 1` in coefficient table header**: for Gamma family the dispersion
  is estimated; the summary print path for that is not tested.
- **`reliability = "unavailable"` branch**: the code has a branch where z/p are set
  to NA when `method_used = "not_computed"`; this branch has no test.
- **`summary()` called on a fit where `std_errors` length mismatches `beta`**: the
  fallback `sqrt(diag(V))` path has no test.

---

## 3. `test-emmeans-glmm.R`

### What is tested

Covers the `emmeans` integration layer (`recover_data.mm_glmm` +
`emm_basis.mm_glmm`) on a single synthetic binomial GLMM
(`y ~ x * g + (1 | subject)`, seed 71). File-level `skip_if_not_installed` guards
for both `emmeans` and `estimability`. Five `test_that` blocks:

1. **Reference grid for categorical factor** — `emmeans(fit, ~ g)` returns
   `emmGrid`, has `g`, `emmean`, `SE` columns, two rows, all finite link-scale
   estimates.

2. **`type='response'` applies linkinv** — checks that the response-scale column
   (`prob`/`response`/`rate`) is in `(0, 1)` for binomial logit.

3. **Response scale matches manual `plogis()`** — compares `em_resp[[prob_col]]`
   to `plogis(s_link$emmean)` at tolerance 1e-10.

4. **Pairwise contrast uses asymptotic z (df=Inf)** — verifies `z.ratio` column
   present and finite.

5. **`emm_basis.mm_glmm` structure** — directly calls `recover_data.mm_glmm` and
   `emm_basis.mm_glmm` and checks: correct named components (`X`, `bhat`, `nbasis`,
   `V`, `dffun`, `dfargs`, `misc`); `bhat` length and `V` dimensions match
   `fixef(fit)` length; `dffun` returns `Inf`.

6. **`misc` carries link info** — checks that `b$misc$tran` or `b$misc$inv.lbl` is
   non-NULL.

### Tolerances and assertions

| Assertion | Tolerance / type |
|---|---|
| Response-scale vs manual | `expect_equal`, 1e-10 |
| `z.ratio` finite | `expect_true(all(is.finite(...)))` |
| Dimensions/structure | `expect_named`, `expect_equal(length(...))` |
| `dffun(…)` | `expect_identical(..., Inf)` |
| Probability bounds | `expect_true(all(prob > 0 & prob < 1))` |

### Skip conditions

- File-level `skip_if_not_installed("emmeans")` and `skip_if_not_installed("estimability")`.
- No `skip_on_cran`.

### What is NOT tested (gaps)

- **Poisson and Gamma families**: all tests use binomial/logit. `emmeans` with
  Poisson/log (rate-scale backtransform) and Gamma/log are not exercised.
- **Non-default links (probit, cloglog)**: `type='response'` backtransform for
  probit/cloglog is not tested.
- **Continuous predictor marginal means**: tests only use `~ g` (categorical). `emmeans(fit, ~ x)` at a grid of x-values is not tested.
- **`emmeans()` via the top-level dispatch (not direct `emm_basis` call)**: tests 5
  and 6 call `recover_data.mm_glmm` and `emm_basis.mm_glmm` directly; the
  standard `emmeans(fit, ~ ...)` dispatch path is tested only in tests 1–4 but those
  don't verify `bhat`/`V` dimensions.
- **`contrast()` estimate value parity**: the pairwise contrast test checks structure
  (columns, finiteness) but does not compare the z-ratio or p-value to an
  independent calculation.
- **Multiplicity adjustment**: no test covers `emmeans(..., adjust = "bonferroni")`
  or similar.
- **`emmeans()` on a revived fit**: no test checks that serialized/revived `mm_glmm`
  objects work correctly through the emmeans dispatch.
- **`emmeans(fit, ~ g * x)` interaction grid**: interaction grids (the common
  `~ g | x` or `~ g * x` form) are not tested.
- **`recover_data` attribute contract**: the test accesses `attr(rd, "terms")` and
  `attr(rd, "xlev")` but does not verify the full `recover_data` contract (e.g.,
  that `attr(rd, "call")` and `attr(rd, "predictors")` are correct).
- **`nbasis` correctness**: `b$nbasis` is checked only for presence in
  `expect_named`; its numeric content is not validated.

---

## 4. Cross-cutting gaps (all three files)

These gaps apply across the three surveyed files and are not covered elsewhere in
the GLMM-specific tests (though some are partially covered in `test-phase4.R`):

- **`predict.mm_glmm` refusal message/class**: `predict.mm_glmm` raises
  `mm_inference_unavailable`; none of the three surveyed files test this. Covered
  in `test-phase4.R` only for the revived case.
- **`ranef.mm_glmm(condVar=TRUE)` refusal**: the `condVar=TRUE` path raises with
  reason `"random_effect_conditional_variance_unavailable_for_glmm"`; not tested in
  any of the three files.
- **`vcov.mm_glmm` matrix properties**: the vcov matrix is used implicitly in
  summary tests, but its symmetry, positive-semi-definiteness, and dimension are
  not directly asserted.
- **`coef.mm_glmm` output shape**: not tested in any of the three files.
- **`model.matrix.mm_glmm`**: not tested.
- **`model.frame.mm_glmm`**: not tested.
- **`formula.mm_glmm`**: not tested.
- **`df.residual.mm_glmm`**: not tested.
- **`ngrps`-equivalent**: no test checks the number of grouping levels via any
  extractor.
- **`glmm()` input validation errors** (`mm_arg_error`): `nAGQ < 1`, `nAGQ` not
  integer-like, `weights`/`subset`/`contrasts` non-NULL — none tested in the three
  surveyed files (partially covered in `test-phase4.R`).
- **Gamma family parity case**: no fixture case for Gamma/log GLMM parity.

---

## 5. Summary table

| Gap | Severity | Classification |
|---|---|---|
| AIC/BIC parity not asserted for GLMM parity cases | minor | test-gap |
| fitted/residuals parity missing from GLMM parity cases | minor | test-gap |
| VarCorr parity missing from GLMM parity cases | minor | test-gap |
| ranef parity missing from GLMM parity cases | minor | test-gap |
| sigma/dispersion parity missing | minor | test-gap |
| No Gamma family parity fixture case | major | test-gap |
| No non-default-link parity cases (probit, cloglog, sqrt) | minor | test-gap |
| `summary(revived_fit, tests='coefficients')` NA-stats path not tested | major | test-gap |
| `vcov_status = "unavailable"` branch in summary not tested | minor | test-gap |
| Poisson/Gamma families in summary Wald-z tests | minor | test-gap |
| `ranef.mm_glmm(condVar=TRUE)` refusal not tested (three files) | minor | test-gap |
| `predict.mm_glmm` refusal not tested in three files | minor | test-gap |
| `emmeans` Poisson/Gamma `type='response'` backtransform not tested | minor | test-gap |
| `emmeans` on revived fit not tested | minor | test-gap |
| `contrast()` estimate/p-value parity not verified | minor | test-gap |
| `emm_basis` `nbasis` content not validated | cosmetic | test-gap |
| `recover_data` full attribute contract not validated | cosmetic | test-gap |
| `coef.mm_glmm` / `model.matrix.mm_glmm` / `formula.mm_glmm` not tested | minor | test-gap |
| `predict.mm_glmm` GLMM prediction is `out-of-scope-by-design` per current Rust contract | — | out-of-scope-by-design |
| Profile-LL CIs for GLMM (`confint(method="profile")`) not tested | — | out-of-scope-by-design |
| `joint_laplace` method fit | — | upstream-blocked |
