# Test Survey: tests-4 — test-marginal.R, test-emmeans.R, test-lincomb.R

Survey date: 2026-05-31  
Files examined:
- `tests/testthat/test-marginal.R`
- `tests/testthat/test-emmeans.R`
- `tests/testthat/test-lincomb.R`

Related source:
- `R/marginal.R` — mm_grid / mm_predictions / mm_means / mm_comparisons
- `R/emmeans.R` — recover_data.mm_lmm/glmm, emm_basis.mm_lmm/glmm
- `R/inference.R` (lines 1874–2075) — mm_lincomb.mm_lmm/mm_glmm/default

---

## test-marginal.R

### What is tested

| Test | Summary |
|---|---|
| `mm_grid() builds a fixed-effect reference grid` | Grid shape, column names, X ncol vs fixef length, trt/block cells crossed with `at=list(x=0)` |
| `mm_predictions() returns contract-shaped population prediction rows` | Full column-contract check, quantity/target/scale/comparison/method fields, uniqueness of labels, finiteness of estimates and CIs |
| `mm_means() averages nuisance fixed-factor cells through contrast()` | Row count equals nlevels(trt), L row equals colMeans of grid$X filtered to that level, estimate matches L %*% fixef exactly (tol 1e-8) |
| `mm_comparisons() returns pairwise differences within by groups` | 6 rows for 3-level trt within 2-level block, quantity/comparison/status fields, label contains " - ", by/specs list columns |
| `weighted LMM covariance feeds marginal standard errors` | Weighted fit produces different vcov than unweighted; SE for predictions and means matches sqrt(diag(L V L')) (tol 1e-8); mm_status/mm_schema_name attributes checked |
| `unsupported marginal quantity requests are typed` | `mm_grid(fit, ~ missing)` → `mm_arg_error`; `mm_comparisons(fit, ~ trt, comparison = "ratio")` → `mm_inference_unavailable` |

### Tolerances used
- Estimate vs L%*%fixef: 1e-8
- SE vs manual sqrt(LVL'): 1e-8

### Skips
None (`skip_on_cran` / `skip_if_not_installed` absent).

### What is NOT tested

1. **`mm_means()` and `mm_comparisons()` with `weights = "proportional"`** — The `mm_cell_weights()` branch for proportional weighting is entirely uncovered. The `equal` path is the only one exercised.

2. **`mm_predictions()` / `mm_means()` / `mm_comparisons()` with inference method `"auto"` or `"satterthwaite"`** — Every call in the suite uses `method = "asymptotic"`. The auto-dispatch and Satterthwaite branches of the marginal functions are not exercised; the df-bearing CI path in `mm_marginal_intervals()` is not reached.

3. **`mm_grid()` with the `by` argument separate from the formula pipe** — The `by` argument is only exercised via the formula `~ trt | block` inside `mm_comparisons()`. Direct use of `mm_grid(fit, ~trt, by = "block")` is not tested.

4. **`mm_predictions()` called without a pre-built grid (via `specs` argument)** — All `mm_predictions()` calls receive an explicit `grid`; the `mm_resolve_grid` path that builds a grid on the fly from `specs` is not exercised for predictions.

5. **`mm_grid()` with a numeric `specs` variable** — All tests use a factor spec. The grid with a numeric displayed variable (≤10 unique values returned as-is) and the `cov.reduce` fallback for wide numeric ranges are not tested.

6. **`mm_grid()` with a custom `cov.reduce` function** — Only the default (`mean`) is used; the error branch for a `cov.reduce` that returns a wrong-length or NA value is not tested.

7. **`mm_grid()` with `at` supplying an out-of-range factor level** — The `mm_grid_values_like()` error path for unknown levels is not triggered.

8. **`mm_comparisons()` on a specs variable with only one level (triggering the "needs at least two means" error)** — The `mm_pairwise_rows()` empty-pairs branch is not tested.

9. **`mm_means()` / `mm_comparisons()` with an explicit pre-built `grid` passed in** — Every call rebuilds the grid internally; the `grid =` argument path is not separately tested.

10. **`print.mm_grid` and `print.mm_marginal_quantity`** — No test verifies the print methods run without error or produce expected output.

11. **Interaction models (e.g., `y ~ trt * x + (1|subject)`)** — All tests use additive fixed effects; interaction term expansion into the reference grid is not exercised.

---

## test-emmeans.R

### What is tested

| Test | Summary |
|---|---|
| `optional emmeans methods reproduce native marginal mean estimates` | `emmeans::emmeans(mm_lmm)` estimates match `mm_means()` (tol 1e-8); SE finite; df infinite (asymptotic); `basis$V` attributes `mm_method`/`mm_status`/`mm_schema_name` checked; init message includes covariance-schema string |
| `optional emmeans pairwise estimates agree with fixed-effect differences` | `contrast(..., method="pairwise")` a-b and a-c estimates match `-fixef["trtb"]` and `-fixef["trtc"]` (tol 1e-8); SE finite |
| `emmeans support methods are exported for conditional registration` | `recover_data.mm_lmm` and `emm_basis.mm_lmm` are in `getNamespaceExports("mixeff")`; both are registered inside the emmeans namespace |

### Tolerances used
- Estimate agreement with `mm_means()`: 1e-8
- Pairwise difference vs fixef: 1e-8

### Skips
`skip_if_not_installed("emmeans")` and `skip_if_not_installed("estimability")` on the first two tests; `skip_if_not_installed("emmeans")` only on the export test. No `skip_on_cran`.

### What is NOT tested

1. **`recover_data.mm_lmm` and `emm_basis.mm_glmm` / `recover_data.mm_glmm` exports** — The export-check test only confirms `mm_lmm` S3 registration. `recover_data.mm_glmm` and `emm_basis.mm_glmm` export membership is not asserted in this file (covered in `test-emmeans-glmm.R` but not here in the LMM-focused file).

2. **`emm_basis.mm_lmm` with `method = "satterthwaite"` or `method = "auto"` producing finite df** — The only call uses `method = "asymptotic"`, so the `dffun` branch that routes through `df_for_contrast()` and returns a finite value is never hit.

3. **`emm_basis.mm_lmm` init message when vcov is unavailable** — The `mm_emmeans_init_messages()` unavailable branch (status != "available") is not exercised; no test constructs a rank-deficient or otherwise degraded fit.

4. **`recover_data.mm_lmm` with an explicit `data =` override** — The `data %||% object$model_frame` branch for the non-NULL data path is not tested.

5. **`emmeans::emmeans(mm_lmm)` pairwise SE matches `mm_comparisons()` SE** — Only the point estimates are checked against fixef differences; no test cross-validates the SE produced by emmeans against the native `mm_comparisons()` standard errors.

6. **`emm_basis.mm_lmm` when emmeans/estimability are absent** — The error path (`mm_inference_unavailable`) for missing optional packages is not triggered.

7. **`emmeans` grid with a numeric covariate at a non-zero value** — All calls use `at = list(x = 0)`; no test varies the numeric covariate to check that `emm_basis` X-matrix varies correctly.

---

## test-lincomb.R

### What is tested

| Test | Summary |
|---|---|
| `mm_lincomb() reproduces hand-rolled Wald z on mm_glmm` | Two-coefficient interaction lincomb against hand_wald; estimate/SE/statistic/p/CI (tol 1e-10–1e-12); statistic_name="z"; df=NA |
| `mm_lincomb() exposes the underlying vcov status as an attribute` | `attr(out, "mm_status")` is a list with status/method/reliability/reason; status = "available" |
| `mm_lincomb() with method='asymptotic' on mm_lmm matches hand-rolled Wald z` | Named numeric weights; estimate/SE/p (tol 1e-10–1e-12); statistic_name="z"; df=NA |
| `mm_lincomb() default on mm_lmm uses Satterthwaite df via df_for_contrast()` | Default method="auto" → statistic_name="t", finite df>0; cross-validated against `contrast()` output (tol 1e-6–1e-10) |
| `mm_lincomb() accepts named list and 1-row data.frame` | Three input forms (vector, list, 1-row df) all produce identical estimate and SE (tol 1e-12) |
| `mm_lincomb() rejects malformed weights` | NULL, unnamed vector, NA value, duplicate names, unknown coefficient name, multi-row df → all `mm_arg_error` |
| `mm_lincomb() rejects unsupported method on mm_glmm` | `method="satterthwaite"` and `method="kenward_roger"` on mm_glmm → `mm_arg_error` |
| `mm_lincomb() default method errors on non-fit input` | `mm_lincomb(list(), ...)` → `mm_arg_error` |
| `mm_lincomb() level argument moves the CI as expected` | 0.99 interval wider than 0.95; estimate and SE identical |

### Tolerances used
- Estimate / SE: 1e-12
- Statistic / p_value: 1e-10
- df (vs contrast()): 1e-6

### Skips
None.

### What is NOT tested

1. **`mm_lincomb(mm_lmm, ..., method = "kenward_roger")` happy path** — `"kenward_roger"` is in `match.arg()` for `mm_lincomb.mm_lmm` and routed through `df_for_contrast(fit, L, method="kenward_roger")`. The only KR appearance in tests is the error test for `mm_glmm`; no test verifies KR produces a finite df and t-statistic on an `mm_lmm`.

2. **`mm_lincomb(mm_lmm, ..., method = "satterthwaite")` with infinite or non-positive df** — The code path where `df` comes back non-finite sets `p_value = NA_real_` and `lower/upper = NA_real_`. No test exercises this fallback.

3. **`mm_lincomb()` on a single-coefficient intercept-only lincomb** — No test sets `weights = c("(Intercept)" = 1)`, though this is the simplest valid case and exercises a different region of the V matrix.

4. **`mm_lincomb()` with `method = "auto"` on mm_glmm** — `method = NULL` is accepted as equivalent to "asymptotic" (the guard allows NULL and "auto"), but no test explicitly passes `method = "auto"` to `mm_lincomb.mm_glmm`.

5. **`mm_lincomb()` output column names are verified** — No test checks `names(out)` against the expected schema columns (`estimate`, `std_error`, `statistic`, `statistic_name`, `df`, `p_value`, `lower`, `upper`, `method`).

6. **`mm_lincomb()` on a revived (de-serialised) fit** — No round-trip test verifies that a fit saved to JSON and revived still produces identical lincomb results.

7. **`mm_lincomb()` with weights that span all fixed-effect terms** — Current tests use 1–2 non-zero coefficients; no test exercises a fully dense weight vector to stress the full V %*% w computation.

8. **`mm_lincomb()` mm_status attribute when vcov is unavailable** — The `mm_lincomb_status_from_vcov()` path for a degraded/unavailable covariance matrix is not covered; no test constructs a singular or rank-deficient fit.
