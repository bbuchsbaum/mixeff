# Test Specification — emmeans / Marginal Means Family

**Family:** emmeans / marginal means  
**Source gap report:** `assessment/gap/emmeans.md`  
**Parity probe:** `assessment/parity/inf-emmeans.md` (cake dataset)  
**In-scope gaps addressed:** 4 (1 partial×3 sub-gaps + 1 in-scope-missing + 1 test-gap promoting 2 latent bugs to verified defects)  
**Total test specs:** 11  
**Written:** 2026-05-31

---

## Classification summary

| Gap from gap report | Classification | Severity | Test specs |
|---|---|---|---|
| `mode=` / `lmer.df=` silently ignored by `emm_basis.mm_lmm` | partial | major | TS-01, TS-02, TS-03 |
| Rank-deficient design: estimable cells lose SE/df; non-estimable cell shows fabricated mean | partial | major | TS-04, TS-05 |
| Test gap: existing tests never exercise `mode=` / rank-deficient paths | test-gap | major | (root cause for TS-01–05) |
| `emm_options(lmer.df=)` global default has no effect (same root as `mode=`) | partial | minor | TS-06 |
| `lmerTest::show_tests()` — no equivalent | in-scope-missing | minor | TS-07 |
| Reference-grid averaging broken for interaction models (cake parity probe) | partial | blocker | TS-08, TS-09, TS-10 |
| `mm_means` / `mm_comparisons` wrong on interaction models (same upstream bug) | partial | blocker | TS-11 |

---

## TS-01 — `mode="asymptotic"` yields `df = Inf` via emmeans

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** `mode=` being silently swallowed by `...`, causing an intended asymptotic z-test to silently return a finite Satterthwaite t-test — a no-silent-surgery violation (PRD/CLAUDE.md).

### Dataset / formula
`sleepstudy` (lme4 built-in, n=180).  
`lmm(Reaction ~ Days + (1 | Subject), sleepstudy)`

### Assertion
```r
test_that("emmeans mode='asymptotic' yields df = Inf for mm_lmm", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  em <- emmeans::emmeans(fit, ~ Days, at = list(Days = c(0, 5, 9)),
                         mode = "asymptotic")
  s  <- as.data.frame(summary(em))

  # All df must be Inf when mode = "asymptotic"
  expect_true(all(is.infinite(s$df)),
              info = paste("Expected Inf df with mode='asymptotic'; got",
                           paste(round(s$df, 2), collapse = ", ")))
})
```

### Tolerance
`df == Inf` (exact, not numeric).

### Notes
- Fix requires `emm_basis.mm_lmm` to accept `mode` and `lmer.df` in its signature (not just `...`), normalize hyphen→underscore, and map to the internal `method=` token.
- Upstream fixture: none needed; pure R bridge fix.

---

## TS-02 — `mode="satterthwaite"` / `mode="kenward-roger"` yield finite df matching `df_for_contrast`

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** `mode=` tokens for finite-df methods being silently ignored (same root as TS-01); also ensures hyphen→underscore normalisation (`"kenward-roger"` → `"kenward_roger"`).

### Dataset / formula
Same as TS-01 (`sleepstudy`, random-intercept model).

### Assertion
```r
test_that("emmeans mode='satterthwaite' and mode='kenward-roger' yield finite df", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  em_satt <- emmeans::emmeans(fit, ~ Days, at = list(Days = 0),
                              mode = "satterthwaite")
  em_kr   <- emmeans::emmeans(fit, ~ Days, at = list(Days = 0),
                              mode = "kenward-roger")

  df_satt <- as.data.frame(summary(em_satt))$df
  df_kr   <- as.data.frame(summary(em_kr))$df

  # Must be finite
  expect_true(all(is.finite(df_satt)),
              info = "mode='satterthwaite' should return finite df")
  expect_true(all(is.finite(df_kr)),
              info = "mode='kenward-roger' (hyphen) should return finite df")

  # Must agree with df_for_contrast to within 0.5 df units
  L <- matrix(c(0, 1), nrow = 1)          # contrast on Days coefficient
  native_satt <- as.numeric(df_for_contrast(fit, L, method = "satterthwaite")[[1]])
  native_kr   <- as.numeric(df_for_contrast(fit, L, method = "kenward_roger")[[1]])

  expect_equal(df_satt, native_satt, tolerance = 0.5)
  expect_equal(df_kr,   native_kr,   tolerance = 0.5)
})
```

### Tolerance
df numeric comparison: absolute tolerance 0.5 (df units).

### Notes
- Verifies hyphen variant `"kenward-roger"` is normalised to underscore `"kenward_roger"` inside `emm_basis.mm_lmm`.
- Also verifies the deprecated alias `lmer.df=` maps correctly (see TS-03).

---

## TS-03 — `lmer.df=` deprecated alias is honoured

**Priority:** P1  
**Kind:** parity-vs-lme4 / error-message (no-silent-surgery)  
**Guards against:** `lmer.df=` (the legacy emmeans global-option style argument) being silently dropped.

### Dataset / formula
Same as TS-01 (`sleepstudy`).

### Assertion
```r
test_that("emmeans lmer.df= alias is honoured by emm_basis.mm_lmm", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  em_inf  <- emmeans::emmeans(fit, ~ Days, at = list(Days = 0),
                              lmer.df = "asymptotic")
  em_fin  <- emmeans::emmeans(fit, ~ Days, at = list(Days = 0),
                              lmer.df = "satterthwaite")

  df_inf <- as.data.frame(summary(em_inf))$df
  df_fin <- as.data.frame(summary(em_fin))$df

  expect_true(all(is.infinite(df_inf)),
              info = "lmer.df='asymptotic' should yield Inf df")
  expect_true(all(is.finite(df_fin)),
              info = "lmer.df='satterthwaite' should yield finite df")
})
```

### Tolerance
Exact `Inf` / finite check.

---

## TS-04 — Rank-deficient design: non-estimable cell is flagged `nonEst`

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / error-message  
**Guards against:** `nbasis = estimability::all.estble` hardcode causing `emm_basis.mm_lmm` to silently return a fabricated numeric mean for a cell that is not estimable (PRD "no fabrication" contract violation).

### Dataset / formula
Synthetic missing-cell `A*B` factorial: factor `A` with levels `{a1, a2, a3}`, factor `B` with levels `{b1, b2}`; cell `(a3, b2)` is empty (zero observations). Random intercept `(1 | subject)`.

```r
mk_missing_cell_fit <- function() {
  set.seed(99L)
  A <- factor(c(rep("a1", 20), rep("a2", 20), rep("a3", 10)))
  B <- factor(c(rep("b1", 10), rep("b2", 10),   # a1
                rep("b1", 10), rep("b2", 10),   # a2
                rep("b1", 10)))                  # a3 b1 only; b2 is missing
  subject <- factor(rep(seq_len(12L), length.out = length(A)))
  y <- rnorm(length(A)) + as.numeric(A) * 0.5 + as.numeric(B) * 0.3
  lmm(y ~ A * B + (1 | subject),
      data.frame(y=y, A=A, B=B, subject=subject),
      control = mm_control(verbose = -1))
}
```

### Assertion
```r
test_that("emmeans marks non-estimable cell nonEst in rank-deficient design", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("estimability")
  fit <- mk_missing_cell_fit()

  em  <- emmeans::emmeans(fit, ~ A * B)
  s   <- as.data.frame(summary(em))

  # The (a3, b2) cell must be flagged nonEst, not a numeric value
  row <- s[s$A == "a3" & s$B == "b2", ]
  expect_equal(nrow(row), 1L)
  expect_true(is.na(row$emmean) || identical(as.character(row$emmean), "NonEst"),
              info = "Non-estimable cell (a3, b2) must be NA / NonEst, not a fabricated numeric")
})
```

### Tolerance
Categorical: cell value must be `NA` or `"NonEst"` (emmeans `nonEst` marker), not a finite number.

### Notes
- Fix requires deriving `nbasis` from the engine's reported rank/null space rather than hardcoding `estimability::all.estble`.
- May require upstream Rust engine to expose the null-space basis, or can be computed in R from the `X_train` column rank.
- **Upstream fixture needed if** the null-space basis must come from the Rust engine's fixed-effect design info.

---

## TS-05 — Rank-deficient design: estimable cells retain valid SE and df

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** V being entirely unavailable for a rank-deficient fit, causing *all* cells — including estimable ones — to lose SE and df (a regression against lme4 behaviour where only the non-estimable cell is `nonEst`).

### Dataset / formula
Same missing-cell fit as TS-04.

### Assertion
```r
test_that("emmeans estimable cells have finite SE/df in rank-deficient design", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("estimability")
  fit <- mk_missing_cell_fit()

  em <- emmeans::emmeans(fit, ~ A * B)
  s  <- as.data.frame(summary(em))

  # All cells EXCEPT (a3, b2) must have finite SE and df
  estimable_rows <- !(s$A == "a3" & s$B == "b2")
  expect_true(all(is.finite(s$SE[estimable_rows])),
              info = "Estimable cells must have finite SE")
  expect_true(all(is.finite(s$df[estimable_rows])),
              info = "Estimable cells must have finite df")
  # Sanity: at least 5 estimable cells (3x2 minus 1)
  expect_gte(sum(estimable_rows), 5L)
})
```

### Tolerance
Finite / non-NA check. No specific numeric tolerance; the SE values themselves are not compared to lme4 numerically here (guarded separately by parity probe if estimates are fixed).

### Notes
Depends on the same fix as TS-04 (proper `nbasis` derivation). Both TS-04 and TS-05 must pass together to certify rank-deficient parity.

---

## TS-06 — `emm_options(lmer.df=)` global default is respected

**Priority:** P2  
**Kind:** integration  
**Guards against:** The `emm_options(lmer.df = "asymptotic")` global setting having no effect on mixeff fits (same root as `mode=` / `lmer.df=` being dropped). This is minor because it is a convenience alias for the session-level default.

### Dataset / formula
`sleepstudy`, same random-intercept model as TS-01.

### Assertion
```r
test_that("emm_options(lmer.df='asymptotic') changes mixeff df to Inf", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  old <- emmeans::emm_options(lmer.df = "asymptotic")
  on.exit(emmeans::emm_options(old), add = TRUE)

  em <- emmeans::emmeans(fit, ~ Days, at = list(Days = 0))
  df <- as.data.frame(summary(em))$df

  expect_true(all(is.infinite(df)),
              info = "Global emm_options(lmer.df='asymptotic') should drive mixeff df to Inf")
})
```

### Tolerance
Exact `Inf` check.

### Notes
This test will pass automatically once `mode=` / `lmer.df=` is read in `emm_basis.mm_lmm` (emmeans passes the resolved option value as `lmer.df=` into `emm_basis`). No additional code change needed beyond TS-01 fix.

---

## TS-07 — `show_tests()` equivalent or documented absence

**Priority:** P2  
**Kind:** error-message (documented divergence)  
**Guards against:** A user calling the lmerTest function `lmerTest::show_tests()` and getting an uninformative error rather than a clear "not implemented" message with a pointer to the mixeff equivalent.

### Assertion
This is a **documentation assurance**, not a purely numeric test. The spec requires either:

(a) A `show_tests.mm_lmm` method that produces output equivalent to `lmerTest::show_tests()` (the design matrix of contrasts underlying an F-test), or  
(b) A clear, structured condition with class `mm_inference_unavailable` and a `message` that explicitly names the mixeff equivalent (e.g., `contrast(fit, L)` or `test_effect(fit, term)`).

```r
test_that("show_tests or its absence produces an informative message", {
  skip_if_not_installed("lmerTest")
  fit <- lmm(Reaction ~ Days + (1 | Subject),
             data = lme4::sleepstudy,
             control = mm_control(verbose = -1))

  # Option A: method exists and returns something printable
  if (existsMethod("show_tests", "mm_lmm") ||
      !is.null(getS3method("show_tests", "mm_lmm", optional = TRUE))) {
    result <- show_tests(fit)
    expect_true(is.data.frame(result) || is.list(result),
                info = "show_tests.mm_lmm should return a data frame or list")
  } else {
    # Option B: calling lmerTest::show_tests on a mixeff fit raises an
    # informative condition, not an opaque R dispatch error
    err <- tryCatch(
      lmerTest::show_tests(fit),
      error = function(e) e,
      warning = function(w) w
    )
    msg <- if (inherits(err, "condition")) conditionMessage(err) else ""
    expect_true(
      nchar(msg) > 0L,
      info = "Dispatching lmerTest::show_tests on mm_lmm should produce a message"
    )
  }
})
```

### Notes
- The gap report classifies `show_tests` as `in-scope-missing` / minor. The minimum bar is Option B (informative message). Implementing Option A is a bonus.
- No upstream fixture needed.

---

## TS-08 — `emmeans(fit, ~ recipe)` marginal means match lme4 on interaction model (cake)

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** The reference-grid averaging failure identified in `assessment/parity/inf-emmeans.md`: marginal means over a nuisance factor being evaluated at only the reference level instead of being averaged over all levels. This is the primary blocker for emmeans parity.

### Dataset / formula
`cake` (lme4 built-in, n=270).  
`lmm(angle ~ recipe * temperature + (1 | recipe:replicate), cake, REML = TRUE)`

### lme4 reference values (from parity probe, 2026-05-31)
```
recipe A: emmean=33.1, SE=1.74, df=42
recipe B: emmean=31.6, SE=1.74, df=42
recipe C: emmean=31.6, SE=1.74, df=42
```
(lme4 uses Kenward-Roger df by default; mixeff is expected to use Satterthwaite or KR — df tolerance is looser.)

### Assertion
```r
test_that("emmeans(mm_lmm) marginal means over interaction match lme4 on cake", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("lme4")
  data("cake", package = "lme4")

  fit_mm <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
                data = cake, REML = TRUE,
                control = mm_control(verbose = -1))

  fit_lme4 <- lme4::lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
                         data = cake, REML = TRUE)

  em_mm   <- as.data.frame(summary(emmeans::emmeans(fit_mm,   ~ recipe)))
  em_lme4 <- as.data.frame(summary(emmeans::emmeans(fit_lme4, ~ recipe)))

  # Estimates must agree to 1e-2 (generous, given df-method differences)
  expect_equal(sort(em_mm$emmean), sort(em_lme4$emmean),
               tolerance = 1e-2,
               info = "Marginal recipe means must average over all temperature levels")

  # SE must agree to 5% relative tolerance
  expect_equal(sort(em_mm$SE), sort(em_lme4$SE),
               tolerance = 0.05,
               info = "Marginal recipe SE must be consistent with full averaging")
})
```

### Tolerances
- estimate: absolute 1e-2 (generous because lme4 uses KR averaging weights; mixeff may use equal weights)
- SE: relative 5%

### Notes
- Root cause (from parity probe): `recover_data.mm_lmm` or the marginal averaging logic in `mm_group_basis` does not expand all levels of nuisance factors when building the reference grid for `~ recipe`. The averaging logic in `mm_group_basis` must include all temperature levels and average over them with equal or proportional weights.
- Both the emmeans bridge *and* `mm_means` share the same bug (confirmed in probe). Fixing `mm_group_basis` / `mm_grid.mm_lmm` will fix both paths.
- This test must pass before TS-09 and TS-10 are meaningful.

---

## TS-09 — `emmeans(fit, ~ temperature)` marginal means match lme4 on cake

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** Same reference-grid averaging failure as TS-08, applied to the numeric (ordered factor) temperature axis.

### Dataset / formula
Same as TS-08.

### lme4 reference values
```
temperature 175: emmean=28.0, SE=1.18, df=77.4
temperature 185: emmean=30.0, SE=1.18, df=77.4
temperature 195: emmean=31.4, SE=1.18, df=77.4
temperature 205: emmean=32.2, SE=1.18, df=77.4
temperature 215: emmean=35.8, SE=1.18, df=77.4
temperature 225: emmean=35.4, SE=1.18, df=77.4
```

### Assertion
```r
test_that("emmeans(mm_lmm) marginal temperature means match lme4 on cake", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("lme4")
  data("cake", package = "lme4")

  fit_mm   <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
                  data = cake, REML = TRUE,
                  control = mm_control(verbose = -1))
  fit_lme4 <- lme4::lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
                         data = cake, REML = TRUE)

  em_mm   <- as.data.frame(summary(emmeans::emmeans(fit_mm,   ~ temperature)))
  em_lme4 <- as.data.frame(summary(emmeans::emmeans(fit_lme4, ~ temperature)))

  expect_equal(em_mm$emmean, em_lme4$emmean, tolerance = 1e-2)
  expect_equal(em_mm$SE,     em_lme4$SE,     tolerance = 0.05)
})
```

### Tolerances
- estimate: 1e-2 absolute
- SE: relative 5%

### Notes
The most egregious failure from the probe is temperature=215: lme4=35.8, mixeff=19.3 (diff=16.5). Any fix that brings all 6 temperature estimates within 0.01 of lme4 certifies this row.

---

## TS-10 — `pairs(emmeans(fit, ~ recipe))` pairwise contrasts match lme4 on cake

**Priority:** P1  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** Pairwise contrasts being derived from wrong marginal means (downstream of TS-08 bug).

### Dataset / formula
Same as TS-08.

### lme4 reference values
```
A - B: estimate=1.478, SE=2.46, df=42, p=0.820
A - C: estimate=1.522, SE=2.46, df=42, p=0.810
B - C: estimate=0.044, SE=2.46, df=42, p=1.000
```

### Assertion
```r
test_that("pairs(emmeans(mm_lmm, ~ recipe)) matches lme4 on cake", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("lme4")
  data("cake", package = "lme4")

  fit_mm   <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
                  data = cake, REML = TRUE,
                  control = mm_control(verbose = -1))
  fit_lme4 <- lme4::lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
                         data = cake, REML = TRUE)

  pw_mm   <- as.data.frame(summary(
    emmeans::pairs(emmeans::emmeans(fit_mm,   ~ recipe))))
  pw_lme4 <- as.data.frame(summary(
    emmeans::pairs(emmeans::emmeans(fit_lme4, ~ recipe))))

  # Align by contrast label
  pw_mm   <- pw_mm[order(pw_mm$contrast), ]
  pw_lme4 <- pw_lme4[order(pw_lme4$contrast), ]

  expect_equal(pw_mm$estimate, pw_lme4$estimate, tolerance = 1e-2)
  expect_equal(pw_mm$SE,       pw_lme4$SE,       tolerance = 0.05)
  # p-values: qualitative agreement (same significance direction)
  expect_equal(sign(pw_mm$estimate), sign(pw_lme4$estimate))
})
```

### Tolerances
- estimate: 1e-2 absolute
- SE: relative 5%
- p-values: sign / direction only (df method differences can shift p-values)

### Notes
Depends on TS-08 (correct marginal means) being fixed first.

---

## TS-11 — `mm_means` / `mm_comparisons` match lme4 emmeans on interaction model (cake)

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** The native `mm_means` / `mm_comparisons` API returning wrong marginal means for interaction models. The parity probe confirmed `mm_means(fit_mm, ~ recipe)` returns the same wrong values as the emmeans bridge (29.1/26.9/27.9 vs lme4's 33.1/31.6/31.6), confirming the bug is in `mm_group_basis` upstream of both surfaces.

### Dataset / formula
Same as TS-08.

### Assertion
```r
test_that("mm_means matches lme4 emmeans marginal means on interaction model", {
  skip_if_not_installed("lme4")
  data("cake", package = "lme4")

  fit_mm   <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
                  data = cake, REML = TRUE,
                  control = mm_control(verbose = -1))
  fit_lme4 <- lme4::lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
                         data = cake, REML = TRUE)

  # Native API
  mm_r <- mm_means(fit_mm, ~ recipe, method = "asymptotic")

  # Reference from lme4+emmeans
  em_lme4 <- as.data.frame(
    summary(emmeans::emmeans(fit_lme4, ~ recipe, mode = "asymptotic")))

  # Align by recipe label
  native_est <- mm_r$table$estimate[order(mm_r$table$label)]
  lme4_est   <- em_lme4$emmean[order(as.character(em_lme4$recipe))]

  expect_equal(native_est, lme4_est, tolerance = 1e-2,
               info = "mm_means must average over all temperature levels, not just reference")

  # Comparisons (pairwise differences)
  comp <- mm_comparisons(fit_mm, ~ recipe, method = "asymptotic")
  expect_true(all(is.finite(comp$table$estimate)),
              info = "mm_comparisons estimates must be finite")
  expect_true(all(is.finite(comp$table$std_error)),
              info = "mm_comparisons SEs must be finite")
})
```

### Tolerances
- estimate: 1e-2 absolute (generous; KR vs Satterthwaite do not affect point estimates, only df)

### Notes
- Root fix is in `mm_group_basis` or the reference-grid construction in `mm_grid.mm_lmm`: nuisance factor columns must be expanded to all levels, not collapsed to the reference level.
- Once this test passes, TS-08 and TS-09 should also pass (same code path).
- `reliability = "low"` on `mm_means` output is acceptable so long as estimates are correct; the test does not assert `reliability = "high"`.

---

## Cross-cutting notes

### Fix dependency graph

```
TS-08, TS-09, TS-11  <-- fix mm_group_basis / reference-grid averaging in mm_grid.mm_lmm
TS-10                <-- depends on TS-08
TS-01, TS-02, TS-03, TS-06  <-- fix emm_basis.mm_lmm: accept mode= / lmer.df= in signature
TS-04, TS-05         <-- fix emm_basis.mm_lmm: derive nbasis from rank null-space
TS-07                <-- independent; documentation / dispatch only
```

### Upstream fixtures required
- TS-04 / TS-05: If the null-space basis for a rank-deficient fit must come from the Rust engine, an upstream fixture exposing the design null space (or the column rank + aliased column indices) is needed. This can alternatively be computed in R from `model.matrix(fit, type="fixed")` column rank — no Rust change required for the R-only path.
- TS-08 / TS-09 / TS-11: No upstream fixture needed. The reference-grid bug is in `mm_grid.mm_lmm` / `mm_group_basis` (R-side logic in `R/marginal.R`).

### Parity tolerance reference (PRD §11)
| Quantity | Tolerance |
|---|---|
| fixef estimate | 1e-4 |
| logLik | 1e-3 |
| sigma | 1e-4 |
| theta | 1e-3 |
| emmeans estimate (this spec) | 1e-2 (relaxed; KR vs Satt weighting) |
| df | 1.0 or sign (finite vs Inf) |

The emmeans estimate tolerance of 1e-2 (rather than the 1e-4 fixef tolerance) is justified because marginal means are averages over a reference grid; small differences in the averaging weights (equal vs proportional, KR correction) can introduce sub-percent-level differences that are scientifically immaterial.
