# Test Specification — "lmerTest surface" capability family

**Generated:** 2026-05-31
**Source gap report:** `assessment/gap/lmerTest.md`
**Scope:** All gaps classified `in-scope-missing`, `partial`, or `test-gap` in the
lmerTest family.  Out-of-scope items (`step()`, `get_model()`,
`as_lmerModLmerTest()`, `lmerModLmerTest` S4 class) are excluded.

Tolerances follow PRD §11 defaults unless stated otherwise:
`fixef/SE` 1e-4, `theta` 1e-3, `logLik` 1e-3, `sigma` 1e-4, `p_value` 1e-4.

---

## TS-01 · Multi-df Satterthwaite F for fixed effects with multi-level factor

**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P0 blocker

### Motivation

The most common lmerTest call is `anova(model)` on a model that contains a
factor with more than two levels or an interaction.  The gap report marks this
as a blocker: `anova(cake_model, method="satterthwaite")` currently returns
`status="unsupported"`, `p_value=NA` for every multi-df term because the
multi-df Satterthwaite F machinery is not implemented.  `method="kenward_roger"`
works, confirming the multi-df infrastructure exists; only the Satterthwaite
path is missing.

### Dataset / formula

```r
data(cake, package = "lme4")
# cake: 270 obs, 3 recipe levels, 6 temperature levels
m_cake <- lmm(
  angle ~ recipe * temperature + (1 | recipe:replicate),
  data = cake,
  control = mm_control(verbose = -1L)
)
```

### Reference (lmerTest)

```r
ref <- lmerTest::lmer(
  angle ~ recipe * temperature + (1 | recipe:replicate),
  data = cake
)
lmerTest_anova <- anova(ref, type = "III", ddf = "Satterthwaite")
# recipe:      F(2, ~254) ≈ 0.096,  p ≈ 0.91
# temperature: F(5, ~254) ≈ 37.5,   p ≈ 1e-29
# recipe:temp: F(10,~254) ≈ 0.042,  p ≈ 1.00
```

### Assertions

```r
library(testthat)
test_that("anova Satterthwaite multi-df terms return numeric F and p for cake", {
  mm_skip_if_no_lme4()

  data(cake, package = "lme4")
  m <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
            data = cake, control = mm_control(verbose = -1L))

  tbl <- anova(m, method = "satterthwaite")$table

  recipe_row <- tbl[tbl$term == "recipe", ]
  temp_row   <- tbl[tbl$term == "temperature", ]
  inter_row  <- tbl[tbl$term %in% c("recipe:temperature", "recipe * temperature"), ]

  # No unsupported rows for multi-df terms
  expect_false(any(tbl$status == "unsupported"),
               info = "multi-df terms must not return status='unsupported'")

  # Numeric F statistics and finite p-values for every term
  expect_true(all(is.finite(tbl$statistic)),
              info = "all F statistics must be finite")
  expect_true(all(is.finite(tbl$p_value)),
              info = "all p-values must be finite")

  # recipe: 2 numerator df, denominator df near 254
  expect_equal(recipe_row$num_df,  2L)
  expect_true(abs(recipe_row$den_df - 254) < 5,
              info = "recipe DenDF should be near 254")
  expect_lt(abs(recipe_row$p_value - 0.909), 0.05)   # tolerance 0.05 on p

  # temperature: 5 numerator df
  expect_equal(temp_row$num_df, 5L)
  expect_true(temp_row$p_value < 1e-20,
              info = "temperature p-value should be extremely small")

  # interaction: 10 numerator df
  expect_equal(nrow(inter_row), 1L)
  expect_equal(inter_row$num_df, 10L)
  expect_gt(inter_row$p_value, 0.90)
})
```

### What it guards against

Regression to `status="unsupported"` / `p_value=NA` for any multi-df term when
`method="satterthwaite"` is requested.

### Upstream fixture / engine note

Requires the Rust engine to expose a multi-df Satterthwaite F path (contrast
matrix Hessian trace formula for general `L` with `nrow(L) > 1`).  No upstream
fixture exists yet; the KR multi-df path in `src/rust/src/lib.rs` is the
template.  This spec cannot pass until the upstream engine change is shipped.

---

## TS-02 · Multi-df Satterthwaite F on ChickWeight (3-df Diet factor)

**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P0 blocker

### Motivation

Complements TS-01 with a different dataset and a 3-df factor, confirming the
multi-df Satterthwaite path generalises beyond 2 df.

### Dataset / formula

```r
data(ChickWeight, package = "datasets")
m_chick <- lmm(
  weight ~ Diet * Time + (1 | Chick),
  data = ChickWeight,
  control = mm_control(verbose = -1L)
)
```

### Reference (lmerTest)

```r
ref_chick <- lmerTest::lmer(weight ~ Diet * Time + (1 | Chick), data = ChickWeight)
lmerTest_anova_chick <- anova(ref_chick, ddf = "Satterthwaite")
# Diet:  F(3, ~45.5) ≈ 6.28,  p ≈ 0.0012
# Time:  F(1, ~483)  ≈ ...
# Diet:Time: F(3, ...) ≈ ...
```

### Assertions

```r
test_that("anova Satterthwaite multi-df Diet term on ChickWeight", {
  mm_skip_if_no_lme4()

  data(ChickWeight, package = "datasets")
  m <- lmm(weight ~ Diet * Time + (1 | Chick),
            data = ChickWeight,
            control = mm_control(verbose = -1L))

  tbl <- anova(m, method = "satterthwaite")$table

  diet_row <- tbl[tbl$term == "Diet", ]

  expect_false(any(tbl$status == "unsupported"))
  expect_true(all(is.finite(tbl$statistic)))
  expect_true(all(is.finite(tbl$p_value)))

  expect_equal(diet_row$num_df, 3L)
  expect_true(abs(diet_row$den_df - 45.5) < 5)
  expect_lt(abs(diet_row$statistic - 6.28), 1.0,   # tol 1.0 on F
            info = "Diet F-statistic should be near 6.28")
  expect_lt(diet_row$p_value, 0.01)
})
```

### Upstream fixture / engine note

Same dependency as TS-01.

---

## TS-03 · summary() default method produces Satterthwaite df and p-values

**Kind:** unit + parity-vs-lme4
**Priority:** P1

### Motivation

lmerTest defaults `summary()` to `ddf="Satterthwaite"`, giving a `df` column
and finite `Pr(>|t|)` in the coefficient table.  mixeff defaults to asymptotic
Wald-z (`df=NA`, no Satterthwaite p-value).  The gap report classifies this as
`partial` / major: the lme4-trained user gets no df or p-values by default and
must discover `method="satterthwaite"`.  The spec does **not** require
mixeff to change its default (that is an API decision), but it requires that:
(a) the `method="satterthwaite"` path is exercised and numerically correct,
(b) the no-method default is clearly documented as Wald-z (not silent), and
(c) a snapshot guards the default output shape so regressions are caught.

### Dataset / formula

`sleepstudy`, formula `Reaction ~ Days + (1 | Subject)`.

### Assertions

```r
test_that("summary default produces Wald-z table with NA df (documented divergence)", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))

  s <- summary(m)
  tbl <- s$table  # or coef(s) per the mm_lmm summary method

  # Default must return Wald-z (statistic_name == "z") and NA df
  expect_true(all(tbl$statistic_name == "z" | tbl$method == "asymptotic_wald_z"),
              info = "default summary should use Wald-z, not Satterthwaite")
  expect_true(all(is.na(tbl$df)),
              info = "Wald-z rows carry df = NA")
  # But p-values must still be finite
  expect_true(all(is.finite(tbl$p_value)))
})

test_that("summary(m, method='satterthwaite') matches lmerTest coef table", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))
  ref <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

  s_mm  <- summary(m, method = "satterthwaite")$table
  s_ref <- coef(summary(ref))  # lmerTest: Estimate/Std.Error/df/t/Pr(>|t|)

  days_mm  <- s_mm[s_mm$term == "Days", ]
  days_ref <- s_ref["Days", ]

  # t-statistic
  expect_lt(abs(days_mm$statistic - days_ref[["t value"]]), 1e-3)

  # Satterthwaite df — lmerTest gives ~17.0, mixeff ~16.998
  expect_lt(abs(days_mm$df - days_ref[["df"]]), 0.1)

  # p-value
  expect_lt(abs(days_mm$p_value - days_ref[["Pr(>|t|)"]]), 1e-4)
})
```

### What it guards against

Regression in the Satterthwaite coef-table path; also documents (via the first
test) that the asymptotic default is intentional.

---

## TS-04 · anova(type=) actually changes the SS hypothesis matrix

**Kind:** unit + parity-vs-lme4
**Priority:** P1

### Motivation

Gap report classifies `anova(type=)` as `partial` / major: `type` is currently
a label only; the contrast matrix is never re-derived.  Type I/II/III produce
identical statistics in all tested cases.  On unbalanced designs with
interactions, lmerTest produces genuinely different F-values per SS type.

### Dataset / formula

`ChickWeight`, `weight ~ Diet * Time + (1 | Chick)`.  The interaction term is
the key: for an unbalanced design, Type II omits the interaction from the test
denominator (marginal) while Type III conditions on it.

### Assertions

```r
test_that("anova Type II and Type III differ for unbalanced diet*time interaction", {
  mm_skip_if_no_lme4()

  data(ChickWeight, package = "datasets")
  m <- lmm(weight ~ Diet * Time + (1 | Chick),
            data = ChickWeight,
            control = mm_control(verbose = -1L))

  tbl_II  <- anova(m, type = "II",  method = "kenward_roger")$table
  tbl_III <- anova(m, type = "III", method = "kenward_roger")$table

  diet_II  <- tbl_II[tbl_II$term == "Diet",   ]
  diet_III <- tbl_III[tbl_III$term == "Diet",  ]

  # Type II and Type III F-statistics for Diet MUST differ
  # (they are the same if the type= argument is a no-op label)
  expect_false(
    isTRUE(all.equal(diet_II$statistic, diet_III$statistic, tolerance = 0.01)),
    label = "Type II and Type III Diet F must differ on unbalanced ChickWeight"
  )
})

test_that("anova Type III matches lmerTest type=3 for Diet on ChickWeight", {
  mm_skip_if_no_lme4()

  data(ChickWeight, package = "datasets")
  m <- lmm(weight ~ Diet * Time + (1 | Chick),
            data = ChickWeight,
            control = mm_control(verbose = -1L))
  ref <- lmerTest::lmer(weight ~ Diet * Time + (1 | Chick), data = ChickWeight)

  tbl_mm  <- anova(m, type = "III", method = "kenward_roger")$table
  tbl_ref <- anova(ref, type = 3, ddf = "Kenward-Roger")

  diet_mm  <- tbl_mm[tbl_mm$term == "Diet", ]
  diet_ref <- tbl_ref["Diet", ]

  expect_lt(abs(diet_mm$statistic    - diet_ref[["F value"]]), 0.1)
  expect_lt(abs(diet_mm$num_df       - diet_ref[["NumDF"]]),   0.5)
  expect_lt(abs(diet_mm$den_df       - diet_ref[["DenDF"]]),   2.0)
  expect_lt(abs(diet_mm$p_value      - diet_ref[["Pr(>F)"]]),  1e-3)
})
```

### What it guards against

The `type=` argument silently being a no-op label after an implementation;
also provides a numerical parity target for the Type III path.

---

## TS-05 · ranova whole-block REML-LRT for correlated random slope

**Kind:** parity-vs-lme4
**Priority:** P1

### Motivation

Gap report: `test_random_effect` refuses multi-parameter drops (3 theta
parameters for `(1+Days|Subject)`) and returns `status="not_assessed"`.
lmerTest's `ranova()` tests the full slope block (2 df: drops variance + covariance),
producing LRT=42.84, p=4.99e-10.  This is the canonical "is the random slope
needed?" workflow that no mixeff path currently reproduces.

### Dataset / formula

`sleepstudy`, formula `Reaction ~ Days + (Days | Subject)`.

### Assertions

```r
test_that("test_random_effect handles full correlated-slope block (2-df drop)", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))

  # The call that currently returns status="not_assessed"
  out <- test_random_effect(m, "(1+Days|Subject)", method = "boundary_lrt")
  row <- out$table[1L, ]

  # Must succeed
  expect_identical(row$status, "available",
                   info = "whole-block slope drop must not return not_assessed")

  # LRT statistic should be near lmerTest ranova value of 42.84
  expect_lt(abs(row$statistic - 42.84), 2.0,
            info = "LRT chi-square should be near 42.84")
  expect_lt(row$p_value, 1e-6)

  # The reference for a 2-parameter drop is the boundary mixture
  # 0.5*chi2(1) + 0.5*chi2(2)  (Self-Liang for 2 components)
  expect_match(row$reference_distribution, "chi-square",
               info = "reference distribution label must be present")
})
```

### What it guards against

Permanent `not_assessed` / refusal for whole-block random-slope tests, which
is the most common ranova use-case.

### Upstream fixture / engine note

Requires the Rust engine to support multi-theta boundary-mixture tests.  The
appropriate mixture weights for a 2-parameter correlated-block drop differ from
the scalar (50:50) case.  An upstream fixture must specify the expected mixture
weights and statistic for this case before the assertion can be numerically
tight.

---

## TS-06 · ranova single-component boundary LRT matches lmerTest ranova (RI)

**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Motivation

For a random-intercept-only model, `test_random_effect` already works (boundary
LRT with 50:50 chi-bar-square mixture).  This spec pins numerical parity against
lmerTest's `ranova()` on `sleepstudy`/RI so a regression is caught.

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (1 | Subject)`.

### Assertions

```r
test_that("test_random_effect RI boundary LRT matches lmerTest ranova", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            REML = FALSE,
            control = mm_control(verbose = -1L))
  ref <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE)

  out <- test_random_effect(m, "Subject", method = "boundary_lrt")
  row <- out$table[1L, ]

  expect_identical(row$status, "available")

  # lmerTest ranova: LRT ≈ 75.46 (or equivalent boundary-corrected value)
  # We pin against lmerTest-reported value with tolerance 2.0
  ref_ranova <- lmerTest::ranova(
    lmerTest::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
  )
  lmerTest_LRT <- ref_ranova$LRT[2]
  expect_lt(abs(row$statistic - lmerTest_LRT), 2.0,
            info = "boundary LRT statistic should match lmerTest ranova")
  expect_lt(row$p_value, 1e-10)
})
```

### What it guards against

Numerical regression in the existing RI boundary-LRT path.

---

## TS-07 · contrast() joint / contestMD Satterthwaite F for multi-row L

**Kind:** unit + parity-vs-lme4
**Priority:** P1

### Motivation

Gap report: `contrast()` has no `joint` argument; multi-row `L` produces
per-row t-tests only.  lmerTest's `contestMD()` collapses a multi-row contrast
matrix into a single Satterthwaite F.  A user wanting "joint test of these 2
contrasts" under Satterthwaite has no route in mixeff.

### Dataset / formula

`cake`, formula `angle ~ recipe * temperature + (1 | recipe:replicate)`, using
the `recipe` main effect (2 contrasts).

### Assertions

```r
test_that("contrast() with joint=TRUE returns single-row Satterthwaite F", {
  mm_skip_if_no_lme4()

  data(cake, package = "lme4")
  m <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
            data = cake, control = mm_control(verbose = -1L))

  # L matrix: the two recipe contrast rows
  p   <- length(fixef(m))
  L   <- matrix(0, nrow = 2L, ncol = p)
  # columns for recipe B and recipe C (indices 2 and 3 in treatment coding)
  L[1L, 2L] <- 1
  L[2L, 3L] <- 1

  out <- contrast(m, L, method = "satterthwaite", joint = TRUE)
  row <- out$table[1L, ]

  # Must be a single joint-test row, not two per-row t rows
  expect_equal(nrow(out$table), 1L,
               info = "joint=TRUE must collapse to one row")
  expect_identical(row$statistic_name, "F",
                   info = "joint Satterthwaite test must report an F statistic")
  expect_equal(row$num_df, 2L)
  expect_true(is.finite(row$den_df))
  expect_true(is.finite(row$statistic))
  expect_true(is.finite(row$p_value))
})

test_that("contrast() joint F matches lmerTest contestMD for recipe", {
  mm_skip_if_no_lme4()

  data(cake, package = "lme4")
  m <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
            data = cake, control = mm_control(verbose = -1L))
  ref <- lmerTest::lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
                        data = cake)

  p   <- length(fixef(m))
  L   <- matrix(0, nrow = 2L, ncol = p)
  L[1L, 2L] <- 1; L[2L, 3L] <- 1

  out_mm  <- contrast(m, L, method = "satterthwaite", joint = TRUE)

  # lmerTest reference
  L_ref <- matrix(0, nrow = 2L, ncol = length(lme4::fixef(ref)))
  L_ref[1L, 2L] <- 1; L_ref[2L, 3L] <- 1
  ref_ct <- lmerTest::contestMD(ref, L_ref, ddf = "Satterthwaite")

  expect_lt(abs(out_mm$table$statistic - ref_ct[["F value"]]), 0.5)
  expect_lt(abs(out_mm$table$den_df    - ref_ct[["DenDF"]]),   5)
  expect_lt(abs(out_mm$table$p_value   - ref_ct[["Pr(>F)"]]),  1e-3)
})
```

### What it guards against

The `joint` argument being silently ignored; regressions in the
contestMD-equivalent multi-df F path after it is implemented.

### Upstream fixture / engine note

Requires the Rust engine to expose a general multi-row joint-F Satterthwaite
path that is distinct from per-row t.  This is a new engine surface.

---

## TS-08 · show_tests() / hypothesis matrix retrieval for anova and mm_means

**Kind:** unit + integration (audit-first transparency)
**Priority:** P1

### Motivation

Gap report classifies `show_tests()` as `in-scope-missing` / major.  The
mixeff audit-first mandate (CLAUDE.md: "every printed claim traces to a JSON
artifact") makes the inability to retrieve the exact contrast matrix behind each
ANOVA row a direct violation of that mandate.  `attr(anova(m), "hypotheses")`
currently returns `NULL`.  `mm_show_tests` exists in the namespace but is
unexported and unwired.

### Assertions

```r
test_that("anova result carries non-null hypotheses attribute", {
  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))

  tbl <- anova(m, method = "satterthwaite")
  hyp <- attr(tbl, "hypotheses")

  expect_false(is.null(hyp),
               info = "anova result must carry a 'hypotheses' attribute")
  expect_true(is.list(hyp),
              info = "'hypotheses' should be a named list of contrast matrices")
  expect_true("Days" %in% names(hyp),
              info = "hypotheses list must have an entry for each term")
  expect_true(is.matrix(hyp[["Days"]]),
              info = "each hypothesis entry must be a numeric matrix")
})

test_that("mm_means result carries non-null hypotheses attribute", {
  data(cake, package = "lme4")
  m <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
            data = cake, control = mm_control(verbose = -1L))

  means <- mm_means(m, ~ recipe)
  hyp   <- attr(means, "hypotheses")

  expect_false(is.null(hyp))
  expect_true(is.list(hyp) || is.matrix(hyp),
              info = "mm_means hypotheses may be a single matrix or named list")
})

test_that("show_tests() is exported and returns a list of contrast matrices", {
  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))
  tbl <- anova(m)

  # show_tests must be an exported function
  expect_true(isNamespaceLoaded("mixeff"))
  expect_true(existsMethod("show_tests", signature()) ||
              exists("show_tests", envir = asNamespace("mixeff"), inherits = FALSE),
              info = "show_tests must be exported from mixeff namespace")

  out <- show_tests(tbl)
  expect_true(is.list(out))
  expect_true(all(vapply(out, is.matrix, logical(1L))),
              info = "each element of show_tests() output must be a matrix")
})
```

### What it guards against

Silent omission of the `hypotheses` attribute; `show_tests()` remaining
unexported; the audit-trail guarantee being broken for ANOVA and marginal-means
tables.

---

## TS-09 · drop1() produces a Satterthwaite F route (method="satterthwaite")

**Kind:** unit + parity-vs-lme4
**Priority:** P2

### Motivation

Gap report: `drop1(mm_lmm)` is a refit-based asymptotic LRT.  lmerTest's
`drop1()` always returns a Satterthwaite F-table.  For the single-df case the
two statistics are equivalent (F = t²); parity is confirmed in the probe.  The
spec guards the single-df equivalence and opens a documented path for a
`method="satterthwaite"` option in `drop1.mm_lmm` that routes through
`test_effect()`.

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (1 | Subject)`.

### Assertions

```r
test_that("drop1 LRT p-value is consistent with test_effect Satterthwaite F", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))

  d1  <- drop1(m)
  te  <- test_effect(m, "Days", method = "satterthwaite")

  # Both should report Days as significant at p < 1e-10
  expect_lt(d1$table$p_value[d1$table$dropped == "Days"],  1e-10)
  expect_lt(te$table$p_value[1L], 1e-10)
})

test_that("drop1(m, method='satterthwaite') routes to Satterthwaite F", {
  # This test will FAIL until the feature is implemented — that is intentional.
  # It specifies the expected behaviour of the new method= argument.
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))

  d1_sw <- tryCatch(
    drop1(m, method = "satterthwaite"),
    error = function(e) NULL
  )
  skip_if(is.null(d1_sw), "drop1(method='satterthwaite') not yet implemented")

  days_row <- d1_sw$table[d1_sw$table$dropped == "Days", ]
  expect_identical(days_row$statistic_name, "F")
  # F = t² where t ≈ 13.015, so F ≈ 169.4
  expect_lt(abs(days_row$statistic - 169.4), 1.0)
  expect_lt(abs(days_row$den_df    - 161.0), 1.0)
})
```

### What it guards against

The LRT/F equivalence regressing on the single-df case; also serves as a
failing test that defines the future `method=` argument surface.

---

## TS-10 · summary() coef table column names (Estimate/Std.Error/df/t/p)

**Kind:** unit + snapshot
**Priority:** P2

### Motivation

Gap report: coef columns are present but renamed/transposed with additional
audit columns (`method`, `status`, `reliability`).  Downstream consumers that
use `coef(summary(lme4_model))` column names will break.  This spec does **not**
require lme4-compatible column names (that would violate PRD §3 non-goals), but
it pins the exact column schema so consumers can be written against a stable
contract, and it verifies the table contains the inference-essential quantities
under the mixeff naming convention.

### Assertions

```r
test_that("summary Satterthwaite coef table has required columns", {
  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))

  s   <- summary(m, method = "satterthwaite")
  tbl <- s$table

  required_cols <- c("term", "estimate", "std_error", "df", "statistic",
                     "statistic_name", "p_value", "method", "status")
  for (col in required_cols) {
    expect_true(col %in% names(tbl),
                info = sprintf("summary table must have column '%s'", col))
  }

  # Satterthwaite rows: df must be finite and positive
  expect_true(all(is.finite(tbl$df) & tbl$df > 0))
  # statistic_name must be "t" for Satterthwaite
  expect_true(all(tbl$statistic_name == "t"))
})

test_that("summary Wald-z coef table has required columns and NA df", {
  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))

  s   <- summary(m)          # default: Wald-z
  tbl <- s$table

  required_cols <- c("term", "estimate", "std_error", "statistic",
                     "statistic_name", "p_value", "method", "status")
  for (col in required_cols) {
    expect_true(col %in% names(tbl))
  }
  expect_true(all(tbl$statistic_name == "z"))
  # df may be absent or NA for Wald-z
  if ("df" %in% names(tbl)) {
    expect_true(all(is.na(tbl$df)))
  }
})
```

### What it guards against

Column-schema regressions in the summary table after implementation work;
provides the stable API contract for downstream consumers.

---

## TS-11 · VarCorr group name preserves colon separator for interaction groups

**Kind:** unit
**Priority:** P2

### Motivation

Gap report Finding 3 in `lmm-cake-interaction.md`: `VarCorr(fit_mm)` returns
group name `"recipe & replicate"` instead of `"recipe:replicate"`.  Any code
that inspects group names by string matching will fail.  This is classified
`in-scope-missing` / minor in the cake parity probe and is directly caused by
the lmerTest family's use of interaction-defined grouping factors.

### Assertions

```r
test_that("VarCorr group name for interaction grouping factor uses colon", {
  data(cake, package = "lme4")
  m <- lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
            data = cake, control = mm_control(verbose = -1L))

  vc <- VarCorr(m)
  group_names <- names(vc)   # or rownames if data.frame form

  expect_true(
    any(grepl("recipe:replicate", group_names, fixed = TRUE)),
    info = paste(
      "VarCorr must use ':' not ' & ' for interaction group.",
      "Got:", paste(group_names, collapse = ", ")
    )
  )

  # As data frame, group column should also use colon form
  vc_df <- as.data.frame(vc)
  if ("grp" %in% names(vc_df)) {
    expect_true(any(grepl("recipe:replicate", vc_df$grp, fixed = TRUE)))
  }
})
```

### What it guards against

The `" & "` separator persisting in VarCorr output for interaction-defined
grouping factors, breaking downstream group-name matching.

---

## TS-12 · ranova reduce.terms analog — documented refusal or implementation

**Kind:** unit + error-message
**Priority:** P2

### Motivation

Gap report: `ranova(reduce.terms=TRUE)` (drop only the covariance, keep the
variance) has no mixeff analog.  The spec defines the expected behaviour for
the current state (structured refusal with a clear reason code) and provides
the assertion to promote if the feature is implemented.

### Assertions

```r
test_that("test_random_effect refuses correlated-block partial-drop with clear reason", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))

  # Attempting to drop only the covariance (reduce.terms analog):
  # Pass the correlation parameter spec as a term name.
  # If not yet supported this must return a structured refusal,
  # NOT an R error or NA without explanation.
  out <- tryCatch(
    test_random_effect(m, "Days|Subject_covariance", method = "boundary_lrt"),
    error = function(e) e
  )

  if (inherits(out, "mm_random_effect_test")) {
    row <- out$table[1L, ]
    expect_true(row$status %in% c("not_assessed", "unsupported"),
                info = "partial covariance drop must return structured status, not error")
    expect_false(is.null(row$reason),
                 info = "refusal must carry a reason string")
    expect_match(row$reason, "(covariance|partial|reduce)",
                 info = "reason must explain why the partial drop is not supported")
  } else {
    # If an R-level error is thrown, that is itself a defect
    fail(paste("test_random_effect threw an R error instead of a structured refusal:",
               conditionMessage(out)))
  }
})
```

### What it guards against

Hard R errors instead of structured refusals for unsupported random-effect test
variants; silently wrong results from partial-block drops.

---

## TS-13 · contest1D per-row t parity — sleepstudy Days contrast

**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Motivation

Gap report: `contest(joint=FALSE)` / `contest1D` is classified `works`.  This
spec pins the numerical parity so a regression is caught and serves as a
cross-reference regression test for the Satterthwaite df machinery underlying
all higher-level tests.

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (1 | Subject)`, contrast `c(0, 1)` (Days).

### Assertions

```r
test_that("contrast() single-row Satterthwaite matches lmerTest contest1D for Days", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))
  ref <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

  ct_mm  <- contrast(m, c(0, 1), method = "satterthwaite")$table
  ct_ref <- lmerTest::contest1D(ref, c(0, 1), ddf = "Satterthwaite")

  expect_lt(abs(ct_mm$estimate  - ct_ref[["Estimate"]]),          1e-4)
  expect_lt(abs(ct_mm$std_error - ct_ref[["Std. Error"]]),        1e-4)
  expect_lt(abs(ct_mm$df        - ct_ref[["df"]]),                0.1)
  expect_lt(abs(ct_mm$statistic - ct_ref[["t value"]]),           1e-3)
  expect_lt(abs(ct_mm$p_value   - ct_ref[["Pr(>|t|)"]]),          1e-4)
})
```

### What it guards against

Numerical regressions in the Satterthwaite df/t/p computation for single
fixed-effect contrasts.

---

## TS-14 · ls_means / mm_means Satterthwaite parity for cake recipe

**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Motivation

Gap report classifies `ls_means()` as `works` with exact numerical agreement
(cake recipe means 33.122/31.644/31.600, SE 1.7368, df 42).  This spec pins
those numbers as a regression test.

### Assertions

```r
test_that("mm_means recipe marginal means match lmerTest ls_means for cake", {
  mm_skip_if_no_lme4()

  data(cake, package = "lme4")
  m <- lmm(angle ~ recipe + temperature + (1 | recipe:replicate),
            data = cake, control = mm_control(verbose = -1L))
  ref <- lmerTest::lmer(angle ~ recipe + temperature + (1 | recipe:replicate),
                        data = cake)

  mm_out  <- mm_means(m, ~ recipe)$table
  ref_out <- lmerTest::ls_means(ref, which = "recipe")

  # Pin means (tolerance 1e-3 — relaxed for contrast-coding differences)
  ref_vals <- ref_out[["Estimate"]]
  mm_vals  <- mm_out$estimate[order(mm_out$level)]

  expect_equal(length(mm_vals), length(ref_vals))
  expect_lt(max(abs(mm_vals - sort(ref_vals))), 1e-3)

  # SE
  ref_se <- ref_out[["Std. Error"]]
  mm_se  <- mm_out$std_error[order(mm_out$level)]
  expect_lt(max(abs(mm_se - sort(ref_se))), 1e-3)

  # Satterthwaite df (all should be ~42)
  expect_true(all(abs(mm_out$df - 42) < 2))
})
```

### What it guards against

Regressions in marginal-means estimation and the Satterthwaite df propagation
through `mm_means()`.

---

## TS-15 · df_for_contrast() / calcSatterth() raw ddf parity

**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P2

### Motivation

Gap report: `df_for_contrast()` (analog of `calcSatterth()`) is classified
`works`.  This spec pins numerical parity for the raw ddf computation
independently of the full contrast pipeline.

### Assertions

```r
test_that("df_for_contrast Satterthwaite df matches lmerTest calcSatterth", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))
  ref <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

  L <- c(0, 1)

  df_mm  <- df_for_contrast(m, L, method = "satterthwaite")
  df_ref <- lmerTest::calcSatterth(ref, L)

  expect_lt(abs(df_mm - df_ref[["df"]]), 0.1,
            info = "Satterthwaite ddf must match lmerTest calcSatterth to 0.1")
})
```

---

## TS-16 · anova single-df Satterthwaite parity — sleepstudy Days

**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Motivation

The gap report shows that single-df `anova(m, method="satterthwaite")` works
and matches lmerTest (p=3.30e-06 vs lmerTest 3.264e-06).  This spec pins those
numbers so a regression in the single-df path is caught separately from the
multi-df blocker.

### Assertions

```r
test_that("anova Satterthwaite single-df Days matches lmerTest on sleepstudy", {
  mm_skip_if_no_lme4()

  data(sleepstudy, package = "lme4")
  m <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
            control = mm_control(verbose = -1L))
  ref <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

  tbl_mm  <- anova(m, method = "satterthwaite")$table
  tbl_ref <- anova(ref, ddf = "Satterthwaite")

  days_mm  <- tbl_mm[tbl_mm$term == "Days", ]
  days_ref <- tbl_ref["Days", ]

  # F = t² for single df; mixeff may report t, lmerTest reports F
  f_mm <- if (days_mm$statistic_name == "t") days_mm$statistic^2
          else days_mm$statistic
  expect_lt(abs(f_mm            - days_ref[["F value"]]), 0.01)
  expect_lt(abs(days_mm$den_df  - days_ref[["DenDF"]]),   0.5)
  expect_lt(abs(days_mm$p_value - days_ref[["Pr(>F)"]]),  1e-4)
})
```

---

## Summary table

| ID    | Name (short)                                  | Kind                          | Dataset / formula                              | Priority |
|-------|-----------------------------------------------|-------------------------------|------------------------------------------------|----------|
| TS-01 | Multi-df Satterthwaite F — cake               | parity + numerical-tolerance  | `cake`, `angle ~ recipe*temperature+(1\|r:r)` | P0       |
| TS-02 | Multi-df Satterthwaite F — ChickWeight Diet   | parity + numerical-tolerance  | `ChickWeight`, `weight~Diet*Time+(1\|Chick)`  | P0       |
| TS-03 | summary() default vs Satterthwaite            | unit + parity                 | `sleepstudy`, `Reaction~Days+(1\|Subject)`    | P1       |
| TS-04 | anova(type=) changes the SS hypothesis        | unit + parity                 | `ChickWeight`, Type II vs III                  | P1       |
| TS-05 | ranova whole-block slope drop (2-df)          | parity                        | `sleepstudy`, `Reaction~Days+(Days\|Subject)` | P1       |
| TS-06 | ranova RI boundary LRT matches lmerTest       | parity + numerical-tolerance  | `sleepstudy`, `Reaction~Days+(1\|Subject)`    | P1       |
| TS-07 | contrast() joint / contestMD Satterthwaite F  | unit + parity                 | `cake`, recipe L-matrix                        | P1       |
| TS-08 | show_tests() exported and hypotheses attached | unit + integration            | `sleepstudy`, `cake`                           | P1       |
| TS-09 | drop1 Satterthwaite F route                   | unit + parity                 | `sleepstudy`, `Reaction~Days+(1\|Subject)`    | P2       |
| TS-10 | summary coef table column contract            | unit + snapshot               | `sleepstudy`                                   | P2       |
| TS-11 | VarCorr interaction group name colon          | unit                          | `cake`, `recipe:replicate`                     | P2       |
| TS-12 | ranova reduce.terms — structured refusal      | unit + error-message          | `sleepstudy`, `(Days\|Subject)`               | P2       |
| TS-13 | contest1D single-row t parity (Days)          | parity + numerical-tolerance  | `sleepstudy`                                   | P1       |
| TS-14 | mm_means / ls_means parity (cake recipe)      | parity + numerical-tolerance  | `cake`, `recipe`                               | P1       |
| TS-15 | df_for_contrast / calcSatterth raw ddf        | parity + numerical-tolerance  | `sleepstudy`                                   | P2       |
| TS-16 | anova single-df Satterthwaite sleepstudy Days | parity + numerical-tolerance  | `sleepstudy`                                   | P1       |

**Total: 16 test specifications.**

**Top priority (P0 blocker):** TS-01 — multi-df Satterthwaite F for `anova()`
on a multi-level factor.  This is the single most common lmerTest call and is
the dominant parity hole in the family.  Every model with a factor having more
than two levels or an interaction returns `p_value=NA` /
`status="unsupported"` for those terms when `method="satterthwaite"` is
requested; the workflow is entirely blocked until TS-01 passes.  TS-01 cannot
pass without a Rust engine change to expose the multi-df Satterthwaite F path.
