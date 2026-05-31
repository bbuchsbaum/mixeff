# Test Specification — Random-Effects Formula Syntax

**Family:** Random-effects formula syntax (`(x|g)`, `(x||g)`, nesting `/`,
interaction `:`, crossing `*`, grouping, `dummy()`, `0+`/`-1` intercept
handling, multiple grouping factors, formula-utility surface).

**Date:** 2026-05-31
**Gap report:** `assessment/gap/formula-re-syntax.md`
**Classification scope:** only `in-scope-missing`, `partial`, and `test-gap`
items are specced here; `works` items are already covered or not a gap;
`out-of-scope-by-design` items are excluded per PRD §3.

---

## How to read this document

Each spec has:

- **Name** — the `test_that()` description.
- **Kind** — one of: `unit`, `parity-vs-lme4`, `integration`, `snapshot`,
  `error-message`, `numerical-tolerance`, `cross-session-revival`,
  `performance`.
- **File** — where the test should live (new file unless noted).
- **Dataset / formula** — exact inputs.
- **Assertion** — what `expect_*` calls must pass, with numeric tolerances
  where applicable (PRD §11: fixef 1e-4, theta 1e-3, logLik 1e-3, sigma 1e-4).
- **Guards against** — the regression scenario this catches.
- **Priority** — P0 (release blocker), P1 (high), P2 (nice-to-have).
- **Upstream fixture / engine change needed** — if yes, notes what the Rust
  engine or upstream crate must expose before the test can go green.

---

## Spec TS-RE-01 — `factor()` on LHS: clear error + documented workaround

**Priority:** P1

**Kind:** error-message

**File:** `tests/testthat/test-parse-formula.R`

**Dataset / formula:**
```r
# No real dataset required for the parse-level check.
mm_parse_formula("y ~ x + (factor(grp) | g)")
```

**Assertion:**

```r
test_that("mm_parse_formula raises mm_formula_error for factor() on LHS of |", {
  expect_error(
    mm_parse_formula("y ~ x + (factor(grp) | g)"),
    class = "mm_formula_error"
  )
  # The error message must name the construct so the user knows what was rejected.
  err <- tryCatch(
    mm_parse_formula("y ~ x + (factor(grp) | g)"),
    mm_formula_error = function(e) e
  )
  expect_match(conditionMessage(err), "factor", ignore.case = TRUE)
})
```

**Companion integration test** (end-to-end workaround must fit):

```r
test_that("precomputed indicator columns substitute for factor() on LHS", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  # Synthetic grouping variable with 3 levels.
  set.seed(42L)
  sleepstudy$cond <- factor(sample(c("a", "b", "c"), nrow(sleepstudy), replace = TRUE))
  # Expand to indicator columns (the documented workaround).
  ind <- model.matrix(~ 0 + cond, data = sleepstudy)
  sleepstudy$cond_b <- ind[, "condb"]
  sleepstudy$cond_c <- ind[, "condc"]

  fit <- lmm(
    Reaction ~ Days + (0 + cond_b + cond_c | Subject),
    data = sleepstudy,
    control = mm_control(verbose = -1)
  )
  expect_s3_class(fit, "mm_lmm")
  expect_true(is.finite(logLik(fit)))
  expect_equal(length(fit$theta), 3L)  # 2x2 diagonal has 3 Cholesky params
})
```

**Guards against:** silent acceptance of `factor()` inside RE terms (which the
engine cannot handle stateless), and against the workaround path silently
breaking (e.g. if column-name handling regresses).

---

## Spec TS-RE-02 — `dummy()` on LHS: clear error + documented workaround

**Priority:** P1

**Kind:** error-message

**File:** `tests/testthat/test-parse-formula.R`

**Dataset / formula:**
```r
mm_parse_formula("distance ~ age + (0 + dummy(Sex, 'Female') | Subject)")
```

**Assertion:**

```r
test_that("mm_parse_formula raises mm_formula_error for dummy() inside RE term", {
  expect_error(
    mm_parse_formula("distance ~ age + (0 + dummy(Sex, 'Female') | Subject)"),
    class = "mm_formula_error"
  )
  err <- tryCatch(
    mm_parse_formula("distance ~ age + (0 + dummy(Sex, 'Female') | Subject)"),
    mm_formula_error = function(e) e
  )
  expect_match(conditionMessage(err), "dummy", ignore.case = TRUE)
})
```

**Companion integration test** (workaround must fit):

```r
test_that("precomputed indicator column substitutes for dummy() inside RE term", {
  skip_if_not_installed("lme4")
  data("Orthodont", package = "nlme")
  Orthodont$female <- as.integer(Orthodont$Sex == "Female")

  fit <- lmm(
    distance ~ age + (0 + female | Subject),
    data = Orthodont,
    control = mm_control(verbose = -1)
  )
  expect_s3_class(fit, "mm_lmm")
  expect_true(is.finite(logLik(fit)))
  expect_length(fit$theta, 1L)  # scalar random slope, no intercept
})
```

**Guards against:** `dummy()` being silently accepted or the workaround path
regressing.

---

## Spec TS-RE-03 — `us()` tag is an alias for bare `(x|g)` syntax

**Priority:** P2

**Kind:** error-message + integration

**File:** `tests/testthat/test-parse-formula.R`

**Rationale:** `us(x|g)` is a minor spelling; the model it requests is fully
supported via `(x|g)`. The error must be informative — it should tell the user
to use `(x|g)` instead.

**Assertion:**

```r
test_that("mm_parse_formula raises mm_formula_error for us() tag with helpful message", {
  err <- tryCatch(
    mm_parse_formula("y ~ x + (us(Days | Subject))"),
    mm_formula_error = function(e) e
  )
  expect_s3_class(err, "mm_formula_error")
  # Message must point toward the bare-bar synonym.
  expect_match(
    conditionMessage(err),
    "Days.*Subject|\\(.*\\|.*\\)",
    perl = TRUE,
    ignore.case = TRUE
  )
})
```

**Guards against:** `us()` being silently accepted in a future refactor or the
error message becoming generic/unhelpful.

---

## Spec TS-RE-04 — `diag()` tag: clear error + `||` synonym documented

**Priority:** P1

**Kind:** error-message

**File:** `tests/testthat/test-parse-formula.R`

**Assertion:**

```r
test_that("mm_parse_formula raises mm_formula_error for diag() tag", {
  expect_error(
    mm_parse_formula("y ~ x + diag(1 + a + b | g)"),
    class = "mm_formula_error"
  )
  err <- tryCatch(
    mm_parse_formula("y ~ x + diag(1 + a + b | g)"),
    mm_formula_error = function(e) e
  )
  # Error must mention diag and ideally mention || as the synonym.
  expect_match(conditionMessage(err), "diag", ignore.case = TRUE)
})
```

**Companion test** — the `||` synonym must produce the same covariance
structure (diagonal) and must fit:

```r
test_that("(a + b || g) fits and produces a diagonal covariance", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  set.seed(1L)
  sleepstudy$x2 <- rnorm(nrow(sleepstudy))

  fit <- lmm(
    Reaction ~ Days + x2 + (Days + x2 || Subject),
    data = sleepstudy,
    control = mm_control(verbose = -1)
  )
  expect_s3_class(fit, "mm_lmm")
  expect_true(is.finite(logLik(fit)))
  vc <- VarCorr(fit)
  # Off-diagonal correlations must be zero (diagonal structure).
  corr_rows <- vc$table[!is.na(vc$table$correlation), ]
  if (nrow(corr_rows) > 0) {
    expect_true(all(abs(corr_rows$correlation) < 1e-10),
                info = "|| should produce zero correlations")
  }
})
```

**Guards against:** `diag()` being accepted silently; `||` covariance structure
silently acquiring non-zero correlations.

---

## Spec TS-RE-05 — `cs()` tag: in-scope-missing, clear error, no silent path

**Priority:** P0

**Kind:** error-message

**File:** `tests/testthat/test-parse-formula.R`

**Rationale:** `cs()` is the only covariance family that is neither reachable
via an alternative syntax nor deferred by PRD §3. A user wanting
compound-symmetric covariance has no path and must receive a clear, actionable
error — not a generic parse failure.

**Assertion:**

```r
test_that("mm_parse_formula raises mm_formula_error for cs() tag (unsupported family)", {
  err <- tryCatch(
    mm_parse_formula("y ~ x + cs(1 + x | g)"),
    mm_formula_error = function(e) e
  )
  expect_s3_class(err, "mm_formula_error")
  # The message should name cs() and indicate it is not currently supported.
  expect_match(conditionMessage(err), "cs", ignore.case = TRUE)
})

test_that("lmm() raises mm_formula_error for cs() in model formula", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  expect_error(
    lmm(Reaction ~ Days + cs(Days | Subject),
        data = sleepstudy,
        control = mm_control(verbose = -1)),
    class = "mm_formula_error"
  )
})
```

**Guards against:** `cs()` being silently ignored (fitting unstructured instead)
or throwing an opaque internal error with no actionable guidance.

**Note:** When `cs()` support is added (upgrading this gap to `works`), these
error tests must be replaced by parity tests against lme4 `cs()` fits.

---

## Spec TS-RE-06 — `nobars()` R utility is exported and correct

**Priority:** P1

**Kind:** unit

**File:** `tests/testthat/test-parse-formula.R` (or a new
`test-formula-utilities.R`)

**Rationale:** `nobars()` is the most commonly needed formula-utility by
downstream tooling (emmeans, lmerTest shims, model-frame construction). Its
absence is a major ergonomic gap.

**Assertion:**

```r
test_that("mm_nobars() removes random-effect terms from a formula", {
  expect_identical(
    mm_nobars(Reaction ~ Days + (Days | Subject)),
    Reaction ~ Days
  )
  expect_identical(
    mm_nobars(y ~ 1 + x + (1 | g) + (1 | h)),
    y ~ 1 + x
  )
  expect_identical(
    mm_nobars(y ~ x + (1 | a/b)),
    y ~ x
  )
})

test_that("mm_nobars() on a formula with no RE terms returns the formula unchanged", {
  expect_identical(
    mm_nobars(y ~ x + z),
    y ~ x + z
  )
})

test_that("mm_nobars() raises an error on non-formula input", {
  expect_error(mm_nobars("y ~ x + (1|g)"), class = "mm_formula_error")
  expect_error(mm_nobars(NULL),             class = "mm_formula_error")
})
```

**Guards against:** `mm_nobars()` not being exported (namespace regression) or
returning the wrong formula after nested/interaction grouping expansion.

**Upstream fixture / engine change needed:** An R-side `mm_nobars()` wrapper
must be implemented (delegating to `reformulas::nobars()` or to a new
`wrap__mm_nobars` Rust binding). No upstream crate change is required if the
R-side delegation is used.

---

## Spec TS-RE-07 — `findbars()` R utility is exported and correct

**Priority:** P1

**Kind:** unit

**File:** `tests/testthat/test-formula-utilities.R`

**Assertion:**

```r
test_that("mm_findbars() extracts all RE sub-formulae from a model formula", {
  bars <- mm_findbars(Reaction ~ Days + (Days | Subject))
  expect_length(bars, 1L)
  # The single bar expression: Days | Subject
  expect_identical(deparse(bars[[1L]]), "Days | Subject")
})

test_that("mm_findbars() returns multiple bars for cross-classified model", {
  bars <- mm_findbars(y ~ x + (1 | subject) + (1 | item))
  expect_length(bars, 2L)
})

test_that("mm_findbars() returns NULL/empty list for a formula with no RE terms", {
  bars <- mm_findbars(y ~ x + z)
  expect_true(length(bars) == 0L)
})

test_that("mm_findbars() expands (1|a/b) to two bars", {
  bars <- mm_findbars(y ~ x + (1 | a/b))
  # After slash expansion: (1|a) + (1|a:b) => two entries
  expect_length(bars, 2L)
})
```

**Guards against:** `mm_findbars()` not being exported or returning incorrect
parse results after slash/interaction expansion.

**Upstream fixture / engine change needed:** R-side implementation delegating
to `reformulas::findbars()`. No Rust change required.

---

## Spec TS-RE-08 — `isNested()` predicate is exported and correct

**Priority:** P2

**Kind:** unit

**File:** `tests/testthat/test-formula-utilities.R`

**Assertion:**

```r
test_that("mm_is_nested() returns TRUE when f1 is nested within f2", {
  skip_if_not_installed("lme4")
  data("Pastes", package = "lme4")
  # batch/cask: cask is nested within batch in the Pastes dataset
  expect_true(mm_is_nested(Pastes$cask, Pastes$batch))
})

test_that("mm_is_nested() returns FALSE for non-nested factors", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  # Days and Subject are not nested
  expect_false(mm_is_nested(sleepstudy$Days, sleepstudy$Subject))
})

test_that("mm_is_nested() raises a typed error on non-factor/non-vector input", {
  expect_error(mm_is_nested(NULL, NULL), class = "mm_formula_error")
})
```

**Guards against:** `mm_is_nested()` not being exported or computing nesting
incorrectly (e.g. using the wrong direction).

**Upstream fixture / engine change needed:** R-side implementation delegating
to `reformulas::isNested()`. No Rust change required.

---

## Spec TS-RE-09 — `expandDoubleVerts()` utility is exported and returns a formula

**Priority:** P2

**Kind:** unit

**File:** `tests/testthat/test-formula-utilities.R`

**Rationale:** The expansion behaviour already exists inside the Rust parser
(confirmed by `mm_parse_formula`), but it is not callable as a
formula→formula utility. Downstream tooling expects to call it on a formula
object.

**Assertion:**

```r
test_that("mm_expand_double_verts() expands || to split single-bar terms", {
  result <- mm_expand_double_verts(y ~ x + (a + b || g))
  # Must return a formula object.
  expect_s3_class(result, "formula")
  deparsed <- deparse(result)
  # Should contain two (or three) single-bar blocks, no remaining ||.
  expect_false(grepl("||", deparsed, fixed = TRUE))
  expect_true(grepl("(1 | g)", deparsed, fixed = TRUE) ||
              grepl("(0 + a | g)", deparsed, fixed = TRUE))
})

test_that("mm_expand_double_verts() leaves a single-bar formula unchanged", {
  f <- y ~ x + (a + b | g)
  result <- mm_expand_double_verts(f)
  expect_identical(deparse(result), deparse(f))
})
```

**Guards against:** `mm_expand_double_verts()` not being exported, or returning
a string instead of a formula.

**Upstream fixture / engine change needed:** R-side implementation. Could
delegate to `reformulas::expandDoubleVerts()`.

---

## Spec TS-RE-10 — Manifest `formula_features` declares crossing and double-bar expansions

**Priority:** P1

**Kind:** unit (manifest accuracy)

**File:** `tests/testthat/test-manifest.R`

**Rationale:** The gap report (§ Notes) notes that `mm_formula_manifest()$formula_features`
does not advertise `crossing_expansion` or `double_bar_split` even though both
work in practice. The manifest is the machine-readable capability record used
by auditors and downstream tooling; under-declaring supported features is a
documentation bug that can cause downstream tools to fall back unnecessarily.

**Assertion:**

```r
test_that("manifest formula_features includes crossing and double-bar transformations", {
  feat <- mm_formula_manifest()$formula_features
  expect_true(
    "crossing_expansion" %in% feat$transformations,
    info = "(1|a*b) => (1|a)+(1|b)+(1|a:b) expansion must be declared"
  )
  expect_true(
    "double_bar_split" %in% feat$transformations,
    info = "(x||g) => (1|g)+(0+x|g) split must be declared"
  )
})

test_that("manifest random_term_forms includes (1|a*b) full crossing", {
  feat <- mm_formula_manifest()$formula_features
  expect_true(
    "(1 | a * b)" %in% feat$random_term_forms ||
    any(grepl("\\*", feat$random_term_forms)),
    info = "Full-crossing syntax (1|a*b) must appear in random_term_forms"
  )
})
```

**Guards against:** Manifest silently dropping newly-advertised capabilities
(preventing version drift between what the engine does and what it claims).

**Upstream fixture / engine change needed:** Yes — the Rust `mm_formula_manifest()`
function in `src/rust/src/lib.rs` must add `"crossing_expansion"` and
`"double_bar_split"` to the `transformations` vector, and add `"(1 | a * b)"`
to `random_term_forms`. This is a metadata-only change with no algorithmic
consequence.

---

## Spec TS-RE-11 — `factor()` on LHS parity: precomputed workaround vs lme4

**Priority:** P1

**Kind:** parity-vs-lme4 + numerical-tolerance

**File:** `tests/testthat/test-parse-formula.R` or a new
`test-formula-re-syntax-parity.R`

**Rationale:** The gap report confirms that precomputed indicator columns fit
correctly, but no parity test exists comparing the fixef and VarCorr output to
lme4's `(factor(grp)|g)` model.

**Dataset / formula:**
```r
# lme4 reference: lmer(Reaction ~ Days + (factor(cond) | Subject), sleepstudy_aug)
# mixeff workaround: lmm(Reaction ~ Days + (0 + cond_b + cond_c | Subject), sleepstudy_aug)
# where sleepstudy_aug has synthetic cond factor with 3 levels.
```

**Assertion:**

```r
test_that("precomputed-indicator workaround matches lme4 factor() RE fit", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  set.seed(7L)
  sleepstudy$cond <- factor(sample(c("a", "b", "c"), nrow(sleepstudy),
                                   replace = TRUE))
  # lme4 reference (intercept + 2 contrasts, treatment coding)
  ref <- suppressMessages(lme4::lmer(
    Reaction ~ Days + (factor(cond) | Subject),
    data = sleepstudy, REML = TRUE
  ))
  # mixeff workaround: precompute contrast columns
  ind <- model.matrix(~ cond, data = sleepstudy)
  sleepstudy$cond_b <- ind[, "condb"]
  sleepstudy$cond_c <- ind[, "condc"]
  fit <- lmm(
    Reaction ~ Days + (1 + cond_b + cond_c | Subject),
    data = sleepstudy, REML = TRUE,
    control = mm_control(verbose = -1)
  )
  # fixef parity (intercept + Days must match; cond contrast coefs may differ
  # in name but values must align)
  expect_equal(
    unname(fixef(fit)[c("(Intercept)", "Days")]),
    unname(lme4::fixef(ref)[c("(Intercept)", "Days")]),
    tolerance = 1e-3,
    info = "fixef Intercept and Days must match lme4 factor() fit"
  )
  # logLik must match within tolerance
  expect_equal(
    as.numeric(logLik(fit)),
    as.numeric(stats::logLik(ref)),
    tolerance = 1e-3,
    info = "logLik must match lme4 factor() fit"
  )
})
```

**Guards against:** Numeric divergence silently creeping in when the workaround
is used; also documents the exact parameterisation alignment.

---

## Spec TS-RE-12 — `||` with factor predictor: correct diagonal (better than lme4 default)

**Priority:** P1

**Kind:** parity-vs-lme4 + integration

**File:** `tests/testthat/test-formula-re-syntax-parity.R`

**Rationale:** The gap report notes that mixeff produces the *correct*
zero-correlation result for `(cond||Subject)` where lme4's default `"split"`
method gets it wrong. This deserves a regression test to ensure the behaviour
is not accidentally degraded.

**Dataset / formula:**
```r
# sleepstudy augmented with 3-level cond factor
# Formula: Reaction ~ Days + (cond || Subject)
```

**Assertion:**

```r
test_that("(factor || g) produces zero correlations (diagonal) for each contrast", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  set.seed(42L)
  sleepstudy$cond <- factor(sample(c("a", "b", "c"), nrow(sleepstudy),
                                   replace = TRUE))
  fit <- lmm(
    Reaction ~ Days + (cond || Subject),
    data = sleepstudy, REML = TRUE,
    control = mm_control(verbose = -1)
  )
  expect_s3_class(fit, "mm_lmm")
  expect_true(is.finite(logLik(fit)))
  # All off-diagonal correlations in VarCorr must be zero.
  vc <- VarCorr(fit)
  corr_rows <- vc$table[!is.na(vc$table$correlation), ]
  if (nrow(corr_rows) > 0) {
    expect_true(
      all(abs(corr_rows$correlation) < 1e-10),
      info = "(cond||Subject) must produce zero off-diagonal correlations"
    )
  }
})
```

**Guards against:** A future refactor introducing spurious correlations in the
double-bar expansion path for factor predictors.

---

## Spec TS-RE-13 — `subbars()` utility: exported and bar→plus substitution correct

**Priority:** P2

**Kind:** unit

**File:** `tests/testthat/test-formula-utilities.R`

**Assertion:**

```r
test_that("mm_subbars() replaces | and || with + producing a model.frame-safe formula", {
  result <- mm_subbars(Reaction ~ Days + (Days | Subject))
  expect_s3_class(result, "formula")
  deparsed <- deparse(result)
  expect_false(grepl("|", deparsed, fixed = TRUE))
  # model.frame() must not error on the result.
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  expect_silent(model.frame(result, data = sleepstudy))
})
```

**Guards against:** `mm_subbars()` not being exported or producing a formula
that still contains unresolved bar operators.

**Upstream fixture / engine change needed:** R-side implementation delegating
to `reformulas::subbars()`.

---

## Spec TS-RE-14 — Interaction-grouping ordering: `(1|a:b)` vs lme4 `b:a`

**Priority:** P2

**Kind:** snapshot + integration

**File:** `tests/testthat/test-formula-re-syntax-parity.R`

**Rationale:** The gap report notes that the nesting expansion `(1|a/b)` renders
as `(1|a:b)` in mixeff (correct a:b ordering) whereas lme4 renders `b:a`.
This is semantically equivalent (factor combination is commutative) but the
difference can confuse users comparing `VarCorr` printouts and ranef row names.
The test documents the divergence as expected.

**Assertion:**

```r
test_that("(1|a/b) expansion uses a:b ordering; lme4 uses b:a — documented divergence", {
  canon <- mm_parse_formula("y ~ x + (1 | a/b)")
  # mixeff canonical form: (1 | a) + (1 | a:b)
  expect_true(grepl("a:b", canon, fixed = TRUE))
  expect_false(grepl("b:a", canon, fixed = TRUE))

  skip_if_not_installed("lme4")
  data("Pastes", package = "lme4")
  fit <- lmm(
    strength ~ (1 | batch/cask),
    data = Pastes, REML = TRUE,
    control = mm_control(verbose = -1)
  )
  ref <- suppressMessages(
    lme4::lmer(strength ~ (1 | batch/cask), data = Pastes, REML = TRUE)
  )
  # logLik must match despite ordering difference in factor labels
  expect_equal(
    as.numeric(logLik(fit)),
    as.numeric(stats::logLik(ref)),
    tolerance = 1e-3,
    info = "logLik must match lme4 despite a:b vs b:a ranef label ordering"
  )
  # fixef must match
  expect_equal(
    unname(fixef(fit)),
    unname(lme4::fixef(ref)),
    tolerance = 1e-4
  )
})
```

**Guards against:** Label ordering change becoming a semantic divergence
(i.e. the factors being combined incorrectly rather than just named differently).

---

## Summary table

| ID | Name (abbreviated) | Classification | Kind | Priority |
|----|---------------------|---------------|------|----------|
| TS-RE-01 | `factor()` on LHS: error + workaround fits | partial | error-message + integration | P1 |
| TS-RE-02 | `dummy()` on LHS: error + workaround fits | partial | error-message + integration | P1 |
| TS-RE-03 | `us()` tag: informative error citing bare-bar synonym | partial | error-message | P2 |
| TS-RE-04 | `diag()` tag: error + `||` synonym fits diagonal | partial | error-message + integration | P1 |
| TS-RE-05 | `cs()` tag: clear error, no silent path | in-scope-missing | error-message | **P0** |
| TS-RE-06 | `mm_nobars()` exported and correct | in-scope-missing | unit | P1 |
| TS-RE-07 | `mm_findbars()` exported and correct | in-scope-missing | unit | P1 |
| TS-RE-08 | `mm_is_nested()` exported and correct | in-scope-missing | unit | P2 |
| TS-RE-09 | `mm_expand_double_verts()` exported and returns formula | partial | unit | P2 |
| TS-RE-10 | Manifest declares crossing + double-bar expansions | test-gap (manifest) | unit | P1 |
| TS-RE-11 | `factor()` workaround parity vs lme4 (numerical) | partial | parity-vs-lme4 + numerical-tolerance | P1 |
| TS-RE-12 | `||` with factor: correct diagonal (regression guard) | works (guard) | parity-vs-lme4 + integration | P1 |
| TS-RE-13 | `mm_subbars()` exported and model.frame-safe | in-scope-missing | unit | P2 |
| TS-RE-14 | `(1|a/b)` ordering divergence vs lme4 documented | test-gap | snapshot + integration | P2 |

**Total specs: 14** (some specs contain companion sub-tests; counting by logical
gap covered).

**Top priority:** TS-RE-05 — `cs()` (compound-symmetric) is the only
covariance family that is in-scope-missing, not reachable via any synonym, and
not deferred by PRD §3. A user requesting `cs()` today receives a generic
parse error with no actionable guidance. The spec requires two tests
(parse-level and end-to-end via `lmm()`) and no upstream engine change — only
a better error message with explicit acknowledgment that CS covariance is not
yet implemented and guidance on the available alternatives (`full` via `(x|g)`
or `diag` via `(x||g)`).

---

## Upstream changes required

| Spec | What is needed | Where |
|------|---------------|-------|
| TS-RE-06 | R-side `mm_nobars()` export, delegating to `reformulas::nobars()` | `R/parse-formula.R` + `NAMESPACE` |
| TS-RE-07 | R-side `mm_findbars()` export, delegating to `reformulas::findbars()` | `R/parse-formula.R` + `NAMESPACE` |
| TS-RE-08 | R-side `mm_is_nested()` export, delegating to `reformulas::isNested()` | `R/parse-formula.R` + `NAMESPACE` |
| TS-RE-09 | R-side `mm_expand_double_verts()` export | `R/parse-formula.R` + `NAMESPACE` |
| TS-RE-10 | Add `"crossing_expansion"` and `"double_bar_split"` to `transformations`; add `"(1 \| a * b)"` to `random_term_forms` in `mm_formula_manifest()` | `src/rust/src/lib.rs` (metadata only) |
| TS-RE-13 | R-side `mm_subbars()` export, delegating to `reformulas::subbars()` | `R/parse-formula.R` + `NAMESPACE` |

All other specs require only R-side test code; no upstream `mixeff-rs` crate
changes are needed.
