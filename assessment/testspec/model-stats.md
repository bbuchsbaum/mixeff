# Test Specification — Model-statistics generics

**Family:** Model-statistics generics — `logLik` (REML/ML), `deviance`, `AIC`,
`BIC`, `nobs`, `df.residual`, `extractAIC`, `REMLcrit`, `isREML`, `refitML`,
`devcomp`, `getME` component coverage.

**Source gap report:** `assessment/gap/model-stats.md`
**Reference:** lme4 2.0.1 / lmerTest 3.2.1
**Date authored:** 2026-05-31

Gaps classified **in-scope-missing**, **partial**, or **test-gap** in the gap
report drive every spec below. Gaps classified **works** or
**out-of-scope-by-design** are excluded.

Standard dataset used throughout unless otherwise stated:

```r
library(lme4); library(mixeff)
sleep_reml <- lmm(Reaction ~ Days + (Days|Subject), sleepstudy,
                  control = mm_control(verbose = -1))          # REML fit
sleep_ml   <- lmm(Reaction ~ Days + (Days|Subject), sleepstudy,
                  REML = FALSE, control = mm_control(verbose = -1))  # ML fit
```

lme4 reference:
```r
ref_reml <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
ref_ml   <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                       REML = FALSE)
```

---

## MS-01 — `logLik` `nall` attribute

**Name:** `logLik-nall-attribute`
**Kind:** parity-vs-lme4
**Priority:** P2
**Gap classification:** partial (minor)
**Source finding:** F1 / summary-table row "logLik `nall` attribute"

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML and ML fits.

### Assertion

```r
test_that("logLik carries nall attribute matching nobs", {
  mm_skip_if_no_lme4()
  ll <- logLik(sleep_reml)
  expect_identical(attr(ll, "nall"), nobs(sleep_reml))
  # lme4 reference value is 180
  expect_identical(attr(ll, "nall"), 180L)

  ll_ml <- logLik(sleep_ml)
  expect_identical(attr(ll_ml, "nall"), nobs(sleep_ml))
})
```

### What it guards against

`bbmle::AICtab`, `AICcmodavg`, and other downstream packages inspect `nall`
when computing small-sample-corrected AIC (AICc). A `NULL` `nall` silently
produces wrong AICc values in those packages.

### Implementation note

`logLik.mm_lmm` in `R/methods-extract.R:291-298` must add
`nall = object$nobs` to the `structure()` call.

---

## MS-02 — `logLik(m, REML=FALSE)` on a REML fit (silent wrong-value bug)

**Name:** `logLik-REML-override-REML-fit`
**Kind:** parity-vs-lme4
**Priority:** P0 (blocker)
**Gap classification:** partial (major) / test-gap
**Source finding:** F1

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.

### Assertion

```r
test_that("logLik(m, REML=FALSE) on a REML fit re-evaluates ML criterion", {
  mm_skip_if_no_lme4()
  # lme4 reference: -875.9929 (ML logLik at REML theta)
  ref_ll_ml <- logLik(ref_reml, REML = FALSE)
  mm_ll_ml  <- logLik(sleep_reml, REML = FALSE)

  # Must differ from the REML logLik
  reml_ll <- logLik(sleep_reml)
  expect_false(isTRUE(all.equal(as.numeric(mm_ll_ml),
                                as.numeric(reml_ll))))

  # Must agree with lme4 to tolerance
  expect_equal(as.numeric(mm_ll_ml), as.numeric(ref_ll_ml),
               tolerance = 1e-3,
               label = "logLik(REML fit, REML=FALSE): ML value at REML theta")

  # Class must be preserved
  expect_s3_class(mm_ll_ml, "logLik")
  expect_identical(attr(mm_ll_ml, "nobs"), nobs(sleep_reml))
})
```

### What it guards against

`logLik.mm_lmm` currently has no `REML` parameter; the argument silently falls
into `...` and is swallowed, returning the wrong (REML) log-likelihood. Any
LRT that calls `logLik(reml_fit, REML=FALSE)` gets a plausible-but-wrong
number with no diagnostic — a direct violation of the "no silent surgery"
contract in CLAUDE.md.

### Upstream fixture / engine change needed

This requires that the Rust engine (or the R wrapper) can re-evaluate the
objective function at stored θ under a different criterion flag. If the engine
does not expose a re-evaluate path, the wrapper must at minimum refuse with a
typed error (`mm_inference_unavailable`) rather than return silently wrong
output. The spec passes under either behavior: correct value (parity) OR a
`mm_inference_unavailable` condition that the test can adjust to assert.

---

## MS-03 — `logLik(m, REML=TRUE)` on an ML fit

**Name:** `logLik-REML-override-ML-fit`
**Kind:** parity-vs-lme4
**Priority:** P1
**Gap classification:** partial (major) / test-gap
**Source finding:** F1

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, ML fit.

### Assertion

```r
test_that("logLik(m, REML=TRUE) on an ML fit re-evaluates REML criterion", {
  mm_skip_if_no_lme4()
  ref_ll_reml <- logLik(ref_ml, REML = TRUE)
  mm_ll_reml  <- logLik(sleep_ml, REML = TRUE)

  # Must differ from the ML logLik
  ml_ll <- logLik(sleep_ml)
  expect_false(isTRUE(all.equal(as.numeric(mm_ll_reml),
                                as.numeric(ml_ll))))

  expect_equal(as.numeric(mm_ll_reml), as.numeric(ref_ll_reml),
               tolerance = 1e-3,
               label = "logLik(ML fit, REML=TRUE): REML value at ML theta")
  expect_s3_class(mm_ll_reml, "logLik")
})
```

### What it guards against

Mirror of MS-02 for the ML→REML direction. Same silent-wrong-value risk.

---

## MS-04 — `deviance(m)` for a REML fit returns ML deviance

**Name:** `deviance-REML-fit-returns-ML`
**Kind:** parity-vs-lme4 / numerical-tolerance
**Priority:** P1
**Gap classification:** partial (major) / test-gap
**Source finding:** F2

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.

### Assertion

```r
test_that("deviance() on a REML fit returns ML deviance, not REMLcrit", {
  mm_skip_if_no_lme4()
  # lme4 documented contract: deviance() always returns the ML deviance
  # lme4 ref: deviance(ref_reml, REML=FALSE) == 1751.986
  ref_dev <- deviance(ref_reml, REML = FALSE)   # 1751.99
  mm_dev  <- deviance(sleep_reml)

  expect_equal(mm_dev, ref_dev, tolerance = 1e-3,
               label = "deviance(REML fit): ML deviance, not REMLcrit")

  # Sanity: ML deviance != REML criterion
  remlcrit_approx <- -2 * as.numeric(logLik(sleep_reml))
  expect_false(isTRUE(all.equal(mm_dev, remlcrit_approx, tolerance = 1e-2)),
               label = "deviance() must not return REMLcrit for REML fit")
})
```

### What it guards against

The gap report shows `m$deviance` for a REML fit stores `1743.628` (the REML
criterion), not `1751.99` (the ML deviance). This means `deviance()` currently
silently returns the wrong quantity for REML fits, violating lme4's documented
contract. The test confirms the fix.

### Upstream fixture / engine change needed

The Rust artifact may need to store both the REML criterion and the ML
deviance in separate fields. The R wrapper must map `deviance.mm_lmm` to the
ML deviance field regardless of how the model was fitted.

---

## MS-05 — `deviance(m, REML=FALSE)` is honoured, not silently ignored

**Name:** `deviance-REML-arg-not-ignored`
**Kind:** parity-vs-lme4
**Priority:** P1
**Gap classification:** partial (major) / test-gap
**Source finding:** F2

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.

### Assertion

```r
test_that("deviance(m, REML=FALSE) returns ML deviance, not a silent no-op", {
  mm_skip_if_no_lme4()
  # When REML=FALSE is explicitly supplied, result must equal ML deviance
  mm_dev_ml <- deviance(sleep_reml, REML = FALSE)
  ref_dev   <- deviance(ref_reml, REML = FALSE)   # 1751.99

  expect_equal(mm_dev_ml, ref_dev, tolerance = 1e-3,
               label = "deviance(m, REML=FALSE): explicit ML deviance")
})

test_that("deviance(m, REML=TRUE) on an ML fit returns REMLcrit or typed refusal", {
  mm_skip_if_no_lme4()
  # Either: returns the REML criterion correctly, OR refuses with a typed error.
  # Silent wrong value is NOT acceptable.
  result <- tryCatch(
    deviance(sleep_ml, REML = TRUE),
    mm_inference_unavailable = function(e) "refused"
  )
  if (!identical(result, "refused")) {
    ref_remlcrit <- deviance(ref_ml, REML = TRUE)
    expect_equal(result, ref_remlcrit, tolerance = 1e-3,
                 label = "deviance(ML fit, REML=TRUE): REML criterion")
  }
})
```

### What it guards against

`deviance.mm_lmm` has no `REML` parameter; the argument is silently ignored.
This means any code doing `deviance(m, REML=FALSE)` gets the wrong answer.
The "no silent surgery" contract requires either correctness or a typed
diagnostic.

---

## MS-06 — `extractAIC` method is defined

**Name:** `extractAIC-method-exists`
**Kind:** unit / integration
**Priority:** P0 (blocker)
**Gap classification:** in-scope-missing (major)
**Source finding:** F3

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML and ML fits.
Also test via `drop1()` or `step()` to confirm base-R selection callers work.

### Assertion

```r
test_that("extractAIC.mm_lmm is defined and returns (edf, AIC) pair", {
  # extractAIC must not throw "no applicable method"
  eaic_reml <- extractAIC(sleep_reml)
  expect_length(eaic_reml, 2L)
  expect_true(is.numeric(eaic_reml))
  # edf is the effective degrees of freedom (= dof)
  expect_equal(eaic_reml[[1L]], sleep_reml$dof,
               label = "extractAIC edf == dof")
  # AIC component equals AIC()
  expect_equal(eaic_reml[[2L]], AIC(sleep_reml), tolerance = 1e-6,
               label = "extractAIC AIC component == AIC()")
})

test_that("extractAIC(m, k=log(nobs)) equals BIC up to sign convention", {
  n <- nobs(sleep_ml)
  eaic_bic <- extractAIC(sleep_ml, k = log(n))
  expect_equal(eaic_bic[[2L]], BIC(sleep_ml), tolerance = 1e-6,
               label = "extractAIC(k=log(n)) matches BIC()")
})

test_that("extractAIC allows drop1() to run without error", {
  mm_skip_if_no_lme4()
  # drop1 needs extractAIC; it must complete without dispatch error
  expect_no_error(
    drop1(sleep_ml, test = "none"),
    message = "drop1() requires extractAIC dispatch"
  )
})
```

### What it guards against

Without `extractAIC.mm_lmm`, `step()` and `drop1()` abort immediately with
"no applicable method for 'extractAIC'". This breaks the base-R model-
selection workflow that all lme4 users rely on.

---

## MS-07 — `REMLcrit` generic is defined

**Name:** `REMLcrit-method-exists`
**Kind:** unit / parity-vs-lme4
**Priority:** P1
**Gap classification:** in-scope-missing (major)
**Source finding:** F3

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.

### Assertion

```r
test_that("REMLcrit() returns the REML criterion for a REML fit", {
  mm_skip_if_no_lme4()
  rc <- REMLcrit(sleep_reml)
  ref_rc <- REMLcrit(ref_reml)   # lme4: 1743.628

  expect_true(is.numeric(rc) && length(rc) == 1L && is.finite(rc))
  expect_equal(rc, ref_rc, tolerance = 1e-3,
               label = "REMLcrit(): REML criterion matches lme4")

  # Value must equal -2 * REML logLik
  expect_equal(rc, -2 * as.numeric(logLik(sleep_reml)), tolerance = 1e-8,
               label = "REMLcrit() == -2 * logLik(REML fit)")
})

test_that("REMLcrit() on an ML fit throws a typed error", {
  expect_error(REMLcrit(sleep_ml),
               class = "mm_inference_unavailable",
               label = "REMLcrit() on ML fit must refuse, not crash")
})
```

### What it guards against

`REMLcrit` is the only lme4 function that exposes the REML criterion as a
named scalar. Without it, programmatic guards like `if (isREML(m))
REMLcrit(m)` fail. Also required by `lmerTest` internals.

---

## MS-08 — `isREML` generic is defined

**Name:** `isREML-method-exists`
**Kind:** unit / parity-vs-lme4
**Priority:** P1
**Gap classification:** in-scope-missing (major)
**Source finding:** F3

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML and ML fits.

### Assertion

```r
test_that("isREML() returns TRUE for REML fit and FALSE for ML fit", {
  expect_true(isREML(sleep_reml))
  expect_false(isREML(sleep_ml))

  # Return value is scalar logical
  expect_identical(length(isREML(sleep_reml)), 1L)
  expect_identical(storage.mode(isREML(sleep_reml)), "logical")
})

test_that("isREML() agrees with internal REML flag", {
  expect_identical(isREML(sleep_reml), isTRUE(sleep_reml$REML))
  expect_identical(isREML(sleep_ml),  isTRUE(sleep_ml$REML))
})
```

### What it guards against

`if (isREML(m))` is the standard lme4 idiom for criterion guards. Without
the generic, code in `lmerTest`, `emmeans`, and user scripts that calls
`isREML()` throws "no applicable method".

---

## MS-09 — `refitML` converts a REML fit to an ML fit

**Name:** `refitML-REML-to-ML`
**Kind:** unit / parity-vs-lme4 / numerical-tolerance
**Priority:** P1
**Gap classification:** in-scope-missing (major)
**Source finding:** F3

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit as input.

### Assertion

```r
test_that("refitML() converts a REML fit to an ML fit", {
  mm_skip_if_no_lme4()
  ml_from_reml <- refitML(sleep_reml)

  # Result is an mm_lmm object
  expect_s3_class(ml_from_reml, "mm_lmm")

  # isREML must be FALSE on the result
  expect_false(isREML(ml_from_reml))

  # logLik must match a direct ML fit to tolerance
  expect_equal(as.numeric(logLik(ml_from_reml)),
               as.numeric(logLik(sleep_ml)),
               tolerance = 1e-3,
               label = "refitML() logLik matches direct ML fit")

  # fixef must match to tolerance
  expect_equal(unname(fixef(ml_from_reml)), unname(fixef(sleep_ml)),
               tolerance = 1e-4,
               label = "refitML() fixef matches direct ML fit")
})

test_that("refitML() on an ML fit is a no-op or returns a valid ML fit", {
  result <- refitML(sleep_ml)
  expect_false(isREML(result))
  expect_equal(as.numeric(logLik(result)), as.numeric(logLik(sleep_ml)),
               tolerance = 1e-6)
})
```

### What it guards against

The standard LRT workflow calls `refitML()` before `anova(m1, m2)`. Without
it, users converting REML fits for LRT must fit twice from scratch.
`refit(m, newresp)` (response-refit) already exists but does not address the
criterion-refit use case.

---

## MS-10 — `getME(m, "devcomp")` is available (PRD §6 commitment)

**Name:** `getME-devcomp`
**Kind:** unit / integration
**Priority:** P1
**Gap classification:** in-scope-missing (major)
**Source finding:** F4 / PRD §6 lines 279-280

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.

### Assertion

```r
test_that("getME(m, 'devcomp') returns a named list with cmp and dims sub-lists", {
  dc <- getME(sleep_reml, "devcomp")

  # Must be a list (or named numeric vector) — not an error
  expect_true(is.list(dc) || is.numeric(dc),
              label = "getME('devcomp') must not throw")

  # Core fields the cmp sub-list carries in lme4
  if (is.list(dc) && !is.null(dc$cmp)) {
    expect_true(all(c("ldL2", "ldRX2", "wrss", "ussq", "pwrss",
                      "drsum", "REML", "n_snb", "ML") %in% names(dc$cmp)),
                label = "devcomp$cmp has expected lme4 fields")
  }

  # At minimum: if it is a list it must have 'cmp' and 'dims'
  if (is.list(dc)) {
    expect_true("cmp"  %in% names(dc), label = "devcomp has 'cmp' field")
    expect_true("dims" %in% names(dc), label = "devcomp has 'dims' field")
  }
})

test_that("getME(m, 'devcomp')$cmp['REML'] is consistent with isREML()", {
  dc <- getME(sleep_reml, "devcomp")
  if (is.list(dc) && !is.null(dc$cmp)) {
    expect_equal(dc$cmp[["REML"]], 1L,
                 label = "devcomp$cmp['REML']==1 for REML fit")
  }
  dc_ml <- getME(sleep_ml, "devcomp")
  if (is.list(dc_ml) && !is.null(dc_ml$cmp)) {
    expect_equal(dc_ml$cmp[["REML"]], 0L,
                 label = "devcomp$cmp['REML']==0 for ML fit")
  }
})
```

### What it guards against

PRD §6 explicitly commits to `getME(m, "devcomp")`. The current default branch
aborts with a structured "not available" error instead of returning a value.
Downstream packages that decompose the deviance components (e.g. checking
weighted residual sum of squares) rely on this name.

### Upstream fixture / engine change needed

The Rust artifact must expose the deviance components (`ldL2`, `ldRX2`,
`wrss`, `REML` flag, `pwrss`, `n_snb`) in a named sub-artifact. The wrapper
`getME.mm_lmm` must read and return them under the `"devcomp"` key.

---

## MS-11 — `getME(m, "Zt")` and `getME(m, "Lambdat")` do not crash

**Name:** `getME-Zt-Lambdat-dispatch-bug`
**Kind:** unit / parity-vs-lme4
**Priority:** P0 (blocker)
**Gap classification:** partial (major) — declared-but-broken
**Source finding:** F4 / `assessment/parity/inf-getME.md`

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.

### Assertion

```r
test_that("getME(m, 'Zt') returns a q x N sparse matrix", {
  Zt <- getME(sleep_reml, "Zt")

  # Must be a Matrix-package sparse object
  expect_true(inherits(Zt, c("sparseMatrix", "Matrix")),
              label = "getME('Zt') returns a sparseMatrix")

  # Dimensions: q x N (36 x 180 for sleepstudy correlated-slope model)
  expect_identical(dim(Zt), c(36L, 180L),
                   label = "getME('Zt') dimensions: q x N")

  # Must equal t(getME(m, 'Z'))
  Z <- getME(sleep_reml, "Z")
  expect_equal(as.matrix(Zt), as.matrix(Matrix::t(Z)), tolerance = 1e-10,
               label = "Zt == t(Z)")
})

test_that("getME(m, 'Lambdat') returns the q x q upper-triangular factor", {
  Lambdat <- getME(sleep_reml, "Lambdat")

  expect_true(inherits(Lambdat, c("sparseMatrix", "Matrix")),
              label = "getME('Lambdat') returns a sparseMatrix")

  expect_identical(dim(Lambdat), c(36L, 36L),
                   label = "getME('Lambdat') dimensions: q x q")

  # Must equal t(Lambda)
  Lambda <- getME(sleep_reml, "Lambda")
  expect_equal(as.matrix(Lambdat), as.matrix(Matrix::t(Lambda)),
               tolerance = 1e-10,
               label = "Lambdat == t(Lambda)")
})
```

### What it guards against

Both `Zt` and `Lambdat` are declared supported names in `getME.mm_lmm` but
crash with `t.default` dispatch because `t()` is called on a `Matrix` sparse
object without the `Matrix::` qualifier. The fix is a two-character change per
line in `R/revive.R` lines 128 and 130. Without the fix any caller that relies
on `Zt` or `Lambdat` (merTools, MixedPower, custom code) gets an error on a
documented name.

---

## MS-12 — `getME(m, "N")`, `"n"`, `"p"`, `"q"` dimension scalars

**Name:** `getME-dimension-scalars`
**Kind:** unit / parity-vs-lme4
**Priority:** P1
**Gap classification:** in-scope-missing (major)
**Source finding:** F4

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.
Expected values: `N=180`, `n=18`, `p=2`, `q=36`.

### Assertion

```r
test_that("getME(m, 'N'/'n'/'p'/'q') return correct dimension scalars", {
  # N: total observations
  expect_equal(getME(sleep_reml, "N"), 180L,
               label = "getME('N') == nobs")

  # n: number of grouping-factor levels (subjects)
  expect_equal(getME(sleep_reml, "n"), 18L,
               label = "getME('n') == n_subjects")

  # p: number of fixed-effect parameters
  expect_equal(getME(sleep_reml, "p"), 2L,
               label = "getME('p') == length(fixef)")

  # q: total random-effect dimensions (2 RE terms × 18 subjects = 36)
  expect_equal(getME(sleep_reml, "q"), 36L,
               label = "getME('q') == ncol(Z)")
})

test_that("getME dimension scalars are consistent with other extractors", {
  expect_equal(getME(sleep_reml, "N"), nobs(sleep_reml))
  expect_equal(getME(sleep_reml, "p"), length(fixef(sleep_reml)))
  expect_equal(getME(sleep_reml, "q"), ncol(getME(sleep_reml, "Z")))
})
```

### What it guards against

`N`, `n`, `p`, `q` are the most commonly used dimension queries in downstream
packages (merTools, MixedPower, glmmTMB compatibility layers). All four
currently throw "component X is not available". They are trivially derivable
from stored fields.

---

## MS-13 — `getME(m, "Gp")` is available (PRD §6 commitment)

**Name:** `getME-Gp`
**Kind:** unit
**Priority:** P1
**Gap classification:** in-scope-missing (major)
**Source finding:** F4 / PRD §6

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.
Expected: integer vector `c(0L, 36L)` (one random-effect group spanning
columns 1–36 of Z).

### Assertion

```r
test_that("getME(m, 'Gp') returns the group-pointer integer vector", {
  Gp <- getME(sleep_reml, "Gp")

  expect_true(is.integer(Gp) || is.numeric(Gp),
              label = "getME('Gp') returns numeric/integer")
  # Length = number of random-effect groups + 1
  expect_equal(length(Gp), 2L,
               label = "getME('Gp') length == n_RE_groups + 1 (sleepstudy has 1)")
  # First element always 0, last element == q
  expect_equal(Gp[[1L]], 0L,
               label = "Gp[1] == 0")
  expect_equal(Gp[[length(Gp)]], getME(sleep_reml, "q"),
               label = "Gp[last] == q")
})
```

### What it guards against

`Gp` is the standard indexer for slicing Z-columns by random-effect group.
PRD §6 commits to it. Without it, tools that use `getME(m, "Gp")` to extract
per-group blocks of `Lambda` or `Z` fall through to the "not available" error.

---

## MS-14 — `getME(m, "lower")` is available (PRD §6 commitment)

**Name:** `getME-lower`
**Kind:** unit
**Priority:** P1
**Gap classification:** in-scope-missing (major)
**Source finding:** F4 / PRD §6

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.
Expected: a numeric vector of lower bounds for each theta element (0 for
diagonal elements, `-Inf` for off-diagonal Cholesky factors).

### Assertion

```r
test_that("getME(m, 'lower') returns theta lower-bound vector", {
  lower <- getME(sleep_reml, "lower")

  expect_true(is.numeric(lower),
              label = "getME('lower') returns numeric")
  expect_equal(length(lower), length(getME(sleep_reml, "theta")),
               label = "getME('lower') length == length(theta)")

  # For a correlated random-slope model (full covariance):
  # theta = (sigma_int, rho, sigma_slope) with lower = c(0, -Inf, 0)
  # Diagonal entries >= 0, off-diagonal entries are unconstrained
  diag_positions <- lower >= 0
  expect_true(any(diag_positions),
              label = "At least one theta has lower bound >= 0 (diagonal SD)")
})
```

### What it guards against

`lower` is used by profiling code, boundary detection, and optimizer restart
logic. PRD §6 commits to it. Without it, any tool that checks whether theta
is on a boundary via the lower bounds gets an error.

---

## MS-15 — `getME(m, "optinfo")` is available (PRD §6 commitment)

**Name:** `getME-optinfo`
**Kind:** unit / integration
**Priority:** P2
**Gap classification:** in-scope-missing (major)
**Source finding:** F4 / PRD §6

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.

### Assertion

```r
test_that("getME(m, 'optinfo') returns optimizer metadata", {
  optinfo <- getME(sleep_reml, "optinfo")

  # Must not throw; must be a list or named vector
  expect_true(is.list(optinfo) || is.character(optinfo),
              label = "getME('optinfo') must not throw")

  # If a list, must carry at least the optimizer name or status
  if (is.list(optinfo)) {
    has_expected <- any(c("optimizer", "conv", "val", "status") %in% names(optinfo))
    expect_true(has_expected,
                label = "optinfo list carries at least one of: optimizer/conv/val/status")
  }
})
```

### What it guards against

PRD §6 commits to `optinfo`. It is used by convergence-checking utilities
(e.g., `allFit` comparisons). Without it, those utilities error on a
documented name.

---

## MS-16 — `getME(m, "sigma")` delegates to `sigma()`

**Name:** `getME-sigma-alias`
**Kind:** unit
**Priority:** P2
**Gap classification:** partial (minor)
**Source finding:** F4 / `assessment/parity/inf-getME.md`

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML fit.

### Assertion

```r
test_that("getME(m, 'sigma') returns the residual SD", {
  s <- getME(sleep_reml, "sigma")
  expect_equal(s, sigma(sleep_reml), tolerance = 1e-10,
               label = "getME('sigma') == sigma()")
  expect_true(is.numeric(s) && length(s) == 1L && s > 0,
              label = "getME('sigma') is a positive scalar")
})
```

### What it guards against

lme4 supports `getME(m, "sigma")` as an alias for `sigma()`. mixeff currently
errors. Downstream packages that route through `getME` for all scalar
parameters (e.g. model-comparison utilities) will hit this error.

---

## MS-17 — `getME(m, "REML")` / `"is_REML"` return criterion flag

**Name:** `getME-REML-flag`
**Kind:** unit
**Priority:** P2
**Gap classification:** partial (minor)
**Source finding:** F4

### Dataset / formula

`sleepstudy`, `Reaction ~ Days + (Days|Subject)`, REML and ML fits.

### Assertion

```r
test_that("getME(m, 'REML') returns 2L for REML fit, 0L for ML fit", {
  # lme4 convention: REML deviance offset; 2L=REML, 0L=ML
  reml_flag <- getME(sleep_reml, "REML")
  ml_flag   <- getME(sleep_ml,   "REML")

  expect_true(reml_flag != 0L,
              label = "getME('REML') non-zero for REML fit")
  expect_equal(ml_flag, 0L,
               label = "getME('REML') == 0 for ML fit")
})

test_that("getME(m, 'is_REML') returns TRUE/FALSE logical", {
  expect_true(getME(sleep_reml, "is_REML"))
  expect_false(getME(sleep_ml,  "is_REML"))
})
```

### What it guards against

Programmatic guards in downstream packages (e.g. `if (getME(m, "is_REML"))`)
currently error in mixeff. The state exists as `m$REML` but is inaccessible
via the lme4 name.

---

## MS-18 — `AIC(m1, m2)` / `BIC(m1, m2)` multi-model route gives actionable error

**Name:** `AIC-BIC-multi-model-refusal-message`
**Kind:** error-message
**Priority:** P2
**Gap classification:** partial (minor) — deliberate routing
**Source finding:** summary table rows "AIC multi-model" / "BIC multi-model"

### Dataset / formula

Two ML fits on `sleepstudy`: one with `Days` slope, one random-intercept only.

### Assertion

```r
test_that("AIC(m1, m2) gives a typed mm_inference_unavailable with compare() hint", {
  m1 <- lmm(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE,
             control = mm_control(verbose = -1))
  m2 <- lmm(Reaction ~ Days + (1|Subject), sleepstudy, REML = FALSE,
             control = mm_control(verbose = -1))

  err <- expect_error(AIC(m1, m2), class = "mm_inference_unavailable")
  expect_match(conditionMessage(err), "compare()", fixed = TRUE,
               label = "AIC multi-model error mentions compare()")
})

test_that("BIC(m1, m2) gives a typed mm_inference_unavailable with compare() hint", {
  m1 <- lmm(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE,
             control = mm_control(verbose = -1))
  m2 <- lmm(Reaction ~ Days + (1|Subject), sleepstudy, REML = FALSE,
             control = mm_control(verbose = -1))

  err <- expect_error(BIC(m1, m2), class = "mm_inference_unavailable")
  expect_match(conditionMessage(err), "compare()", fixed = TRUE,
               label = "BIC multi-model error mentions compare()")
})
```

### What it guards against

The multi-model refusal is intentional (users are routed to `compare()`), but
this test confirms the error is (a) a typed condition and (b) contains the
actionable routing hint. A plain `stop()` without a class or the hint is a
regression.

---

## Summary table

| Spec ID | Name | Kind | Priority | Classification |
|---------|------|------|----------|----------------|
| MS-01 | `logLik-nall-attribute` | parity-vs-lme4 | P2 | partial |
| MS-02 | `logLik-REML-override-REML-fit` | parity-vs-lme4 | **P0** | partial / test-gap |
| MS-03 | `logLik-REML-override-ML-fit` | parity-vs-lme4 | P1 | partial / test-gap |
| MS-04 | `deviance-REML-fit-returns-ML` | parity-vs-lme4 | P1 | partial / test-gap |
| MS-05 | `deviance-REML-arg-not-ignored` | parity-vs-lme4 | P1 | partial / test-gap |
| MS-06 | `extractAIC-method-exists` | unit / integration | **P0** | in-scope-missing |
| MS-07 | `REMLcrit-method-exists` | unit / parity-vs-lme4 | P1 | in-scope-missing |
| MS-08 | `isREML-method-exists` | unit / parity-vs-lme4 | P1 | in-scope-missing |
| MS-09 | `refitML-REML-to-ML` | unit / parity-vs-lme4 | P1 | in-scope-missing |
| MS-10 | `getME-devcomp` | unit / integration | P1 | in-scope-missing (PRD §6) |
| MS-11 | `getME-Zt-Lambdat-dispatch-bug` | unit / parity-vs-lme4 | **P0** | partial (broken) |
| MS-12 | `getME-dimension-scalars` | unit / parity-vs-lme4 | P1 | in-scope-missing |
| MS-13 | `getME-Gp` | unit | P1 | in-scope-missing (PRD §6) |
| MS-14 | `getME-lower` | unit | P1 | in-scope-missing (PRD §6) |
| MS-15 | `getME-optinfo` | unit / integration | P2 | in-scope-missing (PRD §6) |
| MS-16 | `getME-sigma-alias` | unit | P2 | partial (minor) |
| MS-17 | `getME-REML-flag` | unit | P2 | partial (minor) |
| MS-18 | `AIC-BIC-multi-model-refusal-message` | error-message | P2 | partial (minor) |

**Total specs: 18**

**Top priority:** MS-11 (`getME-Zt-Lambdat-dispatch-bug`) — two declared
supported names crash with a method-dispatch error that has a trivial two-line
fix (`t()` → `Matrix::t()` in `R/revive.R` lines 128 and 130); also
co-P0 with MS-02 (`logLik-REML-override-REML-fit`) which silently returns a
wrong number — the worst possible failure mode — and MS-06
(`extractAIC-method-exists`) which blocks base-R `step()`/`drop1()`.

---

## Upstream fixtures / engine changes required

| Spec | Requirement |
|------|-------------|
| MS-02, MS-03 | Rust engine must expose a re-evaluate path: given stored θ, compute the objective under the *opposite* criterion (ML given REML θ, or REML given ML θ). If not implemented, wrapper must refuse with a typed error rather than return silently wrong output. |
| MS-04, MS-05 | Rust artifact must store *both* the REML criterion and the ML deviance in distinct named fields so the R wrapper can separate them. Currently `m$deviance` carries the fit-criterion value. |
| MS-09 (`refitML`) | Wrapper needs a re-fit call with `REML=FALSE`; the existing `refit(m, newresp)` path does not cover this. |
| MS-10 (`devcomp`) | Rust artifact must expose deviance components (`ldL2`, `ldRX2`, `wrss`, `REML` flag, `pwrss`, `n_snb`) in a named sub-artifact. |
| MS-13 (`Gp`), MS-14 (`lower`) | Rust artifact or R-side algebra must produce `Gp` (group-pointer vector) and `lower` (theta lower-bounds vector); both are derivable from the stored random-effect structure without a live Rust handle. |

All other specs (MS-01, MS-06–MS-08, MS-11–MS-12, MS-15–MS-18) can be
implemented in pure R without upstream changes.
