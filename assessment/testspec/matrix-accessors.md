# Test Specification: Matrix & Structure Accessors

**Family:** Matrix & structure accessors (`getME`, `model.matrix`, `model.frame`,
`terms`, `vcov`, `weights`, `ngrps`, `hatvalues`, …)
**Source gap report:** `assessment/gap/matrix-accessors.md`
**Written:** 2026-05-31
**Standard dataset/formula (default unless noted):**
```r
library(mixeff); library(lme4)
data(sleepstudy, package = "lme4")
mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
          control = mm_control(verbose = -1))
fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
```

Tolerance conventions follow `planning/PRD.md` §11:
fixef 1e-4 · theta 1e-3 · logLik 1e-3 · sigma 1e-4 · matrix elements 1e-4.

---

## SPEC-MA-01 — `getME(.,"Zt")` returns the correct transposed sparse matrix

**Priority:** P0 blocker
**Kind:** unit + parity-vs-lme4
**Gap classification:** partial (major) — advertised name errors with cryptic
`"argument is not a matrix"` from `base::t` called on a `dgCMatrix`.
**Root cause identified in:** `R/revive.R` lines 128 and 130 (bare `t()` not
dispatching Matrix S4 method).

### Assertion

```r
test_that("getME(fit, 'Zt') returns transposed Z without error", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("Matrix")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  # Must not throw
  Zt_mf <- getME(mf, "Zt")

  # Class must be a sparse Matrix (dgCMatrix or compatible)
  expect_true(inherits(Zt_mf, "Matrix"))

  # Dimensions: Zt is q×N (transpose of N×q Z)
  Z_mf <- getME(mf, "Z")
  expect_equal(dim(Zt_mf), rev(dim(Z_mf)))

  # Numeric values match lme4's Zt within tolerance 1e-4
  Zt_fm <- lme4::getME(fm, "Zt")
  expect_equal(as.matrix(Zt_mf), as.matrix(Zt_fm), tolerance = 1e-4,
               ignore_attr = TRUE)
})
```

**What it guards against:** regression of the `base::t` dispatch bug; any future
change that re-introduces the `"argument is not a matrix"` error on a sparse
`dgCMatrix`.

**Upstream fixture / engine change required:** None — this is a pure R-side fix
(`Matrix::t()` instead of `t()`).

---

## SPEC-MA-02 — `getME(.,"Lambdat")` returns the correct transposed relative-covariance factor

**Priority:** P0 blocker
**Kind:** unit + parity-vs-lme4
**Gap classification:** partial (major) — same `base::t` dispatch bug as Zt,
triggered on the lazy-computed `Lambda` block-diagonal matrix.

### Assertion

```r
test_that("getME(fit, 'Lambdat') returns transposed Lambda without error", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("Matrix")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  # Must not throw
  Lambdat_mf <- getME(mf, "Lambdat")

  # Class must be a sparse Matrix
  expect_true(inherits(Lambdat_mf, "Matrix"))

  # Dimensions: Lambdat is q×q (transpose of Lambda q×q)
  Lambda_mf <- getME(mf, "Lambda")
  expect_equal(dim(Lambdat_mf), rev(dim(Lambda_mf)))

  # Values match lme4 within tolerance 1e-3 (theta tolerance)
  Lambdat_fm <- lme4::getME(fm, "Lambdat")
  expect_equal(as.matrix(Lambdat_mf), as.matrix(Lambdat_fm), tolerance = 1e-3,
               ignore_attr = TRUE)
})
```

**What it guards against:** same `base::t` dispatch regression as SPEC-MA-01 on
the Lambda branch.

**Upstream fixture / engine change required:** None.

---

## SPEC-MA-03 — `getME(.,"sigma")` is a documented name that must not error

**Priority:** P1
**Kind:** unit + error-message
**Gap classification:** partial (minor) — `getME(mf,"sigma")` throws "not
available" error even though `sigma(mf)` works and lme4 documents `"sigma"` as a
valid `getME` name.

### Assertion

```r
test_that("getME(fit, 'sigma') returns the residual SD scalar", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  # Must not throw
  s <- getME(mf, "sigma")

  expect_true(is.numeric(s) && length(s) == 1L && is.finite(s))
  # Must agree with the sigma() generic
  expect_equal(s, sigma(mf), tolerance = 1e-10)
})
```

**What it guards against:** regression where `getME(.,"sigma")` is wired to the
unavailable fallback despite the value being readily available from `fit$sigma`.

**Upstream fixture / engine change required:** None.

---

## SPEC-MA-04 — `ngrps()` returns a named integer vector matching lme4

**Priority:** P0 blocker
**Kind:** parity-vs-lme4
**Gap classification:** in-scope-missing (major) — `ngrps(mf)` errors; lme4
returns `c(Subject = 18)`.

### Assertion

```r
test_that("ngrps() returns a named integer matching lme4", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  ng_mf <- ngrps(mf)
  ng_fm <- lme4::ngrps(fm)

  # Integer type, named
  expect_true(is.integer(ng_mf))
  expect_true(!is.null(names(ng_mf)))

  # Same group names and counts
  expect_setequal(names(ng_mf), names(ng_fm))
  for (grp in names(ng_fm)) {
    expect_equal(ng_mf[[grp]], ng_fm[[grp]],
                 info = sprintf("ngrps mismatch for group '%s'", grp))
  }
})

test_that("ngrps() works for a crossed two-grouping model", {
  skip_if_not_installed("lme4")
  data(Pastes, package = "lme4")
  mf2 <- lmm(strength ~ (1 | batch) + (1 | sample), data = Pastes,
             control = mm_control(verbose = -1))
  fm2 <- lme4::lmer(strength ~ (1 | batch) + (1 | sample), data = Pastes)

  ng_mf2 <- ngrps(mf2)
  ng_fm2 <- lme4::ngrps(fm2)
  expect_setequal(names(ng_mf2), names(ng_fm2))
  for (grp in names(ng_fm2)) {
    expect_equal(ng_mf2[[grp]], ng_fm2[[grp]],
                 info = sprintf("ngrps mismatch for group '%s' in Pastes", grp))
  }
})
```

**What it guards against:** `ngrps()` being absent or returning wrong counts,
especially for crossed or nested grouping structures.

**Upstream fixture / engine change required:** None — group level counts are
derivable from `getME(fit,"flist")` which already works.

---

## SPEC-MA-05 — `weights()` returns a numeric vector of length N (not NULL) for unweighted fits

**Priority:** P0 blocker
**Kind:** parity-vs-lme4
**Gap classification:** partial (major) — `weights(mf)` returns `NULL`; lme4
returns an all-ones numeric vector of length N.  NULL breaks `weighted.mean`
and broom-family code.

### Assertion

```r
test_that("weights() returns length-N all-ones vector for an unweighted LMM", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  w_mf <- weights(mf)
  w_fm <- weights(fm)   # lme4 returns all-ones numeric(180)

  expect_true(!is.null(w_mf),
              info = "weights() must not return NULL for an unweighted fit")
  expect_true(is.numeric(w_mf))
  expect_equal(length(w_mf), nobs(mf))
  expect_equal(w_mf, w_fm, tolerance = 1e-10,
               info = "weights() should match lme4 all-ones vector")
})

test_that("weighted.mean of fitted values works after weights()", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  w <- weights(mf)
  # Must not error; weighted.mean(x, NULL) errors in base R
  expect_no_error(weighted.mean(fitted(mf), w))
})
```

**What it guards against:** downstream code (`broom.mixed::augment`, manual
`weighted.mean`) crashing on NULL weights returned from an unweighted fit.

**Upstream fixture / engine change required:** None.

---

## SPEC-MA-06 — `terms()` method exists and returns a `terms` object

**Priority:** P1
**Kind:** unit + integration
**Gap classification:** in-scope-missing (major) — `terms(mf)` errors with
`"no terms component nor attribute"`; no `terms.mm_lmm` method.

### Assertion

```r
test_that("terms(fit) returns a terms object without error", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  tt <- terms(mf)

  expect_true(inherits(tt, "terms"))
  # Fixed-effect variables must appear
  vars <- attr(tt, "term.labels")
  expect_true("Days" %in% vars)
})

test_that("terms(fit, fixed.only = TRUE) returns only fixed terms", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  tt_fixed <- terms(mf, fixed.only = TRUE)
  expect_true(inherits(tt_fixed, "terms"))
  # Random-effect bar notation must not appear in fixed-only terms
  labs <- attr(tt_fixed, "term.labels")
  expect_false(any(grepl("|", labs, fixed = TRUE)))
})
```

**What it guards against:** `terms.mm_lmm` being absent, which would break any
code (e.g. `model.matrix` on a new `data.frame`) that calls `terms()` on the
fitted model.

**Note:** emmeans integration already works via dedicated `recover_data` /
`emm_basis`; this spec guards the standalone `terms()` accessor path.

**Upstream fixture / engine change required:** None.

---

## SPEC-MA-07 — `hatvalues()` returns a length-N numeric vector for LMMs

**Priority:** P1
**Kind:** unit + numerical-tolerance (lme4 parity)
**Gap classification:** in-scope-missing (major) — no `hatvalues.mm_lmm` method
at all; `hatvalues(mf)` errors.

### Assertion

```r
test_that("hatvalues() returns a length-N numeric vector", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  hv <- hatvalues(mf)

  expect_true(is.numeric(hv))
  expect_equal(length(hv), nobs(mf))
  # All hat values in (0, 1]
  expect_true(all(hv > 0 & hv <= 1 + 1e-8),
              info = "hat values must be in (0,1]")
  # Sum of hat values = number of model parameters (trace of hat matrix)
  # For an LMM with p fixed parameters the sum is approximately p;
  # a loose check guards against degenerate all-zero returns.
  expect_true(sum(hv) > 1,
              info = "sum of hat values must exceed 1")
})

test_that("hatvalues() numeric values are close to lme4 within 1e-4", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  hv_mf <- hatvalues(mf)
  hv_fm <- hatvalues(fm)

  expect_equal(hv_mf, hv_fm, tolerance = 1e-4,
               ignore_attr = TRUE,
               info = "hatvalues parity failed vs lme4")
})
```

**What it guards against:** leverage diagnostics being completely unavailable;
downstream Cook's-distance and influence workflows failing silently or crashing.

**Note on implementation path:** `hatvalues` for an LMM equals the diagonal of
`H = Z Λ (Z Λ)' (Z Λ (Z Λ)' + σ² I)^{-1}` plus the contribution from X; a
viable shortcut is `diag(X %*% vcov(fit) %*% t(X)) / sigma(fit)^2` combined
with the RE hat matrix.  If the Rust side does not yet expose `L` or `RX`
directly, an R-side approximation using the already-available Z, Lambda, and
vcov matrices is acceptable at this tolerance.  Mark as
`attr(hv, "mm_method") = "r_side_approximation"` if so; the parity test will
catch value-level regressions.

**Upstream fixture / engine change required:** Potentially — if an exact
Cholesky-based implementation is desired, `getME(.,"L")` and/or `getME(.,"RX")`
need to be wired from the Rust FFI.  The parity test is written to tolerate an
R-side approximation so it can be satisfied without blocking on the FFI work.

---

## SPEC-MA-08 — `vcov(fit, correlation = TRUE)` either returns a correlation matrix or raises a typed diagnostic

**Priority:** P1
**Kind:** error-message + parity-vs-lme4
**Gap classification:** partial (major) — `vcov(mf, correlation = TRUE)` silently
ignores the argument and returns the ordinary covariance matrix, violating the
project's no-silent-behavior principle.

### Assertion

```r
test_that("vcov(fit, correlation = TRUE) either honors the flag or raises mm_arg_error", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  # Option A: correlation matrix is implemented
  # Option B: a typed diagnostic is raised instead of silent ignore
  result <- tryCatch(
    vcov(mf, correlation = TRUE),
    mm_arg_error = function(cnd) cnd,
    mm_inference_unavailable = function(cnd) cnd
  )

  if (inherits(result, "condition")) {
    # Typed refusal is acceptable — silent ignore is not
    expect_true(
      inherits(result, "mm_arg_error") || inherits(result, "mm_inference_unavailable"),
      info = "vcov(correlation=TRUE) must raise a typed mm_ condition, not fail silently"
    )
    expect_true(nzchar(conditionMessage(result)),
                info = "typed refusal must carry a non-empty message")
  } else {
    # If implemented, result must be a symmetric matrix with ones on the diagonal
    p <- length(fixef(mf))
    expect_true(is.matrix(result))
    expect_equal(dim(result), c(p, p))
    expect_equal(diag(result), rep(1, p), tolerance = 1e-10,
                 info = "correlation matrix must have 1s on the diagonal")
    expect_true(isSymmetric(result, tol = 1e-8))
  }
})

test_that("vcov(fit, correlation = TRUE) parity vs lme4 when implemented", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  result <- tryCatch(vcov(mf, correlation = TRUE), condition = function(cnd) NULL)
  if (is.null(result)) skip("vcov(correlation=TRUE) not yet implemented")

  expected <- stats::cov2cor(as.matrix(stats::vcov(fm)))
  expect_equal(result, expected, tolerance = 1e-4, ignore_attr = TRUE)
})
```

**What it guards against:** silent ignore of `correlation = TRUE` which gives the
user wrong output (a covariance matrix they believe is a correlation matrix)
without any diagnostic.

**Upstream fixture / engine change required:** None.

---

## SPEC-MA-09 — `getME(.,"u")` and `getME(.,"b")` expose spherical and conditional RE modes

**Priority:** P1
**Kind:** unit + numerical-tolerance
**Gap classification:** in-scope-missing (major) — both names return "not
available"; `ranef()` works (per-term data frames) but the stacked-vector forms
required for algebra (leverage, Cook's D, score equations) are absent.

### Assertion

```r
test_that("getME(fit, 'u') returns a numeric vector of length q (spherical modes)", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  u <- getME(mf, "u")

  expect_true(is.numeric(u))
  q <- ncol(getME(mf, "Z"))     # 36 for the slopes | Subject term
  expect_equal(length(u), q,
               info = "length(u) must equal ncol(Z)")
})

test_that("getME(fit, 'u') numeric values match lme4 within theta tolerance", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  u_mf <- getME(mf, "u")
  u_fm <- lme4::getME(fm, "u")
  expect_equal(u_mf, u_fm, tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("getME(fit, 'b') returns Lambda %*% u", {
  skip_if_not_installed("Matrix")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  u <- getME(mf, "u")
  b <- getME(mf, "b")
  Lambda <- getME(mf, "Lambda")

  expect_equal(length(b), length(u))
  # b = Lambda u, so ||b - Lambda u|| should be tiny
  expect_equal(as.numeric(b), as.numeric(Lambda %*% u), tolerance = 1e-8)
})
```

**What it guards against:** downstream leverage/diagnostics code that needs the
stacked-vector spherical modes; currently only per-group data frames are
available via `ranef()`.

**Upstream fixture / engine change required:** Uncertain — need to confirm whether
the Rust FFI serializes the spherical `u` vector.  If not, this is
upstream-blocked; if yes, it is unwired R-side only.  Track separately in the
mote bead once FFI contract is checked.

---

## SPEC-MA-10 — `getME(.,"Gp")` returns a group-pointer integer vector

**Priority:** P2
**Kind:** unit + parity-vs-lme4
**Gap classification:** in-scope-missing (major) — "not available"; needed to
slice Z/u/Lambda by RE term.

### Assertion

```r
test_that("getME(fit, 'Gp') returns an integer pointer vector of length n_rtrms + 1", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  Gp_mf <- getME(mf, "Gp")
  Gp_fm <- lme4::getME(fm, "Gp")

  expect_true(is.integer(Gp_mf))
  # Length must be number of RE terms + 1
  expect_equal(length(Gp_mf), length(Gp_fm))
  # Values must match lme4 exactly (cumulative column counts)
  expect_equal(Gp_mf, Gp_fm)
})
```

**What it guards against:** inability to slice block-structured matrices by RE
term; required by any user doing manual variance-component or leverage algebra.

**Upstream fixture / engine change required:** None — derivable R-side from
`flist` and per-term sizes already in the artifact.

---

## SPEC-MA-11 — `getME(.,"devcomp")` returns a list with `$cmp` and `$dims` sub-lists

**Priority:** P1
**Kind:** unit + parity-vs-lme4
**Gap classification:** in-scope-missing (major) — "not available"; used widely
by downstream tooling (e.g. `lme4::isSingular`, broom.mixed).

### Assertion

```r
test_that("getME(fit, 'devcomp') returns cmp and dims sub-lists", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  dc_mf <- getME(mf, "devcomp")
  dc_fm <- lme4::getME(fm, "devcomp")

  expect_true(is.list(dc_mf))
  expect_true(all(c("cmp", "dims") %in% names(dc_mf)))

  # Key deviance-component scalars
  cmp_mf <- dc_mf$cmp
  cmp_fm <- dc_fm$cmp
  shared <- intersect(names(cmp_mf), names(cmp_fm))
  expect_true(length(shared) >= 2L,
              info = "devcomp$cmp must share at least 2 named scalars with lme4")
  for (nm in intersect(c("sigmaML", "sigmaREML", "ldL2", "pwrss"), shared)) {
    expect_equal(cmp_mf[[nm]], cmp_fm[[nm]], tolerance = 1e-4,
                 info = sprintf("devcomp$cmp$%s mismatch vs lme4", nm))
  }

  # Dims must include basic integer dimensions
  dims_mf <- dc_mf$dims
  expect_true("n" %in% names(dims_mf) || "N" %in% names(dims_mf),
              info = "devcomp$dims must contain observation count")
})
```

**What it guards against:** downstream tooling that reads `getME(.,"devcomp")`
being silently broken.

**Upstream fixture / engine change required:** None — deviance components (logLik,
sigma, wrss) and dims (N, p, q) are all available R-side.

---

## SPEC-MA-12 — `getME(.,"lower")` returns the theta lower-bounds vector

**Priority:** P2
**Kind:** unit + parity-vs-lme4
**Gap classification:** in-scope-missing (major) — "not available"; lme4's
`isSingular()` reads this vector; mixeff has its own `is_singular()` path but
the raw `lower` vector is unavailable for users.

### Assertion

```r
test_that("getME(fit, 'lower') returns the theta lower-bounds vector", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  lower_mf <- getME(mf, "lower")
  lower_fm <- lme4::getME(fm, "lower")

  expect_true(is.numeric(lower_mf))
  expect_equal(length(lower_mf), length(getME(mf, "theta")))
  # Lower bounds: diagonal entries 0, off-diagonal entries -Inf
  expect_equal(lower_mf, lower_fm,
               info = "lower bounds vector must match lme4")
})
```

**What it guards against:** users manually implementing `isSingular`-style checks
(lower == theta → boundary) being blocked by a missing `lower` vector.

**Upstream fixture / engine change required:** None — lower bounds are determined
by covariance parameterization (diagonal elements ≥ 0, off-diagonal unrestricted)
which is known R-side from the artifact.

---

## SPEC-MA-13 — `getME(.,"N")`, `"n"`, `"p"`, `"q"` return correct dimension scalars

**Priority:** P2
**Kind:** unit + parity-vs-lme4
**Gap classification:** in-scope-missing (minor) — all return "not available";
dimension scalars are trivially derivable but currently not wired.

### Assertion

```r
test_that("getME() dimension scalars N, n, p, q match lme4", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  for (dim_name in c("n", "p", "q")) {
    val_mf <- getME(mf, dim_name)
    val_fm <- lme4::getME(fm, dim_name)
    expect_equal(val_mf, val_fm,
                 info = sprintf("getME(.,'%s') mismatch vs lme4", dim_name))
  }
  # "N" is a documented alias for "n" in lme4
  n_mf <- getME(mf, "n")
  N_mf <- getME(mf, "N")
  expect_equal(n_mf, N_mf,
               info = "'N' and 'n' must return the same observation count")
})
```

**What it guards against:** loops over lme4 recipes that rely on dimension
scalars crashing silently at "not available".

**Upstream fixture / engine change required:** None.

---

## SPEC-MA-14 — `getME(.,"n_rtrms")` and `"n_rfacs"` return RE term and factor counts

**Priority:** P2
**Kind:** unit + parity-vs-lme4
**Gap classification:** in-scope-missing (minor)

### Assertion

```r
test_that("getME(fit, 'n_rtrms') and 'n_rfacs' return correct counts", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  expect_equal(getME(mf, "n_rtrms"), lme4::getME(fm, "n_rtrms"))
  expect_equal(getME(mf, "n_rfacs"), lme4::getME(fm, "n_rfacs"))
})
```

**What it guards against:** term/factor counts being absent for downstream code
iterating over RE terms programmatically.

**Upstream fixture / engine change required:** None — derivable from `flist` and
artifact random term list.

---

## SPEC-MA-15 — `getME(.,"REML")` / `isREML()` return the REML flag

**Priority:** P2
**Kind:** unit + parity-vs-lme4
**Gap classification:** in-scope-missing (minor)

### Assertion

```r
test_that("getME(fit, 'REML') returns the REML logical flag", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf_reml <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                 REML = TRUE, control = mm_control(verbose = -1))
  mf_ml   <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                 REML = FALSE, control = mm_control(verbose = -1))
  fm_reml <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = TRUE)
  fm_ml   <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)

  expect_equal(getME(mf_reml, "REML"), lme4::getME(fm_reml, "REML"))
  expect_equal(getME(mf_ml,   "REML"), lme4::getME(fm_ml,   "REML"))
  expect_true( getME(mf_reml, "REML"))
  expect_false(getME(mf_ml,   "REML"))
})
```

**What it guards against:** code that checks `getME(.,"REML")` to switch between
ML and REML inference paths receiving "not available" instead of a logical.

**Upstream fixture / engine change required:** None — `fit$REML` is stored.

---

## SPEC-MA-16 — `model.matrix(fit, type = "randomListRaw")` returns per-term matrix list

**Priority:** P2
**Kind:** unit
**Gap classification:** in-scope-missing (minor) — `match.arg` rejects the value;
lme4 uses this type for per-term raw matrices.

### Assertion

```r
test_that("model.matrix(fit, type='randomListRaw') returns a named list of matrices", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  raw_list <- model.matrix(mf, type = "randomListRaw")

  expect_true(is.list(raw_list))
  # One element per RE term
  expect_equal(length(raw_list), 1L,  # one `(Days | Subject)` term
               info = "randomListRaw must have one element per RE term")
  # Each element is a matrix with nrow = N
  for (mat in raw_list) {
    expect_true(is.matrix(mat) || inherits(mat, "Matrix"))
    expect_equal(nrow(mat), nobs(mf))
  }
})
```

**What it guards against:** per-term matrix extraction failing; required for
advanced users constructing custom Ztlist-style objects.

**Upstream fixture / engine change required:** None.

---

## SPEC-MA-17 — `getME(.,"Ztlist")` / `"mmList"` return per-term Z' or design-matrix list

**Priority:** P2
**Kind:** unit + parity-vs-lme4
**Gap classification:** in-scope-missing (major) — "not available"; blocks
`model.matrix(type="randomListRaw")` and per-term algebra.

### Assertion

```r
test_that("getME(fit, 'Ztlist') returns a named list of sparse per-term Z' matrices", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  Ztlist_mf <- getME(mf, "Ztlist")
  Ztlist_fm <- lme4::getME(fm, "Ztlist")

  expect_true(is.list(Ztlist_mf))
  expect_equal(length(Ztlist_mf), length(Ztlist_fm),
               info = "Ztlist length must equal number of RE terms")
  for (i in seq_along(Ztlist_fm)) {
    expect_equal(dim(Ztlist_mf[[i]]), dim(Ztlist_fm[[i]]),
                 info = sprintf("Ztlist element %d dimension mismatch", i))
    expect_equal(as.matrix(Ztlist_mf[[i]]), as.matrix(Ztlist_fm[[i]]),
                 tolerance = 1e-4, ignore_attr = TRUE,
                 info = sprintf("Ztlist element %d values mismatch", i))
  }
})
```

**What it guards against:** per-term RE design matrices being inaccessible to
users doing penalty-based algebra or variance-component slice operations.

**Upstream fixture / engine change required:** None — Z is already assembled
R-side from term-level pieces; the list form is just an unaggregated version.

---

## SPEC-MA-18 — Unknown `getME` names produce a typed `mm_arg_error` (regression guard)

**Priority:** P1
**Kind:** error-message
**Gap classification:** partial (minor) — this currently works correctly for
ordinary unknown names, but `Zt`/`Lambdat` throw a cryptic `base::t` error
instead of a typed condition.  This spec guards the correct-path behavior and
prevents regression once SPEC-MA-01/02 are fixed.

### Assertion

```r
test_that("getME() with an unknown name raises a typed mm_arg_error", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  expect_error(getME(mf, "not_a_real_name"), class = "mm_arg_error")
  expect_error(getME(mf, "devfun"),          class = "mm_arg_error")
  expect_error(getME(mf, "par"),             class = "mm_arg_error")
})

test_that("getME(fit, 'Zt') does NOT throw base-R matrix error (regression for MA-01 fix)", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  # After SPEC-MA-01 fix, this must not error with "argument is not a matrix"
  err <- tryCatch(getME(mf, "Zt"), error = function(e) e)
  if (inherits(err, "error")) {
    expect_false(
      grepl("argument is not a matrix", conditionMessage(err), fixed = TRUE),
      info = "getME(.,'Zt') must not produce the base-R 'argument is not a matrix' error"
    )
  }
})
```

**What it guards against:** regression to the cryptic `base::t` error path after
the Matrix S4 dispatch fix; and future contributors inadvertently breaking the
clean-error path for genuinely unsupported names.

**Upstream fixture / engine change required:** None.

---

## SPEC-MA-19 — `getME(.,"ALL")` returns a named list of all available components

**Priority:** P2
**Kind:** unit + parity-vs-lme4
**Gap classification:** in-scope-missing (major) — "not available"; a real lme4
user reflexively calls `getME(fit,"ALL")` for inspection.

### Assertion

```r
test_that("getME(fit, 'ALL') returns a named list with at least the core available components", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  all_comps <- getME(mf, "ALL")

  expect_true(is.list(all_comps))
  expect_true(!is.null(names(all_comps)))

  # At minimum, the currently-working components must appear
  always_present <- c("X", "Z", "Zt", "Lambda", "Lambdat",
                       "theta", "beta", "fixef", "y", "mu", "flist", "cnms")
  for (nm in always_present) {
    expect_true(nm %in% names(all_comps),
                info = sprintf("getME(.,'ALL') must include component '%s'", nm))
    expect_false(is.null(all_comps[[nm]]),
                 info = sprintf("getME(.,'ALL')$%s must not be NULL", nm))
  }
})
```

**What it guards against:** `getME(.,"ALL")` being entirely broken for users
doing exploratory model inspection.

**Upstream fixture / engine change required:** None — `"ALL"` is just a loop
over all registered available names.

---

## SPEC-MA-20 — Cross-session revival: matrix accessors survive `saveRDS` / `readRDS`

**Priority:** P1
**Kind:** cross-session-revival
**Gap classification:** test-gap — the gap report confirms these work in a live
session but the revival path (`revive()`, `lazy_cache` reconstruction) has no
dedicated test for matrix accessors.

### Assertion

```r
test_that("getME() matrix accessors work after saveRDS/readRDS round-trip", {
  data(sleepstudy, package = "lme4")
  mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            control = mm_control(verbose = -1))

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(mf, tmp)
  mf2 <- readRDS(tmp)

  # Lazy cache must rebuild
  expect_no_error(X2   <- getME(mf2, "X"))
  expect_no_error(Z2   <- getME(mf2, "Z"))
  expect_no_error(Zt2  <- getME(mf2, "Zt"))
  expect_no_error(L2   <- getME(mf2, "Lambda"))
  expect_no_error(Lt2  <- getME(mf2, "Lambdat"))

  # Values must agree with live session
  expect_equal(as.matrix(X2),  as.matrix(getME(mf, "X")),  tolerance = 1e-10)
  expect_equal(as.matrix(Z2),  as.matrix(getME(mf, "Z")),  tolerance = 1e-10)
  expect_equal(as.matrix(Zt2), as.matrix(getME(mf, "Zt")), tolerance = 1e-10)
})
```

**What it guards against:** lazy-cache reconstruction silently failing after
deserialization, especially after the `base::t` fix is applied (where the sparse
method must dispatch correctly in a freshly deserialized context too).

**Upstream fixture / engine change required:** None.

---

## Summary table

| Spec ID      | Name                                             | Priority | Kind                        | Blocked on upstream? |
|--------------|--------------------------------------------------|----------|-----------------------------|----------------------|
| SPEC-MA-01   | `getME(.,"Zt")` dispatch bug fix                 | **P0**   | unit + parity-vs-lme4       | No                   |
| SPEC-MA-02   | `getME(.,"Lambdat")` dispatch bug fix            | **P0**   | unit + parity-vs-lme4       | No                   |
| SPEC-MA-03   | `getME(.,"sigma")` wired to scalar              | P1       | unit + error-message        | No                   |
| SPEC-MA-04   | `ngrps()` named integer parity                   | **P0**   | parity-vs-lme4              | No                   |
| SPEC-MA-05   | `weights()` all-ones vector (not NULL)           | **P0**   | parity-vs-lme4              | No                   |
| SPEC-MA-06   | `terms.mm_lmm` method exists                     | P1       | unit + integration          | No                   |
| SPEC-MA-07   | `hatvalues()` method and numeric parity          | P1       | unit + numerical-tolerance  | Possibly (L/RX FFI)  |
| SPEC-MA-08   | `vcov(correlation=TRUE)` typed refusal or impl.  | P1       | error-message + parity      | No                   |
| SPEC-MA-09   | `getME(.,"u")` and `"b"` stacked RE modes        | P1       | unit + numerical-tolerance  | Check FFI contract   |
| SPEC-MA-10   | `getME(.,"Gp")` group-pointer vector             | P2       | unit + parity-vs-lme4       | No                   |
| SPEC-MA-11   | `getME(.,"devcomp")` deviance decomposition      | P1       | unit + parity-vs-lme4       | No                   |
| SPEC-MA-12   | `getME(.,"lower")` theta lower-bounds            | P2       | unit + parity-vs-lme4       | No                   |
| SPEC-MA-13   | `getME(.,"N"/"n"/"p"/"q")` dimension scalars     | P2       | unit + parity-vs-lme4       | No                   |
| SPEC-MA-14   | `getME(.,"n_rtrms"/"n_rfacs")` RE counts         | P2       | unit + parity-vs-lme4       | No                   |
| SPEC-MA-15   | `getME(.,"REML")` / `isREML()` flag              | P2       | unit + parity-vs-lme4       | No                   |
| SPEC-MA-16   | `model.matrix(type="randomListRaw")`             | P2       | unit                        | No                   |
| SPEC-MA-17   | `getME(.,"Ztlist")` per-term Z' list             | P2       | unit + parity-vs-lme4       | No                   |
| SPEC-MA-18   | Unknown-name `getME` typed error + regression    | P1       | error-message               | No                   |
| SPEC-MA-19   | `getME(.,"ALL")` omnibus dump                    | P2       | unit + parity-vs-lme4       | No                   |
| SPEC-MA-20   | Cross-session revival of matrix accessors        | P1       | cross-session-revival       | No                   |

**Count:** 20 test specifications (some specs contain multiple `test_that` blocks;
the count tracks distinct named specs, not the number of `expect_*` calls).

**Top priority:** SPEC-MA-01 (`getME(.,"Zt")` dispatch bug) — it is a P0 blocker
because the name is advertised as available, already returns the wrong error class
(not `mm_arg_error` but a base-R crash), is a one-line fix (`Matrix::t()` instead
of `t()`), and blocking on it means every downstream user who calls `getME` for
the transpose sparse matrices gets a cryptic unhelpful error — the exact failure
mode the "clearer errors" mission is designed to prevent.
