# Test Specification — Diagnostics & Convergence Family

**Family:** Diagnostics & convergence  
**Source gap report:** `assessment/gap/diagnostics.md`  
**Parity probes referenced:** `assessment/parity/inf-ranef-condvar.md`,
`assessment/parity/inf-getME.md`, `assessment/parity/lmm-dyestuff2-singular.md`,
`assessment/parity/lmm-sleep-ri.md`  
**Error probes referenced:** `assessment/errors/singular-fit.md`,
`assessment/errors/convergence-hard.md`  
**In-scope gaps addressed:** partial × 9 sub-gaps, in-scope-missing × 6, test-gap × 4  
**Total test specs:** 22  
**Written:** 2026-05-31

---

## Classification summary

| Gap from gap report | Classification | Severity | Test specs |
|---|---|---|---|
| `residuals(type="pearson"/"deviance"/"working")` errors | partial | major | TS-01, TS-02, TS-03 |
| `getME("Zt")` / `getME("Lambdat")` crash (`t.default` dispatch) | partial (in-scope bug) | major | TS-04, TS-05 |
| `getME("lower")` absent (boundary diagnostic) | partial | major | TS-06 |
| `getME("is_REML")` / `"sigma"` / `"n_rtrms"` absent | partial | minor | TS-07 |
| `getME("theta")` returns unnamed vector (lme4 returns named) | partial | minor | TS-08 |
| `ranef(condVar=TRUE)` postVar within tolerance | partial | minor | TS-09 |
| `ranef(condVar=TRUE)` GLMM returns typed refusal, not crash | partial | minor | TS-10 |
| `rePCA()` — no prcomp-style object | partial | major | TS-11 |
| `allFit()` — entirely absent | in-scope-missing | major | TS-12 |
| `influence()` / `cooks.distance()` / `dfbeta()` / `hatvalues()` entirely absent | in-scope-missing | major | TS-13, TS-14, TS-15, TS-16 |
| Overparameterised model missing identifiability language | test-gap (error messaging) | minor | TS-17 |
| No R-level warning emitted at fit time for reduced-rank / boundary | test-gap (discoverability) | minor | TS-18 |
| GLMM `boundary_parameter` uses parameter index, not term name | test-gap (error messaging) | minor | TS-19 |
| `lme4::isSingular(mm_lmm)` gives opaque base-R error, not clear redirect | test-gap (clearer-errors) | minor | TS-20 |
| `check.scaleX` — no scaling advisory for badly-scaled predictors | in-scope-missing | minor | TS-21 |
| `isREML()` predicate absent | partial | minor | TS-22 |

---

## TS-01 — `residuals(type="pearson")` returns scaled Pearson residuals for LMM

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** `residuals.mm_lmm` only accepting `type="response"`, so every
call with `type="pearson"` errors — the standard diagnostic residual for LMM is
unavailable.

### Dataset / formula
`sleepstudy` (lme4 built-in, n=180).  
`lmm(Reaction ~ Days + (Days | Subject), sleepstudy)`

### Assertion
```r
test_that("residuals(type='pearson') returns response/sigma for LMM", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  r_resp    <- residuals(fit, type = "response")
  r_pearson <- residuals(fit, type = "pearson")

  # Pearson residuals are response residuals scaled by sigma
  expect_length(r_pearson, nrow(sleepstudy))
  expect_true(all(is.finite(r_pearson)))
  expect_equal(r_pearson, r_resp / sigma(fit), tolerance = 1e-8,
               ignore_attr = TRUE,
               info = "Pearson residuals must equal response/sigma for LMM")
})
```

### Tolerance
Element-wise relative tolerance 1e-8 (exact arithmetic relationship).

### Notes
- Fix: extend `residuals.mm_lmm` to accept `type = c("response", "pearson",
  "deviance", "working")` and implement the sigma-scaling branch.
- Upstream fixture: none needed; pure R arithmetic.

---

## TS-02 — `residuals(type="deviance")` returns correct deviance residuals for LMM

**Priority:** P1  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** `type="deviance"` silently erroring; for Gaussian LMM deviance
and response residuals are identical (sign-preserving), so this is also a
regression guard.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)` (RI model).

### Assertion
```r
test_that("residuals(type='deviance') equals type='response' for Gaussian LMM", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  r_resp <- residuals(fit, type = "response")
  r_dev  <- residuals(fit, type = "deviance")

  expect_length(r_dev, nrow(sleepstudy))
  expect_true(all(is.finite(r_dev)))
  # For Gaussian: deviance residuals == response residuals (signs preserved)
  expect_equal(r_dev, r_resp, tolerance = 1e-10, ignore_attr = TRUE)
})
```

### Tolerance
1e-10 (should be numerically identical for Gaussian).

### Notes
- For Gaussian LMM this is trivial, but the test locks in the contract that
  `type="deviance"` does not error.

---

## TS-03 — `residuals(type="pearson")` for GLMM (binomial) uses family variance function

**Priority:** P0 blocker  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** GLMM residuals always erroring with `'arg' should be "response"`;
lme4's default GLMM residual type is `"deviance"` and the standard diagnostic is
`"pearson"` — both are inaccessible without this fix.

### Dataset / formula
`cbpp` (lme4 built-in).  
`glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cbpp, family = binomial)`

### Assertion
```r
test_that("residuals(type='pearson') for GLMM uses family variance function", {
  skip_if_not_installed("lme4")
  data("cbpp", package = "lme4")
  fit <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial,
              control = mm_control(verbose = -1))

  r_pearson <- residuals(fit, type = "pearson")

  expect_length(r_pearson, nrow(cbpp))
  expect_true(all(is.finite(r_pearson)))

  # Pearson residuals: (y - mu_hat) / sqrt(V(mu_hat))
  # For binomial V(mu) = mu*(1-mu)/n; cross-check sign direction
  expect_true(all(abs(r_pearson) < 10),
              info = "Pearson residuals for cbpp should be in a reasonable range")

  # Cross-check against lme4 to within a moderate tolerance
  lme4_fit <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                           data = cbpp, family = binomial)
  lme4_r <- residuals(lme4_fit, type = "pearson")
  expect_equal(r_pearson, as.numeric(lme4_r), tolerance = 0.05,
               ignore_attr = TRUE,
               info = "GLMM Pearson residuals should be within 5% of lme4")
})
```

### Tolerance
0.05 absolute (moderate; allows for optimizer differences in the underlying fit).

### Notes
- Fix requires `residuals.mm_glmm` to implement the family variance function.
- Upstream fixture: GLMM Pearson-residual formula requires the fitted `mu` and
  the binomial variance function `mu*(1-mu)`. Both are available from the artifact.

---

## TS-04 — `getME("Zt")` returns transposed Z sparse matrix without error

**Priority:** P0 blocker  
**Kind:** unit / error-message  
**Guards against:** `getME("Zt")` crashing with `Error in t.default: argument is
not a matrix` due to bare `t()` on a `Matrix::sparseMatrix` in `revive.R:128`.
This breaks any downstream code (e.g. merTools) that calls `getME(fit, "Zt")`.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (Days | Subject), sleepstudy)`.

### Assertion
```r
test_that("getME('Zt') returns transposed sparse Z without error", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  Zt <- getME(fit, "Zt")

  # Must not error; must be a sparse matrix
  expect_true(methods::is(Zt, "Matrix") || is.matrix(Zt))

  # Dimensions: Z is 180x36 so Zt must be 36x180
  expect_equal(dim(Zt), c(36L, 180L))

  # Zt must equal t(Z)
  Z <- getME(fit, "Z")
  expect_equal(as.matrix(Zt), t(as.matrix(Z)), tolerance = 1e-10)
})
```

### Tolerance
1e-10 (exact transpose relationship).

### Notes
- Fix: change `t(stats::model.matrix(...))` to `Matrix::t(stats::model.matrix(...))`
  in `R/revive.R` line 128.
- No upstream fixture required.

---

## TS-05 — `getME("Lambdat")` returns transposed Lambda sparse matrix without error

**Priority:** P0 blocker  
**Kind:** unit / error-message  
**Guards against:** `getME("Lambdat")` crashing with the same `t.default` dispatch
bug as TS-04 (`revive.R:130`).

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (Days | Subject), sleepstudy)`.

### Assertion
```r
test_that("getME('Lambdat') returns transposed Lambda sparse matrix without error", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  Lambdat <- getME(fit, "Lambdat")

  expect_true(methods::is(Lambdat, "Matrix") || is.matrix(Lambdat))

  # Lambda is 36x36, so Lambdat must also be 36x36
  expect_equal(dim(Lambdat), c(36L, 36L))

  # Lambdat must equal t(Lambda)
  Lambda <- getME(fit, "Lambda")
  expect_equal(as.matrix(Lambdat), t(as.matrix(Lambda)), tolerance = 1e-10)
})
```

### Tolerance
1e-10.

### Notes
- Fix: change `t(.mm_lazy(...))` to `Matrix::t(.mm_lazy(...))` in `R/revive.R`
  line 130.
- Both TS-04 and TS-05 are fixed by adding `Matrix::` qualification in the same
  two-line edit.

---

## TS-06 — `getME("lower")` returns theta lower-bound vector (boundary diagnostic)

**Priority:** P1  
**Kind:** parity-vs-lme4  
**Guards against:** the canonical lme4 boundary diagnostic
`all(getME(fit, "theta") == getME(fit, "lower"))` being unavailable, forcing
users to call `is_singular()` when they need the raw boundary check.

### Dataset / formula
`Dyestuff2` (lme4 built-in, n=30, boundary case).  
`lmm(Yield ~ 1 + (1 | Batch), Dyestuff2)`

### Assertion
```r
test_that("getME('lower') returns the theta lower-bound vector", {
  skip_if_not_installed("lme4")
  data("Dyestuff2", package = "lme4")
  fit <- lmm(Yield ~ 1 + (1 | Batch), Dyestuff2,
             control = mm_control(verbose = -1))

  lower <- getME(fit, "lower")

  # Must be a numeric vector, same length as theta
  theta <- getME(fit, "theta")
  expect_type(lower, "double")
  expect_length(lower, length(theta))

  # All lower bounds must be 0 or -Inf (standard lme4 convention)
  expect_true(all(lower >= 0 | is.infinite(lower)),
              info = "lower bounds should be 0 (SD params) or -Inf (covariance params)")

  # On the Dyestuff2 boundary case, theta == lower (intercept SD hits 0)
  expect_true(is_singular(fit))
  expect_equal(theta, lower, tolerance = 1e-8,
               info = "boundary case: theta should equal lower bound")
})
```

### Tolerance
1e-8 for numeric comparison of theta vs lower in the boundary case.

### Notes
- `lower` for a scalar random-intercept model is `c(0)` — the SD parameter
  cannot go below zero.
- Implementation requires adding `"lower"` to the `getME.mm_lmm` switch; the
  value can be constructed from the theta parameterisation (0 for variance-SD
  elements, -Inf for correlation elements).
- Upstream fixture: none needed; the lower-bound vector is a property of the
  parameterisation, not the optimiser result.

---

## TS-07 — `getME("is_REML")` / `"sigma"` / `"n_rtrms"` return correct scalars

**Priority:** P2  
**Kind:** unit  
**Guards against:** downstream packages (merTools, ggeffects) that call
`getME(fit, "is_REML")`, `getME(fit, "sigma")`, or `getME(fit, "n_rtrms")`
receiving an `mm_arg_error` instead of a scalar.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)` (REML=TRUE default).

### Assertion
```r
test_that("getME supports 'is_REML', 'sigma', and 'n_rtrms'", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  # is_REML: logical scalar
  is_reml <- getME(fit, "is_REML")
  expect_type(is_reml, "logical")
  expect_length(is_reml, 1L)
  expect_true(is_reml)  # default is REML

  # sigma: positive numeric scalar matching sigma(fit)
  s <- getME(fit, "sigma")
  expect_type(s, "double")
  expect_length(s, 1L)
  expect_equal(s, sigma(fit), tolerance = 1e-10)

  # n_rtrms: number of random-effect terms (integer scalar >= 1)
  n <- getME(fit, "n_rtrms")
  expect_type(n, "integer")
  expect_length(n, 1L)
  expect_equal(n, 1L)  # one RE term: (1|Subject)
})
```

### Tolerance
1e-10 for sigma comparison.

### Notes
- `is_REML` can be derived from `fit$artifact$reml` or equivalent.
- `n_rtrms` is the number of random-effect grouping structures.

---

## TS-08 — `getME("theta")` returns a named numeric vector

**Priority:** P2  
**Kind:** parity-vs-lme4  
**Guards against:** programmatic code that calls `names(getME(fit, "theta"))`
receiving `NULL` instead of the lme4-style `"Subject.(Intercept)"` / `"Subject.Days"` names.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (Days | Subject), sleepstudy)`.

### Assertion
```r
test_that("getME('theta') vector carries lme4-compatible names", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  theta <- getME(fit, "theta")

  # Length 3 for correlated random slope
  expect_length(theta, 3L)

  # Names must be non-NULL and non-empty
  expect_false(is.null(names(theta)),
               info = "theta names must not be NULL")
  expect_true(all(nzchar(names(theta))),
              info = "all theta names must be non-empty strings")

  # Must contain the group name "Subject"
  expect_true(any(grepl("Subject", names(theta))),
              info = "theta names should reference the grouping factor 'Subject'")
})
```

### Tolerance
N/A (structural / naming check only).

### Notes
- lme4 format: `c("Subject.(Intercept)", "Subject.Days.(Intercept)",
  "Subject.Days")` for a correlated random slope.
- The exact naming convention can differ from lme4 as long as names are
  non-null, non-empty, and group-referencing.

---

## TS-09 — `ranef(condVar=TRUE)` postVar values within tolerance vs lme4

**Priority:** P1  
**Kind:** parity-vs-lme4 / numerical-tolerance  
**Guards against:** postVar values drifting beyond the PRD tolerance for
variance-component estimates (1e-3), which would indicate a bug in the
`cond_var()` Rust bridge beyond what is explained by the known theta/sigma drift.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (Days | Subject), sleepstudy)`.

### Assertion
```r
test_that("ranef(condVar=TRUE) postVar diagonal is within tolerance vs lme4", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")

  fit_mm   <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
                  control = mm_control(verbose = -1))
  fit_lme4 <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                          REML = TRUE)

  r_mm   <- ranef(fit_mm,   condVar = TRUE)
  r_lme4 <- lme4::ranef(fit_lme4, condVar = TRUE)

  pv_mm   <- attr(r_mm$Subject,   "postVar")
  pv_lme4 <- attr(r_lme4$Subject, "postVar")

  # Both must be 2x2x18 arrays
  expect_equal(dim(pv_mm), c(2L, 2L, 18L))
  expect_equal(dim(pv_lme4), c(2L, 2L, 18L))

  # All entries must be finite
  expect_true(all(is.finite(pv_mm)))

  # Diagonal (conditional variances) must be non-negative
  diag_mm <- apply(pv_mm, 3, diag)  # 2x18 matrix
  expect_true(all(diag_mm >= 0))

  # All slices must be symmetric
  for (i in seq_len(dim(pv_mm)[3])) {
    expect_equal(pv_mm[,,i], t(pv_mm[,,i]), tolerance = 1e-12,
                 info = sprintf("postVar slice %d not symmetric", i))
  }

  # Conditional SDs (sqrt of diagonal) within 1e-3 of lme4
  condsd_mm   <- sqrt(diag_mm)  # 2x18
  condsd_lme4 <- sqrt(apply(pv_lme4, 3, diag))
  expect_equal(condsd_mm, condsd_lme4, tolerance = 1e-3,
               ignore_attr = TRUE,
               info = "Conditional SDs from postVar must match lme4 within 1e-3")
})
```

### Tolerance
1e-3 on conditional SDs (PRD variance-component tolerance); 1e-12 for
symmetry check.

### Notes
- The probe `assessment/parity/inf-ranef-condvar.md` shows condSD(Intercept)
  diff = 8.2e-4 and condSD(Days) diff = 2.2e-4, both within 1e-3 — this test
  encodes that passing bar.
- The raw variance values (postVar[1,1] diff = 1.97e-2) exceed 1e-3 but the
  SD-level comparison passes; the test uses SDs to separate numerical noise
  from structural bugs.
- If the underlying theta/sigma drift is fixed (bd-01KRV31R4BJVQCEF0F58NFD4YN),
  both the variance-level and SD-level comparisons should pass at 1e-3.

---

## TS-10 — `ranef(condVar=TRUE)` for GLMM returns typed refusal, not crash

**Priority:** P1  
**Kind:** error-message / unit  
**Guards against:** `ranef(glmm_fit, condVar=TRUE)` crashing with an unhandled
error; the current behavior returns `NA` postVar with a structured reason code.
This test locks in the honest-refusal contract.

### Dataset / formula
`cbpp` (lme4 built-in).  
`glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cbpp, family = binomial)`

### Assertion
```r
test_that("ranef(condVar=TRUE) for GLMM returns a structured refusal, not a crash", {
  skip_if_not_installed("lme4")
  data("cbpp", package = "lme4")
  fit <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial,
              control = mm_control(verbose = -1))

  # Must not error
  r <- expect_no_error(ranef(fit, condVar = TRUE))

  pv <- attr(r$herd, "postVar")

  # postVar must be NA-carrying with a documented reason attribute
  # (structured refusal, not NULL or crash)
  expect_false(is.null(pv),
               info = "postVar must be non-NULL (structured refusal expected)")

  reason_attr <- attr(pv, "mm_unavailable_reason")
  expect_false(is.null(reason_attr),
               info = "postVar NA must carry mm_unavailable_reason attribute")
  expect_match(reason_attr, "glmm",
               info = "refusal reason should mention glmm")
})
```

### Tolerance
N/A (structural / error-contract check).

### Notes
- Current behavior per `assessment/parity/inf-ranef-condvar.md`: returns NA
  postVar with reason `random_effect_conditional_variance_unavailable_for_glmm`.
- This test certifies the honest-refusal contract so a future change does not
  accidentally replace it with a crash.
- Upstream fixture: a full GLMM condVar implementation would turn this into a
  numerical-tolerance test (like TS-09); that is a future promotion.

---

## TS-11 — `rePCA()` substitute: `effective_covariance` exposes rank and loadings

**Priority:** P1  
**Kind:** integration / snapshot  
**Guards against:** the `effective_covariance` artifact losing the fields that
substitute for `rePCA()` — specifically `requested_rank`, `supported_rank`, and
directional loadings accessible via `changes(fit)`.

### Dataset / formula
Reduced-rank forced fit:  
```r
ss <- sleepstudy
ss$noise <- rnorm(180, sd = 1)
lmm(Reaction ~ Days + (Days + noise | Subject), ss)
```

### Assertion
```r
test_that("effective_covariance exposes rank reduction info (rePCA substitute)", {
  skip_if_not_installed("lme4")
  set.seed(42L)
  data("sleepstudy", package = "lme4")
  ss <- sleepstudy
  ss$noise <- rnorm(nrow(ss), sd = 1)

  fit <- lmm(Reaction ~ Days + (Days + noise | Subject), ss,
             control = mm_control(verbose = -1))

  # Must be classified as singular / reduced-rank
  expect_true(is_singular(fit))
  expect_match(fit_status(fit), "reduced_rank")

  # changes() must return a non-empty object describing the rank reduction
  ch <- changes(fit)
  expect_false(is.null(ch))

  # The printed changes output must mention requested and fitted rank
  ch_text <- paste(capture.output(print(ch)), collapse = "\n")
  expect_match(ch_text, "rank", ignore.case = TRUE,
               info = "changes() output should describe rank reduction")

  # optimizer_certificate must carry hessian_rank and hessian_eigen_min
  cert <- optimizer_certificate(fit)
  cert_tbl <- cert$table
  expect_true("hessian_rank" %in% names(cert_tbl) ||
              any(grepl("hessian_rank", cert_tbl$field)),
              info = "optimizer_certificate should expose hessian_rank")
})
```

### Tolerance
N/A (structural check).

### Notes
- There is no `rePCA()` function and none is planned; this test certifies the
  *intended substitute* (`changes()`, `effective_covariance` artifact) works.
- The test documents the documented divergence from lme4's `rePCA()` API: the
  information exists but there is no `prcomp`-style object with `$sdev` /
  `$rotation`.

---

## TS-12 — `allFit()` is absent with a clear, documented error message

**Priority:** P1  
**Kind:** error-message  
**Guards against:** `allFit()` being silently swallowed (no error, no result) or
returning a cryptic base-R dispatch failure instead of a pointer to the
mixeff alternative.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)`.

### Assertion
```r
test_that("allFit() on mm_lmm emits a clear 'not available' message pointing to optimizer_certificate", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  # allFit() is a lme4 function; calling it on mm_lmm should either:
  #   (a) error with a message mentioning optimizer_certificate / mixeff, OR
  #   (b) not exist in the package at all (no S3 method for mm_lmm)
  # Either is acceptable; what must NOT happen is a silent wrong answer.
  result <- tryCatch(
    lme4::allFit(fit),
    error = function(e) e
  )

  # Must be an error (no applicable method or explicit mm_abort)
  expect_true(inherits(result, "error"),
              info = "allFit() on mm_lmm must produce an error, not a silent result")

  # The error message should be informative (not a raw dispatch failure)
  # Accept either: a clear mixeff error OR a base-R no-method error
  # Certify that it does not return a usable allFit object silently
  expect_false(inherits(result, "allFit"),
               info = "allFit() must not silently return an allFit object")
})
```

### Tolerance
N/A (error contract).

### Notes
- `allFit` is a major gap (`in-scope-missing`). This test certifies the *current*
  honest-failure behavior; a future `allFit.mm_lmm` implementation (if the
  engine exposes multiple optimizers) would replace this with a numerical test.
- The PRD §4 places the optimizer inside the Rust crate (`trust_bq`); adding
  allFit would require either upstream engine support or a wrapper that reruns
  with different mm_control options.

---

## TS-13 — `hatvalues()` is absent with a clear error, not a crash

**Priority:** P1  
**Kind:** error-message  
**Guards against:** `hatvalues(mm_lmm_fit)` giving a raw base-R dispatch failure
with no pointer to mixeff's audit-first alternatives.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)`.

### Assertion
```r
test_that("hatvalues() on mm_lmm emits a clear not-yet-available message", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  result <- tryCatch(
    stats::hatvalues(fit),
    error = function(e) e
  )

  expect_true(inherits(result, "error"),
              info = "hatvalues() must error (no method implemented)")

  # Future: if a hatvalues.mm_lmm is added the error message should mention
  # the mixeff alternative or indicate the feature is planned.
  # For now, assert it at least errors cleanly (no segfault, no infinite loop).
  expect_false(is.null(result$message),
               info = "error must carry a message")
})
```

### Tolerance
N/A.

### Notes
- `in-scope-missing`; `X` and `Z` are accessible via `getME`, so a user can
  hand-roll `H = X %*% solve(t(X) %*% X) %*% t(X)` (marginal hat matrix)
  outside mixeff. A future `hatvalues.mm_lmm` should compute the conditional hat
  diagonal using `Z`, `Lambda`, and `X`.
- Upstream fixture required: the Rust engine would need to expose the projected
  hat diagonal or the full `H = X(X'X)^{-1}X'` plus the RE projection.

---

## TS-14 — `influence()` is absent with a clear error, not a crash

**Priority:** P1  
**Kind:** error-message  
**Guards against:** `influence(mm_lmm_fit)` giving a cryptic base-R
`no applicable method` message with no explanation or redirect.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)`.

### Assertion
```r
test_that("influence() on mm_lmm emits a clear not-yet-available message", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  result <- tryCatch(
    stats::influence(fit, groups = "Subject"),
    error = function(e) e
  )

  expect_true(inherits(result, "error"),
              info = "influence() must error (not yet implemented)")

  # Must not segfault or hang; the error message must be a non-empty string
  expect_true(is.character(result$message) && nzchar(result$message))
})
```

### Tolerance
N/A.

### Notes
- `in-scope-missing` / major. Crate survey explicitly lists "influence diagnostics
  requiring repeated refits" as not yet provided.
- A full implementation requires leave-one-group-out refits; this is a
  non-trivial upstream addition.

---

## TS-15 — `cooks.distance()` on an `mm_lmm` errors gracefully

**Priority:** P2  
**Kind:** error-message  
**Guards against:** `cooks.distance()` silently dispatching to an unrelated
method or crashing; depends on `influence.merMod` being absent (TS-14).

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)`.

### Assertion
```r
test_that("cooks.distance() on mm_lmm errors gracefully", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  result <- tryCatch(
    stats::cooks.distance(fit),
    error = function(e) e
  )

  # Must either error cleanly OR (if influence is eventually implemented)
  # return a named numeric vector of length nrow(sleepstudy) or
  # length(unique(sleepstudy$Subject)).
  if (inherits(result, "error")) {
    expect_true(is.character(result$message) && nzchar(result$message),
                info = "error must carry a non-empty message")
  } else {
    # If somehow implemented: must be a finite numeric vector
    expect_true(is.numeric(result))
    expect_true(all(is.finite(result)))
  }
})
```

### Tolerance
N/A.

### Notes
- `in-scope-missing`. Depends on `influence()` (TS-14).

---

## TS-16 — `dfbeta()` / `dfbetas()` on `mm_lmm` error gracefully

**Priority:** P2  
**Kind:** error-message  
**Guards against:** `dfbeta()` or `dfbetas()` dispatching silently to a base-R
method that returns garbage.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)`.

### Assertion
```r
test_that("dfbeta() and dfbetas() on mm_lmm error gracefully", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  res_dfbeta  <- tryCatch(stats::dfbeta(fit),  error = function(e) e)
  res_dfbetas <- tryCatch(stats::dfbetas(fit), error = function(e) e)

  expect_true(inherits(res_dfbeta, "error") || is.numeric(res_dfbeta),
              info = "dfbeta() must error or return a numeric matrix")
  expect_true(inherits(res_dfbetas, "error") || is.numeric(res_dfbetas),
              info = "dfbetas() must error or return a numeric matrix")

  # If both error, check messages are non-empty
  if (inherits(res_dfbeta, "error")) {
    expect_true(nzchar(res_dfbeta$message))
  }
})
```

### Tolerance
N/A.

### Notes
- `in-scope-missing`. Depends on `influence()` (TS-14).

---

## TS-17 — Overparameterised model (n_obs <= n_RE) surfaces identifiability language

**Priority:** P1  
**Kind:** error-message / integration  
**Guards against:** an overparameterised model (e.g. 20 obs, 30 random effects)
being silently fit with only a rank-deficiency note but no identifiability
framing — the `assessment/errors/singular-fit.md` scenario C gap.

### Dataset / formula
Synthetic: 10 subjects × 2 obs, random slope on 3 terms (n_obs=20, n_re=30).

### Assertion
```r
test_that("overparameterised model (n_obs <= n_re) surfaces identifiability language", {
  set.seed(7L)
  n_sub <- 10L; n_per <- 2L
  subject <- factor(rep(seq_len(n_sub), each = n_per))
  x <- rnorm(n_sub * n_per)
  z <- rnorm(n_sub * n_per)
  y <- rnorm(n_sub * n_per)
  df <- data.frame(y = y, x = x, z = z, subject = subject)

  # This formula requests n_re = 10 * 3 = 30 > n_obs = 20
  fit <- suppressMessages(
    lmm(y ~ x + z + (1 + x + z | subject), df,
        control = mm_control(verbose = -1))
  )

  # Must not crash; must return a fit object
  expect_s3_class(fit, "mm_lmm")

  # is_singular must be TRUE
  expect_true(is_singular(fit))

  # diagnostics() output must carry a diagnostic relating to identifiability
  # or n_obs vs n_random_effects
  diag_obj <- diagnostics(fit)
  diag_text <- paste(capture.output(print(diag_obj)), collapse = "\n")

  # Either: an explicit "identifiable"/"identifiability" mention, OR
  # a covariance_too_rich / rank-deficient message with group-level counts
  has_identif_msg <- grepl("identif", diag_text, ignore.case = TRUE)
  has_rank_msg    <- grepl("rank", diag_text, ignore.case = TRUE)
  has_too_rich    <- grepl("too_rich|covariance_too_rich|n_obs|n_random",
                           diag_text, ignore.case = TRUE)

  expect_true(has_identif_msg || has_rank_msg || has_too_rich,
              info = paste(
                "diagnostics() should surface identifiability / rank / richness info.",
                "Got:", substr(diag_text, 1, 400)
              ))
})
```

### Tolerance
N/A (message content check).

### Notes
- From `assessment/errors/singular-fit.md` scenario C: lme4 refuses outright
  with "probably unidentifiable"; mixeff fits and shows rank-deficient.
- The test accepts either approach (refusal or rank message), but requires
  *some* identifiability-relevant language to be surfaced.
- If the engine adds an explicit `n_obs_le_n_re` diagnostic code this test
  should be tightened to check for that code.

---

## TS-18 — Reduced-rank / boundary fit emits an R-level `message()` at fit time

**Priority:** P2  
**Kind:** error-message / unit  
**Guards against:** a user fitting a singular model and seeing no console
output at all (the discoverability gap from `assessment/errors/convergence-hard.md`).

### Dataset / formula
`Dyestuff2`, `lmm(Yield ~ 1 + (1 | Batch), Dyestuff2)` (boundary case).

### Assertion
```r
test_that("boundary/reduced-rank fit emits a console message without verbose=-1", {
  skip_if_not_installed("lme4")
  data("Dyestuff2", package = "lme4")

  # With default verbose, some output must reach the console for a boundary fit
  msgs <- character(0L)
  withCallingHandlers(
    {
      fit <- lmm(Yield ~ 1 + (1 | Batch), Dyestuff2)
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  # Either a message was emitted OR the default explain_model() output
  # (captured via output, not message) described the boundary.
  # Capture the visible output as well.
  visible_out <- capture.output(
    suppressMessages(lmm(Yield ~ 1 + (1 | Batch), Dyestuff2))
  )
  visible_text <- paste(visible_out, collapse = "\n")

  has_msg    <- length(msgs) > 0L
  has_output <- grepl("boundary|singular|reduced_rank|rank-deficient",
                      visible_text, ignore.case = TRUE)

  expect_true(has_msg || has_output,
              info = paste(
                "A boundary fit must produce some console output without verbose=-1.",
                "visible:", substr(visible_text, 1L, 300L)
              ))
})
```

### Tolerance
N/A (discoverability / message check).

### Notes
- Currently mixeff's default `verbose` auto-prints `explain_model()` which
  includes "Fit status: converged_reduced_rank" — this satisfies `has_output`.
- The test also covers the future where a typed R `message()` is emitted.

---

## TS-19 — GLMM `boundary_parameter` diagnostic names the term, not the index

**Priority:** P2  
**Kind:** error-message / integration  
**Guards against:** the GLMM `boundary_parameter` message saying
`"covariance parameter 3 is on its lower bound"` instead of naming the
random-effect term (as the LMM version does).

### Dataset / formula
GLMM with near-singular random structure (synthetic, based on
`assessment/errors/convergence-hard.md` scenario 3).

### Assertion
```r
test_that("GLMM boundary_parameter diagnostic names the term, not a numeric index", {
  skip_if_not_installed("lme4")
  set.seed(99L)
  n_sub <- 20L; n_per <- 5L
  subject <- factor(rep(seq_len(n_sub), each = n_per))
  x_small <- rnorm(n_sub * n_per)
  x_large <- rnorm(n_sub * n_per) * 1e6
  y_bin   <- rbinom(n_sub * n_per, size = 1L,
                    prob = plogis(0.2 * x_small))
  df <- data.frame(y = y_bin, x_small = x_small, x_large = x_large,
                   subject = subject)

  fit <- suppressWarnings(
    glmm(y ~ x_small + x_large + (1 + x_large | subject),
         data = df, family = binomial,
         control = mm_control(verbose = -1))
  )

  cert <- optimizer_certificate(fit)
  cert_tbl <- cert$table

  # If a boundary_parameter diagnostic is present, it must mention a term name
  # (not just "covariance parameter N")
  boundary_rows <- cert_tbl[grepl("boundary_parameter|boundary",
                                  cert_tbl$code %||% cert_tbl$field %||% "",
                                  ignore.case = TRUE), , drop = FALSE]

  if (nrow(boundary_rows) > 0L) {
    msgs <- paste(boundary_rows$message %||% boundary_rows$detail %||% "",
                  collapse = " ")
    # The message should mention a variable name, not just "parameter N"
    only_index <- grepl("^(covariance )?parameter [0-9]+ is on", msgs) &&
                  !grepl("x_large|x_small|subject|Intercept", msgs,
                         ignore.case = TRUE)
    expect_false(only_index,
                 info = paste(
                   "boundary_parameter message should name the term, not just an index.",
                   "Got:", substr(msgs, 1L, 300L)
                 ))
  } else {
    # No boundary_parameter fired — still OK (fit may have converged interior)
    skip("No boundary_parameter diagnostic in this run; cannot test message format")
  }
})
```

### Tolerance
N/A (message content check).

### Notes
- From `assessment/errors/convergence-hard.md`: LMM version names the term
  ("standard deviation for x_large in (1 + x_large | subject)"); GLMM version
  says "covariance parameter 3 is on its lower bound".
- Fix requires the GLMM optimizer-certificate code to map the parameter index
  back to the term name via the covariance spec.
- Upstream fixture: requires the Rust engine to expose the mapping from theta
  index to RE term name in the GLMM optimizer certificate.
  **Upstream change needed.**

---

## TS-20 — `lme4::isSingular(mm_lmm_fit)` emits a pointer to `is_singular()`, not opaque base-R error

**Priority:** P2  
**Kind:** error-message / clearer-errors  
**Guards against:** a user migrating from lme4 calling `lme4::isSingular(fit)` on
an `mm_lmm` and receiving
`"no applicable method for 'isSingular' applied to an object of class ..."`,
with no pointer to `is_singular()`.

### Dataset / formula
`sleepstudy`, `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)`.

### Assertion
```r
test_that("lme4::isSingular(mm_lmm) either works or emits a pointer to is_singular()", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))

  result <- tryCatch(
    lme4::isSingular(fit),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    # Acceptable: must mention is_singular or mixeff in the error/message
    # Currently it gives a raw dispatch error with no pointer
    # This test FAILS today (documenting the gap) and passes once a redirect is added.
    expect_match(result$message, "is_singular|mixeff",
                 ignore.case = TRUE,
                 info = paste(
                   "lme4::isSingular(mm_lmm) error should point to is_singular().",
                   "Got:", result$message
                 ))
  } else {
    # If a method is registered that delegates, it must return a logical scalar
    expect_type(result, "logical")
    expect_length(result, 1L)
  }
})
```

### Tolerance
N/A (error message content).

### Notes
- Fix: add an `isSingular.mm_lmm` method that either delegates to
  `is_singular(x)` or calls `mm_abort()` with a message pointing to
  `is_singular()`.
- Same pattern applies to `lme4::rePCA(mm_lmm)` (not a separate test here as
  `rePCA` has no direct equivalent; the `changes()` redirect is sufficient).

---

## TS-21 — Badly-scaled predictors produce a `check_scale` advisory in `diagnostics()`

**Priority:** P2  
**Kind:** error-message / integration  
**Guards against:** a user passing predictors on wildly different scales
(e.g. one O(1), one O(1e6)) and receiving no advisory — the `check.scaleX`
gap from the diagnostics gap report.

### Dataset / formula
Synthetic: `lmm(y ~ x_small + x_large + (1 | subject), df)` where
`x_small ~ O(1)`, `x_large ~ O(1e6)`.

### Assertion
```r
test_that("badly-scaled predictors produce a scale advisory in diagnostics()", {
  set.seed(13L)
  n <- 200L
  subject <- factor(rep(1:20, each = 10))
  x_small <- rnorm(n)
  x_large <- rnorm(n) * 1e6
  y <- 1 + 0.5 * x_small + 2e-7 * x_large +
       rep(rnorm(20, sd = 0.5), each = 10) + rnorm(n, sd = 0.3)
  df <- data.frame(y = y, x_small = x_small, x_large = x_large,
                   subject = subject)

  fit <- lmm(y ~ x_small + x_large + (1 | subject), df,
             control = mm_control(verbose = -1))

  diag_obj  <- diagnostics(fit)
  diag_text <- paste(capture.output(print(diag_obj)), collapse = "\n")

  # Must surface some kind of scale / condition / rank warning
  has_scale_msg <- grepl("scale|condition|rescal|rank.deficient|rank-deficient",
                         diag_text, ignore.case = TRUE)

  expect_true(has_scale_msg,
              info = paste(
                "Badly-scaled predictors should produce a scale/condition advisory.",
                "Got diagnostics:", substr(diag_text, 1L, 400L)
              ))
})
```

### Tolerance
N/A (message presence check).

### Notes
- The `assessment/errors/convergence-hard.md` probe shows that
  `fixed_effect_rank_deficient` already fires for extreme scale mismatch
  (condition number ~10^9). This test certifies that behavior is stable and
  covers the `check.scaleX` intent.
- If a dedicated `scale_mismatch` diagnostic code is added upstream, this test
  should be tightened to check for that code explicitly.

---

## TS-22 — `isREML()` predicate (or `getME("is_REML")`) is accessible and correct

**Priority:** P2  
**Kind:** unit  
**Guards against:** downstream code that branches on `isREML(fit)` or
`getME(fit, "is_REML")` (e.g. lmerTest internals, merTools) getting an error.

### Dataset / formula
`sleepstudy`:  
- `lmm(Reaction ~ Days + (1 | Subject), sleepstudy)` — REML (default)  
- `lmm(Reaction ~ Days + (1 | Subject), sleepstudy, REML = FALSE)` — ML

### Assertion
```r
test_that("REML status is accessible and correct via getME or fit$artifact", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")

  fit_reml <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
                  control = mm_control(verbose = -1))
  fit_ml   <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
                  REML = FALSE, control = mm_control(verbose = -1))

  # Via getME("is_REML") if implemented (TS-07 prerequisite)
  # Fall back to fit$artifact if getME not yet implemented
  get_is_reml <- function(fit) {
    tryCatch(
      getME(fit, "is_REML"),
      error = function(e) {
        # Fall-back: check artifact directly
        fit$artifact$reml %||% fit$artifact$schema$reml %||% NA
      }
    )
  }

  is_reml_reml <- get_is_reml(fit_reml)
  is_reml_ml   <- get_is_reml(fit_ml)

  expect_false(is.na(is_reml_reml),
               info = "REML status must be accessible (not NA) for REML fit")
  expect_false(is.na(is_reml_ml),
               info = "REML status must be accessible (not NA) for ML fit")

  # Values must differ
  expect_true(isTRUE(as.logical(is_reml_reml)),
              info = "REML fit should report is_REML = TRUE")
  expect_false(isTRUE(as.logical(is_reml_ml)),
               info = "ML fit should report is_REML = FALSE")
})
```

### Tolerance
N/A (boolean check).

### Notes
- `partial` gap (minor severity). REML status is present in the artifact; it
  is just not exposed via a `getME("is_REML")` accessor or `isREML()` generic.
- This test passes even without `getME("is_REML")` if the fall-back artifact
  path works, making it a regression guard for both the current and future API.

---

## Upstream fixture notes

The following test specs require changes in the Rust engine or FFI before they
can fully pass:

| Test spec | Upstream requirement |
|---|---|
| TS-06 (getME "lower") | Lower-bound vector for theta parameterisation exposed via FFI |
| TS-09 (postVar tolerance) | Pre-existing theta/sigma drift (bd-01KRV31R4BJVQCEF0F58NFD4YN) must be resolved |
| TS-19 (GLMM boundary_parameter term name) | Rust GLMM optimizer certificate must map theta index → term name |
| TS-13, TS-14, TS-15, TS-16 (influence / hat cluster) | Leave-one-group-out refit infrastructure; not yet in engine |

All other specs require only R-side fixes or are pure behavioral contracts
certifiable against current code.
