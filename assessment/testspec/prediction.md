# Test Specification — Prediction & Residuals

**Family:** Prediction & residuals (`predict`, `fitted`, `residuals`, `simulate`, `bootMer`)
**Date:** 2026-05-31
**Source gap report:** `assessment/gap/prediction.md`
**Reference implementations:** lme4 2.0.1, lmerTest 3.2.1
**Passing these specs certifies:** mixeff achieves lme4 parity (or an honest documented
divergence) for every in-scope-missing, partial, or test-gap item in this family.

---

## Priority legend

| Code | Meaning |
|------|---------|
| P0   | Blocker — silent wrong answer or primary use-case failure; must pass before any release. |
| P1   | Major capability gap that a typical lme4 user will hit; should pass before GA. |
| P2   | Minor / niche gap; good to have but acceptable as deferred. |

---

## Spec 1 — `predict(random.only=TRUE)` must raise a typed error, not return the full fit

**Name:** `predict_random_only_raises_not_silently_wrong`
**Kind:** error-message
**Priority:** P0

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))
```

### Assertion

```r
test_that("predict(random.only=TRUE) raises mm_arg_error, not silent wrong answer", {
  expect_error(
    stats::predict(fit, random.only = TRUE),
    class = "mm_arg_error"
  )
})
```

The error message must mention `random.only` by name so the user knows why it was
rejected. The current silent-wrong-answer behaviour (returning full conditional fitted
values identical to `predict(fit)`) is the defect being guarded.

### What it guards against

`predict(random.only = TRUE)` returning the full conditional prediction with no
diagnostic. This is the single most serious finding in the gap report: a documented
lme4 argument that silently returns a plausible but wrong answer.

### Upstream fixture / engine change needed?

None. The fix is purely in `predict.mm_lmm`: detect `random.only` in `...` and
dispatch to `mm_abort(..., class = "mm_arg_error")`.

---

## Spec 2 — `residuals(scaled=TRUE)` must raise a typed error, not return unscaled residuals

**Name:** `residuals_scaled_TRUE_raises_not_silently_wrong`
**Kind:** error-message
**Priority:** P0

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))
```

### Assertion

```r
test_that("residuals(scaled=TRUE) raises mm_arg_error, not unscaled residuals", {
  expect_error(
    residuals(fit, scaled = TRUE),
    class = "mm_arg_error"
  )
  # Guard: current (broken) path returns identical to residuals(fit, scaled=FALSE)
  raw <- residuals(fit)
  expect_false(
    tryCatch({
      r <- residuals(fit, scaled = TRUE)
      isTRUE(all.equal(unname(r), unname(raw)))
    }, error = function(e) FALSE),
    info = "scaled=TRUE must not silently return unscaled residuals"
  )
})
```

Alternatively, if `scaled=TRUE` is implemented (Spec 7), this spec becomes a
numerical-tolerance test rather than an error-message test; the two specs are
mutually exclusive depending on implementation choice.

### What it guards against

`residuals(scaled = TRUE)` returning raw residuals with no diagnostic — the same
class of silent-wrong-answer violation as Spec 1.

### Upstream fixture / engine change needed?

None for the error path. For the implemented path (Spec 7) no Rust work is needed:
`scaled = response / sigma`.

---

## Spec 3 — `simulate(newdata=...)` must raise a typed error, not ignore `newdata`

**Name:** `simulate_newdata_raises_not_silently_ignored`
**Kind:** error-message
**Priority:** P0

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))
nd <- sleepstudy[1:5, ]
```

### Assertion

```r
test_that("simulate(newdata=) raises mm_arg_error, not returns 180 training rows", {
  expect_error(
    simulate(fit, newdata = nd),
    class = "mm_arg_error"
  )
  # Guard: current (broken) path returns 180 rows ignoring newdata
  result <- tryCatch(simulate(fit, newdata = nd), error = function(e) e)
  if (!inherits(result, "error")) {
    expect_equal(nrow(result), 5L,
                 info = "if newdata is accepted, row count must match newdata")
  }
})
```

### What it guards against

`simulate(newdata = nd)` silently simulating from the 180-row training set instead of
the 5 supplied rows — a silent wrong-answer violation of the audit-first contract.

### Upstream fixture / engine change needed?

None for the error path. Full newdata-simulate support would require Rust contract
work (out of scope for this spec).

---

## Spec 4 — `predict(newdata=)` with NA in the response column must use level policy, not NaN guard

**Name:** `predict_newdata_na_response_uses_level_policy`
**Kind:** error-message + parity-vs-lme4
**Priority:** P0

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
train <- sleepstudy[sleepstudy$Subject != "372", ]
train$Subject <- droplevels(train$Subject)
# newdata with NA response — the canonical lme4 idiom for prediction
nd <- data.frame(Days = c(0, 5, 9), Subject = factor("372"), Reaction = NA_real_)
fit <- lmm(Reaction ~ Days + (Days | Subject), data = train, REML = TRUE,
           control = mm_control(verbose = -1))
ref <- suppressMessages(suppressWarnings(
  lme4::lmer(Reaction ~ Days + (Days | Subject), data = train, REML = TRUE)
))
```

### Assertions

```r
test_that("predict(allow.new.levels=FALSE) raises mm_inference_unavailable for NA-response newdata", {
  expect_error(
    stats::predict(fit, newdata = nd, re.form = NULL, allow.new.levels = FALSE),
    class = "mm_inference_unavailable"
  )
  # Must NOT raise mm_data_error about NaN in Reaction
  err <- tryCatch(
    stats::predict(fit, newdata = nd, re.form = NULL, allow.new.levels = FALSE),
    error = function(e) e
  )
  expect_false(inherits(err, "mm_data_error"),
               info = "NaN guard must not fire before the new-levels policy check")
})

test_that("predict(allow.new.levels=TRUE) returns population predictions for NA-response newdata", {
  obs <- stats::predict(fit, newdata = nd, re.form = NULL, allow.new.levels = TRUE)
  ref_pred <- stats::predict(ref, newdata = nd, re.form = NULL, allow.new.levels = TRUE)
  expect_true(all(is.finite(obs)))
  expect_equal(length(obs), 3L)
  expect_equal(unname(obs), unname(ref_pred), tolerance = 1e-4)
})
```

### What it guards against

The Rust bridge applying its NaN/Inf guard to the response column before the
new-levels policy check fires, causing a `mm_data_error` ("NaN in Reaction") when
the user legitimately passes `NA` in the response column of newdata — the canonical
lme4 idiom. Documented in `assessment/parity/inf-predict-newdata.md` Finding B.

### Upstream fixture / engine change needed?

The fix is in `mm_predict_conditional_newdata()` (`R/predict.R:188`): strip or replace
the response column before passing newdata to `mm_translate_data`. No Rust changes
required.

---

## Spec 5 — `residuals(type="pearson")` matches lme4 (response / sigma)

**Name:** `residuals_pearson_parity`
**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))
ref <- suppressMessages(suppressWarnings(
  lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
))
```

### Assertion

```r
test_that("residuals(type='pearson') matches lme4 within tolerance", {
  obs <- residuals(fit, type = "pearson")
  exp <- residuals(ref, type = "pearson")
  expect_equal(unname(obs), unname(exp), tolerance = 1e-2)
  # For Gaussian LMM pearson = response / sigma; verify algebraic identity
  expect_equal(unname(obs), unname(residuals(fit)) / sigma(fit), tolerance = 1e-10)
})
```

Tolerance 1e-2 is consistent with the conditional-fitted tolerance (driven by the
known sigma/theta optimizer gap). The algebraic identity `pearson = response/sigma`
must hold to 1e-10 (no Rust round-trip involved).

### What it guards against

`match.arg` hard-error on `type="pearson"` — currently the only allowed type is
`"response"`. This is a trivially computable type that lme4 users routinely request.

### Upstream fixture / engine change needed?

None. Pure R: `out <- object$residuals / sigma(object)`.

---

## Spec 6 — `residuals(type="deviance")` matches lme4 for LMM (equals response)

**Name:** `residuals_deviance_parity_lmm`
**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Dataset / formula

Same as Spec 5.

### Assertion

```r
test_that("residuals(type='deviance') matches lme4 for LMM (deviance == response for Gaussian)", {
  obs <- residuals(fit, type = "deviance")
  exp <- residuals(ref, type = "deviance")
  expect_equal(unname(obs), unname(exp), tolerance = 1e-2)
  # For Gaussian LMM deviance residuals equal response residuals
  expect_equal(unname(obs), unname(residuals(fit, type = "response")),
               tolerance = 1e-10)
})
```

### What it guards against

`match.arg` hard-error on `type="deviance"`. For LMM the answer is identical to the
response residuals and requires no new computation.

### Upstream fixture / engine change needed?

None. Pure R alias.

---

## Spec 7 — `residuals(scaled=TRUE)` matches lme4 (response / sigma)

**Name:** `residuals_scaled_TRUE_parity`
**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Dataset / formula

Same as Spec 5.

### Assertion

```r
test_that("residuals(scaled=TRUE) matches lme4 scaled residuals", {
  obs <- residuals(fit, scaled = TRUE)
  exp <- residuals(ref, scaled = TRUE)
  # lme4 scaled = response / sigma; verify parity
  expect_equal(unname(obs), unname(exp), tolerance = 1e-2)
  # Algebraic identity: must equal response / sigma to machine precision
  expect_equal(unname(obs), unname(residuals(fit)) / sigma(fit), tolerance = 1e-10)
})
```

Note: This spec supersedes Spec 2 if the feature is implemented. If the decision is
to raise an error, retain Spec 2 and mark this spec as deferred.

### What it guards against

Silent return of unscaled residuals when `scaled=TRUE` is passed — the identical
silent-wrong-answer pattern as the `random.only=TRUE` defect.

### Upstream fixture / engine change needed?

None. Pure R: `out <- object$residuals / sigma(object)`.

---

## Spec 8 — `predict(newdata=, na.action=na.pass)` propagates NA rows

**Name:** `predict_newdata_na_action_na_pass`
**Kind:** parity-vs-lme4
**Priority:** P1

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))
ref <- suppressMessages(suppressWarnings(
  lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
))
# Three rows: rows 1 and 3 are valid, row 2 has NA in a predictor
nd <- sleepstudy[c(1, 2, 3), ]
nd$Days[2L] <- NA_real_
```

### Assertion

```r
test_that("predict(newdata, na.action=na.pass) returns NA for NA-predictor rows", {
  obs <- stats::predict(fit, newdata = nd, na.action = stats::na.pass)
  exp <- stats::predict(ref, newdata = nd, na.action = stats::na.pass)
  # lme4 returns c(pred1, NA, pred3) — length 3, NA in position 2
  expect_equal(length(obs), 3L)
  expect_true(is.na(obs[[2L]]))
  expect_true(is.finite(obs[[1L]]))
  expect_true(is.finite(obs[[3L]]))
  # Finite values must match lme4
  expect_equal(obs[c(1, 3)], exp[c(1, 3)], tolerance = 1e-2, ignore_attr = TRUE)
})
```

Currently mixeff hard-errors with `mm_data_error: numeric column 'Days' contains a
non-finite value`. lme4 returns `NA` for that row and preserves row count (na.pass
semantics). This breaks common cross-validation and grid-prediction workflows.

### What it guards against

Silent hard error on a standard lme4 workflow: passing a data frame with NA predictor
values to `predict` when `na.action = na.pass` is requested.

### Upstream fixture / engine change needed?

The population path (`re.form=NA`) in `mm_predict_fixed_only` already calls
`stats::model.frame(..., na.action = stats::na.pass)`, so the fix there is to forward
the NA rows through and emit `NA` predictions for them. The conditional path
(`mm_predict_conditional_newdata`) would need similar treatment before handing off to
the Rust FFI, or by inserting `NA` for rows that were filtered out.

---

## Spec 9 — `predict(se.fit=TRUE)` population path returns finite SEs from `vcov`

**Name:** `predict_se_fit_population_path_finite`
**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))
ref <- suppressMessages(suppressWarnings(
  lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
))
nd <- sleepstudy[1:5, ]
```

### Assertion

```r
test_that("predict(se.fit=TRUE, re.form=NA) returns finite population SEs matching lme4", {
  obs <- stats::predict(fit, newdata = nd, re.form = NA, se.fit = TRUE)
  ref_pred <- stats::predict(ref, newdata = nd, re.form = NA, se.fit = TRUE)

  expect_true(all(is.finite(obs$se.fit)),
              info = "population-path SEs must be finite; vcov(fit) is available")
  expect_equal(unname(obs$fit), unname(ref_pred$fit), tolerance = 1e-6)
  expect_equal(unname(obs$se.fit), unname(ref_pred$se.fit), tolerance = 1e-3)
})

test_that("predict(se.fit=TRUE, re.form=NULL) in-sample still returns NA SEs with reason attr", {
  obs <- stats::predict(fit, se.fit = TRUE)
  expect_true(all(is.na(obs$se.fit)))
  expect_identical(attr(obs, "mm_unavailable_reason"),
                   "prediction_se_unavailable_phase_2")
})
```

The population-path SE is `sqrt(diag(X_new %*% vcov(fit) %*% t(X_new)))` — exact for
LMM. `vcov(fit)` already returns the full 2×2 fixed-effect covariance. This is
computable entirely on the R side; the conditional SE (which requires the joint
FE+RE/prediction variance from Rust) remains `NA` and is out-of-scope-by-design.

### What it guards against

Returning `NA` for population-path SEs when `vcov(fit)` is already available on the R
side — a more conservative stance than necessary that breaks downstream SE-based
workflows (e.g., manual confidence bands, `emmeans` SE propagation).

### Upstream fixture / engine change needed?

None for the population path. No Rust changes required.

---

## Spec 10 — `simulate()` output structure matches lme4 exactly

**Name:** `simulate_output_structure_parity`
**Kind:** snapshot + parity-vs-lme4
**Priority:** P1

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))
ref <- suppressMessages(suppressWarnings(
  lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
))
```

### Assertion

```r
test_that("simulate() output structure matches lme4", {
  obs <- simulate(fit, nsim = 3, seed = 42)
  exp <- simulate(ref, nsim = 3, seed = 42)

  expect_s3_class(obs, "data.frame")
  expect_equal(dim(obs), c(180L, 3L))
  expect_equal(colnames(obs), paste0("sim_", 1:3))
  expect_equal(rownames(obs), rownames(sleepstudy))
  expect_false(is.null(attr(obs, "seed")))

  # Grand mean must be within 3 SE of lme4 grand mean at nsim=3
  # (structure check; moments are verified in Spec 11)
  expect_true(all(sapply(obs, is.numeric)))
  expect_false(anyNA(obs))
})

test_that("simulate(seed=X) is reproducible", {
  s1 <- simulate(fit, nsim = 2, seed = 42)
  s2 <- simulate(fit, nsim = 2, seed = 42)
  expect_identical(s1, s2)
})

test_that("simulate() different seeds produce different draws", {
  s1 <- simulate(fit, nsim = 1, seed = 99)
  s2 <- simulate(fit, nsim = 1, seed = 100)
  expect_false(isTRUE(all.equal(s1, s2)))
})
```

This codifies the structure parity already confirmed in `assessment/parity/inf-simulate.md`.

### What it guards against

Regressions to the simulate output contract (class, dims, colnames, rownames, seed
attribute, no NA).

### Upstream fixture / engine change needed?

None.

---

## Spec 11 — `simulate()` moment parity with lme4 at nsim=200

**Name:** `simulate_moments_parity`
**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Dataset / formula

Same as Spec 10.

### Assertion

```r
test_that("simulate() grand mean matches lme4 within 3 Monte Carlo SEs at nsim=200", {
  set.seed(7)
  obs <- simulate(fit, nsim = 200, seed = 7)
  exp <- simulate(ref, nsim = 200, seed = 7)

  obs_mean <- mean(unlist(obs))
  exp_mean <- mean(unlist(exp))
  fixef_mean <- fixef(fit)[["(Intercept)"]]

  # SE of grand mean ~ sigma_total / sqrt(N * nsim)
  # For sleepstudy: sigma_total ~ 50, N=180, nsim=200; SE ~ 50/sqrt(36000) ~ 0.26
  mc_se <- 50 / sqrt(180 * 200)
  expect_lt(abs(obs_mean - fixef_mean), 3 * mc_se,
            label = "mixeff grand mean within 3 MC SEs of fixef intercept")
  expect_lt(abs(obs_mean - exp_mean), 6 * mc_se,
            label = "mixeff and lme4 grand means within 6 MC SEs of each other")
})
```

Tolerance is Monte Carlo, not numeric. The 200-sim threshold keeps the test runtime
under ~2 s while providing ~0.26 SE resolution on the grand mean.

### What it guards against

Systematic bias in the simulate path — e.g., using a wrong sigma or wrong VarCorr
entries that shifts the grand mean by more than sampling noise.

### Upstream fixture / engine change needed?

None.

---

## Spec 12 — `simulate()→refit()` round-trip preserves model structure

**Name:** `simulate_refit_roundtrip`
**Kind:** integration
**Priority:** P1

### Dataset / formula

Same as Spec 10.

### Assertion

```r
test_that("simulate()->refit() round-trip returns a valid mm_lmm with plausible fixef", {
  sim <- simulate(fit, nsim = 1, seed = 1L)
  refitted <- refit(fit, newresp = sim[[1L]])
  expect_s3_class(refitted, "mm_lmm")
  expect_equal(length(fixef(refitted)), length(fixef(fit)))
  expect_true(all(is.finite(fixef(refitted))))
  expect_true(is.finite(sigma(refitted)))
  # fixef should differ from original (new random data) but be in the same ballpark
  expect_gt(max(abs(fixef(refitted) - fixef(fit))), 1e-6)
})
```

### What it guards against

Breakage of the `refit()` building block used by any future `bootMer` implementation.

### Upstream fixture / engine change needed?

None.

---

## Spec 13 — `bootMer(type="parametric")` equivalent exists and returns valid CIs

**Name:** `bootmer_parametric_basic`
**Kind:** integration
**Priority:** P1

### Dataset / formula

```r
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))
```

### Assertion

```r
test_that("bootMer-equivalent returns parametric bootstrap CIs for fixef", {
  # PRD §10 Phase 4 lists bootMer(type='parametric') as in-scope.
  # This test is provisional: it should be skipped until the function exists.
  skip_if_not(exists("bootMer", asNamespace("mixeff")),
              "bootMer not yet in mixeff namespace")

  boot_out <- bootMer(fit, FUN = fixef, nsim = 50, type = "parametric",
                      seed = 42)
  expect_true(is.list(boot_out) || inherits(boot_out, "boot"))
  if (is.list(boot_out)) {
    # mixeff may return a list with $t (nsim x p matrix)
    expect_equal(ncol(boot_out$t), length(fixef(fit)))
    expect_equal(nrow(boot_out$t), 50L)
  }
  # All bootstrap replicates for intercept should be finite
  t_int <- if (is.list(boot_out)) boot_out$t[, 1L] else boot_out$t[, 1L]
  expect_true(all(is.finite(t_int)))
  # Intercept boot distribution should center near fit$beta[1] within 10 units
  expect_lt(abs(mean(t_int) - fixef(fit)[[1L]]), 10)
})
```

### What it guards against

The Phase 4 `bootMer` deliverable shipping without a basic smoke test. The `skip_if`
guard means the test is a latent spec until the function is implemented; it will
auto-activate when the namespace export appears.

### Upstream fixture / engine change needed?

None (uses existing `simulate()` + `refit()` building blocks).

---

## Spec 14 — `simulate(use.u=TRUE)` raises a typed error, not silent drop

**Name:** `simulate_use_u_raises_not_silently_ignored`
**Kind:** error-message
**Priority:** P2

### Dataset / formula

Same as Spec 10.

### Assertion

```r
test_that("simulate(use.u=TRUE) raises mm_arg_error, not silently ignores the arg", {
  expect_error(
    simulate(fit, use.u = TRUE),
    class = "mm_arg_error"
  )
})
```

`use.u` was deprecated in lme4 in favour of `re.form`; it should be refused loudly so
users know to switch to `re.form = NULL`.

### What it guards against

Silent drop of `use.u` via `...`, which would simulate from the wrong distribution
without warning.

### Upstream fixture / engine change needed?

None. R-side detection of `use.u` in `...`.

---

## Spec 15 — `predict(newparams=)` raises a typed error, not silent drop

**Name:** `predict_newparams_raises_not_silently_ignored`
**Kind:** error-message
**Priority:** P2

### Dataset / formula

Same as Spec 1.

### Assertion

```r
test_that("predict(newparams=list(...)) raises mm_arg_error, not silently ignores", {
  expect_error(
    stats::predict(fit, newparams = list(beta = fixef(fit), theta = c(1, 0, 1))),
    class = "mm_arg_error"
  )
})
```

### What it guards against

Counterfactual predictions returning the fitted model's predictions with no diagnostic
when the user intended a counterfactual.

### Upstream fixture / engine change needed?

None.

---

## Spec 16 — `simulate(newparams=)` raises a typed error, not silent drop

**Name:** `simulate_newparams_raises_not_silently_ignored`
**Kind:** error-message
**Priority:** P2

### Dataset / formula

Same as Spec 10.

### Assertion

```r
test_that("simulate(newparams=list(...)) raises mm_arg_error, not silently ignores", {
  expect_error(
    simulate(fit, newparams = list(beta = fixef(fit), sigma = 1)),
    class = "mm_arg_error"
  )
})
```

### What it guards against

A simulate call under modified parameters returning draws from the original
parameters with no diagnostic.

### Upstream fixture / engine change needed?

None.

---

## Spec 17 — `predict()` conditional parity on Pastes / intercept-only random slope

**Name:** `predict_conditional_pastes_parity`
**Kind:** parity-vs-lme4 + numerical-tolerance
**Priority:** P1

### Dataset / formula

```r
data("Pastes", package = "lme4")
# A model with two nested grouping factors
fit <- lmm(strength ~ 1 + (1 | batch/cask), data = Pastes, REML = TRUE,
           control = mm_control(verbose = -1))
ref <- suppressMessages(suppressWarnings(
  lme4::lmer(strength ~ 1 + (1 | batch/cask), data = Pastes, REML = TRUE)
))
```

### Assertion

```r
test_that("predict() conditional matches lme4 on Pastes (nested grouping)", {
  obs <- stats::predict(fit)
  exp <- stats::predict(ref)
  expect_equal(unname(obs), unname(exp), tolerance = 1e-2)
})

test_that("predict(re.form=NA) population matches lme4 on Pastes", {
  obs <- stats::predict(fit, re.form = NA)
  exp <- stats::predict(ref, re.form = NA)
  expect_equal(unname(obs), unname(exp), tolerance = 1e-6)
})
```

Tolerance 1e-2 allows for the known optimizer theta gap; 1e-6 for fixed-only path.

### What it guards against

Prediction regressions on nested/multi-grouping models beyond the single-group
`sleepstudy` baseline.

### Upstream fixture / engine change needed?

None beyond confirming the Pastes model fits.

---

## Spec 18 — `residuals(type="working")` matches lme4 for LMM (equals response)

**Name:** `residuals_working_parity_lmm`
**Kind:** parity-vs-lme4
**Priority:** P2

### Dataset / formula

Same as Spec 5.

### Assertion

```r
test_that("residuals(type='working') matches lme4 for LMM (working == response for Gaussian)", {
  obs <- residuals(fit, type = "working")
  exp <- residuals(ref, type = "working")
  expect_equal(unname(obs), unname(exp), tolerance = 1e-2)
  # For Gaussian LMM working == response
  expect_equal(unname(obs), unname(residuals(fit, type = "response")),
               tolerance = 1e-10)
})
```

### What it guards against

`match.arg` hard-error for `type="working"` (currently not in the allowed set).

### Upstream fixture / engine change needed?

None.

---

## Spec 19 — `fitted()` / `predict()` agree on in-sample data

**Name:** `fitted_predict_consistency`
**Kind:** unit
**Priority:** P1

### Dataset / formula

Same as Spec 5.

### Assertion

```r
test_that("fitted() and predict() are identical for in-sample data", {
  expect_equal(unname(fitted(fit)), unname(stats::predict(fit)),
               tolerance = 1e-12)
  expect_equal(unname(fitted(fit)), unname(stats::predict(fit, newdata = sleepstudy)),
               tolerance = 1e-2)
})
```

This codifies the existing `mm_expect_prediction_lme4_parity` helper assertion and
guards against any future divergence between the cached `$fitted` field and the
predict dispatch.

### What it guards against

Regression where the `$fitted` cache and the live predict path return different values
(e.g., if `newdata=training` routes through the Rust FFI while in-sample uses the
cache).

### Upstream fixture / engine change needed?

None.

---

## Summary

| # | Name | Kind | Priority | Gap classification | Upstream needed? |
|---|------|------|----------|--------------------|-----------------|
| 1 | `predict_random_only_raises_not_silently_wrong` | error-message | P0 | in-scope-missing (silent wrong) | No |
| 2 | `residuals_scaled_TRUE_raises_not_silently_wrong` | error-message | P0 | in-scope-missing (silent wrong) | No |
| 3 | `simulate_newdata_raises_not_silently_ignored` | error-message | P0 | in-scope-missing (silent wrong) | No |
| 4 | `predict_newdata_na_response_uses_level_policy` | error-message + parity | P0 | in-scope-missing (wrong error class) | No (R-side fix) |
| 5 | `residuals_pearson_parity` | parity + numerical-tolerance | P1 | in-scope-missing | No |
| 6 | `residuals_deviance_parity_lmm` | parity + numerical-tolerance | P1 | in-scope-missing | No |
| 7 | `residuals_scaled_TRUE_parity` | parity + numerical-tolerance | P1 | in-scope-missing | No |
| 8 | `predict_newdata_na_action_na_pass` | parity | P1 | in-scope-missing | R-side only |
| 9 | `predict_se_fit_population_path_finite` | parity + numerical-tolerance | P1 | partial | No |
| 10 | `simulate_output_structure_parity` | snapshot + parity | P1 | test-gap (assurance) | No |
| 11 | `simulate_moments_parity` | parity + numerical-tolerance | P1 | test-gap (assurance) | No |
| 12 | `simulate_refit_roundtrip` | integration | P1 | test-gap (assurance) | No |
| 13 | `bootmer_parametric_basic` | integration | P1 | in-scope-missing (Phase 4) | No |
| 14 | `simulate_use_u_raises_not_silently_ignored` | error-message | P2 | in-scope-missing | No |
| 15 | `predict_newparams_raises_not_silently_ignored` | error-message | P2 | in-scope-missing | No |
| 16 | `simulate_newparams_raises_not_silently_ignored` | error-message | P2 | in-scope-missing | No |
| 17 | `predict_conditional_pastes_parity` | parity + numerical-tolerance | P1 | test-gap | No |
| 18 | `residuals_working_parity_lmm` | parity | P2 | in-scope-missing | No |
| 19 | `fitted_predict_consistency` | unit | P1 | test-gap (assurance) | No |

**Total specs: 19**

**Top priority (P0 blocker):** `predict_random_only_raises_not_silently_wrong` (Spec 1) —
a documented lme4 argument silently returns a numerically plausible but wrong answer
with no diagnostic, directly violating the project's "No silent surgery" contract
(`CLAUDE.md`) and the "clearer errors" mandate (PRD §1).

---

## Implementation notes

### All three P0 silent-wrong-answer bugs share a single root cause

`predict.mm_lmm`, `residuals.mm_lmm`, and `simulate.mm_lmm` each absorb unrecognized
arguments via `...`. The minimal fix for all three P0 specs is to add early detection
of the specific lme4 argument names before `...` is silently discarded:

```r
# In predict.mm_lmm before the target dispatch:
mm_intercept_unsupported_dots(
  ...,
  unsupported = c("random.only", "newparams"),
  context = "predict"
)

# In residuals.mm_lmm:
mm_intercept_unsupported_dots(
  ...,
  unsupported = c("scaled", "newparams"),
  context = "residuals"
)

# In simulate.mm_lmm:
mm_intercept_unsupported_dots(
  ...,
  unsupported = c("newdata", "newparams", "use.u"),
  context = "simulate"
)
```

A single `mm_intercept_unsupported_dots()` helper can implement this pattern once and
be applied to all three methods, covering Specs 1, 2, 3, 14, 15, 16 simultaneously.

### Specs 5, 6, 7, 18 require a single one-line change in `residuals.mm_lmm`

Extend `match.arg` to `type = c("response", "pearson", "deviance", "working")` and
add the trivial R-side computation for each type. No Rust work needed.

### Spec 9 requires adding a `se.fit` handler to the population path only

In `mm_predict_fixed_only()`: after computing `pred <- as.numeric(mm_new %*% beta)`,
add `se <- sqrt(diag(mm_new %*% vcov(fit) %*% t(mm_new)))`. Gate on
`isTRUE(se.fit)` and `identical(target, "population")`.

### Spec 13 is a latent spec (auto-activates when Phase 4 ships `bootMer`)

The `skip_if_not(exists("bootMer", ...))` guard means the test file can be committed
immediately; it becomes active automatically when the function is exported.
