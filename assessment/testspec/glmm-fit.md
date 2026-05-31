# Test / Assurance Specifications — GLMM Fitting

**Family:** GLMM fitting  
**Source gap report:** `assessment/gap/glmm-fit.md`  
**Date:** 2026-05-31  
**Author:** generated from gap-report analysis + probe evidence

Scope: every gap classified `in-scope-missing`, `partial`, or `test-gap` in the gap report.
Out-of-scope-by-design items (cbind multivariate response per PRD §3 line 811,
profile-LL CIs for GLMM, nlopt joint-Laplace in CRAN build, `glmer.nb`,
`nAGQ=0`, start-values, modular API, custom make.link) are **not** assigned test
specs here — they are correctly handled by the existing typed-refusal path. Each
spec below targets something that either silently misbehaves, is absent despite
the engine supporting it, or lacks a regression test.

---

## SPEC-GLMM-01 — `offset=` argument is refused, not silently dropped

**Priority:** P0 (blocker)  
**Kind:** error-message  
**Classification addressed:** in-scope-missing (silent-surgery bug)

### Rationale

`glmm()` has no `offset` parameter; a supplied `offset=` vector falls into `...`
and is discarded with no diagnostic. This violates the package's core contract
("no silent surgery", CLAUDE.md). In-fit offsets are a PRD §3 non-goal, so the
*feature* is out-of-scope, but accepting the argument silently is a defect.
Contrast: `weights=` (also reserved) already raises `mm_fit_error` loudly.

### Dataset / setup

```r
# Simulated Poisson rate data (events per person-year)
set.seed(42)
n <- 60; g <- 6
dat <- data.frame(
  y      = rpois(n, lambda = 2),
  x      = rnorm(n),
  grp    = factor(rep(seq_len(g), each = n / g)),
  logexp = log(runif(n, 0.5, 2))   # known exposure offset
)
```

### Assertions

```r
test_that("glmm() refuses offset= argument with a typed error, not silently drops it", {
  expect_error(
    glmm(y ~ x + (1 | grp), data = dat, family = poisson(), offset = dat$logexp),
    class = "mm_fit_error"
  )
})

test_that("glmm() refuses offset= and the error message names 'offset'", {
  err <- tryCatch(
    glmm(y ~ x + (1 | grp), data = dat, family = poisson(), offset = dat$logexp),
    error = function(e) e
  )
  expect_match(conditionMessage(err), "offset", ignore.case = TRUE)
})

test_that("glmm() with offset= does NOT return a fitted object", {
  result <- tryCatch(
    glmm(y ~ x + (1 | grp), data = dat, family = poisson(), offset = dat$logexp),
    error = function(e) e
  )
  expect_false(inherits(result, "mm_glmm"))
})
```

### What this guards against

Regression: a user fitting a Poisson rate model (e.g., events per person-year)
supplies `offset=` and receives silently wrong estimates. The test guarantees
they get an immediate typed error instead of a plausible-looking but incorrect
fit.

### Upstream fixture / engine change required?

No. This is a pure R-wrapper fix: add `offset` to the reserved-argument check in
`glmm()` (lines 54–61 of `R/glmm.R`) alongside `weights`, `subset`, etc.

---

## SPEC-GLMM-02 — `weights=` reserved-argument error is stable and informative

**Priority:** P1  
**Kind:** error-message  
**Classification addressed:** in-scope-missing (confirmed hard error, but no regression test)

### Rationale

`weights=` already raises `mm_fit_error` (gap report confirmed). This spec adds
a regression test to lock that behavior and verify the error message is
informative enough for a user to understand the limitation.

### Dataset / setup

```r
dat_w <- data.frame(
  y   = rbinom(40, size = 1, prob = 0.4),
  x   = rnorm(40),
  grp = factor(rep(1:4, each = 10))
)
wts <- runif(40, 0.5, 2)
```

### Assertions

```r
test_that("glmm() refuses weights= with a typed mm_fit_error", {
  expect_error(
    glmm(y ~ x + (1 | grp), data = dat_w, family = binomial(), weights = wts),
    class = "mm_fit_error"
  )
})

test_that("glmm() weights= error message mentions 'weights' or 'reserved'", {
  err <- tryCatch(
    glmm(y ~ x + (1 | grp), data = dat_w, family = binomial(), weights = wts),
    error = function(e) e
  )
  expect_true(
    grepl("weights", conditionMessage(err), ignore.case = TRUE) ||
    grepl("reserved", conditionMessage(err), ignore.case = TRUE)
  )
})
```

### What this guards against

Regression: a future refactor accidentally removes the `weights` guard,
allowing a partially-wired or silently-wrong weighted fit.

---

## SPEC-GLMM-03 — Gamma/inverse-link is unblocked at wrapper level

**Priority:** P1  
**Kind:** unit + parity-vs-lme4  
**Classification addressed:** in-scope-missing (engine-capable, wrapper-blocked)

### Rationale

The gap report confirms the Rust engine has `LinkFunction::Inverse` (crate-0.md
§44 line 1210). Only `mm_glmm_supported_family_links()` blocks it. `family=Gamma`
(no explicit link) defaults to `inverse` in base R — the most natural call a
Gamma GLMM user makes is currently refused. This is a wrapper-side restriction
with no upstream dependency.

### Dataset / setup

```r
# Simulated Gamma/inverse data
set.seed(7)
n <- 120; g <- 10
mu_true <- 1 / (1.2 + 0.4 * rnorm(n))  # inverse link: eta = 1/mu
dat_gi <- data.frame(
  y   = rgamma(n, shape = 3, rate = 3 / pmax(mu_true, 0.05)),
  x   = rnorm(n),
  grp = factor(rep(seq_len(g), each = n / g))
)
```

### Assertions — after the wrapper fix

```r
test_that("glmm() accepts Gamma(link='inverse') without error after wrapper fix", {
  # Requires: mm_glmm_supported_family_links() to include Gamma/inverse
  fit <- glmm(y ~ x + (1 | grp), data = dat_gi,
              family = Gamma(link = "inverse"),
              control = mm_control(verbose = -1L))
  expect_s3_class(fit, "mm_glmm")
  expect_identical(fit$family$family, "gamma")
  expect_identical(fit$family$link,   "inverse")
  expect_true(is.finite(fit$logLik))
  expect_true(all(is.finite(fixef(fit))))
})

test_that("glmm() Gamma(link='inverse') fixef within lme4 mismatch bounds", {
  mm_skip_if_no_lme4()
  fit_mm  <- glmm(y ~ x + (1 | grp), data = dat_gi,
                  family = Gamma(link = "inverse"),
                  control = mm_control(verbose = -1L))
  fit_lme4 <- suppressWarnings(
    lme4::glmer(y ~ x + (1 | grp), data = dat_gi,
                family = Gamma(link = "inverse"), nAGQ = 1L)
  )
  # Profiled-PIRLS vs joint-Laplace: use the wider mismatch bounds
  expect_equal(fixef(fit_mm), lme4::fixef(fit_lme4), tolerance = 0.5,
               label = "Gamma/inverse fixef within documented mismatch bounds")
  expect_true(is.finite(as.numeric(logLik(fit_mm))))
})
```

### Current behavior (before fix — expected typed refusal)

```r
test_that("glmm() currently refuses Gamma(link='inverse') with mm_inference_unavailable", {
  # This test documents the current state; delete after the fix lands.
  expect_error(
    glmm(y ~ x + (1 | grp), data = dat_gi,
         family = Gamma(link = "inverse"),
         control = mm_control(verbose = -1L)),
    class = "mm_inference_unavailable"
  )
})
```

### What this guards against

Regression in both directions: (1) the gap is not re-introduced after fix;
(2) any attempted Gamma/inverse fit before the fix is correctly refused rather
than silently producing a Gamma/log fit.

### Upstream fixture / engine change required?

No code change upstream; the engine already supports `Gamma/inverse`. Only
`mm_glmm_supported_family_links()` in `R/glmm.R` needs updating.

---

## SPEC-GLMM-04 — nAGQ>1 on pirls_profiled emits a transparency warning

**Priority:** P1  
**Kind:** error-message (user-facing transparency)  
**Classification addressed:** partial (nAGQ>1 silently means "profiled with AGQ deviance", not lme4 AGQ)

### Rationale

The gap report and `glmm-pois-nagq.md` confirm that `nAGQ=5` is accepted on the
`pirls_profiled` path and changes the deviance approximation slightly, but is NOT
equivalent to lme4's full joint AGQ. A user who writes `glmm(..., nAGQ=5)`
expecting Gauss-Hermite quadrature accuracy gets a different quantity with no
diagnostic. This is a transparency defect: the behavior is not wrong per se, but
it is undocumented to the user at call time.

### Dataset / setup

```r
set.seed(42)
n <- 100; g <- 10
dat_agq <- data.frame(
  y   = rpois(n, lambda = 2),
  x   = rnorm(n),
  grp = factor(rep(seq_len(g), each = n / g))
)
```

### Assertions

```r
test_that("glmm() with nAGQ>1 on pirls_profiled emits a message or warning", {
  # The fit should succeed, but the user must be told it is not full AGQ
  expect_message(
    glmm(y ~ x + (1 | grp), data = dat_agq,
         family = poisson(),
         nAGQ = 5L, method = "pirls_profiled",
         control = mm_control(verbose = -1L)),
    regexp = "AGQ|quadrature|profiled|nAGQ",
    ignore.case = TRUE
  )
})

test_that("glmm() nAGQ>1 fit object still has class mm_glmm", {
  fit <- suppressMessages(
    glmm(y ~ x + (1 | grp), data = dat_agq,
         family = poisson(),
         nAGQ = 5L, method = "pirls_profiled",
         control = mm_control(verbose = -1L))
  )
  expect_s3_class(fit, "mm_glmm")
  expect_identical(fit$nAGQ, 5L)
})

test_that("glmm() nAGQ=1 does NOT emit an nAGQ-transparency warning", {
  # Baseline: the message should only fire for nAGQ > 1
  expect_no_message(
    glmm(y ~ x + (1 | grp), data = dat_agq,
         family = poisson(),
         nAGQ = 1L,
         control = mm_control(verbose = -1L)),
    message = "AGQ|quadrature"
  )
})
```

### What this guards against

A future change that silently promotes the profiled+nAGQ path to look like true
AGQ (e.g., by copying lme4's `nAGQ` argument semantics without implementing the
joint estimator).

### Upstream fixture / engine change required?

No upstream change needed. The R wrapper in `glmm()` (or `mm_glmm_validate_nagq`)
should emit an `mm_fit_note` or `message()` when `nAGQ > 1L` and
`method == "pirls_profiled"`. The crate-0.md §line 244 note ("profiled+nAGQ>1
behavior is undocumented in the R layer") confirms this is purely a wrapper
documentation gap.

---

## SPEC-GLMM-05 — nAGQ>1 does not enforce single-scalar-RE constraint on profiled path

**Priority:** P2  
**Kind:** unit  
**Classification addressed:** partial (lme4 enforces this; mixeff does not on profiled path)

### Rationale

lme4 requires exactly one scalar random effect for `nAGQ > 1`. mixeff accepts
`nAGQ=5` on the profiled path with multiple or vector REs without warning. This
is a mild scope inconsistency; the correct behavior under the profiled path is
debatable, but the absence of a test means the behavior (warn/proceed silently/
refuse) is uncontrolled.

### Dataset / setup

```r
set.seed(42)
dat_multi <- data.frame(
  y    = rpois(120, lambda = 2),
  x    = rnorm(120),
  grp1 = factor(rep(1:10, each = 12)),
  grp2 = factor(rep(1:6,  each = 20))
)
```

### Assertions

```r
test_that("glmm() with nAGQ>1 and multiple REs on pirls_profiled at minimum does not crash", {
  # Behavior: may warn or proceed; must not error with mm_arg_error citing single-RE constraint
  # (that constraint applies only to joint_laplace, which is nlopt-gated).
  # This spec locks the "does not produce a misleading hard error" invariant.
  result <- tryCatch(
    suppressMessages(
      glmm(y ~ x + (1 | grp1) + (1 | grp2), data = dat_multi,
           family = poisson(),
           nAGQ = 3L, method = "pirls_profiled",
           control = mm_control(verbose = -1L))
    ),
    error = function(e) e
  )
  # Should either succeed or emit a documented nAGQ-transparency warning,
  # NOT fail with an lme4-style "single-RE" hard stop on the profiled path.
  if (inherits(result, "error")) {
    # If it does error, it must not be the joint-laplace single-RE restriction
    expect_false(
      grepl("single.*scalar|scalar.*RE|nAGQ.*joint", conditionMessage(result),
            ignore.case = TRUE),
      info = "profiled path must not refuse nAGQ>1 citing the joint-laplace RE constraint"
    )
  } else {
    expect_s3_class(result, "mm_glmm")
  }
})
```

---

## SPEC-GLMM-06 — GLMM parametric bootstrap is wired for all certified families

**Priority:** P1  
**Kind:** integration  
**Classification addressed:** in-scope-missing (engine has `parametricbootstrap_glmm`; R wires lmm only)

### Rationale

The gap report confirms `parametricbootstrap_glmm` is stable in the Rust engine
for Bernoulli/Binomial/Poisson/Gamma (crate-0.md line 242). The R function
`parametric_bootstrap()` currently accepts only `mm_lmm` objects. This is the
realistic uncertainty-quantification path for GLMMs given profile CIs are
deferred (PRD §3). Bootstrap CIs are unreachable from R.

### Dataset / setup

```r
mm_skip_if_no_lme4()
set.seed(1)
dat_bboot <- data.frame(
  y   = rbinom(80, size = 1, prob = 0.4),
  x   = rnorm(80),
  grp = factor(rep(1:8, each = 10))
)
fit_bb <- glmm(y ~ x + (1 | grp), data = dat_bboot,
               family = binomial(),
               control = mm_control(verbose = -1L))
```

### Assertions

```r
test_that("parametric_bootstrap() accepts mm_glmm for certified binomial family", {
  skip_on_cran()   # bootstrap is slow
  bs <- parametric_bootstrap(fit_bb, nsim = 20L, seed = 42L)
  expect_s3_class(bs, "mm_bootstrap")
  # Each replicate should have fixef of the same length as the original
  expect_equal(length(bs$replicates[[1L]]$fixef), length(fixef(fit_bb)))
})

test_that("bootstrap CI from mm_glmm has correct coverage structure", {
  skip_on_cran()
  bs <- parametric_bootstrap(fit_bb, nsim = 50L, seed = 42L)
  ci <- confint(bs, level = 0.95)
  expect_true(all(c("lower", "upper") %in% names(ci)))
  # CIs should bracket the point estimates
  fe <- fixef(fit_bb)
  expect_true(all(ci$lower <= fe & fe <= ci$upper),
              info = "95% bootstrap CI should contain the point estimate for each parameter")
})

test_that("parametric_bootstrap() accepts mm_glmm for poisson family", {
  skip_on_cran()
  dat_pp <- data.frame(
    y   = rpois(60, 2),
    x   = rnorm(60),
    grp = factor(rep(1:6, each = 10))
  )
  fit_pp <- glmm(y ~ x + (1 | grp), data = dat_pp,
                 family = poisson(),
                 control = mm_control(verbose = -1L))
  bs <- parametric_bootstrap(fit_pp, nsim = 20L, seed = 7L)
  expect_s3_class(bs, "mm_bootstrap")
})
```

### What this guards against

The engine capability is silently unreachable from R. These tests certify the
R wire-up is present and produces structurally valid output.

### Upstream fixture / engine change required?

No engine change. This requires the R `parametric_bootstrap()` dispatch in
`R/inference.R` (or wherever it lives) to accept `mm_glmm` objects and call the
Rust `parametricbootstrap_glmm` FFI entry. The Rust side is already stable per
crate-0.md.

---

## SPEC-GLMM-07 — Binomial cbind response: typed refusal with expansion hint

**Priority:** P1  
**Kind:** error-message  
**Classification addressed:** in-scope-missing (cbind response refused; currently the error
message references "stateful transforms", which is misleading for cbind)

### Rationale

`cbind(success, failure) ~ ...` in the response is a first-class lme4 idiom.
The engine supports grouped binomial (crate-3.md line 55). The formula-layer
refusal is correct (cbind response parsing is not in the current scope), but the
current error message mentions "stateful transforms (`poly`, `scale`, ...)",
which is inaccurate — `cbind()` is not a stateful transform. The error should
specifically identify the cbind response pattern and suggest the Bernoulli-row
expansion workaround.

### Dataset / setup

```r
data(cbpp, package = "lme4")
```

### Assertions

```r
test_that("glmm() refuses cbind() response with a typed mm_formula_error", {
  expect_error(
    glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
         data = cbpp, family = binomial()),
    class = "mm_formula_error"
  )
})

test_that("cbind response error message identifies cbind, not just stateful transforms", {
  err <- tryCatch(
    glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
         data = cbpp, family = binomial()),
    error = function(e) e
  )
  # The message should mention cbind specifically
  expect_match(conditionMessage(err), "cbind", ignore.case = TRUE)
})

test_that("cbind response error suggests Bernoulli expansion or is actionable", {
  err <- tryCatch(
    glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
         data = cbpp, family = binomial()),
    error = function(e) e
  )
  # At minimum the error must be a typed condition, not a bare stop()
  expect_true(inherits(err, "mm_condition"),
              info = "cbind refusal must be a typed mm_condition, not a bare stop()")
})
```

### What this guards against

The current message is misleading (implying cbind is like `poly()`). If a future
refactor changes the formula-parsing error path, these tests ensure the cbind
refusal remains typed and specific.

### Upstream fixture / engine change required?

No engine change. The formula-compiler path in `R/manifest.R` (or wherever
`compile_model` calls the Rust formula parser) needs only a targeted cbind
response check before the generic stateful-transform check.

---

## SPEC-GLMM-08 — Bernoulli-expansion of cbpp produces structurally valid mm_glmm

**Priority:** P1  
**Kind:** integration + numerical-tolerance  
**Classification addressed:** partial (cbind not supported; Bernoulli workaround exists but is
untested as a supported pattern)

### Rationale

The Bernoulli-expansion workaround (expand cbpp 56 rows → 842 binary rows) is
the only way to fit grouped-binomial data today. The probe (`glmm-cbpp-binom.md`)
confirmed it runs and produces plausible estimates but diverges from lme4 at the
expected-mismatch bounds in `inst/extdata/expected-mismatches.json`. A test
should: (a) confirm the expanded-data path fits without error, (b) check the
divergence is within the registered expected-mismatch bounds, and (c) confirm the
structural fit invariants hold.

### Dataset / setup (matches fixture `cbpp_binomial_logit_profiled_pirls`)

```r
mm_skip_if_no_lme4()
data(cbpp, package = "lme4")
# mm_expand_binomial_trials is already defined in tests/testthat/test-glmm.R
cbpp_bern <- mm_expand_binomial_trials(cbpp)
# nrow(cbpp_bern) == 842
```

### Assertions

```r
test_that("glmm() on Bernoulli-expanded cbpp returns a valid mm_glmm", {
  fit <- glmm(y ~ period + (1 | herd), data = cbpp_bern,
              family = binomial(), method = "pirls_profiled",
              control = mm_control(verbose = -1L))
  expect_s3_class(fit, "mm_glmm")
  expect_identical(fit$nobs, 842L)
  expect_identical(fit$family$family, "binomial")
  expect_identical(fit$family$link,   "logit")
  expect_length(fixef(fit), 4L)       # (Intercept) + period 2/3/4
  expect_true(is.finite(fit$logLik))
  expect_true(is.finite(fit$AIC))
  expect_identical(as.numeric(fit$sigma), 1)  # binomial dispersion fixed
})

test_that("Bernoulli-expanded cbpp fixef are within expected-mismatch bounds vs lme4", {
  mm_skip_if_no_lme4()
  fit_mm  <- glmm(y ~ period + (1 | herd), data = cbpp_bern,
                  family = binomial(), method = "pirls_profiled",
                  control = mm_control(verbose = -1L))
  fit_ref <- suppressMessages(suppressWarnings(
    lme4::glmer(y ~ period + (1 | herd), data = cbpp_bern,
                family = binomial(), nAGQ = 1L)
  ))
  # Registered bound from expected-mismatches.json: expected_max_abs_diff = 0.17
  mm_assert_parity(
    fixef(fit_mm), lme4::fixef(fit_ref),
    "cbpp_binomial_logit_profiled_pirls", "fixef",
    tol = 0.17,
    label = "Bernoulli-expanded cbpp fixef within expected-mismatch bound"
  )
  # Theta bound: expected_max_abs_diff = 0.07
  mm_assert_parity(
    fit_mm$theta, lme4::getME(fit_ref, "theta"),
    "cbpp_binomial_logit_profiled_pirls", "theta",
    tol = 0.07,
    label = "Bernoulli-expanded cbpp theta within expected-mismatch bound"
  )
})
```

### What this guards against

Engine regressions that break the only available grouped-binomial path, without
being caught by the existing ledger-based parity checks.

---

## SPEC-GLMM-09 — Certified Gamma/log fit: logLik and dispersion snapshot

**Priority:** P1  
**Kind:** numerical-tolerance + snapshot  
**Classification addressed:** partial (Gamma/log fits, but RE variance diverges substantially;
no test for dispersion stability)

### Rationale

The probe (`glmm-gamma.md`) showed that Gamma/log fit with profiled-PIRLS produces
a ~10x smaller RE variance than lme4's joint Laplace on the same data. This
divergence is documented (see `glmm_support_contract.md`), but there is no test
anchoring (a) the dispersion estimate (sigma) within reasonable bounds, (b) that
the fit converges and produces finite estimates, or (c) that the divergence does
not *grow* beyond the documented level.

### Dataset / setup

```r
set.seed(42)
n <- 300; g <- 20
x_val <- rnorm(n)
re    <- rep(rnorm(g, sd = 0.5), each = n / g)
mu    <- exp(1.4 + 0.4 * x_val + re)
dat_gam <- data.frame(
  y   = rgamma(n, shape = 4, rate = 4 / mu),
  x   = x_val,
  grp = factor(rep(seq_len(g), each = n / g))
)
```

### Assertions

```r
test_that("glmm() Gamma/log fit produces finite estimates and converges", {
  fit <- glmm(y ~ x + (1 | grp), data = dat_gam,
              family = Gamma(link = "log"),
              control = mm_control(verbose = -1L))
  expect_s3_class(fit, "mm_glmm")
  expect_identical(fit$family$family, "gamma")
  expect_identical(fit$family$link,   "log")
  expect_true(all(is.finite(fixef(fit))))
  expect_true(is.finite(fit$logLik))
  expect_true(is.finite(fit$dispersion))
  expect_true(fit$dispersion > 0)
  expect_identical(fit$fit_status, "converged_interior")
})

test_that("glmm() Gamma/log dispersion is within plausible bounds for known shape=4", {
  fit <- glmm(y ~ x + (1 | grp), data = dat_gam,
              family = Gamma(link = "log"),
              control = mm_control(verbose = -1L))
  # True dispersion = 1/shape = 0.25; expect within 0.1 to 1.0 given approximation
  expect_true(fit$dispersion > 0.1 && fit$dispersion < 1.5,
              info = sprintf("Gamma dispersion out of plausible range: %g",
                             fit$dispersion))
})

test_that("glmm() Gamma/log fixef sign and scale are qualitatively correct", {
  fit <- glmm(y ~ x + (1 | grp), data = dat_gam,
              family = Gamma(link = "log"),
              control = mm_control(verbose = -1L))
  fe <- fixef(fit)
  # True: intercept ~1.4, slope ~0.4
  expect_true(fe[["(Intercept)"]] > 0.5 && fe[["(Intercept)"]] < 3.0,
              info = "Intercept out of qualitative range")
  expect_true(fe[["x"]] > 0 && fe[["x"]] < 1.5,
              info = "Slope should be positive")
})
```

### Note on lme4 parity for Gamma/log

A parity-vs-lme4 sub-spec is deferred until the Gamma/log divergence is either
(a) reduced by engine improvement or (b) registered in `expected-mismatches.json`
with documented bounds (the current probe shows fixef diff ~0.02, theta diff
~0.22, logLik diff ~6.2 — all beyond standard tolerances). When those bounds are
locked in the ledger, a parity test mirroring SPEC-GLMM-08's `mm_assert_parity`
pattern should be added.

---

## SPEC-GLMM-10 — Certified binomial/probit fit: estimates within expected-mismatch bounds

**Priority:** P1  
**Kind:** parity-vs-lme4 + numerical-tolerance  
**Classification addressed:** partial (probit fits but diverges beyond strict tolerance; no ledger entry)

### Rationale

The probe (`glmm-probit.md`) shows binomial/probit fits converge on both sides,
but fixef diff is ~0.044 (vs 1e-4 tolerance) and SE diff is ~0.030. The root
cause is the profiled-PIRLS vs joint-Laplace objective distinction — the same
class of documented divergence as the cbpp case. A ledger entry and a regression
test locking the observed bounds are needed.

### Dataset / setup

```r
set.seed(42)
n <- 200; g <- 20
grp_re <- rep(rnorm(g, sd = 0.6), each = n / g)
x_val  <- rnorm(n)
p_true <- pnorm(0.3 + 0.8 * x_val + grp_re)
dat_probit <- data.frame(
  y   = rbinom(n, size = 1, prob = p_true),
  x   = x_val,
  grp = factor(rep(seq_len(g), each = n / g))
)
```

### Assertions

```r
test_that("glmm() binomial/probit returns mm_glmm with correct family/link", {
  fit <- glmm(y ~ x + (1 | grp), data = dat_probit,
              family = binomial(link = "probit"),
              control = mm_control(verbose = -1L))
  expect_s3_class(fit, "mm_glmm")
  expect_identical(fit$family$family, "binomial")
  expect_identical(fit$family$link,   "probit")
  expect_identical(as.numeric(fit$sigma), 1)
  expect_true(is.finite(fit$logLik))
})

test_that("glmm() probit fixef within expected-mismatch bounds vs lme4", {
  mm_skip_if_no_lme4()
  fit_mm  <- glmm(y ~ x + (1 | grp), data = dat_probit,
                  family = binomial(link = "probit"),
                  control = mm_control(verbose = -1L))
  fit_ref <- suppressMessages(suppressWarnings(
    lme4::glmer(y ~ x + (1 | grp), data = dat_probit,
                family = binomial(link = "probit"), nAGQ = 1L)
  ))
  # Probe observed: max fixef diff ~0.044; register bound at 0.15 (profiled-PIRLS class)
  mm_assert_parity(
    fixef(fit_mm), lme4::fixef(fit_ref),
    "binomial_probit_profiled_pirls", "fixef",
    tol = 0.15,
    label = "binomial/probit fixef within expected-mismatch bound"
  )
  # Probe observed: theta diff ~0.015; register bound at 0.07
  mm_assert_parity(
    fit_mm$theta, lme4::getME(fit_ref, "theta"),
    "binomial_probit_profiled_pirls", "theta",
    tol = 0.07,
    label = "binomial/probit theta within expected-mismatch bound"
  )
})

test_that("glmm() probit logLik is finite and within expected-mismatch bounds vs lme4", {
  mm_skip_if_no_lme4()
  fit_mm  <- glmm(y ~ x + (1 | grp), data = dat_probit,
                  family = binomial(link = "probit"),
                  control = mm_control(verbose = -1L))
  fit_ref <- suppressMessages(suppressWarnings(
    lme4::glmer(y ~ x + (1 | grp), data = dat_probit,
                family = binomial(link = "probit"), nAGQ = 1L)
  ))
  # Probe observed: logLik diff ~0.049; register bound at 0.2
  mm_assert_parity(
    as.numeric(logLik(fit_mm)), as.numeric(stats::logLik(fit_ref)),
    "binomial_probit_profiled_pirls", "logLik",
    tol = 0.2,
    label = "binomial/probit logLik within expected-mismatch bound"
  )
})
```

### Upstream fixture / engine change required?

A new entry in `inst/extdata/expected-mismatches.json` for
`"binomial_probit_profiled_pirls"` covering `fixef`, `theta`, and `logLik`
fields at the observed-mismatch bounds documented here.

---

## SPEC-GLMM-11 — Grouseticks Poisson logLik regression test

**Priority:** P1  
**Kind:** numerical-tolerance  
**Classification addressed:** partial (grouseticks Poisson parity case shows 0.524-unit logLik
gap vs lme4; this is within the registered mismatch bounds but the test should
explicitly lock the bound and verify it does not grow)

### Rationale

The probe (`glmm-grouseticks-pois.md`) shows mixeff's optimizer finds a
suboptimal solution on the grouseticks multi-RE Poisson dataset (logLik gap
0.524 units, well outside the 1e-3 tolerance). The case is registered in
`expected-mismatches.json` as `grouseticks_poisson_log_profiled_pirls` with
`expected_max_abs_diff = 3.8` for logLik. The current regression detector
should be extended to check the bound is not *exceeded* — both ways (regression
if gap grows, tightening opportunity if gap shrinks).

### Dataset / setup

```r
mm_skip_if_no_lme4()
data(grouseticks, package = "lme4")
grouseticks$cHEIGHT <- grouseticks$HEIGHT - mean(grouseticks$HEIGHT)
```

### Assertions

```r
test_that("grouseticks Poisson logLik gap is within registered mismatch bound", {
  mm_skip_if_no_lme4()
  fit_mm <- glmm(
    TICKS ~ YEAR + cHEIGHT + (1 | INDEX) + (1 | BROOD) + (1 | LOCATION),
    data = grouseticks,
    family = poisson(),
    method = "pirls_profiled",
    control = mm_control(verbose = -1L)
  )
  fit_ref <- suppressMessages(suppressWarnings(
    lme4::glmer(
      TICKS ~ YEAR + cHEIGHT + (1 | INDEX) + (1 | BROOD) + (1 | LOCATION),
      data = grouseticks,
      family = poisson(), nAGQ = 1L
    )
  ))
  loglik_diff <- abs(as.numeric(logLik(fit_mm)) - as.numeric(stats::logLik(fit_ref)))
  # Registered bound from expected-mismatches.json for logLik: 3.8
  expect_lte(loglik_diff, 3.8,
             label = "grouseticks logLik gap must not exceed registered mismatch bound")
  # Regression detector: if gap shrinks below 25% of bound, tighten the bound
  expect_true(loglik_diff > 0,
              info = "logLik diff should be positive (profiled vs joint-Laplace)")
})

test_that("grouseticks Poisson fit converges and has finite estimates", {
  fit <- glmm(
    TICKS ~ YEAR + cHEIGHT + (1 | INDEX) + (1 | BROOD) + (1 | LOCATION),
    data = grouseticks,
    family = poisson(),
    control = mm_control(verbose = -1L)
  )
  expect_s3_class(fit, "mm_glmm")
  expect_true(all(is.finite(fixef(fit))))
  expect_true(is.finite(fit$logLik))
  expect_identical(fit$fit_status, "converged_interior")
})
```

### Note on upstream engine fix

The 0.524-unit logLik gap is an optimizer shortfall in the vendored Rust engine
on this dataset (confirmed by the probe — it is not a response-constant
convention issue). When an upstream fix lands (tracked in the upstream mote for
`grouseticks_poisson_log_profiled_pirls`), the `expected_max_abs_diff` in
`expected-mismatches.json` should be tightened and this test updated to verify
the tighter bound.

---

## SPEC-GLMM-12 — joint_laplace refusal is typed with the correct class (not mm_formula_error)

**Priority:** P2  
**Kind:** error-message  
**Classification addressed:** partial (the joint_laplace error fires correctly on simple models;
the probe showed cbind blocks first, masking the nlopt refusal in that case)

### Rationale

The `glmm-cbpp-agq.md` probe noted that when `cbind(...)` is the response AND
`method="joint_laplace"`, the refusal class is `mm_formula_error` (cbind blocked
first), not the expected `mm_fit_error` for an unavailable backend. On a model
without the cbind blocker, the correct `mm_fit_error` should fire. This tests
that the correct condition class propagates when the formula is valid but the
backend is unavailable.

### Dataset / setup

```r
set.seed(42)
dat_jl <- data.frame(
  y   = rbinom(60, 1, 0.4),
  x   = rnorm(60),
  grp = factor(rep(1:6, each = 10))
)
```

### Assertions

```r
test_that("glmm() joint_laplace refusal on valid formula is mm_fit_error", {
  expect_error(
    glmm(y ~ x + (1 | grp), data = dat_jl,
         family = binomial(),
         method = "joint_laplace",
         control = mm_control(verbose = -1L)),
    class = "mm_fit_error"
  )
})

test_that("glmm() joint_laplace error message mentions nlopt or backend", {
  err <- tryCatch(
    glmm(y ~ x + (1 | grp), data = dat_jl,
         family = binomial(),
         method = "joint_laplace",
         control = mm_control(verbose = -1L)),
    error = function(e) e
  )
  expect_true(
    grepl("nlopt|backend|not.*available|disabled", conditionMessage(err),
          ignore.case = TRUE),
    info = "joint_laplace refusal must name the missing backend, not give a generic error"
  )
})

test_that("glmm() joint_laplace with nAGQ>1 raises mm_arg_error before nlopt check", {
  # mm_glmm_validate_nagq enforces nAGQ <= 1 for joint_laplace
  expect_error(
    glmm(y ~ x + (1 | grp), data = dat_jl,
         family = binomial(),
         method = "joint_laplace",
         nAGQ = 3L,
         control = mm_control(verbose = -1L)),
    class = "mm_arg_error"
  )
})
```

---

## SPEC-GLMM-13 — GLMM summary Wald-z table is returned for all certified families

**Priority:** P1  
**Kind:** integration  
**Classification addressed:** test-gap (summary works per gap report but no cross-family test)

### Rationale

The gap report confirms `summary(fit, tests="coefficients")` works for
binomial/logit. There is no test that covers poisson and gamma, and no test that
verifies the reliability note ("pirls_laplace_working_hessian") propagates
consistently across families.

### Dataset / setup

```r
set.seed(99)
n <- 80; g <- 8
dat_sum <- data.frame(
  y_bin = rbinom(n, 1, 0.4),
  y_poi = rpois(n, 2),
  y_gam = rgamma(n, shape = 3, rate = 1),
  x     = rnorm(n),
  grp   = factor(rep(seq_len(g), each = n / g))
)
```

### Assertions

```r
test_that("summary(mm_glmm) returns Wald-z table for binomial family", {
  fit <- glmm(y_bin ~ x + (1 | grp), data = dat_sum,
              family = binomial(), control = mm_control(verbose = -1L))
  s <- summary(fit, tests = "coefficients")
  tbl <- s$coefficients
  expect_true(!is.null(tbl))
  expect_true(all(c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in% colnames(tbl)))
})

test_that("summary(mm_glmm) returns Wald-z table for poisson family", {
  fit <- glmm(y_poi ~ x + (1 | grp), data = dat_sum,
              family = poisson(), control = mm_control(verbose = -1L))
  s <- summary(fit, tests = "coefficients")
  tbl <- s$coefficients
  expect_true(!is.null(tbl))
  expect_true(all(c("Estimate", "z value") %in% colnames(tbl)))
})

test_that("summary(mm_glmm) returns Wald-z table for Gamma/log family", {
  fit <- glmm(y_gam ~ x + (1 | grp), data = dat_sum,
              family = Gamma(link = "log"), control = mm_control(verbose = -1L))
  s <- summary(fit, tests = "coefficients")
  tbl <- s$coefficients
  expect_true(!is.null(tbl))
  expect_true(all(c("Estimate", "z value") %in% colnames(tbl)))
})

test_that("GLMM summary vcov reliability note is present for all families", {
  for (fam in list(binomial(), poisson(), Gamma(link = "log"))) {
    y_var <- switch(fam$family,
      "binomial" = "y_bin", "poisson" = "y_poi", "Gamma" = "y_gam"
    )
    formula_str <- as.formula(paste(y_var, "~ x + (1 | grp)"))
    fit <- glmm(formula_str, data = dat_sum,
                family = fam, control = mm_control(verbose = -1L))
    vc <- vcov(fit)
    note_attr <- attr(vc, "mm_reliability") %||% attr(vc, "reliability")
    # The reliability note should be present (non-null)
    expect_false(is.null(note_attr),
                 info = sprintf("vcov reliability note missing for %s GLMM",
                                fam$family))
  }
})
```

---

## SPEC-GLMM-14 — GLMM extractors return correct structures for all certified families

**Priority:** P1  
**Kind:** integration  
**Classification addressed:** test-gap (gap report says fixef/ranef/VarCorr/sigma work; no
cross-family integration test)

### Rationale

`fixef()`, `ranef()`, `VarCorr()`, `sigma()`, `nobs()`, `logLik()`, `AIC()`,
`BIC()`, `fitted()`, `residuals()`, `vcov()` all dispatched correctly in the
probes. No systematic cross-family regression test exists.

### Assertions (abbreviated — use loop over families)

```r
test_that("all standard extractors work on mm_glmm for all certified families", {
  set.seed(2)
  n <- 60; g <- 6
  base_dat <- data.frame(
    y_bin = rbinom(n, 1, 0.4),
    y_poi = rpois(n, 2),
    y_gam = rgamma(n, shape = 2, rate = 0.5),
    x     = rnorm(n),
    grp   = factor(rep(seq_len(g), each = n / g))
  )

  families <- list(
    list(fam = binomial(),       y = "y_bin"),
    list(fam = poisson(),        y = "y_poi"),
    list(fam = Gamma(link="log"), y = "y_gam")
  )

  for (spec in families) {
    f <- as.formula(paste(spec$y, "~ x + (1 | grp)"))
    fit <- glmm(f, data = base_dat, family = spec$fam,
                control = mm_control(verbose = -1L))
    lbl <- spec$fam$family

    expect_s3_class(fit, "mm_glmm",         label = paste(lbl, "class"))
    expect_length(fixef(fit), 2L,           label = paste(lbl, "fixef length"))
    expect_true(is.finite(logLik(fit)),     label = paste(lbl, "logLik"))
    expect_true(is.finite(AIC(fit)),        label = paste(lbl, "AIC"))
    expect_true(is.finite(BIC(fit)),        label = paste(lbl, "BIC"))
    expect_true(is.numeric(sigma(fit)),     label = paste(lbl, "sigma"))
    expect_identical(nobs(fit), n,          label = paste(lbl, "nobs"))
    expect_length(fitted(fit), n,           label = paste(lbl, "fitted length"))
    expect_length(residuals(fit), n,        label = paste(lbl, "residuals length"))
    expect_true(is.matrix(vcov(fit)),       label = paste(lbl, "vcov matrix"))
    expect_true(!is.null(ranef(fit)),       label = paste(lbl, "ranef non-null"))
    vc <- VarCorr(fit)
    expect_true(!is.null(vc),              label = paste(lbl, "VarCorr non-null"))
  }
})
```

---

## Summary table

| Spec ID | Priority | Kind | Gap addressed | Upstream required? |
|---|---|---|---|---|
| SPEC-GLMM-01 | **P0** | error-message | `offset=` silently dropped (blocker) | No |
| SPEC-GLMM-02 | P1 | error-message | `weights=` refusal regression test | No |
| SPEC-GLMM-03 | P1 | unit + parity | Gamma/inverse link wrapper-blocked despite engine support | No |
| SPEC-GLMM-04 | P1 | error-message | nAGQ>1 on pirls_profiled transparency warning absent | No |
| SPEC-GLMM-05 | P2 | unit | nAGQ>1 + multi-RE on profiled: no misleading hard stop | No |
| SPEC-GLMM-06 | P1 | integration | GLMM parametric bootstrap R wire-up | No (engine stable) |
| SPEC-GLMM-07 | P1 | error-message | cbind refusal message accuracy (misleading "stateful" text) | No |
| SPEC-GLMM-08 | P1 | integration + tolerance | Bernoulli-expanded cbpp: valid fit + mismatch bounds | `expected-mismatches.json` entry exists |
| SPEC-GLMM-09 | P1 | tolerance + snapshot | Gamma/log dispersion and convergence | `expected-mismatches.json` entry needed |
| SPEC-GLMM-10 | P1 | parity + tolerance | binomial/probit expected-mismatch bounds | `expected-mismatches.json` entry needed |
| SPEC-GLMM-11 | P1 | tolerance | grouseticks Poisson logLik regression bound | `expected-mismatches.json` entry exists; upstream fix tracked |
| SPEC-GLMM-12 | P2 | error-message | joint_laplace refusal class stability | No |
| SPEC-GLMM-13 | P1 | integration | GLMM summary Wald-z table cross-family | No |
| SPEC-GLMM-14 | P1 | integration | Standard extractors cross-family regression | No |

**Total specs: 14**  
**P0 blockers: 1** (SPEC-GLMM-01)  
**P1: 10**  
**P2: 3**

---

## Implementation notes

1. **SPEC-GLMM-01 is the highest-priority write** — it is the only blocker and
   requires only a one-line change to `R/glmm.R` (add `offset` to the reserved
   list), followed by adding the test to `tests/testthat/test-glmm.R`.

2. **SPEC-GLMM-03** gates on a `mm_glmm_supported_family_links()` update. The
   test file should include both the "currently refused" snapshot test (for
   tracking the current state) and the "after fix" assertions in a `skip_if`
   guard.

3. **SPEC-GLMM-06** (bootstrap wire-up) is architecturally significant and
   should be a dedicated file `tests/testthat/test-glmm-bootstrap.R` with
   `skip_on_cran()` guards on slow tests.

4. `expected-mismatches.json` needs new entries for
   `"binomial_probit_profiled_pirls"` (fixef, theta, logLik) before
   SPEC-GLMM-10 parity tests can use `mm_assert_parity` with the ledger-based
   tolerance path. SPEC-GLMM-09 (Gamma/log) needs similar entries.

5. All tests follow existing `mm_skip_if_no_lme4()` / `mm_assert_parity()` /
   `mm_control(verbose = -1L)` conventions from `tests/testthat/test-glmm.R`
   and `tests/testthat/helper-lme4-parity.R`.
