# Test Specification: `simulate` / `refit` / `update` family

Spec date: 2026-05-31  
Gap report: `assessment/gap/simulate-refit-update.md`  
Reference surface: lme4 2.0.1, lmerTest 3.2.1  
Package under test: mixeff (main branch)

Each spec covers one in-scope gap (classified `in-scope-missing`, `partial`, or
`test-gap` in the gap report). Specs are ordered by priority then by logical
grouping. Upstream-fixture requirements are called out explicitly.

---

## Spec SRU-01 — `simulate(newparams=)` must error with a typed refusal

| Field | Value |
|---|---|
| **Name** | `simulate_newparams_must_refuse_not_silently_ignore` |
| **Priority** | **P0 blocker** |
| **Kind** | error-message |
| **Dataset/formula** | `sleepstudy`, `Reaction ~ Days + (Days|Subject)` |
| **Gap classification** | in-scope-missing (contract violation — PRD §3 "no silent surgery") |

### Assertion

`simulate(fit, seed = 42, newparams = list(beta = c(9999, 9999), theta = 10, sigma = 999))`
must **not** silently return output identical to `simulate(fit, seed = 42)`. It must
either:

- **(preferred)** honor `newparams` and return simulated draws centered near `9999`; or  
- raise a condition of class `mm_inference_unavailable` or `mm_arg_error` with a
  non-empty message that names the ignored argument.

The current behavior — returning the default-param simulation with no warning or
error — is a direct violation of the "no silent model surgery" contract.

```r
library(lme4)
library(mixeff)
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))

default_sim <- simulate(fit, seed = 42)
newparams_sim <- simulate(fit, seed = 42,
                          newparams = list(beta = c(9999, 9999),
                                           theta = 10, sigma = 999))

# Contract: output must differ from the default simulation OR an error must
# have been raised. Both are acceptable; silent identity is not.
expect_false(identical(default_sim, newparams_sim))
# Preferred path — if the feature is not yet implemented:
expect_error(
  simulate(fit, seed = 42,
           newparams = list(beta = c(9999, 9999), theta = 10, sigma = 999)),
  class = c("mm_inference_unavailable", "mm_arg_error")
)
```

### What it guards against

Regression of the silent-ignore bug: `newparams` absorbed by `...` with no
effect. This is the primary path for power analysis and simulation studies.

---

## Spec SRU-02 — `simulate(newdata=)` must error with a typed refusal

| Field | Value |
|---|---|
| **Name** | `simulate_newdata_must_refuse_not_silently_ignore` |
| **Priority** | **P0 blocker** |
| **Kind** | error-message |
| **Dataset/formula** | `sleepstudy`, `Reaction ~ Days + (Days|Subject)` |
| **Gap classification** | in-scope-missing (contract violation — PRD §3 "no silent surgery") |

### Assertion

`simulate(fit, newdata = sleepstudy[1:20, ])` must **not** return a 180-row
data frame. It must either:

- honor `newdata` and return a 20-row data frame; or  
- raise `mm_inference_unavailable` or `mm_arg_error` naming `newdata`.

Returning the full-data simulation while accepting `newdata` silently is
asymmetric with `predict.mm_lmm`, which does honor `newdata`, and violates §3.

```r
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))

err <- tryCatch(
  simulate(fit, newdata = sleepstudy[1:20, ], seed = 1),
  mm_inference_unavailable = function(cnd) cnd,
  mm_arg_error = function(cnd) cnd
)
# Must be a typed refusal (not 180 rows of silent output)
expect_s3_class(err, c("mm_inference_unavailable", "mm_arg_error"),
                exact = FALSE)
expect_match(conditionMessage(err), "newdata", fixed = TRUE)
```

### What it guards against

Silent asymmetry between `predict(newdata=)` (honored) and
`simulate(newdata=)` (silently dropped).

---

## Spec SRU-03 — `simulate(use.u=)` must error with a typed refusal

| Field | Value |
|---|---|
| **Name** | `simulate_use_u_must_refuse_not_silently_ignore` |
| **Priority** | P1 |
| **Kind** | error-message |
| **Dataset/formula** | `sleepstudy`, `Reaction ~ Days + (Days|Subject)` |
| **Gap classification** | partial (contract violation — absorbed by `...`) |

### Assertion

`simulate(fit, use.u = TRUE)` and `simulate(fit, use.u = FALSE)` must either:

- map `use.u = TRUE` to `re.form = NULL` (conditional) and
  `use.u = FALSE` to `re.form = ~0` (population-level), matching lme4 behavior;
  or  
- raise a typed condition naming the unrecognized argument.

The current behavior — accepted silently with no effect — is a contract violation.

```r
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))

# Either: typed refusal
err <- tryCatch(
  simulate(fit, use.u = TRUE, seed = 1),
  mm_inference_unavailable = function(cnd) cnd,
  mm_arg_error = function(cnd) cnd
)
expect_s3_class(err, c("mm_inference_unavailable", "mm_arg_error"),
                exact = FALSE)

# Or (preferred path — map to re.form):
# s_true  <- simulate(fit, use.u = TRUE,  seed = 1)   # conditional
# s_false <- simulate(fit, use.u = FALSE, seed = 1)   # population
# expect_false(identical(s_true, s_false))
```

### What it guards against

Silent pass-through of lme4 legacy alias that currently has no effect, allowing
users to believe they are simulating conditionally when they are not.

---

## Spec SRU-04 — `refit()` must accept a single-column data.frame from `simulate()`

| Field | Value |
|---|---|
| **Name** | `refit_accepts_simulate_dataframe_column` |
| **Priority** | **P0 blocker** |
| **Kind** | integration |
| **Dataset/formula** | phase4 synthetic data, `y ~ x + z + (1|subject)` |
| **Gap classification** | partial (major — breaks documented lme4 bootstrap loop) |

### Assertion

The documented lme4 bootstrap idiom `lapply(sims, refit, object = fit)` must
work without modification. This requires `refit.mm_lmm` to accept a
single-column data frame as `newresp`, extracting the column as a numeric
vector.

```r
fit <- lmm(y ~ x + z + (1 | subject), df,
           control = mm_control(verbose = -1))
sims <- simulate(fit, nsim = 3, seed = 99)

# Direct data.frame column (currently fails with mm_arg_error)
ref1 <- refit(fit, sims[, 1, drop = FALSE])  # must NOT error
expect_s3_class(ref1, "mm_lmm")
expect_equal(nobs(ref1), nobs(fit))

# lapply idiom (the documented lme4 pattern)
refits <- lapply(sims, refit, object = fit)  # must NOT error
expect_length(refits, 3L)
expect_true(all(vapply(refits, inherits, logical(1), what = "mm_lmm")))

# Numeric vector path must still work
ref2 <- refit(fit, sims[[1]])
expect_s3_class(ref2, "mm_lmm")
```

### What it guards against

Regression of the data.frame-column rejection; ensures the standard lme4
bootstrap loop works verbatim without requiring users to call `sims[[i]]`
instead of passing `sims` directly.

---

## Spec SRU-05 — `refit(newweights=)` must error with a typed refusal

| Field | Value |
|---|---|
| **Name** | `refit_newweights_must_refuse_not_silently_ignore` |
| **Priority** | P1 |
| **Kind** | error-message |
| **Dataset/formula** | `sleepstudy`, `Reaction ~ Days + (Days|Subject)` |
| **Gap classification** | partial (contract violation — absorbed by `...`) |

### Assertion

`refit(fit, y, newweights = rep(1, nobs(fit)))` must either:

- honor `newweights` and re-fit with the supplied weights; or  
- raise `mm_inference_unavailable` or `mm_arg_error` with a message naming
  `newweights`.

The current behavior — silently proceeding without using the supplied weights —
is a §3 contract violation.

```r
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))
sims <- simulate(fit, nsim = 1, seed = 7)

err <- tryCatch(
  refit(fit, sims[[1]], newweights = rep(1, nobs(fit))),
  mm_inference_unavailable = function(cnd) cnd,
  mm_arg_error = function(cnd) cnd
)
expect_s3_class(err, c("mm_inference_unavailable", "mm_arg_error"),
                exact = FALSE)
expect_match(conditionMessage(err), "newweights", ignore.case = TRUE)
```

### What it guards against

Silent `newweights` drop — user believes they changed observation weights in a
bootstrap loop but actually did not.

---

## Spec SRU-06 — `simulate.mm_glmm` dispatches and returns the correct structure

| Field | Value |
|---|---|
| **Name** | `simulate_mm_glmm_dispatches_and_returns_dataframe` |
| **Priority** | **P0 blocker** |
| **Kind** | integration |
| **Dataset/formula** | cbpp (expanded to binary rows), `y ~ period + (1|herd)`, `family = binomial()` |
| **Gap classification** | in-scope-missing (major — PRD §10 Phase 4 explicitly lists GLMM+simulate) |

### Assertion

`simulate(glmm_fit, nsim = 3, seed = 42)` on an `mm_glmm` object must:

1. not dispatch to `simulate.default` or raise "no applicable method";
2. return a `data.frame` with `nrow == nobs(glmm_fit)` and `ncol == 3`;
3. use column names `sim_1`, `sim_2`, `sim_3`;
4. contain integer or numeric values consistent with the binomial family
   (values in {0, 1} for Bernoulli, or non-negative integers for aggregated
   binomial);
5. carry `attr(out, "seed") == 42`.

```r
pair <- mk_cbpp_glmm_fit()   # from helper or inline
fit  <- pair$fit

sims <- simulate(fit, nsim = 3, seed = 42)

expect_s3_class(sims, "data.frame")
expect_equal(nrow(sims), nobs(fit))
expect_equal(ncol(sims), 3L)
expect_equal(names(sims), c("sim_1", "sim_2", "sim_3"))
expect_identical(attr(sims, "seed"), 42L)  # or 42

# Values must be 0/1 (Bernoulli rows)
expect_true(all(vapply(sims, function(col) all(col %in% c(0L, 1L)), logical(1))))

# Reproducibility
sims2 <- simulate(fit, nsim = 3, seed = 42)
expect_equal(sims, sims2)
```

### Upstream fixture requirement

Requires `mm_simulate_once` (or a new `simulate.mm_glmm` method) to draw from
the GLMM inverse-link / family distribution. The crate `lib.rs` must expose the
fitted linear predictor or conditional mean for GLMM objects, or the R side
must reconstruct it from `fixef` and `ranef`. This is an **upstream-blocked**
spec until the engine exposes a GLMM simulation surface.

---

## Spec SRU-07 — `refit.mm_glmm` dispatches and preserves the GLMM class

| Field | Value |
|---|---|
| **Name** | `refit_mm_glmm_dispatches_and_returns_mm_glmm` |
| **Priority** | **P0 blocker** |
| **Kind** | integration |
| **Dataset/formula** | cbpp (expanded), `y ~ period + (1|herd)`, `family = binomial()` |
| **Gap classification** | in-scope-missing (major — PRD §10 Phase 4) |

### Assertion

`refit(glmm_fit, new_y)` must:

1. dispatch to a `refit.mm_glmm` method (not error with "no applicable method");
2. return an `mm_glmm` object;
3. preserve `$family`, `$method`, and formula from the original fit;
4. record `$refit$source == "refit"`.

```r
pair <- mk_cbpp_glmm_fit()
fit  <- pair$fit
sims <- simulate(fit, nsim = 1, seed = 5)  # requires SRU-06 to pass first

ref <- refit(fit, sims[[1]])
expect_s3_class(ref, "mm_glmm")
expect_identical(ref$family$family, fit$family$family)
expect_identical(ref$family$link,   fit$family$link)
expect_equal(nobs(ref), nobs(fit))
expect_identical(deparse1(formula(ref)), deparse1(formula(fit)))
expect_identical(ref$refit$source, "refit")
```

### Upstream fixture requirement

Depends on SRU-06. Requires the R layer to call `glmm()` with the new response
vector and the stored family, formula, and data. No new Rust surface needed
beyond what `glmm()` already uses.

---

## Spec SRU-08 — Simulate→refit→parametric-bootstrap round-trip for LMM

| Field | Value |
|---|---|
| **Name** | `simulate_refit_bootstrap_roundtrip_lmm` |
| **Priority** | P1 |
| **Kind** | integration |
| **Dataset/formula** | phase4 synthetic, full `y ~ x + z + (1|subject)` vs. reduced `y ~ x + (1|subject)` |
| **Gap classification** | test-gap (happy path exists but SRU-04's data.frame-column fix must be in place) |

### Assertion

The standard parametric-bootstrap workflow using `simulate`, `refit` (with
single-column data.frame), and `parametric_bootstrap` must complete end-to-end
without manual vector extraction.

```r
full    <- lmm(y ~ x + z + (1 | subject), df, REML = FALSE,
               control = mm_control(verbose = -1))
reduced <- lmm(y ~ x     + (1 | subject), df, REML = FALSE,
               control = mm_control(verbose = -1))

# Simulate from the reduced (null) model
sims <- simulate(reduced, nsim = 5, seed = 123)

# refit the full model on each simulated dataset — lme4 idiom
lrt_null <- vapply(seq_len(ncol(sims)), function(i) {
  r <- refit(reduced, sims[, i, drop = FALSE])   # data.frame column form
  f <- refit(full,    sims[, i, drop = FALSE])
  pmax(0, deviance(r) - deviance(f))
}, numeric(1))

expect_length(lrt_null, 5L)
expect_true(all(is.finite(lrt_null)))
expect_true(all(lrt_null >= 0))

# Also verify parametric_bootstrap() agrees in direction
boot <- parametric_bootstrap(reduced, full, nsim = 5, seed = 123)
expect_s3_class(boot, "mm_parametric_bootstrap")
expect_true(boot$status %in% c("available", "not_assessed"))
```

### What it guards against

Ensures the three-component workflow (simulate → refit → bootstrap) works
end-to-end with the lme4-idiomatic column form, after SRU-04 is implemented.

---

## Spec SRU-09 — `simulate` `re.form = ~0` is handled (alias for `NA`)

| Field | Value |
|---|---|
| **Name** | `simulate_re_form_tilde_zero_population_level` |
| **Priority** | P1 |
| **Kind** | parity-vs-lme4 |
| **Dataset/formula** | `sleepstudy`, `Reaction ~ Days + (Days|Subject)` |
| **Gap classification** | partial (minor — `~0` not handled; `NA` works) |

### Assertion

`simulate(fit, re.form = ~0, seed = 1)` must behave identically to
`simulate(fit, re.form = NA, seed = 1)` rather than raising an error. lme4
accepts both `NA` and `~0` as the population-level signal.

```r
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))

s_na  <- simulate(fit, re.form = NA,  seed = 1)
s_tilde0 <- simulate(fit, re.form = ~0, seed = 1)  # currently errors

expect_s3_class(s_tilde0, "data.frame")
expect_equal(dim(s_tilde0), c(180L, 1L))
expect_equal(s_tilde0, s_na)
```

### What it guards against

Users translating lme4 code that uses `re.form = ~0` hit an opaque error;
`mm_prediction_target()` should map `~0` to `"population"`.

---

## Spec SRU-10 — `refitML()` user-facing verb

| Field | Value |
|---|---|
| **Name** | `refitML_exported_verb_available` |
| **Priority** | P2 |
| **Kind** | unit |
| **Dataset/formula** | phase4 synthetic, `y ~ x + z + (1|subject)`, REML = TRUE |
| **Gap classification** | partial (minor — internal behavior exists; no exported verb) |

### Assertion

`refitML(fit)` exported from the mixeff namespace must refit a REML model by
ML and return an `mm_lmm` with `REML == FALSE`.

```r
fit_reml <- lmm(y ~ x + z + (1 | subject), df, REML = TRUE,
                control = mm_control(verbose = -1))
expect_true(isTRUE(fit_reml$REML))

fit_ml <- refitML(fit_reml)

expect_s3_class(fit_ml, "mm_lmm")
expect_false(isTRUE(fit_ml$REML))
expect_equal(nobs(fit_ml), nobs(fit_reml))
expect_identical(deparse1(formula(fit_ml)), deparse1(formula(fit_reml)))
# ML logLik >= REML logLik on the same parameterization is not guaranteed,
# but they should be close
expect_true(abs(as.numeric(logLik(fit_ml)) -
                as.numeric(logLik(fit_reml))) < 5)
```

### What it guards against

Users who expect the standard `refitML(fit)` shorthand (used in lme4/lmerTest
vignettes) for REML→ML conversion when computing LRTs.

---

## Spec SRU-11 — `bootMer`-style general FUN bootstrap refuses with a typed message

| Field | Value |
|---|---|
| **Name** | `bootmer_general_fun_refuses_with_typed_message` |
| **Priority** | P1 |
| **Kind** | error-message |
| **Dataset/formula** | `sleepstudy`, `Reaction ~ Days + (Days|Subject)` |
| **Gap classification** | in-scope-missing (major — no general FUN harness) |

### Assertion

When a user calls `lme4::bootMer(fit, fixef, nsim = 2)` on an `mm_lmm` object,
the failure must be a typed mixeff diagnostic, not the raw S3-dispatch error
"no applicable method for 'isLMM' applied to … 'mm_lmm'". The preferred path is
a native `boot_mer()` / `bootstrap_fun()` equivalent; the minimum acceptable
path is a clear refusal at the mixeff boundary.

```r
library(lme4)
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))

# Option A: mixeff provides its own general FUN bootstrap
# result <- mm_bootstrap(fit, FUN = fixef, nsim = 5, seed = 1)
# expect_true(is.matrix(result) || is.list(result))

# Minimum contract — if not yet implemented, must refuse clearly:
err <- tryCatch(
  mm_bootstrap(fit, FUN = fixef, nsim = 5, seed = 1),
  mm_inference_unavailable = function(cnd) cnd,
  error = function(cnd) cnd
)
# Must not be a raw "no applicable method" dispatch error
if (inherits(err, "error")) {
  expect_false(grepl("no applicable method", conditionMessage(err)))
  expect_s3_class(err, c("mm_inference_unavailable", "mm_arg_error"),
                  exact = FALSE)
}
```

### Note

The native `parametric_bootstrap()` LRT path is already implemented and tested.
This spec concerns the additional general-FUN harness that lme4's `bootMer`
provides.

---

## Spec SRU-12 — `simulate` `newparams` numerical parity with lme4 (post-implementation)

| Field | Value |
|---|---|
| **Name** | `simulate_newparams_numerical_parity_lme4` |
| **Priority** | P1 (blocked on SRU-01 implementation) |
| **Kind** | parity-vs-lme4 / numerical-tolerance |
| **Dataset/formula** | `sleepstudy`, `Reaction ~ Days + (Days|Subject)` |
| **Gap classification** | in-scope-missing (numerical parity spec for when the feature is implemented) |

### Assertion

Once `newparams` is honored, the grand mean of `nsim = 500` simulated responses
under `beta = c(1000, 50), sigma = 10` must be within 3 SE of 1000 (the
intercept, which dominates), and within tolerance of lme4's grand mean under
the same hypothetical parameters.

```r
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))
params <- list(beta = c(1000, 50), theta = 1, sigma = 10)

s_mm <- simulate(fit, nsim = 500, seed = 77, newparams = params)
s_lm <- simulate(lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy),
                 nsim = 500, seed = 77, newparams = params)

grand_mm <- mean(unlist(s_mm))
grand_lm <- mean(unlist(s_lm))

# Mean should be near beta[1] = 1000
expect_true(abs(grand_mm - 1000) < 10)  # 10 = generous; SE ≈ 0.4
# Parity with lme4 grand mean
expect_true(abs(grand_mm - grand_lm) < 5)
```

### Tolerances

| Quantity | Tolerance | Rationale |
|---|---|---|
| `abs(grand_mm - 1000)` | < 10 | 3 SE at nsim=500 ≈ 1.2; 10 is conservative |
| `abs(grand_mm - grand_lm)` | < 5 | Monte Carlo noise |

---

## Spec SRU-13 — `simulate.mm_glmm` moments parity with lme4 (post-SRU-06)

| Field | Value |
|---|---|
| **Name** | `simulate_mm_glmm_moments_parity_lme4` |
| **Priority** | P1 (blocked on SRU-06) |
| **Kind** | parity-vs-lme4 / numerical-tolerance |
| **Dataset/formula** | cbpp expanded, `y ~ period + (1|herd)`, binomial |
| **Gap classification** | in-scope-missing (numerical parity for GLMM simulate) |

### Assertion

Once `simulate.mm_glmm` is implemented, the empirical proportion of 1s across
`nsim = 500` simulated Bernoulli draws must match lme4's empirical proportion
within tolerance.

```r
pair  <- mk_cbpp_glmm_fit()
fit_m <- pair$fit
fit_l <- lme4::glmer(y ~ period + (1 | herd), pair$data,
                     family = binomial())

s_m <- simulate(fit_m, nsim = 500, seed = 88)
s_l <- simulate(fit_l, nsim = 500, seed = 88, re.form = NULL)

prop_m <- mean(unlist(s_m))
prop_l <- mean(unlist(s_l))

# Both should be near the observed incidence rate
observed_rate <- mean(pair$data$y)
expect_true(abs(prop_m - observed_rate) < 0.05)
expect_true(abs(prop_m - prop_l) < 0.05)
```

### Tolerances

| Quantity | Tolerance |
|---|---|
| `abs(prop_m - observed_rate)` | < 0.05 |
| `abs(prop_m - prop_l)` | < 0.05 |

### Upstream fixture requirement

Requires SRU-06 (GLMM simulate surface). May require the crate to expose the
fitted conditional mean (`mu_i`) per observation for GLMM objects, or the R
layer to reconstruct it from `fixef + ranef` via the inverse-link.

---

## Spec SRU-14 — `refit` NA round-trip with `na.action` attr

| Field | Value |
|---|---|
| **Name** | `refit_na_action_attr_shortening_accepted` |
| **Priority** | P2 |
| **Kind** | integration |
| **Dataset/formula** | sleepstudy with 5 NAs inserted, `Reaction ~ Days + (Days|Subject)` |
| **Gap classification** | partial (minor) |

### Assertion

When a model is fit on data with NAs (excluded by `na.action`), and a simulated
response from `simulate()` carries an `na.action` attribute indicating omitted
rows, `refit()` must accept the shortened vector rather than requiring
`length(newresp) == nrow(original_data)`.

```r
ss_na <- sleepstudy
ss_na$Reaction[c(1, 10, 50, 100, 150)] <- NA
fit_na <- lmm(Reaction ~ Days + (Days | Subject), ss_na,
              control = mm_control(verbose = -1))

# simulate() on the NA-excluded fit returns nobs(fit_na) rows
sims <- simulate(fit_na, nsim = 1, seed = 3)
expect_equal(nrow(sims), nobs(fit_na))  # 175, not 180

# refit must accept this shorter vector
ref <- refit(fit_na, sims[[1]])
expect_s3_class(ref, "mm_lmm")
expect_equal(nobs(ref), nobs(fit_na))
```

---

## Spec SRU-15 — `allFit` refuses with a typed mixeff diagnostic (not raw dispatch error)

| Field | Value |
|---|---|
| **Name** | `allfit_refuses_with_typed_message_not_raw_dispatch` |
| **Priority** | P2 |
| **Kind** | error-message |
| **Dataset/formula** | phase4 synthetic, `y ~ x + z + (1|subject)` |
| **Gap classification** | in-scope-missing (minor — no cross-optimizer harness) |

### Assertion

`lme4::allFit(fit)` on an `mm_lmm` object must not produce the raw dispatch
error "no applicable method for 'isGLMM' applied to … 'mm_lmm'". Either a
native `mm_allFit()` equivalent must exist, or the object must carry an `isGLMM`
/ `isLMM` method that returns an informative condition.

```r
fit <- lmm(y ~ x + z + (1 | subject), df,
           control = mm_control(verbose = -1))

err <- tryCatch(
  lme4::allFit(fit),
  error = function(cnd) cnd
)
# Must not be a raw "no applicable method" dispatch error
if (inherits(err, "error")) {
  expect_false(grepl("no applicable method", conditionMessage(err),
                     fixed = TRUE),
               info = paste("raw dispatch error:", conditionMessage(err)))
}
```

---

## Priority Summary

| ID | Name (short) | Priority | Kind | Status |
|---|---|---|---|---|
| SRU-01 | `simulate(newparams=)` typed refusal | **P0** | error-message | in-scope-missing |
| SRU-02 | `simulate(newdata=)` typed refusal | **P0** | error-message | in-scope-missing |
| SRU-04 | `refit()` accepts simulate df column | **P0** | integration | partial |
| SRU-06 | `simulate.mm_glmm` dispatches | **P0** | integration | in-scope-missing |
| SRU-07 | `refit.mm_glmm` dispatches | **P0** | integration | in-scope-missing |
| SRU-03 | `simulate(use.u=)` typed refusal | P1 | error-message | partial |
| SRU-05 | `refit(newweights=)` typed refusal | P1 | error-message | partial |
| SRU-08 | simulate→refit→bootstrap round-trip | P1 | integration | test-gap |
| SRU-09 | `re.form = ~0` population-level alias | P1 | parity-vs-lme4 | partial |
| SRU-11 | `bootMer`-style FUN typed refusal | P1 | error-message | in-scope-missing |
| SRU-12 | `newparams` numerical parity (post-impl) | P1 | numerical-tolerance | blocked on SRU-01 |
| SRU-13 | GLMM simulate moments parity (post-impl) | P1 | numerical-tolerance | blocked on SRU-06 |
| SRU-10 | `refitML()` exported verb | P2 | unit | partial |
| SRU-14 | `refit` NA round-trip | P2 | integration | partial |
| SRU-15 | `allFit` typed refusal | P2 | error-message | in-scope-missing |

Total specs: **15**  
Top priority (single most important): **SRU-01** — `simulate(newparams=)` must
not silently ignore its argument. This is the only spec that is both a P0 blocker
and a direct PRD §3 contract violation in an otherwise-working code path (the
simulate surface is live; the bug is that `newparams` is absorbed by `...`
silently). Fixing it — either by honoring the argument or by raising a typed
refusal — requires a one-line guard at the top of `simulate.mm_lmm` and unblocks
the power-analysis use case for all LMM users.
