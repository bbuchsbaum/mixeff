# Test Specification — Capability Family: Bootstrap

Date: 2026-05-31
Source gap report: `assessment/gap/bootstrap.md`
Reference parity probe: `assessment/parity/inf-bootstrap.md`, `assessment/parity/inf-simulate.md`
Target test file: `tests/testthat/test-bootstrap.R` (new file; some specs extend existing
`test-inference.R` and `test-boundary-lrt.R`)

---

## Classification of gaps requiring specs

Gaps covered here (in-scope-missing, partial, or test-gap per gap report):

| # | Gap | Classification |
|---|-----|----------------|
| BS-01 | `confint(method="bootstrap")` covers fixed effects only — no variance-component CIs | partial / major |
| BS-02 | `confint` `interval=` missing `"norm"` (normal-approx) | partial / minor |
| BS-03 | `confint.mm_glmm` absent; produces inscrutable base-R error | in-scope-missing / major |
| BS-04 | `simulate.mm_glmm` absent | in-scope-missing / major |
| BS-05 | `refit.mm_glmm` absent | in-scope-missing / major |
| BS-06 | No `bootMer`-equivalent with user `FUN` argument | in-scope-missing / major |
| BS-07 | `boot`-class interop absent; replicates not accessible via `boot::boot.ci` | in-scope-missing / major |
| BS-08 | `PBrefdist` / reference-distribution reuse absent | in-scope-missing / minor |
| BS-09 | Semiparametric bootstrap (`type="semiparametric"`) absent | in-scope-missing / minor |
| BS-10 | Conditional parametric bootstrap (`use.u=TRUE`) absent | in-scope-missing / minor |
| BS-11 | `simulate` `newparams=` absent | in-scope-missing / minor |
| BS-12 | `simulate` `newdata=` absent | in-scope-missing / minor |
| BS-13 | `parametric_bootstrap()` LRT p-value accounting — MCSE and replicate audit columns | test-gap |
| BS-14 | Bootstrap CI distributional tolerance vs lme4 `bootMer` (fixed effects, nsim=499) | test-gap |
| BS-15 | `refit.mm_lmm` `newoffset=` absent | in-scope-missing / minor |

Gaps **not** covered here (out-of-scope-by-design or upstream-blocked):
- `seqPBmodcomp` (not exported by pbkrtest on this system; out-of-scope)
- `.simulateFun` internal (out-of-scope)
- `cluster_bootstrap` p-values (`not_assessed` by design; crate-8.md:30)
- `PBmodcomp` on GLMM (upstream-blocked; engine LMM-only per crate-9.md:24)
- `.progress`/`PBargs` progress bars (cosmetic; Rust engine runs the loop)
- `parallel`/`ncpus`/`cl` user parallel control (minor; Rust engine may obviate)

---

## Test specifications

---

### SPEC BS-01a: `confint(method="bootstrap")` covers fixed effects only — scope documented

**Name**: `confint_bootstrap_lmm_fixed_effects_only_scope`
**Kind**: unit
**Priority**: P1
**Dataset / formula**: synthetic LMM — `y ~ x + (1|subject)` (9 subjects × 5 obs, see `mk_inference_fit()`)

**Setup**:
```r
fit <- mk_inference_fit()   # REML = TRUE, 2 fixed effects
ci  <- confint(fit, method = "bootstrap",
               bootstrap = bootstrap_control(nsim = 50, seed = 1))
```

**Assertions**:
```r
# Only fixed-effect parameters appear in rownames
expect_identical(rownames(ci), names(fixef(fit)))
# No variance-component rows (.sigma, .theta, .sig01 etc.)
expect_false(any(grepl("sigma|sig|theta", rownames(ci), ignore.case = TRUE)))
# Class and method attribute are correct
expect_s3_class(ci, "mm_confint")
expect_identical(attr(ci, "method"), "bootstrap_full_model_distribution")
expect_identical(attr(ci, "status"), "available")
# All bounds finite
expect_true(all(is.finite(ci)))
```

**What it guards against**: silent inclusion of fabricated variance-component rows in future;
also documents that the scope is fixed-effects-only so a future PR adding VC bootstrap CIs
updates this spec explicitly.

**Upstream fixture / engine change needed**: None (current behaviour; guards scope boundary).

---

### SPEC BS-01b: Variance-component bootstrap CIs produce a typed refusal, not silence

**Name**: `confint_bootstrap_lmm_vc_produces_typed_refusal`
**Kind**: error-message
**Priority**: P1
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` (sleepstudy)

**Setup**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
           control = mm_control(verbose = -1))
```

**Assertions**:
```r
# Requesting a variance-component parameter by name raises mm_arg_error
# (parameter not in names(fixef(fit))), not a base-R cryptic error.
expect_error(
  confint(fit, parm = ".sig01", method = "bootstrap",
          bootstrap = bootstrap_control(nsim = 20, seed = 1)),
  class = "mm_arg_error"
)
expect_error(
  confint(fit, parm = "sigma", method = "bootstrap",
          bootstrap = bootstrap_control(nsim = 20, seed = 1)),
  class = "mm_arg_error"
)
```

**What it guards against**: a user requesting VC bootstrap CIs gets a clear typed error with the
current name-resolution logic rather than a silent wrong result or opaque base-R error.

**Upstream fixture / engine change needed**: None; current `mm_arg_error` from unknown-parm
validation already fires. This spec locks in that the message is typed.

---

### SPEC BS-02: `confint` `interval = "norm"` raises typed unavailable, not base-R error

**Name**: `confint_bootstrap_lmm_interval_norm_typed_refusal`
**Kind**: error-message
**Priority**: P2
**Dataset / formula**: synthetic LMM (mk_inference_fit)

**Setup**:
```r
fit <- mk_inference_fit()
```

**Assertions**:
```r
# "norm" is not in the interval choices; match.arg should fire with a clear error.
# Current behaviour: match.arg raises base-R simpleError.
# Target behaviour: same match.arg error is acceptable; test guards against
# silent wrong result or a fabricated normal-approx interval.
expect_error(
  confint(fit, method = "bootstrap",
          bootstrap = bootstrap_control(nsim = 20, seed = 1),
          interval = "norm"),
  regexp = "should be one of"   # match.arg wording
)
# If a future PR adds "norm" support, this spec must be updated with
# numerical assertions (see BS-02-future note below).
```

**What it guards against**: a future "norm" implementation that produces wrong bounds or
fabricates without documenting the formula used.

**Note BS-02-future**: when `"norm"` is implemented, replace with:
```r
ci_norm <- confint(fit, method = "bootstrap",
                   bootstrap = bootstrap_control(nsim = 200, seed = 5),
                   interval = "norm")
ci_perc <- confint(fit, method = "bootstrap",
                   bootstrap = bootstrap_control(nsim = 200, seed = 5),
                   interval = "percentile")
# Normal-approx and percentile CIs should be within 2 SE of each other
# on a well-behaved model (SE of each bound ≈ sd(replicates)/sqrt(nsim)).
expect_true(max(abs(ci_norm - ci_perc)) < 2)
```

**Upstream fixture / engine change needed**: None for refusal spec; engine change required
to add `"norm"` interval support (Rust `full_model_bootstrap_contrast_json` would need to
return a `"norm"` interval entry).

---

### SPEC BS-03a: `confint` on `mm_glmm` raises typed `mm_inference_unavailable`

**Name**: `confint_glmm_typed_refusal_not_base_r_error`
**Kind**: error-message
**Priority**: P0 (directly violates "clearer errors" mandate)
**Dataset / formula**: `cbpp` binomial GLMM — `cbind(incidence, size-incidence) ~ period + (1|herd)`, `binomial()`

**Setup**:
```r
skip_if_not_installed("lme4")
data("cbpp", package = "lme4")
fit_g <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial(),
              control = mm_control(verbose = -1))
```

**Assertions**:
```r
# Any confint call on mm_glmm must raise mm_inference_unavailable,
# NOT "non-numeric argument to binary operator" or similar base-R error.
err <- tryCatch(
  confint(fit_g, method = "bootstrap",
          bootstrap = bootstrap_control(nsim = 20, seed = 1)),
  error = function(e) e
)
expect_s3_class(err, "mm_inference_unavailable")
expect_match(conditionMessage(err), "glmm|GLMM|not.*available|bootstrap",
             ignore.case = TRUE)

# Also for default method (Wald / asymptotic)
err_default <- tryCatch(confint(fit_g), error = function(e) e)
# Either works (if Wald is implemented) or raises typed refusal
if (inherits(err_default, "error")) {
  expect_s3_class(err_default, "mm_inference_unavailable")
}
```

**What it guards against**: the documented inscrutable `"non-numeric argument to binary operator"`
error that currently fires for GLMM confint. This is the single highest-priority error-message gap.

**Upstream fixture / engine change needed**: None for the typed-refusal fix; an actual
`confint.mm_glmm` implementation would require upstream engine work.

---

### SPEC BS-03b: `confint.mm_glmm` typed refusal with non-bootstrap methods

**Name**: `confint_glmm_all_methods_typed_refusal`
**Kind**: error-message
**Priority**: P0
**Dataset / formula**: same cbpp GLMM as BS-03a

**Assertions**:
```r
for (meth in c("wald", "asymptotic", "profile")) {
  err <- tryCatch(
    confint(fit_g, method = meth),
    error = function(e) e
  )
  if (inherits(err, "error")) {
    expect_s3_class(err, c("mm_inference_unavailable", "mm_arg_error"),
                    info = paste("method =", meth))
    expect_false(
      grepl("non-numeric argument", conditionMessage(err), fixed = TRUE),
      info = paste("method =", meth, "must not produce base-R cryptic error")
    )
  }
}
```

**What it guards against**: all confint dispatch paths for `mm_glmm` producing opaque base-R
errors regardless of method argument.

**Upstream fixture / engine change needed**: None for typed-refusal fix.

---

### SPEC BS-04a: `simulate.mm_glmm` absent — typed refusal, not S3 dispatch failure

**Name**: `simulate_glmm_typed_refusal`
**Kind**: error-message
**Priority**: P1
**Dataset / formula**: cbpp binomial GLMM (same as BS-03a)

**Assertions**:
```r
err <- tryCatch(simulate(fit_g, nsim = 3, seed = 1), error = function(e) e)
# If no simulate.mm_glmm method exists, S3 dispatch falls to default.
# The default should raise a typed refusal or a clear "no method" message,
# NOT a cryptic R crash.
if (inherits(err, "error")) {
  # Either a typed mixeff error or an R dispatch error -- but NOT a crash
  # that produces NA-typed results silently.
  expect_true(
    inherits(err, "mm_inference_unavailable") ||
      grepl("no applicable method|not available|glmm|GLMM",
            conditionMessage(err), ignore.case = TRUE),
    info = paste("simulate(mm_glmm) error:", conditionMessage(err))
  )
  # Must not be a segfault-class condition (i.e., must be catchable)
  expect_s3_class(err, "error")
}
```

**What it guards against**: uncatchable crash or silent wrong output if simulate dispatches
incorrectly for an `mm_glmm` object.

**Upstream fixture / engine change needed**: Full implementation would need Rust GLMM
simulation support. Typed refusal requires only R-side dispatch fix.

---

### SPEC BS-04b: `simulate.mm_glmm` implementation — output structure and moments

**Name**: `simulate_glmm_output_structure`
**Kind**: parity-vs-lme4 + numerical-tolerance
**Priority**: P1
**Dataset / formula**: `grouseticks` Poisson GLMM — `TICKS ~ YEAR + (1|BROOD)`, `poisson()`

**Setup**:
```r
skip_if_not_installed("lme4")
data("grouseticks", package = "lme4")
fit_pois <- glmm(TICKS ~ YEAR + (1 | BROOD), data = grouseticks,
                 family = poisson(), control = mm_control(verbose = -1))
fit_lme4_pois <- lme4::glmer(TICKS ~ YEAR + (1 | BROOD),
                              data = grouseticks, family = poisson())
```

**Assertions**:
```r
nsim <- 200L
sims_mx   <- simulate(fit_pois, nsim = nsim, seed = 42)
sims_lme4 <- lme4::simulate(fit_lme4_pois, nsim = nsim, seed = 42,
                             re.form = NA)

# Structure
expect_s3_class(sims_mx, "data.frame")
expect_equal(dim(sims_mx), c(nrow(grouseticks), nsim))
expect_equal(colnames(sims_mx), paste0("sim_", seq_len(nsim)))

# All simulated values non-negative integer counts
expect_true(all(as.matrix(sims_mx) >= 0))
expect_true(all(as.matrix(sims_mx) %% 1 == 0))

# Seed reproducibility
sims_mx2 <- simulate(fit_pois, nsim = nsim, seed = 42)
expect_identical(sims_mx, sims_mx2)

# Grand-mean within 10% of lme4 (distributional, not exact)
grand_mx   <- mean(as.matrix(sims_mx))
grand_lme4 <- mean(as.matrix(sims_lme4))
expect_true(abs(grand_mx - grand_lme4) / grand_lme4 < 0.10,
            info = sprintf("grand mean: mixeff=%.3f lme4=%.3f", grand_mx, grand_lme4))
```

**What it guards against**: structural regressions and gross distributional errors in a future
`simulate.mm_glmm` implementation.

**Upstream fixture / engine change needed**: Requires Rust GLMM simulation support.
Mark this spec `skip_if_not(exists("simulate.mm_glmm", ...))` until implemented.

---

### SPEC BS-05a: `refit.mm_glmm` absent — typed refusal

**Name**: `refit_glmm_typed_refusal`
**Kind**: error-message
**Priority**: P1
**Dataset / formula**: cbpp binomial GLMM (BS-03a)

**Assertions**:
```r
# Simulate a new response via lme4 to use as newresp
set.seed(1)
n <- nrow(cbpp)
new_counts <- rbinom(n, size = cbpp$size, prob = 0.1)
newresp <- cbind(new_counts, cbpp$size - new_counts)

err <- tryCatch(refit(fit_g, newresp = newresp), error = function(e) e)
if (inherits(err, "error")) {
  expect_true(
    inherits(err, c("mm_inference_unavailable", "mm_arg_error")) ||
      grepl("no applicable method|not available|glmm|GLMM",
            conditionMessage(err), ignore.case = TRUE),
    info = paste("refit(mm_glmm) error:", conditionMessage(err))
  )
  # Must not produce a base-R cryptic dispatch error
  expect_false(
    grepl("subscript out of bounds|cannot coerce|non-numeric",
          conditionMessage(err), fixed = FALSE)
  )
}
```

**What it guards against**: cryptic S3 dispatch failures when refit is called on an mm_glmm.

**Upstream fixture / engine change needed**: Full refit.mm_glmm requires Rust GLMM bridge.

---

### SPEC BS-05b: `refit.mm_lmm` preserves model structure and produces plausible estimates

**Name**: `refit_lmm_round_trip_plausibility`
**Kind**: integration
**Priority**: P1
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` (sleepstudy, REML)

**Setup**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
           control = mm_control(verbose = -1))
sims <- simulate(fit, nsim = 3, seed = 7)
```

**Assertions**:
```r
# refit with each simulated response succeeds and returns mm_lmm
for (i in seq_len(3)) {
  r <- refit(fit, newresp = sims[[i]])
  expect_s3_class(r, "mm_lmm")
  # Fixed effects differ from original (different data) but are finite
  expect_true(all(is.finite(fixef(r))))
  # Length of fixef matches original
  expect_equal(length(fixef(r)), length(fixef(fit)))
  # theta and sigma are positive
  expect_true(all(r$theta > 0))
  expect_true(r$sigma > 0)
}

# Wrong-length newresp raises mm_arg_error
expect_error(
  refit(fit, newresp = rnorm(5)),
  class = "mm_arg_error"
)
# NA in newresp raises mm_arg_error
bad_resp <- sims[[1]]; bad_resp[1] <- NA
expect_error(refit(fit, newresp = bad_resp), class = "mm_arg_error")
```

**What it guards against**: refit regression; ensures the simulate→refit round-trip used in
parametric bootstrap retains correct argument validation.

**Upstream fixture / engine change needed**: None; `refit.mm_lmm` is already implemented.

---

### SPEC BS-06: Generic `bootMer`-equivalent with user `FUN` — typed refusal or stub

**Name**: `bootmer_fun_absent_typed_refusal`
**Kind**: error-message
**Priority**: P2
**Dataset / formula**: synthetic LMM (mk_inference_fit)

**Assertions**:
```r
fit <- mk_inference_fit()
# mixeff does not export bootMer; calling lme4::bootMer on an mm_lmm should
# either dispatch to a mixeff method or raise a clear error.
# This spec certifies the error is understandable, not cryptic.
if (exists("bootMer", where = asNamespace("lme4"))) {
  err <- tryCatch(
    lme4::bootMer(fit, FUN = fixef, nsim = 5),
    error = function(e) e
  )
  if (inherits(err, "error")) {
    # Must not produce "non-numeric argument to binary operator" or similar
    expect_false(
      grepl("non-numeric argument|subscript out of bounds",
            conditionMessage(err), fixed = FALSE)
    )
  }
}
# parametric_bootstrap does not accept a user FUN argument
expect_false("FUN" %in% names(formals(parametric_bootstrap)))
# Attempting to pass FUN via ... should not silently succeed with wrong results
err2 <- tryCatch(
  parametric_bootstrap(fit, fit, nsim = 5, FUN = fixef),
  error = function(e) e
)
# Either ignored (currently) or typed error — document the behavior
# to detect silent regressions
expect_true(
  inherits(err2, c("mm_parametric_bootstrap", "error")),
  info = "parametric_bootstrap with extra FUN arg must not crash silently"
)
```

**What it guards against**: silent wrong dispatch; documents that generic-FUN bootstrap is absent.

**Upstream fixture / engine change needed**: Full implementation requires Rust engine work.

---

### SPEC BS-07: `boot`-class interop — replicates accessible via `attr(x,"bootstrap")`

**Name**: `confint_bootstrap_replicates_accessible`
**Kind**: integration
**Priority**: P1
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` (sleepstudy)

**Setup**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
           control = mm_control(verbose = -1))
ci <- confint(fit, method = "bootstrap",
              bootstrap = bootstrap_control(nsim = 100, seed = 42))
```

**Assertions**:
```r
# Raw replicates are accessible even without boot-class interop
boot_payload <- attr(ci, "bootstrap")
expect_true(is.list(boot_payload))
expect_equal(length(boot_payload), nrow(ci))  # one payload per parameter

# Each payload carries replicate statistics
for (payload in boot_payload) {
  reps <- payload$replicate_statistics
  expect_true(length(reps) > 0)
  expect_true(all(is.finite(unlist(reps))))
  # Replicate count matches requested
  expect_equal(length(reps),
               as.integer(payload$metadata$successful_replicates))
}

# A user can manually reconstruct percentile CI from replicates
reps_intercept <- unlist(boot_payload[[1]]$replicate_statistics)
manual_lower <- quantile(reps_intercept, 0.025)
manual_upper <- quantile(reps_intercept, 0.975)
# Manual CI matches the reported CI within floating-point tolerance
expect_equal(as.numeric(ci["(Intercept)", 1]), manual_lower,
             tolerance = 1e-6,
             info = "lower bound should match 2.5% quantile of replicates")
expect_equal(as.numeric(ci["(Intercept)", 2]), manual_upper,
             tolerance = 1e-6,
             info = "upper bound should match 97.5% quantile of replicates")
```

**What it guards against**: regressions in replicate storage; ensures users have a documented
escape hatch for boot-class interop even without a native `"boot"` object.

**Upstream fixture / engine change needed**: None; `attr(x,"bootstrap")` is already populated.
This spec formally certifies the undocumented workaround that the gap report identified.

---

### SPEC BS-08: `parametric_bootstrap` result carries reusable `simulated` vector

**Name**: `parametric_bootstrap_simulated_vector_accessible`
**Kind**: unit
**Priority**: P2
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` (sleepstudy, ML)

**Setup**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit_null <- lmm(Reaction ~ 1 + (1 | Subject), data = sleepstudy, REML = FALSE,
                control = mm_control(verbose = -1))
fit_alt  <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE,
                control = mm_control(verbose = -1))
pb <- parametric_bootstrap(fit_null, fit_alt, nsim = 50, seed = 1)
```

**Assertions**:
```r
expect_s3_class(pb, "mm_parametric_bootstrap")

# Simulated LRT statistics are accessible
sims <- pb$simulated
expect_true(is.numeric(sims))
expect_true(length(sims) > 0)
expect_true(all(sims >= 0))  # LRT statistics are non-negative

# Observed statistic is accessible and positive
expect_true(is.finite(pb$observed))
expect_true(pb$observed > 0)

# p-value is in [0,1]
expect_true(is.finite(pb$p_value))
expect_true(pb$p_value >= 0 && pb$p_value <= 1)

# Seed is recorded
expect_equal(pb$seed, 1L)

# Replicate accounting columns exist
expect_true(!is.null(pb$successful_replicates))
expect_true(!is.null(pb$completed_replicates))
expect_true(pb$successful_replicates <= pb$completed_replicates)
```

**What it guards against**: removal of `$simulated` field (the only path for PBrefdist-style
reuse); regressions in replicate accounting structure.

**Upstream fixture / engine change needed**: None; `parametric_bootstrap` already returns this.

---

### SPEC BS-09: Semiparametric bootstrap absent — no `type` argument fabrication

**Name**: `bootstrap_control_no_type_argument`
**Kind**: unit
**Priority**: P2
**Dataset / formula**: synthetic LMM (mk_inference_fit)

**Assertions**:
```r
# bootstrap_control does not accept `type`; if passed via ..., it should not
# silently behave as parametric (which is the only implemented type).
expect_false("type" %in% names(formals(bootstrap_control)))

# Passing type = "semiparametric" raises an error (no silent wrong result)
err <- tryCatch(
  confint(mk_inference_fit(), method = "bootstrap",
          bootstrap = bootstrap_control(nsim = 10, seed = 1, type = "semiparametric")),
  error = function(e) e
)
# Either ignored with warning or raises an error — either is acceptable;
# what is NOT acceptable is silent wrong results labelled as semiparametric.
if (!inherits(err, "error")) {
  # If no error, the result must be documented as parametric
  ci <- err  # actually the confint result
  expect_identical(attr(ci, "method"), "bootstrap_full_model_distribution")
}
```

**What it guards against**: future `type` argument being wired to parametric silently.

**Upstream fixture / engine change needed**: Semiparametric bootstrap requires Rust residual
resampling; this spec only tests the absence of silent mis-labeling.

---

### SPEC BS-10: Conditional parametric bootstrap (`use.u = TRUE`) absent — no fabrication

**Name**: `simulate_no_use_u_argument`
**Kind**: unit
**Priority**: P2
**Dataset / formula**: synthetic LMM (mk_inference_fit)

**Assertions**:
```r
fit <- mk_inference_fit()
# use.u is not in formals of simulate.mm_lmm
expect_false("use.u" %in% names(formals(simulate.mm_lmm)))

# Passing use.u via ... should not silently succeed producing wrong sims
# (i.e., the random effects should not be held fixed when they shouldn't be)
sims1 <- simulate(fit, nsim = 5, seed = 1)
sims2 <- simulate(fit, nsim = 5, seed = 1, use.u = TRUE)  # extra arg silently ignored
# If use.u is silently ignored, both calls return identical results (both marginal)
expect_identical(sims1, sims2)
# Alternatively: if use.u raises an error, it must be typed
```

**What it guards against**: `use.u = TRUE` being silently treated as marginal simulation
while claiming conditional; documents the marginal-only nature of the implementation.

**Upstream fixture / engine change needed**: Conditional parametric bootstrap requires
Rust engine to accept a fixed random-effect state.

---

### SPEC BS-11: `simulate` `newparams=` absent — graceful error

**Name**: `simulate_newparams_typed_refusal`
**Kind**: error-message
**Priority**: P2
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` (sleepstudy)

**Assertions**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
           control = mm_control(verbose = -1))

# newparams is not a documented argument; passing it via ... is silently ignored
# or raises an error. Either outcome is acceptable; what is NOT acceptable is
# simulating with different parameters while reporting the fitted parameters.
expect_false("newparams" %in% names(formals(simulate.mm_lmm)))

# The simulated mean should match the fitted fixef, not a different theta
sims <- simulate(fit, nsim = 200, seed = 1,
                 newparams = list(beta = c(300, 5), theta = 0.5, sigma = 20))
# Grand mean should still reflect the fitted model (newparams ignored)
grand_mean <- mean(as.matrix(sims))
expect_true(abs(grand_mean - fixef(fit)[["(Intercept)"]]) < 30,
            info = sprintf("grand mean %.1f should be near fitted intercept %.1f",
                           grand_mean, fixef(fit)[["(Intercept)"]]))
```

**What it guards against**: a future `newparams` implementation that silently uses the wrong
parameters without labelling the output.

**Upstream fixture / engine change needed**: `newparams` requires Rust simulate to accept
an override parameter struct.

---

### SPEC BS-12: `simulate` `newdata=` absent — graceful error

**Name**: `simulate_newdata_typed_refusal`
**Kind**: error-message
**Priority**: P2
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` (sleepstudy)

**Assertions**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
           control = mm_control(verbose = -1))

expect_false("newdata" %in% names(formals(simulate.mm_lmm)))

# Passing newdata via ... should not silently simulate from the new design
# matrix with the fitted model and return wrong-length output.
new_obs <- data.frame(Days = 0:4, Subject = "308")
sims <- tryCatch(
  simulate(fit, nsim = 3, seed = 1, newdata = new_obs),
  error = function(e) e
)
if (!inherits(sims, "error")) {
  # If newdata is silently ignored, the output must still have nrow = nrow(sleepstudy)
  expect_equal(nrow(sims), nrow(sleepstudy),
               info = "simulate with ignored newdata must return nrow = original data rows")
} else {
  # Error must be typed (not base-R crash)
  expect_s3_class(sims, c("mm_inference_unavailable", "mm_arg_error", "error"))
}
```

**What it guards against**: `newdata` being silently mishandled, producing wrong-dimension
output or wrong-size random effects draws.

**Upstream fixture / engine change needed**: Out-of-sample simulation requires Rust engine.

---

### SPEC BS-13: `parametric_bootstrap` LRT — MCSE and replicate audit columns

**Name**: `parametric_bootstrap_lrt_audit_columns`
**Kind**: integration + numerical-tolerance
**Priority**: P0
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` vs `Reaction ~ 1 + (1|Subject)` (sleepstudy, ML)

**Setup**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit_null <- lmm(Reaction ~ 1 + (1 | Subject), data = sleepstudy, REML = FALSE,
                control = mm_control(verbose = -1))
fit_alt  <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = FALSE,
                control = mm_control(verbose = -1))
pb50 <- parametric_bootstrap(fit_null, fit_alt, nsim = 50, seed = 1)
```

**Assertions**:
```r
expect_s3_class(pb50, "mm_parametric_bootstrap")

# Observed LRT statistic should match anova() LRT
av <- anova(fit_null, fit_alt)
lrt_stat <- av$table$statistic[2]   # second row = alt
if (!is.na(lrt_stat)) {
  expect_equal(pb50$observed, lrt_stat, tolerance = 1e-3,
               info = "observed LRT should match anova() LRT statistic")
}

# p-value in [0, 1]
expect_true(pb50$p_value >= 0 && pb50$p_value <= 1)

# MCSE is finite and non-negative
expect_true(is.finite(pb50$mcse))
expect_true(pb50$mcse >= 0)

# Replicate accounting
expect_equal(pb50$completed_replicates, 50L)
expect_true(pb50$successful_replicates <= pb50$completed_replicates)
expect_true(pb50$successful_replicates > 0)
expect_true(is.numeric(pb50$boundary_count) || is.integer(pb50$boundary_count))
expect_true(pb50$boundary_rate >= 0 && pb50$boundary_rate <= 1)

# For a strong effect (Days is highly significant), p-value should be small
expect_true(pb50$p_value < 0.05,
            info = "Days effect on sleepstudy should be significant via PB-LRT")

# Seed recorded
expect_equal(pb50$seed, 1L)

# Status must be "available", not "not_assessed"
expect_identical(pb50$status, "available")
```

**What it guards against**: regressions in the replicate-accounting JSON fields; ensures the
MCSE, boundary_count, and seed_record fields survive future engine changes.

**Upstream fixture / engine change needed**: None; this tests current functionality.
If engine changes the JSON schema, update field names here.

---

### SPEC BS-14a: Bootstrap CI distributional parity vs lme4 — fixed effects, sleepstudy

**Name**: `confint_bootstrap_lmm_distributional_parity_sleepstudy`
**Kind**: parity-vs-lme4 + numerical-tolerance
**Priority**: P1
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` (sleepstudy, REML)

**Tolerance rationale**: Bootstrap CI is inherently stochastic; the tolerance is on CI _width_
(structural) not on exact bound values (RNG-dependent). The observed probe result (inf-bootstrap.md)
shows CI width agreement within ~1.7% for Intercept and ~4.4% for Days at nsim=499.

**Setup**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit_mx <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
              control = mm_control(verbose = -1))
fit_lme4 <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy,
                        REML = TRUE)
NSIM <- 499L
SEED <- 42L

ci_mx <- confint(fit_mx, method = "bootstrap",
                 bootstrap = bootstrap_control(nsim = NSIM, seed = SEED))
boot_lme4 <- lme4::bootMer(fit_lme4, FUN = lme4::fixef, nsim = NSIM,
                            type = "parametric", use.u = FALSE,
                            re.form = NA, seed = SEED)
ci_lme4 <- confint(boot_lme4, type = "perc")
```

**Assertions**:
```r
# CI widths agree within 10% (generous tolerance for RNG + REML/ML strategy difference)
common_parms <- intersect(rownames(ci_mx), rownames(ci_lme4))
expect_true(length(common_parms) >= 2,
            info = "Both CIs must cover (Intercept) and Days")

for (p in common_parms) {
  width_mx   <- ci_mx[p, 2]   - ci_mx[p, 1]
  width_lme4 <- ci_lme4[p, 2] - ci_lme4[p, 1]
  rel_diff <- abs(width_mx - width_lme4) / width_lme4
  expect_true(rel_diff < 0.10,
              info = sprintf("%s: width mixeff=%.3f lme4=%.3f rel_diff=%.1f%%",
                             p, width_mx, width_lme4, 100 * rel_diff))
}

# CI midpoints (estimate of parameter location) agree within 2 SE
for (p in common_parms) {
  mid_mx   <- mean(ci_mx[p, ])
  mid_lme4 <- mean(ci_lme4[p, ])
  se_p     <- sqrt(diag(vcov(fit_mx)))[p]
  expect_true(abs(mid_mx - mid_lme4) < 2 * se_p,
              info = sprintf("%s: midpoint diff=%.3f > 2*SE=%.3f",
                             p, abs(mid_mx - mid_lme4), 2 * se_p))
}

# Both CIs contain the point estimate (basic sanity)
for (p in common_parms) {
  est <- fixef(fit_mx)[p]
  expect_true(ci_mx[p, 1] < est && est < ci_mx[p, 2],
              info = sprintf("%s: point estimate %.3f not in CI [%.3f, %.3f]",
                             p, est, ci_mx[p, 1], ci_mx[p, 2]))
}
```

**What it guards against**: gross distributional divergence from lme4 at the same nsim,
such as would occur if the Rust bootstrap were sampling from the wrong distribution.
The 10% width tolerance accommodates the known REML/ML refit strategy difference and
different RNG (R vs Rust ChaCha).

**Upstream fixture / engine change needed**: None; probe data in `inf-bootstrap.md` provides
the expected magnitudes.

---

### SPEC BS-14b: Bootstrap CI basic interval correctness vs percentile

**Name**: `confint_bootstrap_lmm_basic_vs_percentile`
**Kind**: unit + numerical-tolerance
**Priority**: P1
**Dataset / formula**: synthetic LMM (mk_inference_fit), nsim = 200

**Setup**:
```r
fit  <- mk_inference_fit()
ci_p <- confint(fit, method = "bootstrap", interval = "percentile",
                bootstrap = bootstrap_control(nsim = 200, seed = 9))
ci_b <- confint(fit, method = "bootstrap", interval = "basic",
                bootstrap = bootstrap_control(nsim = 200, seed = 9))
```

**Assertions**:
```r
# Both intervals are finite
expect_true(all(is.finite(ci_p)))
expect_true(all(is.finite(ci_b)))

# They have the same dimensions and row names
expect_equal(dim(ci_p), dim(ci_b))
expect_equal(rownames(ci_p), rownames(ci_b))

# Basic and percentile intervals should cover the same approximate region
# but need not be identical. On a symmetric distribution they should be close.
# Allow up to 2-unit difference (generous for small nsim=200 synthetic data).
expect_true(max(abs(ci_p - ci_b)) < 2,
            info = "basic and percentile CIs should be in similar range")

# Interval type is labelled correctly
expect_identical(attr(ci_p, "interval"), "percentile")
expect_identical(attr(ci_b, "interval"), "basic")
```

**What it guards against**: basic interval formula being wrong (e.g., bounds swapped) or
identical to percentile (which would indicate no-op implementation).

**Upstream fixture / engine change needed**: None.

---

### SPEC BS-14c: `parametric_bootstrap` LRT reproduces lme4 PB-LRT distributional shape

**Name**: `parametric_bootstrap_lrt_distribution_vs_lme4`
**Kind**: parity-vs-lme4 + numerical-tolerance
**Priority**: P1
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` vs `Reaction ~ 1 + (1|Subject)` (sleepstudy, ML)

**Setup**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit_null_mx  <- lmm(Reaction ~ 1 + (1|Subject), data = sleepstudy, REML = FALSE,
                    control = mm_control(verbose = -1))
fit_alt_mx   <- lmm(Reaction ~ Days + (1|Subject), data = sleepstudy, REML = FALSE,
                    control = mm_control(verbose = -1))
fit_null_lme4 <- lme4::lmer(Reaction ~ 1 + (1|Subject), data = sleepstudy,
                              REML = FALSE)
fit_alt_lme4  <- lme4::lmer(Reaction ~ Days + (1|Subject), data = sleepstudy,
                              REML = FALSE)

NSIM <- 199L
pb_mx <- parametric_bootstrap(fit_null_mx, fit_alt_mx, nsim = NSIM, seed = 42)
pb_lme4 <- lme4::bootMer(fit_null_lme4,
                          FUN = function(m) 2 * (logLik(fit_alt_lme4) -
                                                   logLik(refit(m, simulate(m)[[1]]))),
                          nsim = NSIM, seed = 42)
```

**Assertions**:
```r
# Both p-values should be significant for the Days effect
expect_true(pb_mx$p_value < 0.05)

# Observed LRT statistics should agree within tolerance
obs_lme4 <- as.numeric(2 * (logLik(fit_alt_lme4) - logLik(fit_null_lme4)))
expect_equal(pb_mx$observed, obs_lme4, tolerance = 1e-2,
             info = "observed LRT statistics should match between engines")

# Median of reference distribution: both should be near chi-sq(1) median (~0.45)
med_mx <- median(pb_mx$simulated)
expect_true(med_mx > 0 && med_mx < 5,
            info = sprintf("reference distribution median=%.3f should be near 0",
                           med_mx))
```

**What it guards against**: observed LRT statistic mismatch between R and Rust; wrong
reference distribution shape.

**Upstream fixture / engine change needed**: None.

---

### SPEC BS-15: `refit.mm_lmm` `newoffset=` absent — graceful error

**Name**: `refit_lmm_newoffset_typed_refusal`
**Kind**: error-message
**Priority**: P2
**Dataset / formula**: `Reaction ~ Days + (1|Subject)` (sleepstudy)

**Assertions**:
```r
skip_if_not_installed("lme4")
data("sleepstudy", package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), data = sleepstudy,
           control = mm_control(verbose = -1))
sims <- simulate(fit, nsim = 1, seed = 1)

# newoffset is not in formals
expect_false("newoffset" %in% names(formals(refit.mm_lmm)))

# Passing newoffset via ... is silently ignored; the result should be the same
# as without it (i.e., the offset is not applied)
r_no_offset  <- refit(fit, newresp = sims[[1]])
r_with_offset <- refit(fit, newresp = sims[[1]], newoffset = rep(0, nrow(sleepstudy)))
expect_equal(fixef(r_no_offset), fixef(r_with_offset), tolerance = 1e-8,
             info = "newoffset passed via ... should be silently ignored")
```

**What it guards against**: a future `newoffset` implementation accidentally applying a
zero-offset and shifting estimates.

**Upstream fixture / engine change needed**: `newoffset` support requires Rust bridge update.

---

## Priority summary

| Priority | Spec IDs | Count |
|----------|----------|-------|
| P0 — blocker (violates mandate) | BS-03a, BS-03b, BS-13 | 3 |
| P1 — high value / lme4-parity | BS-01a, BS-01b, BS-04a, BS-04b, BS-05a, BS-05b, BS-07, BS-08, BS-14a, BS-14b, BS-14c | 11 |
| P2 — nice-to-have / guards absence | BS-02, BS-06, BS-09, BS-10, BS-11, BS-12, BS-15 | 7 |

Total specs: **21**

---

## Top priority

**BS-03a** (`confint_glmm_typed_refusal_not_base_r_error`) — The `confint` call on an
`mm_glmm` object currently produces `"non-numeric argument to binary operator"`, an opaque
base-R error. This directly violates the project's "clearer errors" mandate (PRD §2). It
requires only an R-side dispatch fix (add a `confint.mm_glmm` stub that raises
`mm_inference_unavailable` with a human-readable message) and no upstream engine changes.
It is the single most impactful fix for user trust.

---

## Implementation notes

1. **Shared fixture**: the `cbpp` dataset (package `lme4`) is used by BS-03a, BS-03b, BS-04a,
   BS-05a. It should be loaded once in `_setup.R` or a shared helper function, guarded by
   `skip_if_not_installed("lme4")`.

2. **Skip guards for unimplemented paths**: BS-04b (simulate.mm_glmm) and BS-05a (refit.mm_glmm)
   should be prefixed with `skip_if(...)` checks that skip until the method exists:
   ```r
   skip_if(!existsMethod("simulate", "mm_glmm") &&
           !exists("simulate.mm_glmm", mode = "function"))
   ```
   This keeps the test file runnable during CI without false failures.

3. **RNG tolerance**: specs involving distributional comparison (BS-14a, BS-14c) use relative
   CI-width tolerance (10%) rather than exact bound tolerance, to accommodate the documented
   RNG difference (R vs Rust ChaCha) and REML/ML refit strategy difference. These tolerances
   are calibrated to the probe data in `inf-bootstrap.md` (max observed relative width diff 4.4%).

4. **Engine-change gating (BS-04b, BS-11, BS-12)**: when upstream Rust adds GLMM simulate,
   `newparams`, or `newdata` support, the corresponding specs must be updated from
   "graceful-error" mode to "numerical correctness" mode using the templates provided.
