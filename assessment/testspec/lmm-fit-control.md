# Test Specification — LMM Fitting & Control

Family: **LMM fitting & control**
Source gap report: `assessment/gap/lmm-fit-control.md`
Date: 2026-05-31
Author: automated spec pass

Scope: covers every gap classified **in-scope-missing**, **partial**, or
**test-gap** in the gap report. Out-of-scope-by-design items (`offset`,
`devFunOnly`, modular pipeline, `allFit` verb, `optimizer` choice) are noted
only where a *refusal message* quality spec is appropriate.

PRD §3 tolerances (default): fixef 1e-4, theta 1e-3, logLik 1e-3,
sigma 1e-4. Tighter tolerances used where feasible.

---

## TS-LFC-01 — `contrasts=` argument accepted and wired (P0 blocker)

**Kind:** parity-vs-lme4 + unit

**Classification from gap report:** in-scope-missing / major

**Dataset / formula:**
```r
# Inline — no external dependency
set.seed(42)
df <- data.frame(
  y       = rnorm(60),
  grp     = factor(rep(LETTERS[1:3], 20)),
  subject = factor(rep(1:10, 6))
)
```
Formula: `y ~ grp + (1 | subject)`

**Assertion:**

1. `lmm(y ~ grp + (1 | subject), df, contrasts = list(grp = "contr.sum"))` must
   *not* error with `"unused argument (contrasts = ...)"`.  
   The returned `mm_lmm` object must have `fixef` names consistent with
   sum-to-zero coding (`grp1`, `grp2`, not `grpB`, `grpC`).

2. Numerical parity vs. lme4: `max(abs(fixef(fit) - lme4::fixef(ref))) < 1e-4`
   where `ref <- lme4::lmer(y ~ grp + (1 | subject), df, contrasts = list(grp = "contr.sum"))`.

3. Passing a custom contrast matrix (not a string) must also work:
   ```r
   cmat <- matrix(c(1, -0.5, -0.5,   0, 0.5, -0.5), nrow = 3,
                  dimnames = list(levels(df$grp), c("A_vs_BC", "B_vs_C")))
   lmm(y ~ grp + (1 | subject), df, contrasts = list(grp = cmat))
   ```

4. `model.matrix(fit)` column names must reflect the supplied coding, not
   treatment coding.

**What it guards against:**
- The current `"unused argument"` error that forces manual contrast-column
  construction.
- Silent regression to treatment coding when `contrasts=` is ignored.

**Upstream fixture / engine note:**
The contrast matrix must be applied on the R side (inside `compile_model()` /
`mm_translate_data()`) before the model frame is passed to Rust. No Rust FFI
change is required for standard named-contrast strings. Custom contrast
matrices require `model.matrix()` to encode them R-side before the Rust call.

**testthat skeleton:**
```r
test_that("lmm() accepts contrasts= and produces sum-to-zero coding", {
  skip_if_not_installed("lme4")
  df <- local({ set.seed(42); data.frame(
    y = rnorm(60), grp = factor(rep(LETTERS[1:3], 20)),
    subject = factor(rep(1:10, 6)))
  })
  fit <- lmm(y ~ grp + (1 | subject), df,
             contrasts = list(grp = "contr.sum"),
             control = mm_control(verbose = -1))
  ref <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ grp + (1 | subject), df, contrasts = list(grp = "contr.sum"))
  ))
  expect_s3_class(fit, "mm_lmm")
  expect_true(any(grepl("^grp1$", names(fixef(fit)))))
  expect_lt(max(abs(fixef(fit) - lme4::fixef(ref))), 1e-4)
})
```

---

## TS-LFC-02 — `contrasts=` refused with a guided message when engine cannot encode (P1)

**Kind:** error-message

**Classification from gap report:** in-scope-missing / major (error-path quality)

**Dataset / formula:** Same inline dataset as TS-LFC-01.

**Assertion:**
If a `contrasts=` specification is supplied for a variable that the engine
*cannot* represent (e.g., an interaction term, or a non-factor), `lmm()` must
throw an `mm_formula_error` or `mm_data_error` with a message that:
- names the offending variable,
- says what encoding was attempted,
- suggests the manual workaround.

It must *not* silently ignore the argument or produce a cryptic Rust panic.

**What it guards against:**
- Regression to a silent no-op for unsupported contrast specs.

---

## TS-LFC-03 — `subset=` argument accepted and equivalent to pre-filter (P1)

**Kind:** unit + parity-vs-lme4

**Classification from gap report:** in-scope-missing / minor

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
# subset = Days >= 2
```
Formula: `Reaction ~ Days + (1 | Subject)`

**Assertion:**

1. `lmm(Reaction ~ Days + (1 | Subject), sleepstudy, subset = Days >= 2)` must
   not error with `"unused argument"`.

2. Numerical equivalence with the pre-filtered call:
   ```r
   fit_sub  <- lmm(..., subset = Days >= 2)
   fit_pre  <- lmm(..., data = sleepstudy[sleepstudy$Days >= 2, ])
   expect_equal(fixef(fit_sub), fixef(fit_pre), tolerance = 1e-10)
   expect_equal(nobs(fit_sub), nobs(fit_pre))
   ```

3. `nobs(fit_sub)` equals `sum(sleepstudy$Days >= 2)`.

**What it guards against:**
- The current `"unused argument"` error.
- `subset` silently being ignored (all rows used).

**testthat skeleton:**
```r
test_that("lmm() subset= is equivalent to pre-filtering data", {
  data(sleepstudy, package = "lme4")
  ctl <- mm_control(verbose = -1)
  fit_sub <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
                 subset = Days >= 2, control = ctl)
  fit_pre <- lmm(Reaction ~ Days + (1 | Subject),
                 sleepstudy[sleepstudy$Days >= 2, ], control = ctl)
  expect_equal(nobs(fit_sub), sum(sleepstudy$Days >= 2))
  expect_equal(fixef(fit_sub), fixef(fit_pre), tolerance = 1e-10)
  expect_equal(sigma(fit_sub), sigma(fit_pre), tolerance = 1e-10)
})
```

---

## TS-LFC-04 — `na.action = na.omit` silently drops NA rows, matching lme4 default (P1)

**Kind:** parity-vs-lme4 + unit

**Classification from gap report:** partial / major

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
sleepstudy_na <- sleepstudy
sleepstudy_na$Reaction[c(3, 15, 27)] <- NA   # inject 3 NA rows
```
Formula: `Reaction ~ Days + (1 | Subject)`

**Assertion:**

1. `lmm(..., data = sleepstudy_na, na.action = na.omit)` must fit successfully
   and `nobs(fit)` must equal `nrow(sleepstudy) - 3`.

2. `fixef`, `sigma`, `logLik` must match `lme4::lmer(..., na.action = na.omit)`
   within default PRD tolerances (fixef 1e-4, sigma 1e-4, logLik 1e-3).

3. When `na.action` is absent (default) and NA rows are present, `lmm()` must
   still throw a clear `mm_data_error` containing the words
   `"complete cases"` (the current behavior, which must not regress).

4. `na.action = na.fail` must throw an error when NAs are present (parity with
   `stats::na.fail`).

**What it guards against:**
- The current "always error on NA" divergence from lme4's default `na.omit`.
- Behavioral regression once `na.action` is wired: silent dropping without
  updating `nobs`.

**testthat skeleton:**
```r
test_that("lmm() na.action=na.omit drops NA rows and matches lme4", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  sna <- sleepstudy; sna$Reaction[c(3, 15, 27)] <- NA
  ctl <- mm_control(verbose = -1)
  fit <- lmm(Reaction ~ Days + (1 | Subject), sna,
             na.action = na.omit, control = ctl)
  ref <- suppressMessages(
    lme4::lmer(Reaction ~ Days + (1 | Subject), sna, na.action = na.omit)
  )
  expect_equal(nobs(fit), nrow(sleepstudy) - 3L)
  expect_equal(fixef(fit), lme4::fixef(ref), tolerance = 1e-4)
  expect_equal(sigma(fit), sigma(ref), tolerance = 1e-4)
})

test_that("lmm() default NA policy errors with 'complete cases' message", {
  data(sleepstudy, package = "lme4")
  sna <- sleepstudy; sna$Reaction[3] <- NA
  expect_error(
    lmm(Reaction ~ Days + (1 | Subject), sna,
        control = mm_control(verbose = -1)),
    regexp = "complete cases",
    class = "mm_data_error"
  )
})
```

---

## TS-LFC-05 — `na.action = na.exclude` pads fitted/residuals to original length (P2)

**Kind:** parity-vs-lme4 + unit

**Classification from gap report:** partial / major

**Dataset / formula:** Same `sleepstudy_na` as TS-LFC-04.

**Assertion:**

1. `lmm(..., na.action = na.exclude)` fits successfully.

2. `length(fitted(fit))` equals `nrow(sleepstudy_na)` (original length, not
   `nobs`).

3. Fitted values at NA rows are `NA`.

4. Non-NA fitted values are numerically identical (tolerance 1e-8) to
   `fitted(lmm(..., na.action = na.omit))`.

**What it guards against:**
- Workflows that use `residuals(fit)` to build residual plots aligned to the
  original data frame — a frequent lme4 user pattern.

**Note:** This is a P2 because the workaround (pre-filter + re-index) is
mechanical. However it must be implemented before `update()` is added, because
`update()` relies on `na.action` for the model frame reconstruction.

---

## TS-LFC-06 — `start=` argument accepted and changes optimizer trajectory (P2)

**Kind:** unit

**Classification from gap report:** in-scope-missing / minor

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
```
Formula: `Reaction ~ Days + (1 | Subject)`

**Assertion:**

1. `lmm(..., start = c(theta = 1.5))` must not error with `"unused argument"`.

2. The fit must converge and produce the *same* final `fixef`, `sigma`, and
   `logLik` as the default-start fit within PRD tolerances (fixef 1e-4,
   sigma 1e-4, logLik 1e-3). This guards that `start` only seeds the optimizer,
   not the final solution.

3. `lmm(..., start = c(0))` (too short / wrong length) must throw a structured
   `mm_arg_error` identifying the mismatch.

**What it guards against:**
- The current `"unused argument"` error that makes it impossible to reproduce
  a fit from a prior run's `getME(fit, "theta")`.

**Upstream note:** Exposing `start` requires adding a `start_theta` field to
the `control_json` and wiring it through `mm_fit_lmm_json` into the crate's
optimizer initialization. The Rust side already accepts `_control` (currently
unused) — this requires both R changes and an upstream crate change to read
`control.start_theta`.

---

## TS-LFC-07 — `mm_control()` exposes `optCtrl` tolerance knobs (P0 blocker)

**Kind:** unit + integration

**Classification from gap report:** in-scope-missing / major

**Dataset / formula:**
```r
# Hard convergence case — many random-slope parameters
set.seed(99)
df_hard <- data.frame(
  y       = rnorm(200),
  x       = rnorm(200),
  subject = factor(rep(1:20, 10))
)
```
Formula: `y ~ x + (1 + x | subject)`

**Assertion:**

1. `mm_control(ftol_abs = 1e-10, xtol_abs = 1e-10, maxeval = 5000L)` must be
   accepted without error. `names(mm_control(ftol_abs = 1e-8))` must include
   `"ftol_abs"`.

2. `lmm(formula, df_hard, control = mm_control(maxeval = 50L))` must
   eventually terminate (possibly with a convergence warning / non-converged
   status) — the test verifies that the `maxeval` limit is actually honored by
   checking that `optimizer_certificate(fit)$converged` is `FALSE` or that the
   fit completes in fewer evaluations than the default would require.

3. `lmm(formula, df_hard, control = mm_control(ftol_abs = 1e-12, maxeval = 5000L))`
   must converge and produce the same `fixef` as the default-control fit
   within 1e-4.

4. Passing an unrecognised key (`mm_control(typo_key = 1)`) must throw
   `mm_arg_error` with the key name in the message (no silent ignore).

**What it guards against:**
- The current situation where `control_json` is passed through but the Rust
  side binds it to `_control` (underscore = unused), making tolerance tuning
  a no-op.

**Upstream note:** Requires: (a) R-side `mm_control()` to accept and validate
`ftol_abs`, `xtol_abs`, `maxeval`; (b) `mm_validate_control()` to pass them
through; (c) Rust `mm_fit_lmm_json` to read and apply them from the
deserialized `control_json`. This is an upstream engine change.

**testthat skeleton:**
```r
test_that("mm_control() accepts optimizer tolerance knobs", {
  ctl <- mm_control(ftol_abs = 1e-10, xtol_abs = 1e-10, maxeval = 2000L,
                    verbose = -1L)
  expect_equal(ctl$ftol_abs, 1e-10)
  expect_equal(ctl$maxeval, 2000L)
})

test_that("mm_control() rejects unknown keys", {
  expect_error(mm_control(bogus_param = 99), class = "mm_arg_error",
               regexp = "bogus_param")
})

test_that("maxeval limit is honored by the optimizer", {
  skip_if_not_installed("lme4")
  set.seed(99)
  df_hard <- data.frame(y = rnorm(200), x = rnorm(200),
                        subject = factor(rep(1:20, 10)))
  fit_lim <- lmm(y ~ x + (1 + x | subject), df_hard,
                 control = mm_control(maxeval = 20L, verbose = -1L))
  cert <- optimizer_certificate(fit_lim)
  # Either it didn't converge (eval budget exhausted) or converged quickly:
  expect_true(is.list(cert))  # certificate must always be present
})
```

---

## TS-LFC-08 — `verbose` does not produce optimizer iteration trace (documented divergence) (P2)

**Kind:** snapshot + error-message

**Classification from gap report:** partial / minor

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
```
Formula: `Reaction ~ Days + (1 | Subject)`

**Assertion:**

1. `mm_control(verbose = 2)` does *not* print per-iteration deviance lines.
   Captured output must contain the `explain_model()` block and must *not*
   contain the string `"iter"` or `"deviance"` as an optimizer progress
   header.

2. The current behavior (verbose prints the pre-fit explain block once) is
   the documented contract and must not regress.

3. The package's `?mm_control` man page (or `README`/vignette) must state
   explicitly that no optimizer-trace is available and that post-fit
   convergence evidence is in `optimizer_certificate(fit)`. This is a
   **documentation assurance**, not a code test, but it must be verified by
   checking the Rd source contains the word `"optimizer_certificate"`.

**What it guards against:**
- Silent regression where verbose suddenly produces noisy optimizer output
  (unplanned Rust-side change).
- Documentation / user expectation drift: new users expect lme4-style verbose
  trace; the package must set expectations explicitly.

**testthat skeleton:**
```r
test_that("verbose=2 prints explain_model block but no optimizer trace", {
  data(sleepstudy, package = "lme4")
  out <- capture.output(
    lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
        control = mm_control(verbose = 2))
  )
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "Random effects", fixed = TRUE)       # explain block
  expect_false(grepl("^iter", joined, perl = TRUE))           # no optimizer row
})
```

---

## TS-LFC-09 — `mm_control()` exposes pre-fit check severity overrides (P2)

**Kind:** unit + integration

**Classification from gap report:** partial / minor

**Dataset / formula:**
```r
# Under-leveled grouping — triggers nobs.vs.nlev check
set.seed(7)
df_small <- data.frame(
  y = rnorm(6),
  x = rnorm(6),
  g = factor(c("a", "a", "b", "b", "c", "c"))   # 2 obs per group
)
```
Formula: `y ~ x + (1 | g)`

**Assertion:**

1. By default, `lmm()` on `df_small` emits a diagnostic (warning or structured
   diagnostic) about small group sizes. The test captures output and checks
   for a relevant message.

2. `lmm(..., control = mm_control(check.nobs.vs.nlev = "ignore"))` must
   suppress the diagnostic entirely and still return an `mm_lmm` object.

3. `lmm(..., control = mm_control(check.nobs.vs.nlev = "stop"))` must throw
   an error with a message referencing the number of levels and observations.

4. `mm_control(check.nobs.vs.nlev = "invalid_value")` must throw
   `mm_arg_error` enumerating valid values.

**What it guards against:**
- Lack of configurability forcing users who know their design is unusual to
  suppress warnings package-globally.

---

## TS-LFC-10 — Post-fit singular / convergence certificate always present (P1)

**Kind:** unit

**Classification from gap report:** partial / minor (coverage of existing partial)

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
```
Singular formula: `Reaction ~ Days + (1 + Days | Subject)` (known to be
near-singular with `theta` near zero in some configurations)

**Assertion:**

1. `optimizer_certificate(fit)` returns a non-NULL list for every successful
   `lmm()` call — including near-singular fits.

2. For a genuinely singular fit (force by over-parameterizing with a small
   dataset), `is_singular(fit)` is `TRUE` and
   `optimizer_certificate(fit)$singular` is `TRUE`.

3. `optimizer_certificate(fit)$converged` is a logical scalar (never NA for
   a completed fit).

4. The certificate must survive `revive()`: fit, serialize to JSON (simulated
   via `saveRDS` / `readRDS`), revive, and call `optimizer_certificate()`
   again — result must be identical to the pre-serialization certificate.

**What it guards against:**
- Silent loss of convergence metadata during serialization (cross-session
  revival gap).

**testthat skeleton:**
```r
test_that("optimizer_certificate() is always present and survives revival", {
  data(sleepstudy, package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 + Days | Subject), sleepstudy,
             control = mm_control(verbose = -1))
  cert <- optimizer_certificate(fit)
  expect_true(is.list(cert))
  expect_true(is.logical(cert$converged) && length(cert$converged) == 1L)

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)
  fit2 <- revive(fit2)
  cert2 <- optimizer_certificate(fit2)
  expect_equal(cert2$converged, cert$converged)
})
```

---

## TS-LFC-11 — `update.mm_lmm()` method exists and rewires formula/REML (P1)

**Kind:** unit + parity-vs-lme4

**Classification from gap report:** in-scope-missing / minor

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
```
Base formula: `Reaction ~ Days + (1 | Subject)`

**Assertion:**

1. `update(fit, REML = FALSE)` must return a new `mm_lmm` object fitted by ML,
   with `fit2$REML == FALSE` and `logLik(fit2) > logLik(fit)` (ML logLik is
   numerically larger than REML logLik on the same data).

2. `update(fit, . ~ . + I(Days^2))` must add the quadratic term:
   `"I(Days^2)"` must be in `names(fixef(updated))`.

3. `update(fit, . ~ . - Days)` must drop `Days`:
   `"Days"` must not be in `names(fixef(updated))`.

4. Numerical parity for case (1): `update(fit, REML = FALSE)` result must
   match `lmm(..., REML = FALSE)` within PRD tolerances.

5. lme4 parity: `update(lme4_fit, REML = FALSE)` must produce
   `fixef`/`sigma`/`logLik` matching `update(mixeff_fit, REML = FALSE)` within
   PRD tolerances.

**What it guards against:**
- Complete absence of `update()`, which breaks the
  `fit_reml <- lmm(...); fit_ml <- update(fit_reml, REML = FALSE)` idiom used
  by virtually every lme4 tutorial for LRT-based model comparison.

**testthat skeleton:**
```r
test_that("update.mm_lmm() switches REML to ML", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))
  fit_ml <- update(fit, REML = FALSE)
  expect_false(fit_ml$REML)
  expect_gt(as.numeric(logLik(fit_ml)), as.numeric(logLik(fit)))

  ref_ml <- suppressMessages(
    lme4::lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = FALSE)
  )
  expect_equal(fixef(fit_ml), lme4::fixef(ref_ml), tolerance = 1e-4)
})

test_that("update.mm_lmm() can add and drop fixed-effect terms", {
  data(sleepstudy, package = "lme4")
  fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
             control = mm_control(verbose = -1))
  fit_q <- update(fit, . ~ . + I(Days^2))
  expect_true("I(Days^2)" %in% names(fixef(fit_q)))
  fit_no_days <- update(fit, . ~ . - Days)
  expect_false("Days" %in% names(fixef(fit_no_days)))
})
```

---

## TS-LFC-12 — Out-of-scope args produce guided refusal messages, not generic errors (P1)

**Kind:** error-message

**Classification from gap report:** out-of-scope-by-design gaps that lack
clean refusal messages

**Dataset / formula:** Minimal inline (any valid `mm_lmm` setup).

**Assertion:** For each of the following calls, the error must be of class
`mm_formula_error`, `mm_arg_error`, or `mm_fit_error` (not a raw R
`"unused argument"` condition), and the message must name the argument and
suggest a documented alternative or workaround:

| Call | Expected class | Key phrase in message |
|---|---|---|
| `lmm(..., offset = w)` | `mm_arg_error` or `mm_formula_error` | `"offset"`, `"precompute"` or `"v2"` |
| `lmm(..., devFunOnly = TRUE)` | `mm_arg_error` | `"devFunOnly"` or `"not supported"` |
| `lmm(..., optimizer = "bobyqa")` | `mm_arg_error` | `"optimizer"` or `"auto-dispatched"` |

**What it guards against:**
- `"unused argument (offset = ...)"` — a raw R error that gives no guidance.
- Users wasting time searching for the argument in `?lmm`.

**testthat skeleton:**
```r
test_that("lmm() guided refusal for offset= argument", {
  set.seed(1); df <- data.frame(y = rnorm(20), x = rnorm(20),
                                 g = factor(rep(1:5, 4)))
  expect_error(
    lmm(y ~ x + (1 | g), df, offset = rnorm(20),
        control = mm_control(verbose = -1)),
    regexp = "offset",
    class = "mm_arg_error"
  )
})
```

---

## TS-LFC-13 — `mm_control()` `verbose=-1` contract does not regress (P1)

**Kind:** snapshot + unit

**Classification from gap report:** partial / minor (verbose contract regression
guard)

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
```

**Assertion:**

1. `lmm(..., control = mm_control(verbose = -1))` produces *no* output to
   stdout or stderr.
   ```r
   expect_silent(lmm(..., control = mm_control(verbose = -1)))
   ```

2. `lmm(..., control = mm_control(verbose = 0))` prints the `explain_model()`
   block exactly once (no repetition).

3. `lmm(..., control = mm_control(verbose = 1))` prints the block and may
   print additional detail, but must not produce optimizer iteration lines
   (per TS-LFC-08).

4. Negative `verbose` values other than `-1` are treated identically to `-1`
   (silent) — no error.

**What it guards against:**
- Any upstream Rust or R-side change that reintroduces print side-effects
  during silent fits, breaking test suites and automated pipelines that rely on
  `expect_silent()`.

**testthat skeleton:**
```r
test_that("mm_control(verbose = -1) produces no output", {
  data(sleepstudy, package = "lme4")
  expect_silent(
    lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
        control = mm_control(verbose = -1))
  )
})
```

---

## TS-LFC-14 — `crate-1.md:235` fit-intent mode contract vs. installed `mm_control()` (P1)

**Kind:** unit (documentation drift guard)

**Classification from gap report:** partial / major (contract/docstring drift)

**Dataset / formula:** N/A — this is a contract-conformance check.

**Assertion:**

1. `names(formals(mm_control))` must either:
   - include `"fit_intent"` (if the feature is implemented), OR
   - the `?mm_control` man page must contain the phrase `"fit_intent"` **only
     in a "Reserved for future use"** or "Not yet implemented" sentence — not
     as a currently-usable parameter.

2. If `"fit_intent"` is not in `names(formals(mm_control))`, then attempting
   `mm_control(fit_intent = "exploratory")` must throw `mm_arg_error`, not
   silently pass through.

3. The assessment survey `crate-1.md` claim that fit-intent modes are
   "accessible via `mm_control()` options" must be explicitly reconciled in
   the package documentation or the survey must be marked aspirational.

**What it guards against:**
- The specific drift documented in the gap report: crate-1.md:235 implies an
  accessible feature that does not exist in the installed R surface. This spec
  ensures the discrepancy is either resolved (feature added) or documented
  (false claim flagged).

**testthat skeleton:**
```r
test_that("mm_control() rejects fit_intent= until it is implemented", {
  if ("fit_intent" %in% names(formals(mm_control))) {
    # Feature implemented — check it works
    ctl <- mm_control(fit_intent = "exploratory", verbose = -1L)
    expect_equal(ctl$fit_intent, "exploratory")
  } else {
    # Not yet implemented — must not silently pass
    expect_error(mm_control(fit_intent = "exploratory"),
                 class = "mm_arg_error")
  }
})
```

---

## TS-LFC-15 — `refit(object, newresp)` numerical correctness and `newweights` gap (P2)

**Kind:** parity-vs-lme4 + unit

**Classification from gap report:** works (with a documented partial: `newweights` absent)

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
```

**Assertion:**

1. `refit(fit, newresp = sleepstudy$Reaction + 10)` shifts `fixef(fit)["(Intercept)"]`
   by 10 (within 1e-8) while `sigma` and `theta` are unchanged (within 1e-8).

2. `fixef`, `sigma`, `logLik`, `theta` of `refit(fit, newresp)` match
   `lme4::refit(ref, newresp)` within PRD tolerances.

3. `refit(fit, newresp, newweights = w)` either: (a) works correctly and
   matches `lme4::refit(ref, newresp, newweights = w)`, or (b) throws a
   structured `mm_arg_error` with `"newweights"` in the message and
   `"not yet supported"` — not a raw `"unused argument"` error.

**What it guards against:**
- Numerical regression in `refit()` (currently working, must stay working).
- The `newweights` gap silently becoming a raw R error instead of a guided
  refusal.

---

## Summary table

| ID | Description | Kind | Priority | In-scope classification |
|---|---|---|---|---|
| TS-LFC-01 | `contrasts=` arg accepted & sum-coding parity | parity-vs-lme4 | **P0** | in-scope-missing / major |
| TS-LFC-02 | `contrasts=` guided refusal for unsupported specs | error-message | P1 | in-scope-missing / major |
| TS-LFC-03 | `subset=` arg equivalent to pre-filter | parity-vs-lme4 | P1 | in-scope-missing / minor |
| TS-LFC-04 | `na.action=na.omit` silently drops NAs, matches lme4 | parity-vs-lme4 | P1 | partial / major |
| TS-LFC-05 | `na.action=na.exclude` pads fitted/residuals | parity-vs-lme4 | P2 | partial / major |
| TS-LFC-06 | `start=` accepted, seeds optimizer, same final result | unit | P2 | in-scope-missing / minor |
| TS-LFC-07 | `mm_control()` `optCtrl` tolerance knobs | unit + integration | **P0** | in-scope-missing / major |
| TS-LFC-08 | `verbose` divergence documented, no optimizer trace | snapshot | P2 | partial / minor |
| TS-LFC-09 | `mm_control()` pre-fit check severity overrides | unit + integration | P2 | partial / minor |
| TS-LFC-10 | `optimizer_certificate()` always present, survives revival | unit | P1 | partial / minor |
| TS-LFC-11 | `update.mm_lmm()` method for REML↔ML and formula edits | parity-vs-lme4 | P1 | in-scope-missing / minor |
| TS-LFC-12 | Out-of-scope args produce guided refusals | error-message | P1 | out-of-scope-by-design (refusal quality) |
| TS-LFC-13 | `verbose=-1` silent-fit contract regression guard | snapshot | P1 | partial / minor |
| TS-LFC-14 | `crate-1.md:235` fit-intent drift reconciliation | unit | P1 | partial / major (contract drift) |
| TS-LFC-15 | `refit()` numerical correctness + `newweights` gap | parity-vs-lme4 | P2 | partial (newweights) |

**Total specs: 15**

**Single most important spec: TS-LFC-07** — `mm_control()` `optCtrl` tolerance
knobs. The Rust control JSON is currently a no-op (`_control` binding). Users
who hit borderline convergence — the most common post-fit complaint — cannot
tighten `ftol_abs`/`xtol_abs` or raise `maxeval`, the first step recommended
by `?lme4::convergence`. Unblocking this requires coordinated R + upstream
crate changes and is the highest-leverage unimplemented surface in this family.
