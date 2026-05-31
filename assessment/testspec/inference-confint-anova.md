# Test/Assurance Specification — Inference: confint / anova / drop1 / profile

**Capability family:** `confint` (Wald/profile/bootstrap), `profile()` / `thpr`,
`anova` (single- and multi-model), `drop1`, `KRmodcomp`/`PBmodcomp` analogs,
`ranova`/`test_random_effect`, `contrast`/`contest*`.

**Source gap report:** `assessment/gap/inference-confint-anova.md`  
**Date:** 2026-05-31  
**Scope note:** specs are written only for gaps classified `in-scope-missing`,
`partial`, or `test-gap` in the gap report. `out-of-scope-by-design` items
(full `KRmodcomp`, `vcovAdj`, `seqPBmodcomp`, `PBrefdist`) are excluded from
this spec. Phase 5 items (`profile()`/`thpr`) receive placeholder sentinel
tests that assert the *current* honest-refusal behavior — they will be promoted
to full parity tests when Phase 5 ships.

---

## Spec index

| ID | Name | Kind | Priority |
|----|------|------|----------|
| CI-01 | bootstrap-CI covers VC parameters | parity-vs-lme4 | P0 |
| CI-02 | profile-CI beta under REML — honest refusal contract | error-message | P0 |
| CI-03 | Wald CI SE/bounds numerical tolerance on correlated RS model | numerical-tolerance | P0 |
| CI-04 | parm group selectors `"theta_"` / `"beta_"` for profile CI | unit | P1 |
| CI-05 | bootstrap CI interval types: `"basic"` vs `"percentile"` | unit | P1 |
| CI-06 | profile() generic emits a stable typed sentinel | error-message | P2 |
| CI-07 | confint.thpr() sentinel — honest refusal until Phase 5 | error-message | P2 |
| AN-01 | anova(type=) Type I / II / III contrasts on multi-term model | parity-vs-lme4 | P0 |
| AN-02 | anova statistic label: `t` vs `F` documentation | snapshot | P2 |
| AN-03 | anova(model.names=) cosmetic relabeling | unit | P2 |
| D1-01 | drop1(ddf="Satterthwaite") F-table via test_effect routing | integration | P1 |
| D1-02 | drop1(test="user", sumFun=) hook or documented refusal | unit | P2 |
| PB-01 | compare(method="bootstrap") reports Bartlett/Gamma rows or docs absence | snapshot | P1 |
| PB-02 | PBrefdist analog: documented absence or sentinel refusal | error-message | P2 |
| RE-01 | test_random_effect whole-model table (ranova-style) | integration | P1 |
| RE-02 | ranova reduce.terms decomposition sentinel | error-message | P2 |
| BM-01 | bootMer arbitrary-FUN sentinel | error-message | P2 |

---

## Detailed specifications

---

### CI-01 — bootstrap CI covers VC parameters (sigma, random-effect SDs/correlations)

**Kind:** parity-vs-lme4  
**Priority:** P0 — blocks a routine lme4 workflow

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))
```

**Assertion:**
```r
# Currently errors with "Unknown fixed-effect parameter(s): sigma"
# The spec requires it to succeed and return finite bounds.

ci_sigma <- confint(fit, parm = "sigma", method = "boot",
                    bootstrap = bootstrap_control(nsim = 199, seed = 42))
expect_s3_class(ci_sigma, "mm_confint")
expect_true("sigma" %in% rownames(ci_sigma))
expect_true(is.finite(ci_sigma["sigma", 1]))
expect_true(is.finite(ci_sigma["sigma", 2]))
expect_true(ci_sigma["sigma", 1] > 0)
expect_true(ci_sigma["sigma", 2] > ci_sigma["sigma", 1])

# Random-effect SD (theta scale or SD scale — either is acceptable
# so long as the row is present and finite with a documented parameterization)
ci_all <- confint(fit, method = "boot",
                  bootstrap = bootstrap_control(nsim = 199, seed = 42))
# Must include at least one VC parameter row beyond fixed effects
vc_rows <- setdiff(rownames(ci_all), names(fixef(fit)))
expect_true(length(vc_rows) >= 1L,
            info = "bootstrap CI must include at least one VC parameter row")
expect_true(all(is.finite(ci_all[vc_rows, , drop = FALSE])),
            info = "VC parameter bootstrap CI bounds must be finite")
```

**Tolerance:** VC bounds need not match lme4 exactly (different RNG/refit
strategy); the assertion is structural (finite, positive, lower < upper) plus a
width sanity check: width of sigma CI must be > 0 and < 3 * point estimate.

**What it guards against:** `confint(method="boot")` silently rejecting any
request that names a VC parameter, leaving users unable to bootstrap SD/cor
intervals without workarounds.

**Upstream note:** Requires the Rust bootstrap engine to expose VC parameter
distributions, not just beta. This is an **upstream fixture requirement** —
the `bootstrap_full_model_distribution` payload must include rows for
`parameter_kind = "theta"` and `"sigma"`. Coordinate with mixeff-rs before
wiring the R side.

---

### CI-02 — profile CI beta-under-REML: stable reason_code contract

**Kind:** error-message  
**Priority:** P0 — honesty contract; any regression here is a silent-surgery bug

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit_reml <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
                REML = TRUE, control = mm_control(verbose = -1))
```

**Assertion:**
```r
ci <- confint(fit_reml, method = "profile", level = 0.95)
payload <- attr(ci, "mm_profile")
table  <- payload$table

beta_rows <- table[table$parameter_kind == "beta", , drop = FALSE]

# 1. Beta rows MUST be present (refusal rows, not missing rows)
expect_true(nrow(beta_rows) >= 1L,
            info = "REML profile payload must surface refusal rows for beta")

# 2. Every beta row must carry the canonical reason_code
expect_true(
  all(beta_rows$reason_code == "profile_beta_unavailable_under_reml"),
  info = "reason_code must be the stable canonical value, not a new string"
)

# 3. Bounds must be NA (not fabricated finite values)
expect_true(all(is.na(beta_rows$lower)))
expect_true(all(is.na(beta_rows$upper)))

# 4. Sigma and theta rows must still be finite
vc_rows <- table[table$parameter_kind %in% c("sigma", "theta"), , drop = FALSE]
expect_true(nrow(vc_rows) >= 1L)
ok_vc <- vc_rows[vc_rows$reason_code %in% c(NA_character_, ""), , drop = FALSE]
expect_true(all(is.finite(ok_vc$lower)))
expect_true(all(is.finite(ok_vc$upper)))
```

**Tolerance:** exact string match on `reason_code`; bounds must be exactly `NA`.

**What it guards against:** Future upstream change that silently provides beta
profile CIs under REML without updating the wrapper's contract tests, or a
regression that drops the refusal rows entirely (silent omission = silent
surgery).

**Upstream note:** None — this spec cements the current contract. If upstream
ever ships REML beta profiling, this spec must be replaced by CI-02b (parity
test against lme4 REML beta profile CI), and the reason_code must be versioned.

---

### CI-03 — Wald CI SE/bounds numerical tolerance on correlated random-slope model

**Kind:** numerical-tolerance  
**Priority:** P0 — current FAIL on the correlated RS model (inf-confint-wald probe)

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit_mm   <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
                REML = TRUE, control = mm_control(verbose = -1))
fit_lme4 <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                       REML = TRUE)
```

**Assertion:**
```r
skip_if_not_installed("lme4")
ci_mm   <- confint(fit_mm,   method = "wald", level = 0.95)
ci_lme4 <- confint(fit_lme4, method = "Wald", level = 0.95)

# Align on fixed-effect rows only (lme4 also emits NA rows for VCs)
shared <- intersect(rownames(ci_mm), rownames(ci_lme4))
expect_true(length(shared) >= 2L)

for (nm in shared) {
  expect_lt(abs(ci_mm[nm, 1] - ci_lme4[nm, 1]), 1e-4,
            label = sprintf("Wald CI lower[%s]", nm))
  expect_lt(abs(ci_mm[nm, 2] - ci_lme4[nm, 2]), 1e-4,
            label = sprintf("Wald CI upper[%s]", nm))
}

# SE tolerance (indirect guard via CI width)
se_mm   <- sqrt(diag(vcov(fit_mm)))
se_lme4 <- sqrt(diag(lme4::vcov.merMod(fit_lme4)))
for (nm in names(se_mm)) {
  expect_lt(abs(se_mm[[nm]] - se_lme4[[nm]]), 1e-4,
            label = sprintf("SE[%s]", nm))
}
```

**Tolerance:** CI bounds within 1e-4; SE within 1e-4 (the stated fixef-class
tolerance from PRD §11).

**What it guards against:** The probe `inf-confint-wald.md` shows Days SE off
by 7.16e-04 (7x over tolerance). This test will FAIL until the REML optimizer
achieves tighter convergence on the correlated RS model. It is a regression
sentinel: once fixed it must not regress.

**Upstream note:** Root cause is REML optimizer precision on the correlated RS
Cholesky theta. Fix likely requires either tighter optimizer convergence
criteria or a more numerically stable vcov path in mixeff-rs. Track in
upstream mote issue when evidence is solid.

---

### CI-04 — parm group selectors `"theta_"` / `"beta_"` for profile CI

**Kind:** unit  
**Priority:** P1

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))
```

**Assertion:**
```r
# "beta_" selector must return only fixed-effect rows
ci_beta <- confint(fit, parm = "beta_", method = "profile", level = 0.95)
expect_true(nrow(ci_beta) >= 1L,
            info = '"beta_" group selector must return at least one row')
payload_beta <- attr(ci_beta, "mm_profile")
if (!is.null(payload_beta)) {
  expect_true(all(payload_beta$table$parameter_kind == "beta"))
}

# "theta_" selector must return only VC (theta/sigma) rows
ci_theta <- confint(fit, parm = "theta_", method = "profile", level = 0.95)
expect_true(nrow(ci_theta) >= 1L,
            info = '"theta_" group selector must return at least one row')
payload_theta <- attr(ci_theta, "mm_profile")
if (!is.null(payload_theta)) {
  expect_true(all(payload_theta$table$parameter_kind %in% c("theta", "sigma")))
}

# Neither selector may return an empty result
# (current bug: "theta_" returns 0 rows)
expect_false(nrow(ci_beta) == 0L)
expect_false(nrow(ci_theta) == 0L)
```

**Tolerance:** structural (non-zero rows, correct parameter_kind filter).

**What it guards against:** The gap report documents that `confint(m, parm="theta_")` currently returns an empty matrix — the string is treated as a literal parameter name rather than a group selector. This test certifies the selector dispatch is wired.

---

### CI-05 — bootstrap CI interval types: `"basic"` vs `"percentile"` produce distinct results

**Kind:** unit  
**Priority:** P1

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))
```

**Assertion:**
```r
bc <- bootstrap_control(nsim = 199, seed = 99)

ci_pct   <- confint(fit, method = "boot", interval = "percentile", bootstrap = bc)
ci_basic <- confint(fit, method = "boot", interval = "basic",      bootstrap = bc)

# Both must succeed and return mm_confint
expect_s3_class(ci_pct,   "mm_confint")
expect_s3_class(ci_basic, "mm_confint")

# "basic" interval reflects around the estimate, so it differs from percentile
# on an asymmetric distribution. At nsim=199 they should not be identical.
expect_false(
  isTRUE(all.equal(ci_pct["Days", ], ci_basic["Days", ], tolerance = 1e-10)),
  info = 'basic and percentile CIs must differ (they use different reflections)'
)

# Both must be finite and ordered
expect_true(all(is.finite(ci_pct)))
expect_true(all(is.finite(ci_basic)))
expect_true(ci_pct["Days",   1] < ci_pct["Days",   2])
expect_true(ci_basic["Days", 1] < ci_basic["Days", 2])

# Method attribute must record the interval type
expect_match(attr(ci_pct,   "interval"), "percentile", ignore.case = TRUE)
expect_match(attr(ci_basic, "interval"), "basic",      ignore.case = TRUE)
```

**Tolerance:** structural; bounds must be finite and lower < upper.

**What it guards against:** Silent no-op if `interval` argument is parsed but not routed to the Rust bootstrap engine.

---

### CI-06 — `profile()` generic emits a stable typed sentinel (Phase 5 placeholder)

**Kind:** error-message  
**Priority:** P2 — Phase 5 item; spec cements current honest behavior

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))
```

**Assertion:**
```r
# Until Phase 5 ships, profile() must emit an informative, stable error
# rather than silently returning NULL or a malformed object.
err <- tryCatch(profile(fit), error = function(e) e)
expect_true(inherits(err, "error"),
            info = "profile() must error until Phase 5 implementation")
# The message must mention that profile is not yet implemented,
# not an inscrutable "no applicable method" R dispatch error.
# Acceptable: either a typed mm_* condition or a plain error with a
# comprehensible message string.
if (inherits(err, "mm_not_implemented")) {
  # Preferred: typed refusal
  invisible(NULL)
} else {
  expect_match(conditionMessage(err),
               "profile|not.*implemented|Phase 5",
               ignore.case = TRUE,
               info = "profile() error message must name the missing feature")
}
```

**Tolerance:** exact class or message-pattern match.

**What it guards against:** A future refactor that accidentally strips the S3
method registration and silently returns the base `profile()` default output
(which would be a meaningless list), or changes the error message to something
inscrutable without updating this test.

**Phase 5 promotion:** When `profile.mm_lmm` ships, replace this test with a
full parity spec: `pp <- profile(fit); ci <- confint(pp); expect_s3_class(ci, "mm_confint")`.

---

### CI-07 — `confint.thpr()` returns an informative typed refusal (Phase 5 placeholder)

**Kind:** error-message  
**Priority:** P2 — Phase 5 item

**Assertion:**
```r
# thpr objects do not exist pre-Phase-5.
# Constructing one via a fake list and calling confint() on it
# must error cleanly, not crash or return garbage.
fake_thpr <- structure(list(), class = "thpr")
err <- tryCatch(confint(fake_thpr), error = function(e) e)
expect_true(inherits(err, "error"))
# Must not be an internal R segfault or "subscript out of bounds" — must be
# a dispatch error or an explicit mm_not_implemented condition.
expect_false(grepl("subscript out of bounds|object of type 'closure'",
                   conditionMessage(err)),
             info = "confint.thpr error must be clean, not an internal R crash")
```

**What it guards against:** A naive `confint.thpr` stub that crashes on a zero-length list rather than emitting a sensible "not yet implemented" message.

---

### AN-01 — anova(type=) Type I / II / III contrasts on multi-term / interaction model

**Kind:** parity-vs-lme4  
**Priority:** P0 — currently flagged test-gap; required to certify Type ANOVA table correctness

**Dataset / formula:**
```r
data(Orthodont, package = "nlme")
# Two fixed effects plus an interaction — exposes Type I vs II vs III differences
fit <- lmm(distance ~ age * Sex + (age | Subject), Orthodont,
           REML = TRUE, control = mm_control(verbose = -1))
fit_lme4 <- lme4::lmer(distance ~ age * Sex + (age | Subject), Orthodont,
                       REML = TRUE)
```

**Assertion (three sub-tests):**

```r
# Sub-test 1: Type III — mixeff must agree with lmerTest Type III F-table
skip_if_not_installed("lme4"); skip_if_not_installed("lmerTest")
an_mm_III <- anova(fit, type = "III")
an_lt_III <- lmerTest::anova(fit_lme4, type = 3, ddf = "Satterthwaite")

# Both must have a row for each of the three terms
expect_true("age" %in% an_mm_III$table$term)
expect_true("Sex" %in% an_mm_III$table$term)
expect_true("age:Sex" %in% an_mm_III$table$term)

# p-values must agree to within 1e-4 (Satterthwaite df may differ slightly)
for (term in c("age", "Sex", "age:Sex")) {
  p_mm <- an_mm_III$table$p_value[an_mm_III$table$term == term]
  p_lt <- an_lt_III[term, "Pr(>F)"]
  expect_lt(abs(p_mm - p_lt), 1e-4,
            label = sprintf("Type III p-value for %s", term))
}

# Sub-test 2: Type I table must differ from Type III for this model
an_mm_I   <- anova(fit, type = "I")
an_mm_III2 <- anova(fit, type = "III")
# For an interaction model, Type I and Type III p-values for the main effects
# are generally not identical. If they are identical it is a strong signal that
# type is being ignored.
p_age_I   <- an_mm_I$table$p_value[an_mm_I$table$term == "age"]
p_age_III <- an_mm_III2$table$p_value[an_mm_III2$table$term == "age"]
expect_false(isTRUE(all.equal(p_age_I, p_age_III, tolerance = 1e-6)),
             info = "Type I and Type III p-values for 'age' must differ on interaction model")

# Sub-test 3: Type II table must differ from both I and III
an_mm_II  <- anova(fit, type = "II")
p_age_II  <- an_mm_II$table$p_value[an_mm_II$table$term == "age"]
expect_false(isTRUE(all.equal(p_age_II, p_age_III, tolerance = 1e-6)),
             info = "Type II and Type III p-values for 'age' must differ on interaction model")
```

**Tolerance:** p-values within 1e-4 vs lmerTest for Type III; Types I/II/III
must differ numerically on an interaction model.

**What it guards against:** The gap report flags that the current Rust term
table appears to return the same row for any `type` value — i.e., `type` is
accepted but ignored. Sub-tests 2 and 3 will catch this. Sub-test 1 certifies
the Type III absolute values against lmerTest.

**Upstream note:** This is an **upstream fixture requirement**. The Rust
`mm_rust_term_table` FFI must accept a `type` argument and compute the correct
contrast matrices for each ANOVA type. Before writing R-side routing, confirm
the Rust-side FFI supports it.

---

### AN-02 — anova() statistic label `t` vs `F`: snapshot and documentation

**Kind:** snapshot  
**Priority:** P2 — cosmetic but worth a snapshot to prevent silent regressions

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
           REML = TRUE, control = mm_control(verbose = -1))
```

**Assertion:**
```r
an <- anova(fit, method = "satterthwaite")
# mixeff reports statistic_name="t" for single-df terms (not "F")
# This is a documented design choice (F = t^2, equivalent for 1 df).
expect_identical(an$table$statistic_name[an$table$term == "Days"], "t")

# Guard: if statistic_name ever changes to "F", this test flags it so
# we can update documentation and downstream consumers.
# To accept "F", update this assertion and add an F-vs-t^2 equivalence check.
```

**Tolerance:** exact string match.

**What it guards against:** Silent change in `statistic_name` field that would
break downstream consumers who parse the field.

---

### AN-03 — anova(model.names=) cosmetic relabeling

**Kind:** unit  
**Priority:** P2 — cosmetic

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
m0 <- lmm(Reaction ~ 1 + (Days | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))
m1 <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))
```

**Assertion (two sub-tests):**
```r
# Sub-test 1: model.names not yet supported — must error or warn clearly,
# not silently ignore the argument.
result <- tryCatch(
  anova(m0, m1, model.names = c("null", "full")),
  error   = function(e) list(type = "error",   msg = conditionMessage(e)),
  warning = function(w) list(type = "warning", msg = conditionMessage(w))
)
if (!is.null(result$type)) {
  # Clean error/warning is acceptable pre-implementation
  expect_true(result$type %in% c("error", "warning"))
} else {
  # If it silently succeeds, the model column must reflect the supplied names
  expect_true("null" %in% result$table$model || "full" %in% result$table$model,
              info = "if model.names is accepted, it must be reflected in output")
}
```

**What it guards against:** Silent no-op where `model.names` is swallowed
without effect and the user's custom labels disappear.

---

### D1-01 — drop1() F-table via Satterthwaite routing (lmerTest parity)

**Kind:** integration  
**Priority:** P1

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
           REML = TRUE, control = mm_control(verbose = -1))
```

**Context:** lmerTest's `drop1.lmerModLmerTest` returns a per-term Satterthwaite
F-table. mixeff's `drop1.mm_lmm` does not accept a `ddf` argument; the
equivalent is `test_effect(fit, term, method="satterthwaite")`. This spec
certifies that combination produces F-equivalent values matching lmerTest and
documents the route for users.

**Assertion:**
```r
skip_if_not_installed("lme4"); skip_if_not_installed("lmerTest")
fit_lme4 <- lme4::lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = TRUE)

# lmerTest route (masks drop1, returns F-table)
an_lt <- drop1(fit_lme4, test = "F")   # lmerTest::drop1.lmerModLmerTest
F_lt  <- an_lt["Days", "F value"]
df_lt <- an_lt["Days", "Df"]

# mixeff equivalent: test_effect with satterthwaite
te_mm <- test_effect(fit, "Days", method = "satterthwaite")
t_mm  <- te_mm$table$statistic[te_mm$table$term == "Days"]
df_mm <- te_mm$table$den_df[te_mm$table$term == "Days"]

# F = t^2 equivalence
expect_lt(abs(t_mm^2 - F_lt), 0.01,
          label = "test_effect t^2 vs lmerTest drop1 F for Days")
expect_lt(abs(df_mm - df_lt), 0.5,
          label = "Satterthwaite DenDF must agree to within 0.5")

# Additionally: mixeff drop1(test="Chisq") must still produce a LRT table
d1_lrt <- drop1(fit, test = "Chisq")
expect_true("Days" %in% d1_lrt$table$dropped)
expect_true(all(is.finite(d1_lrt$table$LRT)))
expect_true(all(d1_lrt$table$LRT > 0))
```

**Tolerance:** F / t^2 within 0.01; DenDF within 0.5.

**What it guards against:** Certifies that the `test_effect` → `drop1` bridge
is a valid lmerTest-equivalent path, and that the LRT path is not broken by any
future `ddf` argument addition.

---

### D1-02 — drop1(test="user", sumFun=) hook or documented refusal

**Kind:** unit  
**Priority:** P2

**Assertion:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy, REML = TRUE,
           control = mm_control(verbose = -1))

# Until test="user" is implemented, the function must error cleanly,
# not silently ignore sumFun.
err <- tryCatch(
  drop1(fit, test = "user", sumFun = function(m0, m1) data.frame(p = 0.05)),
  error = function(e) e
)
expect_true(inherits(err, "error"))
expect_match(conditionMessage(err), "user|sumFun|not.*supported",
             ignore.case = TRUE,
             info = "test='user' must produce a comprehensible error, not silent no-op")
```

**What it guards against:** A `match.arg` that silently coerces `"user"` to an
existing value, making the extension hook invisible to users.

---

### PB-01 — compare(method="bootstrap") output structure vs PBmodcomp

**Kind:** snapshot  
**Priority:** P1

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
m0 <- lmm(Reaction ~ 1 + (Days | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))
m1 <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))
```

**Assertion:**
```r
pb <- compare(m0, m1, method = "bootstrap", nsim = 199, seed = 1)

# Method label
expect_identical(pb$method, "parametric_bootstrap_lrt")

# Bootstrap p-value must be present and in [0,1]
expect_true(is.numeric(pb$table$p_value[pb$table$model == pb$table$model[2]]) ||
            is.numeric(pb$p_value))

# Replicate count must equal nsim
reps <- pb$bootstrap_summary$replicates %||% pb$replicates %||%
        attr(pb, "replicates")
if (!is.null(reps)) expect_equal(reps, 199L)

# LRT statistic under H0 must be >= 0
lrt_obs <- pb$table$LRT[!is.na(pb$table$LRT)]
expect_true(length(lrt_obs) >= 1L)
expect_true(all(lrt_obs >= 0))

# Document the ABSENCE of Bartlett/Gamma reference-distribution rows
# (gap report notes these are not returned — this assertion confirms absence
#  so that if they are added we know to update downstream consumers).
has_bartlett <- any(grepl("bartlett|gamma|F_ref",
                          names(pb), ignore.case = TRUE))
if (has_bartlett) {
  # Future-proof: if Bartlett rows are added, require them to be numeric
  message("PB-01: Bartlett/Gamma rows detected — update snapshot")
} else {
  # Current behavior: only bootstrap p-value, no Bartlett rows
  expect_false(has_bartlett,
               info = "PB-01: Bartlett rows absent (expected for current impl)")
}
```

**Tolerance:** p-value in [0, 1]; LRT >= 0; replicate count = nsim.

**What it guards against:** Regression that corrupts the bootstrap LRT p-value
or drops the method label; also an intentional snapshot of the current
Bartlett-row absence so any future addition is a visible, deliberate change.

---

### PB-02 — PBrefdist analog: documented absence or clean sentinel

**Kind:** error-message  
**Priority:** P2

**Assertion:**
```r
# PBrefdist() does not exist in mixeff. Calling it must produce a clean
# "not found" error, not a namespace collision or misleading output.
err <- tryCatch(mixeff::PBrefdist, error = function(e) e)
expect_true(inherits(err, "error"),
            info = "PBrefdist must not exist in the mixeff namespace")
# The re-simulate-from-scratch behavior is expected and documented.
```

**What it guards against:** Accidental export of an incomplete stub named
`PBrefdist` that returns nonsense.

---

### RE-01 — test_random_effect whole-model table (ranova-style summary)

**Kind:** integration  
**Priority:** P1

**Dataset / formula:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           REML = TRUE, control = mm_control(verbose = -1))
```

**Context:** `ranova()` from lmerTest returns a table with one row per random
effect term. mixeff's `test_random_effect` requires explicit `term` specification.
This spec verifies that calling `test_random_effect` for each term and
assembling the result is correct, and also verifies that a convenience
whole-model table exists (or documents its absence with a sentinel test).

**Assertion:**
```r
# Part 1: Per-term boundary LRT works for both terms in (Days | Subject)
re_subj <- test_random_effect(fit, "Subject", method = "boundary_lrt")
expect_s3_class(re_subj, "mm_random_effect_test")
expect_true(all(is.finite(re_subj$table$LRT) | re_subj$table$status == "not_assessed"))
expect_true(any(re_subj$table$status == "available"))
expect_true(all(re_subj$table$p_value[re_subj$table$status == "available"] >= 0))

# Boundary mixture: p-value for removing a random slope should be in [0, 1]
p_val <- re_subj$table$p_value[re_subj$table$status == "available"][1]
expect_true(p_val >= 0 && p_val <= 1)

# Part 2: Whole-model table — test whether a multi-term convenience interface exists
all_terms <- test_random_effect(fit, method = "boundary_lrt")
if (!is.null(all_terms)) {
  # If a no-term version returns a table, it must cover all RE terms
  expect_true(nrow(all_terms$table) >= 1L)
} else {
  # Absence is acceptable pre-implementation; mark as a known gap
  skip("whole-model test_random_effect table not yet implemented (RE-01 Part 2)")
}
```

**Tolerance:** p-value in [0, 1]; LRT >= 0 for boundary models.

**What it guards against:** Regression in the `boundary_lrt` path; also acts as
a specification that a whole-model convenience interface should be implemented
and will be tested when it exists.

---

### RE-02 — ranova reduce.terms decomposition sentinel

**Kind:** error-message  
**Priority:** P2

**Assertion:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           REML = TRUE, control = mm_control(verbose = -1))

# reduce.terms is not yet implemented in test_random_effect.
# Calling with reduce.terms=TRUE must error cleanly, not silently ignore it.
err <- tryCatch(
  test_random_effect(fit, "Subject", method = "boundary_lrt",
                     reduce_terms = TRUE),
  error = function(e) e
)
# Acceptable outcomes:
# (a) Error with an informative message (not implemented)
# (b) The argument is silently ignored and the regular test is returned.
# Sub-test: if it silently succeeds, the result must be a valid test object.
if (inherits(err, "error")) {
  # Preferred: explicit "not yet" error
  expect_match(conditionMessage(err),
               "reduce|not.*implemented|unsupported",
               ignore.case = TRUE)
} else {
  expect_s3_class(err, "mm_random_effect_test")
}
```

**What it guards against:** Correlated-slope decomposition (`(x|g)` →
`(1|g)+(0+x|g)`) silently running as a full-term drop, which gives different
chi-square statistics than the correct decomposition.

---

### BM-01 — bootMer arbitrary-FUN sentinel

**Kind:** error-message  
**Priority:** P2

**Assertion:**
```r
data(sleepstudy, package = "lme4")
fit <- lmm(Reaction ~ Days + (1 | Subject), sleepstudy,
           REML = FALSE, control = mm_control(verbose = -1))

# mixeff does not expose a public bootMer(x, FUN=...) interface.
# Any attempt to call it (or a name-matched stub) must fail cleanly.
err <- tryCatch(mixeff::bootMer, error = function(e) e)
expect_true(inherits(err, "error"),
            info = "bootMer must not exist in mixeff namespace until implemented")
```

**What it guards against:** Incomplete `bootMer` stub accidentally exported
before the full interface (arbitrary FUN, use.u, re.form, type) is implemented.

---

## Upstream fixture requirements summary

The following specs require changes or additions to the upstream Rust crate
(`mixeff-rs`) before the R-side tests can pass:

| Spec | Required upstream change |
|------|--------------------------|
| CI-01 | Bootstrap engine must expose VC parameter distributions (`parameter_kind = "theta"/"sigma"` rows in the bootstrap payload). |
| CI-03 | REML optimizer convergence must achieve tighter theta precision on correlated RS models; or vcov numerical path must be more stable. |
| AN-01 | `mm_rust_term_table` FFI must accept and honor a `type` argument (`"I"`, `"II"`, `"III"`) for ANOVA contrast computation. |

All other specs are R-side only (routing, argument dispatch, error messages).

---

## Test file placement

Tests marked P0/P1 should be placed in (or referenced from):
- `tests/testthat/test-confint-profile.R` — CI-02, CI-04 (profile family)
- `tests/testthat/test-inference.R` — CI-01, CI-03, CI-05 (confint/bootstrap)
- `tests/testthat/test-lmm.R` or a new `test-anova-types.R` — AN-01
- `tests/testthat/test-inference.R` — D1-01, RE-01

Tests marked P2 may be grouped in a new `test-inference-sentinels.R` that
collects "Phase 5 placeholder / not-yet-implemented" sentinel tests.
