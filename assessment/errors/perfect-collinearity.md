# Error-message quality probe: perfect collinearity

**Scenario:** Two perfectly collinear fixed-effect predictors (`x2 = 2 * x1`,
`cor(x1, x2) = 1.0`).  
**Date:** 2026-05-31  
**Probe script:** `assessment/errors/probe-perfect-collinearity.R`

---

## 1. Setup

```r
set.seed(42)
df <- data.frame(
  subject = factor(rep(seq_len(10), each = 4)),
  x1      = rnorm(40),
  y       = rnorm(40)
)
df$x2 <- 2 * df$x1   # perfect linear dependence
```

---

## 2. lme4 behaviour

### lme4::lmer (LMM)

lmer **silently fits** a rank-deficient model.  It drops `x2` and returns
estimates for `(Intercept)` and `x1` only, announcing this via two
**informational messages** (not errors, not typed conditions):

```
fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
boundary (singular) fit: see help('isSingular')
```

The returned `fixef()` contains only two terms; `x2` has been silently
aliased away.  `isSingular()` returns `TRUE`.

**Key lme4 characteristics:**
- Does not error; does not warn; signals via `message()`.
- Silently drops the aliased column from the returned coefficient vector.
- The user has no typed condition to `tryCatch`; the only hook is intercepting
  the messages.
- Does not name which pair of predictors is collinear, only that a column was
  dropped.
- No guidance on how to identify or fix the formula.

### lme4::glmer (Gaussian shortcut)

Identical messages to lmer, plus a deprecation notice.

### lme4::glmer (Binomial/logit)

Same messages; fits, drops `x2` silently.

---

## 3. mixeff behaviour

### mixeff::lmm — rank-deficient data (x2 = 2*x1)

**mixeff does not error.** It fits, returns a full `mm_lmm` object, and
communicates the problem through several layers:

#### 3a. Fit status on the returned object

```
fit_status: "converged_reduced_rank"
```

The enum value itself names the situation explicitly: the model converged on
a reduced-rank design.

#### 3b. Default print output

```
Linear mixed model fit by REML
Formula: y ~ x1 + x2 + (1 | subject)
Fit status: converged_reduced_rank
...
Fixed effects:
(Intercept)          x1          x2
  0.0793242   0.0000000  -0.0131722

Fitted covariance state:
The fitted covariance matrix is rank-deficient.
  r0: requested rank 1; fitted effective rank 0.
Use changes(fit) to see which dimension was unsupported.
Use random_options(spec, group = subject) to inspect lower-dimensional
  covariance choices.
Audit verbs: audit(), diagnostics(), inference_table(), model_report()
```

The print method surfaces the fit status and the covariance reduction. It
retains all three named coefficients in the output (zero-pinning `x1` and
keeping `x2` — which column is zeroed depends on the Rust solver's pivot
choice, not on which was listed first).

#### 3c. diagnostics()

```
Diagnostics:
                        code severity  stage         affected_terms
 fixed_effect_rank_deficient  warning  design_audit  x1
                  scope_note     info  design_audit  r0
          boundary_parameter     info  certification (1 | subject)
          covariance_reduced     info  certification (1 | subject)

Messages:
  fixed_effect_rank_deficient: fixed-effect formula is rank-deficient
    (rank 2 of 3); some requested coefficients are not separately estimable
    from the observed data
  scope_note: `x1` varies within `subject`, so a `subject`-level slope is
    structurally possible
  scope_note: `x2` varies within `subject`, so a `subject`-level slope is
    structurally possible
  boundary_parameter: standard deviation for intercept in (1 | subject)
    is on its lower bound
  covariance_reduced: fitted covariance for (1 | subject) has effective
    rank 0 of requested rank 1
```

`fixed_effect_rank_deficient` is a named diagnostic code with severity
`warning`, reporting the exact rank (2 of 3) and naming the affected term.

#### 3d. summary() — inference table

`x1` is marked `not_estimable` with reason
`"contrast touches aliased or non-finite coefficient directions"`.
`(Intercept)` and `x2` are marked `p_value_unavailable` / reliability `low`
with reason `"standard error is unavailable"`.
The method column says `not_computed` for `x1`.

#### 3e. inference_table()

Same information in tabular form; `status` column clearly distinguishes
`not_estimable` from `p_value_unavailable`.

### mixeff::glmm — rank-deficient data (binomial/logit, x2 = 2*x1)

Same pattern: fits without error, `(Intercept)` and `x1` or `x2` are
retained; the aliased direction is zeroed.  The print, diagnostics, and
inference table surfaces match lmm behaviour.

---

## 4. Side-by-side comparison

| Dimension | lme4::lmer | mixeff::lmm |
|---|---|---|
| Signals error? | No | No |
| Signals warning? | No (message only) | `diagnostics()` code severity=warning |
| Silent coefficient drop | Yes — x2 vanishes from fixef() | No — all three named; aliased one is 0 |
| Typed condition catchable | No | `fit_status = "converged_reduced_rank"` on object; `fixed_effect_rank_deficient` diagnostic code |
| Names the collinear term | No ("1 column" dropped) | Yes — `affected_terms: x1` (the zeroed coefficient) |
| Gives rank info | No | Yes — "rank 2 of 3" |
| Inference refusal | Silent (coefficient missing) | Explicit: `not_estimable`, `p_value_unavailable`, with reason strings |
| Actionable follow-up | None | `changes(fit)`, `random_options()`, `diagnostics()`, `inference_table()` |
| Audit trail | None | Named diagnostic codes, fit_status enum, covariance reduction report |

---

## 5. Issues found

### 5a. Silent wrong answer on coefficient attribution (needs-work)

mixeff returns all three coefficient names with a zero in the aliased slot
(which column is zeroed depends on the Rust solver's column-pivoting order,
not user intent). For `x2 = 2 * x1`:

```
fixef: (Intercept)=0.0793, x1=0.0000, x2=-0.0066
```

lme4 drops `x2` (the later column). mixeff zeros `x1` (the earlier one).
The zeroed coefficient is the wrong one to present to the user — `x2` was the
redundant predictor, not `x1`. This is a **silent wrong attribution** issue:
the user will read the output and believe `x1` contributes nothing, when in
fact it is `x2` that is linearly dependent. The `diagnostics()` `affected_terms`
field names `x1` as the affected term, but the human-visible message says
"some coefficients are not separately estimable" — it does not say "x2 is a
linear function of x1, so x2 is the aliased predictor." This is improvable
but not a crash.

### 5b. Covariance singularity conflated with fixed-effect collinearity (test-gap)

The printed `Fitted covariance state` block and the `covariance_reduced`
diagnostic both appear in the same output as the `fixed_effect_rank_deficient`
diagnostic. Both issues co-occur here because the collinear design also forces
a boundary random-effect solution. Users may conflate the two. No bug, but
the separation in `diagnostics()` (different `stage` values: `design_audit`
vs `certification`) is useful and worth documenting.

### 5c. No upfront refusal / pre-fit warning (needs-work)

mixeff's `compile_model()` / `explain_model()` phase has access to the design
matrix before handing off to Rust. A pre-fit rank check at that stage could
surface `fixed_effect_rank_deficient` as an informational message *before* the
Rust fit runs, matching mixeff's "audit-first" promise. Currently the
collinearity is only detectable post-fit via `diagnostics()` or `fit_status`.

---

## 6. Classification

| Issue | Classification |
|---|---|
| mixeff fits without aborting (same as lme4) | works — both fit rank-deficient |
| fit_status "converged_reduced_rank" on object | works |
| diagnostics() surfaces `fixed_effect_rank_deficient` with rank and term | works |
| summary/inference_table marks aliased term not_estimable | works |
| Wrong column zeroed (x1 instead of x2) without explaining why | needs-work |
| No pre-fit collinearity warning from compile_model | needs-work |
| No typed error raised on perfect collinearity | in-scope-missing / needs-work |

---

## 7. Verdict

**needs-work**

mixeff is materially **better** than lme4 on this scenario:

- It names the rank deficiency (`fixed_effect_rank_deficient`, rank 2 of 3).
- It preserves all coefficient names rather than silently dropping one.
- It provides typed, catchable diagnostics and a structured inference refusal
  (`not_estimable`) rather than a silent message.
- It surfaces actionable follow-up verbs.

However, two gaps keep it from "good":

1. The zeroed coefficient is the *first* collinear column, not the *redundant*
   one — the user would benefit from an explicit message such as "x2 is a
   linear function of x1; x2 has been constrained to zero." The current wording
   ("some requested coefficients are not separately estimable") doesn't identify
   the direction of the dependency.
2. The collinearity is detected only post-fit. An `explain_model()` / pre-fit
   rank check would match the audit-first design and catch the problem before
   the expensive Rust fit.

Neither gap is a crash, panic, or silent wrong answer at the fit level —
but the attribution ambiguity (which column is aliased and *why*) means a user
could be misled about their formula.
