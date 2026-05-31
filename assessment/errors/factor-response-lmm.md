# Error-message quality: factor response passed to lmm / glmm

**Scenario:** `factor-response-lmm`
**Date probed:** 2026-05-31
**mixeff version:** installed from repo (main branch)
**lme4 version:** 2.0.1 / lmerTest 3.2.1 (R 4.5.2)

---

## Setup

```r
set.seed(42)
n <- 40
# 3-level factor response
df <- data.frame(
  y       = factor(sample(c("low", "mid", "high"), n, replace = TRUE)),
  x       = rnorm(n),
  subject = factor(rep(seq_len(10), each = 4))
)
# binary factor response (for glmer/glmm binomial)
df2 <- df
df2$y <- factor(sample(c("no", "yes"), n, replace = TRUE))
```

---

## Verbatim messages

### lme4::lmer (3-level factor response)

```
Class:  simpleError, error, condition
Message:
response must be numeric
```

### lme4::glmer (gaussian, factor response)

First emits a deprecation warning (`calling glmer() with family=gaussian ... as a shortcut to lmer() is deprecated`), then delegates to lmer and produces:

```
Class:  simpleError, error, condition
Message:
response must be numeric
```

### lme4::glmer (binomial, binary factor response)

lme4 **succeeds** — it converts the binary factor to 0/1 silently. Returns a `glmerMod` object (with a singular-fit warning). No error is raised.

### mixeff::lmm (3-level factor response)

```
Class:  mm_fit_error, mm_condition, rlang_error, error, condition
Message:
failed to construct LMM: Invalid argument: Response 'y' not found or not numeric
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct LMM: Invalid argument: Response 'y' not found or not numeric
```

### mixeff::lmm (binary factor response)

Identical to the 3-level case:

```
Class:  mm_fit_error, mm_condition, rlang_error, error, condition
Message:
failed to construct LMM: Invalid argument: Response 'y' not found or not numeric
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct LMM: Invalid argument: Response 'y' not found or not numeric
```

### mixeff::glmm (binomial, binary factor response)

Prints the full design-audit block first (random-effects explanation, design notes), then errors:

```
[... full design audit output ...]

Class:  mm_fit_error, mm_condition, rlang_error, error, condition
Message:
failed to construct GLMM: Invalid argument: Response 'y' not found or not numeric
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct GLMM: Invalid argument: Response 'y' not found or not numeric
```

### mixeff::glmm (poisson, 3-level factor response)

Same structure: design audit prints, then identical bridge error.

---

## Analysis

| Axis | lme4::lmer | mixeff::lmm | mixeff::glmm |
|---|---|---|---|
| Error type | `simpleError` (untyped) | `mm_fit_error` (typed, but wrong class) | `mm_fit_error` (typed, but wrong class) |
| Message text | `response must be numeric` | `failed to construct LMM: Invalid argument: Response 'y' not found or not numeric` | `failed to construct GLMM: Invalid argument: Response 'y' not found or not numeric` |
| Names the variable? | No | Yes — `'y'` | Yes — `'y'` |
| States the type constraint? | Yes — "must be numeric" | Ambiguously — "not found or not numeric" conflates two cases | Same |
| Actionable remedy? | No | No | No |
| Pre-fit check (R side)? | Yes | No — escapes to Rust FFI | No — escapes to Rust FFI |
| Design audit printed before error? | N/A | No | Yes — confusing partial-success impression |
| Condition class correct? | N/A | No — data-type error misclassed as `mm_fit_error` | No — data-type error misclassed as `mm_fit_error` |

### Specific defects

**Defect 1 — Ambiguous "not found or not numeric" message (bridge string leak)**

The Rust error string `Response 'y' not found or not numeric` leaks an internal
disjunction that is confusing in this context. The variable `y` IS found in the
data; what failed is that its type is `factor` / categorical. A user reading
"not found or not numeric" will reasonably wonder whether they mistyped the
variable name. The check that fires here is a type check, not a name lookup, so
the message should say only "not numeric" (or better: "must be numeric; got
factor").

**Defect 2 — Wrong condition class (`mm_fit_error` instead of `mm_data_error`)**

A factor response is a data preparation / type-mismatch problem, not a
model-construction or optimization failure. The correct class is `mm_data_error`,
consistent with how `mm_translate_data` handles other unsupported types (it raises
`mm_data_error` for `Date`, `POSIXct`, etc.). Misfiling this as `mm_fit_error`
means code that catches `mm_data_error` to detect input problems will silently
miss this case.

**Defect 3 — Check fires at Rust FFI rather than in R pre-flight**

`mm_translate_data` already iterates all columns and classifies them. A factor
column is silently promoted to categorical (correct for a grouping variable), but
there is no check that the *response* column (the LHS of the formula) ends up in
`numeric_columns`. The R side knows the response name from the formula; adding a
post-translate assertion before the `.Call` would catch this with a clean
`mm_data_error` before any Rust boundary is crossed. Currently the check leaks
through to the Rust layer, producing the ambiguous bridge message above.

**Defect 4 — Design audit prints before the error in glmm paths**

For `glmm()`, the full `explain_model()` output is printed before the type error
fires. This creates a misleading partial-success impression: the user sees a
detailed analysis of their random-effects structure and then an abrupt error
message. Because the error is about the response type — a fundamental input
problem — the design audit should never have printed. The type check must precede
`compile_model()`.

### Comparison to lme4

lme4's `response must be numeric` is terse and untyped, but it states the
violated constraint clearly and unambiguously in four words. A user immediately
understands the problem.

mixeff's current message has two compounding problems: the "not found **or** not
numeric" phrasing introduces false ambiguity (was the variable missing?), and
the `mm_fit_error` class misleads programmatic callers. The Rust bridge error
string is worse than lme4's in the clarity dimension, despite mixeff having more
context available (it could name the type: `factor`).

Note: lme4::glmer with binomial family and a binary factor **succeeds** by
silently coercing the factor to 0/1. This is lme4's "convenient but
non-transparent" behavior. mixeff's refusal here is consistent with the
audit-first no-silent-surgery design — but the *message* must clearly explain
what the user should do (e.g., `as.integer(y == "yes")`).

---

## Required fixes

1. **Add a response-type check in R before the FFI call** (in `lmm()` and
   `glmm()`): after `compile_model()` / `mm_translate_data()`, verify that the
   LHS variable is in `spec_data$numeric_columns`. If it is in
   `categorical_values` instead, raise:

   ```
   mm_data_error: Response variable `y` is a factor (categorical), but lmm()
   requires a numeric response. Convert to numeric before fitting, e.g.:
   data$y <- as.numeric(data$y)
   ```

   This fires before any Rust call, is typed correctly as `mm_data_error`, and
   provides a concrete remedy.

2. **Fix the condition class** — change from `mm_fit_error` to `mm_data_error`
   wherever a response-type mismatch is reported.

3. **For glmm()** — move the response-type check before `compile_model()` so the
   design audit is never printed when the response is fundamentally wrong.

---

## Verdict

**bug** — The message leaks an ambiguous Rust internal string ("not found or
not numeric") that is less clear than lme4's terse `response must be numeric`.
The condition class is wrong (`mm_fit_error` instead of `mm_data_error`). For
glmm paths, the design audit prints before the error, creating a misleading
partial-success impression. The fix is a pre-FFI R-side type check on the
response variable.

**mixeff quality vs lme4:** `worse-than-lme4` on message clarity (ambiguous
disjunction vs. direct constraint statement), though mixeff correctly names the
variable. The class mismatch and premature audit output are additional regressions
beyond the message text itself.
