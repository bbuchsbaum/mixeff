# Error-message quality probe: contrast-nonexistent-coef

**Scenario:** Calling `contrast()`, `test_effect()`, or `confint()` with a
coefficient or term name that does not exist in the fitted model.

**Probe script:** `assessment/errors/probe-contrast-nonexistent-coef.R`

---

## Setup

```r
set.seed(42)
dat <- data.frame(
  y       = rnorm(80, mean = 10, sd = 2),
  x       = rnorm(80),
  subject = factor(rep(1:20, each = 4))
)
dat$y_bin <- as.integer(dat$y > 10)
```

Fitted models:

- `lme4::lmer(y ~ x + (1 | subject), data = dat)` → fixed effects: `(Intercept)`, `x`
- `lme4::glmer(y_bin ~ x + (1 | subject), data = dat, family = binomial())`
- `mixeff::lmm(y ~ x + (1 | subject), data = dat)`
- `mixeff::glmm(y_bin ~ x + (1 | subject), data = dat, family = binomial())`

---

## lme4 behavior

### lme4::lmer / lme4::glmer — subscripting fixef() with a nonexistent name

lme4 has **no first-class `contrast()` verb with named-coefficient lookup**.
The canonical post-fit contrast path is `multcomp::glht()`.  The natural naive
attempt — subscripting `fixef()` — produces:

```
lme4::fixef(fit_lme4)["NONEXISTENT"]   # → NA (no error, no warning)
lme4::fixef(fit_glmer)["NONEXISTENT"]  # → NA (no error, no warning)
```

**Result: silent `NA`.** No error, no warning, no message. The user receives a
named `NA` and has no indication that the requested coefficient does not exist.
This is the classic "inscrutable lme4 failure mode" — the wrong answer arrives
quietly.

---

## mixeff behavior

### 1. `contrast(fit_mm, L)` — wrong column count (analog of wrong coef specification)

`contrast()` takes a numeric matrix `L`; the column dimension must equal the
number of fixed effects.  Passing a 3-column matrix for a 2-fixed-effect model:

```
ERROR class: mm_arg_error, mm_condition, rlang_error, error, condition
ERROR message: `L` must have 2 column(s), one for each fixed effect.
```

**Assessment:** Clear, typed, actionable. States exactly what is wrong (column
count mismatch) and what is expected (2 columns, one per fixed effect).

### 2. `test_effect(fit_mm, "NONEXISTENT_TERM")` — nonexistent term name

```
ERROR class: mm_arg_error, mm_condition, rlang_error, error, condition
ERROR message: Unknown fixed-effect term(s): NONEXISTENT_TERM.
```

**Assessment:** Clear and direct. Names the unknown term(s) explicitly.
Does not list the available terms in the error message itself — a minor gap
(user must call `mm_fixed_effect_terms(fit)` or read the model summary to
discover what names are valid). However the class `mm_arg_error` signals the
problem type unambiguously, and the available terms are queryable.

### 3. `confint(fit_mm, parm = "NONEXISTENT_COEF")` — nonexistent parameter name

```
ERROR class: mm_arg_error, mm_condition, rlang_error, error, condition
ERROR message: Unknown fixed-effect parameter(s): NONEXISTENT_COEF.
```

**Assessment:** Clear, typed, actionable. Identical pattern to `test_effect`.
Same minor gap: available parameter names not listed inline.

### 4. `contrast(fit_glmm, L)` — no `contrast.mm_glmm` method registered

```
ERROR class: simpleError, error, condition
ERROR message: no applicable method for 'contrast' applied to an object of
class "c('mm_glmm', 'mm_fit', 'mm_compiled')"
```

**Assessment:** This is an **untyped R dispatch error** — a `simpleError`, not
an `mm_arg_error` or `mm_condition`. The message is R's generic S3 dispatch
failure. It tells the user there is no `contrast()` method for `mm_glmm`, but
does not explain whether GLMM contrasts are out-of-scope by design, planned for
a future version, or simply missing. Per the PRD `§3` non-goals, GLMM
post-fit inference via `contrast()` is not listed as certified in v1, but it is
also not explicitly excluded. This gap is a **needs-work** issue: the message
should be either a typed `mm_not_implemented` condition with an explanation, or
a documented statement that post-fit GLMM contrasts are deferred.

### 5. Named `L` with a wrong rowname (semantic label only)

```r
L <- matrix(c(0, 1), nrow = 1)
rownames(L) <- "NONEXISTENT_COEF_LABEL"
colnames(L) <- c("(Intercept)", "x")
contrast(fit_mm, L)
```

Result: **no error**. The contrast is computed correctly (the label is just a
display name for the row). The returned `contrast` column in the table uses
the user-supplied rowname as-is. This is correct behavior: `L` is a numeric
hypothesis matrix, not a lookup by name.

---

## Comparison summary

| Scenario | lme4 | mixeff | Winner |
|---|---|---|---|
| Subscript fixef() with nonexistent name | Silent `NA` — no error | N/A (different API) | mixeff (explicit typed errors) |
| `test_effect("NONEXISTENT")` | N/A (no `test_effect`) | Typed `mm_arg_error`: "Unknown fixed-effect term(s): NONEXISTENT_TERM." | mixeff clearly better |
| `confint(parm="NONEXISTENT")` | Silently drops or subscript-NA | Typed `mm_arg_error`: "Unknown fixed-effect parameter(s): NONEXISTENT_COEF." | mixeff clearly better |
| Wrong-dimension `L` to `contrast()` | N/A | Typed `mm_arg_error`: "`L` must have 2 column(s), one for each fixed effect." | mixeff clearly better |
| `contrast()` on `mm_glmm` | N/A | Untyped `simpleError`: S3 dispatch failure | **needs-work** |

---

## Verdict

**`lmm()` path: GOOD.** All three named-coef error scenarios (`test_effect`,
`confint`, wrong-dim `L`) produce typed `mm_arg_error` conditions with clear,
actionable messages that are substantially better than lme4's silent-`NA`
behavior.

**`glmm()` path: NEEDS-WORK.** The absence of `contrast.mm_glmm` produces a
raw R S3 dispatch error (`simpleError`) with no guidance. It should either:
(a) register a stub method that raises a typed `mm_not_implemented` or
`mm_inference_unavailable` condition explaining the scope boundary, or
(b) be documented as an explicit non-goal with a pointer to the GLMM
asymptotic inference path (`summary(fit, tests="coefficients")`).

**Minor improvement opportunity (both paths):** Error messages for unknown
term/parameter names could list the available names inline (e.g.
"Unknown fixed-effect term(s): NONEXISTENT_TERM. Available terms: 1, x."),
saving the user a follow-up query.

---

## Raw output

```
=== lme4::lmer — accessing nonexistent fixed-effect term ===
Fixed-effect names: (Intercept), x

--- lme4 fixef()['NONEXISTENT'] ---
Result (no error): NA

=== lme4::glmer — accessing nonexistent fixed-effect term ===
Fixed-effect names: (Intercept), x

--- lme4 glmer fixef()['NONEXISTENT'] ---
Result (no error): NA

=== mixeff::lmm — contrast() with incorrect column count ===
ERROR class: mm_arg_error, mm_condition, rlang_error, error, condition
ERROR message: `L` must have 2 column(s), one for each fixed effect.

=== mixeff::lmm — test_effect() with nonexistent term ===
Available fixed-effect terms: 1, x

--- mixeff test_effect(fit, 'NONEXISTENT_TERM') ---
ERROR class: mm_arg_error, mm_condition, rlang_error, error, condition
ERROR message: Unknown fixed-effect term(s): NONEXISTENT_TERM.

=== mixeff::lmm — confint() with nonexistent parm name ===
ERROR class: mm_arg_error, mm_condition, rlang_error, error, condition
ERROR message: Unknown fixed-effect parameter(s): NONEXISTENT_COEF.

=== mixeff::glmm — contrast() dispatched to glmm (no method registered) ===
ERROR class: simpleError, error, condition
ERROR message: no applicable method for 'contrast' applied to an object of
class "c('mm_glmm', 'mm_fit', 'mm_compiled')"

=== mixeff::lmm — contrast() with named L, wrong rowname (semantic test) ===
Result (no error): returned, contrast is valid numerically
Contrast rowname used: NONEXISTENT_COEF_LABEL
  contrast: NONEXISTENT_COEF_LABEL  estimate: 0.262  p_value: 0.322  status: available
```
