# Error-message quality: response variable is a character vector

**Scenario:** `character-response`
**Date probed:** 2026-05-31
**mixeff version:** installed from repo (main branch)
**lme4 version:** 2.0.1 / lmerTest 3.2.1 (R 4.5.2)

---

## Setup

```r
set.seed(42)
n <- 60
df <- data.frame(
  y_char = sample(c("low", "mid", "high"), n, replace = TRUE),  # 3-level character
  y_bin  = sample(c("yes", "no"), n, replace = TRUE),           # binary character
  x      = rnorm(n),
  subj   = factor(rep(1:10, each = 6))
)
```

---

## Verbatim messages

### lme4::lmer (y = character vector)

```
Class:  simpleError, error, condition
Message:
response must be numeric
```

### lme4::glmer (y = 3-level character vector, binomial family)

```
Class:  simpleError, error, condition
Message:
response must be numeric or factor
```

### lme4::glmer (y = binary character "yes"/"no", binomial family)

```
Class:  simpleError, error, condition
Message:
response must be numeric or factor
```

### mixeff::lmm (y = 3-level character vector)

```
Class:  mm_fit_error, mm_condition, rlang_error, error, condition
Message:
failed to construct LMM: Invalid argument: Response 'y_char' not found or not numeric
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct LMM: Invalid argument: Response 'y_char' not found or not numeric
Extra fields on condition: formula
```

### mixeff::glmm (y = 3-level character, binomial family)

```
Class:  mm_fit_error, mm_condition, rlang_error, error, condition
Message:
failed to construct GLMM: Invalid argument: Response 'y_char' not found or not numeric
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct GLMM: Invalid argument: Response 'y_char' not found or not numeric
Extra fields on condition: formula, metadata, spec
```

### mixeff::glmm (y = binary character "yes"/"no", binomial family)

```
Class:  mm_fit_error, mm_condition, rlang_error, error, condition
Message:
failed to construct GLMM: Invalid argument: Response 'y_bin' not found or not numeric
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct GLMM: Invalid argument: Response 'y_bin' not found or not numeric
Extra fields on condition: formula, metadata, spec
```

---

## Root-cause trace

`data-translate.R` (`mm_translate_data()`) routes character columns through the
categorical path — they are passed to Rust as `categorical_values` /
`categorical_levels`, not in `numeric_columns`. When `compile_model()` /
`lmm()` / `glmm()` then call the Rust FFI, the engine looks for the response
variable in the numeric column set, does not find it there, and raises:

> `Invalid argument: Response '<name>' not found or not numeric`

This message propagates back through `mm_abort_from_bridge()` as an
`mm_fit_error`. There is **no R-side pre-fit check** that detects a
character (or factor) response before the FFI boundary is crossed.

---

## Analysis

| Axis | lme4::lmer | lme4::glmer | mixeff::lmm / glmm |
|---|---|---|---|
| Error type | `simpleError` (untyped) | `simpleError` (untyped) | `mm_fit_error` (typed) |
| Message text | `response must be numeric` | `response must be numeric or factor` | `failed to construct LMM: Invalid argument: Response 'y_char' not found or not numeric` |
| Names the column? | No | No | Yes — `'y_char'` |
| Correct type constraint stated? | Yes — "numeric" | Yes — "numeric or factor" | Misleading — "not found or not numeric" implies column is absent |
| Actionable? | Partially — user knows to coerce | Partially — user knows to coerce | Partially — but confusing phrasing undermines it |
| Typed / machine-catchable? | No | No | Yes — `mm_fit_error` |
| Fires before FFI boundary? | Yes (R-side model.frame) | Yes (R-side model.frame) | No — reaches Rust before being caught |

**lme4** produces short, accurate messages: `"response must be numeric"` and
`"response must be numeric or factor"`. They do not name the offending column,
but the constraint is stated plainly. They fire at the R level (model-frame
construction) before any C call.

**mixeff** names the column (`'y_char'`), which is better than lme4 on that
dimension. However the message text has two problems:

1. **"not found or not numeric"** — the disjunctive phrasing implies the
   column might be missing from `data`. The user knows it is present; the
   "not found" branch refers to an internal Rust dataframe concept (the column
   was not placed into the numeric slab because it is categorical). This is an
   internal implementation detail leaking into the user-facing message.

2. **"failed to construct LMM: Invalid argument: ..."** — the outer wrapper
   phrase `failed to construct LMM` is generic. Together with the Rust
   sub-message it reads as a multi-layer stack fragment rather than a clean
   diagnosis.

The typed class (`mm_fit_error`) and the attached `formula` field are genuine
advantages over lme4, but the message text itself is comparable-but-improvable:
it names the column (better than lme4) yet confuses the user with "not found"
(worse than lme4's brevity).

A clean R-side pre-fit guard in `lmm()` / `glmm()` (or in `compile_model()`)
would catch this before any FFI crossing and could produce:

> `mm_data_error: Response variable 'y_char' must be numeric; got character.
>  Coerce with as.numeric() or recode to 0/1 before fitting.`

That message would be unambiguously clearer than lme4's.

---

## Verdict

**needs-work** — mixeff's message is typed and names the column, which are
genuine improvements over lme4. However the "not found or not numeric" phrasing
leaks an internal Rust dataframe concept ("not in the numeric slab") and reads as
though the column is missing, which is misleading. The error also fires at the
FFI boundary (as an `mm_fit_error`) rather than at the R level (where it should
be an `mm_data_error` fired before the Rust call). Both problems are fixable
with a single R-side guard.

**mixeff_quality: comparable** — the message names the column (better than lme4)
but the "not found or" phrase partially undermines clarity (worse than lme4's
concise constraint statement). Net: roughly comparable, with a clear path to
"clearer-than-lme4."

---

## Recommended fix

Add a response-type guard in `lmm()` (and the corresponding path in `glmm()`)
immediately after `compile_model()` returns, before the Rust fit call:

```r
resp_col <- spec$model_frame[[all.vars(formula[[2L]])]]
if (!is.numeric(resp_col)) {
  mm_abort(
    message = sprintf(
      "Response variable `%s` must be numeric; got %s.\n%s",
      deparse(formula[[2L]]),
      class(resp_col)[[1L]],
      if (is.character(resp_col))
        "Coerce with as.numeric() or recode to 0/1 before fitting."
      else if (is.factor(resp_col))
        "Coerce with as.numeric(as.character(.)) or use contrasts explicitly."
      else
        "Only numeric response variables are supported by lmm()."
    ),
    class = "mm_data_error",
    input = resp_col
  )
}
```

This fires as `mm_data_error` (R-level, before any FFI crossing), names the
column with correct context, states the constraint plainly, and offers a
concrete remedy — all without leaking Rust internals.
