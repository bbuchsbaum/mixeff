# Error-message quality: response variable not in data

**Scenario:** `response-missing`
**Date probed:** 2026-05-31
**mixeff version:** installed from repo (main branch)
**lme4 version:** 2.0.1 (R 4.5.2)

---

## Setup

```r
df <- data.frame(
  x       = rnorm(40),
  subject = factor(rep(seq_len(10), each = 4))
)
# 'y' is named in the formula but deliberately absent from df
```

---

## Verbatim messages

### lme4::lmer

```
Class:  simpleError, error, condition
Message:
object 'y' not found
```

### lme4::glmer

```
Class:  simpleError, error, condition
Message:
object 'y' not found
```

### mixeff::lmm

```
Class:  mm_data_error, mm_condition, rlang_error, error, condition
Message:
Variable(s) named in `formula` not found in `data`: `y`.
```

### mixeff::glmm

```
Class:  mm_data_error, mm_condition, rlang_error, error, condition
Message:
Variable(s) named in `formula` not found in `data`: `y`.
```

---

## Analysis

| Axis | lme4 | mixeff |
|---|---|---|
| Error type | `simpleError` (untyped) | `mm_data_error` (typed, machine-catchable) |
| Message text | `object 'y' not found` — inherited from R's `model.frame()` / `eval` machinery; no mention of `formula`, `data`, or what was expected | `Variable(s) named in \`formula\` not found in \`data\`: \`y\`.` — names the relationship, names the variable |
| Actionable? | No — user must infer that `y` is the formula response and that it should be a column | Yes — clearly says the variable is named in the formula but absent from data |
| Typed / catchable? | No | Yes — `mm_data_error` allows `tryCatch(..., mm_data_error = ...)` |
| Names both artefacts (formula + data)? | No | Yes |

**lme4's message** (`object 'y' not found`) leaks R's internal `eval` context. A user who does not know R internals will not understand that "object 'y'" means "column y is not in the data frame you supplied." The message does not distinguish between a missing data column, a missing R variable in the global env, and a typo in the formula.

**mixeff's message** is specific and audit-first: it tells the user exactly which name is absent, where it was expected (the formula), and where it was looked for (the data). The typed condition class `mm_data_error` allows programmatic handling. Both `lmm()` and `glmm()` produce the same quality of message because both delegate formula resolution to `compile_model()`, which validates variable presence before any Rust call.

---

## Verdict

**clearer-than-lme4** — mixeff's `mm_data_error` message is materially clearer: it names the missing variable, explains why it is missing (named in formula but absent from data), and is typed. lme4's `object 'y' not found` is an opaque leaked eval error.

**Overall verdict: good** — no defect; this scenario is handled well and is a concrete illustration of the "less inscrutable errors" promise.

---

## Possible improvement (not a bug)

The message could additionally mention *which side of the formula* the variable appears on, e.g. "response variable `y` named in formula LHS not found in `data`." This would make it marginally more helpful for the common user mistake of naming a response that doesn't exist, vs. a predictor typo. This is a polish item, not a defect.
