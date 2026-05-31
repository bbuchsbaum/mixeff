# Error-message quality assessment: `misspelled-var`

**Scenario:** formula references a column not present in `data` — either a
misspelled fixed-effect predictor (`xx` instead of `x`) or a misspelled
grouping variable (`gg` instead of `g`).

**Date probed:** 2026-05-31  
**mixeff version:** installed at /Users/bbuchsbaum/code/mixeff  
**lme4 version:** 2.0.1  

---

## Verbatim messages

### lme4::lmer — misspelled fixed predictor (`xx`)

```
class: simpleError, error, condition
message: object 'xx' not found
```

### lme4::glmer (binomial) — misspelled fixed predictor (`xx`)

```
class: simpleError, error, condition
message: object 'xx' not found
```

### lme4::lmer — misspelled grouping variable (`gg`)

```
class: simpleError, error, condition
message: object 'gg' not found
```

---

### mixeff::lmm — misspelled fixed predictor (`xx`)

```
class: mm_data_error, mm_condition, rlang_error, error, condition
message: Variable(s) named in `formula` not found in `data`: `xx`.
```

### mixeff::lmm — misspelled grouping variable (`gg`)

```
class: mm_data_error, mm_condition, rlang_error, error, condition
message: Variable(s) named in `formula` not found in `data`: `gg`.
```

### mixeff::glmm (binomial) — misspelled fixed predictor (`xx`)

```
class: mm_data_error, mm_condition, rlang_error, error, condition
message: Variable(s) named in `formula` not found in `data`: `xx`.
```

### mixeff::glmm (binomial) — misspelled grouping variable (`gg`)

```
class: mm_data_error, mm_condition, rlang_error, error, condition
message: Variable(s) named in `formula` not found in `data`: `gg`.
```

---

## Side note: glmm with unsupported family

When `glmm()` is called with `family = gaussian()` (not in the certified
contract), the family-validation guard fires *before* the variable check,
producing an unrelated `mm_inference_unavailable` error. This is correct
priority ordering (validate arguments before touching data), but it means the
misspelled-variable check is masked in that case. Not a bug — just a note that
lme4's `glmer(..., family=gaussian())` issues a deprecation *warning* instead
of an error, so it happens to reach model-frame evaluation and then surfaces
the `object 'xx' not found` error. This divergence is consistent with
mixeff's stricter upfront contract enforcement.

---

## Assessment

| Dimension | lme4 | mixeff |
|---|---|---|
| Error class | `simpleError` (untyped) | `mm_data_error` + `mm_condition` (typed) |
| Message | `object 'xx' not found` | `Variable(s) named in \`formula\` not found in \`data\`: \`xx\`.` |
| Names the missing variable | Yes (bare name) | Yes (backtick-quoted) |
| Names the *source* of the problem | No — says "object" not "column in data" | Yes — explicitly says "named in \`formula\`" and "not found in \`data\`" |
| Handles grouping var same as fixed | Yes | Yes |
| Machine-catchable by class | No (untyped `simpleError`) | Yes (`mm_data_error`) |
| Audit-first: no fabrication | n/a | Correct — fails before fitting |

### Summary

mixeff's message is **clearer than lme4's** on this scenario. lme4 emits a
bare R evaluation error (`object 'xx' not found`) that does not tell the user
the variable was expected from the data frame or that the formula is the
source. mixeff fires early in `compile_model()` (lines 92–103 of
`R/compile.R`), explicitly names the formula as the source and the data frame
as the lookup target, quotes the offending variable name(s), handles multiple
missing vars in a single message, and emits a typed `mm_data_error` condition
that callers can catch programmatically. No crash, no panic, no silent wrong
result.

**Verdict: good.** Meets the "clearer than lme4" promise for this class of
input error.

### Source location

`R/compile.R` lines 92–103 — `all.vars(formula)` diff against `names(data)`.
