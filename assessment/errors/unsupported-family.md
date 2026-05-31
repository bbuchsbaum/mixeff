# Error-message quality probe: `unsupported-family`

**Date**: 2026-05-31  
**mixeff version**: installed (main branch)  
**lme4 version**: 2.0.1  
**Probe script**: `/tmp/probe-unsupported-family.R`

---

## Summary verdict

**good** — mixeff's unsupported-family error path is clearer than lme4's in every
tested scenario. The messages are typed, named, and include a machine-readable
supported-family table as a condition field. No panics, no stack-trace-only
errors, no silent wrong answers.

---

## Scenarios tested

Five scenarios were exercised for both `lme4::glmer` and `mixeff::glmm`.

---

### Scenario A — fully unsupported family (`inverse.gaussian`)

**lme4::glmer**
```
class: simpleError, error, condition
message: positive values only are allowed for the 'inverse.gaussian' family
```
lme4 rejects at the data-validation level (negative/zero response values), not
at the family-dispatch level. The message is a runtime data check, not a family
support signal; it would not fire with all-positive data even though
`inverse.gaussian` is still unsupported by mixeff.

**mixeff::glmm**
```
class: mm_inference_unavailable, mm_condition, rlang_error, error, condition
message: GLMM family/link `inverse.gaussian/1/mu^2` is outside the certified upstream contract.
reason_code: unsupported_glmm_family_link
```
The condition fires before any data is touched, names the family and link
explicitly, carries `reason_code = "unsupported_glmm_family_link"`, and attaches
a `$supported` data frame (see below). Strictly clearer and earlier than lme4.

---

### Scenario B — supported family, unsupported link (`binomial/cauchit`)

**lme4::glmer**
```
class: glmerMod
(no error — fit succeeded)
```
lme4 accepts `binomial(link="cauchit")` and returns a fit object. mixeff refuses
this because `cauchit` is outside the certified Rust engine contract.

**mixeff::glmm**
```
class: mm_inference_unavailable, mm_condition, rlang_error, error, condition
message: GLMM family/link `binomial/cauchit` is outside the certified upstream contract.
reason_code: unsupported_glmm_family_link
```
This is a deliberate scope-narrowing (PRD §3 non-goals / certified surface). The
error is clear and typed. The divergence from lme4 is by design, not a defect.

---

### Scenario C — non-family object passed (character string `"banana"`)

**lme4::glmer**
```
class: simpleError, error, condition
message: object 'banana' of mode 'function' was not found
```
lme4 tries to evaluate the string as a function name and fails with a bare R
lookup error that does not mention `family`.

**mixeff::glmm**
```
class: mm_arg_error, mm_condition, rlang_error, error, condition
message: `family` must be an R family object or family constructor.
```
Clearly names the argument, states the requirement. Better than lme4.

---

### Scenario D — missing `family` argument entirely

**lme4::glmer**
```
class: simpleWarning, warning, condition
message: calling glmer() with family=gaussian (identity link) as a shortcut to lmer()
         is deprecated; please call lmer() directly
```
lme4 silently falls back to `gaussian` (deprecated path), issuing only a warning.
The call proceeds; no hard stop.

**mixeff::glmm**
```
class: mm_arg_error, mm_condition, rlang_error, error, condition
message: `family` is required for `glmm()`.
```
Hard error, names the missing argument. Per audit-first design this is the
correct behaviour.

---

### Scenario E — Poisson with unsupported link (`poisson/identity`)

**lme4::glmer**
```
class: std::runtime_error, C++Error, error, condition
message: PIRLS loop resulted in NaN value
```
lme4 accepts the link, attempts numerical optimisation, and fails with an opaque
C++ runtime error that exposes the internal algorithm name (PIRLS). The root
cause (non-canonical link causing numerical instability) is invisible.

**mixeff::glmm**
```
class: mm_inference_unavailable, mm_condition, rlang_error, error, condition
message: GLMM family/link `poisson/identity` is outside the certified upstream contract.
reason_code: unsupported_glmm_family_link
```
Refuses before any computation, names the cause. Clearly better than lme4's
cryptic NaN message.

---

## Attached condition fields (Scenario B, representative)

The `mm_inference_unavailable` condition carries:

```
$reason_code: "unsupported_glmm_family_link"
$family:      "binomial"
$link:        "cauchit"
$supported:
    family    link
1 binomial   logit
2 binomial  probit
3 binomial cloglog
4  poisson     log
5  poisson    sqrt
6    Gamma     log
```

A caller can catch `mm_inference_unavailable`, read `$reason_code` to confirm
the cause, and display `$supported` without re-parsing the message string.

---

## Classification

| Scenario | Gap classification |
|---|---|
| A — inverse.gaussian | `works` (clearer-than-lme4; lme4 msg is a data check, not family-support) |
| B — binomial/cauchit | `out-of-scope-by-design` (certified surface; lme4 silently fits) |
| C — non-family object | `works` (clearer than lme4's lookup error) |
| D — missing family | `works` (hard error vs lme4 deprecated silent fallback) |
| E — poisson/identity | `works` (clearer than lme4's C++ PIRLS NaN) |

No panics, no stack-trace-only errors, no silent wrong results.

---

## Potential improvement (needs-work note, not a bug)

The message text `"is outside the certified upstream contract"` is accurate but
slightly abstract. Adding one sentence — e.g. "The certified families are:
binomial (logit, probit, cloglog), poisson (log, sqrt), Gamma (log). Use
`mm_glmm_supported_family_link_table()` for a programmatic list." — directly in
the message string (not just as a condition field) would make the error
self-contained without requiring the caller to inspect `$supported`. This is a
polish item, not a defect.
