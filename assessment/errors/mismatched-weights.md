# Error-message quality: mismatched-weights

**Scenario:** `weights` vector has the wrong length (too short, too long, or
zero-length) when calling `lmm()` or `glmm()`.

**Probe script:** `assessment/errors/probe-mismatched-weights.R`

---

## Verbatim messages

### lme4 (all three bad-length cases — lmer and glmer give the same message)

```
variable lengths differ (found for '(weights)')
```

Condition class: base `simpleError` (no typed class hierarchy).

---

### mixeff::lmm — too-short weights (length 20, need 40)

```
`weights` must be a finite positive numeric vector with one value per row in `data`.
```

Condition class: `mm_data_error, mm_condition, rlang_error, error, condition`

---

### mixeff::lmm — too-long weights (length 80, need 40)

```
`weights` must be a finite positive numeric vector with one value per row in `data`.
```

Same class hierarchy as above.

---

### mixeff::lmm — zero-length weights

```
`weights` must be a finite positive numeric vector with one value per row in `data`.
```

Same class hierarchy as above.

---

### mixeff::glmm — any non-NULL weights (reserved parameter)

```
`random`, `weights`, `subset`, custom `na.action`, and `contrasts` are reserved for the fitted GLMM bridge.
```

Condition class: `mm_fit_error, mm_condition, rlang_error, error, condition`

---

## Comparison and assessment

| Dimension | lme4 | mixeff::lmm | mixeff::glmm |
|---|---|---|---|
| Message text | `variable lengths differ (found for '(weights)')` | `` `weights` must be a finite positive numeric vector with one value per row in `data`. `` | Reserved-param error listing all reserved args |
| Names the bad parameter | Yes, indirectly via `'(weights)'` | Yes, explicitly via `` `weights` `` | Yes, explicitly |
| States what is required | No | Yes — "finite positive numeric vector with one value per row in `data`" | N/A (feature not yet wired) |
| Typed condition class | No — base `simpleError` | `mm_data_error` | `mm_fit_error` |
| Actionable? | Barely — user must infer that weights length != nrow(data) | Clear and direct | Partial — message explains reserved status but not when weights will be supported |
| Silent wrong answer? | No | No | No |
| Crash / panic / stack trace? | No | No | No |

### Summary

**mixeff::lmm is clearly better than lme4** for this scenario. lme4 emits a
terse, model-frame-level message (`variable lengths differ`) that does not name
the constraint or explain what is acceptable. mixeff::lmm states exactly what
is required: a finite positive numeric vector with one entry per row. It also
carries a typed `mm_data_error` class that is catchable programmatically, which
lme4's untyped base error is not.

One minor gap: the message does not echo the *observed* length vs the
*expected* length (e.g. "got 20, expected 40"). Adding that detail would make
it marginally more actionable but the current message is already unambiguously
clearer than lme4.

**mixeff::glmm** raises `mm_fit_error` with a clear "reserved" message for any
non-NULL weights. This is honest and non-silent but is categorically different
from the lmm path — weights are not validated for length because the argument
is not yet wired through. The message does not tell the user *when* weights
will be supported or *what* they should do in the meantime (use `lmm()` if
applicable). This is a minor UX gap but not a bug given the stated PRD
out-of-scope status.

### Verdict

- `lmm()`: **good** — clearer than lme4, typed, actionable.
- `glmm()`: **needs-work** — honest reserved-param error, but the message
  conflates several reserved args without guidance on alternatives; does not
  validate length at all (not a length-check path). Not a bug (weights are
  out-of-scope for GLMM in this build), but the message could be improved to
  say "weights are not yet supported for glmm(); use lmm() for weighted
  LMM fits."
