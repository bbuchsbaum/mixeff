# Error-message quality probe: no-random-effect

**Scenario:** formula has no `(.|.)` random-effect term (e.g. `y ~ x`).
**Question:** does mixeff route clearly, or crash/confuse?

---

## Verbatim messages

### lme4::lmer — `y ~ x`

```
Error class : simpleError, error, condition
Message     : No random effects terms specified in formula
```

### lme4::glmer — `y ~ x`

```
Error class : simpleError, error, condition
Message     : No random effects terms specified in formula
```

lme4 raises a plain `simpleError` with a terse, undifferentiated message. It
names the symptom but gives no hint about what to do next (e.g. "use lm()" or
"add a (1|group) term"). The condition class carries no structural information.

---

### mixeff::lmm — `y ~ x`

```
Error class : mm_fit_error, mm_condition, rlang_error, error, condition
Message     : failed to construct LMM: No random effects in formula: this is not a mixed model
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct LMM: No random effects in formula: this is not a mixed model
```

### mixeff::glmm — `y ~ x`, `family = binomial()`

```
Error class : mm_fit_error, mm_condition, rlang_error, error, condition
Message     : failed to construct GLMM: No random effects in formula: this is not a mixed model
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct GLMM: No random effects in formula: this is not a mixed model
```

### mixeff::glmm — `y ~ x`, `family = gaussian()` (unsupported family, masks the random-effect issue)

```
Error class : mm_inference_unavailable, mm_condition, rlang_error, error, condition
Message     : GLMM family/link `gaussian/identity` is outside the certified upstream contract.
```

The gaussian/identity case fires the family-validation guard in `glmm()` before
formula compilation is reached, so the missing-random-effect problem is never
surfaced. This is expected ordering (bad family is caught first) but worth
noting: a user writing `glmm(y ~ x, family = gaussian())` gets only a family
error, not a family + random-effect error.

---

## Analysis

### What mixeff does well

1. **Typed condition class.** `mm_fit_error` is a structured subclass of
   `mm_condition`. Downstream code can catch it specifically or generically,
   which lme4's `simpleError` does not permit.

2. **Semantically correct.** "No random effects in formula: this is not a mixed
   model" correctly names both the structural deficit and its implication.

3. **No panic / segfault.** The Rust engine returns a tagged error string
   (`mm_fit_error: …`); the R bridge routes it to the right typed condition via
   `mm_split_tagged_error` / `mm_abort_from_bridge`. No stack-trace-only crash.

4. **Audit-first contract respected.** The error is raised, not a fabricated
   fit. No silent wrong answer is returned.

### What needs work

1. **Noisy "Caused by" chain.** The rlang parent chain repeats the same message
   verbatim with a `doTryCatch()` frame attached:

   ```
   failed to construct LMM: No random effects in formula: this is not a mixed model
   Caused by error in `doTryCatch()`:
   ! mm_fit_error: failed to construct LMM: …
   ```

   The duplication adds noise without information. The `doTryCatch()` frame is
   an internal implementation detail that should not be visible to users. This
   is a `needs-work` cosmetic defect: the parent condition should be suppressed
   or the bridge should not re-wrap the already-tagged string as a parent.

2. **No actionable next step.** Neither lme4 nor mixeff tells the user "add a
   `(1 | group)` term or use `lm()` instead." mixeff's message is no worse than
   lme4 on this axis, but it could be better given the project's stated goal of
   fewer inscrutable errors. A `cnd_footer` or a `hint` field suggesting `lm()`
   would be a meaningful improvement.

3. **"failed to construct LMM:" prefix is redundant.** The error is already
   classed as `mm_fit_error` and already says "this is not a mixed model." The
   `"failed to construct LMM: "` wrapper conveys nothing new and adds visual
   clutter.

---

## Classification

| Axis | Rating |
|------|--------|
| Panic / crash / silent wrong result | No — clean typed error |
| Worse than lme4 | No — comparable or slightly better (typed class) |
| Actionable hint present | No — neither tool gives one |
| Noisy parent chain | Yes — `doTryCatch()` frame bleeds through |
| Overall verdict | **needs-work** (cosmetic noise, no hint; not a bug) |

The core promise — "no fabricated inference, clear typed diagnostic" — is met.
The `doTryCatch()` parent-chain leakage and missing actionable hint are the two
concrete improvement targets.
