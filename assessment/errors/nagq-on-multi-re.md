# Error-message quality: nAGQ > 1 on multiple / non-scalar random effects

**Scenario:** `nAGQ > 1` requested on a GLMM with either (a) multiple grouping
factors or (b) a single vector (non-scalar) random-effect term.  lme4/glmer
refuses both configurations.

Probe script: `/tmp/nagq_probe.R`
Date: 2026-05-31

---

## lme4 baseline

Both configurations produce the same `simpleError`:

```
nAGQ > 1 is only available for models with a single, scalar random-effects term
```

The message is brief and accurate but gives no actionable fix (no suggestion to
use `nAGQ = 1`, no explanation of what "scalar" means in this context).

---

## mixeff behaviour

### Scenario C — two grouping factors, nAGQ = 5

**Class chain:** `mm_fit_error, mm_condition, rlang_error, error, condition`

**Verbatim message:**
```
failed to fit GLMM: Invalid argument: n_agq = 5 > 1 requires exactly one scalar
random-effects term; this model has 2 term(s) with vsizes [1, 1]
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to fit GLMM: Invalid argument: n_agq = 5 > 1 requires
exactly one scalar random-effects term; this model has 2 term(s) with vsizes [1, 1]
```

The core text — *"n_agq = 5 > 1 requires exactly one scalar random-effects term;
this model has 2 term(s) with vsizes [1, 1]"* — is fully diagnostic: it names
the rule, the requested value, and the actual model structure. The class
`mm_fit_error` is a typed condition that callers can catch programmatically.

The duplicated "Caused by error in `doTryCatch()`" tail is cosmetic noise from
the rlang error chain and does not obscure the core message.

### Scenario D — single vector RE (random slope), nAGQ = 5

**Class chain:** `mm_fit_error, mm_condition, rlang_error, error, condition`

**Verbatim message:**
```
failed to fit GLMM: Invalid argument: n_agq = 5 > 1 requires exactly one scalar
random-effects term; this model has 1 term(s) with vsizes [2]
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to fit GLMM: Invalid argument: n_agq = 5 > 1 requires
exactly one scalar random-effects term; this model has 1 term(s) with vsizes [2]
```

Correctly distinguishes the *non-scalar* case from the *multiple-term* case:
`vsizes [2]` makes clear this is a two-dimensional (intercept + slope) RE block,
not just "too many groups". This is strictly MORE informative than lme4, which
conflates both failure modes in one message.

### Scenario E — single scalar RE, nAGQ = 5 (control / baseline)

**Outcome:** SUCCESS — the model fits without error, confirming nAGQ > 1 is
accepted when the precondition is satisfied. No spurious refusal.

### Scenario F — lmm() with multiple grouping factors (nAGQ N/A)

**Outcome:** SUCCESS — `lmm()` does not accept an `nAGQ` argument; multiple REs
are handled correctly by the LMM path, which is unrelated to quadrature.

---

## Comparison

| Dimension                     | lme4/glmer                                      | mixeff/glmm                                               |
|-------------------------------|------------------------------------------------|-----------------------------------------------------------|
| Error class                   | `simpleError` (untyped)                        | `mm_fit_error` (typed, catchable)                        |
| Names the rule                | yes ("single, scalar random-effects term")     | yes ("exactly one scalar random-effects term")            |
| Reports requested nAGQ value  | no                                             | yes (`n_agq = 5 > 1`)                                    |
| Reports actual model structure| no                                             | yes (`2 term(s) with vsizes [1, 1]` / `vsizes [2]`)      |
| Distinguishes multi-RE vs vector-RE | no (same message for both)              | yes (`vsizes [1, 1]` vs `vsizes [2]`)                    |
| Actionable fix suggested      | no                                             | no (minor gap: could say "use nAGQ = 1")                 |
| Duplicate message tail        | no                                             | yes (cosmetic rlang chain noise, not blocking)            |

---

## Assessment

**Verdict: good** — mixeff's message is **clearer than lme4's** on this scenario.
It is typed (`mm_fit_error`), names the violated constraint, reports the
requested value, reports the actual model topology (number of terms and their
vector sizes), and correctly distinguishes the two failure modes that lme4
conflates into one message.

Two minor improvement opportunities (neither a bug):

1. **No actionable fix suggestion.** Neither lme4 nor mixeff suggests the
   remedy (`nAGQ = 1L`). Adding "use `nAGQ = 1L` to fit this model with the
   Laplace approximation" would be a genuine UX win.

2. **Duplicate rlang chain tail.** The repeated `"Caused by error in
   doTryCatch(): ! mm_fit_error: ..."` line is cosmetic but slightly clutters
   the output. It comes from `mm_abort_from_bridge` wrapping the bridge
   condition inside rlang's chained error format. Not a functional problem.

Neither issue rises to "needs-work" for the parity audit; both are polish items.

**Classification:** `works` (message quality exceeds lme4 baseline)
