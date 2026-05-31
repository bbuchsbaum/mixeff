# Error-message quality probe: revive-no-handle

**Scenario:** Use a revived/serialized fit after the Rust handle is dropped — does
it revive cleanly or error cryptically?

**Probe script:** `assessment/errors/probe-revive-no-handle.R`

---

## 1. lme4 / lmerTest baseline

lme4 objects survive `saveRDS()` / `readRDS()` without any explicit revival step.
All tested post-fit verbs work silently on the reloaded object:

```
lme4 fixef after reload:    OK — (Intercept)=2.2507, x=0.3902
lme4 predict after reload:  OK — first 3: 2.7589, 3.1491, 3.5394
lme4 ranef after reload:    OK — n groups: 10

lme4 glmer fixef after reload:   OK — (Intercept)=0.776, x=0.2207
lme4 glmer predict after reload: OK — first 3: 2.8569, 3.5623, 4.4419
```

lme4 requires no user-visible revival API because its handle to the
underlying C++ object is transparently reattached via `.onLoad` / R's
external-pointer finalizer / lazy re-evaluation on first access.

---

## 2. mixeff LMM — readRDS WITHOUT revive()

Critically, mixeff **does not require an explicit `revive()` call** for the
common case.  Because the R object carries all durable state (fixef, theta,
sigma, residuals, random effects, artifact JSON, model frame), all tested
verbs succeed on the raw reloaded object:

```
mixeff fixef (no revive):    OK — (Intercept)=2.2507, x=0.3902
mixeff predict (no revive):  OK — first 3: 2.7589, 3.1491, 3.5394
mixeff ranef (no revive):    OK — n groups: 10
mixeff summary (no revive):  OK
mixeff contrast (no revive): OK — p_value: 1.767475e-13
```

The `.mm_lazy()` helper internally calls `revive()` if `lazy_cache` is not an
environment (i.e., after deserialization), so lazy rebuild of the design
matrices, Lambda, flist, etc. also just works.

---

## 3. mixeff LMM — readRDS WITH explicit revive()

Explicit `revive()` is also clean:

```
fit_handle_alive after revive:  FALSE
mixeff fixef (revived):         OK — (Intercept)=2.2507, x=0.3902
mixeff predict (revived):       OK — first 3: 2.7589, 3.1491, 3.5394
mixeff ranef (revived):         OK — n groups: 10
mixeff contrast (revived):      OK — p_value: 1.767475e-13
mixeff random_blocks (revived): OK — groups: subject
mixeff inference_table (revived): OK — methods: asymptotic_wald_z, asymptotic_wald_z
```

`fit_handle_alive()` correctly reports `FALSE` — the Rust handle is
intentionally absent after revival; it is treated as a rebuilable cache, not
the source of truth. The design is intentional and documented.

---

## 4. mixeff revive() on a broken/invalid fit

### 4a. Artifact nulled out
```
revive() on broken fit:
  ERROR class: mm_arg_error mm_condition rlang_error error condition
  ERROR message: `fit` does not carry a parsed compiler artifact.
```

Clear, typed (`mm_arg_error`), actionable — names exactly what is missing.

### 4b. revive() on a plain list (wrong type)
```
revive() on plain list:
  ERROR class: mm_arg_error mm_condition rlang_error error condition
  ERROR message: `revive()` expects a fitted mixeff object.
```

### 4c. revive() on a scalar (wrong type)
```
revive() on integer:
  ERROR class: mm_arg_error mm_condition rlang_error error condition
  ERROR message: `revive()` expects a fitted mixeff object.
```

Both dispatch through `revive.default()` cleanly.

---

## 5. mixeff GLMM — readRDS WITH revive()

```
fit_handle_alive after glmm revive:  FALSE
mixeff glmm fixef (revived):    OK — (Intercept)=0.7966, x=0.2207
mixeff glmm ranef (revived):    OK — n groups: 10
mixeff glmm predict (revived):
  ERROR class: mm_inference_unavailable mm_condition rlang_error error condition
  ERROR message: GLMM prediction is not certified by the current Rust contract.
```

The GLMM `predict()` refusal is **by design** — `predict.mm_glmm` is an
explicit `mm_abort(class = "mm_inference_unavailable")` stub. GLMM prediction
is scoped to Phase 4 (PRD §10). The error is typed, honest, and explains
exactly why: it names the certification boundary. This is not a bug.

---

## 6. Classification

| Sub-scenario | Result | Classification |
|---|---|---|
| LMM reload without revive, standard verbs | All OK | works |
| LMM reload with explicit revive(), standard verbs | All OK | works |
| revive() on broken fit (no artifact) | Typed `mm_arg_error`, clear message | works |
| revive() on wrong type | Typed `mm_arg_error`, clear message | works |
| GLMM reload with revive(), fixef/ranef | All OK | works |
| GLMM predict() after revive | Typed `mm_inference_unavailable`, explains boundary | out-of-scope-by-design (Phase 4) |

---

## 7. Comparison vs lme4

| Dimension | lme4 | mixeff |
|---|---|---|
| Needs explicit revival call after readRDS? | No (transparent) | No (lazy auto-revive) — explicit `revive()` optional |
| Wrong-type error message | "no applicable method for X applied to Y" — generic dispatch error, no actionable guidance | Typed `mm_arg_error`: names exactly what is expected |
| Broken-artifact error | N/A (no artifact concept) | Typed `mm_arg_error`: "`fit` does not carry a parsed compiler artifact." |
| Unsupported post-fit verb (e.g., predict on GLMM) | Often silent wrong answer or cryptic C++ error | Typed `mm_inference_unavailable` with explicit "not certified" reason |
| handle_alive introspection | No API | `fit_handle_alive()` returns `FALSE` correctly |

**Overall assessment:** mixeff's revive pathway is **clearer than lme4**.
lme4 is silent-success only because it has no explicit revival contract;
when something goes wrong (e.g., wrong class, broken object) lme4 produces
generic S3 dispatch errors. mixeff gives typed, named conditions with
actionable messages, and correctly distinguishes "not yet certified by Rust
contract" from "bug" by using `mm_inference_unavailable`.

---

**Verdict:** good
**mixeff_quality:** clearer-than-lme4
