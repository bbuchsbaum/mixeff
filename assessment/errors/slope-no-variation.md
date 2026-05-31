# Error-message quality probe: slope-no-variation

**Scenario:** A random slope is requested for a predictor that does not vary
within any grouping level — a between-subjects variable.  The slope is
structurally unidentifiable from the data.

**Dataset:** 20 subjects × 4 observations; `treatment` is a between-subjects
factor (every subject sees only one level), so `(1 + treatment | subject)` has
zero within-subject variation to estimate the slope from.  Also tested with a
numeric between-subjects predictor (`score = as.numeric(subject)`).

Probe script: `assessment/errors/probe-slope-no-variation.R`  
Run date: 2026-05-31

---

## Verbatim output by surface

### 1. lme4::lmer — factor between-subjects slope

```
boundary (singular) fit: see help('isSingular')
  error  :  (none)
  warning:  (none)   [warning was printed directly, not catchable via tryCatch]
  result :  lmerMod
```

lme4 **fits silently** and prints a bare singularity advisory. The message
`boundary (singular) fit: see help('isSingular')` is emitted directly to the
console (not as an R warning condition catchable by `tryCatch`), gives no
indication that the slope is between-subjects or structurally unidentifiable,
and does not name the problematic term.

### 2. lme4::glmer — factor between-subjects slope (binary)

```
boundary (singular) fit: see help('isSingular')
  error  :  (none)
  warning:  (none)
  result :  glmerMod
```

Identical behaviour to lmer: silent fit + same opaque singularity note.

### 3. mixeff::lmm — factor between-subjects slope

```
  error  :  (none)
  warning:  (none)
  result :  mm_lmm
  diagnostics:
                    code severity         stage            affected_terms
   formula_canonicalized     info  design_audit (1 + treatment | subject)
      structural_refusal     info  design_audit (1 + treatment | subject)
      boundary_parameter     info certification (1 + treatment | subject)
      covariance_reduced     info certification (1 + treatment | subject)

Messages:
  formula_canonicalized: random-effect basis was expanded into optimizer columns
  structural_refusal: `treatment: B` does not vary within `subject`, so a
        `subject`-level `treatment: B` slope cannot be estimated from this design
  boundary_parameter: standard deviation for treatment: B in (1 + treatment |
        subject) is on its lower bound
  covariance_reduced: fitted covariance for (1 + treatment | subject) has
        effective rank 1 of requested rank 2
  fit_status: converged_reduced_rank
```

### 4. mixeff::glmm — factor between-subjects slope (binary)

```
  error  :  (none)
  warning:  (none)
  result :  mm_glmm
  diagnostics:
                    code severity         stage            affected_terms
   formula_canonicalized     info  design_audit (1 + treatment | subject)
      structural_refusal     info  design_audit (1 + treatment | subject)
      boundary_parameter     info certification    covariance parameter 1
      boundary_parameter     info certification
  
Messages:
  formula_canonicalized: random-effect basis was expanded into optimizer columns
  structural_refusal: `treatment: B` does not vary within `subject`, so a
        `subject`-level `treatment: B` slope cannot be estimated from this design
  boundary_parameter: covariance parameter 1 is on its lower bound
  boundary_parameter: GLMM covariance state classified as ValidZeroVariance
  fit_status: converged_boundary
```

### 5. mixeff::lmm — numeric between-subjects slope

```
  error  :  (none)
  warning:  (none)
  result :  mm_lmm
  diagnostics:
                       code severity         stage        affected_terms
         structural_refusal     info  design_audit (1 + score | subject)
   optimizer_nonconvergence  warning certification
         boundary_parameter     info certification (1 + score | subject)

Messages:
  structural_refusal: `score` does not vary within `subject`, so a
        `subject`-level `score` slope cannot be estimated from this design
  optimizer_nonconvergence: optimizer stopped before an acceptable convergence
        criterion with return code 'MAXEVAL_REACHED'
  boundary_parameter: standard deviation for score in (1 + score | subject)
        is on its lower bound
```

### 6. lme4::lmer — numeric between-subjects slope

```
boundary (singular) fit: see help('isSingular')
  error  :  (none)
  warning:  (none)
  result :  lmerMod
```

Same opaque singularity note, no mention of between-subjects structure.

---

## Assessment

### lme4 behaviour

lme4 silently returns a fitted object with a singularity boundary and emits
only `boundary (singular) fit: see help('isSingular')` — directly to the
console, not as a catchable R condition.  The message:

- Does not name the slope or predictor involved.
- Does not explain *why* the singularity occurred (between-subjects structure,
  zero within-group variation).
- Does not indicate the fit is untrustworthy beyond a pointer to `?isSingular`.
- Is not actionable without the user independently diagnosing the design.

### mixeff behaviour

mixeff also does **not** error or refuse the fit in this scenario; it fits and
returns an `mm_lmm` / `mm_glmm` object.  However it emits two clearly
named diagnostics at `design_audit` stage:

1. **`structural_refusal`** (severity: info) — names the exact term and
   grouping: "`treatment: B` does not vary within `subject`, so a
   `subject`-level `treatment: B` slope cannot be estimated from this design."
   This is a precise, actionable message that locates the problem structurally.

2. **`covariance_reduced`** (for lmm) / **`boundary_parameter`** (for both) —
   post-fit evidence that the optimizer also hit the boundary, consistent with
   the pre-fit audit.

The `fit_status` field is `converged_reduced_rank` (lmm) or
`converged_boundary` (glmm), not the bare "converged", so downstream code can
gate on this.

### Comparative verdict

| Criterion                          | lme4           | mixeff         |
|------------------------------------|----------------|----------------|
| Names the problematic term         | No             | Yes            |
| Explains why (between-subjects)    | No             | Yes            |
| Catchable as a structured condition| No (console)   | Yes (diagnostics table) |
| fit_status reflects the problem    | No             | Yes            |
| Refuses/errors vs silent           | Silent fit     | Silent fit     |

mixeff's message quality is **materially clearer than lme4's** on this
scenario.  The `structural_refusal` diagnostic names the term, names the group,
and states the structural reason in plain English.

### Remaining gap (needs-work, not a bug)

mixeff does **not** raise an R error or warning at fit time; it silently fits
and surfaces the structural refusal only via `diagnostics()`.  A user who does
not inspect `diagnostics()` or `fit_status()` will receive an `mm_lmm` object
with a `converged_reduced_rank` status and might not notice the issue.

Three improvements would make this scenario fully "best-in-class":

1. **Emit an R warning** (a catchable `mm_fit_warning` condition) at the end of
   `lmm()` / `glmm()` whenever any `structural_refusal` diagnostic is present
   in the artifact, so the user sees it at the REPL without inspecting
   diagnostics manually.
2. **Surface the structural refusal in the default `print.mm_lmm()` output**
   rather than only in `diagnostics()`.
3. For the numeric slope case the optimizer also hits `MAXEVAL_REACHED`
   (`optimizer_nonconvergence`) in addition to the structural refusal.  The
   nonconvergence message is a red herring caused by the structurally
   unidentifiable term; ideally the engine would short-circuit and not try to
   fit a slope the audit already declared unestimable.

**Classification:** `needs-work` — mixeff is clearly better than lme4 on this
scenario but does not yet surface the structural refusal to the user without
an explicit `diagnostics()` call.
