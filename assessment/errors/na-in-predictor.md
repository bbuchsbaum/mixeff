# Error-message quality probe: `na-in-predictor`

**Scenario:** NAs in a predictor (and response, and grouping variable) with
default `na.action`.

**Probe script:** `assessment/errors/probe-na-in-predictor.R`  
**Date run:** 2026-05-31

---

## Scenario A — NA in fixed-effect predictor (`x` has 3 NAs)

### lme4::lmer
```
(succeeded — lme4 silently dropped NA rows via na.omit)
nobs: 37 (expected 37 after NA drop)
```
lme4 applies `na.action = na.omit` silently. No message, no warning, no
indication that 3 rows were discarded. The user has no way to know from the fit
object alone that data were dropped without inspecting `attr(model.frame, "na.action")`.

### lme4::glmer (binomial)
```
(succeeded — lme4 silently dropped NA rows via na.omit)
nobs: 37 (expected 37 after NA drop)
```
Same silent-drop behaviour.

### mixeff::compile_model
```
class: mm_data_error, mm_condition, rlang_error, error, condition
message:
  Missing values in design variable(s): `x` (3 NA). mixeff requires
  complete cases; pass na.omit(data) explicitly before fitting.
```

### mixeff::lmm
```
class: mm_data_error, mm_condition, rlang_error, error, condition
message:
  Missing values in design variable(s): `x` (3 NA). mixeff requires
  complete cases; pass na.omit(data) explicitly before fitting.
```

### mixeff::glmm (binomial)
```
class: mm_data_error, mm_condition, rlang_error, error, condition
message:
  Missing values in design variable(s): `x` (3 NA). mixeff requires
  complete cases; pass na.omit(data) explicitly before fitting.
```

---

## Scenario B — NA in response variable (`y` has 2 NAs)

### lme4::lmer
```
(succeeded — lme4 silently dropped NA rows)
nobs: 38 (expected 38 after NA drop)
```
Silent drop, same as Scenario A.

### mixeff::lmm
```
class: mm_data_error, mm_condition, rlang_error, error, condition
message:
  Missing values in design variable(s): `y` (2 NA). mixeff requires
  complete cases; pass na.omit(data) explicitly before fitting.
```
mixeff catches NAs in the response too — `mm_check_no_na` checks all variables
named in the formula including the LHS.

---

## Scenario C — NA in grouping variable (`g` has 2 NAs)

### lme4::lmer
```
(succeeded — lme4 silently dropped NA rows)
nobs: 38 (expected 38 after NA drop)
```

### mixeff::lmm
```
class: mm_data_error, mm_condition, rlang_error, error, condition
message:
  Missing values in design variable(s): `g` (2 NA). mixeff requires
  complete cases; pass na.omit(data) explicitly before fitting.
```
Grouping variable NAs are caught at the same stage, with the same typed
condition.

---

## Scenario D — `glmm()` custom `na.action` refusal

### mixeff::glmm (na.action = na.pass)
```
class: mm_fit_error, mm_condition, rlang_error, error, condition
message:
  `random`, `weights`, `subset`, custom `na.action`, and `contrasts` are
  reserved for the fitted GLMM bridge.
```
The refusal message is accurate but lumps `na.action` in with four other
reserved arguments. A caller who only passed `na.action = na.pass` gets a
message that looks like a generic "reserved features" wall rather than a
targeted "na.action is not yet implemented" note. This is a minor usability
gap, not a bug.

---

## Assessment

### Verdict: `good` (with one minor `needs-work` sub-item)

| Dimension | lme4 | mixeff |
|---|---|---|
| NA in predictor | Silent success — rows dropped without any message | Typed `mm_data_error`; names column and NA count; gives exact remediation |
| NA in response | Silent success | Typed `mm_data_error`; same quality |
| NA in grouping var | Silent success | Typed `mm_data_error`; same quality |
| Typed / catchable | No condition raised | `mm_data_error` inheriting `mm_condition`; programmatically catchable |
| Actionable remediation | None | "pass `na.omit(data)` explicitly before fitting" — exact fix given |
| Fabrication risk | N/A (succeeds silently) | None — refuses rather than guessing |

### mixeff's message quality
The message is clear, actionable, and non-inscrutable:
- Names the offending variable(s) with backtick quoting.
- Gives the exact NA count per column.
- Provides the exact corrective action (`na.omit(data)`).
- Raises a typed condition (`mm_data_error`) that callers can `tryCatch` generically
  via `mm_condition` or specifically via `mm_data_error`.
- Does not panic, produce a stack trace, or return a silently wrong result.

This is strictly better than lme4's behaviour (silent row-drop with no
user-facing signal).

### One minor `needs-work` item (Scenario D)
`glmm(na.action = na.pass)` raises `mm_fit_error` with a message that lists
five reserved arguments together. A caller who only passed `na.action` sees
"`random`, `weights`, `subset`, custom `na.action`, and `contrasts` are
reserved..." — it is technically accurate but not targeted. A focused message
such as "Custom `na.action` is not yet supported; mixeff always requires
complete cases. Use `na.omit(data)` before calling `glmm()`." would be more
helpful. This is cosmetic and does not rise to a bug.

### Classification
`works` / `good` for the primary NA-in-predictor path.  
`needs-work` (cosmetic) for the custom `na.action` refusal wording in `glmm()`.

### PRD §3 non-goal check
No out-of-scope issues flagged. The no-silent-surgery principle (PRD §8.1) is
correctly enforced: mixeff refuses rather than dropping rows silently.
