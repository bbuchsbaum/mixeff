# Error-message quality probe: malformed-formula

**Scenario:** Syntactically broken or semantically invalid formula strings and objects
passed to `lmm()` and `glmm()`.

**Probe script:** `assessment/errors/probe-malformed-formula.R`
**Date run:** 2026-05-31

---

## Raw output (verbatim)

### Scenario 1 — Empty formula string (`""`)

| caller | message | class |
|---|---|---|
| lme4::lmer | `parsing result not of length one, but 0` | simpleError |
| lme4::glmer | `parsing result not of length one, but 0` | simpleError |
| mixeff::lmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |
| mixeff::glmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |

### Scenario 2 — Garbled string (double tilde: `"y ~~ x + (1 | subject)"`)

| caller | message | class |
|---|---|---|
| lme4::lmer | `invalid model formula` | simpleError |
| lme4::glmer | `invalid model formula` | simpleError |
| mixeff::lmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |
| mixeff::glmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |

### Scenario 3 — Missing `|` in random term: `"y ~ x + (1  subject)"`

| caller | message | class |
|---|---|---|
| lme4::lmer | `<text>:1:13: unexpected symbol\n1: y ~ x + (1  subject\n                ^` | simpleError |
| lme4::glmer | same as lmer | simpleError |
| mixeff::lmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |
| mixeff::glmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |

### Scenario 4 — Empty grouping factor: `"y ~ x + (1 | )"`

| caller | message | class |
|---|---|---|
| lme4::lmer | `<text>:1:14: unexpected ')'\n1: y ~ x + (1 | )\n                 ^` | simpleError |
| lme4::glmer | same as lmer | simpleError |
| mixeff::lmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |
| mixeff::glmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |

### Scenario 5 — Valid R formula but no random effects: `y ~ x`

| caller | message | class |
|---|---|---|
| lme4::lmer | `No random effects terms specified in formula` | simpleError |
| lme4::glmer | `No random effects terms specified in formula` | simpleError |
| mixeff::lmm | `failed to construct LMM: No random effects in formula: this is not a mixed model\nCaused by error in \`doTryCatch()\`:\n! mm_fit_error: failed to construct LMM: No random effects in formula: this is not a mixed model` | mm_fit_error |
| mixeff::glmm | `failed to construct GLMM: No random effects in formula: this is not a mixed model\nCaused by error in \`doTryCatch()\`:\n! mm_fit_error: failed to construct GLMM: No random effects in formula: this is not a mixed model` | mm_fit_error |

### Scenario 6 — Non-formula, non-character argument (numeric `42`)

| caller | message | class |
|---|---|---|
| lme4::lmer | `invalid formula` | simpleError |
| lme4::glmer | `invalid formula` | simpleError |
| mixeff::lmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |
| mixeff::glmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |

### Scenario 7 — Formula as `NA`

| caller | message | class |
|---|---|---|
| lme4::lmer | `invalid formula` | simpleError |
| mixeff::lmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |
| mixeff::glmm | `` `formula` must be a two-sided R formula (lhs ~ rhs). `` | mm_formula_error |

---

## Assessment

### Scenarios 1–4 and 6–7: GOOD — clearly better than lme4

For every purely syntactic malformation (empty string, double tilde, missing pipe,
empty grouping factor, wrong type, NA), mixeff produces:

- A typed `mm_formula_error` condition (catchable by class).
- A plain-English message: `` `formula` must be a two-sided R formula (lhs ~ rhs). ``
- No stack trace, no Rust internals, no cryptic position markers.

lme4's corresponding messages (`"parsing result not of length one, but 0"`,
`"invalid model formula"`, `"<text>:1:13: unexpected symbol"`) are parse-level
fragments that expose R parser internals and give the caller nothing actionable.

mixeff's message is more actionable: it names the expected shape (`lhs ~ rhs`),
and the typed class lets callers catch-and-handle without string matching.

**One note:** Scenarios 3 and 4 involve strings that are syntactically unparseable
as R formulas (`y ~ x + (1  subject)`, `y ~ x + (1 | )`). R's `as.formula()`
would fail at parse time, but lme4 emits R parser position markers because it
calls `parse()` internally. mixeff's `compile_model()` path goes through
`mm_coerce_formula_string()` first; when passed a character string, it validates
structure without calling `parse()`, so it catches the error at a higher level
and emits the same typed message. The Rust parser would give a more specific
"missing grouping variable" message for Scenario 3 if the string were parseable;
instead R-level coercion fires first. This is acceptable: the user still gets a
typed, catchable, plain-English error.

### Scenario 5 — No random effects: NEEDS-WORK

**Issue: internal infrastructure leaks into `conditionMessage()`.**

The full text of `conditionMessage(e)` for `lmm(y ~ x, df)` is:

```
failed to construct LMM: No random effects in formula: this is not a mixed model
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct LMM: No random effects in formula: this is not a mixed model
```

Problems:

1. The first line is the clean, user-facing message.
2. Lines 2–3 are rlang's chained-error formatting leaking through: `Caused by
   error in \`doTryCatch()\`:` is a call-stack frame from the bridge tryCatch,
   not meaningful to users.
3. Line 3 repeats the raw Rust-tagged string (`mm_fit_error: failed to construct
   LMM: ...`) — the `mm_<tag>: ` prefix that `mm_split_tagged_error()` is
   supposed to strip is still present in the chained parent condition and
   surfaced by rlang's formatter.

By comparison, lme4's message is a clean one-liner: `No random effects terms
specified in formula`.

**Root cause:** `mm_abort_from_bridge()` sets `parent = cnd` (the raw bridge
error). rlang appends the parent condition's formatted message to the child
message, including the raw `! mm_fit_error: ...` line from the bridge error.
The `mm_split_tagged_error()` stripping only affects the new condition's
`message` field, not the chained parent.

**Fix direction:** Either pass `parent = NULL` when the child message already
contains the full user-facing content, or strip the tag from the parent
condition's message before chaining it. The parent condition carries no
additional diagnostic value for this error path.

**Class is correct:** `mm_fit_error` is the right typed condition for this case
(distinct from `mm_formula_error`, which is for syntactic/structural formula
problems). The semantic diagnosis ("no random effects") is accurate.

lme4 comparison: lme4 says "No random effects terms specified in formula".
mixeff says "No random effects in formula: this is not a mixed model" — the
additional clause ("this is not a mixed model") is a useful clarification. The
core content is comparable-to-better, but the trailing infrastructure noise
degrades the overall quality.

---

## Verdict summary

| Scenario | mixeff vs lme4 | verdict |
|---|---|---|
| Empty string | Clearly better (typed, plain English vs cryptic length error) | good |
| Garbled string | Clearly better (typed vs "invalid model formula") | good |
| Missing pipe in random term | Better (typed vs R parser position marker) | good |
| Empty grouping factor | Better (typed vs R parser position marker) | good |
| No random effects (`y ~ x`) | Core message accurate; `conditionMessage()` leaks `doTryCatch()` + raw Rust tag | needs-work |
| Non-formula argument | Better (typed vs "invalid formula") | good |
| NA as formula | Better (typed vs "invalid formula") | good |

**Overall verdict: needs-work** (one specific case — no-random-effects —
has a `conditionMessage()` that leaks internal call-frame text and the raw
Rust tag; all syntactic/structural malformations are handled clearly).

**Classification:** `in-scope-missing` for the leaking parent condition in
the no-random-effects path. The fix is small and localized to
`mm_abort_from_bridge()` or its call site in `lmm()`/`glmm()`.
