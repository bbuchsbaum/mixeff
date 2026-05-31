# Gap Report — Random-Effects Formula Syntax

**Family:** Random-effects formula syntax — `(x|g)`, `(x||g)`, nesting `/`,
interaction `:`, crossing `*`, grouping, `dummy()`, `0+`/`-1` intercept
handling, multiple grouping factors, and the `findbars`/`nobars`/`subbars`/
`expandDoubleVerts` utility surface.

**Date:** 2026-05-31
**Reference:** `assessment/survey/lme4-formula-re-syntax.md`
**mixeff entry points:** `R/parse-formula.R` (`mm_parse_formula`),
`R/compile.R` (`compile_model`), `R/manifest.R` (`mm_formula_manifest`),
`R/fit-lmm.R` (`lmm`). Formula parsing/expansion is delegated entirely to the
Rust engine via `wrap__mm_parse_formula` / `wrap__mm_compile_model_json`; there
is no R-side formula AST.

## Method

All status calls below were confirmed by running `library(mixeff)` against
`mixeff` 0.x as installed, comparing to `lme4` 2.0.1 where relevant. Parse
round-trips used `mm_parse_formula()` (canonical `Display` string); end-to-end
fits used `lmm()` on `sleepstudy` (augmented with synthetic `item`/`site`/`cond`
columns). Scripts: `/tmp/re_test.R`, `/tmp/re_test2.R`, `/tmp/re_test3.R`.

## Summary

The **core lme4 RE grammar is fully and faithfully supported**: scalar/slope
intercepts, explicit/implicit/suppressed intercepts, multi-slope terms,
`||`, nesting `/`, interaction `:`, crossing `*`, and multiple cross-classified
grouping factors all parse to correct canonical forms and fit end-to-end.
Parse-time errors for unsupported constructs are clear and typed
(`mm_formula_error`) — meeting the "clearer errors" bar.

The gaps are concentrated in two areas:

1. **Inline transform / covariance-tag wrappers on RE terms**
   (`factor()` on LHS, `dummy()`, and the lme4 ≥2.0 `us()`/`diag()`/`cs()`/
   `ar1()` tags). The engine rejects all `foo(...)` constructs outside a small
   stateless subset. `ar1()` (and structured residual covariance generally) is
   out-of-scope by PRD §3; `factor()`/`dummy()`/`diag()`/`cs()` are at least
   **partial** gaps that a real lme4 user writing RE terms will hit.

2. **The programmer-facing formula-utility surface**
   (`findbars`/`nobars`/`subbars`/`expandDoubleVerts`/`isNested`/`mkReTrms`/
   `splitForm`/`reOnly`/`randint`/`noSpecials`) is entirely absent from the
   mixeff R namespace. These are public reformulas/lme4 functions that
   downstream tooling and power users rely on.

## Capability table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `(1\|g)` scalar random intercept | works | works | — | `mm_parse_formula("y~x+(1\|g)")` → `y ~ 1 + x + (1 \| g)`; `lmm` fits. |
| `(x\|g)` correlated intercept+slope (default `us`) | works | works | — | `(x\|g)` → canonical `(1 + x \| g)`, `cov="full"`, theta=3; `lmm(Reaction~Days+(Days\|Subject))` fits. |
| `(1+x\|g)` explicit intercept | works | works | — | round-trips identically to `(x\|g)`. |
| `(0+x\|g)` / `(-1+x\|g)` slope, no intercept | works | works | — | both → `(0 + x \| g)`; `lmm(...+(0+Days\|Subject))` fits, prints `random_slope_without_intercept` note. |
| `(a+b\|g)` multiple correlated slopes | works | works | — | `(a+b\|g)` → `(1 + a + b \| g)`. |
| `(x\|\|g)` uncorrelated (continuous) | works | works | — | `(Days\|\|Subject)` splits to `(1\|Subject)+(0+Days\|Subject)`; VarCorr std devs match lme4 to ~3 dp (Int 25.0513, Days 5.9885 vs 5.9882; resid 25.565). Cross-block note printed. |
| `(a+b\|\|g)` multi-slope uncorrelated | works | works | — | `(a+b\|\|g)` → `(1 + a + b \|\| g)` canonical. |
| `(x\|\|g)` with a **factor** predictor | works (better than lme4) | works | — | `(cond\|\|Subject)`: mixeff produces a single block diagonalized over contrast columns (`cond:b`,`cond:c`, +0 corr) — i.e. the *correct* zero-correlation-over-contrasts result that lme4's default `"split"` method gets wrong (lme4 needs `diag_special`). No defect. |
| `(x\|g)+(1\|g)` two blocks, same g | works | works | — | `(0+Days\|Subject)+(1\|Subject)` fits; prints "separate random-effect blocks fix the covariance … to zero." |
| `(1\|a/b)` nesting | works | works | — | → `(1\|a) + (1\|a:b)`; `lmm(...+(1\|site/Subject))` fits with `syntax_expansion` note. |
| `(1\|a/b/c)` 3-level nesting | works | works | — | → `(1\|a) + (1\|a:b) + (1\|a:b:c)`. |
| `(1\|a:b)` interaction grouping | works | works | — | → `(1\|a:b)`, single block; fits. Note: lme4 reference renders nesting as `b:a`; mixeff uses `a:b` ordering — semantically equivalent factor combination. |
| `(1\|a*b)` full crossing | works | works | — | → `(1\|a)+(1\|b)+(1\|a:b)`; fits, prints `crossing_likely_unintended` guidance. |
| Multiple grouping factors `(1\|s)+(1\|item)` (cross-classified) | works | works | — | fits; independent scalar blocks, prints independence covariance note. |
| `y ~ 0 + x + (1\|g)` fixed-intercept suppression | works | works | — | `0+`/`-1+` both drop the fixed `1`: canonical `y ~ x + (1 \| g)`. |
| Clear errors on unsupported syntax | works | works | — | unsupported `foo(...)` raises typed `mm_formula_error` naming the construct, position, and the allowed stateless subset — meets the "clearer errors" bar. |
| `factor(grp)` on LHS of `\|` | errors (precompute workaround) | partial | major | `(factor(grp)\|g)` → `mm_formula_error` ("`factor(...)` not in stateless subset"). Workaround: precompute indicator columns and use `(0+col\|g)` — confirmed to fit. A real lme4 user writing `(factor(grp)\|g)` is rejected, not silently handled. Error is clear but no automatic in-formula support. |
| `dummy(f, levels)` selective dummy coding | errors (precompute workaround) | partial | major | `(0+dummy(Sex,"Female")\|g)` → `mm_formula_error` ("`dummy(...)` not in stateless subset"). Same precompute workaround applies. lme4-documented idiom is unavailable in-formula. |
| `us(x\|g)` explicit unstructured tag | errors | partial | minor | `us(x\|g)` → `mm_formula_error`. The *structure itself* is fully supported via bare `(x\|g)` (`cov="full"`); only the lme4≥2.0 named-tag spelling is unrecognized. User reaches the same model without the tag. |
| `diag(x\|g)` / `diag(...,hom=TRUE)` tag | errors (synonym path exists) | partial | major | `diag(1+a+b\|g)` → `mm_formula_error`. The diagonal structure is reachable via the `\|\|` synonym (mixeff `cov="diag"`/split blocks), which mixeff documents as "same model, different font" (PRD §9.5). But the explicit `diag()` tag and `hom=TRUE` homogeneous-variance variant are not parseable. |
| `cs(x\|g)` compound-symmetric tag | in-scope-missing | in-scope-missing | major | `cs(x\|g)` → `mm_formula_error`. No `cov="cs"`/compound-symmetric family in mixeff's covariance vocabulary (`full`/`diag`/`scalar` only). Not listed in PRD §3 non-goals, so not deferred-by-design; genuinely absent. A user wanting CS covariance has no path. |
| `ar1(x\|g)` AR(1) residual/structured tag | out-of-scope-by-design | out-of-scope-by-design | minor | `ar1(x\|g)` → `mm_formula_error`. PRD §3 explicitly defers "residual covariance structures (AR(1), spatial)" to post-v0. Classify as deferred, not a defect. |
| `findbars(term)` | absent | in-scope-missing | major | not in `getNamespaceExports("mixeff")`. No R-side bar extraction; all parsing is in Rust. Programmatic formula manipulation that lme4/reformulas users expect is unavailable. |
| `findbars_x(...)` (with expansion control) | absent | in-scope-missing | minor | absent; mixeff has no R formula AST to operate on. |
| `nobars(term)` / `nobars_` | absent | in-scope-missing | major | absent. `nobars(y~x+(Days\|Subject))` → fixed-only formula is a very common downstream call (e.g. building fixed design); not exposed. |
| `subbars(term)` | absent | in-scope-missing | minor | absent. Used to build `model.frame()`-safe formulas; mixeff builds its own model frame in `compile_model` so internal need is met, but the public helper is missing. |
| `expandDoubleVerts(term)` | partial (internal only) | partial | minor | No exported `expandDoubleVerts`. The expansion **is performed** inside the Rust parser (`(Days\|\|Subject)` canonicalizes to split blocks, visible via `mm_parse_formula`), so the behavior exists but is not callable as a formula→formula utility. |
| `isNested(f1,f2)` | absent | in-scope-missing | minor | absent. mixeff does data-driven nesting/support analysis inside the design audit (e.g. `repeated_unit_unmodeled` notes), but no standalone predicate. |
| `mkReTrms(bars,fr,...)` | out-of-scope-by-design | out-of-scope-by-design | minor | absent. lme4's sparse-matrix builder is replaced by the Rust engine's compiled artifact; exposing `Zt`/`Lambdat`/`cnms` is the `getME()` compat surface (separate family), not a formula-syntax obligation. PRD §8.2 lists `getME` read-only compat, not `mkReTrms`. |
| `splitForm(formula,...)` | absent | in-scope-missing | minor | absent. Special-term parser; mixeff parses specials in Rust and rejects unknown wrappers. No R-callable splitForm. |
| `reOnly(f, response, ...)` | absent | in-scope-missing | minor | absent. RE-only extraction not exposed. |
| `expandAllGrpVar(bb)` | absent | out-of-scope-by-design | cosmetic | internal reformulas helper; expansion happens in Rust. No user-facing obligation. |
| `randint(formula)` | absent | in-scope-missing | cosmetic | absent. Minor convenience; rarely used directly. |
| `noSpecials(term,...)` | absent | in-scope-missing | minor | absent companion to `splitForm`. |
| `getDoublevertDefault()` / `options(lme4.doublevert.default)` | n/a (different model) | out-of-scope-by-design | minor | mixeff has no split-vs-diag_special toggle: for factor predictors it already does the diag-special (correct) thing (see `(cond\|\|Subject)` row), so the lme4 option that exists to fix the buggy default is moot. No equivalent needed. |
| `lmerTest::ranova()` / `rand()` (RE-term LRT reduction) | covered elsewhere | (out of family) | — | Random-effect LRT reduction lives in the inference family (`test_random_effect`, see `man/test_random_effect.Rd` and `tests/testthat/test-boundary-lrt.R`), not in this formula-syntax family. Not scored here. |

## Notes and caveats

- **No silent surgery confirmed.** Every expansion (`/`, `*`, `\|\|`, split
  blocks) is reported back to the user through `explain_model()`/audit notes
  (`syntax_expansion`, `crossing_likely_unintended`, cross-block covariance
  notes). This is the PRD §9.5 contract and it holds in practice — the wrapper
  does *not* hide the transformation.

- **Manifest under-declares the surface.** `mm_formula_manifest()$formula_features`
  lists `random_term_forms` of 8 entries and `transformations` of only
  `implicit_intercept`, `nested_grouping_expansion`, `interaction_grouping`. It
  does **not** advertise `crossing_expansion` or `double_bar_split` even though
  both demonstrably work. Minor: the machine-readable capability record lags the
  actual (broader) capability. Worth reconciling so `audit()` reports the true
  surface.

- **`cs()` is the one clear in-scope-missing covariance family.** `us`(=full),
  `diag`, `scalar` are all reachable; compound-symmetric has no path and is not
  named in PRD §3 non-goals. This is the strongest covariance-family gap a user
  could hit.

- **Inline `factor()`/`dummy()` are "partial" not "missing"** because a correct
  model is reachable by precomputing indicator columns (confirmed to fit), and
  the rejection error is explicit. But the ergonomic lme4 idioms are
  unavailable, so a real user is impeded — major.
