# Gap Report — LMM Fitting & Control

Family: **LMM fitting & control** (lmer REML vs ML, lmerControl, optimizer
choice, weights, offset, subset, na.action, contrasts arg, start, devFunOnly)

Reference: `assessment/survey/lme4-lmm-fit-control.md`
Standard: "everything lme4 does, faster, with clearer errors."
Date: 2026-05-31 | lme4 2.0.1, lmerTest 3.2.1, mixeff (installed)

mixeff entry point for this family is `lmm()`:

```r
lmm(formula, data, REML = TRUE, weights = NULL, control = mm_control())
mm_control(verbose = 0L)   # only `verbose`
```

(`R/fit-lmm.R:38`, `R/mm-control.R:14`.) The Rust FFI `mm_fit_lmm_json`
(`src/rust/src/lib.rs:345`) takes only `formula_string`, `REML`, the
column-oriented model frame, `weights`, and `control_json`. The control JSON is
parsed into a binding named `_control` (underscore = unused;
`src/rust/src/lib.rs:355`), so beyond R-side `verbose` it is effectively a
no-op passthrough today.

---

## Capability table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `REML = TRUE/FALSE` (REML vs ML fit) | Works; values match lme4 to display precision | works | — | `lmm(y~x+(1\|g), REML=TRUE)` logLik -77.52001 == `lmer` -77.52001; `REML=FALSE` logLik -75.15027 distinct. fixef/sigma bit-match lme4. |
| `weights` (prior case weights) | Works; matches lme4 | works | — | `lmm(..., weights=w)` fixef 0.1178093/-0.008828 and sigma 0.9877065 all match `lmer(..., weights=w)`. Validated finite/positive/length in `mm_lmm_weights` (`R/fit-lmm.R:126`). |
| Default factor contrasts (treatment coding) in formula | Works | works | — | `lmm(y~fx+(1\|g))` yields `(Intercept)`, `fx: b`, `fx: c` (treatment coding). Factors handled via the compiler/model-frame path. |
| `contrasts =` argument (Helmert/sum/custom without altering globals) | Absent — `lmm()` has no `contrasts` arg; passing it errors `unused argument`. User must precompute contrast columns by hand. The in-formula message explicitly bans `factor()` and other stateful transforms. | in-scope-missing | major | `lmm(y~factor(x)+(1\|g), contrasts=list(...))` -> `unused argument (contrasts = ...)`. Real lme4 users routinely pass `contrasts=` for sum-to-zero/Helmert coding; here the only path is manual column construction. Not listed in PRD §3 non-goals. |
| `offset` argument | Absent — `unused argument (offset = ...)` | out-of-scope-by-design | major | PRD §3 explicitly defers "in-fit offsets" to not-v0. But this is a real lme4 feature a user will hit; classified out-of-scope per spec, severity major because the gap is user-visible and common (exposure models). |
| `offset()` term inside formula | Refused with a clear structured diagnostic | out-of-scope-by-design | major | `lmm(y~x+offset(w)+(1\|g))` -> `mm_formula_error`: offset is "not in the engine's stateless transform subset ... must be precomputed ... or handled by the host wrapper." Consistent with PRD §3 offset deferral; refusal is explicit (good), but capability is missing. |
| `subset` argument | Absent — `unused argument (subset = ...)` | in-scope-missing | minor | `lmm(..., subset=(x>0))` errors. Trivially worked around with `data[subset, ]`, but lme4 supports it and `update()`/`drop1()` rely on it. Not in PRD §3. Minor because the workaround is one line and lossless. |
| `na.action` argument + NA handling policy | Partial — no `na.action` arg (`unused argument`); NAs always error and demand the user pre-filter. lme4 default silently `na.omit`s; `na.exclude`/`na.pass` unsupported. | partial | major | NA in response -> clear message: "mixeff requires complete cases; pass na.omit(data) explicitly before fitting." Message is clear (meets "clearer errors") but it is a behavioral divergence: lme4 fits by dropping NAs by default. `na.exclude` (pads fitted/residuals to original length) has no analog. Not listed in PRD §3, so the *missing na.action surface* is in-scope; complete-cases-only is a deliberate audit-first stance (no silent surgery) — hence partial, not pure out-of-scope. |
| `start` (theta starting values) | Absent — `unused argument (start = ...)` | in-scope-missing | minor | `lmm(..., start=c(1))` errors. Optimizer is auto-dispatched by the crate (crate-0.md:95 "Callers do not choose"); no start hook exposed. Users supply `start` to reproduce/seed difficult fits. Not in PRD §3. Minor: the crate's auto-dispatch + restart machinery reduces the need, but reproduction-from-prior-fit is impossible. |
| `verbose` (optimizer/PIRLS trace, levels 0/1/2) | Partial — `mm_control(verbose=)` controls only the pre-fit `explain_model()` printout; `-1` suppresses it, `>=0` prints it once. No optimizer iteration trace at any level. | partial | minor | `mm_control(verbose=2)` prints the random-effects explanation, not an optimizer trace. lme4 `verbose=1/2` prints per-iteration deviance / PIRLS steps for convergence diagnosis. mixeff surfaces convergence post hoc via the optimizer certificate instead. |
| `control = lmerControl(...)` object | Partial — `mm_control()` exists but exposes only `verbose`; no optimizer, tolerance, or check.* surface. The Rust side ignores the control JSON. | partial | major | `names(formals(mm_control))` == `"verbose"`. None of lmerControl's ~25 knobs (optimizer, optCtrl tolerances, restart_edge, boundary.tol, calc.derivs, check.* pre/post-fit) have an analog. Survey crate-1.md:235 claims fit-intent modes are "accessible via mm_control() options" — **not true of the shipped `mm_control()`**, which is verbose-only. |
| `optimizer` choice (nloptwrap/bobyqa/Nelder_Mead/nlminbwrap/optimx) | Absent by design — crate auto-dispatches; caller cannot choose | out-of-scope-by-design | minor | crate-0.md:92-95: "Callers do not choose" the optimizer; scalar native vs TrustBQ/NLopt by structure. Reasonable for an opinionated engine; lme4's optimizer menu has no analog. The chosen optimizer + convergence evidence are surfaced via `optimizer_certificate()` (exported). |
| `optCtrl` tolerances (ftol_abs, xtol_abs, maxeval, ...) | Absent — no tolerance/iteration controls exposed | in-scope-missing | major | The `?convergence` workflow's first recommended step is tightening `ftol_abs`/`xtol_abs` to 1e-8 and raising `maxeval`. mixeff offers no knob; a user who hits a borderline convergence warning cannot tighten tolerances or extend evaluations. Not in PRD §3. |
| `restart_edge`, `boundary.tol`, `calc.derivs`, `use.last.params` | Absent as user controls | partial | minor | No R surface. Crate runs its own bounded restart / `verify_convergence()` internally (crate-1.md:134) but it is not configurable from R. Singular/boundary detection is surfaced (`is_singular()` exported) and boundary handling exists, so behavior is covered even though the *controls* are not. |
| Pre-fit `check.*` (nobs.vs.nlev, nlev.gtr.1, rankX, scaleX, ...) configurability | Partial — checks happen and emit structured diagnostics, but action/tolerance are not user-configurable (no `check.* = "ignore"/"warning"/"stop"`) | partial | minor | mixeff's audit-first design performs design checks (rank, support, level counts) and reports them, but you cannot downgrade a `"stop"` to `"warning"` the way `lmerControl(check.nobs.vs.nlev="warning")` allows. The diagnostics are clearer than lme4's; the configurability is absent. |
| Post-fit convergence checks (check.conv.grad/singular/hess) configurability | Partial — singular fits and convergence are reported via certificate; tolerances not user-settable | partial | minor | `is_singular()` exported; `optimizer_certificate()` surfaces stop evidence + gradient/Hessian context. No `.makeCC()`-style tolerance/action configuration. |
| `devFunOnly = TRUE` (return deviance closure) | Absent | out-of-scope-by-design | minor | No analog; mixeff never returns a live R deviance closure (the native handle is a rebuildable cache, source-of-truth is JSON — CLAUDE.md / PRD §3 "no live-handle assumptions"). Manual optimization / numDeriv gradient checks are not supported. Niche; aligns with the serialize-first architecture. |
| `lFormula` / `mkLmerDevfun` / `optimizeLmer` / `mkMerMod` / `mkReTrms` modular pipeline | Absent | out-of-scope-by-design | minor | mixeff exposes `compile_model()`/`explain_model()` for the parse+validate stage but not a four-stage devfun/optimize/package pipeline for custom optimizers. Architecturally incompatible with the Rust engine boundary. Niche power-user surface. |
| `allFit()` multi-optimizer robustness | Absent (no R verb) | out-of-scope-by-design | minor | crate-1.md:134 `verify_convergence()` does bounded restart / jittered-start / alternate-optimizer checks internally and stores results in `fit$artifact`; not exposed as an `allFit`-style R function (crate-2.md:262 marks it "Not exposed to R ... partial"). The robustness *evidence* exists; the interactive verb does not. |
| `update.merMod()` | Absent | in-scope-missing | minor | No `update()` method for `mm_lmm`. Common lme4 idiom for REML<->ML swaps and tweaking formulas. Workaround: re-call `lmm()`. Not in PRD §3. |
| `refit(object, newresp)` | Works (LMM) | works | — | `refit()` exported (`R/simulate.R:13`), re-fits same formula/REML/weights with a new response; validates length + no-NA. `newweights` arg not supported (lme4 allows it). |
| `verbose` -> optimizer iteration trace for live convergence debugging | Absent | partial | minor | See `verbose` row; covered by post-fit certificate instead of live trace. |
| lmerTest auto-Satterthwaite at fit time | N/A here (covered in inference family); fit returns `mm_lmm` with df/inference computed on demand | works | — | Out of this family's scope; noted for completeness. |

---

## Summary judgement

REML/ML and weighted fitting are the load-bearing capabilities of this family,
and both are **fully working with bit-level numerical parity to lme4** (well
inside PRD §3 tolerances). Default factor contrasts work. `refit()` works.

The shipped control surface is thin. `mm_control()` exposes only `verbose`, and
that `verbose` controls the pre-fit explanation print, **not** an optimizer
trace. The Rust FFI accepts a `control_json` but ignores it (binds to
`_control`). Several lmer arguments a routine user reaches for are simply
`unused argument` errors today rather than clear, guided refusals:
`contrasts`, `subset`, `na.action`, `start`, `offset`.

Classification split:
- **out-of-scope-by-design** (PRD §3 / architecture): `offset` arg + in-formula
  `offset()` (explicit §3 deferral), `optimizer` choice, `devFunOnly`, modular
  `lFormula`/`mkLmerDevfun`/... pipeline, `allFit` verb.
- **in-scope-missing** (not in §3, real users hit them): `contrasts` arg
  (major), `optCtrl` tolerances/`maxeval` (major), `subset` (minor), `start`
  (minor), `update()` (minor).
- **partial**: `control`/`lmerControl` surface (major — only verbose),
  `na.action`/NA policy (major — complete-cases-only divergence from lme4's
  default na.omit, no na.exclude), `verbose` trace (minor), pre/post-fit
  `check.*` configurability (minor), restart/boundary controls (minor).

The two biggest practical gaps for a real lme4 migrant are: (1) **no
`contrasts=` argument** — forcing manual contrast-column construction, since
even `factor()` in-formula is refused; and (2) **no convergence-tolerance
controls** (`optCtrl`/`maxeval`), which is the very first lever lme4's
`?convergence` workflow tells users to pull. Both are in-scope-missing/major.

Note: `na.action` is a behavioral divergence — mixeff refuses NA-containing
data with a clear message ("clearer errors" is met) rather than silently
dropping rows, which is the audit-first stance ("no silent surgery"). But the
absence of any `na.action` surface (especially `na.exclude` for padded
fitted/residuals) is a genuine partial gap, not pure design choice.

Also flagged as a **documentation/contract drift**: survey crate-1.md:235 states
fit-intent modes (`design_compiled`/`as_specified`/`exploratory`/`predictive`)
are "accessible via `mm_control()` options." The installed `mm_control()` has
**only** `verbose` — those modes are not reachable from R. Either the docstring
is aspirational or the wiring is missing.
