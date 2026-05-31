# Gap Report — `ranef` / `condVar` & Plotting

Family: ranef condVar & plotting (`ranef(condVar=TRUE)` postVar attr,
`as.data.frame.ranef.mer` with `condsd`, `dotplot`/caterpillar, `qqmath`).

Date: 2026-05-31
mixeff: installed snapshot on branch `main`
lme4 2.0.1, lattice loaded.

Bottom line: the **numeric core** of the family — `ranef(condVar=TRUE)` for
LMMs producing a real `p×p×k` `postVar` array that matches lme4 within
tolerance — **works**. Everything built *on top of* that array for a real
lme4 user (the long-format `as.data.frame` with `condsd`, the
`with conditional variances` print line, and the `dotplot` / `qqmath` /
`plot` caterpillar & diagnostic plots) is **missing**, and in the
`as.data.frame` case actively produces wrong output by silently falling
through to a base method. For GLMMs the `postVar` itself is a typed-refusal
stub. The `drop` and `whichel` arguments are silently ignored.

---

## Capability matrix

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `ranef(condVar=TRUE)` LMM: attach `p×p×k` `postVar` array | works | works | — | `ranef.mm_lmm` routes through Rust `wrap__mm_lmm_cond_var_json` (`R/methods-extract.R:42-61, 144-212`). Live: `dim 2x2x18`, all finite, symmetric, PSD. |
| `postVar` numeric parity vs lme4 (PRD theta tol 1e-3 / sigma 1e-4) | works | works | — | condsd 12.07167 vs lme4 12.07086 (rel 2.4e-4); Days slope condsd 2.305062 vs 2.304839; intercept-only model max abs diff 3.95e-6. `test-ranef-condvar.R` (5 assertions) all pass. |
| `postVar` caching (no refit on repeated calls) | works | works | — | cached on `fit$lazy_cache` via `.mm_lazy(fit,"cond_var",...)` (`methods-extract.R:144-146`); test "caches across repeated calls" passes. |
| `ranef(condVar=FALSE)` default (no `postVar`) | partial | partial | minor | mixeff default IS `condVar=FALSE` (lme4 default is `TRUE`). No `postVar` attached, no `mm_unavailable_reason` — correct. Default differs from lme4 but is the safe direction; documented in `@param condVar`. |
| `ranef(condVar=TRUE)` GLMM | partial (typed refusal) | partial | major | `ranef.mm_glmm` returns `postVar` array of all-`NA` with `attr "mm_unavailable_reason" = "..._for_glmm"` (`methods-extract.R:65-78`). Honest stub, not a crash, but a real lme4 GLMM user gets no conditional variances. |
| `as.data.frame(ranef(...))` → long DF (`grpvar,term,grp,condval,condsd`) | **missing / wrong output** | in-scope-missing | **blocker** | No `as.data.frame.mm_ranef` method. Live call falls through to base `as.data.frame.list`, yielding a WIDE frame with mangled names `Subject..Intercept.`, `Subject.Days`, NO `grpvar/term/grp/condval`, and NO `condsd` even with `condVar=TRUE`. Silently wrong, not an error. This is the canonical ggplot2 caterpillar workflow input. |
| `print.ranef.mer` "with conditional variances for …" line | partial | partial | minor | `print.mm_ranef` (`R/methods-print.R:147-158`) prints the per-group frames but never emits the `with conditional variances for "Subject"` reminder line, even when a real `postVar` is attached. |
| `dotplot.ranef.mer` caterpillar (lattice) | missing | in-scope-missing | **major** | No `dotplot.mm_ranef`; `getS3method("dotplot","mm_ranef")` → none. Live `dotplot(re)` errors `cannot xtfrm data frames` (dispatches to default). The single most common RE visualization is absent. |
| `qqmath.ranef.mer` Q-Q plot w/ CI bars (lattice) | missing | in-scope-missing | **major** | No `qqmath.mm_ranef`. Live: `no applicable method for 'qqmath' applied to … "mm_ranef"`. |
| `plot.ranef.mer` (qqmath/xyplot/splom dispatch) | missing | in-scope-missing | major | No `plot.mm_ranef`. Live `plot(re)` errors `'x' is a list, but does not have components 'x' and 'y'`. |
| `transf` arg on `dotplot`/`qqmath` (e.g. `exp` for OR) | missing | in-scope-missing | major | Depends on the absent plot methods. |
| `level` arg (CI width for error bars) | missing | in-scope-missing | major | Depends on the absent plot methods. |
| `drop=TRUE` (simplify single-col groups to named vectors) | missing (silent) | in-scope-missing | minor | `ranef.mm_lmm`/`ranef.mm_glmm` accept only `condVar` then `...`; `drop=TRUE` is swallowed by `...` and ignored — no simplification, no warning. Violates "no silent surgery" spirit (silently drops a requested transformation). |
| `whichel` (subset of grouping factors) | missing (silent) | in-scope-missing | minor | Same: swallowed by `...`, ignored, no error. |
| `postVar` (deprecated synonym for `condVar`) | missing | partial | cosmetic | Not accepted; would be ignored via `...`. Deprecated in lme4 anyway; low value. |
| Multi-term-same-group block-diagonal `postVar` (`(1\|g)+(0+t\|g)`) | works | works | — | `mm_merge_block_diag_postvar` (`methods-extract.R:208-239`) assembles block-diagonal array; covered by the design comment and exercised by the merge path. |

---

## Detail and evidence

### What genuinely works (the hard part)

`ranef(condVar=TRUE)` for an LMM is the numerically difficult piece, and it
is implemented end-to-end: `ranef.mm_lmm` calls `mm_cond_var_postvars` →
`mm_compute_cond_var_postvars` → Rust `wrap__mm_lmm_cond_var_json`
(`R/methods-extract.R:42-212`), validates the schema header
`mixeff.lmm_cond_var v1`, builds a `p×p×k` array with correct dimnames
(slope × slope × level), aligns it to the `ranef` frame, and caches it.
Parity is well within the PRD bar:

```
mixeff condsd[1] (Intercept): 12.07167   lme4: 12.07086   (rel 6.7e-5)
mixeff condsd[1] Days:         2.305062   lme4:  2.304839
intercept-only (1|Subject):   max abs diff 3.95e-6
```

The aggregate `max(abs(diff))` across the full array is 0.0197, but that is
the absolute scale of the off-diagonal Intercept×Days covariance entry; the
relative error is 2.4e-4 and the `tolerance=1e-3` assertion in
`test-ranef-condvar.R` passes. All 5 tests in that file pass.

So the family's foundation is sound. The gaps are everything a user reaches
for *after* getting the array.

### `as.data.frame` — the worst gap (silently wrong)

There is no `as.data.frame.mm_ranef`. Because `mm_ranef` inherits from
`list`, `as.data.frame()` dispatches to the base list method and returns a
wide frame:

```
    Subject..Intercept. Subject.Days
308            2.249080    9.2010998
...
has condsd col?: FALSE
```

vs lme4's contract:

```
   grpvar        term grp    condval   condsd
1 Subject (Intercept) 308   2.258551 12.07086
```

This is the canonical input to the lme4 ggplot2 caterpillar recipe
(`lme4-ranef-condvar.md` §"Typical ggplot2 workflow"). A user copying that
recipe gets neither the `condval`/`grp` columns it references nor the
`condsd` needed for error bars — and gets **no error**, just wrong columns.
This is a "no silent surgery" violation: the postVar that mixeff went to the
trouble of computing is dropped on the floor at the `as.data.frame` step.
Severity blocker because it breaks the most common downstream workflow
without signalling failure.

### Plot methods — entirely absent

`dotplot`, `qqmath`, and `plot` have no `mm_ranef` methods (NAMESPACE
exports only `print.mm_ranef`, `ranef.mm_glmm`, `ranef.mm_lmm`). Live calls
all error (dispatch to lattice/base defaults). For a package whose standard
is "everything lme4 does, faster, with clearer errors," the caterpillar plot
(`dotplot`) and the random-effect Q-Q plot (`qqmath`) are core diagnostic
output, and they are missing. The errors are also *not* clear mixeff
diagnostics — they are raw lattice/base errors (`cannot xtfrm data frames`,
`no applicable method for 'qqmath'`), which is itself off-contract.

### `print` reminder line

`print.mm_ranef` omits the lme4 `with conditional variances for "Subject"`
line. Minor, but it is the only textual cue that uncertainty is attached, so
its absence is a small parity gap.

### `drop` / `whichel` silently ignored

Both are accepted only through `...` and discarded with no effect and no
warning. A user passing `drop=TRUE` expecting a named vector gets a
one-column data frame instead, silently.

---

## PRD scope check

PRD §3 non-goals were reviewed. None of these capabilities are listed as
deferred or out-of-scope: the §3 deferrals concern profile-LL GLMM CIs,
multivariate `cbind`, Kenward-Roger beyond scalar, AR(1)/spatial residuals,
`nlmer`, `I()/poly()/splines`, GAM smooths, and the model-selection engine.
`ranef`, conditional variances, `as.data.frame`, and the lattice plot
methods are *not* excluded. Therefore the missing pieces here are
**in-scope-missing**, not out-of-scope-by-design. The GLMM-condVar refusal
is the one defensible deferral (conditional variances for GLMM modes are
genuinely harder and adjacent to the deferred GLMM-inference surface), but it
is currently a `partial` typed refusal, not a documented non-goal.

## Test coverage note

`test-ranef-condvar.R` covers the LMM `postVar` array thoroughly (shape,
finiteness, symmetry, PSD, lme4 parity, caching, FALSE-default). There are
**no** tests for `as.data.frame.mm_ranef`, the print reminder line,
`dotplot`/`qqmath`/`plot`, or `drop`/`whichel` — consistent with those
surfaces not existing yet (test-gap follows the implementation gap).
