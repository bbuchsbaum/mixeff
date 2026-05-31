# Gap Report — Model-statistics generics

**Family:** Model-statistics generics — `logLik` (REML/ML), `deviance`, `AIC`,
`BIC`, `nobs`, `df.residual`, `extractAIC`, `REMLcrit`, `devcomp`, component
extraction via `getME`.

**Reference:** lme4 2.0.1 / lmerTest 3.2.1, `assessment/survey/lme4-model-stats.md`.
**mixeff source:** `R/methods-extract.R` (logLik/deviance/AIC/BIC/nobs/df.residual),
`R/revive.R` (getME), `R/fit-lmm.R` / `R/glmm.R` (stored fit fields).
**Date:** 2026-05-31. All statuses below are confirmed by running
`library(mixeff); library(lme4)` on `sleepstudy` (`Reaction ~ Days + (Days|Subject)`),
both REML and ML fits, unless noted.

The bar for this family: "everything lme4 does, faster, with clearer errors."
Numeric agreement on the optimized criterion is excellent. The gaps are in
(a) the REML/ML *override* semantics of `logLik`/`deviance`, (b) entirely
missing generics (`extractAIC`, `isREML`, `REMLcrit`, `refitML`, `devcomp`),
and (c) `getME` coverage that is far narrower than both lme4 and the mixeff
PRD's own promise (PRD §6, lines 279-280).

---

## Summary table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `logLik(m)` value (REML fit) | works | works | — | mixeff `-871.8141` = lme4 `-871.8141`; class `logLik`, `df=6`, `nobs=180` |
| `logLik(m)` value (ML fit) | works | works | — | mixeff = lme4 `-875.99` |
| `logLik` `df` attribute | works | works | — | `attr(ll,"df")==6` matches lme4 |
| `logLik` `nobs` attribute | works | works | — | `attr(ll,"nobs")==180` |
| `logLik` `nall` attribute | **missing** | partial | minor | `attr(logLik(m),"nall")` is `NULL` in mixeff; lme4 returns `180`. `bbmle`/`AICcmodavg` inspect `nall` |
| `logLik(m, REML=FALSE)` override on a REML fit | **wrong value (ignored)** | partial | **major** | mixeff returns `-871.8141` (the REML logLik); lme4 returns `-875.9929` (ML logLik re-evaluated at REML θ). `REML` is swallowed by `...` — `logLik.mm_lmm <- function(object, ...)` has no `REML` arg (`R/methods-extract.R:291`) |
| `logLik(m, REML=TRUE)` override on an ML fit | **wrong value (ignored)** | partial | **major** | mixeff returns the ML logLik unchanged; lme4 re-evaluates the REML criterion at ML θ |
| `deviance(m)` (ML deviance, lme4 contract) | **returns REML criterion for REML fits** | partial | **major** | `m$deviance` for a REML fit stores `1743.628` = `REMLcrit`, not the ML deviance `1751.99`. lme4 `deviance(m, REML=FALSE)` = `1751.99`. mixeff matches lme4's *deprecated* default-arg behavior (which also returns REMLcrit + a deprecation warning) but diverges from the documented "always ML deviance" contract |
| `deviance(m)` (ML fit) | works | works | — | mixeff `1751.939` = lme4 `1751.939` |
| `deviance(m, REML=FALSE)` argument | **silently ignored** | partial | **major** | mixeff returns `1743.628` regardless; lme4 returns `1751.986`. `deviance.mm_lmm <- function(object, ...)` (`R/methods-extract.R:306`) has no `REML` arg and ignores it without a diagnostic — violates the "no silent surgery" contract |
| `AIC(m)` single model (REML fit) | works | works | — | mixeff `1755.628` = lme4 `1755.628`; uses fit criterion (REML) like lme4 |
| `AIC(m)` single model (ML fit) | works | works | — | mixeff `1763.939` = lme4 `1763.939` |
| `AIC(..., k=...)` | works | works | — | `k` honored: `AIC.mm_lmm(object,...,k=2)` computes `-2*logLik + k*dof` |
| `AIC(m1, m2, ...)` multi-model | **refuses** | partial | minor | mixeff aborts `mm_inference_unavailable` directing user to `compare()`. lme4 returns a `df`/`AIC` data frame (with REML-mismatch warning). Clear error, deliberate routing — but a real lme4 user who types `AIC(m1,m2)` hits a wall |
| `BIC(m)` single model | works | works | — | mixeff `1774.786` = lme4 `1774.786` |
| `BIC(m1, m2, ...)` multi-model | **refuses** | partial | minor | same routing-to-`compare()` abort as `AIC` |
| `extractAIC(fit)` | **missing** | in-scope-missing | **major** | `extractAIC(m)` → "no applicable method for 'extractAIC'". No `extractAIC.mm_lmm` in `R/` or NAMESPACE. Breaks `step()`/`drop1()` and any base-R selection caller |
| `nobs(m)` | works | works | — | mixeff `180` = lme4 `180` |
| `df.residual(m)` | works | works | — | mixeff `174` = lme4 `174` (N − p for this model; survey doc's "178" is for a 2-fixef-only example). Note mixeff returns the stored `df_residual`, not recomputed — value agrees here |
| `REMLcrit(m)` | **missing** | in-scope-missing | **major** | lme4's `REMLcrit(m)` errors on mixeff (`no applicable method for 'isLMM'`). No generic/method defined. The value *exists* (`m$deviance` for REML fits IS the REML criterion) but there is no `REMLcrit()` accessor. This is the only lme4 function exposing the REML criterion as a scalar |
| `isREML(m)` | **missing** | in-scope-missing | major | `isREML(m)` → "no applicable method". State is available as `m$REML` (logical) but the lme4 generic is unimplemented; programmatic guards (`if (isREML(m)) ...`) break |
| `refitML(m)` | **missing** | in-scope-missing | major | `refitML(m)` → "no applicable method". No method. Users converting REML→ML before LRT cannot. (`refit(m, newresp)` IS implemented — `R/simulate.R:19` — but that is response-refit, not criterion-refit) |
| `devcomp(m)` (lme4 deprecated → `getME(.,"devcomp")`) | **missing** | in-scope-missing | major | `devcomp(m)` errors (`is(x,"merMod")`); `getME(m,"devcomp")` → "component `devcomp` is not available". **PRD §6 (lines 279-280) explicitly lists `devcomp` as in-scope for `getME`.** Not delivered |
| `getME(m, "X")` | works | works | — | returns `180×2` dense matrix |
| `getME(m, "Z")` | works | works | — | returns `180×36` sparse matrix |
| `getME(m, "Zt")` | **broken** | partial | major | `getME(m,"Zt")` → error "argument is not a matrix". The `Zt = t(model.matrix(object,"random"))` branch (`R/revive.R:128`) calls `t()` on a representation that is not a base matrix, so the documented `Zt` name throws instead of returning the q×N transpose |
| `getME(m, "theta")` | works | works | — | length-3 vector |
| `getME(m, "beta")` / `"fixef"` | works | works | — | length-2 |
| `getME(m, "Lambda")` / `"Lambdat"` | works | works | — | `36×36` |
| `getME(m, "y")` / `"mu"` | works | works | — | length-180 |
| `getME(m, "flist")` / `"cnms"` | works | works | — | length-1 lists |
| `getME(m, "sigma")` | **missing** | partial | minor | not a `getME` name; `sigma(m)` works as its own generic, but `getME(m,"sigma")` errors. lme4 supports it |
| `getME(m, "N"/"n"/"p"/"q")` dimension scalars | **missing** | in-scope-missing | major | all error in mixeff. Common lme4 idiom (`getME(m,"p")`, `getME(m,"N")`). Trivially derivable from stored fields |
| `getME(m, "Gp")` | **missing** | in-scope-missing | major | errors. **PRD §6 lists `Gp` as in-scope** |
| `getME(m, "lower")` | **missing** | in-scope-missing | major | errors. **PRD §6 lists `lower` as in-scope** |
| `getME(m, "optinfo")` | **missing** | in-scope-missing | major | errors. **PRD §6 lists `optinfo` as in-scope** |
| `getME(m, "REML"/"is_REML")` | **missing** | partial | minor | both error; state exists as `m$REML`. lme4 returns `2L`/`TRUE` |
| `getME(m, "Ztlist"/"mmList"/"u"/"b"/"L"/"RX"/"RZX"/"A"/"Lind"/"Tlist"/...)` | **missing** | partial | minor | full lme4 internal-Cholesky/optimizer surface absent. PRD §6 scopes `getME` to a named subset ("X, Z, theta, Lambda, cnms, flist, Gp, lower, devcomp, optinfo"), so deep solver internals beyond that list are by-design narrower |
| `llikAIC(object)` | **missing** | out-of-scope-by-design | cosmetic | lme4-internal helper; mixeff builds its AIC table in `R/methods-summary.R`. Not a user-facing surface; no PRD commitment |
| `devfun2()` / `varianceProf()` | **missing** | out-of-scope-by-design | minor | profile-likelihood machinery. PRD §3 (line 40) defers profile-likelihood CIs; v0 explicitly will not ship them. Cite PRD §3 |
| `refit(m, newresp)` | works | works | — | implemented `R/simulate.R:19` (separate family, listed in survey §4); validates `length(newresp)==nobs` |

---

## Detailed findings (load-bearing)

### F1 — `logLik(m, REML=...)` override is silently ignored (major, partial)

`R/methods-extract.R:291`:
```r
logLik.mm_lmm <- function(object, ...) {
  structure(object$logLik, df = object$dof, nobs = object$nobs, class = "logLik")
}
```
There is no `REML` parameter. lme4's contract (survey §1, §9.1) is that
`logLik(m_reml, REML=FALSE)` *re-evaluates* the ML log-likelihood at the stored
REML θ (returns `-875.99` vs the REML `-871.81`). mixeff returns the stored
optimized criterion unconditionally, so `logLik(m, REML=FALSE)` and
`logLik(m, REML=TRUE)` are both wrong whenever the requested criterion differs
from the fit criterion. This is exactly the kind of cross-criterion query that
LRT/AIC-under-ML workflows make. Worse than a clear refusal: it returns a
plausible-but-wrong number with no diagnostic — a "silent surgery" violation
of CLAUDE.md's stated contract.

### F2 — `deviance()` returns the REML criterion for REML fits, and ignores `REML=` (major, partial)

`m$deviance` for a REML fit stores `1743.628`, which equals `REMLcrit`, not the
ML deviance `1751.99`. lme4's documented contract (survey §1) is that
`deviance()` always returns the ML deviance; `deviance(m_reml, REML=FALSE)`
gives `1751.986`. mixeff's `deviance.mm_lmm` has no `REML` argument and returns
the stored value regardless, so it (a) disagrees with lme4's documented ML-
deviance contract for REML fits and (b) silently ignores `REML=FALSE`. The
underlying Rust fit appears to store the REML criterion in the `deviance` slot
for REML fits; the wrapper does not separate "REML criterion" from "ML
deviance." Note lme4 2.0.1's *default*-arg `deviance(m_reml)` also returns
REMLcrit (with a deprecation warning), so the no-arg call coincidentally
matches — the divergence is real only under `REML=FALSE` and against the
documented contract.

### F3 — Missing generics: `extractAIC`, `REMLcrit`, `isREML`, `refitML`, `devcomp` (major, in-scope-missing)

None are defined (confirmed by NAMESPACE scan + runtime "no applicable method").
`extractAIC` breaks base-R `step()`/`drop1()`. `REMLcrit` is the only scalar
accessor for the REML criterion (the value exists internally but is
unreachable via the lme4 name). `isREML` has the state available (`m$REML`) but
no generic. `refitML` (criterion conversion) is absent though response-`refit`
exists. `devcomp` is absent AND is explicitly promised by PRD §6.

### F4 — `getME` is much narrower than lme4 and than the PRD promise; `Zt` is broken (mixed)

Implemented names: `X, Z, Zt(broken), Lambda, Lambdat, theta, beta, fixef, y,
mu, flist, cnms` (`R/revive.R:114-147`). 

- `Zt` throws "argument is not a matrix" — a declared-but-broken branch (major).
- PRD §6 (lines 279-280) commits to `Gp`, `lower`, `devcomp`, `optinfo` —
  all absent (in-scope-missing, major each).
- Common dimension scalars `N/n/p/q/sigma` and state `REML/is_REML` are absent
  though trivially derivable from stored fields.
- Deep solver internals (`L`, `RX`, `RZX`, `A`, `Ztlist`, `u`, `b`, `Lind`,
  `Tlist`, `devfun`, ...) are outside the PRD's named subset → narrower by
  design, but undocumented as a deliberate boundary in the help page.

`getME.mm_lmm`'s default branch DOES emit a clear, structured error
(`mm_abort` / `mm_arg_error`), which satisfies the "clearer errors" goal for
unsupported names — but the PRD-promised names and the broken `Zt` are defects,
not honest refusals.

---

## Test-gaps observed

No tests exercise: the `REML=` override on `logLik`/`deviance` (the F1/F2
bugs), the `nall` attribute, `getME("Zt")`, `getME("devcomp")` (PRD-promised),
`extractAIC`, `isREML`, `REMLcrit`, `refitML`. `tests/testthat/test-lmm.R:58-68`
only checks that extractors echo stored fields (`deviance(fit)==fit$deviance`),
which cannot catch F2 because it tests against the same wrong stored value.
The parity tutorial tests only assert `as.numeric(logLik(...))` equality (the
no-override path, which works). Classification: **test-gap** for all of the
above — the bugs would have been caught by an lme4 cross-check on the override
arguments.
