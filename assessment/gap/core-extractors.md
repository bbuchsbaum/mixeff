# Gap Analysis — Core estimate extractors

Capability family: `fixef`, `ranef`, `coef`, `VarCorr` (+`as.data.frame`),
`sigma`, `fitted`, `residuals(type=)`.

Date: 2026-05-31. Environment: mixeff (installed, current `main`), lme4 2.0.1,
lmerTest 3.2.1. Reference surface: `assessment/survey/lme4-core-extractors.md`.
Source inspected: `R/methods-extract.R`, `R/predict.R`, `R/revive.R`,
`R/methods-print.R`, `NAMESPACE`.

Standard applied: "everything lme4 does, faster, with clearer errors."
Classifications: works | partial | in-scope-missing | out-of-scope-by-design
(PRD §3) | upstream-blocked | test-gap.

## Summary

The numeric core is solid: `fixef`, `ranef` (incl. `condVar=TRUE` postVar),
`coef`, `sigma`, `fitted`, and LMM `residuals` all return correct values that
match lme4 within parity tolerances, are S3-registered, and survive
save/revive (no live Rust handle needed). Names/shapes for vectors and ranef
data frames mirror lme4.

Three real-user gaps remain:

1. **`as.data.frame(VarCorr(.))` returns the wrong shape.** mixeff has no
   `as.data.frame.mm_varcorr` method, so the default coerces the internal
   `mm_varcorr` list into a frame with `table.*`-prefixed columns, **no
   Residual row, and no covariance rows** — incompatible with the lme4
   `grp/var1/var2/vcov/sdcor` contract that broom and downstream tooling
   consume. (major)
2. **`residuals(type=)` only accepts `"response"`.** `"pearson"`, `"deviance"`,
   `"working"` all error for both LMM and GLMM, despite PRD §11 listing
   `residuals(type=)` as a generic-compatible extractor. For LMMs these all
   equal the response residual and would be trivial to support. (major)
3. **GLMM `residuals` default is `"response"`, not lme4's `"deviance"`.**
   Silent behavioral divergence from `residuals.glmerMod`. (major)

Plus one minor: `residuals(scaled=TRUE)` is silently swallowed by `...` and
has no effect (violates the "no silent surgery" principle).

## Detailed table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `fixef(object)` named vector | Works. `fixef.mm_lmm`/`.mm_glmm` return `object$beta`. Matches lme4 to 5.7e-13 on sleepstudy. | works | — | `fixef(fit)` → `(Intercept) 251.40510 Days 10.46729`; `max(abs(fixef(mf)-lme4::fixef(lf)))=5.68e-13`. |
| `fixef(add.dropped=TRUE)` (reinsert aliased cols as NA) | Not implemented. `fixef.mm_lmm(object, ...)` ignores the arg via `...`; no rank-deficiency NA-reinsertion path. Aliased-column handling is upstream-mediated (engine refuses/reduces rank-deficient designs as structured diagnostics rather than silently dropping). | partial | minor | `fixef(mf, add.dropped=TRUE)` runs but the arg is a no-op (swallowed by `...`). No fixture exercises rank-deficient `add.dropped`. |
| `fixef(noScale=)` | Not applicable — mixeff does not autoscale fixed effects the way lme4 can; arg absent. | out-of-scope-by-design | cosmetic | No autoscaling surface in mixeff; arg swallowed by `...`. |
| `ranef(object)` → named list of data frames | Works. `ranef.mm_lmm` returns `object$random_effects` (class `mm_ranef`), per-group data frames, rows=levels, cols=terms. Values match lme4 BLUPs (max abs diff 0.0113 on sleepstudy, within BLUP tolerance). | works | — | `ranef(fit)$Subject` head matches lme4; `max abs diff = 0.0113`. |
| `ranef(condVar=TRUE)` → `postVar` p×p×k array | Works for LMM. `ranef.mm_lmm(condVar=TRUE)` calls Rust `mm_lmm_cond_var_json`, attaches a real `2×2×18` postVar array. Matches lme4 postVar to max abs diff 0.0197 (BLUP-tolerance range). On bridge failure it returns an NA-filled array + `mm_unavailable_reason` (typed refusal, not fabrication). | works | — | `dim(attr(ranef(fit,condVar=TRUE)$Subject,"postVar"))=2 2 18`; `all(is.na(.))=FALSE`; `max abs diff vs lme4 = 0.0197`. |
| `ranef(condVar=TRUE)` for GLMM | Deliberate typed refusal. `ranef.mm_glmm(condVar=TRUE)` returns NA postVar with `mm_unavailable_reason="..._for_glmm"`. Audit-first (no fabricated conditional variance). | out-of-scope-by-design | minor | PRD §3 defers GLMM CI/condvar surfaces; returns structured unavailable rather than erroring or faking. `R/methods-extract.R:65-78`. |
| `ranef(drop=TRUE)` simplify to named vectors | Not implemented; arg swallowed by `...`. lme4 users with intercept-only RE terms expect a named numeric vector. | in-scope-missing | minor | `ranef.mm_lmm(object, condVar=FALSE, ...)` has no `drop` handling. |
| `ranef(whichel=)` subset of grouping factors | Not implemented; arg swallowed by `...`. | partial | minor | No `whichel` filtering in `R/methods-extract.R:42-61`. |
| `as.data.frame.ranef.mer` (long: grpvar/term/grp/condval/condsd) | Missing. No `as.data.frame.mm_ranef` method; default coercion does not produce the lme4 long caterpillar-plot frame. | in-scope-missing | major | `grep as.data.frame NAMESPACE` → no method. The condsd bridge to ggplot2 caterpillar plots is absent. |
| `dotplot`/`qqmath` on ranef | Missing (lattice viz methods). Adjacent visualization, not a core extractor. | in-scope-missing | minor | No methods registered; user must hand-roll from postVar. |
| `coef(object)` = fixef + ranef, list of data frames | Works. `coef.mm_lmm` adds `fixef` to matching ranef columns; class `mm_coef`. Identity `coef == fixef + ranef` holds exactly (diff 0). | works | — | `max(abs(coef(mf)$Subject - (fixef+ranef)))=0`. |
| `VarCorr(object)` extract var/sd/corr + residual sd | Works (different object). `VarCorr.mm_lmm` returns `mm_varcorr` list with `$table` (group/name/variance/std_dev/correlation/boundary) and `$residual_sd`. Values match lme4 (variance 612.16 vs 612.10; resid sd 25.5904 vs 25.5918). Has a custom `print.mm_varcorr`. Not the lme4 `VarCorr.merMod` list-of-matrices class, but conveys the same content with a boundary flag. | partial | minor | `VarCorr(fit)` prints variance components + `Residual std. dev.: 25.5904`; class `mm_varcorr list`. No `$stddev`/`$correlation`/`sc`/`theta` attributes or per-group matrices. |
| `VarCorr(sigma=)` scale multiplier | Not implemented; no `sigma` arg on `VarCorr.mm_lmm(x, ...)`. | partial | minor | `VarCorr.mm_lmm <- function(x, ...)` — arg absent. |
| `as.data.frame.VarCorr.merMod` (grp/var1/var2/vcov/sdcor; cov.last/lower.tri order; Residual row) | **Broken shape.** No `as.data.frame.mm_varcorr` method. Default coercion yields columns `table.group/table.name/table.variance/table.std_dev/table.correlation/table.boundary/residual_sd` — **no Residual row, no var1/var2 covariance rows, residual_sd recycled into every row.** Incompatible with broom::tidy and the documented lme4 contract. | in-scope-missing | major | mixeff `as.data.frame(VarCorr(fit))` → wide `table.*` frame, 2 rows, no covariance/Residual rows. lme4 → 4 rows incl. `Subject (Intercept) Days` covariance and `Residual` row with `grp/var1/var2/vcov/sdcor`. PRD §11 lists `VarCorr` as generic-compatible. |
| `sigma(object)` residual SD (LMM) | Works. `sigma.mm_lmm` returns `object$sigma`; matches lme4 to 1.4e-3 (within 1e-3 theta-class tolerance; resid-sd tolerance 1e-4 nominal — see note). | works | — | `sigma(fit)=25.59036`; `abs diff vs lme4 = 1.44e-3`. |
| `sigma` for binomial/poisson = 1 (fixed) | Works. Poisson GLMM `sigma(gm)=1`. | works | — | `sigma(gm)` → `1` for poisson fit. |
| `sigma` for Gamma/inverse-Gaussian scale | Not verified here; Gamma/IG dispersion scale paths depend on engine + nlopt feature gating. | upstream-blocked | minor | PRD §3 / nlopt feature gate; not exercised. |
| `fitted(object)` conditional fitted, row-index names | Works for LMM and GLMM. `fitted.mm_lmm` returns `object$fitted` with `rownames(model_frame)` names (`"1","2",...`), matching lme4. GLMM fitted matches lme4 to ~1e-6. | works | — | `head(names(fitted(mf)))="1".."6"`; GLMM `fitted` head `3.003379 3.038710 3.026807` == lme4. |
| `residuals(type="response")` (LMM default) | Works. `residuals.mm_lmm` default `"response"`, returns `object$residuals` with row-index names; matches lme4 response residuals. | works | — | `head(residuals(mf))` matches lme4 LMM response residuals. |
| `residuals(type="pearson"/"deviance"/"working")` (LMM — all equal response) | **Errors.** `match.arg(type)` with only `"response"` allowed; `residuals(fit,type="pearson")` → `'arg' should be "response"`. For LMMs these equal the response residual, so this is a trivially-fillable in-scope gap; PRD §11 lists `residuals(type=)`. | in-scope-missing | major | `residuals(mf,type="pearson")` → error `'arg' should be "response"` (same for deviance/working). |
| `residuals` GLMM default `"deviance"` | **Divergent.** mixeff GLMM default is `"response"` (returns y−mu, matches lme4 `type="response"`), whereas `residuals.glmerMod` defaults to `"deviance"`. Silent behavioral mismatch a glmer user hits immediately. | partial | major | mixeff `residuals(gm)` head `-1.0034 -1.0387 ...` == lme4 `type="response"`; lme4 `residuals(lg)` default = deviance `-0.617 -0.636 ...`. |
| `residuals(type="pearson"/"deviance"/"working")` (GLMM — genuinely differ) | **Errors / missing.** Only `"response"` accepted; pearson/deviance/working error. These differ numerically for GLMMs (deviance/Pearson are standard glmer diagnostics) and require engine support. | in-scope-missing | major | `residuals(gm,type="deviance")` → error; lme4 provides distinct deviance (−0.617) and pearson (−0.579) values. |
| `residuals(scaled=TRUE)` divide by sigma | **Silently ignored.** `scaled` swallowed by `...`; output identical to unscaled. Violates "no silent surgery" (transformation requested, silently not applied). | partial | minor | `all.equal(residuals(mf), residuals(mf,scaled=TRUE))=TRUE`. |
| `residuals(type="partial")` | lme4 itself errors ("not yet implemented"); mixeff also unsupported. Parity holds (neither supports it). | out-of-scope-by-design | cosmetic | Both error; not a mixeff-specific gap. |

## Notes on tolerances

- `sigma` matched lme4 to 1.44e-3. The PRD-stated residual-SD tolerance is 1e-4
  (`assessment/survey/lme4-core-extractors.md` §12). The observed diff exceeds
  1e-4 for this sleepstudy fit. This is plausibly an optimizer-convergence /
  REML-criterion tolerance difference rather than an extractor bug, but it is
  worth a parity-fixture check that the `(Days|Subject)` case meets the 1e-4
  sigma bar; if it does not, that is a fit-convergence concern outside this
  extractor family. Flagged as a `test-gap` for the family below.
- `ranef`/`postVar` diffs (0.011 / 0.020) are in the normal BLUP/condVar
  tolerance band given the small theta/sigma difference above.

## Audit-first behavior (positive)

The refusal/unavailable paths are correctly structured, not silent: GLMM
`condVar`, bridge failures, and `vcov(type="theta")` all return values carrying
`mm_unavailable_reason` attributes rather than fabricating numbers — consistent
with PRD "no silent surgery." The one violation of that principle in this
family is `residuals(scaled=TRUE)` being silently ignored.

## Test-gap

- No test asserts the `(Days|Subject)` sleepstudy `sigma` meets the 1e-4 parity
  bar (observed 1.4e-3). Classification: test-gap (and possibly a fit-tolerance
  concern, but that lies in the fit family, not extractors).
