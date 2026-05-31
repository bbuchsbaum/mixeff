# Gap Analysis — Diagnostics & Convergence

**Family:** Diagnostics & convergence (`isSingular`, convergence warnings/checks
(`checkConv`), `rePCA`, `hatvalues`, `influence.merMod`, `cooks.distance`,
`dfbeta`/`dfbetas`, `allFit`, theta/Hessian checks, conditional variances,
residual diagnostics).

**Reference:** `assessment/survey/lme4-diagnostics.md` (lme4 2.0.1 / lmerTest 3.2.1).

**Standard:** everything lme4 does, faster, with clearer errors. Gaps that a real
lme4 diagnostics user would hit and that are in scope are at least **major**.

**Evidence basis:** live runs against installed `mixeff` + `lme4` on `sleepstudy`
(non-singular `Reaction ~ Days + (Days|Subject)` and a forced reduced-rank
`(Days + noise | Subject)` fit). Source: `R/revive.R`, `R/diagnostics.R`,
`R/methods-extract.R`, `R/predict.R`, `R/methods-print.R`, `NAMESPACE`.

---

## Summary table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `isSingular(object, tol)` | works (via own `is_singular()`) | works | — | `is_singular(fm)` returns `FALSE`; on reduced-rank fit returns `TRUE`. Own generic `is_singular.mm_lmm` (R/revive.R:554) reading `fit_status` + `effective_covariance$status`. `lme4::isSingular(fm)` errors (no method) — see Cross-API note. |
| `getSingTol()` / `options(lme4.singular.tolerance)` | partial | partial | minor | `is_singular(x, tol=)` accepts `tol` but it is "Reserved for compatibility" and **ignored** (R/revive.R:542,554); singularity is status-driven, not eigenvalue/tol-driven. No `getSingTol()` and no `lme4.singular.tolerance` option. |
| `rePCA(x)` (prcomp per grouping factor: SDs + rotation) | partial | partial | major | No `rePCA()` function/method (`getS3method` FALSE; not in NAMESPACE). The *purpose* (find near-zero variance dimensions) is served by the audit-first `effective_covariance` artifact: `requested_rank`, `supported_rank`, `directions` (with loadings), `unsupported_directions`, surfaced via `changes(fit)` and the print "Fitted covariance state" block. But there is no prcomp-style object (no orthogonal SDs / rotation matrix a user can index), so workflows that call `rePCA(fm)$Subject$sdev` have no port. |
| `checkConv(...)` (gradient/Hessian post-opt checks) | partial | partial | major | No `checkConv()` and no `lmerControl(check.conv.*)` knobs. The engine ships an `OptimizerCertificate`: `optimizer_certificate(fm)` exposes `status` (`converged_interior`/`converged_boundary`/`converged_reduced_rank`), `free_gradient_norm`, `projected_gradient_norm`, `hessian_eigen_min`, `hessian_rank`, `information_rank`. This *is* a derivative-backed convergence check, surfaced honestly. But there is no user-tunable check (no `.makeCC`, no `tol`/`action`), no per-check messages list, and no way to set custom gradient/Hessian tolerances. |
| `convergence` help-topic workflow (tighten optCtrl, restart, numDeriv recheck, allFit) | partial | partial | minor | Optimizer is engine-side (`trust_bq`); `mm_control()` exists but the documented lme4 remedies (xtol_abs/ftol_abs tuning, restart_edge, numDeriv recheck) are not user-exposed verbs. No narrative convergence vignette mapping the lme4 recovery steps. |
| `lmerControl(check.conv.grad/singular/hess, restart_edge, boundary.tol, calc.derivs, check.nobs.vs.*, check.rankX, check.scaleX)` | partial | partial | major | None of the `check.*` control toggles exist. Identifiability/rank/scale checks are instead performed engine-side as design diagnostics (`covariance_too_rich`, `effective_covariance` rank reduction) and surfaced via `diagnostics(fit)` / `audit_design()` / `explain_model()` — non-silent by design. But a user cannot set `action`/`tol` per check, cannot turn checks off, and cannot request `check.scaleX`-style fixed-effect scaling warnings. |
| `.makeCC(action, tol, relTol)` | in-scope-missing | in-scope-missing | minor | No equivalent; depends on a `check.conv.*` framework that does not exist. |
| `fm@optinfo` (`$derivs$gradient`, `$derivs$Hessian`, `$conv`, `$feval`, `$message`, `$val`) | partial | partial | major | No `optinfo` slot. `optimizer_certificate(fm)` covers `status`, gradient *norms* (not the full gradient vector), `hessian_eigen_min`/`hessian_rank` (not the full Hessian matrix), and `iterations` (~`feval`). The **full gradient vector and Hessian matrix are not exposed**, so the standard lme4 recheck `numDeriv::hessian(getME(fm,'devfun'), theta)` vs `fm@optinfo$derivs$Hessian` cannot be reproduced, and `Matrix::rcond(Hessian)` cannot be computed (only eigen-min is given). |
| `allFit(object, ...)` (multi-optimizer gold-standard) | in-scope-missing | in-scope-missing | major | No `allFit`; no `mm_lmm` method, not in NAMESPACE. Engine uses a single optimizer (`trust_bq`). The canonical lme4 false-positive convergence check (refit with all optimizers and compare fixef/llik/theta) has no port. |
| `getME(object, name)` diagnostic names: `theta`, `lower`, `devfun`, `devarg`, `par`, `ST`, `L`, `Lambda`, `Lambdat` | partial | partial | major | `getME(fm, "theta"/"beta"/"fixef"/"X"/"Z"/"Zt"/"Lambda"/"Lambdat"/"y"/"mu"/"flist"/"cnms")` all OK. **`"lower"`, `"devfun"`, `"devarg"`, `"par"`, `"ST"`, `"L"` all error** ("component not available"). The boundary test `theta == getME(fm,'lower')` (lme4's documented singularity diagnostic) cannot be done; `getME(fm,'devfun')` for external gradient/Hessian recheck is unavailable. |
| `devfun2(fm, useSc, scale)` (SD/correlation-scale deviance fn) | in-scope-missing | in-scope-missing | minor | No equivalent. Used for manual profiling / interpretable Hessian; absent. |
| `profile(fitted, which, ...)` / `varianceProf` / `logProf` / `confint(method="profile")` | out-of-scope-by-design (LMM is in scope; was deferred) | partial | major | PRD §3 lists "profile-likelihood CIs" as a v0 non-goal, but the untracked `tests/testthat/test-confint-profile.R` and `R/compare.R`/`R/inference.R` changes suggest profile CIs are being added for LMM. No `profile.mm_lmm` method exists today (`getS3method` FALSE). `varianceProf`/`logProf` (profile transforms) and the `thpr` object are absent. Classify the GLMM case as out-of-scope (PRD §3 defers profile-LL for GLMM); the LMM case is at most partial/in-progress. |
| `influence(model, groups, ...)` (deletion influence, leave-one-group-out) | in-scope-missing | in-scope-missing | major | No `influence` method (`getS3method` FALSE; `influence(fm)` errors "no applicable method"). Crate survey (crate-1.md:262) explicitly lists "Influence diagnostics requiring repeated refits" as not yet provided. lme4 users doing per-subject deletion influence have no port. |
| `cooks.distance(influence.merMod)` | in-scope-missing | in-scope-missing | major | Depends on `influence.merMod`; no `cooks.distance.mm_*` method. Absent. |
| `dfbeta(influence.merMod, which)` / `dfbetas` | in-scope-missing | in-scope-missing | major | No methods; depend on `influence`. Absent. |
| `hatvalues(model, fullHatMatrix)` (leverage diagonal / full hat matrix) | in-scope-missing | in-scope-missing | major | No `hatvalues.mm_lmm` (`getS3method` FALSE; `hatvalues(fm)` errors). Leverage diagnostics for LMMs are entirely absent; X and Z are available via `getME`/`model.matrix` so a user could hand-roll, but there is no supported verb. |
| `residuals(object, type, scaled)` — `type` in {response, pearson, deviance, working} | partial | partial | major | `residuals.mm_lmm/glmm` (R/predict.R:109) hard-codes `type = c("response")`; `residuals(fm, type="pearson"/"deviance"/"working")` **all error** with `'arg' should be "response"`. `scaled=TRUE` is *silently ignored* (returns same length but the `scaled` arg is swallowed by `...`, not applied — verify). For GLMM this is a real defect: lme4's GLMM residual default is `"deviance"`, and Pearson residuals are the standard residual diagnostic — neither is available. |
| `plot(merMod)` residual diagnostic plot | in-scope-missing | in-scope-missing | minor | No `plot.mm_lmm`; no residual-vs-fitted diagnostic plot. (Plotting is thin in lme4 too; minor.) |
| `fortify.merMod` (add resid/hat/cooks columns for ggplot) | in-scope-missing | in-scope-missing | cosmetic | Deprecated in lme4 in favour of broom.mixed; no port. |
| `dotplot.ranef.mer` (caterpillar plot w/ PI from condVar) | in-scope-missing | in-scope-missing | minor | No dotplot method, but condVar (postVar) is available, so a user can build it. |
| `ranef(object, condVar=TRUE)` (postVar conditional variances) — LMM | works | works | — | `ranef(fm, condVar=TRUE)` attaches a real `postVar` array (2x2x18 for `(Days|Subject)`), routed through Rust `cond_var()` bridge (R/methods-extract.R:144). Note **default is `condVar=FALSE`** (lme4 default is `TRUE`); minor API divergence. |
| `ranef(condVar=TRUE)` — GLMM | partial | partial | minor | Returns `NA` postVar with structured reason `random_effect_conditional_variance_unavailable_for_glmm` (R/methods-extract.R:65) — honest refusal, not a crash. lme4 does provide GLMM condVar, so this is a real gap, but flagged honestly. |
| `VarCorr(x)` (singularity inspection: zero var / ±1 corr) | works | works | — | `VarCorr.mm_lmm`/`mm_glmm` exported (NAMESPACE). Provides the variance/correlation view used to eyeball singularity. |
| `isLMM`/`isGLMM`/`isNLMM`/`isREML` predicates | partial | partial | minor | Not exported as such; model type is carried by S3 class (`mm_lmm`/`mm_glmm`) and `inherits()`. `isREML` has no direct accessor, though REML status is in the artifact/print. lme4 diagnostic code that branches on `isLMM(x)` will not run unmodified. |
| `Nelder_Mead`/`nloptwrap`/`nlminbwrap` exported optimizers | out-of-scope-by-design | out-of-scope-by-design | minor | Optimization is engine-side (`trust_bq` in Rust). PRD architecture (§4) places the optimizer inside the crate; exposing R-level optimizer functions is not a goal. Not a parity defect. |
| `check.rankX` (fixed-effect design rank) | partial | partial | minor | No `check.rankX` toggle, but rank-deficiency of the design is handled non-silently engine-side via design diagnostics rather than column-dropping options. The lme4 `"message+drop.cols"` behavior (auto-drop, multiple action levels) is not user-configurable. |
| `check.scaleX` (badly-scaled fixed columns) | in-scope-missing | in-scope-missing | minor | No fixed-effect scaling warning surfaced; lme4 warns by default. A user with predictors on wildly different scales gets no advisory. |
| `check.nobs.vs.nlev` / `nRE` / `nlev.gtr.1` / `nlev.gtreq.5` (RE identifiability) | partial | partial | minor | RE support is assessed engine-side (`support: sufficient/too_rich`, `min rows/group`, `group levels`) and surfaced in `random_blocks()` / print. The specific lme4 `"stop"` on `nlevels <= 1` and `"warning"` on `< 5 levels` are not reproduced as distinct, configurable checks, but the data is present. |

---

## Notes on classification

### Audit-first design substitutes (PRD §3, CLAUDE.md "no silent surgery")

`mixeff` deliberately replaces several lme4 diagnostic verbs with a single
audit-first channel:

- **Singularity / rank** → `is_singular()`, the `effective_covariance` artifact
  (`requested_rank` vs `supported_rank`, per-direction loadings,
  `unsupported_directions`), `changes()`, and the print "Fitted covariance state"
  block. This covers `isSingular` (works) and the *intent* of `rePCA` (partial:
  the information exists, the prcomp object does not).
- **Convergence / Hessian** → `optimizer_certificate()` with derivative-backed
  `free_gradient_norm`, `projected_gradient_norm`, `hessian_eigen_min`,
  `hessian_rank`. This is a genuine post-optimization KKT/Hessian check
  (cf. crate-1.md:135,146) and is *better surfaced* than lme4's terse warnings.
  It is classified **partial** vs `checkConv`/`optinfo` only because the full
  gradient vector and Hessian matrix are not handed back, and there are no
  user-tunable check thresholds.

These substitutes mean the family is **not** "missing" wholesale — the core
singularity and convergence questions are answered, often more clearly than
lme4. The genuine in-scope gaps are the influence/leverage cluster and residual
types.

### Cross-API behavior (clearer-errors standard)

Calling `lme4::isSingular(fm)`, `lme4::rePCA(fm)`, `hatvalues(fm)`,
`influence(fm)` on an `mm_lmm` produces the generic R error
`no applicable method for '...' applied to an object of class
"c('mm_lmm','mm_fit','mm_compiled')"`. That is not a *clear* mixeff error — a
user migrating from lme4 gets a base-R dispatch failure, not a pointer to the
mixeff equivalent (`is_singular()`, `optimizer_certificate()`, `changes()`).
This is a minor but real "clearer-errors" miss for the absent verbs.

### Highest-impact in-scope gaps a real lme4 user will hit

1. **Influence / leverage cluster** (`influence.merMod`, `cooks.distance`,
   `dfbeta`, `dfbetas`, `hatvalues`) — entirely absent, in-scope, **major**.
   Per-group deletion influence is a routine LMM diagnostic.
2. **Residual types** (`pearson`, `deviance`, `working`) — only `response` is
   supported; for GLMM the lme4 default (`deviance`) and the standard diagnostic
   (`pearson`) both error. **major**.
3. **Full gradient/Hessian access** for external recheck (`getME(.,"devfun")`,
   `optinfo$derivs$Hessian`, `rcond`) — only summary scalars are exposed.
   **major**.
4. **`allFit`** multi-optimizer gold-standard check — absent, **major**.
5. **`getME` boundary diagnostic** (`"lower"`, so `theta == lower`) — absent,
   **major**.

### Out-of-scope / deferred (do not count as defects)

- Profile-likelihood CIs for **GLMM** (PRD §3 explicit non-goal). LMM profile CIs
  appear to be in progress (untracked `test-confint-profile.R`).
- R-level optimizer functions (`Nelder_Mead`, `nloptwrap`, `nlminbwrap`):
  optimization lives in the Rust crate by design (PRD §4).
- `fortify.merMod`: deprecated upstream.

---

## Verification commands used

```r
library(mixeff); library(lme4); data(sleepstudy)
fm <- lmm(Reaction ~ Days + (Days|Subject), data=sleepstudy)
fit_status(fm)                       # "converged_interior"
is_singular(fm)                      # FALSE
optimizer_certificate(fm)            # status, gradient norms, hessian_eigen_min/rank
getME(fm, "theta")                   # OK; "lower"/"devfun"/"ST"/"L"/"par" -> error
ranef(fm, condVar=TRUE)              # postVar 2x2x18 present (LMM)
residuals(fm, type="pearson")        # ERROR: 'arg' should be "response"
# reduced-rank:
ss <- sleepstudy; ss$noise <- rnorm(180)
fs <- lmm(Reaction ~ Days + (Days+noise|Subject), data=ss)
is_singular(fs)                      # TRUE; fit_status "converged_reduced_rank"
changes(fs)                          # effective covariance rank 2 < requested 3
# absent verbs (all "no applicable method"):
lme4::isSingular(fm); lme4::rePCA(fm); hatvalues(fm); influence(fm)
```
