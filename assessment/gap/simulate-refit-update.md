# Gap Analysis: `simulate` / `refit` / `update` family

Survey date: 2026-05-31
mixeff branch: main
lme4 2.0.1, lmerTest 3.2.1 (installed)
Reference surface: `assessment/survey/lme4-simulate-refit-update.md`

## Scope context

- PRD §10 Phase 4 explicitly lists this family as **in scope**: "Phase 4 — GLMM
  + simulate + refit + compare … `glmm()`, `simulate`, `refit`, `compare`,
  multi-model `anova`, `drop1`, parametric bootstrap."
- PRD §3 non-goals relevant here: **"No silent model surgery; reductions and
  refusals always cross the boundary."** This makes any silently-ignored
  argument a contract violation, not a missing feature. PRD §3 does NOT list
  `newparams`, `newdata`-simulation, `bootMer`, `allFit`, `refitML`,
  `simulate.formula`, or `cluster.rand` as explicit non-goals.
- `cbind(y1, y2) ~` multivariate / 2-column binomial response is out of scope
  per §3 (multivariate cross-outcome) and is rejected by the engine's stateless
  transform subset.

## mixeff inventory (R/)

- `R/simulate.R`: `refit()` generic + `refit.mm_lmm`; `simulate.mm_lmm`.
- `R/compare.R`: `parametric_bootstrap(null, alternative, nsim, seed)` +
  `print.mm_parametric_bootstrap`; `bootstrap_control()` in `R/inference.R`.
- NAMESPACE exports: `refit`, `simulate.mm_lmm`, `parametric_bootstrap`,
  `bootstrap_control`, `importFrom(stats, update)`. There is **no**
  `update.mm_lmm`, `simulate.mm_glmm`, `refit.mm_glmm`, `refitML`, `bootMer`, or
  `allFit`.

## Capability table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `simulate.merMod` basic (`nsim`, `seed`) | works | works | — | `simulate(fit, nsim=2, seed=1)` returns 180×2 data.frame, cols `sim_1`,`sim_2`, `attr(,"seed")` set. RNG state restored via `mm_with_seed()` (`R/simulate.R:189`). |
| `simulate` `re.form = NULL` (conditional, default) | works | works | — | Default `re.form=NULL` → conditional draw via `mm_simulate_random_mean()`. |
| `simulate` `re.form = NA` / `~0` (population) | partial | partial | minor | `re.form=NA` works (population mean + residual). But `~0` (formula form) is not handled — only literal `NULL`/`NA` map to targets via `mm_prediction_target()`. lme4 accepts `~0`. |
| `simulate` partial `re.form` (e.g. `~(Days\|Subject)`) | errors (typed) | in-scope-missing | major | `simulate(fit, re.form=~(Days\|Subject))` → `mm_inference_unavailable`: "`re.form` requests beyond NULL and NA are not available for simulation." Clear refusal (good), but the feature a real user needs (condition on some REs, simulate others) is absent. |
| `simulate` `newparams` (hypothetical θ/β/σ) | **silently ignored** | in-scope-missing (contract violation) | **blocker** | `simulate(fit, seed=42, newparams=list(beta=c(9999,9999),theta=10,sigma=999))` returns output **identical** to default-param simulate (mean 300.63 both). Absorbed by `...`, never read. This is the primary power-analysis path and is silently dropped — direct violation of PRD §3 "no silent surgery". |
| `simulate` `newdata` (out-of-sample simulation) | **silently ignored** | in-scope-missing (contract violation) | **blocker** | `simulate(fit, newdata=sleepstudy[1:20,])` returns **180** rows, not 20. `newdata` absorbed by `...` and never used. Note: `predict.mm_lmm` *does* honor `newdata`, so the asymmetry is a real surprise. |
| `simulate` `use.u` (legacy alias for `re.form`) | **silently ignored** | partial (contract violation) | major | `simulate(fit, use.u=TRUE)` accepted with no error and no effect; `use.u` is absorbed by `...`. lme4 maps `use.u=TRUE`↔`re.form=NULL`, `FALSE`↔`re.form=~0`. mixeff ignores it. |
| `simulate` `family` override | missing | in-scope-missing | minor | Not handled (LMM Gaussian only); see GLMM rows. |
| `simulate` `cluster.rand` (non-normal RE draws) | missing | in-scope-missing | minor | `mm_rmvnorm()` hardcodes Gaussian RE draws; no hook to swap the cluster generator. Real users doing RE-robustness studies hit this. Not a §3 non-goal. |
| `simulate` `allow.new.levels` | missing (moot) | in-scope-missing | minor | Since `newdata` is ignored, `allow.new.levels` for simulation is also absent. (`predict.mm_lmm` supports it for prediction.) |
| `simulate.mm_glmm` (GLMM responses) | **absent** | in-scope-missing | **major** | `getS3method("simulate","mm_glmm")` → none. `simulate(glmm_fit)` → "no applicable method for 'simulate' applied to … 'mm_glmm'". PRD §10 Phase 4 pairs GLMM with simulate explicitly; crate-4.md:190 notes "`simulate.mm_lmm` — only wired for LMM". |
| `simulate.formula` / `.simulateFun` (no fitted model) | absent | in-scope-missing | major | No `simulate` formula method. The fit-free data-generating path for power analysis is unavailable. Not a §3 non-goal. |
| `simulate` `cond.sim`, `weights`, `offset` (formula path) | absent | in-scope-missing | minor | Only reachable via the absent `simulate.formula` path. |
| `refit.merMod` numeric `newresp` | works | works | — | `refit(fit, sims[[1]])` returns `mm_lmm`; re-runs `lmm()` with stored model frame + REML setting (`R/simulate.R:19`). Records `$refit` provenance. |
| `refit` accepting a `simulate()` 1-col data.frame directly | errors | partial | major | lme4 idiom `lapply(sims, refit, object=fit)` works because `refit.merMod` accepts a single-column df. `refit(fit, sims[,1,drop=FALSE])` → "`newresp` must be a numeric vector …". Users must `sims[[1]]` / `unlist`. Breaks the documented lme4 bootstrap loop verbatim. |
| `refit` NA round-trip (`na.action` attr on newresp) | not handled | partial | minor | `refit.mm_lmm` rejects any `newresp` with `anyNA` and requires `length == nobs`; no `na.action`-attr shortening path. lme4 supports passing `simulate()` output from an NA-excluded fit. |
| `refit` `newweights` | **silently ignored** | partial (contract violation) | major | `refit(fit, y, newweights=rep(1,180))` runs with no error; `newweights` absorbed by `...`. `refit.mm_lmm` always reuses `object$weights`. lme4 lets you change weights. Silent drop = §3 violation. |
| `refit` `rename.response` | missing | in-scope-missing | cosmetic | Not supported. Low real-world impact. |
| `refit` `maxit` (GLMM inner loop) | missing (moot) | in-scope-missing | minor | No `refit.mm_glmm` at all. |
| `refit.mm_glmm` | **absent** | in-scope-missing | **major** | `getS3method("refit","mm_glmm")` → none; `refit(glmm_fit, y)` → "no applicable method". PRD §10 Phase 4 lists GLMM+refit. |
| `refitML` (REML→ML refit) | absent (as a verb) | partial | minor | No exported `refitML`. The capability exists internally: `compare()`/multi-model `anova()` refit REML→ML for fixed-effect comparison (`R/compare.R`, test-phase4.R:205), and `boundary_lrt`/`parametric_bootstrap` refuse REML with a typed reason. No standalone user-facing `refitML(fit)`. |
| `update.merMod` formula./args | works (via `stats::update.default`) | works | — | No `update.mm_lmm` method exists (`getS3method` → none); `update(fit, .~.-Days)` succeeds through `stats::update.default` re-running the stored `lmm()` call. Stored `$call` is `lmm(formula=…, data=…)`. |
| `update` `evaluate = FALSE` | works | works | — | `update(fit, .~.-Days, evaluate=FALSE)` returns the unevaluated `lmm(...)` call. |
| `update` changing `data`/`REML`/`control`/`weights` | works (via default) | works | — | Inherited from `stats::update.default` since all are `lmm()` args. Not separately tested. |
| `update.lmerModLmerTest` class preservation | N/A (out-of-scope-by-design) | out-of-scope-by-design | — | mixeff has a single `mm_lmm` class with Satterthwaite/inference always available; there is no two-tier `lmerMod`→`lmerModLmerTest` upgrade to preserve. lme4 §8 behaviour is moot here. |
| `bootMer` (FUN-driven (semi)parametric bootstrap) | absent | in-scope-missing | major | `bootMer(fit, fixef, nsim=2)` → "no applicable method for 'isLMM' applied to … 'mm_lmm'" (mixeff objects don't satisfy lme4's `isLMM`). No native equivalent for arbitrary `FUN`. `parametric_bootstrap()` only does the nested-model LRT, and `compare(method="bootstrap")` / `confint(method="bootstrap")` cover specific statistics — not a general `FUN` harness. `type="semiparametric"` (residual resampling) absent. Not a §3 non-goal. |
| `bootMer` `parallel`/`ncpus`/`cl` | absent (moot) | in-scope-missing | minor | No `bootMer`; bootstrap parallelism is engine-internal. |
| `allFit` (multi-optimizer convergence check) | absent | in-scope-missing | minor | `allFit(fit)` → "no applicable method for 'isGLMM' applied to … 'mm_lmm'". No native cross-optimizer robustness harness. Arguably tension with §3 "no model-selection" but allFit is diagnostic, not selection. |
| `parametric_bootstrap()` nested-model LRT (mixeff extra) | works | works | — | `parametric_bootstrap(reduced, full, nsim=2, seed=101)` returns `mm_parametric_bootstrap`; routes through certified Rust LRT; refuses REML pairs and mismatched data with typed reasons (test-phase4.R:284–365). This is the mixeff-native counterpart to the bootMer-for-LRT idiom. |

## Key findings

### Blockers / majors (in-scope, lme4 users will hit)

1. **`simulate(newparams=)` silently ignored** (blocker). The canonical power-
   analysis path produces output identical to default-parameter simulation with
   no error or warning. Either honor `newparams` or refuse it with a typed
   `mm_inference_unavailable`/`mm_arg_error`. Current behaviour violates PRD §3
   "no silent model surgery". Repro:
   `simulate(fit, seed=42, newparams=list(beta=c(9999,9999),theta=10,sigma=999))`
   equals `simulate(fit, seed=42)`.

2. **`simulate(newdata=)` silently ignored** (blocker). Returns full-data-length
   vectors regardless. Asymmetric with `predict.mm_lmm`, which honors `newdata`.

3. **No `simulate.mm_glmm` / `refit.mm_glmm`** (major). PRD §10 Phase 4 pairs
   GLMM with simulate+refit; both methods are absent and error with the raw
   S3-dispatch message ("no applicable method"), not a typed mixeff diagnostic.

4. **`refit` rejects the `simulate()` data.frame column form** (major). The
   documented lme4 loop `lapply(sims, refit, object=fit)` fails; users must
   extract `sims[[i]]`.

5. **`refit(newweights=)` and `simulate(use.u=)` silently ignored** (major,
   contract violations). Absorbed by `...` with no effect.

6. **No general `bootMer(FUN=)` harness** (major). Only the LRT-specific
   `parametric_bootstrap()` and method-keyed `compare`/`confint` bootstraps
   exist; arbitrary user statistics and semiparametric (residual-resample)
   bootstrap are unavailable.

### Partials / minors

- `simulate` `re.form`: only `NULL`/`NA` literals; `~0` and partial formulas
  unsupported (partial refusal is typed and clear).
- `refitML` exists only as internal compare/anova behaviour, not a user verb.
- `refit` NA round-trip, `rename.response`, `cluster.rand`, `family` override,
  `allFit`: in-scope-missing, lower frequency.

### Works

- `simulate.mm_lmm` (LMM, `nsim`/`seed`/conditional/population), `refit.mm_lmm`
  (numeric newresp), `update` (via `stats::update.default`, incl.
  `evaluate=FALSE`), and the mixeff-native `parametric_bootstrap()` LRT.

### Out-of-scope-by-design

- `update.lmerModLmerTest` class-upgrade behaviour (mixeff has one `mm_lmm`
  class; inference always available — §3 design).
- `cbind(y1,y2)~` 2-column binomial simulate (multivariate, §3).

### Test gaps

- test-phase4.R:189 covers `simulate`/`refit` happy paths and the LRT bootstrap,
  but **no test asserts** that `newparams`/`newdata`/`use.u`/`newweights` are
  rejected (or honored). The silent-ignore bugs are unguarded — adding refusal
  tests would have surfaced findings 1, 2, 5.

## Evidence files
- Source: `R/simulate.R`, `R/compare.R` (`parametric_bootstrap`),
  `R/inference.R` (`bootstrap_control`), NAMESPACE.
- Repro scripts run live against installed `mixeff` + `lme4` 2.0.1 (see table
  rows for exact calls and observed output).
- Tests: `tests/testthat/test-phase4.R:189,205,284–365`.
