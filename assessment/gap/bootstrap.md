# Gap Analysis — Capability Family: Bootstrap

Date: 2026-05-31
Reference surface: `assessment/survey/lme4-bootstrap.md`
Packages: mixeff (installed), lme4 2.0.1, lmerTest 3.2.1, pbkrtest, boot.

## Scope of this family

lme4/pbkrtest bootstrap surface: `bootMer()` (generic FUN-based parametric &
semiparametric bootstrap), `confint(method="boot")`, `simulate.merMod()`,
`.simulateFun()`, `refit()`, `PBmodcomp()`/`PBrefdist()`/`seqPBmodcomp()`,
`boot`-class integration (`boot.ci`, `plot.boot`), `.progress`/`PBargs`,
`parallel`/`ncpus`/`cl`, `cluster.rand`, `newparams`.

mixeff exports in this family (confirmed via `ls("package:mixeff")`):
`bootstrap_control`, `parametric_bootstrap`, `refit`, plus S3 methods
`simulate.mm_lmm`, `confint.mm_lmm` (with `method="bootstrap"`), and the
internal `test_effect(method=c("bootstrap","bootstrap_lrt","cluster_bootstrap"))`
/ `contrast(method="bootstrap")` paths. No `bootMer`, no `PBmodcomp`.

---

## Capability-by-capability table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `bootMer(x, FUN, ...)` — generic bootstrap of an arbitrary user `FUN(merMod)` statistic | **Absent.** No `bootMer`; no exported function accepts a user `FUN`. Bootstrap is only available through fixed-purpose entry points (fixed-effect CIs, term tests, LRT). | in-scope-missing | major | `exists("bootMer", where="package:mixeff")` → FALSE. `formals(parametric_bootstrap)`/`bootstrap_control`/`confint.mm_lmm` have no `FUN` argument. A user wanting to bootstrap, e.g., an R²/ICC/derived quantity has no path. |
| `confint(merMod, method="boot")` — bootstrap CIs (default scope = variance components + fixed effects) | **Partial.** `confint(fit, method="bootstrap", bootstrap=bootstrap_control(...))` works but covers **fixed effects only**. lme4's default bootstrap CI also covers `.sig01`,`.sig02`,`.sig03`,`.sigma` (variance/correlation/residual). | partial | major | mixeff: `confint(fit, method="bootstrap")` rownames = `(Intercept)`,`Days` only. lme4: `confint(lmer(...), method="boot")` rownames = `.sig01 .sig02 .sig03 .sigma (Intercept) Days`. Variance-component CIs are reachable via `method="profile"` (theta1..theta3, sigma) but NOT via bootstrap. |
| `confint` `boot.type=` (`perc`,`basic`,`norm`) | **Partial.** mixeff `interval=` accepts only `"percentile"` and `"basic"`. `"norm"` (normal-approx) is missing; `"bca"`/`"stud"` also missing (lme4 lacks bca/stud too). | partial | minor | `confint(fit, method="bootstrap", interval="norm")` → error `'arg' should be one of "percentile","basic"`. |
| `confint(method="boot")` on **GLMM** | **Broken / absent.** No `confint.mm_glmm`; dispatch falls through to default `confint`, then errors. | in-scope-missing | major | `glmm(y~x+(1|g), binomial())` then `confint(g, method="bootstrap")` → `"non-numeric argument to binary operator"`; `confint(g)` same. `exists("confint.mm_glmm", ns)` → FALSE. This is an inscrutable error, violating the "clearer errors" mandate. |
| `simulate.merMod()` — response simulation, `nsim`, `seed`, `re.form` (NULL/NA), `use.u` | **Partial.** `simulate.mm_lmm` works for `nsim`/`seed`/`re.form` (NULL→conditional, NA→population). Returns a data frame of sims. Missing: `use.u`, `newdata`, `newparams`, `family`, `cluster.rand`, `na.action`. | partial | major | `simulate(fit, nsim=3, seed=1)` returns 180×3 df with `mm_method="r_side_gaussian_parametric"`. `formals(simulate.mm_lmm)` = object,nsim,seed,re.form,... — no `newdata`/`newparams`/`use.u`. `re.form` beyond NULL/NA refused with `mm_inference_unavailable`. |
| `simulate.merMod()` on **GLMM** | **Absent.** No `simulate.mm_glmm`. | in-scope-missing | major | `exists("simulate.mm_glmm", ns)` → FALSE; `simulate(g)` for an `mm_glmm` does not dispatch to a mixeff method. |
| `simulate` `newparams=` — simulate from non-fitted theta/beta/sigma | **Absent.** | in-scope-missing | minor | Not in `simulate.mm_lmm` formals. |
| `simulate` `cluster.rand=` — non-normal random-effect simulation (sensitivity analyses) | **Absent.** RE draws are hard-wired Gaussian (`mm_rmvnorm`). | in-scope-missing | minor | `mm_simulate_random_mean` uses `mm_rmvnorm` only; no hook to substitute a generator. |
| `simulate` `newdata=` — out-of-sample simulation | **Absent.** | in-scope-missing | minor | Not in formals; sims always use stored `model_frame`. |
| `.simulateFun()` — extended programmatic simulation engine | **Absent (internal-only in lme4).** | out-of-scope-by-design | cosmetic | lme4 exports it but it is a power-user internal; mixeff's `mm_simulate_once` is the analog (not exported). Not a user-facing gap. |
| `refit(object, newresp)` (LMM) | **Works.** `refit.mm_lmm` refits via `lmm()` with stored model frame and REML setting; validates `newresp` length/NA. | works | — | `refit(fit, newresp=...)` returns a new `mm_lmm`; rejects wrong-length/NA response with a clear `mm_arg_error`. |
| `refit(object, newresp)` on **GLMM** | **Absent.** No `refit.mm_glmm`. | in-scope-missing | major | `refit` is generic but only `refit.mm_lmm` exists; `refit(g, ...)` for `mm_glmm` has no method. |
| `refit` `newoffset=` | **Absent.** | in-scope-missing | minor | `formals(refit)` = object,newresp,...; no `newoffset`. |
| Parametric bootstrap LRT for nested model comparison (`PBmodcomp` `PBtest` column) | **Works (LMM, ML).** `parametric_bootstrap(null, alt, nsim, seed)` and `compare(method="bootstrap")` return a certified p-value with replicate accounting (successful/completed/boundary/MCSE/seed). | works | — | `parametric_bootstrap(null, fit, nsim=50, seed=1)` → status available, observed 23.54, p.value 0, MCSE 0, 50/50 reps. Backed by Rust `mm_bootstrap_lrt_json`. |
| `PBmodcomp` additional columns: `Bartlett`, `Gamma`, `F` approximations | **Absent.** mixeff returns the raw bootstrap (`PBtest`-equivalent) p-value only. | in-scope-missing | minor | `print.mm_parametric_bootstrap` reports observed/p.value/MCSE/reps only; no Bartlett/Gamma/F. These are pbkrtest conveniences for accelerating/smoothing the same reference distribution. |
| `PBmodcomp` flexible `smallModel` (term string / `~.-term` / contrast matrix L) | **Partial.** `parametric_bootstrap()` requires two fitted `mm_lmm` objects; `test_effect(method="bootstrap_lrt", term=)` covers the drop-a-term case. No formula/L form. | partial | minor | `parametric_bootstrap` formals require fitted `null`,`alternative`; non-`mm_lmm` → `mm_arg_error`. Term-drop covered by `test_effect`. |
| `PBmodcomp` on **GLMM** | **Absent (engine LMM-only).** | upstream-blocked | major | Crate contract (`crate-9.md:24`, `crate-2.md:153`): `stats::parametric_bootstrap_lrt` is **LMM only**. GLMM parametric-bootstrap LRT not provided by the engine. A real lmer→glmer user comparing nested GLMMs via PB has no route. |
| `PBrefdist()` — pre-compute & reuse reference distribution across comparisons | **Absent.** `parametric_bootstrap` re-simulates each call; the `simulated` replicate vector is on the result but cannot be fed back as `ref`. | in-scope-missing | minor | `formals(parametric_bootstrap)` has no `ref`. Result object stores `simulated` but no API to reuse it. |
| `seqPBmodcomp()` — adaptive early-stopping bootstrap | **Absent.** (Also not exported in installed pbkrtest.) | out-of-scope-by-design | cosmetic | Not exported by pbkrtest on this system; an efficiency convenience, not a capability gap. |
| **Semiparametric bootstrap** (`type="semiparametric"`, residual resampling) | **Absent.** No `type=` argument anywhere; all bootstrap paths are parametric. | in-scope-missing | minor | No `type` in `bootstrap_control`/`parametric_bootstrap`/`confint`. lme4's own semiparametric is "experimental/partial," so the practical loss is limited, but it is a named lme4 capability that is fully absent. |
| `type="parametric", use.u=TRUE` — conditional parametric bootstrap (u fixed) | **Absent.** No `use.u`/conditional-bootstrap toggle. | in-scope-missing | minor | No `use.u` in any bootstrap entry point. (Marginal parametric bootstrap is the implemented mode.) |
| `boot`-class return / `boot::boot.ci`, `plot.boot`, `as.data.frame` interop | **Absent.** `confint(method="bootstrap")` returns `c("mm_confint","matrix")`, not a `"boot"` object. `parametric_bootstrap` returns `mm_parametric_bootstrap`. | in-scope-missing | major | `class(confint(fit, method="bootstrap"))` = `mm_confint`,`matrix`. No `boot.ci`/`plot.boot` interop. lme4 users routinely pipe `bootMer()` output into `boot::boot.ci(type=c("norm","basic","perc"))` with transforms `h`/`hinv` — entirely unavailable. (mixeff does expose the raw replicate payload via `attr(x,"bootstrap")` and `result$simulated`, so a determined user can reconstruct CIs manually.) |
| `.progress` (`txt`/`tk`/`win`) progress bars | **Absent.** No `.progress` argument. | in-scope-missing | minor | Not in any formals. For long bootstraps (nsim≥1000) lme4 users rely on this; the Rust engine runs the loop, so a progress callback would need plumbing. |
| `PBargs` (progress-bar constructor args) | **Absent.** | in-scope-missing | cosmetic | Tied to `.progress`; not present. |
| `parallel`/`ncpus`/`cl` (multicore/snow parallel bootstrap) | **Absent at the R surface.** No `parallel`/`ncpus`/`cl` arguments; parallelism (if any) is internal to the Rust engine and not user-controllable. | in-scope-missing | minor | Not in `bootstrap_control` or any bootstrap entry point. Bootstrap speed is the headline selling point ("faster"), so the Rust engine may obviate user parallel control — but there is no documented knob and no `cl` for reproducible cluster setups. |
| `seed` control for reproducibility | **Works.** `bootstrap_control(seed=)`, `parametric_bootstrap(seed=)`, `simulate(seed=)` all honor and record the seed; engine returns `seed_record`. | works | — | `bootstrap_control(nsim=50, seed=1)` → result reports `seed: 1`; reproducible across runs. `mm_with_seed` preserves/restores `.Random.seed`. |
| Replicate failure accounting (`bootFail`, `boot.fail.msgs`) | **Works (better than lme4).** Engine returns `successful_replicates`/`completed_replicates`/`failed_refits`/`boundary_count`/`boundary_rate`; `failed_refit_policy` is user-controllable (`exclude`/`count_extreme`/`abort`). | works | — | `confint(...,method="bootstrap")` print shows requested/successful/failed_refits/boundary_rate per parameter. This is a structured-diagnostic improvement over lme4's attributes (consistent with PRD audit-first design). |
| `cluster_bootstrap` (estimator distribution) | **Partial by design.** `test_effect(method="cluster_bootstrap")` returns an unavailable table with a stable reason code; p-values `not_assessed` in schema 1.0.0. | out-of-scope-by-design | minor | `crate-8.md:30`: cluster bootstrap is an estimator-distribution target; term p-values deliberately `not_assessed`. Documented refusal, not a defect. |

---

## Classification summary

- **works**: `refit` (LMM), `parametric_bootstrap` LRT (LMM/ML), seed control,
  replicate failure accounting. Core LRT-bootstrap use case and reproducibility
  are solid and arguably better-instrumented than lme4.
- **partial**: `confint(method="boot")` (fixed effects only, no variance-component
  CIs), `boot.type`/`interval` (no `norm`), `simulate` (no use.u/newdata/newparams/
  cluster.rand), flexible `smallModel`.
- **in-scope-missing (the real lme4-user pain points)**: generic `bootMer`+`FUN`;
  bootstrap CIs for variance components; GLMM `confint`/`refit`/`simulate`;
  `boot`-class interop (`boot.ci`/`plot.boot`); semiparametric & conditional
  (`use.u`) bootstrap; `PBrefdist` reuse; `.progress`/`PBargs`; user parallel
  control.
- **upstream-blocked**: parametric-bootstrap LRT for GLMM (engine is LMM-only per
  crate contract).
- **out-of-scope-by-design**: `.simulateFun` (internal), `seqPBmodcomp`,
  `cluster_bootstrap` p-values (documented `not_assessed`).

## Notable observations

1. **GLMM is the biggest hole.** `confint`, `refit`, `simulate`, and PB-LRT all
   work for `mm_lmm` and are absent or broken for `mm_glmm`. The GLMM
   `confint(method="bootstrap")` path produces `"non-numeric argument to binary
   operator"` — an inscrutable base-R error that directly violates the project's
   "clearer errors" mandate. This should at minimum be a clean
   `mm_inference_unavailable` refusal.
2. **Variance-component bootstrap CIs are the second hole.** lme4's
   `confint(method="boot")` *defaults* to including `.sig01/.sigma`; mixeff's
   bootstrap CI silently covers fixed effects only. Profile-LL CIs fill this for
   LMM (theta/sigma), so a workaround exists, but the bootstrap surface is not at
   parity. (Note: PRD §3 lists profile-likelihood CIs as a v0 non-goal, yet
   `confint(method="profile")` *is* implemented and certified — the crate contract
   `crate-8.md:159` lists profile-LL CI as release-blocking/certified, so §3's
   bullet appears stale rather than the feature being missing.)
3. **No `bootMer`/FUN generality.** The single most flexible lme4 bootstrap
   primitive — bootstrap any statistic — has no analog. Everything in mixeff is a
   fixed-purpose endpoint.
4. **No `boot`-class interop.** lme4 leans on the `boot` ecosystem
   (`boot.ci`, transforms, plots). mixeff returns its own classes; raw replicates
   are recoverable via `attr(x,"bootstrap")` / `result$simulated` but require
   manual reconstruction.

## Verdict

Parity is **partial**. The certified LMM LRT-bootstrap and reproducibility/
diagnostics are at or above lme4 quality, but a real lme4 bootstrap user hits
multiple hard walls: no generic `bootMer`, no variance-component bootstrap CIs,
no GLMM bootstrap (and an inscrutable error there), no `boot`-class interop, and
no semiparametric variant.
