# Gap Analysis — Inference: confint / anova / drop1 / profile

**Capability family:** `confint` (profile/Wald/boot, `parm`), `profile()`, `anova`
(single & multi-model REML-refit + LRT), `drop1`, `KRmodcomp`, `PBmodcomp`,
plus the pbkrtest/lmerTest companions (`ranova`, `contest*`, `SATmodcomp`,
`vcovAdj`, `show_tests`, `bootMer`, `PBrefdist`).

**Date:** 2026-05-31
**Reference:** `assessment/survey/lme4-inference-confint-anova.md`
**Environment:** mixeff (installed), lme4 2.0.1, lmerTest 3.2.1, pbkrtest, emmeans 2.0.3.
All reproductions run against `lme4::sleepstudy`, `Reaction ~ Days + (Days | Subject)`.

## Verdict

**Partial.** The statistically load-bearing surface is present and behaves
correctly where implemented: Wald/profile/bootstrap CIs for fixed effects,
single- and multi-model `anova`, `drop1`, asymptotic LRT with REML→ML auto-refit,
a parametric-bootstrap LRT (the `PBmodcomp` analog), a scalar Kenward-Roger F-test,
and a boundary-aware random-effect LRT (the `ranova` analog). The honesty contract
is honored throughout — every uncertified path emits a structured refusal rather
than a fabricated number, which is exactly the design intent.

But several capabilities a real lme4 user routinely reaches for are missing or
narrower than lme4:
- No `profile()` generic / `thpr` object (no reusable profile, no `xyplot`/`splom`,
  no `logProf`/`varianceProf`). Phase 5 scope per PRD §10.
- Profile and bootstrap CIs do **not** cover the full lme4 parameter set: bootstrap
  CIs are fixed-effects-only; profile CIs omit beta under REML by contract.
- No `parm` group selectors (`"theta_"`, `"beta_"`), no `boot.type="norm"`,
  no `signames`, no `boot.scale="vcov"`.
- No public Kenward-Roger / parametric-bootstrap model-comparison verbs matching
  `KRmodcomp`/`PBmodcomp` signatures (only a scalar-term KR F and a `compare()`
  bootstrap LRT exist), and KR is documented as out-of-scope-beyond-scalar.

## Capability table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `confint(method="Wald")` (fixed effects) | works | works | — | `confint(m, method="wald")` returns `2.5%/97.5%` for `(Intercept)`, `Days`; `method="asymptotic"` is an accepted synonym. Carries `status: not_certified_by_rust_inference_contract` (honesty marker, by design). `R/inference.R:844-896`. |
| `confint(method="Wald")` for VC params returns `NA` | works (different shape) | works | minor | lme4 returns `NA` rows for VC params; mixeff simply restricts Wald to fixed effects and errors on VC names (`confint(m, parm="sigma")` → "Unknown fixed-effect parameter(s): sigma"). Behavior is defensible but not row-for-row identical. |
| `confint(method="profile")` — fixed effects | partial | partial | major | Under **ML** beta profile CIs are returned (`(Intercept)`, `Days` both finite). Under **REML** beta rows are explicit typed refusals (`reason_code=profile_beta_unavailable_under_reml`, lower/upper `NA`) by upstream contract. lme4 profiles beta under REML. `R/inference.R:1646-1712`, `mm_lmm_profile_confint_json` FFI. |
| `confint(method="profile")` — VC params (theta/sigma) | works | works | — | Returns finite CIs for `theta1..theta3`, `sigma`; near-boundary `theta2` correctly returns `NA` lower with finite upper. `method: profile_likelihood / status: available`. |
| `confint(method="boot")` — fixed effects | works | works | — | `confint(m, method="boot")` and `="bootstrap"` both run a full-model parametric bootstrap (999 reps default), report percentile CIs plus a run summary (successful/failed refits, boundary rate) and the honesty note that distributions "do not certify fixed-effect hypothesis-test p-values". `R/inference.R:1274-1302`. |
| `confint(method="boot")` — VC params (SD/cor) | in-scope-missing | in-scope-missing | major | lme4 default bootstrap CI covers fixed effects **and** all VC parameters on sd/cor scale. mixeff rejects VC param names: `confint(m, parm="sigma", method="boot")` → "Unknown fixed-effect parameter(s): sigma." A real user expecting bootstrap CIs on random-effect SDs hits a wall. |
| `confint(parm=)` integer/name subsetting (fixed effects) | works | works | — | `confint(m, parm="Days")` and integer indices map through `names(object$beta)`. `R/inference.R:867-882`. |
| `confint(parm="theta_"/"beta_")` group selectors | in-scope-missing | in-scope-missing | minor | `confint(m, parm="theta_", method="profile")` returns an **empty** matrix (the literal string is not recognized as a group selector; it matches no row). lme4 documents `"theta_"`/`"beta_"` as first-class group selectors. |
| `confint(... boot.type=, boot.scale=, FUN=, signames=, zeta=)` | in-scope-missing | partial | minor | mixeff `confint.mm_lmm` supports `interval=c("percentile","basic")` only; no `"norm"`, no `boot.scale="vcov"`, no user `FUN`, no `signames`, no `zeta` cutoff override. `R/inference.R:844-848`. |
| `profile.merMod()` generic → `thpr` object | in-scope-missing | in-scope-missing | major | `profile(m)` → "no applicable method for 'profile' applied to an object of class mm_lmm". No `thpr` object, so the documented heavy-use workflow `pp<-profile(m); confint(pp, level=...)` is impossible; every `confint(method="profile")` re-profiles from scratch. PRD §10 schedules `profile()` in **Phase 5** (not yet shipped). |
| `confint.thpr()` (CI from stored profile) | in-scope-missing | in-scope-missing | major | Depends on `thpr`; absent (see above). The intentional profile/confint separation that lme4 users exploit for repeated levels is unavailable. |
| `log()/logProf()`, `varianceProf()`, `xyplot.thpr()`, `splom.thpr()` | in-scope-missing | in-scope-missing | minor | No profile object means none of the profile transform/plot utilities exist. Visualization of boundary/non-normal VC profiles is unavailable. |
| `anova()` single-model — sequential deviance table | works (richer) | works | — | `anova(m)` returns a Type III fixed-effect table with `num_df/den_df/statistic/p_value/method`; default method `auto`→Satterthwaite, so mixeff actually gives **more** than lme4's no-p-value table (it bakes in the lmerTest-style result). `R/compare.R:267-311`. |
| `anova()` single-model `type=` (I/II/III) | partial | partial | major | `type` argument is accepted (`c("III","II","I")`, default III) and stamped on the output, **but** the printed/returned table is identical across types — the term-level F is the same single-`Days` row regardless of `type`. For a one-term model that is correct; there is no evidence the Rust term table actually recomputes Type I vs II vs III contrasts for multi-term/interaction models. Needs a multi-term repro to fully certify; flagged as partial/test-gap. `R/compare.R:269-302`. |
| `anova(m1, m2, ...)` multi-model LRT | works | works | — | `anova(m0, m)` delegates to `compare()`, auto-refits REML→ML (logLik −887.74/−875.97 are the ML values), returns `df/logLik/AIC/BIC/delta_df/LRT/p_value`, marks model 1 `reference_model`. `R/compare.R:273-282`, `R/compare.R:572-652`. |
| `anova(refit=TRUE)` REML→ML auto-refit | works | works | — | Default `refit_for_comparison="auto"` refits REML fits to ML before comparison; `refit_for_comparison="error"` refuses with a structured `mm_inference_unavailable`. `R/compare.R:488-508`. Surfaces the refit (no silent surgery), exceeding lme4's silent behavior. |
| `anova(model.names=)` | in-scope-missing | partial | cosmetic | No `model.names` argument; models are auto-labeled `m1/m2` with formula strings. Cosmetic relabeling only. |
| `anova(ddf="Satterthwaite")` (lmerTest F-table) | works | works | — | `anova(m, method="satterthwaite")` (the default `auto` resolves to Satterthwaite) returns df + p-value. NOTE: statistic is reported as `t` for single-df terms, not `F` (`statistic_name="t"`); equivalent (F = t²) but not label-identical to lmerTest. |
| `anova(ddf="Kenward-Roger")` (lmerTest F-table) | partial | partial / out-of-scope-by-design | major | `anova(m, method="kenward_roger")` returns a KR result (`den_df=17`, `statistic_name="t"`, `reliability=moderate`). Works for scalar/single-term tests. PRD §3 lists Kenward-Roger as a v0 non-goal and the task guardrails defer "KR beyond scalar"; multi-df / general KR F is **out-of-scope-by-design**, scalar KR is **partial/works**. At singularity it refuses (`kenward_roger_unavailable_at_boundary`, `R/inference.R:1111-1118`). |
| `anova(ddf="lme4")` fallback | partial | partial | minor | `method="none"` produces an unavailable-effect table (no p-values) which is the spiritual equivalent, but there is no literal `ddf="lme4"` passthrough to a raw deviance table. |
| `drop1(test="none")` (AIC-only) | works | works | — | `drop1(m, test="none")` returns `df/logLik/AIC/BIC` per dropped term, `method="none"`. Refits each reduced fixed-effect model, preserves random effects, reports reduced formulas. `R/compare.R:339-404`. |
| `drop1(test="Chisq")` (asymptotic LRT) | works | works | — | `drop1(m, test="Chisq")` adds `LRT` + `p_value` via `pchisq` on the deviance difference (`Days`: LRT 23.54, p 1.23e-6). |
| `drop1(test="user", sumFun=)` (KR/PB hook) | in-scope-missing | in-scope-missing | minor | mixeff `drop1` exposes only `test=c("none","Chisq")`; no `sumFun`/`test="user"` hook to plug KR or parametric bootstrap per term. lme4's documented extensibility point is absent. `R/compare.R:339-343`. |
| `drop1(ddf="Satterthwaite"/"Kenward-Roger")` (lmerTest F-drop1) | in-scope-missing | in-scope-missing | major | mixeff `drop1` has no `ddf` argument and only does LRT/AIC, not the per-term Satterthwaite/KR F-test that lmerTest's `drop1.lmerModLmerTest` provides. A user wanting F-based single-term deletion must instead use `test_effect(fit, term, method=)` (per-term, not a `drop1` table). |
| `KRmodcomp(largeModel, smallModel, betaH, L)` | partial / out-of-scope-by-design | out-of-scope-by-design | major | No `KRmodcomp` verb. Scalar KR is reachable via `test_effect(fit, term, method="kenward_roger")` / `anova(method="kenward_roger")` (single-term F, confirmed working). The full pbkrtest interface — restriction matrix `L`, `betaH`, large-vs-small model comparison, multi-df F — is not provided; KR is a PRD §3 non-goal (deferred beyond scalar). `getKR()` analog: none. |
| `PBmodcomp(largeModel, smallModel, nsim, ref, seed)` | partial | partial | major | The parametric-bootstrap LRT analog **works**: `compare(m0, m, method="bootstrap", nsim=199, seed=1)` returns `method="parametric_bootstrap_lrt"`, a bootstrap p-value, replicate count, and MCSE (`199/199 replicates, MCSE=0`). Gaps vs `PBmodcomp`: (1) only the bootstrap p-value is reported — no Bartlett/Gamma/F reference-distribution rows; (2) no `ref=` reuse (no `PBrefdist` analog); (3) the verb is `compare()`/`parametric_bootstrap()`, not a `PBmodcomp(large,small)` signature; (4) the bootstrap CI/LRT machinery is fixed-effect / model-pair oriented. `R/compare.R:61-90`, exported `parametric_bootstrap`. |
| `seqPBmodcomp()` (early-stop sequential PB) | in-scope-missing | in-scope-missing | minor | No sequential / early-stop parametric bootstrap. Efficiency-only feature. |
| `SATmodcomp()` (Satterthwaite model comparison) | partial | partial | minor | No standalone `SATmodcomp(large, small)` verb; Satterthwaite is available per-term (`test_effect`/`anova(method="satterthwaite")`) and for contrasts (`contrast(..., method="satterthwaite")`, confirmed), but not as a two-model nested comparison. |
| `X2modcomp()` (chi-square model comparison) | works (different verb) | works | minor | The asymptotic-LRT chi-square comparison is exactly what `anova(m0,m)` / `compare(method="lrt")` deliver; just not under the `X2modcomp` name/framework. |
| `vcovAdj()` (KR-adjusted covariance) | in-scope-missing | out-of-scope-by-design | minor | No public KR-adjusted fixed-effect covariance accessor; follows from KR-beyond-scalar being a PRD §3 non-goal. |
| `PBrefdist()` (precompute PB reference) | in-scope-missing | in-scope-missing | minor | No reusable reference-distribution object; each `compare(method="bootstrap")` re-simulates. |
| `ranova()` (random-effect LRT, boundary correction) | works (different verb) | works | minor | `test_random_effect(fit, term, method="boundary_lrt")` provides the boundary-aware (Self-Liang 50:50 mixture) random-effect variance-component LRT with REML→ML auto-refit — the statistical core of `ranova`. `R/inference.R:287-353`. Gaps: not named `ranova`, no whole-model table over all RE terms in one call, no `reduce.terms` (correlated-slope decomposition to `(1|g)+(0+x|g)`). |
| `ranova(reduce.terms=TRUE)` correlated-slope reduction | in-scope-missing | in-scope-missing | minor | `test_random_effect` drops a whole term; it does not first decompose `(x|g)` into independent components for separate tests as `ranova(reduce.terms=TRUE)` does. |
| `contest()/contest1D()/contestMD()` (contrast tests w/ Satterthwaite/KR df) | works (different verb) | works | minor | `contrast(fit, L, method=c("satterthwaite","kenward_roger",...))` tests arbitrary linear contrasts with the requested df method (confirmed: `contrast(m, rbind(c(0,1)), method="satterthwaite")` → estimate, df 16.98, statistic, p). Covers the `contest*` functionality under a different name; single-df statistic labeled `t` rather than `F`. `R/inference.R` `contrast`. |
| `show_tests()` (inspect internal L matrices) | in-scope-missing | in-scope-missing | minor | No `show_tests` analog to reveal the per-term contrast matrices used in the ANOVA decomposition. Reduces auditability of Type II/III tables (mild tension with the audit-first ethos, but lme4-niche). |
| `bootMer()` (low-level PB engine) | partial | partial | minor | No public `bootMer(x, FUN, ...)` accepting an arbitrary statistic; the bootstrap engine is exposed only through `bootstrap_control()` + `confint(method="boot")` / `parametric_bootstrap()` / `contrast(method="bootstrap")`. No `use.u`/`re.form`/`type="semiparametric"` knobs. |

## Notes on classification

- **Profile/bootstrap parameter coverage (the headline gaps).** Two distinct
  holes a real lme4 user hits immediately: (a) `confint(method="boot")` refuses
  VC parameters, and (b) there is no `profile()`/`thpr` object so profiles can't
  be reused or plotted. Both are in-scope-missing (Phase 5 per PRD §10) and rated
  **major** — not blockers, because the fixed-effect paths work and VC profile CIs
  are available via `confint(method="profile")`.

- **Kenward-Roger.** PRD §3 explicitly lists Kenward-Roger as a v0 non-goal and the
  task guardrails defer KR beyond scalar. Scalar/single-term KR F is implemented and
  works (`test_effect`/`anova(method="kenward_roger")`), so the *full* `KRmodcomp`
  interface (restriction matrices, `betaH`, multi-df F, `vcovAdj`) is classified
  **out-of-scope-by-design**, while the scalar path is **partial/works**.

- **REML beta profile refusal** is an upstream-contract decision (the Rust engine
  does not profile beta under REML), surfaced as an explicit typed refusal rather
  than silently. Under ML, beta profile CIs are returned. Classified **partial**
  (the refusal is honest and documented, but it is narrower than lme4, which
  profiles beta under REML).

- **Honesty contract upheld.** Every uncertified path (`bootstrap_not_run`,
  `kenward_roger_unavailable_at_boundary`, `profile_beta_unavailable_under_reml`,
  Wald `not_certified_by_rust_inference_contract`) crosses the boundary as a
  structured diagnostic with a `reason_code`. This matches the PRD's "no silent
  surgery" mandate and is a clear improvement over lme4's silent approximations.

- **Naming divergence.** Much of pbkrtest/lmerTest is present under mixeff verbs
  (`compare` ⊇ `anova`-multimodel/`X2modcomp`/`PBmodcomp`-LRT,
  `test_random_effect` ≈ `ranova`, `contrast` ≈ `contest*`,
  `test_effect`/`anova(method=)` ≈ `SATmodcomp`/scalar-`KRmodcomp`). A user porting
  literal `KRmodcomp`/`PBmodcomp`/`ranova`/`contest` code will not find those names.

- **Single-df statistic labeling.** Satterthwaite/KR single-df tests report
  `statistic_name="t"` rather than lmerTest's `F` (with F = t²). Numerically
  equivalent; cosmetic mismatch worth noting for users comparing tables.

- **Type II/III not exercised here.** `anova(type=)` is accepted but a one-term
  model cannot demonstrate distinct Type I/II/III behavior; whether the Rust term
  table truly recomputes the marginal contrasts for multi-term/interaction models
  is unverified — flagged partial/test-gap, recommend a multi-term repro.

## Reproductions (abbreviated)

```r
library(mixeff); data(sleepstudy, package="lme4")
m  <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy, REML=TRUE)
m0 <- lmm(Reaction ~ 1 + (Days | Subject), sleepstudy, REML=TRUE)

confint(m, method="wald")                 # fixed effects, honesty-marked
confint(m, method="profile")              # theta/sigma finite; beta NA under REML
confint(m, method="boot")                 # fixed effects only
confint(m, parm="sigma", method="boot")   # ERROR: Unknown fixed-effect parameter(s): sigma
confint(m, parm="theta_", method="profile") # empty matrix (no group selector)
profile(m)                                # ERROR: no applicable method for 'profile'
anova(m)                                  # Type III, Satterthwaite p-value
anova(m0, m)                              # multi-model LRT, REML->ML auto-refit
drop1(m, test="Chisq")                    # per-term asymptotic LRT
test_effect(m, "Days", method="kenward_roger")     # scalar KR F (label "t")
test_random_effect(m, "Subject")          # ranova analog, boundary mixture
compare(m0, m, method="bootstrap", nsim=199, seed=1) # PBmodcomp-LRT analog
contrast(m, rbind(c(0,1)), method="satterthwaite")   # contest analog
```
