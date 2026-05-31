# mixeff — Canonical Capability & CRAN-Readiness Assessment

**Date:** 2026-05-31 · **Reference stack:** lme4 2.0.1, lmerTest 3.2.1, emmeans 2.0.3, pbkrtest 0.5.5, mixeff 0.1.0
**Standard applied:** *"Do everything lme4/lmerTest does, but faster and with clearer (non-inscrutable) errors."*
**Evidence base:** 17 capability-gap reports (`assessment/gap/`), 30 empirical parity probes (`assessment/parity/`), 21 error-quality probes (`assessment/errors/`), 15 test-coverage specifications + 8 test-inventory surveys (`assessment/testspec/`, `assessment/survey/`), the live build/test baseline (`assessment/build-test-baseline.log`), and a completeness-critic pass over packaging/FFI/CRAN surfaces.

---

## 1. Executive summary

mixeff today is a **strong, audit-first LMM engine with a real speed advantage and genuinely clearer input-validation errors** — but it is **not yet "everything lme4 does," and it is not yet CRAN-credible** because of a recurring contract violation and an entirely unexamined packaging/FFI risk surface.

**What is true and verified:**
- **LMMs are at statistical parity with lme4.** Across every Gaussian topology probed (random intercept, correlated/uncorrelated slopes, nested, crossed, interactions, weights, ML, singular boundary), fixef and logLik/AIC/BIC match to ~1e-6 or better — well inside the PRD tolerances (`parity/lmm-*.md`). The only LMM divergences are sub-percent optimizer drift in SE/sigma/ranef on small/correlated models, which **vanishes at large N** (`parity/speed-scaling.md`).
- **The speed promise holds.** mixeff is faster on every timed operation, and the advantage *grows with N*: ~8×→11×→13× on fits from N=1k→50k, ~37× on parametric bootstrap, ~3× on profile CIs (`parity/speed-scaling.md`, `parity/inf-bootstrap.md`, `parity/inf-confint-profile.md`).
- **The audit-first refusal machinery is real and frequently clearer than lme4** — typed `mm_*` conditions replace lme4's silent NA-row drops, silent `NA` returns, opaque C++ PIRLS panics, and parser fragments on a broad class of input-validation failures (`errors/`: 8 "clearer-than-lme4" verdicts).

**What undercuts the goal — be candid:**
1. **A cross-cutting "silent-surgery" defect class.** Multiple documented lme4 arguments (`simulate(newparams=/newdata=)`, `predict(random.only=)`, `residuals(scaled=)`, `logLik/deviance(REML=)`, `glmm(offset=)`, non-default factor `contrasts()`, `emmeans(mode=/lmer.df=)`) are silently swallowed by `...` and return **plausible-but-wrong numbers with no signal**. This violates the package's own foundational contract (CLAUDE.md "no silent surgery"; PRD §3, §8.1) and is the project's clearest reputational risk.
2. **One outright wrong-answer blocker:** `emmeans()` / native `mm_means` on **interaction models** returns marginal means evaluated at the reference level instead of averaging (lme4 35.8 vs mixeff 19.3), with sign-flipped contrasts and only a non-committal `reliability="low"` flag (`parity/inf-emmeans.md`).
3. **GLMMs are second-class and sometimes genuinely under-converged.** Beyond the *documented* profiled-PIRLS-vs-joint-Laplace convention gap, two canonical benchmarks (grouseticks Poisson, probit) converge to a *suboptimal* optimum — a real correctness concern, likely an upstream PIRLS termination-tolerance issue (`parity/glmm-grouseticks-pois.md`, `parity/glmm-probit.md`, `parity/speed-glmm.md`).
4. **lmerTest parity has a dominant hole:** multi-df Satterthwaite F is unimplemented, so `anova()` leaves every multi-level factor/interaction at `p=NA` unless the user switches to KR (`gap/lmerTest.md`).
5. **The CRAN-build / FFI-lifecycle surface is completely unaudited** — no `R CMD check --as-cran`, no offline-vendor-build validation, and an unexamined `panic="abort"` configuration that may turn a Rust panic into an R-session *abort* (auto-reject) rather than a catchable condition (completeness-critic §A1).

**Verdict:** The numeric core and speed thesis are proven. The "clearer errors" thesis is mostly delivered but is *broken on degenerate-data structure*. The package is **not release-ready**: it currently fails its own non-silent-surgery contract in ~9 documented places, has 14 failing tests, and has never been run through the CRAN gate.

---

## 2. Current build/test health (baseline log)

Source: `assessment/build-test-baseline.log` (R CMD INSTALL + full non-stopping test run, 2026-05-31 13:44 EDT).

- **Build:** clean `R CMD INSTALL` **succeeded** — Rust compiled, vignettes built. (No `R CMD check --as-cran` was run; see §7.)
- **Tests:** `pass=2394 fail=14 warn=0 skip=17` across **44 test files**.
- **The 14 failures cluster in 9 files**, and per project memory (MEMORY.md, mote `bd-01KRV31R4BJVQCEF0F58NFD4YN`) **predate the recent mixeff-rs re-pin** — they are tracked, not new regressions:

| File | Failing context | Nature |
|---|---|---|
| `test-glmm.R` | GLMM cases match lme4 within ledger bounds (8 failures) | `cbpp_binomial_logit_profiled_pirls`: theta 0.625 vs 0.642, logLik -277.534 vs -277.502, deviance 555.069 vs 555.004 — exceed the recorded `expected-mismatches.json` bounds |
| `test-bw-lme-tutorial.R` | Bodo Winter core outputs (3) + LRT (1) | optimizer drift exceeding tutorial tolerances |
| `test-emmeans.R` | emmeans reproduces native marginal means | the interaction-mean bug (§5) |
| `test-lme4-parity.R` | random-effect modes where labels align | ranef/label drift |
| `test-parity-scoreboard.R` | scoreboard artifact emission | downstream of parity drift |
| `test-phase2-revive.R` | revive | — |
| `test-phase4.R` | glmm family/link surface vs contract | — |

**Reading:** the build is healthy; the failures are the *known* GLMM-divergence and LMM optimizer-drift cases bumping against ledger tolerances after the re-pin. They are not silent — the parity harness catches them — but several (the GLMM `cbpp` and grouseticks shortfalls) represent genuine under-convergence (§4), not pure convention. The ledger (`inst/extdata/expected-mismatches.json`) needs re-baselining *and* the underlying GLMM optimizer issue triaged before these can honestly be reclassified as `expected_mismatch`.

---

## 3. Capability gap matrix (ranked)

Classification rubric: **in-scope-missing | partial | out-of-scope-by-design (PRD §3) | upstream-blocked | works**.

### 3a. Highest-priority pattern — silent-surgery argument swallowing (cross-cutting, BLOCKER)

Mostly cheap R-side fixes (detect recognized-but-unsupported arg in `...`, refuse with a typed condition).

| Defect | Wrong behavior observed | Source |
|---|---|---|
| `simulate(newparams=)` | ignored; identical to default-param output (breaks power analysis) | `gap/simulate-refit-update.md:40` |
| `simulate(newdata=)` | ignored; returns 180 training rows not 20 | `gap/simulate-refit-update.md:41`, `gap/prediction.md:71` |
| `predict(random.only=TRUE)` | ignored; returns full conditional prediction | `gap/prediction.md:51` |
| `residuals(scaled=TRUE)` | ignored; returns unscaled residuals | `gap/core-extractors.md:67` |
| `logLik/deviance(REML=FALSE)` | ignored; returns fit-criterion value — **numerically wrong** under requested criterion | `gap/model-stats.md:32,36` |
| `glmm(offset=)` | silently dropped (Poisson rate models silently wrong) | `gap/glmm-fit.md:64` |
| non-default factor `contrasts()` | silently ignored; wrong fixef/SE/logLik, no warning; `lmm()` has no `contrasts=` arg | `parity/lmm-contrasts.md` |
| `emmeans(mode=/lmer.df=)` | ignored; asymptotic request silently yields finite df | `gap/emmeans.md:53` |
| `vcov(correlation=TRUE)`, `refit(newweights=)`, `simulate(use.u=)` | ignored silently | `gap/matrix-accessors.md:69`, `gap/simulate-refit-update.md:42,52` |

**Single recommended fix:** a shared `mm_intercept_unsupported_dots()` helper across all S3 methods, raising typed `mm_arg_error`/`mm_inference_unavailable`. Closes 6+ specs at once (`gap/prediction.md:103`, `gap/simulate-refit-update.md:74`).

### 3b. Consolidated family matrix

| Family | Core works? | Biggest in-scope-missing / partial | Out-of-scope (PRD §3) | Upstream-blocked |
|---|---|---|---|---|
| Core extractors (`gap/core-extractors.md`) | yes | `as.data.frame(VarCorr)` wrong shape (breaks broom); `residuals(type=)` only "response"; GLMM resid default diverges | GLMM condVar refusal | Gamma/IG dispersion scale (nlopt) |
| Matrix/structure accessors (`gap/matrix-accessors.md`) | partial | `getME("Zt"/"Lambdat")` **crash** (`Matrix::t()` 1-line fix, `R/revive.R:128,130`); `ngrps()` errors; `weights()` returns NULL; `terms()`/`hatvalues()` absent | `getME("devfun")`; `fortify` | none |
| Model stats (`gap/model-stats.md`) | yes (values) | `logLik/deviance(REML=)` wrong; `extractAIC`/`REMLcrit`/`isREML`/`refitML`/`devcomp` missing (`devcomp` **promised by PRD §6**) | `llikAIC` | none |
| Diagnostics/convergence (`gap/diagnostics.md`) | partial | influence/leverage cluster absent (`influence`,`cooks.distance`,`dfbeta(s)`,`hatvalues`); `rePCA`; `allFit` | R-level optimizer menu; `fortify` | profile-LL CI for GLMM |
| lmerTest surface (`gap/lmerTest.md`) | partial | **multi-df Satterthwaite F unimplemented** (anova p=NA); `anova(type=)` label-only; `summary()` defaults Wald-z not Satterthwaite; `ranova` differs; no `contestMD`/`joint` | `step()`,`get_model()`,`as_lmerModLmerTest` | none (scalar+multi-df KR now work) |
| Inference confint/anova (`gap/inference-confint-anova.md`) | yes (fixed-effect) | bootstrap CIs exclude variance comps; **no `profile()`/`thpr` object** (Phase 5); no `parm` group selectors; `drop1` no `ddf` F-table | KR beyond scalar; profile object (Phase 5) | REML beta-profile refusal (engine contract) |
| emmeans (`gap/emmeans.md`) | yes (broad) | `mode=`/`lmer.df=` swallowed; rank-deficient designs lose SE/df + show fabricated mean; **interaction-mean bug** | `qdrg`,`pbkrtest.limit` | none |
| Prediction/residuals (`gap/prediction.md`) | yes (point) | 3 silent-surgery bugs (§3a); NA-in-newdata hard-errors; `pearson`/`deviance` resid; population `se.fit` | conditional `se.fit`/interval; GLMM predict SE | GLMM SE (nlopt) |
| Bootstrap (`gap/bootstrap.md`) | yes (LMM LRT) | no generic `bootMer(FUN=)`; no VC bootstrap CIs; **`confint(mm_glmm,method="bootstrap")` cryptic base-R error** | `cluster_bootstrap` p-values | GLMM PB-LRT (engine LMM-only) |
| GLMM fitting (`gap/glmm-fit.md`) | yes (certified matrix) | `offset=` dropped; `cbind(succ,fail)` rejected; `weights=` unwired; Gamma default link blocked **despite engine support**; nAGQ>1 != true AGQ | `nAGQ=0`, `glmer.nb`, custom links | `joint_laplace`/true AGQ (nlopt-gated off for CRAN) |
| ranef condVar/plots (`gap/ranef-condvar.md`) | yes (postVar) | `as.data.frame(ranef)` wrong (breaks caterpillar); `dotplot`/`qqmath`/`plot` absent (**not §3-deferred**); `drop`/`whichel` ignored | — | — |
| simulate/refit/update (`gap/simulate-refit-update.md`) | yes (LMM) | silent `newparams`/`newdata` (§3a); no `simulate.mm_glmm`/`refit.mm_glmm`; `refit` rejects `simulate()` df form | `update.lmerModLmerTest` | GLMM simulate (engine) |
| Formula RE-syntax (`gap/formula-re-syntax.md`) | yes (full grammar) | `cs()` compound-symmetric no path (**not §3**); inline `factor()`/`dummy()` rejected; `findbars`/`nobars`/`subbars`/`isNested` absent | `ar1()`/spatial; `mkReTrms` | — |
| Datasets/utils (`gap/datasets-utils.md`) | n/a | `rePCA`; modular `lFormula`/`mkLmerDevfun`; `findbars`/`nobars`/`ngrps`/`isNested` | modular API; `allFit`; `checkConv`; `step` | — |

### 3c. Ranked in-scope gaps blocking "everything lme4 does"

1. **Silent-surgery cluster (§3a).** *Blocker.* Cheap fix, highest stakes.
2. **Multi-df Satterthwaite F unimplemented in `anova()`.** *Blocker.* Dominant lmerTest hole (`gap/lmerTest.md:45,76`). Needs upstream Rust.
3. **emmeans interaction-mean wrong answer + GLMM second-class status.** *Blocker / major* (`parity/inf-emmeans.md`, `gap/glmm-fit.md`, `gap/bootstrap.md:30`).
4. **Broken day-one accessors:** `getME("Zt"/"Lambdat")` crash, `ngrps()` errors, `weights()` NULL, `as.data.frame(VarCorr/ranef)` wrong shape (`gap/matrix-accessors.md:30`, `gap/core-extractors.md:58`, `gap/ranef-condvar.md:31`).
5. **Influence/leverage diagnostic cluster absent** (`gap/diagnostics.md:36-39`).
6. **`anova(type=)` is a label only** — Type II/III never recompute (`gap/lmerTest.md:47`).
7. **Missing model-stats generics** (`extractAIC` breaks base `step`/`drop1`; `devcomp` promised by PRD §6).
8. **Bootstrap breadth** (no `bootMer(FUN=)`, no VC bootstrap CIs).
9. **No `profile()`/`thpr` object** (Phase 5) — note `confint(method="profile")` itself *works* for LMM.
10. **ranef plotting absent (not §3-deferred); cheap residual paths absent.**

### 3d. Genuinely out-of-scope (PRD §3 — do NOT count as defects)

`step()`/`get_model()`/model-selection engine; `cbind(y1,y2)~` multivariate; KR **beyond scalar** (`vcovAdj`, multi-df `KRmodcomp` — note §3's KR bullet is **stale**: scalar KR works and multi-df KR F now works for `anova`, per `gap/lmerTest.md:68`); `ar1()`/spatial covariance; `nlmer`; `I()`/`poly()`/splines/GAM; modular `lFormula`/`mkLmerDevfun`; R-level optimizer menu; **profile-LL CIs for GLMM**. Two **stale-PRD notes** to fix: profile-LL CI for **LMM** is implemented/working, and KR is now in-scope-working.

### 3e. Upstream-blocked (engine / nlopt)

GLMM parametric-bootstrap LRT (engine LMM-only); `joint_laplace`/true AGQ and GLMM SE/profile (nlopt feature-gated off for CRAN — **asserted, not yet verified against the actual `--as-cran` feature set**, §7); REML beta-profile refusal (engine contract); Gamma/IG dispersion scale. **Caveat:** several "GLMM" gaps are **wrapper-side, not upstream** — the crate survey certifies grouped binomial, observation weights, Gamma+inverse, and GLMM parametric bootstrap; those are in-scope-missing wiring.

---

## 4. Empirical parity & speed evidence

30 probes under `assessment/parity/`, each backed by a runnable `*-probe.R`. Tolerances: fixef 1e-4, theta 1e-3, logLik 1e-3, sigma 1e-4.

### LMM parity (15 probes) — within tolerance
All classic topologies pass; **fixef and logLik/AIC/BIC within ~1e-6**. Recurring sub-tolerance pattern (correlated-slope sleepstudy):

| Quantity | Typical max abs diff | Tol | Status |
|---|---|---|---|
| fixef | 1e-12–1e-13 | 1e-4 | within |
| logLik/AIC/BIC | 4e-6–8e-6 | 1e-3 | within |
| theta | 1e-5–1.5e-4 | 1e-3 | within |
| sigma | 1.3e-3–1.4e-3 | 1e-4 | **~14x over** |
| SE (Days) | up to 7.2e-4 | 1e-4 | **~7x over** |
| ranef/fitted | up to 1.1e-2 | 1e-4 | **over** |

Root cause (single hypothesis, attributed convergently): the Rust trust-region REML optimizer stopping at a slightly different theta/sigma on a flat likelihood ridge (logLik agrees to ~4e-6, same basin). Severity **minor**; **vanishes at N=50,000** (`parity/speed-scaling.md`). *Caveat (critic §D3): all probes share one optimizer; thread-count nondeterminism was never excluded — see §7.* The recurrent `vcov[1,1] EXCEEDS-TOL` flags are a **probe-script artifact** (1e-8 applied to raw vcov vs the 1e-4 spec), not a defect.

### Inference surface (10 probes) — mostly works; two bugs, one blocker
- **Works/within-tol:** `anova`/LRT, `drop1` LRT (+Satterthwaite F=t^2, KR available), `confint(method="profile")` (**LMM profile CIs computed** — ML within ~1e-3; honest typed REML-beta refusal `profile_beta_unavailable_under_reml`), parametric bootstrap (CI within ~0.5% width), `simulate` (exact seed reproducibility). Typed refusals clean.
- **`inf-confint-wald` — divergent (major):** SE drift (7.2e-4) propagates to CI bounds; same REML root cause.
- **`inf-getME` — two confirmed bugs:** `getME("Zt")`/`getME("Lambdat")` **crash** (`"argument is not a matrix"`, bare `t()` on a sparse `Matrix`; fix `Matrix::t()` at `R/revive.R:128,130`); core components (X, Z, theta, beta, Lambda) match exactly.
- **`inf-emmeans` — BLOCKER:** interaction-model marginal means evaluated at reference level, not averaged (temperature=215: 19.3 vs lme4 35.8; contrast sign reversals). Native `mm_means`/`mm_comparisons` show the same wrong values -> root is reference-grid population in `recover_data.mm_lmm` (**suspected** missing `xlev`, unconfirmed — critic §D5). `reliability="low"` does not signal wrongness. **Silent wrong answer for a core advertised feature.**
- **`inf-predict-newdata` — major bug:** population predictions exact (5.7e-13); `allow.new.levels=TRUE` with canonical NA-response newdata **fails** with a misleading `mm_data_error` (response column validated before stripping, `R/predict.R:~188`). Dummy-response workaround reproduces lme4 to 1e-12.
- **`lmm-contrasts` — silent wrong answer (major):** non-default `contrasts()` ignored, no warning.
- **`lmm-offset` — in-scope-missing (major):** typed `mm_formula_error` refusal (honest) but no offset pathway.

### GLMM parity (6 probes) — divergence by design, plus genuine under-convergence
GLMMs use **profiled fast-PIRLS** (`docs/glmm_support_contract.md`, pre-registered in `inst/extdata/expected-mismatches.json`) — not lme4's joint Laplace, so coefficient/logLik gaps are expected and classed `documented_divergence`.
- **Documented convention (expected):** `cbind(success,failure)~` refused (Phase-4 deferred); Bernoulli-expanded workaround diverges within registered bounds; AGQ `nAGQ>1` is profiled-PIRLS-with-AGQ-deviance, not joint AGQ; `glmm-gamma` RE variance ~10x smaller than lme4 (0.023 vs 0.246) **by design but emitted silently** (transparency gap, major).
- **Genuine under-convergence (HIGH PRIORITY, beyond the documented gap):** `glmm-grouseticks-pois` lands **0.524 logLik short** (~500x tol), SEs inflated ~4x; `glmm-probit` fixef ~5% low, SE **underestimated ~11–21%** (anti-conservative); `speed-glmm` (N=10k, 5 seeds) lme4 *always* achieves higher logLik, fixef gap ~3.5–4.3e-3 while theta is identical -> **fixed-effect subproblem under-tightened**, and `mm_control()` exposes no tolerance knob. **A minimal reproducer is a named, deferred deliverable** (CLAUDE.md prerequisite before filing the upstream mote issue).

### Speed — promise confirmed
| Operation | lme4 | mixeff | Speedup |
|---|---|---|---|
| `speed-scaling` N=1k/10k/50k | 0.032/0.146/0.783 s | 0.004/0.014/0.062 s | **8x->11x->13x** |
| GLMM binomial N=10k | 0.188 s | 0.055 s | 3.4x |
| Parametric bootstrap nsim=499 | 1.907 s | 0.051 s | **37x** |
| ML profile CI | 1.003 s | 0.337 s | 3.0x |
| `anova`/LRT fit+compare | 0.034 s | 0.007 s | 4.8x |
| `simulate` | 1.9 ms/call | 0.25 ms/call | 7.6x |

The ratio *improves* with size. **Honest qualification (critic §D1):** for GLMM, "faster" must be stated as *"faster where it converges to the same optimum, which for GLMM it sometimes does not"* — `speed-glmm` is faster *and* converging to a worse optimum. A faster wrong answer is not parity.

---

## 5. Error-message quality — promise kept / broken

21 probes (`assessment/errors/`). Tally: **8 clearer-than-lme4, 8 needs-work, 4 bug/worse-than-lme4**, plus scoped-by-design refusals.

### Promise delivered (clearer than lme4)
R-side input validation firing *before compute*, typed catchable `mm_*` conditions where lme4 leaks base-R/C++/parser fragments or silently drops data: missing/misspelled variables (`misspelled-var.md`, `response-missing.md`), NA in design / all-NA response (`na-in-predictor.md`, `all-na-response.md` — lme4 *silently drops*), mismatched weights length (`mismatched-weights.md`), malformed formulas (`malformed-formula.md`), unsupported family/link (`unsupported-family.md` — typed refusal + `$supported` table vs lme4's silent-fit / C++ NaN), predict() guards (`newdata-missing-col.md`, `new-levels-predict.md`), nonexistent-coef lookup (`contrast-nonexistent-coef.md` — lme4 returns silent `NA`), `nAGQ>1` topology (`nagq-on-multi-re.md`), convergence/singular diagnostics (`convergence-hard.md`, `singular-fit.md`, `too-few-groups.md` — typed diagnostics vs lme4's uncatchable `message()`), revive (`revive-no-handle.md`).

### Promise broken (bugs — must fix)
1. **Silent wrong answer — empty data to `glmm()`** (`empty-data.md`): zero-row frame returns a structurally valid fit with fabricated `fixef=0`, `logLik=0`. Fix: `nrow==0` guard -> `mm_data_error`.
2. **Silent wrong answer — single-obs-per-group `lmm()`** (`single-obs-per-group.md`): lme4 hard-errors; mixeff fits silently with `is_singular=FALSE`/`boundary=FALSE` (both misleading). The `mm_not_identifiable` class exists but is unused.
3. **Rust-panic leak — empty data to `lmm()`** (`empty-data.md`): `mm_bridge_error: Matrix index out of bounds`. **Re-severitize:** with `panic="abort"` (`tools/config.R:74`) this may **abort the R session**, not throw — an automatic CRAN reject (critic §A1/§D2). Probe before classifying as mere phrasing.
4. **Structural misdiagnosis — n=1 row** (`one-row.md`): mixeff reports "Constant response" when the true fault is structural (n=1). Worse than lme4's accurate pre-fit guard. Fix: pre-fit `nrow>=2` + >=2-level guards in `R/fit-lmm.R`/`R/glmm.R`.

### Worse text / needs-work
Factor/character response leaks `Response 'y' not found or not numeric` ("not found" wrongly implies a typo), wrong condition class, audit block prints before error (`factor-response-lmm.md`, `character-response.md`); Rust-enum leak in fix hints (`Use NewReLevels::Population`, `new-levels-predict.md`); `doTryCatch()`/raw-`mm_*:`-tag leakage in chained conditions (`malformed-formula.md`, `no-random-effect.md`); ambiguous grouping-not-factor message; no REPL signal for reduced-rank/boundary fits; wrong-column zeroing under collinearity; GLMM `boundary_parameter` index not term name; reserved-arg wall for `glmm()`; duplicated diagnostic rows; GLMM `contrast()` untyped dispatch error.

**Root cause across failures:** absence of **pre-fit R-side structural guards** (`nrow>=2`, >=2 levels/grouping factor, response numeric, n_groups<n_obs) — exactly the audit-first checks PRD §8.1 mandates — plus systemic **bridge-condition formatting** leakage.

---

## 6. The assurance backlog (P0 -> P2)

230 specs across 15 families (`assessment/testspec/`) vs 44 existing test files (`assessment/survey/tests-{0..7}.md`). P0 = release blocker; P1 = typical user hits it; P2 = hardening.

### P0 — release blockers
- **P0.1** Silent-wrong-answer in `predict`/`residuals`/`simulate` — **zero coverage today**. One shared `mm_intercept_unsupported_dots()` helper closes 6+ specs.
- **P0.2** `getME("Zt"/"Lambdat")` crash — `Matrix::t()` fix at `R/revive.R:128,130` + regression + revival tests.
- **P0.3** `confint(mm_glmm)` cryptic `"non-numeric argument to binary operator"` — R-side typed `mm_inference_unavailable` stub; register `confint.mm_glmm` in NAMESPACE.
- **P0.4** Absent/broken accessors breaking downstream tooling: `as.data.frame(VarCorr)`, `as.data.frame(ranef, condVar=TRUE)` (discards working Rust `postVar` — no-silent-surgery violation), `ngrps()`, `weights()`, `extractAIC.mm_lmm`, `residuals(type=...)` `match.arg` widening.
- **P0.5** Silent-wrong-value model stats: `logLik(REML=FALSE)` returns REML value; `deviance()` on REML fit returns REMLcrit; `predict(newdata, NA-response)` NaN guard fires before new-levels policy.
- **P0.6** `anova(type=)` no-op label; multi-df Satterthwaite returns `NA`. **Both need upstream Rust FFI changes.**
- **P0.7** `glmm(offset=)` silently dropped — one-line R fix.
- **P0.8** `simulate.mm_glmm`/`refit.mm_glmm` absent (PRD §10 Phase 4); `refit()` rejects single-column df. SRU-06/07 upstream-blocked; SRU-04 R-side.

### P1 — high value (parity credibility core)
- **No ML (`REML=FALSE`) LMM parity fixture** — all 7 classic fixtures are REML (`tests-0.md` gap #1).
- **`vcov(fit)` never compared to `vcov(ref)`** — only self-consistency checked (`tests-0.md` gap #2).
- **No numeric lme4/lmerTest parity on inference quantities** (p, t/F, df checked only for finiteness — `tests-2.md` G11).
- **GLMM parity is 2 cases** — no Gamma/probit/cloglog/sqrt, no fitted/resid/VarCorr/ranef/AIC parity. **Re-baseline the cbpp ledger entry currently failing (§2).**
- **`ranef()` BLUP parity absent from all 8 tutorial files.**
- **Boundary/singular discipline:** `is_singular()` on a non-singular fit never asserted FALSE; dyestuff2 surfacing untested.
- **lme4-name typed redirects** for `isSingular`/`allFit`/`bootMer`/`influence`/`hatvalues`/`cooks.distance`.
- **Cross-session revival** for every new R method to guard the JSON-source-of-truth contract (PRD §8).
- **emmeans interaction bug regression** + `mode=`/`lmer.df=` swallow tests.
- **Fit control:** `contrasts=`, `subset=`, `na.action`, `mm_control(ftol_abs/maxeval)`, `update.mm_lmm()`.

### P2 — hardening
Argument validation (`bootstrap_control(nsim<=0)`, `confint(level=1.5)`, `confint(parm="nonexistent")`); plotting dispatch; `ranef(drop=/whichel=)` swallow; GLMM summary across families; Phase-5 typed-refusal sentinels (`profile()`, `confint.thpr`, `PBrefdist`, `bootMer`); VarCorr `" & "` vs `":"` naming; schema/manifest coverage.

**Caveat:** several specs intentionally pin *current* divergences as documented locks (e.g. `ranef-condvar.md` TS-RC-10 postVar at 5e-2; cake-interaction VarCorr naming). These are test-gaps to *lock*, not bugs to fix.

**Highest-leverage:** (1) `mm_intercept_unsupported_dots()` helper; (2) `Matrix::t()` fix + tests; (3) widen `residuals` `match.arg` + add `as.data.frame.mm_*` + `ngrps`/`weights`/`extractAIC`/`terms`; (4) ML-fit + `vcov`-vs-lme4 parity assertions; (5) coordinate the 3 upstream Rust items via mote.

---

## 7. CRAN-readiness checklist

**The single largest gap in the assessment itself: nearly the entire CRAN-build / FFI-lifecycle surface is unexamined** (completeness-critic §A). mixeff is an R-over-Rust FFI package; the audit largely treated it as pure-R.

| Gate | Status | Action |
|---|---|---|
| `R CMD check --as-cran` clean | **NOT RUN** | Run on a clean offline build; inventory every NOTE/WARNING/ERROR. Use `cran-prepare`/`r-cmd-check`. **Gates everything.** |
| Tests pass | **14 failing** (§2) | Re-baseline GLMM ledger + triage under-convergence; resolve emmeans/contrasts/tutorial failures. |
| `panic="abort"` => no reachable input aborts R | **UNVERIFIED** (`tools/config.R:74`) | Probe empty/n=1/all-NA: does R *abort* or *throw*? Abort = auto-reject. Re-severitize `empty-data.md`. |
| Offline vendor build validates | **UNVERIFIED** | Run `tools/check-offline-install.R`, `check-vendor-drift.R`, `check-license-note.R`; CRAN builds offline. |
| Actual CRAN feature set verified | **ASSUMED** | Resolve `@CRAN_FLAGS@`; confirm nlopt gating; re-run GLMM-divergence probes against the `--as-cran`-flagged build. |
| `SystemRequirements: Cargo>=1.78.0` graceful degradation | **UNVERIFIED** | Probe build with Rust absent/old. |
| Cross-platform (Win/macOS-ARM/Linux) | **UNVERIFIED** | `Makevars.win.in`, `mixeff-win.def`, `@TARGET@`/`@PROFILE@`; `.github/` CI new/untracked. |
| `inst/LICENSE.note` covers bundled-Rust copyright | **UNREVIEWED** | CRAN requires bundled-source disclosure. |
| Rd examples run / `\dontrun`->`\donttest` | **NOT RUN** | Only 3/9 Rd files run unguarded; 6 use `\dontrun` (likely a NOTE). |
| spell-check / URL-check / `checkRd` | **NOT RUN** | `devtools::spell_check`, URL check. |
| NAMESPACE <-> `man/` reconciliation | **NOT DONE** | Check exported-but-undocumented / documented-but-unexported drift. |
| Determinism / threading | **UNVERIFIED** | Fit one model 50x and `RAYON_NUM_THREADS=1` vs default — is theta bit-identical? Settles flat-ridge-vs-nondeterminism (§4). |
| R-seed bridging | **UNVERIFIED** | Does `set.seed(1); simulate(fit)` reproduce from R's global `.Random.seed` and restore it? Seed flows to a *Rust* RNG (`R/inference.R`). |

---

## 8. Top 10 prioritized next actions

1. **Run `R CMD check --as-cran` on a clean offline build** and inventory every NOTE/WARNING/ERROR. Highest priority — gates everything (§7).
2. **Probe `panic="abort"`:** enumerate empty/n=1/all-NA inputs; record abort-vs-throw. Re-severitize the error-quality "bugs" (§5, §7).
3. **Eliminate the silent-surgery cluster** with the shared `mm_intercept_unsupported_dots()` helper across all S3 methods + regression tests (§3a, P0.1).
4. **Fix the emmeans interaction-mean blocker** — confirm the suspected `recover_data.mm_lmm`/`xlev` root cause, patch, add regression (§4, P1).
5. **Triage GLMM under-convergence**: isolate a minimal grouseticks/probit reproducer (CLAUDE.md prerequisite), file the upstream `mixeff-rs` mote issue, re-baseline the failing `cbpp` ledger entry (§2, §4).
6. **Two-line `Matrix::t()` fix** at `R/revive.R:128,130` + `getME("Zt"/"Lambdat")` tests across revival (§3b, P0.2).
7. **Widen `residuals` `match.arg`; add `as.data.frame.mm_varcorr`/`as.data.frame.mm_ranef`, `ngrps`, `weights`, `extractAIC`, `terms`** — pure-R, restores broom/merTools/`step` interop (P0.4).
8. **Add an ML-fit LMM parity case and a `vcov`-vs-lme4 parity assertion** — closes the two structural parity holes most undermining the PRD §11 claim (P1).
9. **Coordinate the three upstream Rust items via mote** — multi-df Satterthwaite F, `anova(type=)` SS recompute, GLMM simulate surface (gate P0.6/P0.8).
10. **Validate offline vendor build + verify the actual CRAN feature set, then re-run GLMM probes against that build**; settle determinism (`RAYON_NUM_THREADS`) and R-seed bridging; run Rd examples / spell / URL / NAMESPACE-vs-`man` reconciliation (§7).

---

### Resolved internal contradictions (for the record)
- **LMM profiling:** `confint(method="profile")` *is implemented and works* for LMM (ML CIs within ~1e-3, typed REML-beta refusal). What is *absent* is the `profile()` generic / reusable `thpr` object (PRD §10 Phase 5). Both source sections are correct about different objects.
- **GLMM speed verdict:** faster, but honestly only "faster where it converges to the same optimum"; GLMM sometimes does not (§4).
- **Stale PRD §3:** KR is now in-scope-working (scalar + multi-df anova F); only KR *beyond scalar* (`vcovAdj`) remains deferred.

*Evidence: `assessment/build-test-baseline.log`, `assessment/gap/*.md` (17), `assessment/parity/*.md` (30), `assessment/errors/*.md` (21), `assessment/testspec/*.md` (15), `assessment/survey/tests-{0..7}.md` (8), `inst/extdata/expected-mismatches.json`, `tools/config.R`, `R/revive.R`, `R/predict.R`, `R/inference.R`, `R/glmm.R`.*
