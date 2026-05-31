# Gap Report — emmeans / Marginal Means Family

**Family:** emmeans / marginal means (recover_data, emm_basis, emmeans, contrast,
pairs, joint_tests, ref_grid, lsmeans, emtrends, eff_size, df/adjust handling).

**Reference:** `assessment/survey/lme4-emmeans.md` (emmeans 2.0.3, lme4 2.0.1,
lmerTest 3.2.1, all installed and exercised live).

**mixeff source:** `R/emmeans.R` (the bridge), `R/marginal.R` (native
`mm_grid`/`mm_means`/`mm_comparisons`/`mm_predictions`), `R/inference.R`
(`contrast`, `df_for_contrast`, `test_effect`, `test_random_effect`), `R/zzz.R`
(`.emm_register` hook).

**Method of assessment:** Ran live `library(mixeff); library(emmeans); library(lme4)`
comparisons on `sleepstudy` (LMM) and an expanded-Bernoulli `cbpp` (GLMM). All
status calls below are backed by executed output, not source reading alone.

---

## Summary verdict

The emmeans bridge is real and broad: `emmeans()`, `ref_grid()`, `lsmeans()`,
`contrast()` (all built-in `.emmc` methods), `pairs()`, `joint_tests()`,
`emtrends()`, `eff_size()`, `summary()`/`confint()`/`test()` with `infer=` and
`adjust=`, `at=`, `weights=`, and GLMM `type="response"` back-transforms all work
and match lme4 numerically (e.g., `joint_tests` F-ratios agree to 3 decimals;
GLMM odds ratios and probabilities back-transform correctly). The native
lmerTest-style surface (`mm_means`, `mm_comparisons`, `test_effect`,
`test_random_effect`, `df_for_contrast`) provides contract-preserving
equivalents of `ls_means`/`difflsmeans`/`contest*`/`ranova`/`calcSatterth`.

Two concrete defects a real lme4 user would hit:

1. **`mode=` / `lmer.df=` are silently ignored** by `emm_basis.mm_lmm`. The
   bridge uses its own `method=` argument; the emmeans-standard `mode=` /
   `lmer.df=` arguments never reach it, so `emmeans(fit, ~x, mode="asymptotic")`
   silently returns finite Satterthwaite-ish df instead of `Inf`. This is both a
   parity gap and a no-silent-surgery violation (PRD/CLAUDE.md).

2. **Rank-deficient (missing-cell) designs lose SEs/df for *all* cells**, and a
   non-estimable cell is shown with a fabricated numeric mean rather than
   `nonEst`. lme4 marks only the offending cell `nonEst` and keeps valid
   estimates+SEs+df for every estimable cell.

---

## Capability table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `recover_data.merMod` (plumbing) | `recover_data.mm_lmm` / `.mm_glmm` present + registered via `.emm_register` on emmeans load (`R/zzz.R:49-52`, `R/emmeans.R:27,101`) | works | — | `emmeans(fit, ~grp)` builds grid without user action; `ref_grid(fit)` returns proper grid. |
| `emm_basis.merMod` — X, bhat, V, dffun | Implemented (`R/emmeans.R:47-97`); V from `mixedmodels.fixed_effect_covariance_matrix`; per-contrast dffun via `df_for_contrast` | works | — | `emmeans(fit,~grp)` SE=9.63, df=21.8 match lme4 KR/Satt to 3 sig figs. |
| `mode=` / `lmer.df=` arg (`"kenward-roger"`/`"satterthwaite"`/`"asymptotic"`) | **Silently ignored.** Bridge exposes `method=` instead; emmeans passes `mode=`, which is swallowed by `...`. | partial | **major** | `summary(emmeans(fit,~grp,mode="asymptotic"))$df` → 21.77 (should be Inf). `mode="kenward-roger"`, `mode="satterthwaite"`, `lmer.df="asymptotic"` ALL → 21.77. By contrast `method="asymptotic"`→Inf, `method="satterthwaite"`→21.77. lme4: KR 21.767, Satt 21.766, asymptotic Inf. Violates no-silent-surgery. |
| Mode name `kenward-roger` (hyphen) vs `kenward_roger` (underscore) | Underscore is the mixeff `method=` token; hyphenated emmeans `mode=` token is ignored entirely (see above) | partial | major | Same repro as above — neither hyphen nor underscore `mode=` changes df. |
| Automatic df fallback chain (KR→Satt→asymptotic) | mixeff `method="auto"` resolves internally (KR/Satt available; auto≈Satt here) but is not driven by emmeans `mode=` | partial | minor | `df_for_contrast(fit,L,method="auto")`=159.2 equals `"satterthwaite"`; `"kenward_roger"`=159.2; `"asymptotic"`=NA (not Inf). |
| `pbkrtest.limit` / `lmerTest.limit` size guards | No emmeans-side size guard; df always computed by the Rust engine regardless of N | out-of-scope-by-design | minor | The limits exist in lme4 because KR/Satt are expensive in R; mixeff computes df in Rust. No analogous guard, and none needed. Document the difference. |
| `ref_grid()` (incl. `at=`, `cov.reduce`) | works | works | — | `ref_grid(fit)` reduces Days to mean 4.5; `emmeans(fit,~Days,at=list(Days=c(0,5,9)))` honored. |
| `emmeans()` core specs (`~a`, `~a*b`, `~a|b`, `pairwise~a`) | works | works | — | `emmeans(fit,~grp)` matches lme4 means; `~Days` averages over grp. |
| `emmeans(weights=)` (`equal`/`proportional`/…) | Accepted and passed through (handled by emmeans on the linfct) | works | — | `emmeans(fit,~grp,weights="proportional")` runs; equal cell sizes make it equal here. |
| `lsmeans()` alias | works | works | — | `lsmeans(fit,~grp)` returns lsmean column, identical numbers. |
| `contrast()` — all built-in `.emmc` methods | works (pairwise, dunnett, poly, eff, etc. via emmeans machinery) | works | — | `contrast(em,"dunnett")`, `contrast(em,"poly")` both produce correct contrasts with tukey/dunnettx adjustment notes. |
| `pairs()` / `reverse=` | works | works | — | `pairs(emmeans(fit,~grp))` → 3 tukey-adjusted comparisons, df=159. |
| `joint_tests()` | works, numerically matches lme4 | works | — | mixeff: Days F=167.668 df2=159; grp F=0.172 df2=159.21. lme4: Days F=167.669; grp F=0.172. |
| `emtrends()` | works | works | — | `emtrends(fit,~grp,var="Days")` → trend 10.5, SE 0.808, df 159. |
| `eff_size()` | works (requires `sigma()` + `df.residual()`, both implemented) | works | — | `eff_size(em,sigma=sigma(fit),edf=df.residual(fit))` returns effect sizes; `sigma(fit)`=31.15. |
| `summary.emmGrid` `infer=`, `adjust=`, `level=` | works | works | — | `summary(em,infer=c(TRUE,TRUE))` adds t.ratio/p.value; tukey/dunnettx adjustments applied. |
| `confint.emmGrid` / `test.emmGrid` | works (inherited from emmeans on the grid) | works | — | CIs present in all summaries; `test()`/`joint=` route through emmeans. |
| `update`/`predict`/`plot`/`vcov`/`coef`/`as.data.frame`/`as.list`/`rbind` emmGrid methods | works (emmeans-owned, operate on the populated grid) | works | — | `as.data.frame(pairs(em))` used in tests; these are emmeans methods on a valid emmGrid. |
| GLMM `emm_basis.mm_glmm` — asymptotic z, link labels | works | works | — | `emmeans(fg,~period)` df=Inf on logit scale; init label "Results are given on the logit scale". V from `pirls_laplace_working_hessian`. |
| GLMM `type="response"` back-transform | works | works | — | `summary(emmeans(fg,~period,type="response"))` → prob column with back-transformed CIs; `pairs(...,type="response")` → odds.ratio with log-odds-ratio note. |
| Non-estimable cells via `estimability::nonest.basis()` (→ `nonEst`) | mixeff hardcodes `nbasis = estimability::all.estble` (`R/emmeans.R:91,162`); never flags `nonEst`. Engine drops aliased columns at fit time instead. | partial | **major** | Missing-cell factorial `A*B` (a3 never b2): lme4 marks `a3 b2` `nonEst`, gives valid estimates/SE/df for the other 5 cells. mixeff: V reported unavailable (`rank_deficient_fixed_effects`), so ALL 6 cells get `SE=NA, df=Inf, CL=NA`, and the non-estimable `a3 b2` cell shows a fabricated numeric mean (282). Two problems: (a) estimable cells lose inference; (b) numeric mean shown for a non-estimable cell. |
| Clear diagnostics for rank-deficient designs | mixeff emits structured reasons (`fixed_effect_rank_deficient`, `fixed_effect_empty_cell`) — clearer than lme4's terse "dropping 1 column" | works | — | mixeff prints which cell (A=a3,B=b2) and rank (5 of 6). This is the "clearer errors" win. |
| `lmerTest::ls_means()` / `lsmeansLT()` | Native equivalent: `mm_means(fit, ~factor)` (contract-preserving, status/reason fields) + emmeans `emmeans()` | works | — | `mm_means` exported (`NAMESPACE:152`); also `emmeans(fit,~grp)` covers it. |
| `lmerTest::difflsmeans()` | Native equivalent: `mm_comparisons()` / `pairs(emmeans(...))` | works | — | `mm_comparisons` exported (`NAMESPACE:145`); pairwise differences with Satt df. |
| `lmerTest::contest()` / `contest1D()` / `contestMD()` | Native equivalent: `contrast(fit, L)` (single/multi-row L) + `test_effect(fit, term)` for term-level F | works | — | `contrast`/`test_effect` exported (`NAMESPACE:129,169`); produce estimate/SE/df/t and F-tests. |
| `lmerTest::calcSatterth()` | Native equivalent: `df_for_contrast(fit, L, method="satterthwaite")` | works | — | exported (`NAMESPACE:130`); returns 159.2 for grp contrast. |
| `lmerTest::ranova()` | Native equivalent: `test_random_effect(fit, term)` (boundary-aware LRT) | works | — | `test_random_effect` exported (`NAMESPACE:170`); LRT route for dropping RE terms. |
| `lmerTest::show_tests()` | No equivalent (design-matrix-of-test diagnostic) | in-scope-missing | minor | Niche diagnostic; lme4 users rarely call it directly. Not in PRD non-goals. |
| `lmerTest::step()` | No equivalent — automated backward elimination | out-of-scope-by-design | minor | PRD §3 non-goals: no model-selection/recommendation engine (`recommend_model`/`auto_random_effects`). `step()` is exactly that. |
| `qdrg()` | N/A (emmeans-side helper for non-model objects); mixeff fits are real model objects with a bridge | out-of-scope-by-design | cosmetic | Not needed — `emm_basis.mm_lmm` is the proper hook. |
| `add_grouping()`, `add_submodels()`, `pwpm()`, `pwpp()`, `emmip()`, `as.glht()`, `regrid()` | emmeans-owned methods that operate on a populated `emmGrid`; inherited for free once the bridge produces a valid grid | works | — | These dispatch on `emmGrid`, not on `mm_lmm`; they work on any grid the bridge builds (e.g., `regrid` used implicitly by `type="response"`). Not separately re-tested but architecturally inherited. |
| `emm_options()` / `get_emm_option()` (`lmer.df`, limits) | emmeans-global; `lmer.df` default has no effect because `mode=` is ignored (see row above) | partial | minor | Setting `emm_options(lmer.df="asymptotic")` would not change mixeff df — same root cause as the `mode=` defect. |

---

## Detailed findings

### Finding 1 — `mode=` / `lmer.df=` silently ignored (major)

`emm_basis.mm_lmm` (`R/emmeans.R:47-51`) declares its df selector as `method =
c("auto","satterthwaite","kenward_roger","asymptotic","none")`. emmeans, however,
passes the user's df request as `mode=` (and the deprecated alias `lmer.df=`),
not `method=`. Those arguments fall into `...` and are dropped. Net effect: the
documented emmeans df controls are inert.

Repro (live):
```
mode="asymptotic"     -> df 21.77   (lme4: Inf)
mode="kenward-roger"  -> df 21.77
mode="satterthwaite"  -> df 21.77
lmer.df="asymptotic"  -> df 21.77
method="asymptotic"   -> df Inf      (mixeff's own arg DOES work)
method="satterthwaite"-> df 21.77
```
This is a parity defect (lme4 users will type `mode=`) and a no-silent-surgery
violation: a requested asymptotic z-test silently becomes a finite-df t-test.
Fix is small: accept `mode`/`lmer.df` in the signature, normalize hyphen→
underscore, and map to the internal `method`. The current tests
(`tests/testthat/test-emmeans.R:24,55`) pass `method=` directly, which is why
this was never caught — see Finding 3.

### Finding 2 — rank-deficient designs (major)

`emm_basis.mm_lmm`/`.mm_glmm` always set `nbasis = estimability::all.estble`
(`R/emmeans.R:91,162`), asserting every linear function is estimable. For a
full-rank fit this is fine. For a missing-cell design the engine reports the
fixed-effect covariance as unavailable (`mm_status="unavailable",
reason="rank_deficient_fixed_effects"`), so `V` carries no usable numbers and
emmeans renders `SE=NA, df=Inf, CL=NA` for every cell — including the estimable
ones. Simultaneously, because `nbasis` claims full estimability, the genuinely
non-estimable cell is shown with a fabricated point estimate rather than `nonEst`.

lme4 behavior on the same data: only the offending cell is `nonEst`; the five
estimable cells get correct estimates, SEs, and KR df. mixeff's structured
diagnostics (which cell, what rank) are genuinely clearer, but the loss of
inference on estimable cells and the fabricated mean on the non-estimable cell
are regressions against both lme4 and the package's own "no fabrication"
contract. Correct fix: derive `nbasis` from the real rank/null space (an
`estimability::nonest.basis`-equivalent built from the engine's design rank) so
non-estimable rows are flagged `nonEst` while estimable rows keep their SE/df.

### Finding 3 — test gap on `mode=` (test-gap, major root)

`tests/testthat/test-emmeans.R` and `test-emmeans-glmm.R` exercise the bridge
only with `method="asymptotic"` passed directly. They never call emmeans with
the user-facing `mode=`/`lmer.df=`, and never test a rank-deficient grid. That
is why Findings 1 and 2 are latent. Recommended: add tests asserting (a)
`mode="asymptotic"` yields `df=Inf`, `mode="kenward-roger"`/`"satterthwaite"`
yield finite df matching `df_for_contrast`, and (b) a missing-cell factorial
flags `nonEst` on the empty cell while preserving SE/df on estimable cells.

### What is solidly at or above parity

- `joint_tests` matches lme4 F-ratios to 3 decimals (Days 167.668 vs 167.669).
- GLMM `type="response"` produces correct probabilities and odds ratios with
  proper scale annotations and `df=Inf`.
- `emtrends`, `eff_size`, all built-in `contrast` methods, `at=`, `weights=`,
  `lsmeans`, `ref_grid` all work and read identically to lme4 output.
- The native surface (`mm_means`/`mm_comparisons`/`contrast`/`test_effect`/
  `df_for_contrast`/`test_random_effect`) covers the lmerTest functions with the
  added benefit of row-level status/reason fields the lme4 stack lacks.

### Out-of-scope (PRD §3) — not defects

- `lmerTest::step()` = automated model selection → PRD §3 non-goals
  (no recommendation/selection engine).
- `qdrg()` = grid-from-raw-components helper → unnecessary; mixeff has a proper
  `emm_basis` hook.
- `pbkrtest.limit`/`lmerTest.limit` = R-side cost guards → not applicable; df is
  computed in Rust.
