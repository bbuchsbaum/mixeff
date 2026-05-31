# Gap Analysis — "lmerTest surface" capability family

**Date:** 2026-05-31
**Reference:** `assessment/survey/lme4-lmerTest.md` (lmerTest 3.2.1, lme4 2.0.1)
**mixeff source:** `R/inference.R`, `R/compare.R`, `R/marginal.R`, `R/methods-summary.R`, `NAMESPACE`
**Verification:** Live `Rscript` comparisons against `lmerTest::lmer()` on `sleepstudy`, `cake`, `ChickWeight`.

---

## Orientation

mixeff does **not** re-expose lmerTest's verb names. The lmerTest surface is
re-expressed as a different, audit-first API:

| lmerTest verb | mixeff analog |
|---------------|---------------|
| `summary(ddf=)` coef table | `summary(mm_lmm, method=)` / `inference_table()` |
| `anova(type=, ddf=)` | `anova(mm_lmm, type=, method=)` |
| `drop1()` (Satterthwaite F) | `drop1(mm_lmm)` (LRT) + `test_effect()` |
| `ranova()` / `rand()` | `test_random_effect()` (boundary LRT) |
| `contest()` / `contestMD()` / `contest1D()` | `contrast()` + `test_effect()` |
| `ls_means()` / `lsmeansLT()` | `mm_means()` |
| `difflsmeans()` | `mm_comparisons()` |
| `calcSatterth()` | `df_for_contrast()` |
| `show_tests()` | internal `mm_show_tests` (not exported) |
| `step()`, `get_model()` | (none — out of scope §3) |
| `as_lmerModLmerTest()` | (none — no lme4 interop, by design §3) |

Numerically, where computed, mixeff Satterthwaite/KR results match lmerTest to
the documented tolerances (verified `sleepstudy` Days df 16.98 vs 17.00; `cake`
ls_means estimates/SE/df identical; KR multi-df F for `cake`/`ChickWeight`
matches lmerTest exactly).

---

## Capability table

| lme4/lmerTest capability | mixeff status | classification | severity | Evidence / repro |
|---|---|---|---|---|
| `summary(ddf="Satterthwaite")` coef table with `df` + `Pr(>\|t\|)` | partial | partial | major | `summary(mm_lmm)` **defaults to asymptotic Wald-z with `df=NA`**, not Satterthwaite. lmerTest defaults to Satterthwaite. User must call `summary(m, method="satterthwaite")` to get df+t p-values, which then match (df 16.998 vs lmerTest 17.000). A lme4 user expecting `summary()` to show df/p will not get them by default. |
| `summary(ddf="Kenward-Roger")` | works | works | — | `summary(m, method="kenward_roger")` resolves KR. KR scalar df matches (df 17). |
| `summary(ddf="lme4")` (no p-values) | works | works | — | `summary(m, tests="none")` / `method="none"` gives the no-p-value path. |
| `summary` coef columns `Estimate/Std.Error/df/t/Pr` | partial | partial | minor | Columns present but renamed/transposed into an audit table (`method`, `status`, `reliability` columns added; statistic labeled `t value`/`z value`). Not column-compatible with `coef(summary())` consumers. |
| `anova(type="III", ddf="Satterthwaite")` single-df term | works | works | — | `anova(m, method="satterthwaite")` on `sleepstudy` Days: p=3.30e-06 vs lmerTest 3.264e-06 (within tol). Reported as `t` statistic, not `F` (t²=F). |
| `anova(...)` **multi-df** term under Satterthwaite (factor >2 levels, interaction) | **in-scope-missing** | **in-scope-missing** | **blocker** | `anova(cake_model, method="satterthwaite")` returns `status="unsupported"`, `p_value=NA`, `reason="multi-df fixed-effect contrast tests are not implemented in this scaffold"` for `recipe` (2 df) and `recipe:temp` (2 df). Same for `ChickWeight` `Diet` (3 df). lmerTest gives F for all. **Any model with a multi-level factor or interaction gets no Satterthwaite p-value for those terms — the lmerTest default workflow.** |
| `anova(...)` multi-df term under Kenward-Roger | works | works | — | `anova(cake, method="kenward_roger")`: recipe F=0.0957/2df DenDF=254.0, recipe:temp F=0.0417/2df — matches lmerTest exactly. `ChickWeight` Diet F=6.27/3df DenDF=45.5 vs lmerTest 6.275/46.03. So the multi-df machinery exists for KR but not Satterthwaite. |
| `anova(type=)` actually changes the SS hypothesis (I/II/III) | partial | partial | major | `type` is a **label only**: `anova.mm_lmm` sets `table$type <- type` (R/compare.R:300) and never forwards `type` to the Rust term builder (`mm_rust_term_table` takes no type arg). Type II and Type III return identical statistics in all tested cases (`ChickWeight` Diet III vs II byte-identical: F=1.71→unsupported, Time identical). lmerTest derives genuinely different contrast matrices per type. mixeff is effectively always one fixed parameterization. |
| `anova(m1, m2, ...)` model-comparison LRT | works | works | — | `anova(mm_lmm, other)` dispatches to `compare()` (R/compare.R:274). Provides LRT for nested models. |
| `drop1(ddf="Satterthwaite")` single-term F-table | partial | partial | major | `drop1(mm_lmm)` is a **refit-based LRT** (`statistic_name`/`LRT`, chi-square p), not lmerTest's Satterthwaite single-term-deletion F. Different statistic, different p (sleepstudy Days: mixeff LRT chi² p=1.23e-06 vs lmerTest F p=3.26e-06). For the Satterthwaite/KR drop-one-term F, the user must instead use `anova()`/`test_effect()`. No `ddf`/`force_get_contrasts` args. |
| `ranova()` / `rand()` REML-LRT for random terms | partial | partial | major | `test_random_effect()` is the analog but differs substantively: (1) it uses the **Self-Liang 50:50 boundary mixture** reference, not ranova's plain (conservative) chi-square; (2) it tests **exactly one variance/covariance parameter** and refuses multi-parameter drops — `test_random_effect(sleepstudy_slope_model, "(1+Days\|Subject)")` returns `status="not_assessed"`, `p=NA` (theta_parameters=3, can't drop the whole slope block), whereas `ranova` reports the full `Days in (Days\|Subject)` 2-df test (LRT=42.84, p=4.99e-10). (3) No table over all random terms at once; one term per call. So the common "is the random slope needed?" ranova row is **not reproducible**. |
| `ranova(reduce.terms=)` toggle | in-scope-missing | in-scope-missing | minor | No `reduce.terms` analog; `test_random_effect` only drops one parameter and refuses correlated-block reductions. |
| `contest(joint=TRUE)` multi-df F (contestMD) | partial | partial | major | `contrast()` has **no `joint` argument** — it always returns per-row t-tests. Multi-row `L` under satterthwaite yields per-row t, never a collapsed F. The joint/contestMD path exists only inside term `anova()` and (per above) is unsupported under Satterthwaite for multi-df. A user wanting "joint test of these 2 contrasts" via Satterthwaite F has no route. |
| `contest(joint=FALSE)` per-row t (contest1D) | works | works | — | `contrast(m, L, method="satterthwaite")` gives per-row Estimate/SE/df/t/p; matches lmerTest contest1D semantics (sleepstudy Days contrast df 16.98, p 3.30e-06). |
| `contest1D()` single-df t | works | works | — | Covered by `contrast()` single-row; df/t/p match. |
| `contestMD()` multi-df F | partial | partial | major | Same as contest joint: KR multi-df F works through `test_effect()`/`anova()`; Satterthwaite multi-df F unsupported; no general contestMD over arbitrary `L` (only term-keyed families). `rhs` vector supported per-row. |
| `contest(confint=, level=)` CI columns | works | works | — | `contrast()` rows carry estimate/SE/df; CIs available via `mm_means`/`mm_predictions` (`conf_low`/`conf_high`) and `confint()`. |
| `contest(check_estimability=)` | works | works | — | `estimability(mm_lmm, L)` exposes the upstream estimability certificate (status/rank/requested_rank/reason). Always assessed in `contrast()` rows too. |
| `calcSatterth(model, L)` raw ddf | works | works | — | `df_for_contrast(m, L, method="satterthwaite")` returns per-contrast Satterthwaite df with method/reason attrs; df matches `contrast()` df. |
| `ls_means()` / `lsmeansLT()` | works | works | — | `mm_means(m, ~factor)` matches lmerTest exactly: `cake` recipe means 33.122/31.644/31.600, SE 1.7368, df 42, identical t/p. Equal-weight averaging matches LS-means semantics; also offers `weights="proportional"`. |
| `difflsmeans()` (pairwise) | works | works | — | `mm_comparisons(m, ~factor)` produces pairwise differences of marginal means with the same Satterthwaite df machinery (routes through `contrast()`). |
| `ls_means(which=)` factor subset | works | works | — | `mm_means(m, ~ recipe)` selects the displayed dimension; `by=`/`at=` extend it. |
| `ls_means(ddf="Kenward-Roger")` | works | works | — | `mm_means(..., method="kenward_roger")` supported (same KR path). |
| `show_tests()` (anova/ls_means hypothesis matrices) | in-scope-missing | in-scope-missing | major | No exported `show_tests`. An internal `mm_show_tests` symbol exists in the namespace but is **not exported** and not wired to `mm_anova`/`mm_means`: `attr(anova(m), "hypotheses")` is `NULL`. Given mixeff's audit-first mandate (every claim traces to an artifact), the inability to inspect the exact contrast matrix behind each ANOVA/ls_means row is a notable transparency gap. The per-row `details$contrast_family` exposes rank/family metadata but not the L matrix itself. |
| `step()` backward elimination | out-of-scope-by-design | out-of-scope-by-design | — | PRD §3 (lines 44–50): "No model selection or random-effects recommendation engine… never rank, select, or substitute models." `step()` is automatic model selection; deliberately excluded. R9 risk register asserts no "recommend"/"should"/"drop the random slope" language. |
| `get_model()` (final model from step) | out-of-scope-by-design | out-of-scope-by-design | — | Only meaningful with `step()`; same §3 exclusion. |
| `as_lmerModLmerTest(model)` upgrade lme4 fit | out-of-scope-by-design | out-of-scope-by-design | — | PRD §3 line 43: mixeff does not mask/interoperate with `lme4::lmer`. There is no `lmerMod` ingestion path; mixeff fits its own `mm_lmm` via the Rust engine. No coercion from a foreign object class. |
| `lmerModLmerTest` class / cached Satterthwaite slots | out-of-scope-by-design | out-of-scope-by-design | — | mixeff's `mm_lmm` is JSON-artifact-backed (CLAUDE.md: "JSON artifacts are the source of truth"); Satterthwaite machinery lives in Rust, computed on demand, not cached as S4 slots. Not a defect — different architecture. |
| Kenward-Roger beyond scalar | works (within scope) | works | — | KR multi-df F is implemented and matches lmerTest (`cake`, `ChickWeight`). Note PRD §3 historically listed KR as a v0 non-goal (line 41) but it is now implemented; survey scope note allows KR. |
| `pbkrtest` soft-dependency for KR | works | works | — | mixeff KR is native (Rust), no `pbkrtest` dependency; parity tests gate on `mm_skip_if_no_pbkrtest()` only for the reference side. |
| Clearer errors than lme4 | works | works | — | Refusals are structured (`status`, `reason`, `reason_code`, `reliability`) rather than silent. e.g. multi-df Satterthwaite returns an explicit `"...not implemented in this scaffold"` reason rather than a number — honest, though the underlying capability is still missing. |

---

## Highest-impact gaps (a real lmerTest user will hit these)

1. **Multi-df fixed-effect F under Satterthwaite is unimplemented** (blocker).
   `anova(model)` — lmerTest's single most common call — silently leaves every
   multi-level factor and interaction term with `p_value = NA` /
   `status = "unsupported"` unless the user switches to `method="kenward_roger"`.
   Only single-df terms get a Satterthwaite p. Reason string: *"multi-df
   fixed-effect contrast tests are not implemented in this scaffold."* This is
   the dominant parity hole in the family.

2. **`anova(type=)` does not change the hypothesis** (major). Type I/II/III are
   labels; the contrast matrix is never re-derived per SS type. Results are
   correct only for designs where the types coincide; users on unbalanced
   designs with interactions will get a mislabeled (single-parameterization)
   answer rather than the requested SS type.

3. **`summary()` default is Wald-z, not Satterthwaite** (major). df/p columns
   are absent by default; the lme4-trained user must opt in with `method=`.

4. **`ranova` semantics differ** (major). `test_random_effect` is a
   one-parameter boundary-mixture test, not ranova's whole-block REML-LRT;
   the canonical "is the random slope needed?" 2-df row cannot be reproduced.

5. **`drop1` is an LRT, not a Satterthwaite F-table** (major), and
   **`contrast()` has no `joint`/contestMD F** (major).

6. **`show_tests()` is unexported / unwired** (major for an audit-first
   package): hypothesis matrices behind ANOVA/ls_means rows are not retrievable.

Out-of-scope by PRD §3: `step()`, `get_model()`, `as_lmerModLmerTest()`, the
`lmerModLmerTest` S4 class and cached slots.
