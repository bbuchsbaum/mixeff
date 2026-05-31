# Gap Analysis — Datasets & utilities

**Family:** Datasets & utilities (shipped datasets; `mkReTrms`/`mkNewReTrms`;
`dummy`; `GHrule`/`GQdk`/`GQN`; modular fitting `lFormula`/`glFormula`/
`mkLmerDevfun`/`optimizeLmer`/`mkMerMod`; formula-manipulation utilities;
simulation templates; diagnostics utilities `isSingular`/`getSingTol`/`rePCA`/
`ngrps`; convergence utilities `allFit`/`checkConv`; ~35 internal/advanced
exports; lmerTest utilities).

**Reference:** `assessment/survey/lme4-datasets-utils.md` (lme4 2.0.1,
lmerTest 3.2.1).

**Date:** 2026-05-31. **Verdict:** Largely missing by design; the few overlaps
(`is_singular`, `getME`, `mm_formula_manifest`, `mm_parse_formula`) are partial
or serve a different purpose than the lme4 originals.

---

## Evidence baseline

- `NAMESPACE` exports no dataset, no `findbars`/`nobars`/`subbars`/`reOnly`/
  `isNested`, no `mkReTrms`/`mkNewReTrms`, no `dummy`, no `GHrule`/`GQdk`/`GQN`,
  no `lFormula`/`glFormula`/`mkLmerDevfun`/`optimizeLmer`/`mkMerMod`/
  `updateGlmerDevfun`, no `isSingular`/`getSingTol`/`rePCA`/`ngrps`, no
  `allFit`/`checkConv`, no simulation templates, no lmerTest `step`/
  `as_lmerModLmerTest`/`show_tests`/`get_model`.
- `data(package="mixeff")$results` is empty → **no shipped datasets**.
  `exists(<fn>, asNamespace("mixeff"))` is `FALSE` for every formula/RE/modular/
  quadrature/dummy utility listed above.
- mixeff overlaps that DO exist: `is_singular()` (S3), `getME()` (S3 + default),
  `mm_parse_formula()`, `mm_formula_manifest()` (zero-arg capability report),
  `VarCorr`, `model.matrix`, `model.frame`.

All `lmm()` examples below run against `lme4::sleepstudy` (mixeff carries no copy).

---

## Capability table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| **Shipped datasets** (`sleepstudy`, `cbpp`, `cbpp2`, `Pastes`, `Penicillin`, `Dyestuff`, `Dyestuff2`, `cake`, `grouseticks`, `grouseticks_agg`, `VerbAgg`, `InstEval`, `Arabidopsis`, `salamander`, `schizophrenia`, `toenail`, `gopherdat2`, `culcita*`) | missing | out-of-scope-by-design | minor | `data(package="mixeff")$results` empty. mixeff is a thin wrapper; users load data from `lme4`/`lmerTest`, which the env already provides. Tests/vignettes use `lme4::sleepstudy` directly. PRD §3 frames mixeff as a non-drop-in wrapper, not a dataset republisher. A real user incurs only `library(lme4)`; no functional loss. |
| **lmerTest datasets** (`carrots`, `ham`, `TVbo`) | missing | out-of-scope-by-design | minor | Same as above; available via `library(lmerTest)`. |
| `findbars(term)` | missing | in-scope-missing | minor | Not exported/defined. mixeff parses formulas internally (Rust) and surfaces structure via `mm_parse_formula()` + the random-effects manifest, but offers no programmatic `findbars`-style RE-term extractor for arbitrary formulas. A user scripting against RE terms (common in package code) has no replacement. |
| `nobars(term)` | missing | in-scope-missing | minor | Not present. No fixed-effects-only formula extractor. |
| `subbars(term)` | missing | in-scope-missing | minor | Not present. |
| `expandDoubleVerts(term)` | missing | partial (capability exists internally) | minor | `||` is a supported formula operator (see formula-re-syntax family / `mm_formula_manifest()$formula_features`), so mixeff *expands* `||` during compile, but exposes no standalone `expandDoubleVerts` inspector returning the expanded formula. |
| `reOnly(f, ...)` | missing | in-scope-missing | minor | Not present. |
| `isNested(f1, f2)` | missing | in-scope-missing | minor | Not present. Note: mixeff's design-support layer (`audit_design`, random-effects manifest "support" line) describes nesting/support qualitatively but offers no boolean `isNested(a,b)` predicate over two factors. A user verifying design structure programmatically has no direct analog. |
| `mkReTrms(bars, fr, ...)` | missing | out-of-scope-by-design | minor | RE-term construction (Zt/theta/Lind/Gp/flist) happens inside the Rust crate; not exposed as an R constructor. This is the core hook of lme4's *modular fitting API*, which mixeff does not offer (single-call `lmm()`/`glmm()` only). No PRD §3 line names `mkReTrms` explicitly, but the modular pipeline as a whole is absent by architecture (Rust owns assembly). Closest surface: `getME(fit,"theta")` works; `getME(fit,"flist")` works. `getME(fit,"Zt")` is **broken** (see below). |
| `mkNewReTrms(object, newdata, ...)` | missing | partial | minor | No exported analog, but the *capability* (RE construction for new data / new levels) is delivered internally by `predict.mm_lmm`/`predict.mm_glmm`. Users get prediction at new levels through `predict()`, not through a callable `mkNewReTrms`. |
| **Modular LMM pipeline** `lFormula` / `mkLmerDevfun` / `optimizeLmer` | missing | out-of-scope-by-design | major | None exist. There is no way to obtain a raw deviance function, intercept the parsed model frame before optimization, or run a custom optimizer. mixeff's `compile_model()` + `lmm()` are a single fused path with the optimizer inside Rust. The PRD's audit-first design (every fit step crosses the boundary as a structured artifact) is the substitute philosophy, but a user who relies on `mkLmerDevfun` to profile/visualize the likelihood surface or swap optimizers has no path. Marked **major** because power users and package authors genuinely use this. Not literally enumerated in PRD §3, but consistent with "not a drop-in replacement." |
| **Modular GLMM pipeline** `glFormula` / `mkGlmerDevfun` / `optimizeGlmer` / `updateGlmerDevfun` | missing | out-of-scope-by-design | major | None exist. Same architecture: GLMM fitting is fused inside Rust (`glmm()`); some Laplace/AGQ paths are additionally `nlopt`-feature-gated for CRAN. No user-accessible nAGQ-staged deviance function. |
| `mkMerMod(rho, opt, ...)` | missing | out-of-scope-by-design | minor | mixeff returns `mm_lmm`/`mm_glmm` objects, not `merMod`; there is no modular pipeline whose results would need assembling, so no analog. |
| `GHrule(ord)` | missing | in-scope-missing | minor | Not exported/defined. No way to inspect/verify the Gauss-Hermite rule mixeff (Rust) uses for AGQ. A user verifying quadrature has no R-side rule generator. |
| `GQdk(d, k)` / `GQN` | missing | in-scope-missing | cosmetic | Not present; rarely used directly. |
| `dummy(f, levelsToKeep)` | missing | in-scope-missing | minor | Not exported/defined. The lme4 helper is used to build contrast-specific RE design columns, e.g. `(0 + dummy(Sex,"Female") | Subject)`. Because mixeff cannot consume `dummy()` inside a formula, that idiom is unavailable; whether mixeff supports an equivalent factor-slope encoding is a formula-syntax question outside this family, but the named utility is absent. |
| `mkDataTemplate` / `mkParsTemplate` (simulation templates) | missing | out-of-scope-by-design | minor | Not present. mixeff has `simulate.mm_lmm` (simulate-refit family) but no template-from-formula generators. Power/sample-size scaffolding is not a stated v1 goal. |
| `isSingular(object, tol)` | **works (renamed)** | works | — | `is_singular(fit)` is exported and returns a logical. `is_singular(lmm(Reaction~Days+(Days|Subject), sleepstudy))` → `FALSE`. Name differs from lme4 (`is_singular` vs `isSingular`); a user typing `isSingular(fit)` gets "could not find function". Functionally equivalent for the core query. |
| `getSingTol()` | missing | partial | minor | No exported tolerance accessor. `is_singular()` presumably uses a fixed/internal threshold; users cannot inspect or override it the way `getSingTol()`/`options()` allow. |
| `rePCA(x)` | missing | in-scope-missing | major | No analog. lme4's `rePCA` is the standard tool (Bates et al. 2015) for diagnosing *which* RE components are near-zero, not just whether the fit is singular. mixeff's `is_singular()` gives only the boolean; the diagnostic richness (component SDs, rotation) is absent. Real lme4 users diagnosing over-parameterized RE structures will hit this. The PRD's RE-guidance layer (§9.5, "describe consequences / what each option assumes is zero") is adjacent but is passive/qualitative, not a PCA of the fitted covariance. |
| `ngrps(object)` | missing | partial | minor | No exported `ngrps`. The number of grouping-factor levels IS surfaced qualitatively in the RE manifest ("group levels: 18") and `getME(fit,"flist")` returns the grouping factors (so `sapply(getME(fit,"flist"), nlevels)` is a workaround), but there is no direct `ngrps(fit)` returning the named integer vector. |
| `allFit(object, ...)` | missing | out-of-scope-by-design | minor | No multi-optimizer refit harness. mixeff has a single Rust optimizer plus an `optimizer_certificate()` artifact attesting convergence; the audit-first design replaces "re-fit with every optimizer and compare" with a structured convergence certificate. Robustness-across-optimizers checking is not available. |
| `checkConv(...)` | missing | out-of-scope-by-design | cosmetic | Internal lme4 convergence checker; mixeff's equivalent is the `optimizer_certificate` / `diagnostics` artifact. |
| `getME(object, name)` | **partial** | partial | major | Exported with an `mm_lmm` method. `getME(fit,"theta")` → correct theta vector; `getME(fit,"flist")` → `"Subject"`. **But** `getME(fit,"Zt")` ERRORS: `Error in t.default(stats::model.matrix(object, type = "random")): argument is not a matrix`. So the general accessor covers some of the ~30 lme4 names but at least one documented modular name (`Zt`) is broken. This is a real defect a user would hit. (Cross-ref matrix-accessors family.) |
| `isLMM` / `isGLMM` / `isNLMM` / `isREML` | missing | partial | minor | No exported predicates. Class is testable via `inherits(fit,"mm_lmm")`/`inherits(fit,"mm_glmm")`; REML state is recoverable from `mm_control`/fit metadata but not via an `isREML()` predicate. NLMM is out of scope (nlmer non-goal, PRD §3). |
| `getL`, `getReCovs`, `anyStructured`, `isNewMerMod`, `forceNewMerMod`, `devcomp`, `llikAIC` | missing | out-of-scope-by-design | cosmetic | Low-level `merMod`-internal accessors; mixeff has no `merMod`. `VarCorr(fit)` covers the RE-covariance need; `logLik`/`AIC`/`BIC` cover `llikAIC`. |
| `mkVarCorr`, `formatVC`, `cov2sdcor`/`sdcor2cov`, `Cv_to_Sv` family, `mlist2vec`/`vec2mlist`/`vec2STlist` | missing | out-of-scope-by-design | cosmetic | Covariance-parameterization plumbing for lme4's theta encoding. mixeff exposes its own `parameterization()` / `theta_map` artifacts instead; these specific converters are package-author internals. |
| `varianceProf` / `logProf` / `golden` | missing | out-of-scope-by-design | cosmetic | Profile-CI plumbing; profile-likelihood CIs for GLMM are deferred (PRD §3) and these are lme4-profile-object internals. |
| `Nelder_Mead`/`NelderMead`, `nlminbwrap`, `nloptwrap`, `rePos`, `merPredD`, `lmerResp`/`glmResp`/`lmResp`/`nlsResp`, `glmFamily` | missing | out-of-scope-by-design | cosmetic | Optimizer wrappers and C++-backed reference classes specific to lme4's R/C++ bridge. mixeff's optimizer lives in Rust; these have no meaning here. |
| `namedList`, `factorize`, `lme4_testlevel` | missing | out-of-scope-by-design | cosmetic | Minor lme4 internals. |
| `as_lmerModLmerTest(model)` | missing | out-of-scope-by-design | minor | Converts a plain `lmerMod` to an lmerTest object. mixeff fits are not `lmerMod`; Satterthwaite/KR inference is built into mixeff's own inference layer (inference family), so no conversion step is needed/possible. |
| `get_model(step_object)` | missing | out-of-scope-by-design | minor | Tied to lmerTest `step()` backward elimination — a model-selection routine explicitly disclaimed (PRD §3: "No model selection ... random-effects recommendation engine"). |
| `show_tests(object)` | missing | partial | minor | lmerTest's contrast-matrix display for ANOVA tests. mixeff has rich estimability/contrast surfaces (`estimability()`, `contrast()`, `df_for_contrast()`, `show_tests`-adjacent reporting), but no single `show_tests`-named verb dumping the L-matrices for an ANOVA table. (Cross-ref inference family.) |

---

## mixeff utilities with no direct lme4 counterpart (context, not gaps)

- `mm_parse_formula(formula)` → returns canonicalized formula string (e.g.
  `"Reaction ~ 1 + Days + (1 + Days | Subject)"`). Loosely overlaps the
  `findbars`/`nobars` *intent* (formula introspection) but returns a normalized
  string, not extractable term lists.
- `mm_formula_manifest()` → **zero-argument** machine-readable capability/schema
  self-report (`mixeff_rust_version`, `crate_version`, `schema_versions`,
  `formula_features`, `capabilities`). Not an lme4 analog; it is the audit-first
  design's "what does this build know how to do?" answer. (Confirmed: it takes
  no formula/data args — `mm_formula_manifest(formula, data)` errors with
  "unused arguments".)
- Random-effects manifest / `explain_model` printout surfaces grouping levels,
  theta count, covariance type, and design support — qualitatively covering the
  *intent* of `ngrps`, `isNested`, and partially `rePCA`, but as prose, not as
  the structured return values lme4 programmers consume.

---

## Severity rationale

- The standard is "everything lme4 does, faster, with clearer errors." For this
  family most items are deliberately out of scope: mixeff does not republish
  datasets (trivially available from lme4/lmerTest) and does not expose lme4's
  R/C++ modular pipeline (fitting is fused inside Rust by architecture).
- **Major** gaps a real user hits: (1) the entire user-accessible **modular
  fitting pipeline** (`lFormula`/`mkLmerDevfun`/`optimizeLmer` and GLMM peers)
  — power users and package authors rely on it; (2) **`rePCA`** — the canonical
  RE-singularity *diagnosis* tool, only the boolean `is_singular()` is offered;
  (3) **`getME(fit,"Zt")` is broken** — a documented accessor name that errors.
- **In-scope-missing minor** items (`findbars`/`nobars`/`subbars`/`reOnly`/
  `isNested`/`dummy`/`GHrule`/`getSingTol`) are formula/diagnostic conveniences a
  scripting user could reasonably expect; each has a workaround or is rarely used.
- Naming: `is_singular` (vs `isSingular`) means copy-pasted lme4 code fails with
  "could not find function" rather than a clearer redirect — a small clarity gap
  against the project's "clearer errors" promise.
