# GLMM Fitting — Gap Report

**Family:** GLMM fitting (glmer families, links, nAGQ/Laplace/AGQ, cbind binomial, weights, offset, dispersion)
**Reference:** `assessment/survey/lme4-glmm-fit.md` (lme4 2.0.1)
**mixeff entry point:** `glmm()` (`R/glmm.R`), FFI `wrap__mm_fit_glmm_json` (`src/rust/src/lib.rs:459`)
**Date:** 2026-05-31
**Standard:** everything lme4 does, faster, with clearer errors.

All findings below were confirmed by running `library(mixeff); library(lme4)` comparisons,
not by source reading alone. Repro scripts were run in `/tmp`.

---

## Executive summary

`glmm()` fits the **certified family/link matrix only**: binomial (logit/probit/cloglog),
poisson (log/sqrt), Gamma (log). For 0/1 Bernoulli, Poisson, Gamma it returns a fitted
`mm_glmm` with fixef, theta, dispersion, logLik, deviance, AIC/BIC, ranef, VarCorr, and a
Wald-z summary table. The numerical engine is the upstream **profiled fast-PIRLS** path, which
is *deliberately not* lme4's joint-Laplace fit — divergences are documented in
`inst/extdata/expected-mismatches.json` and disclaimed by PRD §3 (no bit-exact reproduction).

The two material defects a real `glmer` user would hit are:

1. **`offset` is silently dropped** (blocker / silent-surgery contract violation). The `offset`
   argument is not in the `glmm()` signature, so it falls into `...` and is discarded with no
   error. lme4 honors it; mixeff returns the no-offset fit. This violates CLAUDE.md "no silent
   surgery." In-fit offsets are a PRD §3 non-goal, so the *capability* is out of scope — but the
   *silent acceptance* is a bug regardless.
2. **`cbind(success, failure)` binomial responses are rejected** as an unsupported formula
   construct (`mm_formula_error`), with no documented host-side expansion path. This is a
   first-class lme4 idiom (the canonical cbpp example) and there is no `weights=`-denominator
   alternative either (weights is also rejected).

Everything else is either a clean typed refusal (good clear-error behavior) or an
out-of-scope-by-design deferral that is correctly surfaced.

---

## Capability-by-capability

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `glmer()` top-level fit | works (profiled fast-PIRLS) | works | — | `glmm(y~x+(1|g), family=binomial)` returns `mm_glmm`; fixef/theta/logLik/AIC populated. |
| binomial, logit link | works | works | — | Bernoulli 0/1 fit; fixef (-0.098, 0.932) vs lme4 (-0.108, 1.001). |
| binomial, probit link | works | works | — | `family=binomial(link="probit")` fits; fixef (0.417, 0.444). |
| binomial, cloglog link | works | works | — | `family=binomial(link="cloglog")` fits OK. |
| binomial, cauchit link | in-scope-missing → refused | partial | minor | Not in `mm_glmm_supported_family_links()`; refused with typed `mm_inference_unavailable`. lme4 supports cauchit. Clear error, but capability absent. |
| binomial, log link | in-scope-missing → refused | partial | minor | Same: refused, not in matrix. lme4 lists `log` as valid for binomial. |
| poisson, log link | works | works | — | Poisson fit; fixef (0.557, 0.517) vs lme4 (0.529, 0.517); theta 0.840 vs 0.840. |
| poisson, sqrt link | works | works | — | Allowed by matrix and engine (`LinkFunction::Sqrt`). |
| poisson, identity link | in-scope-missing → refused | partial | minor | Not in wrapper matrix; refused. lme4 supports it. |
| Gamma, log link | works | works | — | `Gamma(link="log")` fits; fixef + dispersion (0.63) returned. |
| Gamma, inverse link (canonical/default) | engine-capable, wrapper-blocked | in-scope-missing | major | `family=Gamma` (default inverse) → `mm_inference_unavailable` "Gamma/inverse outside certified contract." Crate survey crate-0.md §44 confirms the **engine certifies Gamma+Inverse** (`LinkFunction::Inverse` exists, line 1210); only `mm_glmm_supported_family_links()` blocks it. This is the *default* Gamma link in base R, so a naive `family=Gamma` user is refused for a link the engine supports. |
| Gamma, identity link | in-scope-missing → refused | partial | minor | Refused; lme4 supports. |
| inverse.gaussian family | in-scope-missing → refused | partial | minor | `family=inverse.gaussian` → typed refusal "inverse.gaussian/1/mu^2 outside certified contract." Engine has `Family::InverseGaussian` label (line 1244) but `glmm_family()` (line 1193) does not map it. Rare family; clear error. |
| gaussian via glmm (identity) | refused | out-of-scope-by-design | cosmetic | `family=gaussian` → typed refusal. Users route to `lmm()`. Reasonable. |
| `negative.binomial(theta)` fixed-theta family | in-scope-missing → refused | partial | minor | `family=MASS::negative.binomial(2)` → typed `mm_inference_unavailable`. No NB family in engine. |
| `glmer.nb()` (estimated theta) | absent | out-of-scope-by-design | minor | No `glmer.nb` analog. NB not in engine family enum. Overdispersed counts unsupported. Not in Phase 4 scope (PRD §10 line 806). Clear path: no verb exists, so no silent failure. |
| quasi / quasibinomial / quasipoisson | refused (clearer than lme4) | works | — | `family=quasipoisson` → typed refusal. lme4 also rejects quasi; mixeff message is clearer/typed. |
| custom `make.link()` link objects | refused | out-of-scope-by-design | cosmetic | Only the fixed link names map. Refused via matrix gate. lme4 passes arbitrary links through. Niche. |
| **`cbind(success, failure)` binomial response** | **rejected, no alternative** | **in-scope-missing** | **major** | `glmm(cbind(incidence,size-incidence)~period+(1|herd), family=binomial)` → `mm_formula_error`: "cbind(...) not in stateless transform subset." No host-side expansion, and `weights=` denominator path is also blocked. Grouped binomial (the canonical cbpp idiom) is unreachable except by manual Bernoulli expansion. Crate survey crate-3.md line 55 says engine supports "Binomial (grouped); trial weights supported" — so the gate is wrapper-side. |
| `weights` (prior/observation weights) | reserved → hard refusal | in-scope-missing | major | `weights=` in `glmm()` → `mm_fit_error` "reserved for the fitted GLMM bridge" (`R/glmm.R:54`). Crate survey says engine supports observation weights (crate-1.md line 39). Wrapper has not wired them. Hard error (good — not silent), but capability absent. Blocks grouped-binomial-via-weights and any weighted GLM. |
| **`offset` argument** | **SILENTLY DROPPED** | **in-scope-missing (silent-surgery bug)** | **blocker** | `glmm(...)` has no `offset` param; `offset=rep(3,n)` lands in `...` and is discarded. Repro: lme4 intercept with offset=3 → -2.47; mixeff intercept identical to no-offset (0.557). `fit$offset` is NULL. Violates CLAUDE.md "no silent surgery." In-fit offsets are PRD §3 non-goal, so the *feature* is out-of-scope, but **silently ignoring a supplied offset instead of refusing it is a defect.** Contrast: `weights` (also reserved) errors loudly; `offset` does not. |
| `offset()` in formula | refused (typed) | out-of-scope-by-design | minor | `y ~ x + offset(logexp) + (1|g)` → `mm_formula_error` "offset(...) not in stateless transform subset." PRD §3 non-goal (in-fit offsets, line 42). This path at least errors loudly, unlike the `offset=` argument. |
| nAGQ = 1 (Laplace) | works (profiled approx, not true Laplace) | partial | major | Default. Returns `method="pirls_profiled"`, not lme4 joint Laplace. Estimates diverge by documented tolerances (see fast-PIRLS row). |
| nAGQ = 0 (PQL-like fast stage) | absent | out-of-scope-by-design | minor | No `nAGQ=0` exposure; validator requires `nAGQ >= 1` (`mm_glmm_validate_nagq`, line 219). The profiled path *is* the fast path. Reasonable. |
| nAGQ > 1 (adaptive Gauss-Hermite) | accepted but NOT true AGQ | partial | major | `nAGQ=10` accepted on `pirls_profiled`; estimates shift slightly (logLik -187.52→-187.29) but this is **not lme4 AGQ** — it is profiled-PIRLS AGQ metadata. The true joint-AGQ path requires `method="joint_laplace"` + nlopt, which is disabled. Crate-0.md line 244 flags that profiled+nAGQ>1 behavior is undocumented in the R layer. A user expecting Gauss-Hermite quadrature accuracy gets a different quantity with no warning. |
| nAGQ > 1 single-scalar-RE constraint | not enforced on profiled path | partial | minor | lme4 requires exactly one scalar RE for nAGQ>1; mixeff accepts nAGQ>1 on profiled path regardless of RE structure (no `invalid_agq_request`). The constraint only bites on the (disabled) joint path. |
| `method="joint_laplace"` (true Laplace) | refused (nlopt disabled) | upstream-blocked | major | `method="joint_laplace"` → `mm_fit_error` "requires the upstream nlopt backend, disabled in this vendored build." Wired but unreachable in CRAN build. PRD §5.2/§10 R1: nlopt feature-gated for CRAN; available under `--features nlopt` for R-Universe/dev. Clear typed refusal. This is *the* path to lme4-parity GLMM fits, so its absence in the default build is significant. |
| true joint AGQ (`joint_laplace` + nAGQ>1) | refused (nlopt disabled) | upstream-blocked | major | Same refusal path; `mm_arg_error` if nAGQ>1 + joint, else nlopt error. Engine certifies it only for culcitalogreg rows (crate-1.md line 48). |
| dispersion: binomial/poisson fixed at 1 | works | works | — | `sigma()` returns dispersion; binomial/poisson canonical fixed. |
| dispersion: Gamma/inverse.gaussian estimated | partial | partial | minor | Gamma+log returns `dispersion` (0.63 in test); `fit$sigma == fit$dispersion`. inverse.gaussian unsupported, so its dispersion is unreachable. |
| observation-level RE for overdispersion `(1|obs)` | works (mechanism) | works | — | Standard RE term; supported via formula. Standard glmer overdispersion workaround available. |
| `glmerControl()` optimizer/convergence knobs | partial (mm_control) | out-of-scope-by-design | minor | mixeff uses `mm_control()`, not glmerControl. Optimizer/check knobs differ by design (audit-first, different engine). Not a 1:1 surface. |
| two-stage optimization (nAGQ0initStep) | n/a | out-of-scope-by-design | cosmetic | Engine-internal; not exposed. |
| `start`/`mustart`/`etastart` initialization | absent | out-of-scope-by-design | minor | No start-value passthrough in `glmm()`. Not in Phase 4 scope. |
| `subset`, `na.action`, `contrasts` | reserved → hard refusal | in-scope-missing | minor | All in the `R/glmm.R:54` reserved-reject list; loud `mm_fit_error`. lme4 supports all three. Loud refusal is acceptable per "no silent surgery," but parity gap. |
| `devFunOnly` | absent | out-of-scope-by-design | cosmetic | Low-level; not exposed. |
| modular API (glFormula/mkGlmerDevfun/optimizeGlmer) | absent | out-of-scope-by-design | cosmetic | Internal lme4 surface; not a mixeff goal. |
| `GHrule`/`GQN`/`GQdk` quadrature exports | absent | out-of-scope-by-design | cosmetic | Internal; not exposed. |
| `summary()` Wald-z table for GLMM | works | works | — | `summary(fit, tests="coefficients")` returns z value / Pr(>z) table with reliability note ("pirls_laplace_working_hessian"). |
| `fixef`/`ranef`/`VarCorr`/`sigma`/`vcov` | works | works | — | All populated in `mm_glmm` (see `R/glmm.R:106-149`). Covered by other family reports (extractors). |
| `predict`/`fitted`/`residuals` | works | works | — | `fitted`/`residuals` stored; `predict.mm_glmm` exists (`R/predict.R`). Covered by prediction report. |
| `anova`/`drop1` LRT | partial | partial | minor | Out of this family's strict scope (inference report); GLMM LRT support partial. |
| `confint(method="profile")` GLMM | absent | out-of-scope-by-design | minor | PRD §3: profile-LL CIs deferred to v2 for GLMM. |
| `bootMer` / parametric bootstrap GLMM | partial (engine has it; R wires lmm only) | in-scope-missing | major | Crate-0.md line 242: `parametricbootstrap_glmm` is **stable in the engine** (Bernoulli/Binomial/Poisson/Gamma) but R `parametric_bootstrap()` only accepts `mm_lmm`. GLMM bootstrap CIs unreachable from R. (Belongs primarily to the bootstrap report; noted here because it is the realistic GLMM uncertainty path given profile CIs are deferred.) |
| `isSingular`/`rePCA` for GLMM | partial | partial | minor | Diagnostics surface exists; coverage in diagnostics report. |
| `isGLMM`/`isLMM` predicates | works (class-based) | works | — | `inherits(fit,"mm_glmm")`; no `isGLMM()` verb but class is queryable. |

---

## Severity rationale

- **Blocker:** `offset=` silently dropped — violates the package's core "no silent surgery"
  contract (CLAUDE.md, PRD §3 line 35). Any Poisson rate model (events per person-year, the
  most common offset use) gives a wrong fit with no warning.
- **Major:** cbind binomial unreachable; `weights` unwired; Gamma/inverse (default link)
  wrapper-blocked despite engine support; nAGQ>1 silently means something other than AGQ;
  joint_laplace (the lme4-parity path) disabled in CRAN build; GLMM parametric bootstrap
  unwired. Each is something a routine `glmer` user hits.
- **Minor/cosmetic:** missing links (cauchit, identity), NB family, start values, modular API,
  control-knob parity — all niche or correctly surfaced as typed refusals.

## Classification notes vs PRD §3

- **In-fit offsets are an explicit PRD §3 non-goal (line 42).** So the offset *feature* is
  out-of-scope-by-design. The defect is not the absence of the feature; it is that the
  `offset=` **argument is silently swallowed** rather than refused. The formula-based
  `offset()` correctly errors; the argument path does not. Classified the silent-drop as
  in-scope-missing/blocker because the contract is "every refusal crosses the boundary."
- **No bit-exact reproduction (PRD §3 line 37)** covers the fast-PIRLS vs joint-Laplace
  divergence, which is documented in `inst/extdata/expected-mismatches.json`
  (`cbpp_binomial_logit_profiled_pirls`). That divergence is therefore *not* a defect, but the
  fact that the only lme4-parity path (joint_laplace) is nlopt-gated off makes the divergence
  the default behavior — worth the user's attention.
- **nlopt gating (PRD §5.2 R1)** makes joint_laplace/AGQ upstream-blocked in the CRAN build,
  available under `--features nlopt` for R-Universe/dev.

## Confirmed clear-error wins (the "clearer errors" half of the standard)

- quasi families: typed `mm_inference_unavailable` (lme4's is a bare stop).
- unsupported family/link: typed refusal carrying a `supported` table.
- joint_laplace without nlopt: typed `mm_fit_error` naming the missing backend, not a silent fallback.
- cbind / formula offset: typed `mm_formula_error` listing the allowed transform subset.

## Test-gap observations

- `tests/testthat/test-glmm.R` exercises ledger-bounded parity for the certified matrix but
  there is **no test asserting that a supplied `offset=` is either honored or refused** — the
  silent-drop is untested. Recommend a regression test (test-gap, major).
- No test for `weights=`, `cbind`, or `nAGQ>1`-vs-`nAGQ=1` divergence on the profiled path.
