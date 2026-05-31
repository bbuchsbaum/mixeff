# Gap Report: Matrix & Structure Accessors

**Family:** Matrix & structure accessors (`getME`, `model.matrix`, `model.frame`,
`terms`, `vcov`, `weights`, `ngrps`, `fitted`, `hatvalues`, `sigma`, `nobs`, `coef`, …)
**Reference:** `assessment/survey/lme4-matrix-accessors.md` (lme4 2.0.1, lmerTest 3.2.1)
**Assessed:** 2026-05-31 against installed `mixeff` (fit fn is `lmm()` / `glmm()`,
classes `mm_lmm`/`mm_glmm` inheriting `mm_fit`, `mm_compiled`)
**Standard:** "everything lme4 does, faster, with clearer errors."

## How this was tested

```r
library(mixeff); library(lme4); library(emmeans)
data(sleepstudy, package = "lme4")
mf <- lmm(Reaction ~ Days + (Days | Subject), data = sleepstudy)   # mixeff
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)         # lme4
```

Every status below was confirmed by running the accessor on `mf`, not by reading
source alone.

---

## Summary table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `getME(., "X")` | works | works | — | `getME(mf,"X")` → `matrix` (180×2), correct dimnames |
| `getME(., "Z")` | works | works | — | `getME(mf,"Z")` → `dgCMatrix` (180×36), has colnames |
| `getME(., "Zt")` | **errors (advertised but broken)** | partial | **major** | `getME(mf,"Zt")` → `Error: argument is not a matrix`. Switch arm calls bare `t()` on a `dgCMatrix`; resolves to `base::t` (Matrix S4 `t` not dispatched in pkg ns). `Matrix::t(Z)` succeeds. Root cause in `R/revive.R:128`. |
| `getME(., "Lambda")` | works | works | — | `getME(mf,"Lambda")` → `dtCMatrix` |
| `getME(., "Lambdat")` | **errors (advertised but broken)** | partial | **major** | Same bug as Zt: `R/revive.R:130` `t(.mm_lazy(...))` → `Error: argument is not a matrix` |
| `getME(., "theta")` | works | works | — | `getME(mf,"theta")` → numeric (3) |
| `getME(., "beta")` / `"fixef"` | works | works | — | numeric vector returned |
| `getME(., "y")` | works | works | — | numeric (180) |
| `getME(., "mu")` | works | works | — | numeric (180), = `fitted(mf)` |
| `getME(., "flist")` | works | works | — | `mm_flist` object |
| `getME(., "cnms")` | works | works | — | `mm_cnms` object |
| `getME(., "sigma")` | **errors** | partial | minor | `getME(mf,"sigma")` → not-available error, though `sigma(mf)` works (25.x). lme4 docs call `getME(.,"sigma")` a secondary path, but it is a documented name; absent here. |
| `getME(., "u")` (spherical RE mode) | missing | in-scope-missing | major | `getME(mf,"u")` → "not available". Conditional modes exist internally (ranef works) but the spherical `u` is not exposed. |
| `getME(., "b")` (RE mode = Λu) | missing | in-scope-missing | major | `getME(mf,"b")` → "not available". Stacked `b` vector not exposed (though `ranef()` gives per-term tables). |
| `getME(., "Gp")` | missing | in-scope-missing | major | "not available" — groups pointer; needed to slice Z/u/Lambda by term. |
| `getME(., "Tp")` | missing | in-scope-missing | minor | theta pointer not exposed. |
| `getME(., "L")` (sparse Cholesky) | missing | in-scope-missing | major | "not available" — needed for leverage / log-det / hatvalues. |
| `getME(., "Lind")` | missing | in-scope-missing | minor | theta→Lambda index map not exposed. |
| `getME(., "Tlist")` / `"ST"` | missing | in-scope-missing | minor | template / S-T factorisation not exposed. |
| `getME(., "A")` (= ZtΛ') | missing | in-scope-missing | minor | scaled sparse model matrix not exposed. |
| `getME(., "RX")` | missing | in-scope-missing | major | Cholesky of profiled X'X; the basis lme4 uses for `vcov`. mixeff returns vcov from a payload instead (works), but `RX` itself is unavailable for users doing manual SE / projection math. |
| `getME(., "RZX")` | missing | in-scope-missing | minor | cross-term block not exposed. |
| `getME(., "Ztlist")` / `"mmList"` | missing | in-scope-missing | major | per-term Z' blocks / raw per-term model matrices not exposed; also blocks `model.matrix(type="randomListRaw")`. |
| `getME(., "Tlist")`, `"p_i"`, `"l_i"`, `"q_i"`, `"m_i"`, `"k"`, `"m"` | missing | in-scope-missing | minor | term-dimension scalars/vectors not exposed. |
| `getME(., "N"/"n"/"p"/"q")` | missing | in-scope-missing | minor | basic dimension scalars not exposed (recoverable via `nobs`, `length(fixef)`, `dim(getME(.,"Z"))`, but not by name). |
| `getME(., "n_rtrms"/"n_rfacs")` | missing | in-scope-missing | minor | RE term / factor counts not exposed. |
| `getME(., "REML"/"is_REML")` | missing | in-scope-missing | minor | REML flag not exposed via getME (REML status may be elsewhere; `isREML` generic not provided either). |
| `getME(., "devcomp")` (`$cmp`,`$dims`) | missing | in-scope-missing | major | "not available". Deviance decomposition (ldL2, wrss, pwrss, sigmaML/REML) and dims vector unavailable. Used widely by downstream tooling. |
| `getME(., "offset")` | missing | in-scope-missing | minor | offset vector not exposed. |
| `getME(., "lower")` | missing | in-scope-missing | major | theta lower bounds not exposed; lme4's `isSingular()` reads this. mixeff has `is_singular()` (own path) so user-facing singularity works, but the `lower` vector itself is unavailable. |
| `getME(., "par"/"devarg"/"devfun")` | missing | in-scope-missing / partial | minor | optimizer par / deviance closure not exposed. `devfun` (profiling closure) is conceptually out-of-scope-by-design given the audit-first JSON model, but no clear error explains its absence. |
| `getME(., "ALL")` | missing | in-scope-missing | major | lme4's omnibus dump unsupported; `getME(mf,"ALL")` → "not available". A real lme4 user reflexively calls this. |
| `getME()` clear-error quality | partial | partial | minor | Unknown names give a clean typed `mm_arg_error` ("`getME()` component `X` is not available."), which **meets** the clearer-errors bar — but `Zt`/`Lambdat` throw the cryptic base-R `"argument is not a matrix"` instead, violating it. |
| `model.matrix(., type="fixed")` | works | works | — | `matrix` 180×2 |
| `model.matrix(., type="random")` | works | works | — | `dgCMatrix` 180×36 with colnames |
| `model.matrix(., type="randomListRaw")` | missing | in-scope-missing | minor | `match.arg` rejects: `Error: 'arg' should be one of "fixed", "random"`. Per-term raw matrices unsupported. |
| `model.matrix(noScale=)` | n/a | out-of-scope-by-design | cosmetic | autoscaling/noScale not part of mixeff's API surface; no scaling round-trip to undo. |
| `model.frame(., fixed.only=)` | partial | partial | minor | `model.frame(mf)` → `data.frame` (works). But signature is `model.frame.mm_lmm <- function(formula, ...)`; `fixed.only` argument is not honored (ignored, no fixed-only subsetting). |
| `terms(., fixed.only/random.only)` | **missing** | in-scope-missing | **major** | `terms(mf)` → `Error in terms.default: no terms component nor attribute`. No `terms.mm_lmm` method. emmeans still works (`emmeans(mf, ~Days)` returns an `emmGrid`) via dedicated `recover_data`/`emm_basis`, so the gap is the standalone `terms()` accessor, not emmeans integration. |
| `formula(., fixed.only/random.only)` | partial | partial | minor | `formula(mf)` → `formula` (works, full formula). `fixed.only`/`random.only` arguments appear unsupported (not in this family's core but worth noting). |
| `fixef(., add.dropped/noScale)` | works (core) | works | — | returns named numeric vector (named, unlike `getME(.,"beta")`). `add.dropped`/`noScale` flags not assessed here (rank-deficiency family). |
| `vcov(., correlation/full/use.hessian/noScale)` | partial | partial | major | `vcov(mf)` → p×p `matrix` with correct dimnames (works). But `correlation=TRUE` is silently ignored (`attr(.,"factors")` is NULL — no correlation matrix attached). `full=TRUE` (joint (q+p)×(q+p)) unsupported. `use.hessian` unsupported. mixeff adds a non-lme4 `type=c("fixed","theta")` arg. Silent ignore of `correlation` is a clear-error violation. |
| `sigma(.)` | works | works | — | `sigma(mf)` → scalar matching lme4 (25.x) |
| `fitted(.)` | works | works | — | numeric (180) |
| `residuals(., type/scaled)` | partial | partial | minor | `residuals(mf)` works (numeric 180). `type`/`scaled` coverage belongs to the diagnostics family; core path works. |
| `hatvalues(., fullHatMatrix)` | **missing** | in-scope-missing | **major** | `hatvalues(mf)` → `Error: no applicable method for 'hatvalues'`. No method. Leverage diagnostics (Cook's distance workflows) unavailable. Well-defined for LMMs, so in scope. |
| `weights(., type)` | partial | partial | major | `weights(mf)` → **`NULL`**. lme4 returns an all-ones numeric vector of length N for an unweighted fit (`head(weights(fm))` = `1 1 1 1 1 1`). Returning NULL breaks downstream code (e.g. `weighted.mean`, broom) that expects a numeric vector. `type=c("prior","working")` unsupported. |
| `ngrps(.)` | **missing** | in-scope-missing | **major** | `ngrps(mf)` → `Error: Cannot extract the number of groups from this object`. lme4 returns named integer `c(Subject = 18)`. Core RE-structure summary; commonly used. |
| `nobs(.)` | works | works | — | `nobs(mf)` → integer 180 |
| `coef(.)` | works | works | — | `coef(mf)` → `mm_coef` (β + b_i per level) |
| `getData()` / `fortify.merMod()` | missing | in-scope-missing / out-of-scope | minor | no `getData`/`fortify`; broom.mixed-style augmentation not provided. `getData` is a thin convenience (low severity); `fortify` is ggplot2 glue (arguably out-of-scope but no replacement documented). |

---

## Classification notes & PRD alignment

- **None of these gaps are covered by PRD §3 non-goals.** §3 defers profile-LL
  GLMM CIs, multivariate `cbind()`, Kenward-Roger beyond scalar, AR(1)/spatial,
  nlmer, `I()`/`poly()`/splines/GAM, and the model-selection engine. Matrix and
  structure accessors are not among the deferred items, so missing accessors are
  **in-scope-missing**, not out-of-scope.
- The two genuinely defensible out-of-scope items are `getME(.,"devfun")`
  (a live deviance closure conflicts with the audit-first/JSON-artifact model)
  and `fortify`/`noScale` (ggplot2 glue / scaling round-trip mixeff doesn't do).
  Even for these, mixeff offers no *clear typed error explaining the omission* —
  it just falls through to the generic "not available", which is acceptable but
  not the deliberate diagnostic the project's "no silent surgery" principle aims
  for.
- **Nothing here is upstream-blocked with evidence.** The crate already exposes
  X, Z, Lambda, theta, conditional-variance arrays, and the fixed-effect
  covariance payload, so the FE/RE design and Cholesky machinery exist Rust-side.
  Whether `L`, `RX`, `devcomp`, `Gp`, `u`/`b` are cheaply serializable from the
  current FFI would need a crate-contract check before promoting any to
  upstream-blocked; on current evidence they are unwired R-side, i.e.
  in-scope-missing.

## Highest-impact gaps a real lme4 user hits immediately

1. **`getME(.,"Zt")` / `getME(.,"Lambdat")` throw `"argument is not a matrix"`**
   (partial, major) — advertised names that error with a cryptic base-R message.
   Bug, not a missing feature. Fix: use `Matrix::t()` (or `methods::t`) in
   `R/revive.R:128,130`.
2. **`ngrps()` errors** (in-scope-missing, major) — trivial to derive from
   `flist`; a near-universal RE summary call.
3. **`weights()` returns `NULL`** (partial, major) — should return all-ones
   numeric of length N for unweighted fits to match lme4 and not break
   downstream consumers.
4. **`terms()` errors** (in-scope-missing, major) — no `terms.mm_lmm`; emmeans
   is wired separately so this is the bare accessor gap.
5. **`hatvalues()` has no method** (in-scope-missing, major) — leverage
   diagnostics unavailable for LMMs.
6. **`getME(.,"ALL")` and ~30 component names unsupported** (in-scope-missing,
   mostly major/minor) — including `u`, `b`, `Gp`, `L`, `RX`, `devcomp`,
   `lower`, `Ztlist`/`mmList`, dimension scalars.
7. **`vcov(correlation=TRUE)` silently ignored** (partial, major) — violates the
   project's no-silent-behavior principle; either honor it or raise a typed
   diagnostic.

## What genuinely works (parity met)

`getME` for {X, Z, Lambda, theta, beta, fixef, y, mu, flist, cnms};
`model.matrix(type=fixed|random)` with correct dimnames; `model.frame`;
`vcov` (point covariance, correct dimnames); `sigma` (matches lme4);
`fitted`; `residuals` (core); `nobs`; `coef`; `formula` (full); `fixef` (named).
Unknown-name `getME` errors are clean typed `mm_arg_error` conditions — the
clearer-errors bar is met *except* for the Zt/Lambdat base-R leak.
