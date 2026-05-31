# lme4 / lmerTest — Random-Effects Formula Syntax: Reference Survey

**Date:** 2026-05-31  
**Scope:** Every user-facing function, operator, argument, and behaviour that
lme4 (≥ 2.0) and lmerTest expose for specifying the *random-effects* part of a
mixed-model formula.  This document is a pure reference against which mixeff
will be evaluated; it does not assess mixeff.

---

## 1. The Bar Operator — core syntax

The `|` ("bar") operator is the fundamental delimiter of a random-effects term
inside a mixed-model formula.

```
(lhs | grouping_factor)
```

| Syntax | What it expresses | Why a user needs it |
|--------|-------------------|---------------------|
| `(1 \| g)` | Random intercept per level of `g`; scalar variance, no slope. | Baseline random effect for grouped/clustered data. |
| `(x \| g)` | Correlated random intercept **and** random slope for `x` per `g`; unstructured 2×2 covariance matrix (default `us` structure). | Models per-group variation in both baseline and the effect of `x`. |
| `(1 + x \| g)` | Identical to `(x \| g)` — the `1` makes the implicit intercept explicit. | Explicit style preferred by some users; equivalent to `(x \| g)`. |
| `(0 + x \| g)` or `(-1 + x \| g)` | Random slope for `x` with **no** random intercept; only the slope varies by group. | Suppresses the random intercept; useful when intercept variation is not meaningful or when avoiding identification issues. |
| `(x \| g) + (1 \| g)` | Two separate blocks for the same grouping factor — treated as distinct covariance terms. | Allows different covariance structures for different parts of the RE for the same grouping factor. |

### 1.1 Implicit intercept rule

Inside a random-effects term, R's standard formula rules apply: `1` is added
implicitly.  `(x|g)` therefore expands to `(1+x|g)` — a random intercept plus
a random slope.  Suppress the intercept with `0+` or `-1+` on the left of `|`.

---

## 2. The Double-Bar Operator `||`

```
(lhs || g)
```

Specifies **multiple uncorrelated random effects** for the same grouping factor.
lme4 ≥ 2.0 provides two expansion methods:

### 2.1 "split" method (default, backward-compatible)

`(1 + x || g)` expands to `(1|g) + (0+x|g)` — separate scalar blocks.  Works
correctly for **continuous** predictors only.  For factors, `(0+f||g)` does not
produce a diagonal covariance over contrasts; see §3 for `diag()`.

```r
expandDoubleVerts(y ~ x + (a+b||g))
# => y ~ x + ((1|g) + (0+a|g) + (0+b|g))
```

### 2.2 "diag_special" method (recommended for factors)

Set `options(lme4.doublevert.default = "diag_special")` (or call
`getDoublevertDefault()`).  Then `(1+a+b||g)` becomes `diag(1+a+b|g)` — a
single block with a diagonal (zero-correlation) covariance matrix.  Matches
`afex::mixed()` and `glmmTMB` behaviour.

### 2.3 `getDoublevertDefault()` / `options(lme4.doublevert.default)`

| Function / option | Purpose |
|-------------------|---------|
| `getDoublevertDefault()` | Returns `"split"` (default) or `"diag_special"`. |
| `options(lme4.doublevert.default = "diag_special")` | Sets the global double-bar expansion method for the session. |

---

## 3. Covariance-Structure Tags (lme4 ≥ 2.0)

Named function wrappers that precede the RE term and control the shape of the
covariance matrix.  All use the `tag(lhs | g)` syntax.

| Tag | Full syntax | Covariance constraint | Notes |
|-----|-------------|----------------------|-------|
| `us` | `us(x \| g)` | Unstructured positive-semidefinite (default). | `(x\|g)` and `us(x\|g)` are identical; the tag is explicit-only. |
| `diag` | `diag(x \| g)` | Diagonal — all off-diagonal covariances fixed to zero. | `diag(x\|g, hom=TRUE)` imposes equal variances on the diagonal (homogeneous). |
| `cs` | `cs(x \| g)` | Compound-symmetric — one shared variance, one shared covariance. | `cs(x\|g, hom=TRUE)` for homogeneous variances. |
| `ar1` | `ar1(x \| g)` | AR(1) — exponentially decaying correlation. Homogeneous by default. | `ar1(x\|g, hom=FALSE)` for heterogeneous variances. |

These are accessed via `splitForm()` / `findbars_x()` from the `reformulas`
package.  `VarCorr()` output labels the structure (e.g., `(cs)`, `(ar1)`).

---

## 4. Grouping-Factor Operators

Operators applied **to the right of `|`** that create composite or hierarchical
grouping factors.

| Syntax | Expansion | Semantics | Why a user needs it |
|--------|-----------|-----------|---------------------|
| `(1 \| a/b)` | `(1 \| b:a) + (1 \| a)` | Nesting: `b` nested within `a`.  Adds a random intercept for each unique `a`, plus one for each `b`-within-`a` combination. | Standard hierarchical data (e.g., students within schools). |
| `(1 \| a/b/c)` | `(1 \| c:b:a) + (1 \| b:a) + (1 \| a)` | Three-level nesting, expanding to three intercept blocks. | Three-level hierarchies (e.g., measurements within subjects within sites). |
| `(1 \| a:b)` | Single block on the interaction factor `a:b`. | Crossed grouping by the combination of `a` and `b`. | Explicitly model the cross-classified combination without all main effects. |
| `(1 \| a*b)` | `(1 \| a) + (1 \| b) + (1 \| a:b)` | Full crossing — main effects plus interaction, three separate intercept blocks. | Simultaneously model row effects, column effects, and cell effects in a two-way cross-classification. |

lme4 uses `isNested(f1, f2)` (from `reformulas`) to check empirically whether
every level of `f1` occurs within exactly one level of `f2`.

---

## 5. Multiple Grouping Factors (Cross-Classified)

```r
y ~ x + (1|subject) + (1|item)
```

Multiple `(... | g)` terms separated by `+` each generate an independent
covariance block.  Groups need not be nested — this is the standard
cross-classified random effects model.

- Any number of grouping factors may be included.
- Each may carry its own covariance structure tag.
- `mkReTrms()` collects all terms and orders blocks by decreasing number of
  groups (controlled by `reorder.terms`).

---

## 6. Intercept-Suppression in the Fixed Part

Standard R formula rules govern the fixed part of the model.  They affect how
the *random-effects design matrix* is interpreted when the same predictor
appears on both sides.

| Fixed-part syntax | Effect |
|-------------------|--------|
| `y ~ x + (1\|g)` | Fixed intercept present; random intercept is per-group deviation from it. |
| `y ~ 0 + x + (1\|g)` | Fixed intercept suppressed; all group-level variation absorbed by the random intercept (unusual; rarely meaningful). |
| `y ~ -1 + x + (1\|g)` | Same as `0+`. |

---

## 7. `dummy()` — Selective Dummy Coding

```r
dummy(f, levelsToKeep)
```

Wrapper around `model.matrix` that returns `nlevels(f)-1` columns (or a
specified subset of levels) as a numeric indicator matrix.  Designed for use
**inside** a random-effects term when the user wants to allow different random
slopes for specific levels of a factor without the usual treatment-contrast
reparameterisation.

```r
lmer(distance ~ age + (age|Subject) + (0+dummy(Sex,"Female")|Subject),
     data=Orthodont)
```

The column named `dummy(Sex, "Female")` appears as a random effect column in
the `cnms` slot of the fitted model.

---

## 8. `factor()` Inside a Random-Effects Term

lme4 permits `factor(x)` on the LHS of `|`:

```r
(factor(grp) | g)    # unstructured covariance over treatment contrasts
(0+factor(grp) | g)  # one variance per level (indicator parameterisation)
```

The two forms differ in parameterisation:
- `(factor(grp)|g)` — intercept + contrasts, treatment-coded (default contrasts).
- `(0+factor(grp)|g)` — one indicator per level, all variances free.

The `||` double-bar does **not** produce the intended diagonal covariance for
factors under the default `"split"` method; use `diag(factor(grp)|g)` or set
`"diag_special"` mode.

---

## 9. Formula-Processing Utility Functions

These are the programmer-facing tools that expose the formula-parsing pipeline.
They now live in the `reformulas` package; lme4 re-exports them as deprecated
shims.

### 9.1 `findbars(term)` / `findbars_x(...)`

| Function | Signature | What it returns | User-facing purpose |
|----------|-----------|-----------------|---------------------|
| `findbars(term)` | `term`: formula or language | List of `lhs \| rhs` expressions, one per RE block, after slash-expansion. | Extract all random-effects sub-formulae from a full model formula; foundation for programmatic formula manipulation. |
| `findbars_x(term, debug, specials, default.special, target, expand_doublevert_method)` | Same + expansion control | List of RE expressions with structure tags applied. | Richer version used internally; `expand_doublevert_method` = `"split"` or `"diag_special"` controls `||` handling. |

Note: `(1|a/b)` → `[[1 | b:a], [1 | a]]`; slash expansion happens inside
`findbars`.

### 9.2 `nobars(term)` / `nobars_(term)`

Removes all random-effects terms from a formula, returning the fixed-effects
formula.

```r
nobars(Reaction ~ Days + (Days|Subject))
# => Reaction ~ Days
```

Used by lme4 to build the fixed-effects model matrix and by downstream tools
that need the fixed specification.

### 9.3 `subbars(term)`

Replaces `|` and `||` with `+` throughout the formula, producing a formula
suitable for `model.frame()` (which does not understand the bar operators).

```r
subbars(Reaction ~ Days + (Days|Subject))
# => Reaction ~ Days + (Days + Subject)
```

### 9.4 `expandDoubleVerts(term)`

Expands all `||` terms in a formula to separate `|` terms (the "split" method).

```r
expandDoubleVerts(y ~ x + (a+b||g))
# => y ~ x + ((1|g) + (0+a|g) + (0+b|g))
```

### 9.5 `isNested(f1, f2)`

Predicate: does every level of factor `f1` occur with exactly one level of
factor `f2`?  Returns `TRUE`/`FALSE`.  Used to check structural nesting in
data before fitting.

### 9.6 `mkReTrms(bars, fr, drop.unused.levels, reorder.terms, reorder.vars, calc.lambdat, sparse)`

Builds the full random-effects model matrices from `findbars` output and a
model frame.  Returns: `Zt` (sparse transpose RE design matrix), `Ztlist`,
`Lambdat` (sparse relative covariance factor), `Lind`, `theta`, `lower`,
`flist`, `cnms`, `Gp`, `nl`, `ord`.  Central to the lme4 computational
machinery; also used by glmmTMB and other packages.

### 9.7 `splitForm(formula, defaultTerm, allowFixedOnly, allowNoSpecials, debug, specials)`

Parses a formula into `fixedFormula`, `reTrmFormulas`, `reTrmAddArgs`, and
`reTrmClasses`.  Handles "special" RE terms of the form `foo(x|g)` (`diag`,
`cs`, `ar1`, `us`, `rr`, etc.).

### 9.8 `reOnly(f, response, bracket, doublevert_split)`

Extracts only the random-effects part of a formula, optionally including the
response.

```r
reOnly(y ~ x + (1|g))           # => ~(us(1 | g))
reOnly(y ~ x + (1|g), response=TRUE)  # => y ~ (us(1 | g))
```

### 9.9 `expandAllGrpVar(bb)`

Expands a list of bare `1|f` grouping expressions, applying slash/interaction
logic.  Used internally after `findbars`.

### 9.10 `randint(formula)`

Adds a `(1|g)` random-intercept term to a formula if it does not already
contain one.  Minor convenience for programmatic formula building.

### 9.11 `noSpecials(term, delete, debug, specials)`

Removes or identifies "special" RE-term wrappers (`diag(...)`, `cs(...)`, etc.)
from a formula.  Companion to `splitForm`.

---

## 10. Modular Pipeline Functions

These expose the internal formula → model-matrix pipeline for advanced users.

| Function | Signature | Purpose |
|----------|-----------|---------|
| `lFormula(formula, data, REML, ...)` | Full lmer args | Parses formula and builds the model frame + RE terms structure without fitting. Returns a list with `fr`, `X`, `reTrms`, `REML`, `formula`, `wmsgs`. |
| `glFormula(formula, data, family, ...)` | Full glmer args | Same for GLMMs. |
| `mkReTrms(bars, fr, ...)` | See §9.6 | Builds sparse RE matrices from parsed bars. |
| `mkNewReTrms(object, newdata, ...)` | fitted model + new data | Rebuilds RE design matrices for prediction. |

---

## 11. lmerTest Additions Relevant to RE Syntax

lmerTest wraps lme4's `lmer()` and adds inference tools.  It does not extend
the formula syntax per se but provides functions that operate on the random
structure.

| Function | Signature | Purpose |
|----------|-----------|---------|
| `ranova(model, reduce.terms=TRUE, ...)` | `lmerModLmerTest` | ANOVA-like table of random-effect term deletions/reductions via LRT.  Alias: `rand()`. |
| `rand(model, ...)` | Same | Alias for `ranova`. |

### `ranova` reduction rules

| RE term form | `reduce.terms=TRUE` action |
|--------------|---------------------------|
| `(f1 + f2 \| gr)` | Reduced to `(f2\|gr)` and `(f1\|gr)` separately; LRT vs original. |
| `(f1 \| gr)` | Reduced to `(1\|gr)`. |
| `(1 \| gr)` | Simply removed (cannot reduce further). |
| `(0 + f1 \| gr)` or `(-1 + f1 \| gr)` | Reduced to `(1\|gr)`. |
| `(1 \| gr1/gr2)` | Auto-expanded to `(1\|gr2:gr1) + (1\|gr1)` via `findbars_x`. |
| Structured tags (`diag`, `cs`, `ar1`) | `reduce.terms` ignored; term is always removed. |

---

## 12. VarCorr — Variance-Covariance Extraction

```r
VarCorr(x, sigma=1, ...)
as.data.frame(VarCorr(x), order=c("cov.last","lower.tri"))
```

`VarCorr` extracts all estimated variance/covariance components.  The result is
a list of matrices (one per grouping factor), each carrying `stddev` and
`correlation` attributes.  Key attributes on the object:

- `"theta"` — covariance parameters on the Cholesky scale.
- `"rho"` — correlation parameters for structured terms.
- `"profpar"` — profiling parameterisation of theta.

`as.data.frame` produces a tidy row-per-parameter table, column ordering
controlled by `order`.

---

## 13. ranef — Conditional Modes

```r
ranef(object, condVar=TRUE, drop=FALSE, whichel=names(ans), postVar=FALSE, ...)
```

Returns a list of data frames (one per grouping factor) of conditional modes
(= conditional means for LMMs).  With `condVar=TRUE` attaches a `"postVar"`
array of posterior variances — key for caterpillar plots and uncertainty
quantification.  S3 methods: `dotplot.ranef.mer`, `qqmath.ranef.mer`,
`as.data.frame.ranef.mer`.

---

## 14. `dummy()` — Detail

```r
dummy(f, levelsToKeep)
```

Arguments:

| Argument | Type | Effect |
|----------|------|--------|
| `f` | any object coercible to `factor` | Input factor. |
| `levelsToKeep` | character vector (optional) | Subset of levels to keep; default drops first level. |

Returns a numeric model matrix.  Intended use: `(0+dummy(Sex,"Female")|g)`
produces a single indicator column as a RE predictor without the full
treatment-contrast matrix.

---

## 15. Summary: Supported Formula Building-Blocks

| Category | Syntax forms |
|----------|-------------|
| **Intercept** | `1` (implicit or explicit), `0`, `-1` (suppress) |
| **Single slope** | `(x\|g)`, `(1+x\|g)`, `(0+x\|g)`, `(-1+x\|g)` |
| **Multiple slopes** | `(a+b\|g)`, `(1+a+b\|g)` |
| **Double-bar uncorrelated** | `(x\|\|g)`, `(a+b\|\|g)` |
| **Nesting** | `(1\|a/b)`, `(1\|a/b/c)` |
| **Interaction grouping** | `(1\|a:b)`, `(x\|a:b)` |
| **Full crossing** | `(1\|a*b)` |
| **Multiple grouping factors** | `(1\|a) + (1\|b)` (cross-classified) |
| **Covariance tags** | `us()`, `diag()`, `diag(hom=TRUE)`, `cs()`, `cs(hom=TRUE)`, `ar1()`, `ar1(hom=FALSE)` |
| **Inline transforms** | `dummy(f, levels)`, `factor(x)` on LHS of `\|` |
| **Formula utilities** | `findbars`, `findbars_x`, `nobars`, `subbars`, `expandDoubleVerts`, `isNested`, `mkReTrms`, `splitForm`, `reOnly`, `expandAllGrpVar`, `randint`, `noSpecials` |
| **lmerTest RE tools** | `ranova()` / `rand()` |

---

## 16. Authoritative Sources Used

- `help("lmer", package="lme4")` — formula argument documentation and covariance structure tags.
- `help("findbars", package="reformulas")` — findbars_x full signature.
- `help("expandDoubleVerts", package="reformulas")`.
- `help("nobars", package="reformulas")`.
- `help("subbars", package="reformulas")`.
- `help("isNested", package="reformulas")`.
- `help("mkReTrms", package="reformulas")`.
- `help("splitForm", package="reformulas")`.
- `help("dummy", package="lme4")`.
- `help("getDoublevertDefault", package="lme4")`.
- `help("ranef", package="lme4")`.
- `help("VarCorr", package="lme4")`.
- `help("lFormula", package="lme4")`.
- `help("ranova", package="lmerTest")`.
- Live `Rscript` verification of all formula patterns cited above.
- `reformulas` package exports: `getNamespaceExports("reformulas")`.
