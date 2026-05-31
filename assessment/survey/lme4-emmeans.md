# lme4 / lmerTest — emmeans / Marginal Means Surface Reference

**Purpose:** Exhaustive reference of every user-facing function, argument, and
behavior that lme4 (via emmeans) and lmerTest offer in the "emmeans / marginal
means" capability family. This is the reference baseline against which mixeff
is judged. It does not assess mixeff; it documents what lme4/lmerTest provide.

Packages surveyed (installed versions):

| Package   | Version |
|-----------|---------|
| lme4      | 2.0.1   |
| lmerTest  | 3.2.1   |
| emmeans   | 2.0.3   |
| pbkrtest  | (optional, Kenward-Roger) |

---

## 1. Integration Architecture

### 1.1 `recover_data.merMod` (emmeans internal, registered S3)

**What it does:** Reconstructs the model-frame data needed to build the
reference grid from a fitted `merMod` (lme4 model object).  Refuses nonlinear
mixed models (`nlmer`) with an explicit error.

**Signature:**
```r
recover_data(object, ...)
```

**How called:** Automatically by `ref_grid()` and `emmeans()`.

**Why users rely on it:** This is the invisible plumbing that lets emmeans
operate on lme4 objects without any user action.

---

### 1.2 `emm_basis.merMod` (emmeans internal, registered S3)

**What it does:** Constructs the linear-functional basis (design matrix `X`,
fixed-effect estimates `bhat`, covariance `V`, df function `dffun`) from a
`merMod` object.  Dispatches on `isLMM` vs `isGLMM` and selects the df method.

**Signature:**
```r
emm_basis(object, trms, xlev, grid,
          mode            = get_emm_option("lmer.df"),
          lmer.df,                               # deprecated alias for mode
          disable.pbkrtest = get_emm_option("disable.pbkrtest"),
          pbkrtest.limit   = get_emm_option("pbkrtest.limit"),
          disable.lmerTest = get_emm_option("disable.lmerTest"),
          lmerTest.limit   = get_emm_option("lmerTest.limit"),
          options, ...)
```

**`mode` / `lmer.df` valid values:**

| Value | Df calculation |
|-------|---------------|
| `"kenward-roger"` | Kenward-Roger adjusted V + `pbkrtest::Lb_ddf()`. Default. |
| `"satterthwaite"` | Satterthwaite via `lmerTest::calcSatterth()`. |
| `"asymptotic"` | df = Inf; z-test. |

**Automatic fallback chain:** kenward-roger → satterthwaite → asymptotic when
the required package is unavailable or the dataset exceeds the size limit.

**Size limits (rows in model):**
- `pbkrtest.limit` default **3000** — above this KR is skipped with a note.
- `lmerTest.limit` default **3000** — above this Satterthwaite is skipped.

**Non-estimable cells:** Uses `estimability::nonest.basis()` on a rank-deficient
model matrix so non-estimable contrasts surface as `nonEst` in the output table
rather than silently producing garbage numbers.

**GLMM behavior:** df = Inf (asymptotic z); link labels forwarded to emmeans via
`.std.link.labels()` so `type = "response"` back-transforms automatically.

---

## 2. Core emmeans Functions

### 2.1 `ref_grid()`

**What it does:** Builds an `emmGrid` reference grid from a fitted model by
constructing the Cartesian product of all predictor levels/values.

**Signature:**
```r
ref_grid(object,
         at,                              # named list — override predictor values
         cov.reduce   = mean,             # how to reduce numeric covariates
         cov.keep     = get_emm_option("cov.keep"),  # keep covariates with ≤N levels
         mult.names,                      # for multivariate responses
         mult.levs,
         options      = get_emm_option("ref_grid"),
         data,                            # override training data
         df,                              # override degrees of freedom
         type,                            # prediction type ("response","link",…)
         regrid,                          # transform the grid
         nesting,                         # declare nesting structure
         offset,
         sigma,                           # residual SD for effect sizes / PIs
         counterfactuals,
         nuisance     = character(0),     # nuisance variables (averaged out)
         non.nuisance,
         wt.nuis      = "equal",
         rg.limit     = get_emm_option("rg.limit"),
         ...)
```

**Why users rely on it:** Provides a reproducible, inspectable reference point
before computing marginal means; allows fine-grained control over what the
"average" covariate value is; supports `at=` to create custom prediction
profiles.

**Return value:** An `emmGrid` S4 object with slots:
`model.info`, `roles`, `grid`, `levels`, `matlevs`, `linfct`, `bhat`, `nbasis`,
`V`, `dffun`, `dfargs`, `misc`, `post.beta`.

---

### 2.2 `emmeans()`

**What it does:** Computes estimated marginal means (least-squares means) for
factor levels or combinations, averaging over or conditioning on other variables.
Calls `ref_grid()` internally if needed.

**Signature:**
```r
emmeans(object, specs,
        by         = NULL,         # conditioning variable(s)
        fac.reduce = function(coefs) apply(coefs, 2, mean),
        contr,                     # optional contrast specification
        options    = get_emm_option("emmeans"),
        weights,                   # cell-weighting: "equal", "proportional", "outer", "cells", "flat", numeric vector
        offset,
        ...,
        tran)                      # transformation override
```

**Key `specs` forms:**

| Form | Meaning |
|------|---------|
| `~ trt` | Marginal means over `trt` |
| `~ trt * time` | Cell means for all combinations |
| `~ trt \| time` | Means for `trt` within each level of `time` |
| `pairwise ~ trt` | Means + pairwise contrasts in one call |
| `trt.vs.ctrl ~ trt` | Means + treatment-vs-control contrasts |

**`mode` / `lmer.df` (passed through `...` to `emm_basis`):**
Values `"kenward-roger"`, `"satterthwaite"`, `"asymptotic"` — see §1.2.

**`weights` values:**
`"equal"` (default in most contexts), `"proportional"`, `"outer"`, `"cells"`,
`"flat"`, numeric vector.

**Return value:** `emmGrid` S4 object (same class as `ref_grid()` output).

---

### 2.3 `contrast()` / `contrast.emmGrid()`

**What it does:** Computes arbitrary linear contrasts from an `emmGrid`.

**Signature:**
```r
contrast(object,
         method        = "eff",         # contrast method or named list of L vectors
         interaction   = FALSE,         # interaction contrasts
         by,
         offset        = NULL,
         scale         = NULL,
         name          = "contrast",
         options       = get_emm_option("contrast"),
         type,
         adjust,                        # p-value adjustment
         simple,                        # "each factor" for interaction grids
         combine       = FALSE,
         ratios        = TRUE,
         parens,
         enhance.levels = TRUE,
         wts,
         ...)
```

**Built-in contrast methods (`.emmc` functions):**

| Method | Description |
|--------|-------------|
| `"pairwise"` | All pairwise differences (A-B, A-C, B-C, …) |
| `"revpairwise"` | Reversed pairwise (B-A, C-A, …) |
| `"tukey"` | Tukey HSD comparisons |
| `"dunnett"` | Dunnett's test vs. a control |
| `"trt.vs.ctrl"` | Each treatment vs. last level (control) |
| `"trt.vs.ctrl1"` | Each treatment vs. first level |
| `"trt.vs.ctrlk"` | Each treatment vs. last level |
| `"eff"` | Each mean vs. grand mean |
| `"del.eff"` | Delete-one effect contrasts |
| `"helmert"` | Helmert contrasts |
| `"poly"` | Polynomial contrasts |
| `"opoly"` | Orthogonal polynomial contrasts |
| `"consec"` | Consecutive differences |
| `"mean_chg"` | Mean change (before/after style) |
| `"identity"` | Identity (keep rows as-is) |

**Custom contrasts:** Pass a named list of numeric vectors as `method`.

---

### 2.4 `pairs()` / `pairs.emmGrid()`

**What it does:** Shortcut for `contrast(object, method="pairwise")`.

**Signature:**
```r
pairs(x, reverse = FALSE, ...)
```

**`reverse = TRUE`:** Computes B-A, C-A, … instead of A-B, A-C, …

---

### 2.5 `joint_tests()`

**What it does:** Omnibus F (or chi-square) tests for each model term, analogous
to `anova()` output but computed through the emmeans linear-contrast machinery.
Uses the reference grid to construct the tests so the df denominator is the same
Satterthwaite/KR df used for individual contrasts.

**Signature:**
```r
joint_tests(object,
            by        = NULL,
            show0df   = FALSE,
            showconf  = TRUE,
            cov.reduce = make.meanint(1),
            ...)
```

**Return:** Data frame with columns `model term`, `df1`, `df2`, `F.ratio`,
`p.value`.

---

### 2.6 `summary.emmGrid()`

**What it does:** Produces a formatted data frame of estimates, SEs, df,
confidence intervals, and (optionally) hypothesis-test p-values.

**Signature:**
```r
summary(object,
        infer,                       # c(CI=TRUE/FALSE, test=TRUE/FALSE)
        level,                       # confidence level; default 0.95
        adjust,                      # p-value adjustment method
        by,
        cross.adjust = "none",
        type,                        # "response", "link", "mu", "unlink", "log", "scale"
        df,                          # override df
        calc,                        # additional calculated columns
        null,                        # null hypothesis value
        delta,                       # equivalence margin
        side,                        # 0 (two-tailed), 1 (>), -1 (<)
        frequentist,
        bias.adjust = get_emm_option("back.bias.adj"),
        sigma,
        ...)
```

**`adjust` values (p-value MCC):**
Any method accepted by `p.adjust()` plus emmeans-specific:
`"tukey"`, `"scheffe"`, `"sidak"`, `"bonferroni"`, `"dunnettx"`, `"mvt"`,
`"none"`, `"fdr"`, `"holm"`, `"hochberg"`, `"hommel"`, `"BH"`, `"BY"`.

**`type` values:**
`"response"` (back-transform), `"link"` (linear predictor scale),
`"mu"` (mean function), `"unlink"` (without link transform),
`"log"`, `"scale"`.

---

### 2.7 `confint.emmGrid()`

**What it does:** Extracts confidence intervals; thin wrapper around `summary()`
with `infer = c(TRUE, FALSE)`.

**Signature:**
```r
confint(object, parm, level = 0.95, ...)
```

---

### 2.8 `test.emmGrid()`

**What it does:** Performs hypothesis tests; thin wrapper around `summary()` with
`infer = c(FALSE, TRUE)`. With `joint = TRUE` performs joint F-test.

**Signature:**
```r
test(object,
     null    = 0,
     joint   = FALSE,
     verbose = FALSE,
     rows,
     by,
     status  = FALSE,
     ...)
```

---

### 2.9 `update.emmGrid()`

**What it does:** Modifies slots of an `emmGrid` (e.g., change `adjust`, `level`,
`by`, `tran`, `type`, `df`, `sigma`, `infer`).

**Signature:**
```r
update(object, ..., silent = FALSE)
```

---

### 2.10 `predict.emmGrid()`

**What it does:** Returns point predictions, with optional confidence or
prediction intervals.

**Signature:**
```r
predict(object,
        type,
        interval = c("none", "confidence", "prediction"),
        level    = 0.95,
        bias.adjust = get_emm_option("back.bias.adj"),
        sigma,
        ...)
```

---

### 2.11 `plot.emmGrid()`

**What it does:** Plots marginal means with confidence intervals, optional
prediction intervals, and optionally overlays pairwise comparison arrows.

**Signature:**
```r
plot(x, y,
     type,
     CIs         = TRUE,
     PIs         = FALSE,
     comparisons = FALSE,       # overlay pairwise comparison arrows
     colors,
     alpha       = 0.05,
     adjust      = "tukey",
     int.adjust  = "none",
     intervals,
     ...)
```

---

### 2.12 `vcov.emmGrid()`

**What it does:** Returns the covariance matrix for the linear functions
(contrasts or means) in the grid.

**Signature:** `vcov(object, ...)`

---

### 2.13 `coef.emmGrid()`

**What it does:** Returns the linear-functional coefficient matrix (the `linfct`
slot, i.e., the contrast coefficient matrix mapping fixed effects to displayed
quantities). Returns `NULL` when the basis is too complex.

---

### 2.14 `as.data.frame.emmGrid()` / `as.data.frame.summary_emm()`

**What it does:** Converts an `emmGrid` or its summary to a plain data frame for
further manipulation.

---

### 2.15 `as.list.emmGrid()` / `as.emm_list()`

**What it does:** Converts an `emmGrid` to an `emm_list` for multi-component
results (e.g., `pairwise ~ trt` returns both means and contrasts as a list).

---

### 2.16 `rbind.emm_list()`

**What it does:** Pools multiple contrast sets from an `emm_list` and applies
a single multi-comparison adjustment across all of them.

**Signature:**
```r
rbind(..., which, adjust = "bonferroni")
```

---

### 2.17 `emtrends()`

**What it does:** Computes estimated marginal slopes (partial derivatives) of the
response with respect to a continuous variable, evaluated at reference-grid
values. Used to compare slopes across factor levels.

**Signature:**
```r
emtrends(object, specs, var,
         delta.var   = 0.001 * rng,   # finite-difference step size
         max.degree  = 1,              # polynomial degree
         ...)
```

---

### 2.18 `eff_size()`

**What it does:** Converts pairwise contrasts to standardized effect sizes
(Cohen's d style) using a pooled sigma and effective df.

**Signature:**
```r
eff_size(object, sigma, edf,
         method = "pairwise",
         ...)
```

---

### 2.19 `add_grouping()`

**What it does:** Adds a new factor to an `emmGrid` by mapping existing levels to
a coarser grouping, enabling nested or hierarchical marginal mean computations.

**Signature:**
```r
add_grouping(object, newname, refname, newlevs, ...)
```

---

### 2.20 `regrid()`

**What it does:** Re-grids (transforms) an existing `emmGrid` to a different
scale (e.g., log, response, sqrt) so that subsequent operations work on that
scale.  Used to obtain back-transformed CIs that respect the transformation.

**Signature:** Accessible via `update(object, regrid = ...)` or direct call.

---

### 2.21 `pwpm()`

**What it does:** Produces a compact pairwise P-value matrix (upper triangle:
p-values; diagonal: estimates; lower triangle: differences).

**Signature:**
```r
pwpm(emm,
     by,
     reverse = FALSE,
     pvals   = TRUE,
     means   = TRUE,
     diffs   = TRUE,
     flip    = FALSE,
     digits,
     ...)
```

---

### 2.22 `pwpp()`

**What it does:** Pairwise P-value plot — a significance-based graphical summary
of pairwise comparisons.

**Signature:**
```r
pwpp(emm,
     method   = "pairwise",
     by,
     sort     = TRUE,
     values   = TRUE,
     rows     = ".",
     xlab, ylab, xsub = "",
     plim     = numeric(0),
     add.space = 0,
     aes,
     ...)
```

---

### 2.23 `emmip()`

**What it does:** Interaction-profile plot (trace lines) for marginal means.

**Signature:**
```r
emmip(object, formula, ...)
```

---

### 2.24 `qdrg()`

**What it does:** "Quick-and-dirty reference grid" — constructs an emmGrid from
raw model components (coefficients, vcov, df) without a formal model object.
Useful for non-standard models.

**Signature:**
```r
qdrg(formula, data, coef, vcov, df, mcmc, object, subset,
     weights, contrasts, link, qr, ordinal, ...)
```

---

### 2.25 `add_submodels()`

**What it does:** Appends sub-model grids (e.g., from reduced models) to an
existing `emmGrid` for simultaneous multi-model display.

**Signature:**
```r
add_submodels(object, ..., newname = "model")
```

---

### 2.26 `as.glht()` (via multcomp)

**What it does:** Converts an `emmGrid` contrast set to a `glht` object
(multcomp package) for use with `summary.glht()`, `confint.glht()`, and
simultaneous inference machinery.

**Signature:** `as.glht(object, ...)` — requires **multcomp**.

---

### 2.27 `emm_options()` / `get_emm_option()`

**What it does:** Gets and sets package-wide options controlling default behavior.

**Signature:**
```r
emm_options(..., disable)
```

**Key options for mixed-model users:**

| Option | Default | Effect |
|--------|---------|--------|
| `lmer.df` | `"kenward-roger"` | Default df method for `merMod` |
| `pbkrtest.limit` | `3000` | Max N for Kenward-Roger |
| `lmerTest.limit` | `3000` | Max N for Satterthwaite |
| `disable.pbkrtest` | `FALSE` | Force-disable KR |
| `disable.lmerTest` | `FALSE` | Force-disable Satterthwaite |
| `cov.keep` | `"2"` | Covariates with ≤2 unique values kept as factors |
| `back.bias.adj` | `FALSE` | Back-transform bias correction |
| `save.ref_grid` | `FALSE` | Cache last ref_grid in session |
| `rg.limit` | `10000` | Max reference-grid rows before warning |

---

## 3. lmerTest-Specific Functions

lmerTest provides its own marginal-means surface that predates emmeans
integration. All functions below operate on `lmerModLmerTest` objects (produced
by `lmerTest::lmer()` or `as_lmerModLmerTest()`). Internally they delegate to
emmeans but expose a simpler, opinionated interface.

### 3.1 `lmerTest::ls_means()` / `lsmeansLT()`

**What it does:** Computes least-squares (marginal) means with Satterthwaite df
for all fixed-effect factors. `lsmeansLT()` is a deprecated alias.

**Signature:**
```r
ls_means(model, ...)
```

**Output:** A formatted table with columns `Estimate`, `Std. Error`, `df`,
`t value`, `lower`, `upper`, `Pr(>|t|)`. Confidence level 95%, df method
always Satterthwaite.

---

### 3.2 `lmerTest::difflsmeans()`

**What it does:** Computes all pairwise differences of least-squares means with
Satterthwaite df.

**Signature:**
```r
difflsmeans(model, ...)
```

**Output:** Formatted table with `Estimate`, `Std. Error`, `df`, `t value`,
`lower`, `upper`, `Pr(>|t|)`.

---

### 3.3 `lmerTest::contest()`

**What it does:** Tests a general linear hypothesis `L β = 0` via F-test
(multi-row `L`) or t-test (single-row `L`).  Df computed by Satterthwaite.
Dispatches to `contest1D()` or `contestMD()`.

**Signature:**
```r
contest(model, L, ...)
```

---

### 3.4 `lmerTest::contest1D()`

**What it does:** Single-df contrast test (t-test) with Satterthwaite df.

**Signature:**
```r
contest1D(model, L, ...)
```

**Output:** `Estimate`, `Std. Error`, `df`, `t value`, `Pr(>|t|)`.

---

### 3.5 `lmerTest::contestMD()`

**What it does:** Multi-df F-test for a hypothesis `L β = 0` where `L` is a
matrix. Satterthwaite denominator df.

**Signature:**
```r
contestMD(model, L, ...)
```

**Output:** `Sum Sq`, `Mean Sq`, `NumDF`, `DenDF`, `F value`, `Pr(>F)`.

---

### 3.6 `lmerTest::ranova()`

**What it does:** ANOVA-like table for random-effect terms — likelihood-ratio
tests for dropping each random-effect term. Related to marginal quantities
because it helps decide whether a grouping structure is needed at all.

**Signature:**
```r
ranova(model, reduce.terms = TRUE, ...)
```

**Output:** Table with `npar`, `logLik`, `AIC`, `LRT`, `Df`, `Pr(>Chisq)`.

---

### 3.7 `lmerTest::show_tests()`

**What it does:** Shows how the design matrix underlying an anova test is built;
diagnostic function for understanding what hypothesis is being tested.

**Signature:**
```r
show_tests(object, ...)
```

---

### 3.8 `lmerTest::calcSatterth()`

**What it does:** Low-level Satterthwaite df calculation for a single contrast
vector `L` and a fitted `lmerModLmerTest`.

**Signature:**
```r
calcSatterth(model, L)
```

**Return:** List with element `denom` (the denominator df scalar).

---

### 3.9 `lmerTest::step()`

**What it does:** Automated backward-elimination model selection operating on
both random and fixed effects, using Satterthwaite tests. Produces a sequence
of models and summary tables.

**Signature:**
```r
step(object, ...)
```

---

## 4. Non-Estimability Handling

lme4/emmeans handles rank-deficient models (e.g., a missing cell in a factorial
design) through `estimability::nonest.basis()`. Non-estimable contrasts appear
as `nonEst` in the output table (not `NA`, not an error). This applies to:

- `emmeans()` output
- `contrast()` output
- `joint_tests()` output

---

## 5. `emmGrid` S4 Class and Its Methods

The `emmGrid` class is the return type for `ref_grid()`, `emmeans()`, and
`contrast()`. User-accessible S3/S4 methods:

| Method | Description |
|--------|-------------|
| `[` | Subset the grid |
| `+` | Combine two grids |
| `as.data.frame` | Convert summary to data frame |
| `as.glht` | Convert to multcomp `glht` |
| `as.list` | Convert to `emm_list` |
| `coef` | Return linear-functional coefficient matrix |
| `confint` | Confidence intervals |
| `contrast` | Linear contrasts |
| `head` / `tail` | First/last rows of grid |
| `levels<-` | Replace level names |
| `pairs` | Pairwise contrasts |
| `plot` | Plot means with CIs/PIs/comparison arrows |
| `predict` | Point predictions and intervals |
| `print` / `show` | Display the grid |
| `rbind` | Pool multiple grids |
| `str` | Structural summary |
| `subset` | Subset by expression |
| `summary` | Formatted table of estimates, SEs, df, CIs, p-values |
| `test` | Hypothesis tests |
| `update` | Modify grid options |
| `vcov` | Covariance matrix for contrasts |
| `weights` | Cell weights |

---

## 6. The `emm_list` Class

When `emmeans()` is called with a combined spec like `pairwise ~ trt`, it returns
an `emm_list` (a named list of `emmGrid` objects). Extra methods:

| Method | Description |
|--------|-------------|
| `confint.glht_list` | CIs for list elements |
| `contrast.emm_list` | Apply contrast to list elements |
| `pairs.emm_list` | Pairwise contrasts |
| `plot.emm_list` | Plot all list elements |
| `rbind.emm_list` | Pool with joint adjustment |
| `summary.emm_list` | Summarize all elements |
| `as.data.frame.summary_eml` | Convert pooled summary to data frame |

---

## 7. df Handling Flow for lme4 LMMs

```
emmeans() / ref_grid()
  └─ emm_basis.merMod()
       ├─ isLMM → select mode
       │     ├─ "kenward-roger" → pbkrtest::Lb_ddf   (default, N ≤ 3000)
       │     ├─ "satterthwaite" → lmerTest::calcSatterth   (N ≤ 3000)
       │     └─ "asymptotic"   → df = Inf
       └─ isGLMM → df = Inf always
```

---

## 8. lsmeans / Legacy Aliases

`lsmeans()` was an earlier package superseded by emmeans. The emmeans package
provides `lsmeans` as a deprecated alias for `emmeans()` with identical
behavior.

---

## 9. Parity Assessment Reference

Against this surface, the key capabilities to check in mixeff are:

1. `recover_data.mm_lmm` + `recover_data.mm_glmm` — bridge plumbing
2. `emm_basis.mm_lmm` — df method, mode argument, V matrix, nbasis
3. `emm_basis.mm_glmm` — asymptotic df, link labels for `type="response"`
4. S3 registration via `emmeans::.emm_register()` (lme4 style hook)
5. `mode=` argument naming (`kenward-roger` vs `kenward_roger`)
6. `satterthwaite` explicit mode producing correct per-contrast df
7. Non-estimable cell detection via `estimability::nonest.basis()`
8. `lmerTest::ls_means()` / `difflsmeans()` equivalents (native surface)
9. `lmerTest::contest()` / `contest1D()` / `contestMD()` equivalents
10. `lmerTest::ranova()` equivalent (random-effect LRT)
11. `eff_size()` — requires `sigma.mm_lmm`
12. `emtrends()` — slope estimates across factor levels
13. `joint_tests()` — omnibus F-tests through emmeans machinery
14. `add_grouping()` — hierarchical groupings of factor levels
15. `regrid()` / back-transformation for GLMMs
16. `pbkrtest.limit` / `lmerTest.limit` size guards
17. `at=` covariate override for reference grid
18. `weights=` cell weighting (`"proportional"` vs `"equal"`)

---

*Generated: 2026-05-31. Sources: installed package introspection (emmeans 2.0.3,
lme4 2.0.1, lmerTest 3.2.1) via `Rscript` + `getAnywhere()` + live testing.*
