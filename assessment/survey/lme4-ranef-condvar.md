# lme4 `ranef` / `condVar` & Plotting — Public Surface Survey

Survey date: 2026-05-31  
lme4 version: 2.0.1 (installed); lattice required for plot methods.  
Scope: every user-facing function, argument, and behavior in the
`ranef(condVar=TRUE)` / `postVar` / `as.data.frame.ranef.mer` /
`dotplot` / `qqmath` / `plot` / `print` family.

---

## 1. `ranef.merMod` — extracting conditional modes

```r
ranef(object, condVar = TRUE, drop = FALSE,
      whichel = names(ans), postVar = FALSE, ...)
```

Applies to both `lmerMod` and `glmerMod` objects (class `merMod`).

| Argument | Type | What it does / why users rely on it |
|---|---|---|
| `object` | `merMod` | The fitted LMM or GLMM object. |
| `condVar` | logical, default **`TRUE`** | When `TRUE`, attaches the conditional variance-covariance arrays to each element of the returned list as the `"postVar"` attribute. This is the core argument of the family — it is what enables uncertainty display in caterpillar plots and `as.data.frame` `condsd` column. |
| `drop` | logical, default `FALSE` | When `TRUE`, any grouping factor whose data frame has only a single column (typically just an intercept, `"(Intercept)"`) is simplified to a named numeric vector instead of a one-column data frame. Useful when working programmatically. |
| `whichel` | character vector | Subset of grouping factor names to return. Defaults to all grouping factors. Lets callers avoid computing/returning unused groups in models with many crossed random effects. |
| `postVar` | logical, default `FALSE` | Deprecated synonym for `condVar`. Retained for back-compatibility; triggers a warning in recent lme4. |

### Return value: `ranef.mer` object

`ranef()` returns an S3 object of class `"ranef.mer"` which is a named
**list of data frames**, one per grouping factor.

- Each data frame has **k rows** (one per level of the grouping factor)
  and **j columns** (one per random-effect term, named after the term,
  e.g. `"(Intercept)"`, `"Days"`).
- Row names are the levels of the grouping factor.

### The `"postVar"` attribute

When `condVar = TRUE`, each data frame in the list carries an attribute
named `"postVar"` (the name is a historical artifact; lme4 documentation
notes it may eventually be renamed to `"condVar"`).

**Single RE term per grouping factor** (the common case,
e.g. `(1 + Days | Subject)`):

- The attribute is a **3-D array of dimension `j × j × k`** where `j`
  is the number of random-effect columns and `k` is the number of
  grouping levels.
- The `i`-th face `pv[,,i]` is the **conditional variance-covariance
  matrix** for the `i`-th level (a symmetric, positive-semi-definite `j×j`
  matrix).  The diagonal entries are the conditional variances; off-diagonal
  entries are conditional covariances between RE terms for that level.
- For an LMM with a single grouping factor this array gives the diagonal
  blocks of the full block-diagonal conditional variance-covariance of the
  entire random-effects vector.  With multiple grouping factors the blocks
  are no longer block-diagonal globally, but each face is still the
  correct diagonal block.

**Multiple RE terms for the same grouping factor** (the separated form,
e.g. `(1|f) + (0+x|f)`):

- The attribute is a **list of arrays**, one per term, each of dimension
  `1 × 1 × k` (since each term is scalar in this parameterisation).

**Numerical meaning**: each diagonal element `pv[i, i, k]` is
`Var(b_{ik} | y, θ̂)` — the posterior variance of random-effect
component `i` for level `k`, conditional on the data and the estimated
variance-component parameters.  Its square root is the "conditional
standard deviation" (`condsd`) used in uncertainty visualisation.

---

## 2. `as.data.frame.ranef.mer` — long-format conversion

```r
as.data.frame(x, ...)
```

Converts a `ranef.mer` list into a single long-format data frame.
Used heavily to feed ggplot2 caterpillar plots.

| Output column | Description |
|---|---|
| `grpvar` | Name of the grouping factor (character). |
| `term` | Name of the random-effect term, e.g. `"(Intercept)"` or `"Days"` (character). |
| `grp` | Level of the grouping factor (factor/character). |
| `condval` | Value of the conditional mode (numeric). |
| `condsd` | Conditional standard deviation, i.e. `sqrt(postVar[i,i,k])` (numeric). **Only present when `condVar = TRUE` was used in the preceding `ranef()` call.** |

Key behaviors:
- All grouping factors are stacked into a single data frame (rows for all
  `grpvar × term × grp` combinations).
- When `condVar` was `FALSE`, the `condsd` column is absent — the function
  gracefully handles missing `"postVar"` attribute.
- Internally uses the `lme4:::asDf0` helper which also supports a `"condVar"`
  attribute name alias (for forward compatibility with any eventual rename).

---

## 3. `print.ranef.mer`

```r
print(x, ...)
```

Prints the unclassed list of data frames (one per grouping factor) and,
when any element carries a non-null `"postVar"` attribute, appends the
line:

```
with conditional variances for "Subject"
```

This is the only textual indication to the user that uncertainty estimates
are attached — it serves as a reminder that `as.data.frame()` or `dotplot()`
will have `condsd` / error bars available.

---

## 4. `dotplot.ranef.mer` — caterpillar plot (lattice)

```r
## S3 method for class 'ranef.mer'
dotplot(x, data, main = TRUE, transf = I, level = 0.95, ...)
```

Requires the **lattice** package to be loaded (the generic lives in lattice).

| Argument | Description |
|---|---|
| `x` | A `ranef.mer` object (from `ranef(condVar=TRUE)`). |
| `data` | Required by the generic but not used; pass nothing or `NULL`. |
| `main` | Logical; whether to add the grouping factor name as plot title. |
| `transf` | A transformation function applied to the RE values and their error bars before plotting. Default `I` (identity). Useful for e.g. `exp` on logistic regression models to display odds ratios. |
| `level` | Confidence level for the error bars. Default `0.95`. Translates to `±qnorm((1+level)/2) × condsd`. |
| `...` | Passed through to the lattice `dotplot()` call. Supports all lattice panel arguments including `lty`, `lwd`, `col`, `lty.v`, `lwd.v`, `col.line.v` (vertical reference line at zero), `lty.h`, `lwd.h`, `col.line.h` (horizontal grid lines), and `scales` (e.g. `scales=list(x=list(relation="free"))` for independent axes across panels). |

**Return value**: a **named list of `trellis` objects**, one per grouping
factor.  The list is named by the grouping factor names (e.g.
`list(Subject = <trellis>)`).  Individual panels within each trellis are
the RE terms (e.g. `(Intercept)` and `Days`).  RE levels are ordered by
their intercept value on the y-axis (caterpillar ordering).

**Behaviour when `postVar` is absent**: the error-bar segments are simply
not drawn (the `se` variable is `NULL`), but the point estimates are
still plotted.

---

## 5. `qqmath.ranef.mer` — Q-Q plot against standard normal (lattice)

```r
## S3 method for class 'ranef.mer'
qqmath(x, data, main = TRUE, level = 0.95, ...)
```

| Argument | Description |
|---|---|
| `x` | A `ranef.mer` object. |
| `data` | Required by generic; not used. |
| `main` | Add grouping factor name as title. |
| `level` | Confidence level for error bars (same `±qnorm` construction as `dotplot`). |
| `...` | Lattice arguments passed through. |

**Behaviour when `postVar` is present**: ranks each level within each RE
term against `qnorm((rank − 0.5) / n)` and draws ±1.96 (at 0.95 level)
horizontal error bar segments. Panels are faceted by RE term with free
x-scales (`scales = list(x = list(relation = "free"))`).

**Behaviour when `postVar` is absent**: falls back to `qqmath(~values|ind,
data = stack(xt))` without error bars and with free y-scales.

**Return value**: a **named list of `trellis` objects**, one per grouping
factor (same structure as `dotplot`).

---

## 6. `plot.ranef.mer` — generic S3 plot dispatch (lattice)

```r
## S3 method for class 'ranef.mer'
plot(x, y, ...)
```

Dispatches to different lattice plot types depending on the number of
RE columns per grouping factor:

| `ncol(x[[i]])` | Plot type |
|---|---|
| 1 | `qqmath(~col, x[[i]], ...)` — Q-Q plot of the single RE |
| 2 | `xyplot(col1 ~ col2, x[[i]], ...)` — scatter of two RE terms |
| ≥ 3 | `splom(~x[[i]], ...)` — scatterplot matrix |

Does **not** use the `"postVar"` attribute; no error bars.  Primarily a
quick diagnostic tool.

**Return value**: a list of `trellis` objects, one per grouping factor.

---

## 7. Summary of the full surface

| Function / method | Key args | Needs condVar/postVar? | What user gets |
|---|---|---|---|
| `ranef(object, condVar=TRUE)` | `condVar`, `drop`, `whichel`, `postVar` (deprecated) | produces it | `ranef.mer` list of data frames with `"postVar"` attr |
| `ranef(object, condVar=FALSE)` | — | no | `ranef.mer` list without `"postVar"` attr |
| `attr(re[[g]], "postVar")` | — | yes | `j×j×k` array of conditional cov matrices |
| `as.data.frame(re)` | — | optional | long data frame; `condsd` col present iff `condVar=TRUE` |
| `print(re)` | — | optional | prints data frames + "with conditional variances" line |
| `dotplot(re, transf, level, ...)` | `transf`, `level`, lattice `...` | yes (bars omitted otherwise) | named list of `trellis` caterpillar plots |
| `qqmath(re, level, ...)` | `level`, lattice `...` | yes (falls back otherwise) | named list of `trellis` Q-Q plots with CI bars |
| `plot(re, ...)` | lattice `...` | no | list of `trellis` scatter/QQ/splom plots |

### Argument defaults (lme4 vs common user expectation)

- lme4 defaults `condVar = TRUE` in `ranef.merMod`.  This is intentional:
  the common use case (plotting caterpillars) needs the postVar, so lme4
  computes it unless explicitly disabled.
- The `"postVar"` attribute name is the stored name on the data frame, even
  though the argument is called `condVar`.  `asDf0` also checks for a
  `"condVar"` attribute name (forward-compat alias).

### Typical ggplot2 workflow enabled by `as.data.frame`

```r
dd <- as.data.frame(ranef(fit, condVar = TRUE))
ggplot(dd, aes(y = grp, x = condval)) +
  geom_point() +
  facet_wrap(~term, scales = "free_x") +
  geom_errorbarh(aes(xmin = condval - 2*condsd,
                     xmax = condval + 2*condsd), height = 0)
```

This is the canonical lme4 vignette pattern; `condsd` must be present in
the data frame for error bars to render.
