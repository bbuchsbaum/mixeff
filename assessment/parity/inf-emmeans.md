# Parity Assessment: inf-emmeans

**Cell:** inf-emmeans  
**Dataset:** cake (lme4 built-in, n=270)  
**Formula:** `angle ~ recipe*temperature + (1|recipe:replicate)`  
**Focus:** `emmeans()` + `pairs()` vs lme4/lmerTest — estimates, SE, df, p  
**Date:** 2026-05-31  
**Probe script:** `assessment/parity/inf-emmeans-probe.R`

---

## Environment

| Package   | Version |
|-----------|---------|
| lme4      | 2.0.1   |
| lmerTest  | 3.2.1   |
| emmeans   | 2.0.3   |
| mixeff    | 0.1.0   |

---

## Raw output (abbreviated)

### lme4 fit
- Wall-clock: 0.032 s
- Converged (no warnings)

### mixeff fit
- Wall-clock: 0.013 s
- `fit_status: converged_interior`

### lme4 `emmeans(~ recipe)` (marginal, averaged over temperature)
```
 recipe emmean   SE df lower.CL upper.CL
 A        33.1 1.74 42     29.6     36.6
 B        31.6 1.74 42     28.1     35.1
 C        31.6 1.74 42     28.1     35.1
Degrees-of-freedom method: kenward-roger
```

### mixeff `emmeans(~ recipe)`
```
 recipe emmean   SE   df lower.CL upper.CL
 A        29.1 2.04 77.4     25.1     33.2
 B        26.9 2.04 77.4     22.8     30.9
 C        27.9 2.04 77.4     23.9     32.0
Degrees-of-freedom method: mixeff auto
```

**Observed:** mixeff returns (29.1, 26.9, 27.9) instead of (33.1, 31.6, 31.6).  
**Diagnosis:** (29.1, 26.9, 27.9) exactly matches the lme4 `recipe*temperature` cell means at temperature=175 (the first temperature level). The emmeans bridge is not averaging over temperature levels — it is evaluating only the first level of the temperature grid, i.e., the reference grid is not being populated with all temperature levels.

### lme4 `emmeans(~ temperature)` (marginal, averaged over recipe)
```
 temperature emmean   SE   df
 175           28.0 1.18 77.4
 185           30.0 1.18 77.4
 195           31.4 1.18 77.4
 205           32.2 1.18 77.4
 215           35.8 1.18 77.4
 225           35.4 1.18 77.4
```

### mixeff `emmeans(~ temperature)`
```
 temperature emmean   SE  df
 175           26.3 1.80 219
 185           25.3 1.86 227
 195           28.6 1.98 238
 205           33.6 1.74 211
 215           19.3 2.07 244
 225           34.8 1.39 133
```

**Diagnosis:** Wildly wrong estimates (max abs diff ~16.5 on temperature=215) and incorrect SE/df. These resemble single-cell (one recipe level) predictions rather than averages over recipe — consistent with the same reference-grid population bug.

### lme4 `pairs(emmeans(~ recipe))`
```
 A - B:  est=1.478  SE=2.46  df=42  p=0.820
 A - C:  est=1.522  SE=2.46  df=42  p=0.810
 B - C:  est=0.044  SE=2.46  df=42  p=1.000
```

### mixeff `pairs(emmeans(~ recipe))`
```
 A - B:  est=2.267  SE=2.882  df=77.4  p=0.713
 A - C:  est=1.200  SE=2.882  df=77.4  p=0.909
 B - C:  est=-1.067  SE=2.882  df=77.4  p=0.927
```

**Diagnosis:** Contrast estimates are wrong (derived from incorrect marginal means); SE and df are also wrong. p-values are qualitatively similar (all non-significant) but quantitatively divergent.

### lme4 `pairs(emmeans(~ temperature))`
```
 175-185:  est=-1.978  SE=0.954  df=210  p=0.305
 175-195:  est=-3.444  SE=0.954  df=210  p=0.005
 ...
```

### mixeff `pairs(emmeans(~ temperature))`
```
 175-185:  est=1.049  SE=2.025  df=210  p=0.605
 175-195:  est=-2.288  SE=2.234  df=210  p=0.307
 ...
```

df matches (210 vs 210, within-tol), but estimates and SE diverge severely (max abs diff ~17.9 on estimates, ~1.42 on SE). Sign reversals occur (e.g. 175-185: lme4=-1.978 vs mixeff=+1.049).

---

## Numerical comparison table

| Quantity | Max Abs Diff | Tolerance | Status |
|---|---|---|---|
| recipe emmean: estimate | 4.78 | 1e-4 | **EXCEEDS-TOL** |
| recipe emmean: SE | 0.301 | 1e-4 | **EXCEEDS-TOL** |
| recipe emmean: df | 35.4 | 1.0 | **EXCEEDS-TOL** |
| temperature emmean: estimate | 16.51 | 1e-4 | **EXCEEDS-TOL** |
| temperature emmean: SE | 0.893 | 1e-4 | **EXCEEDS-TOL** |
| temperature emmean: df | 166.4 | 1.0 | **EXCEEDS-TOL** |
| recipe\*temp emmean: estimate | 18.42 | 1e-4 | **EXCEEDS-TOL** |
| recipe\*temp emmean: SE | 1.642 | 1e-4 | **EXCEEDS-TOL** |
| recipe\*temp emmean: df | 168.97 | 1.0 | **EXCEEDS-TOL** |
| pairs(recipe): estimate | 1.111 | 1e-4 | **EXCEEDS-TOL** |
| pairs(recipe): SE | 0.426 | 1e-4 | **EXCEEDS-TOL** |
| pairs(recipe): df | 35.4 | 1.0 | **EXCEEDS-TOL** |
| pairs(recipe): p.value | 0.108 | 1e-4 | **EXCEEDS-TOL** |
| pairs(temp): estimate | 17.91 | 1e-4 | **EXCEEDS-TOL** |
| pairs(temp): SE | 1.421 | 1e-4 | **EXCEEDS-TOL** |
| pairs(temp): **df** | 5.99e-4 | 1.0 | **WITHIN-TOL** |
| pairs(temp): p.value | 0.996 | 1e-4 | **EXCEEDS-TOL** |

Only one quantity passes tolerance: `pairs(~ temperature)` df (Satterthwaite df=210 matches for the contrasts, since both use the same residual df path for within-temperature contrasts).

---

## Root cause analysis

### Primary bug: reference grid population in `recover_data.mm_lmm`

The `emmeans` bridge constructs its reference grid from the call + model frame. When `emmeans` builds a grid for `~ recipe`, it must populate the grid with **all levels of nuisance predictors** (here: temperature) and then average over them. The fact that mixeff's `emmeans(~ recipe)` returns values equal to the lme4 cell means at temperature=175 (the first/reference level of `temperature`) shows that the reference grid contains only one temperature level (the reference level encoded as 0 in the contrast columns), rather than all six levels.

The most likely mechanism: `recover_data.mm_lmm` passes `frame = object$model_frame`, but `object$model_frame` may not carry the `xlev` attribute (factor levels) in the form emmeans requires to reconstruct all levels of `temperature`. Without the full factor level information, emmeans falls back to a single-row grid (the reference level), so all predictions are evaluated at temperature=175 rather than averaged over {175, 185, 195, 205, 215, 225}.

This is also consistent with the `~ recipe * temperature` interaction grid showing wildly varying SE and df per cell (ranging from df=87 to df=246 per cell), which would not occur in a properly balanced grid.

### Secondary issue: df method mismatch

- lme4 via emmeans defaults to **Kenward-Roger** (df=42 for recipe, df=77.4 for temperature).
- mixeff via emmeans uses **Satterthwaite** (df=77.4 for recipe marginal means — which would be reasonable for the temperature marginal means, not recipe).

This is a secondary concern compared to the reference grid bug, since incorrect estimates dominate. Once estimates are correct, the df comparison can be re-evaluated.

### Native API (`mm_means` / `mm_comparisons`) also divergent

`mm_means(fit_mm, ~ recipe)` returns (29.13, 26.87, 27.93) — the same wrong values as the emmeans bridge — with df=77.4. This confirms the bug is not in the emmeans bridge itself but upstream in the marginal means computation. The native mm_means API has the same incorrect reference-grid averaging behavior.

`mm_means(fit_mm, ~ temperature)` returns (26.3, 25.3, 28.6, 33.6, 19.3, 34.8) — very different from lme4's (28.0, 30.0, 31.4, 32.2, 35.8, 35.4). Temperature=215 is especially egregious (19.3 vs 35.8, diff=16.5).

`mm_comparisons(fit_mm, ~ recipe)` returns contrasts (−2.27, −1.20, +1.07) derived from the wrong means.

---

## mm_means reliability flag

Both `mm_means` outputs carry `reliability = "low"` — mixeff is correctly flagging that the covariance source (`fixed_co....`) has limited reliability. This is an honest diagnostic, but the larger problem is that the estimates themselves are wrong, not just their uncertainty quantification.

---

## Speed

| Operation | lme4 mean/call | mixeff mean/call | Ratio (mm/lme4) |
|---|---|---|---|
| `emmeans(fit, ~ recipe)` | 0.032 s | 0.016 s | 0.50x (2x faster) |
| Model fit | 0.032 s | 0.013 s | ~0.4x (2.5x faster) |

Speed is not the concern here.

---

## Classification

**Outcome: divergent** — the emmeans bridge and native mm_means/mm_comparisons API produce systematically wrong marginal means for models with interaction terms. The reference-grid averaging over nuisance predictors (non-focal factors in the `specs`) is broken. This is an **in-scope-missing / partial** capability (the API exists and runs without error, but produces wrong numerical results).

**Severity: blocker** — emmeans support is advertised and the bridge runs without error, but all marginal means and contrasts for interaction models are wrong. A user relying on `emmeans(fit_mm, ~ recipe)` for a `recipe*temperature` model would obtain values corresponding only to the reference temperature level rather than population-averaged means, with no warning that the averaging failed. The `reliability = "low"` flag is present but does not indicate the estimates themselves are wrong. This silently produces incorrect scientific conclusions.

---

## What passes

- `pairs(~ temperature)` **df**: WITHIN-TOL (both give df=210 for within-temperature contrasts)
- Model fit itself converges correctly (fixef, theta, sigma all within tolerance per lmm-cake-interaction probe)
- mixeff is faster (2x on emmeans calls)
- The emmeans bridge does not error; the `mm_means`/`mm_comparisons` API does not error

## What fails

- All emmeans marginal means for any marginal (recipe, temperature, recipe*temperature) — estimates, SE, df
- All pairwise contrasts from those means — estimates, SE, p-values
- Native `mm_means` / `mm_comparisons` show the same wrong estimates (bug is upstream of the emmeans bridge)
- df method: mixeff uses Satterthwaite; lme4/emmeans defaults to Kenward-Roger (minor secondary issue)
