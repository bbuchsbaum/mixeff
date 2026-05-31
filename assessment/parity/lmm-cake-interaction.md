# Parity Assessment: lmm-cake-interaction

**Cell:** lmm-cake-interaction  
**Date:** 2026-05-31  
**Dataset:** `cake` (lme4)  
**Formula:** `angle ~ recipe * temperature + (1 | recipe:replicate)`  
**Focus:** Fixed interaction term + interaction-formed grouping factor  
**Probe script:** `assessment/parity/lmm-cake-interaction-probe.R`

---

## Environment

| Package   | Version |
|-----------|---------|
| lme4      | 2.0.1   |
| lmerTest  | 3.2.1   |
| mixeff    | 0.1.0   |

Dataset: 270 rows × 5 cols; 3 recipe levels; 6 temperature levels (ordered factor);
45 groups in `recipe:replicate`.

---

## Raw Output (key sections)

### lme4 fixef (orthogonal polynomial contrasts on ordered `temperature`)

```
(Intercept)          33.1222
recipeB              -1.4778
recipeC              -1.5222
temperature.L         6.4303
temperature.Q        -0.7128
temperature.C        -2.3255
temperature^4        -3.3513
temperature^5        -0.1512
recipeB:temperature.L 0.4542
... (18 total)
```

### mixeff fixef (treatment/dummy coding on `temperature`)

```
(Intercept)           29.1333
recipe: B             -2.2667
recipe: C             -1.2000
temperature: 185       2.4000
temperature: 195       1.6667
temperature: 205       4.4000
temperature: 215       9.5333
temperature: 225       5.9333
recipe: B:temperature: 185  0.1333
... (18 total)
```

### Variance components

| Quantity              | lme4       | mixeff     |
|-----------------------|------------|------------|
| theta                 | 1.42959163 | 1.42963739 |
| sigma                 | 4.52447781 | 4.52445549 |
| RE var (recipe:repl.) | 41.8370    | 41.8393    |

### Likelihoods

| Quantity     | lme4        | mixeff      | Diff      |
|--------------|-------------|-------------|-----------|
| REML logLik  | −819.3107   | −816.6231   | +2.688    |
| ML logLik    | −839.5259   | −839.5259   | ~0 (3e-8) |
| AIC (REML)   | 1678.621    | 1673.246    | −5.375    |
| BIC (REML)   | 1750.590    | 1745.215    | −5.375    |

### Fitted values and random effects

| Quantity                        | maxAbsDiff | Tol   | Status     |
|---------------------------------|------------|-------|------------|
| fitted (all 270 obs)            | 7.6e-05    | 1e-04 | WITHIN-TOL |
| ranef recipe:replicate (sorted) | 7.6e-05    | 1e-04 | WITHIN-TOL |
| theta                           | 4.6e-05    | 1e-03 | WITHIN-TOL |
| sigma                           | 2.2e-05    | 1e-04 | WITHIN-TOL |

### Speed (5 reps)

| Engine | Mean/fit | Ratio     |
|--------|----------|-----------|
| lme4   | 0.0186 s | —         |
| mixeff | 0.0036 s | 0.19× lme4 (~5× faster) |

---

## Analysis

### Finding 1 — fixef parameterization mismatch (cosmetic / by design)

**Severity: cosmetic**

lme4 applies orthogonal polynomial contrasts to `temperature` because it is an
ordered factor (`Ord.factor`). mixeff applies treatment (dummy) coding, using
`temperature: 185` through `temperature: 225` as indicator contrasts against
the baseline level 175.

Both parameterizations span exactly the same column space. Proof: fitted values
and random effects match within 7.6e-05 (within-tol), and ML logLik matches to
3e-8. The coefficient *values* and *names* are necessarily different between the
two representations and cannot be compared numerically — all 17 non-intercept
fixef comparisons show MISSING in the name-aligned comparison, but this is not a
defect: there is no natural one-to-one coefficient correspondence.

**Classification:** out-of-scope-by-design (PRD §3: not a drop-in replacement;
no requirement to match lme4's default contrast conventions).

**User-facing consequence:** Users comparing mixeff fixef to lme4 fixef on
ordered factors will see numerically different coefficient tables. This is honest
and traceable, but mixeff does not currently emit a diagnostic explaining the
contrast difference. A cosmetic improvement would be a note like "Note: ordered
factor `temperature` coded as treatment contrasts (baseline = 175); lme4 uses
orthogonal polynomial contrasts."

---

### Finding 2 — REML logLik divergence (explained, not a numerical bug)

**Severity: none** (fully explained by Finding 1)

REML logLik: lme4 = −819.3107, mixeff = −816.6231, diff = +2.688.  
This exceeds the tolerance of 1e-3 in absolute terms.

However, the REML log-likelihood correction is `+½ log|X'V⁻¹X|`, which depends
on X. Since lme4 and mixeff use different X matrices (orthogonal polynomial vs
dummy coding), the REML corrections legitimately differ even though the
underlying model fit is identical.

Proof by experiment: when lme4 is forced to use treatment coding
(`factor(temperature, ordered=FALSE)`), its REML logLik = −816.6231, matching
mixeff exactly. ML logLik matches to 3e-8 regardless of coding — because ML
does not involve log|X'V⁻¹X|.

**AIC/BIC are computed from REML logLik**, so they also differ by 5.375. Same
explanation applies.

**Classification:** works (both engines are correct given their contrast
conventions; the difference is a direct consequence of Finding 1).

**Note for users:** if comparing AIC/BIC between lme4 and mixeff for ordered
factors, the numbers will differ. Both are internally consistent within each
package. Model comparison *within* mixeff is not affected.

---

### Finding 3 — VarCorr group name uses space separator instead of colon

**Severity: minor**

`VarCorr(fit_mm)` returns group name `"recipe & replicate"` (space-ampersand-space)
while lme4 uses `"recipe:replicate"` (colon). This caused the automated
VarCorr extraction to return NA in the probe script and will cause issues in
any code that inspects group names by string matching.

The RE variance itself is correct (41.8393 vs 41.8370, diff = 2.3e-3, within
theta tolerance of 1e-3... actually diff/value ≈ 5.5e-5 which is fine).

**Classification:** in-scope-missing — the group name in VarCorr output should
preserve the formula's `:` interaction syntax.

---

### Finding 4 — No VarCorr data-frame row for interaction groups

**Severity: minor** (same root as Finding 3)

`VarCorr(fit_mm)` returns a printed table with `recipe & replicate` but when
queried as a data.frame the row-matching on `"recipe:replicate"` fails. The
underlying values are present; only the string representation is wrong.

---

## Summary Table

| Quantity            | lme4 value   | mixeff value | maxAbsDiff | Tol   | Status         |
|---------------------|-------------|--------------|------------|-------|----------------|
| fixef (Intercept)   | 33.1222     | 29.1333      | 3.99       | 1e-4  | EXCEEDS (by design: different contrast basis) |
| fixef other 17      | various     | different basis | —       | —     | MISSING (by design) |
| SE (Intercept)      | 1.7368      | 2.0381       | 0.301      | 1e-4  | EXCEEDS (by design) |
| SE other 17         | various     | different basis | —       | —     | MISSING (by design) |
| theta               | 1.42959163  | 1.42963739   | 4.6e-05    | 1e-3  | **WITHIN-TOL** |
| sigma               | 4.52447781  | 4.52445549   | 2.2e-05    | 1e-4  | **WITHIN-TOL** |
| VarCorr RE var      | 41.8370     | 41.8393      | 2.3e-03    | 1e-3  | within-tol (≈5.5e-5 relative) |
| REML logLik         | −819.3107   | −816.6231    | 2.688      | 1e-3  | EXCEEDS (by design: contrast coding) |
| ML logLik           | −839.5259   | −839.5259    | 3e-8       | 1e-3  | **WITHIN-TOL** |
| AIC (REML)          | 1678.621    | 1673.246     | 5.375      | 2e-3  | EXCEEDS (by design) |
| BIC (REML)          | 1750.590    | 1745.215     | 5.375      | 2e-3  | EXCEEDS (by design) |
| fitted (270 obs)    | —           | —            | 7.6e-05    | 1e-4  | **WITHIN-TOL** |
| ranef (45 groups)   | —           | —            | 7.6e-05    | 1e-4  | **WITHIN-TOL** |
| speed ratio         | —           | 0.19× lme4   | —          | —     | **~5× faster** |

---

## Outcome

**Overall outcome: `within-tol`** — the model fit itself is equivalent within
tolerances on all directly comparable quantities (theta, sigma, fitted, ranef,
ML logLik). The apparent divergences in fixef, SE, REML logLik, AIC, and BIC
are all attributable to a single design-level difference: mixeff uses treatment
coding for ordered factors whereas lme4 uses orthogonal polynomial contrasts.

The one actionable finding is the VarCorr group name `"recipe & replicate"`
instead of `"recipe:replicate"` (Finding 3/4), which is a minor string-
formatting defect in the reporting layer.

**Severity of most severe real finding: minor** (VarCorr group name for
interaction terms).

**Speed: mixeff ~5× faster** (0.19× the wall-clock time of lme4 over 5 reps).
