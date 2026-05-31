# lme4 / lmerTest Public Surface: Matrix & Structure Accessors

**Survey date:** 2026-05-31  
**lme4 version:** 2.0.1  
**lmerTest version:** 3.2.1 (no additions to this family)  
**Purpose:** Exhaustive reference against which mixeff is judged. Does NOT assess mixeff — that is a separate step.

---

## 1. `getME(object, name, ...)`

**Package:** lme4  
**Class:** `merMod` (S3 method `getME.merMod`)  
**Source doc:** `?getME`

The single omnibus extractor. Pass `name = "ALL"` to return every component as a named list. Accepts a character vector of names; returns a named list when `length(name) > 1`.

### Complete name catalogue

| Name | Return type | What it is / why a user needs it |
|------|-------------|----------------------------------|
| `"X"` | `matrix` (dense) | Fixed-effects design matrix (N × p). Needed for manual SE computation, leverage, projection. |
| `"Z"` | `dgCMatrix` (sparse) | Random-effects design matrix (N × q). Full stacked matrix across all grouping factors. |
| `"Zt"` | `dgCMatrix` | Transpose of Z (q × N). lme4 stores Z' internally; the direct accessor avoids an extra transpose. |
| `"Ztlist"` | named list of `dgCMatrix` | Z' split per random-effects term. Needed when per-term block structure matters (e.g. custom variance diagnostics). |
| `"mmList"` | list of matrices | Raw (un-transposed, un-scaled) model matrices for each RE term before Lambda scaling. |
| `"y"` | numeric vector (n) | Response vector as seen by the optimizer. |
| `"mu"` | numeric vector (N) | Conditional mean of the response (= fitted values for LMMs). |
| `"u"` | numeric vector (q) | Conditional mode of the *spherical* random effects U (before Lambda scaling). |
| `"b"` | numeric vector (q) | Conditional mode of the actual random effects b = Lambda · u. |
| `"Gp"` | integer vector (k+1) | Groups pointer: `Gp[i]` is the first index in q-space of the i-th RE term; `Gp[k+1] = q`. Used to slice Z/Lambda/u by term. |
| `"Tp"` | integer vector (k+1) | Theta pointer: analogous to Gp but for the theta vector. Needed to map theta elements back to terms. |
| `"L"` | `dCHMsimpl` (sparse Cholesky) | Cholesky factor of the penalized RE system. Used for computing leverage, log-determinant, and solving the mixed-model equations. |
| `"Lambda"` | `dgCMatrix` | Relative covariance factor: the lower-triangular block-diagonal matrix such that Var(b) = σ² Λ Λ'. |
| `"Lambdat"` | `dgCMatrix` | Transpose of Lambda (= Λ'). lme4 stores this directly. |
| `"Lind"` | integer vector | Index mapping: position i in `nnz(Lambda)` is filled by `theta[Lind[i]]`. Needed for understanding/manipulating the Cholesky parameterisation. |
| `"Tlist"` | list of matrices | Template matrices from which Lambda blocks are generated; one per RE term. |
| `"A"` | `dgCMatrix` | Scaled sparse model matrix A = Zt %*% Lambdat (q × N). Intermediate in the PIRLS/REML computation. |
| `"RX"` | dense matrix (p × p) | Cholesky factor of X'X in the profiled system. Used to compute Var(β̂) = σ² (RX'RX)⁻¹ without re-inverting the full system. |
| `"RZX"` | dense matrix (q × p) | Cross-term of the full Cholesky: the off-diagonal block connecting RE and FE in the joint system. |
| `"sigma"` | scalar | Residual standard deviation. (Prefer `sigma(object)` over `getME(., "sigma")`.) |
| `"flist"` | named list of factors | The grouping variables (factors) for each RE term. Used to reconstruct level sets, sort random effects, or build design matrices. |
| `"fixef"` | named numeric vector | Fixed-effects estimates (same as `fixef(object)` but without names in some contexts). |
| `"beta"` | numeric vector | Same as fixef but unnamed. |
| `"theta"` | numeric vector | Lower-triangular Cholesky entries of relative RE covariance matrices (column-major). The optimization parameter vector. |
| `"ST"` | list of matrices | S and T factors in the TSST' factorisation of relative variance matrices; diagonal = S entries, off-diagonal = T entries. Alternative view of theta. |
| `"par"` | numeric vector | Concatenation of theta (and beta for GLMMs with nAGQ > 0). The full deviance-function argument. |
| `"REML"` | integer | 0 = ML-fitted; positive integer = REML-fitted. |
| `"is_REML"` | logical | Same result as `isREML(object)`. |
| `"n_rtrms"` | integer | Number of random-effects *terms* (k). |
| `"n_rfacs"` | integer | Number of distinct grouping *factors*. |
| `"N"` | integer | Rows of X (= n for LMMs; may differ for offset models). |
| `"n"` | integer | Length of response vector y. |
| `"p"` | integer | Columns of X (number of fixed-effect parameters). |
| `"q"` | integer | Columns of Z (total number of RE coefficients). |
| `"p_i"` | integer vector | Number of columns of each raw mmList matrix (columns per RE term before grouping). |
| `"l_i"` | integer vector | Number of levels of each grouping factor. |
| `"q_i"` | integer vector | Number of columns of each term-wise Z block (= p_i × l_i). |
| `"k"` | integer | Number of random-effects terms (= n_rtrms). |
| `"m_i"` | integer vector | Number of covariance parameters in each RE term. |
| `"m"` | integer | Total number of covariance parameters (= length of theta). |
| `"cnms"` | named list | "Component names": maps grouping-factor name → character vector of slope names for that factor's RE block. |
| `"devcomp"` | list: `$cmp` (numeric), `$dims` (integer) | Full deviance decomposition and dimension vector from the optimizer. Contains `ldL2`, `ldRX2`, `wrss`, `ussq`, `pwrss`, `REML`, `dev`, `sigmaML`, `sigmaREML`; and dims `N`, `n`, `p`, `q`, `nth`, `nAGQ`, `reTrms`, `REML`, `GLMM`, `NLMM`. |
| `"offset"` | numeric vector | Model offset (length 0 when none specified). |
| `"lower"` | numeric vector | Lower bounds on theta: 0 for diagonal Cholesky entries, -Inf for off-diagonal. Used by `isSingular()` to detect boundary fits. |
| `"devfun"` | function | The deviance/REML criterion function as a closure. Used for profiling, bootstrap, and custom optimisation. |
| `"devarg"` | numeric vector | Current argument to devfun (= par). |
| `"glmer.nb.theta"` | scalar | Negative-binomial dispersion parameter; only for `glmer.nb()` fits. |

---

## 2. `model.matrix(object, type, noScale, ...)`

**Method:** `model.matrix.merMod`  
**Key argument:** `type = c("fixed", "random", "randomListRaw")`

| `type` value | Return | Use case |
|---|---|---|
| `"fixed"` (default) | Dense `matrix` (N × p) | Standard fixed-effect design matrix; same as `getME(., "X")`. |
| `"random"` | `dgCMatrix` (N × q) | Full stacked sparse random-effect design matrix; same as `getME(., "Z")`. |
| `"randomListRaw"` | Named list of matrices | Per-term raw (un-Lambda-scaled) RE design matrices; same as `getME(., "mmList")`. |

`noScale` (logical): when autoscaling was applied, `noScale = TRUE` returns parameters on the original scale.

---

## 3. `model.frame(formula, fixed.only, ...)`

**Method:** `model.frame.merMod`  
**Key argument:** `fixed.only = FALSE`

Returns the stored `frame` slot — the `data.frame` containing every variable required to parse the model formula. When `fixed.only = FALSE` (default), grouping factors and other random-effect variables are included. Users need this to reconstruct covariates, check NA handling, or feed downstream tools.

---

## 4. `terms(x, fixed.only, random.only, ...)`

**Method:** `terms.merMod`  
**Key arguments:** `fixed.only = TRUE` (default), `random.only = FALSE`

Returns the `terms` object describing the model structure. `fixed.only = TRUE` gives only fixed-effect terms (the lm-style terms object, needed for `model.matrix` rebuilding, emmeans, etc.). `fixed.only = FALSE` gives all terms including random-effect syntax. Mutually exclusive with `random.only = TRUE`.

---

## 5. `formula(x, fixed.only, random.only, ...)`

**Method:** `formula.merMod`

Returns the model formula. `fixed.only = FALSE` (default) returns the full formula including RE syntax. `fixed.only = TRUE` strips RE terms. `random.only = TRUE` returns only the RE portion. Users need the full formula for `update()`, refitting, and documentation.

---

## 6. `fixef(object, add.dropped, noScale, ...)`

**Method:** `fixef.merMod`  
**Key arguments:**
- `add.dropped = FALSE`: if `TRUE`, re-inserts `NA` placeholders for rank-dropped columns so the result always has length = ncol(full design matrix).
- `noScale`: return parameters on original (unscaled) scale.

Returns a **named numeric vector** of fixed-effect estimates. Identical to `getME(., "beta")` but with names. Users rely on names for downstream use in emmeans, contrasts, reporting.

---

## 7. `ranef(object, condVar, drop, whichel, postVar, ...)`

**Method:** `ranef.merMod`  
**Key arguments:**
- `condVar = TRUE`: attach the conditional variance-covariance array as attribute `"postVar"` on each data frame. This is the primary mechanism for caterpillar plots and random-effect uncertainty intervals.
- `drop = FALSE`: if `TRUE`, single-column data frames become named vectors.
- `whichel`: character vector of grouping factor names to include.
- `postVar`: deprecated synonym for `condVar`.

Returns an object of class `"ranef.mer"` (list of data frames, one per grouping factor). The `"postVar"` attribute is a **3-D array** (p × p × n_levels) of conditional variance-covariance matrices. Supports `as.data.frame()` (long format with `grpvar`, `term`, `grp`, `condval`, `condsd` columns), `dotplot()`, and `qqmath()` for visualisation.

---

## 8. `VarCorr(x, sigma, ...)`

**Method:** `VarCorr.merMod`  
**Key argument:** `sigma = 1` (multiplier for SDs)

Returns an object of class `"VarCorr.merMod"`: a list of variance-covariance matrices (one per grouping factor), each carrying attributes `"stddev"` (SDs) and `"correlation"` (correlation matrix), plus residual SD as attribute `"sc"`.

Additional attributes: `"theta"` (raw Cholesky parameters), `"rho"` (structured-covariance parameterisation), `"profpar"` (profiling parameterisation).

`as.data.frame()` method produces columns: `grp`, `var1`, `var2` (NA for variances), `vcov`, `sdcor`.

`print()` accepts `comp = c("Std.Dev.", "Variance")` to control display.

---

## 9. `vcov(object, correlation, sigm, use.hessian, full, noScale, ...)`

**Method:** `vcov.merMod`  
**Key arguments:**
- `correlation = TRUE`: attaches a correlation matrix as `attr(result, "factors")$correlation`.
- `use.hessian = NULL`: auto-selects finite-difference Hessian (exact for GLMMs with nAGQ > 0) vs. RX-based inversion (exact for LMMs, approximate for GLMMs). `FALSE` uses the RX-based method for backward compatibility.
- `full = FALSE`: if `TRUE`, returns the **joint** covariance matrix of both conditional modes (b) and fixed-effect parameters — a sparse matrix of dimension (q + p) × (q + p). Slow for large models.
- `noScale`: return on original scale.

Returns a **p × p** covariance matrix (or (q+p) × (q+p) sparse when `full = TRUE`). Row/column names are the fixed-effect parameter names.

---

## 10. `sigma(object, ...)`

**Method:** `sigma.merMod`

Returns the **residual standard deviation** (scalar). For Gaussian: √(residual variance). For Gamma/InvGauss: √(1/shape). For Binomial/Poisson: always 1 (no scale parameter). The canonical way to get σ; `getME(., "sigma")` is a secondary path.

---

## 11. `fitted(object, ...)`

**Method:** `fitted.merMod`

Returns a numeric vector of length N: fitted values conditional on the estimated random-effect modes. For GLMMs, these are on the response scale (not linear predictor). For more flexible access (marginal fits, new data), use `predict.merMod`.

---

## 12. `residuals(object, type, scaled, ...)`

**Method:** `residuals.merMod`  
**Key arguments:**
- `type`: `"response"` (default for LMMs) or `"deviance"` (default for GLMMs); also `"pearson"`, `"working"`. Note: `"partial"` is not yet implemented.
- `scaled = FALSE`: if `TRUE`, divides by residual SD.

Returns a numeric vector of length N. Pearson residuals differ from `residuals.lme` in that lme4 does not scale by the estimated SD unless `scaled = TRUE`.

---

## 13. `hatvalues(model, fullHatMatrix, ...)`

**Method:** `hatvalues.merMod`  
**Key argument:** `fullHatMatrix = FALSE`

Returns the **diagonal of the hat (projection) matrix** H as a numeric vector of length N. If `fullHatMatrix = TRUE`, returns the full N × N matrix H (memory-intensive). Note: meaningful only for LMMs; the hat-matrix concept is not well-defined for GLMMs.

---

## 14. `weights(object, type, ...)`

**Method:** `weights.merMod`  
**Key argument:** `type = c("prior", "working")`

- `"prior"`: the observation weights supplied by the user at fit time (all 1 if none supplied).
- `"working"`: the final PIRLS weights at convergence (relevant for GLMMs; equals prior weights for LMMs).

---

## 15. `ngrps(object, ...)`

**Method:** `ngrps.merMod`  
**Also:** `ngrps.factor`

Returns a **named integer vector**: one element per distinct grouping factor, giving the number of levels. Users need this to assess the effective sample size for random effects and check degrees of freedom.

---

## 16. `nobs(object, ...)`

**Method:** `nobs.merMod`

Returns the number of observations N (scalar). Needed for AIC/BIC penalty and degrees-of-freedom calculations.

---

## 17. `coef(object, ...)`

**Method:** `coef.merMod`

Returns a list of data frames of class `"coef.mer"`, one per grouping factor. Each row is a level of that factor; columns are fixed-effect names; values are the sum of fixed and random effects (β + b_i). Users use this for subject-specific (conditional) predictions and inspection of BLUPs on the fixed-effect scale.

---

## 18. `getData(object)` / `fortify.merMod(model, data, ...)`

**Function:** `getData.merMod` (S3 method)  
**Function:** `fortify.merMod` (plain function, not an S3 method — avoids ggplot2 dependency)

`getData()` retrieves the original `data` argument from the formula environment. Minimal implementation; does not handle `na.action` or `subset` specially.

`fortify.merMod()` augments the data frame with columns: `.fitted`, `.resid`, `.fixed` (fixed-part only), `.mu`, `.u`, `.offset`, `.sqrtXwt`, `.sqrtrwt`, `.weights`, `.wtres`. Used by ggplot2-based diagnostic plots; broom.mixed is now preferred.

---

## 19. `isSingular(object, tol)` / `getSingTol()`

**Function:** `isSingular` (generic, exported)

Tests whether any random-effects covariance matrix is (near-)singular. Internally checks whether `min(getME(object, "theta")[getME(object, "lower") == 0]) < tol`. Users need this to decide whether to simplify the RE structure, interpret estimates with caution, or reject the model.

`rePCA(object)` provides deeper singularity diagnostics: PCA of each RE covariance block, showing effective dimensionality and the mapping from variance components to orthogonal directions.

---

## 20. `isREML(object)` / `isLMM(object)` / `isGLMM(object)` / `isNLMM(object)`

Type-predicate functions. Return logical scalars. Used to branch on model class before applying class-specific extractors (e.g. `hatvalues` on LMMs only, `sigma = 1` for Binomial GLMMs).

---

## 21. `deviance(object, REML, ...)` / `REMLcrit(object)` / `logLik(object, REML, ...)`

These are standard S3 methods on `merMod`.

| Function | Returns | Notes |
|---|---|---|
| `deviance(object)` | scalar | −2·log-likelihood (absolute deviance for LMMs; conditional-relative for GLMMs). |
| `REMLcrit(object)` | scalar | REML criterion at optimum. Convenience; equivalent to `deviance(object, REML=TRUE)` for REML fits. |
| `logLik(object, REML=NULL)` | `logLik` object with df/nobs attrs | Log-likelihood. For GLMMs, proportional to the true log-density when using Laplace approximation. |

---

## 22. `devcomp` (via `getME`)

Provides a comprehensive named numeric vector (`$cmp`) and integer dimension vector (`$dims`) summarising the fitted model state. Key `$cmp` entries:

| Name | Meaning |
|---|---|
| `ldL2` | 2 × log|L| (log-det of Cholesky) |
| `ldRX2` | 2 × log|RX| |
| `wrss` | Weighted residual sum of squares |
| `ussq` | ‖u‖² |
| `pwrss` | Penalized WRSS = wrss + ussq |
| `REML` | REML criterion (REML fits only) |
| `dev` | Deviance (ML fits only) |
| `sigmaML` | σ̂ (ML) |
| `sigmaREML` | σ̂ (REML) |

Key `$dims` entries: `N`, `n`, `p`, `q`, `nth` (= m = length of theta), `nAGQ`, `reTrms`, `REML`, `GLMM`, `NLMM`, `useSc`.

---

## 23. Summary of the full surface

| Function/method | Signature highlights | Returns |
|---|---|---|
| `getME(object, name)` | 47 named components + "ALL" | component or named list |
| `model.matrix(object, type, noScale)` | type ∈ {fixed, random, randomListRaw} | matrix / dgCMatrix / list |
| `model.frame(formula, fixed.only)` | fixed.only logical | data.frame |
| `terms(x, fixed.only, random.only)` | two logicals | terms object |
| `formula(x, fixed.only, random.only)` | two logicals | formula |
| `fixef(object, add.dropped, noScale)` | add.dropped for rank-deficient | named numeric vector |
| `ranef(object, condVar, drop, whichel)` | condVar for postVar arrays | ranef.mer list of data.frames |
| `VarCorr(x, sigma)` | sigma multiplier; as.data.frame method | VarCorr.merMod list |
| `vcov(object, correlation, use.hessian, full, noScale)` | full for joint (q+p)×(q+p) matrix | matrix / sparseMatrix |
| `sigma(object)` | — | scalar |
| `fitted(object)` | — | numeric vector (N) |
| `residuals(object, type, scaled)` | 5 types; scaled flag | numeric vector (N) |
| `hatvalues(model, fullHatMatrix)` | optionally returns full H | numeric vector (N) or matrix |
| `weights(object, type)` | prior / working | numeric vector |
| `ngrps(object)` | — | named integer vector |
| `nobs(object)` | — | scalar integer |
| `coef(object)` | — | coef.mer list of data.frames |
| `getData(object)` | — | data.frame |
| `fortify.merMod(model, data)` | augments data with fit columns | data.frame |
| `isSingular(object, tol)` | tol default via getSingTol() | logical |
| `isREML / isLMM / isGLMM / isNLMM` | — | logical |
| `deviance / REMLcrit / logLik` | REML argument on deviance/logLik | scalar / logLik |
| `getME(., "devcomp")` | — | list($cmp, $dims) |

---

## lmerTest additions (this family)

lmerTest 3.2.1 does not add new matrix or structure accessors. It overloads `summary()` to attach Satterthwaite df and p-values to the fixed-effects coefficient table, and provides `lmerTest::contest()` / `lmerTest::ranova()` for inference, but these are inference functions, not matrix/structure accessors.

---

*End of survey.*
