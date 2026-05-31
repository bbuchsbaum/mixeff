# lme4/lmerTest Datasets & Utilities — Reference Survey

**Surveyed packages:** lme4 2.0.1, lmerTest 3.2.1, reformulas (re-exported shims)
**Survey date:** 2026-05-31
**Purpose:** Exhaustive reference of the public surface for the "Datasets & utilities"
capability family. mixeff parity is assessed separately; this document only records
what lme4/lmerTest offer and why users rely on each item.

---

## 1. Shipped Datasets

### 1.1 lme4 Datasets

#### `sleepstudy`
- **Format:** data.frame, 180 rows × 3 cols (`Reaction`, `Days`, `Subject`)
- **Variables:** `Reaction` (numeric, ms); `Days` (numeric, 0–9); `Subject` (factor, 18 levels)
- **Design:** Belenky et al. (2003) sleep deprivation study. 18 subjects measured over 10 days
  of progressive sleep restriction (3 h/night). Days 0–1 = adaptation/training, Day 2 =
  baseline, Days 3–9 = deprivation.
- **Typical use:** The canonical lme4 example for a random-slope LMM:
  `lmer(Reaction ~ Days + (Days|Subject), sleepstudy)`. Used in almost every lme4 vignette,
  tutorial, and test suite as the go-to continuous-predictor balanced design.

#### `cbpp`
- **Format:** data.frame, 56 rows × 4 cols (`herd`, `incidence`, `size`, `period`)
- **Variables:** `herd` (factor, 15 levels); `incidence` (count); `size` (denominator count);
  `period` (factor, 4 levels)
- **Design:** Contagious bovine pleuropneumonia in cattle herds across 4 survey periods.
  Binomial outcome (incidence/size per herd-period cell).
- **Typical use:** Standard GLMM binomial example:
  `glmer(cbind(incidence, size - incidence) ~ period + (1|herd), cbpp, binomial)`.
  Used to demonstrate aggregated binomial response and overdispersion diagnostics.

#### `cbpp2`
- **Format:** data.frame, 56 rows × 6 cols (`herd`, `treatment`, `avg_size`, `period`, `size`,
  `incidence`)
- **Design:** Extended version of `cbpp` adding a treatment factor and average size variable.
  Useful for demonstrating more complex GLMM designs with a treatment effect.

#### `Pastes`
- **Format:** data.frame, 60 rows × 4 cols (`strength`, `batch`, `cask`, `sample`)
- **Variables:** `strength` (numeric); `batch` (factor, 10 levels A–J); `cask` (factor, 3
  levels a–c); `sample` (factor, 30 levels, derived as batch:cask interaction)
- **Design:** Paste strength from batches subdivided into casks. Classic nested/hierarchical
  design: cask nested within batch.
- **Typical use:** Demonstrates nested random effects `(1|batch/cask)` or equivalently
  `(1|batch) + (1|sample)`. Used in the lme4 vignette for LMM anatomy; `isNested` example.

#### `Penicillin`
- **Format:** data.frame, 144 rows × 3 cols (`diameter`, `plate`, `sample`)
- **Variables:** `diameter` (numeric, zone of inhibition in mm); `plate` (factor, 24 levels);
  `sample` (factor, 6 levels)
- **Design:** Penicillin concentration assay using a crossed incomplete block design: 6 samples
  assayed on 24 plates.
- **Typical use:** Illustrates a fully crossed two-way random effects model:
  `lmer(diameter ~ (1|plate) + (1|sample), Penicillin)`. Used in lme4 vignette for balanced
  crossed designs.

#### `Dyestuff`
- **Format:** data.frame, 30 rows × 2 cols (`Batch`, `Yield`)
- **Variables:** `Batch` (factor, 6 levels A–F); `Yield` (numeric)
- **Design:** Dye yield from 6 batches, 5 observations per batch. Classic one-way random
  effects.
- **Typical use:** The simplest possible LMM example: `lmer(Yield ~ 1 + (1|Batch), Dyestuff)`.
  Used in lme4 vignette as the first-pass example; demonstrates variance components and REML
  vs ML.

#### `Dyestuff2`
- **Format:** data.frame, 30 rows × 2 cols (`Batch`, `Yield`)
- **Variables:** Same structure as `Dyestuff` but with smaller between-batch variability.
- **Design:** Companion dataset to `Dyestuff` where the estimated batch variance is near zero,
  demonstrating a boundary/singular fit.
- **Typical use:** Demonstrates `isSingular()`, boundary effects, and singular covariance
  matrix handling.

#### `cake`
- **Format:** data.frame, 270 rows × 5 cols (`replicate`, `recipe`, `temperature`, `angle`,
  `temp`)
- **Variables:** `replicate` (factor, 15 levels); `recipe` (factor, 3 levels A–C);
  `temperature` (ordered factor, 6 levels 175–225°C); `angle` (integer, breaking angle);
  `temp` (numeric, temperature as continuous)
- **Design:** Chocolate cake baking experiment. 15 replicates × 3 recipes × 6 temperatures
  = 270 observations. Response is breaking angle of cake.
- **Typical use:** Demonstrates LMM with both fixed factorial effects and a continuous
  covariate alongside random intercepts: `lmer(angle ~ recipe * temp + (1|replicate), cake)`.

#### `grouseticks`
- **Format:** data.frame, 403 rows × 7 cols (`INDEX`, `TICKS`, `BROOD`, `HEIGHT`, `YEAR`,
  `LOCATION`, `cHEIGHT`)
- **Variables:** `TICKS` (count, tick burden per grouse chick); `BROOD` (factor, 118 levels);
  `HEIGHT` (numeric, altitude); `YEAR` (factor, 95/96/97); `LOCATION` (factor, 63 levels);
  `cHEIGHT` (centered altitude)
- **Design:** Tick infestation of red grouse chicks from Scotland; overdispersed counts nested
  within broods and locations.
- **Typical use:** Poisson/negative binomial GLMM with nested random effects and an
  overdispersion random effect. `grouseticks_agg` is an aggregated companion (118 rows with
  brood-level means and variances).

#### `grouseticks_agg`
- **Format:** data.frame, 118 rows × 7 cols (`BROOD`, `meanTICKS`, `varTICKS`, `HEIGHT`,
  `YEAR`, `LOCATION`, `cHEIGHT`)
- **Typical use:** Aggregated version of `grouseticks` for brood-level analyses.

#### `VerbAgg`
- **Format:** data.frame, 7,584 rows × 9 cols (`Anger`, `Gender`, `item`, `resp`, `id`,
  `btype`, `situ`, `mode`, `r2`)
- **Variables:** `resp` (ordered factor 3 levels: no/perhaps/yes); `r2` (binary factor N/Y);
  `id` (factor, 316 subjects); `item` (factor, 24 items); `Anger` (numeric trait score);
  `Gender` (factor); `btype` (behavior type: curse/scold/shout); `situ` (other/self);
  `mode` (want/do)
- **Design:** Verbal aggression item-response study (Vansteelandt 2000). Each of 316 subjects
  responded to 24 scenarios; both subject and item are random effects. Cross-classified design.
- **Typical use:** Demonstrates cross-classified GLMMs with binary or ordinal outcomes:
  `glmer(r2 ~ (Anger + Gender + btype + situ)^2 + (1|id) + (1|item), VerbAgg, binomial)`.
  The largest standard lme4 example dataset; stresses large-scale GLMM fitting.

#### `InstEval`
- **Format:** data.frame, 73,421 rows × 7 cols (`s`, `d`, `studage`, `lectage`, `service`,
  `dept`, `y`)
- **Variables:** `s` (factor, 2,972 students); `d` (factor, 1,128 instructors); `y` (integer
  1–5, evaluation score); `studage`/`lectage` (ordered factors); `service`/`dept` (factors)
- **Design:** ETH Zürich student evaluations of lecturers. Cross-classified: student × instructor
  with large and unbalanced design.
- **Typical use:** Large-scale cross-classified LMM benchmark:
  `lmer(y ~ service + (1|s) + (1|d), InstEval)`. Used to demonstrate performance and
  scalability of lme4 with large n and high-dimensional random effects.

#### `Arabidopsis`
- **Format:** data.frame, 625 rows × 8 cols (`reg`, `popu`, `gen`, `rack`, `nutrient`, `amd`,
  `status`, `total.fruits`)
- **Variables:** `total.fruits` (integer count); `reg` (region factor, 3 levels); `popu`
  (population factor, 9 levels); `gen` (integer, genotype); `rack`/`nutrient`/`amd`/`status`
  (experimental factors)
- **Design:** Banta et al. (2010) Arabidopsis thaliana fruit count study; Poisson/NB GLMM
  with nested random effects (genotype within population within region).
- **Typical use:** Demonstrates negative binomial GLMM with complex nested random structure
  and zero-inflation considerations.

#### `salamander`
- **Format:** data.frame, 360 rows × 8 cols (`Season`, `Experiment`, `TypeM`, `TypeF`,
  `Cross`, `Male`, `Female`, `Mate`)
- **Variables:** `Mate` (binary 0/1); `TypeM`/`TypeF` (factor, R=Rough Butt/W=White Side);
  `Cross` (factor, RR/RW/WR/WW); `Male`/`Female` (integer IDs)
- **Design:** McCullagh & Nelder (1989) salamander mating experiment. Binary response with
  both male and female as crossed random effects.
- **Typical use:** Classic cross-classified GLMM with two random effects in a logistic model.

#### `schizophrenia`
- **Format:** data.frame, 1,603 rows × 9 cols (`id`, `imps79`, `imps79b`, `imps79o`, `int`,
  `TxDrug`, `Week`, `SqrtWeek`, `TxSWeek`)
- **Variables:** `imps79` (numeric, IMPS score); `imps79b` (binary); `imps79o` (ordinal);
  `TxDrug` (treatment indicator); `Week`/`SqrtWeek`/`TxSWeek` (time variables)
- **Design:** Davis (2002) schizophrenia clinical trial longitudinal data. Repeated measures
  on subjects; three response coding variants included.
- **Typical use:** Demonstrates longitudinal GLMM (binary/ordinal outcomes) with time ×
  treatment interactions.

#### `toenail`
- **Format:** data.frame, 1,908 rows × 5 cols (`patientID`, `outcome`, `treatment`, `time`,
  `visit`)
- **Variables:** `outcome` (factor 2 levels, none-or-mild/moderate-to-severe toenail
  separation); `treatment` (factor, itraconazole vs comparator); `time` (numeric); `visit`
  (integer)
- **Design:** De Backer et al. (1998) onychomycosis clinical trial; binary longitudinal
  outcome.
- **Typical use:** Binary repeated-measures GLMM with unbalanced visit times.

#### `gopherdat2`
- **Format:** data.frame, 30 rows × 6 cols (`Site`, `year`, `shells`, `Area`, `density`,
  `prev`)
- **Variables:** `shells` (count of fresh gopher tortoise shell middens); `Area` (site area);
  `density` (tortoise density); `prev` (prevalence of previous year); `Site`/`year`
- **Design:** Bolker et al. Poisson GLMM with offset for site area and site-level random
  effects; small dataset with overdispersion concern.
- **Typical use:** Demonstrates Poisson GLMM with `offset()` term.

#### `culcitalogreg` / `culcitalvolume`
- **Format:** data.frames, 80 and 50 rows respectively
- **Design:** Culcita (crown-of-thorns starfish) predation experiments; blocked design.
- **Typical use:** Logistic GLMM (`culcitalogreg`) and continuous response (`culcitalvolume`)
  with block random effects. Used in the lme4 "Fitting GLMMs with lme4" vignette.

---

### 1.2 lmerTest Datasets

#### `carrots`
- **Format:** data.frame; consumers (103) × products (12) panel study
- **Variables:** `Consumer` (factor); `Product` (factor, 12 levels); `Preference`,
  `Sweetness`, `Bitterness`, `Crispness` (7-point scales); `sens1`, `sens2` (PCA scores);
  demographic variables (`Frequency`, `Gender`, `Age`, `Homesize`, `Work`, `Income`)
- **Design:** Danish consumer preference mapping of 12 carrot varieties harvested 1996.
  Sensory panel PCA scores serve as predictors for consumer preference.
- **Typical use:** Demonstrates lmerTest random-slope model and Satterthwaite F-tests:
  `lmer(Preference ~ sens2 + Homesize + (1 + sens2|Consumer), carrots)`.

#### `ham`
- **Format:** data.frame; 81 consumers × 4 products × 2 information levels
- **Variables:** `Consumer` (factor, 81 levels); `Product` (factor, 4 levels); `Information`
  (factor, 2 levels); `Informed.liking` (numeric); `Gender` (factor); `Age` (numeric)
- **Design:** Næs et al. (2010) conjoint study of dry-cured ham. Consumers rated 4 products
  under true/false information conditions.
- **Typical use:** Demonstrates lmerTest `step()` backward elimination and `anova()` with
  Satterthwaite or Kenward-Roger denominator df:
  `lmer(Informed.liking ~ Product*Information + (1|Consumer), ham)`.

#### `TVbo`
- **Format:** data.frame; 8 assessors × 12 TV sets (3 TVset × 4 Picture combinations)
- **Variables:** `Assessor` (factor, 8 levels); `TVset` (factor, 3 levels); `Picture` (factor,
  4 levels); 15 numeric sensory response variables (Coloursaturation, Colourbalance, Noise,
  Depth, Sharpness, Lightlevel, Contrast, Sharpnessofmovement, Flickeringstationary,
  Flickeringmovement, Distortion, Dimglasseffect, Cutting, Flossyedges, Elasticeffect)
- **Design:** Bang & Olufsen TV sensory assessment. Trained panel of 8 evaluating 12 products
  on 15 attributes.
- **Typical use:** Demonstrates `ranova()` (random-effects ANOVA table) and multi-response
  lmerTest analyses with crossed assessor × product random effects.

---

## 2. Formula Processing Utilities

These functions manipulate mixed-model formulas. As of lme4 ≥ 1.1-30, the
implementations live in the **`reformulas`** package; lme4 exports shim wrappers
that warn and delegate. Users calling them via lme4 still see them in the lme4
namespace.

### `findbars(term)`
- Traverses a formula and returns a list of all random-effects terms (the
  expressions between `(` and `|`). Users rely on this to programmatically
  extract the random-effects specification from a formula.

### `nobars(term)`
- Returns the fixed-effects portion of a mixed model formula, stripping all
  random-effects terms. Useful for extracting the FE design matrix formula or
  comparing to a corresponding fixed-effects model.

### `subbars(term)`
- Substitutes `|` with `+` in all random-effects terms, returning a plain
  formula suitable for `model.frame` construction. Internal utility but exposed
  for users building custom parsers.

### `expandDoubleVerts(term)`
- Expands the `||` (double-bar) notation, which specifies independent random
  effects (zero correlations), into a sum of separate `(0+x|g)` terms. Allows
  users to inspect what `||` expands to, or to use it programmatically.

### `reOnly(f, response=FALSE, bracket=TRUE, doublevert_split=TRUE)`
- Returns a formula containing only the random-effects terms (i.e., a formula
  whose RHS consists solely of the bar-containing terms). Used when the user
  wants a standalone RE formula, e.g. for `simulate()` with `re.form=`.

### `isNested(f1, f2)`
- Tests whether every level of factor `f1` occurs with exactly one level of
  factor `f2`. Returns `TRUE` if `f1` is nested within `f2`. Used to verify
  design structure before choosing a nested vs crossed random-effects
  specification.

---

## 3. Random-Effects Term Construction

### `mkReTrms(bars, fr, drop.unused.levels=TRUE, reorder.terms=TRUE, reorder.vars=FALSE, calc.lambdat=TRUE, sparse=NULL)`
*(now in `reformulas`, re-exported by lme4)*

- Takes the list of bar-expressions from `findbars()` and a model frame `fr`,
  and returns a list describing the random-effects structure: the transposed
  relative covariance factor `Zt`, the covariance template `theta`, grouping
  factor levels, etc. This is the core internal constructor for all lme4 random
  effects; it is exposed so that advanced users and package authors can build
  deviance functions manually (part of the modular fitting API).

Key return components:
- `Zt` — transposed RE design matrix (sparse)
- `theta` — starting covariance parameters
- `Lind` — mapping from `theta` to lower Cholesky factor
- `flist` — list of grouping factors
- `Gp` — group pointers into Zt

### `mkNewReTrms(object, newdata, re.form=NULL, na.action=na.pass, allow.new.levels=FALSE, sparse=...)`
- Constructs random-effects terms for new data from an existing fitted model.
  The primary vehicle for prediction at new levels and for constructing the
  prediction random-effects matrix. Users rely on it indirectly through
  `predict.merMod()` and `simulate.merMod()`.

---

## 4. Modular Fitting Functions

These functions expose each step of the `lmer`/`glmer` fitting pipeline
individually, letting users intercept, modify, or replace any step.
The full pipeline is documented under `?modular`.

### 4.1 LMM Modular Pipeline

#### `lFormula(formula, data=NULL, REML=TRUE, subset, weights, na.action, offset, contrasts=NULL, control=lmerControl(), ...)`
- **Step 1: parse and validate.** Parses the formula and data, checks for
  errors, constructs the model frame (`fr`), fixed-effects design matrix (`X`),
  and random-effects terms list (`reTrms`). Returns a named list with components
  `fr`, `X`, `reTrms`, `REML`, `formula`, `wmsgs` (warnings).
- Users call this to inspect what lme4 derived from their formula/data before
  any numerical computation begins; also enables custom pre-processing.

#### `mkLmerDevfun(fr, X, reTrms, REML=TRUE, start=NULL, verbose=0, control=lmerControl(), ...)`
- **Step 2: build deviance function.** Takes the output of `lFormula` (minus
  `formula`) and returns a compiled objective function (an R closure over a C++
  environment) that computes the REML or ML deviance as a function of the
  covariance parameter vector `theta`.
- Users call this to obtain a raw deviance function for custom optimization,
  profiling, or visualization of the likelihood surface.

#### `optimizeLmer(devfun, optimizer=..., restart_edge=..., boundary.tol=..., start=NULL, verbose=0L, control=list(), ...)`
- **Step 3: optimize.** Minimizes the deviance function from `mkLmerDevfun`
  over `theta`. Returns a list with `par` (optimal theta), `fval` (deviance),
  `conv` (convergence code), and `message`. Uses bobyqa by default.
- Users call this to run custom optimizers or to re-run the optimization with
  different starting values or control settings without rebuilding the deviance
  function.

### 4.2 GLMM Modular Pipeline

#### `glFormula(formula, data=NULL, family=gaussian, subset, weights, na.action, offset, contrasts=NULL, start, mustart, etastart, control=glmerControl(), ...)`
- **Step 1 (GLMM): parse and validate.** Analogous to `lFormula` for GLMMs.
  Additionally handles `family`, `mustart`, and `etastart`. Returns the same
  component list.

#### `mkGlmerDevfun(fr, X, reTrms, family, nAGQ=..., verbose=0L, maxit=100L, control=glmerControl(), ...)`
- **Step 2 (GLMM): build deviance function.** Builds a penalized iteratively
  reweighted least-squares (PIRLS) objective. `nAGQ` controls the quadrature
  approximation: `nAGQ=0` uses the Laplace approximation (fastest); `nAGQ≥1`
  uses adaptive Gauss-Hermite quadrature (more accurate, scalar RE only).

#### `optimizeGlmer(devfun, optimizer=..., restart_edge=FALSE, boundary.tol=..., verbose=0L, control=list(), nAGQ=..., stage, start=NULL, ...)`
- **Step 3 (GLMM): optimize.** Two-stage optimization: Stage 1 (`nAGQ=0`)
  optimizes over `theta` only (fast Laplace); Stage 2 (`nAGQ≥1`) jointly
  optimizes over `theta` and `beta` with full quadrature. The `stage` argument
  selects which stage to run.

#### `updateGlmerDevfun(devfun, reTrms, nAGQ=1L)`
- **Between stages (GLMM).** Upgrades a Stage-1 GLMM deviance function (nAGQ=0)
  to a Stage-2 deviance function (nAGQ≥1) using the current `theta` estimates.
  Required when running the two-stage GLMM pipeline manually.

### 4.3 Shared Final Step

#### `mkMerMod(rho, opt, reTrms, fr, mc, lme4conv=NULL)`
- **Final step.** Takes the C++ deviance function environment (`rho`), the
  optimizer result (`opt`), the RE terms list, model frame, and the original
  `match.call()` object, and assembles a fully populated `merMod` (or
  `lmerMod`/`glmerMod`) S4 object. This is the object returned by `lmer`/`glmer`.
- Users call this when they have run the modular pipeline manually and need to
  convert results into a standard `merMod` object for use with all downstream
  methods (summary, anova, predict, etc.).

---

## 5. Quadrature Rules

### `GHrule(ord, asMatrix=TRUE)`
- Computes a univariate Gauss-Hermite quadrature rule of order `ord` (1–100).
  Returns `ord` rows with columns `z` (node positions), `w` (weights), `ldnorm`
  (log normal density at nodes). The rule integrates exactly polynomials of
  degree up to 2k-1 times the standard normal density.
- Users rely on this to understand/verify the quadrature approximation used in
  `glmer(nAGQ=k)`, or to implement custom Gauss-Hermite-based quadrature in
  other functions.

### `GQdk(d=1L, k=1L)`
- Alternative interface to multivariate Gauss-Hermite quadrature: returns nodes
  and weights for `d`-dimensional integration with `k` points per dimension
  (a sparse grid). Less commonly used directly; underlies multivariate RE
  integration.

### `GQN`
- Internal object (not a function with standard args) containing pre-computed
  Gauss-Hermite quadrature nodes/weights tables. Exposed but primarily for
  internal use.

---

## 6. Dummy Variable Utility

### `dummy(f, levelsToKeep)`
- **Arguments:** `f` — a factor (or coercible object); `levelsToKeep` — optional
  character vector of factor levels to include.
- Returns a `model.matrix` with `nlevels(f)-1` dummy columns (first level dropped
  by default), or only the specified levels if `levelsToKeep` is supplied.
- Users rely on this when constructing random-effects design matrices for
  specific factor contrasts without triggering `model.matrix`'s intercept-column
  behavior. Example use: `(0 + dummy(Sex, "Female") | Subject)` for a
  sex-specific random slope with no correlation to the intercept.

---

## 7. Simulation / Template Utilities

### `mkDataTemplate(formula, data, nGrps=2, nPerGrp=1, rfunc=NULL, ...)`
- Generates a template data frame matching the structure implied by `formula`
  and `data`, with `nGrps` levels per grouping factor and `nPerGrp` observations
  per group. `rfunc` can supply a covariate-generating function.
- Users call this to create synthetic balanced data frames as a skeleton for
  `simulate.merMod()`, especially when exploring power or sample-size questions.
  Described as EXPERIMENTAL in docs.

### `mkParsTemplate(formula, data)`
- Returns a named list with template parameter vectors (`beta`, `theta`, `sigma`)
  compatible with the model implied by `formula` and `data`. Intended to be filled
  in by the user and passed to `simulate.merMod(newparams=...)`.
- Used together with `mkDataTemplate` to simulate data from a model with
  user-specified parameters rather than estimated ones.

---

## 8. Model Diagnostics Utilities

### `isSingular(object, tol=getSingTol())`
- Tests whether a fitted `merMod` is near-singular: any random-effects
  covariance matrix has variance components (in the PCA sense) close to zero.
  `tol` defaults to `getSingTol()` (currently 1e-4).
- Users call this as a post-fit check; models returning `TRUE` may be
  overparameterized and should be simplified.

### `getSingTol()`
- Returns the tolerance used by `isSingular()`. Exposed so users can inspect
  or change the threshold globally via `options()`.

### `rePCA(x)`
- Performs PCA on the random-effects variance-covariance matrices of a fitted
  `merMod`. Returns a `prcomplist` giving standard deviations of orthogonal
  variance components and rotation matrices.
- More informative than `isSingular` alone: shows _which_ components are near
  zero and the mapping from model parameters to principal directions. Used to
  diagnose overfitting in complex random-effects structures (Bates et al. 2015).

### `ngrps(object, ...)`
- Returns the number of levels of each grouping factor for a `merMod` object
  (or a scalar for a plain `factor`). For models with multiple grouping factors
  returns a named integer vector.
- Used to quickly inspect the random-effects sample sizes (number of clusters
  per random effect).

---

## 9. Convergence / Optimizer Utilities

### `allFit(object, meth.tab=NULL, data=NULL, verbose=TRUE, show.meth.tab=FALSE, maxfun=1e5, parallel=c("no","multicore","snow"), ncpus=1L, cl=NULL, catch.errs=TRUE, start_from_mle=TRUE)`
- Re-fits a model with all available optimizers (bobyqa, Nelder_Mead, nlminbwrap,
  plus wrappers via `optimx` and `dfoptim::nmkb`) and compares results. Returns
  a list of fitted models with summary methods.
- Users call this when convergence warnings appear, to verify whether the
  parameter estimates are robust across optimizers or whether a genuine
  convergence failure occurred. `parallel` enables multi-core execution.
  Requires `optimx` and `dfoptim` for full optimizer coverage.

### `checkConv(attr, val, conv, control, lbound)`
- Internal convergence checker exposed in the namespace. Evaluates optimizer
  convergence codes, gradient magnitudes, and Hessian eigenvalues against
  control thresholds. Rarely called directly by users but accessible for
  diagnostics.

---

## 10. Internal / Advanced Utilities (Exported but Low-Level)

These are exported and thus part of the public API, but primarily intended for
package authors extending lme4 or for advanced internal use.

| Function | Brief description |
|---|---|
| `getME(object, name)` | Generic accessor for named internal components of a `merMod` (e.g. `"theta"`, `"beta"`, `"Zt"`, `"Lambda"`, `"flist"` — about 30 named slots). Superseded in part by named extractors but remains the most general accessor. |
| `lmerControl(...)` / `glmerControl(...)` | Build control lists for `lmer`/`glmer`; covered in the LMM-fit-control survey. Listed here because `lFormula`/`glFormula` accept them and they govern the modular pipeline. |
| `isLMM(object)` / `isGLMM(object)` / `isNLMM(object)` / `isREML(object)` | Predicate functions testing model class membership. Used to branch on model type in generic functions. |
| `isNewMerMod(object)` | Tests whether an object is a `merMod` (new-style lme4 ≥ 1.0). |
| `anyStructured(object)` | Tests whether a model has any structured (non-diagonal) covariance components. |
| `getL(object)` | Returns the sparse lower Cholesky factor of the random-effects precision matrix. |
| `getReCovs(object, ...)` | Returns the list of RE covariance matrices. |
| `mkVarCorr(sc, cnms, nc, theta, nms)` | Constructs a `VarCorr`-style list from raw components. Package-author utility. |
| `formatVC(varcor, ...)` | Formats a `VarCorr` object as a character matrix for printing. |
| `cov2sdcor(V)` / `sdcor2cov(M)` | Convert between covariance matrix and SD/correlation parameterization. |
| `Cv_to_Sv(Cv)` / `Sv_to_Cv(Sv)` / `Cv_to_Vv(Cv)` / `Vv_to_Cv(Vv)` | Conversions between correlation, SD, and variance parameterizations (capitalized variants). |
| `mlist2vec(L)` / `vec2mlist(v, ...)` / `vec2STlist(v, ...)` | Convert between list-of-matrices and vector parameterization of covariance components. |
| `devcomp(x, ...)` | Extracts the deviance components (objective, parameters, convergence info) from a `merMod`. |
| `llikAIC(object)` | Returns log-likelihood and AIC for a `merMod` in a single call. |
| `varianceProf(x, ...)` | Converts a profile object to variance-scale parameterization. |
| `logProf(x, ...)` | Converts a profile to log-scale. |
| `golden(f, ...)` | Golden-section search utility used internally for profile CI computation. |
| `Nelder_Mead(fn, par, ...)` / `NelderMead` | Nelder-Mead optimizer wrappers available for user use. `NelderMead` is the reference class; `Nelder_Mead` is the function wrapper. |
| `nlminbwrap(fn, par, ...)` | Wrapper for `nlminb` conforming to lme4's optimizer interface (accepts `lower`/`upper`, returns `par`/`fval`/`conv`). |
| `nloptwrap(fn, par, ...)` | Wrapper for `nloptr` optimizers conforming to lme4's optimizer interface. |
| `rePos` | Reference class for random-effects positions. |
| `merPredD` | Reference class wrapping the C++ `merPredD` object (dense random effects predictor). |
| `lmerResp` / `glmResp` / `lmResp` / `nlsResp` | Reference classes for model response objects. Low-level building blocks of the deviance function. |
| `glmFamily(family)` | Reference class wrapping a GLM family object for use in C++. |
| `namedList(...)` | Utility to create a named list where names are taken from variable names. Minor convenience function. |
| `factorize(val, frloc, ...)` | Coerces variables in a model frame to factors. |
| `lme4_testlevel()` | Returns the current test level (integer 0–3) for conditional test execution in lme4's own test suite. Exposed for external package test suites that want to match lme4's verbosity. |
| `forceNewMerMod(object)` | Coerces an old-style `mer` object to a new-style `merMod`. Legacy compatibility. |

---

## 11. lmerTest Utility Functions

These are distinct from the inference-focused lmerTest functions (anova, ranova,
etc.) but still qualify as utilities in the sense of model manipulation and
dataset access.

### `as_lmerModLmerTest(model, tol=1e-8, ...)`
- Converts an existing `lmerMod` (from plain lme4 `lmer`) to a
  `lmerModLmerTest` object, enabling lmerTest inference methods on a model
  that was originally fitted with lme4.
- Users rely on this when they fit a model with lme4 first (e.g. for speed or
  compatibility reasons) and then want lmerTest-style Satterthwaite/KR
  p-values without refitting.

### `get_model(step_object)`
- Extracts the final model from a `step` object returned by lmerTest's
  `step()` backward-elimination procedure.
- Users call this after `step(fm)` to retrieve the selected final model as a
  standard `lmerModLmerTest` object for further inference.

### `show_tests(object, fractions=FALSE, ...)`
- Displays the contrast matrices used for the F-tests in lmerTest's `anova()`
  output, showing precisely which linear combinations of fixed effects are
  being tested.
- Users call this to understand and verify the estimability and test structure
  for type I/II/III ANOVA tables; essential for understanding non-orthogonal
  designs.

---

## 12. Summary of Coverage by Sub-Category

| Sub-category | Key items | Count |
|---|---|---|
| lme4 datasets | sleepstudy, cbpp, cbpp2, Pastes, Penicillin, Dyestuff, Dyestuff2, cake, grouseticks, grouseticks_agg, VerbAgg, InstEval, Arabidopsis, salamander, schizophrenia, toenail, gopherdat2, culcitalogreg, culcitalvolume | 19 |
| lmerTest datasets | carrots, ham, TVbo | 3 |
| Formula manipulation | findbars, nobars, subbars, expandDoubleVerts, reOnly, isNested | 6 |
| RE term construction | mkReTrms, mkNewReTrms | 2 |
| Modular fitting (LMM) | lFormula, mkLmerDevfun, optimizeLmer | 3 |
| Modular fitting (GLMM) | glFormula, mkGlmerDevfun, optimizeGlmer, updateGlmerDevfun | 4 |
| Final assembly | mkMerMod | 1 |
| Quadrature | GHrule, GQdk, GQN | 3 |
| Dummy variables | dummy | 1 |
| Simulation templates | mkDataTemplate, mkParsTemplate | 2 |
| Model diagnostics | isSingular, getSingTol, rePCA, ngrps | 4 |
| Convergence utilities | allFit, checkConv | 2 |
| Internal/advanced exports | getME, isLMM, isGLMM, isNLMM, isREML, isNewMerMod, anyStructured, getL, getReCovs, mkVarCorr, formatVC, cov2sdcor, sdcor2cov, Cv_to_Sv et al., mlist2vec, vec2mlist, vec2STlist, devcomp, llikAIC, varianceProf, logProf, golden, Nelder_Mead, NelderMead, nlminbwrap, nloptwrap, rePos, merPredD, lmerResp/glmResp/lmResp/nlsResp, glmFamily, namedList, factorize, lme4_testlevel, forceNewMerMod | ~35 |
| lmerTest utilities | as_lmerModLmerTest, get_model, show_tests | 3 |

**Total items documented: ~93**
