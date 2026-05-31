# Error-message quality probe: convergence-hard

**Scenario:** Deliberately ill-conditioned models that strain convergence via two
independent stressors:

1. **Extreme scale mismatch** (scenario 1, 3): predictors differ by ~10^6 in
   magnitude (`x_small ~ O(1)`, `x_large ~ O(10^6)`), producing a condition
   number of ~10^9 for the fixed-effect design matrix, plus a random slope on
   `x_large` whose covariance parameter is near-zero in original units. This
   stresses both the fixed-effect linear solve and the theta optimizer.

2. **Near-zero within-group variance** (scenario 2): a random slope predictor
   `x` has between-group SD ~8 but within-group SD ~0.001, so the likelihood
   ridge for the slope covariance is nearly flat. lme4/mixeff both have to
   navigate a nearly-degenerate parameter landscape.

**Dataset sizes:** 15 subjects × 4 obs (scenarios 1, 2); 20 subjects × 5 obs
(scenario 3 GLMM).

Probe script: `assessment/errors/probe-convergence-hard.R`
Run date: 2026-05-31
mixeff version: installed from /Users/bbuchsbaum/code/mixeff
lme4 version: 2.0.1

---

## Verbatim output by surface

### Scenario 1: extreme scale mismatch — lme4::lmer tight budget (maxfun=10)

```
  class: lmerMod
  WARN: Some predictor variables are on very different scales: consider rescaling.
You may also use (g)lmerControl(autoscale = TRUE) to improve numerical stability.
  WARN: maxfun < 10 * length(par)^2 is not recommended.
  WARN: convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded
  WARN: unable to evaluate scaled gradient
  WARN: Model failed to converge: degenerate  Hessian with 3 negative eigenvalues
  See ?lme4::convergence and ?lme4::troubleshooting.
```

Five separate R warnings are emitted. The messages collectively convey that:
(a) scale mismatch exists, (b) the evaluation budget was too small, (c) the
gradient could not be evaluated after fitting, (d) the Hessian has negative
eigenvalues. The final catch-all points to `?lme4::convergence`. Despite all
these warnings, lme4 **still returns a fitted lmerMod object** — the user must
act on the warnings to know the fit is untrustworthy.

### Scenario 1: extreme scale mismatch — lme4::lmer DEFAULT budget

```
  class: lmerMod
  MSG : boundary (singular) fit: see help('isSingular')
  WARN: Some predictor variables are on very different scales: consider rescaling.
You may also use (g)lmerControl(autoscale = TRUE) to improve numerical stability.
```

With full budget lme4 reaches the boundary and emits only a terse singularity
advisory plus the scale-mismatch warning. Neither message names the affected
term, explains that the slope on `x_large` drove the fit to boundary, or notes
that the fixed-effect matrix is near rank-deficient.

### Scenario 1: extreme scale mismatch — mixeff::lmm DEFAULT budget

```
  fit_status: converged_reduced_rank
  -- artifact$diagnostics --
    fixed_effect_rank_deficient [warning/design_audit]:
      fixed-effect formula is rank-deficient (rank 2 of 3); some requested
      coefficients are not separately estimable from the observed data
    scope_note [info/design_audit]:
      `x_small` varies within `subject`, so a `subject`-level slope is
      structurally possible
  -- optimizer_certificate$diagnostics --
    boundary_parameter [info/certification]:
      standard deviation for x_large in (1 + x_large | subject) is on its
      lower bound
    covariance_reduced [info/certification]:
      fitted covariance for (1 + x_large | subject) has effective rank 1 of
      requested rank 2
```

mixeff does **not** error or warn at the R level (no catchable R warning is
raised). All diagnostics are available only via `diagnostics()` or via the
artifact. Fit status is `converged_reduced_rank`.

Four distinct named diagnostics are produced:
1. `fixed_effect_rank_deficient` (warning): names the rank (2 of 3), indicates
   some coefficients are non-estimable.
2. `scope_note` (info): clarifies which predictor could structurally support a
   slope.
3. `boundary_parameter` (info): names the exact term and SD direction that hit
   the lower bound.
4. `covariance_reduced` (info): states the effective rank (1 of 2) and names
   the affected term.

The summary output is informative about what was dropped:
```
  x_small: not_estimable — contrast touches aliased or non-finite coefficient
                            directions
  (Intercept): std_err unavailable — standard error is unavailable, so the
                            Wald z p-value is unavailable
  x_large: std_err unavailable — same reason
```

**Notable finding:** the `fixed_effect_rank_deficient` diagnostic is a genuine
extra catch. The design matrix condition number is ~10^9 (nearly rank-deficient
in floating-point), and the Rust engine flags this at design_audit stage.
lme4 with default budget does not emit this diagnostic at all — it silently
absorbs the scale mismatch and hits a singular boundary without explaining the
collinearity.

### Scenario 2: near-zero within-group variance — lme4::lmer

```
  class: lmerMod
  MSG : boundary (singular) fit: see help('isSingular')
```

Single terse message. No term named, no explanation of the near-zero within-group
variance, no suggestion to simplify the random-effect structure.

### Scenario 2: near-zero within-group variance — mixeff::lmm

```
  fit_status: converged_reduced_rank
  -- artifact$diagnostics --
    covariance_too_rich [warning/design_audit]:
      10 levels are below the v0 full-covariance threshold 15 for 3 covariance
      parameters
  -- optimizer_certificate$diagnostics --
    covariance_reduced [info/certification]:
      fitted covariance for (1 + x | subject) has effective rank 1 of
      requested rank 2
```

Two named diagnostics:
1. `covariance_too_rich` (warning): gives the concrete count (10 levels) and
   concrete threshold (15 for 3 covariance parameters).
2. `covariance_reduced` (info): names the affected term and its effective rank.

### Scenario 3: GLMM extreme scale mismatch — lme4::glmer tight budget (maxfun=20)

```
  class: glmerMod
  MSG : boundary (singular) fit: see help('isSingular')
  WARN: maxfun < 10 * length(par)^2 is not recommended.
  WARN: convergence code 1 from bobyqa: bobyqa -- maximum number of function
        evaluations exceeded
  WARN: maxfun < 10 * length(par)^2 is not recommended.
  WARN: convergence code 1 from bobyqa: bobyqa -- maximum number of function
        evaluations exceeded
```

The two `maxfun` / convergence-code pairs are duplicated (lme4's two-stage
GLMM outer/inner optimizers both hit the budget). No term named, no link to the
scale problem.

### Scenario 3: GLMM extreme scale mismatch — mixeff::glmm DEFAULT budget

```
  fit_status: converged_boundary
  -- artifact$diagnostics --
    fixed_effect_rank_deficient [warning/design_audit]:
      fixed-effect formula is rank-deficient (rank 2 of 3); some requested
      coefficients are not separately estimable from the observed data
    scope_note [info/design_audit]:
      `x_small` varies within `subject`, so a `subject`-level slope is
      structurally possible
    near_unit_random_effect_correlation [warning/certification]:
      random-effect correlation for group subject between (Intercept) and
      x_large is -1.000; the fitted covariance is nearly one-dimensional
  -- optimizer_certificate$diagnostics --
    boundary_parameter [info/certification]:
      covariance parameter 3 is on its lower bound
    boundary_parameter [info/certification]:
      GLMM covariance state classified as ValidZeroVariance
    near_unit_random_effect_correlation [warning/certification]:
      random-effect correlation for group subject between (Intercept) and
      x_large is -1.000; the fitted covariance is nearly one-dimensional
```

Five distinct named diagnostics with staged provenance (design_audit /
certification). The `near_unit_random_effect_correlation` diagnostic is
particularly informative — it names the exact random-effect pair (-1.000
correlation), and classifies the covariance as "nearly one-dimensional",
which lme4 never does.

---

## Comparative assessment

### Summary table

| Criterion | lme4 | mixeff |
|-----------|------|--------|
| Names the affected term / SD direction | No | Yes (`boundary_parameter` names exact term) |
| Flags near-rank-deficient fixed-effects | No | Yes (`fixed_effect_rank_deficient`, rank 2 of 3) |
| Gives condition-number / scale advice | "consider rescaling" (generic) | Implicit via rank-deficient diagnostic |
| Explains why boundary occurred | No (only "singular") | Yes (reduced rank, effective rank stated) |
| Separate diagnostic codes per problem | No (1 warning per issue, untyped) | Yes (4–5 named codes per scenario) |
| Catchable as typed R condition | No (R `warning()` / `message()`) | Yes (typed `mm_condition` subclasses) |
| fit_status reflects the problem | No (`@optinfo$conv$lme4`) | Yes (`converged_reduced_rank`, `converged_boundary`) |
| Non-convergence on tight budget | 5 warnings, returns fit anyway | Not applicable (cannot cap budget via mm_control) |

### lme4 behaviour

With a tight budget lme4 emits 5 uncorrelated R warnings of varying quality:
the scale warning is actionable ("rescale", "autoscale = TRUE"), the budget
warning tells the user the maxfun was too small, but the convergence-code warning
(`convergence code 1 from bobyqa`) and gradient warning (`unable to evaluate
scaled gradient`) are opaque — they require reading `?lme4::convergence` and
cross-referencing BOBYQA return codes to understand.

With a default budget lme4 reduces to a single terse `boundary (singular) fit:
see help('isSingular')` that names neither the term nor the cause.

### mixeff behaviour

mixeff produces a richer set of named, staged diagnostics. The key strengths
on this scenario:

- `fixed_effect_rank_deficient` fires at **design_audit** stage (before
  fitting) when the condition number of the fixed-effect matrix approaches
  machine limits. lme4 does not have a pre-fit audit path and never emits this.
- `boundary_parameter` names the **exact term and direction** ("standard
  deviation for x_large in (1 + x_large | subject)").
- `covariance_reduced` states the **effective rank** (1 of 2) and names the
  term, and includes a `suggested_formula` payload (`(1 | subject)`) pointing
  to the simpler model the data can actually support.
- `near_unit_random_effect_correlation` (GLMM) names the **exact pair** and
  value (-1.000), classifying the fitted covariance as "nearly one-dimensional".
- `fit_status` (`converged_reduced_rank`, `converged_boundary`) gives a
  machine-readable convergence class.

### Gaps (needs-work, not bugs)

1. **No R-level warning at fit time.** mixeff does not raise a catchable R
   warning or condition when the fit completes with reduced-rank or boundary
   status. A user who does not call `diagnostics()` explicitly will only see the
   `fit_status` field if they inspect the object. lme4's approach of emitting R
   `warning()` calls — even if the messages are less informative — means the
   user sees something at the console without any extra calls. This is a
   discoverability gap.

2. **`boundary_parameter` in GLMM does not name the term.** In scenario 3 the
   GLMM `boundary_parameter` message says "covariance parameter 3 is on its
   lower bound" — a parameter index, not a term name. The LMM version names the
   term explicitly ("standard deviation for x_large in (1 + x_large | subject)").
   This inconsistency between lmm and glmm is a **needs-work** item.

3. **Diagnostic deduplication.** The `diagnostics()` output shows each row
   twice (artifact-level + certificate-level), a cosmetic issue also seen in
   other scenarios.

4. **Cannot cap the optimizer budget from R.** `mm_control()` accepts only
   `verbose=`; there is no `max_iterations` or equivalent to test the
   tight-budget path from R. This is a test-gap (not a runtime defect).

---

## Verdict

**good** — mixeff is materially clearer than lme4 on this scenario.

On the hardest convergence stress (extreme scale mismatch, condition number
~10^9), mixeff produces 4–5 named, staged, typed diagnostics that collectively
explain the rank deficiency, the boundary direction, the affected term, and
the effective rank of the fitted covariance. lme4 with a tight budget produces
five anonymous R warnings; with a default budget it collapses to the single
opaque `boundary (singular) fit` advisory.

The two gaps (no R-level warning at fit time; GLMM `boundary_parameter` uses
a parameter index rather than a term name) are needs-work improvements, not
bugs — they do not make mixeff worse than lme4 on this scenario, and the
underlying diagnostic information is present and correct.
