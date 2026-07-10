# Package index

## Fitting models

Construct and fit linear and generalized linear mixed models.

- [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) : Fit
  a linear mixed-effects model

- [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) :
  Fit a generalized linear mixed model

- [`mm_negative_binomial()`](https://bbuchsbaum.github.io/mixeff/reference/mm_negative_binomial.md)
  :

  Negative-binomial family for
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)

- [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  : Compile a mixed-effects model spec without fitting

- [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
  : Control mixeff fitting behavior

## Extractors (lme4-compatible)

The familiar lme4 accessor surface.

- [`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`coef(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`coef(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`sigma(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`sigma(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`logLik(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`logLik(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`deviance(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`deviance(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`AIC(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`AIC(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`BIC(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`BIC(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`nobs(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`nobs(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`df.residual(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`df.residual(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`formula(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`formula(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`model.frame(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`model.frame(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`ngrps()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`weights(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`weights(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`extractAIC(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`extractAIC(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`terms(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`terms(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`as.data.frame(`*`<mm_varcorr>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`as.data.frame(`*`<mm_ranef>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`model.matrix(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`model.matrix(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`vcov(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  [`vcov(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  : Extract components from a fitted mixeff LMM
- [`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md) :
  Extract low-level model components
- [`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md)
  : Test whether a fit is singular or reduced-rank
- [`predict(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/predict.mm_lmm.md)
  [`fitted(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/predict.mm_lmm.md)
  [`residuals(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/predict.mm_lmm.md)
  [`fitted(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/predict.mm_lmm.md)
  [`residuals(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/predict.mm_lmm.md)
  : Predict from a fitted mixeff LMM
- [`predict(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/predict.mm_glmm.md)
  : Predict from a fitted mixeff GLMM
- [`update(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/update.mm.md)
  [`update(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/update.mm.md)
  : Update and re-fit a mixeff model
- [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md) :
  Refit a mixeff LMM with a new response
- [`simulate(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/simulate.mm_lmm.md)
  : Simulate from a mixeff LMM

## Inference

Contrasts, term tests, degrees of freedom, intervals, and bootstrap.

- [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
  : Contrast fixed effects
- [`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md)
  : Test a fixed-effect term
- [`test_random_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_random_effect.md)
  : Test a random-effect variance component
- [`anova(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/anova.mm_glmm.md)
  : Analysis of deviance for GLMMs
- [`drop1(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/drop1.mm_lmm.md)
  : Drop one fixed-effect term at a time
- [`drop1(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/drop1.mm_glmm.md)
  : Drop one fixed-effect term at a time from a GLMM
- [`confint(`*`<mm_glmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/confint.mm_glmm.md)
  : Confidence intervals for fixed effects of a mixeff GLMM
- [`profile(`*`<mm_lmm>`*`)`](https://bbuchsbaum.github.io/mixeff/reference/profile.mm_lmm.md)
  : Profile a fitted linear mixed model
- [`df_for_contrast()`](https://bbuchsbaum.github.io/mixeff/reference/df_for_contrast.md)
  : Degrees of freedom for a contrast
- [`estimability()`](https://bbuchsbaum.github.io/mixeff/reference/estimability.md)
  : Assess contrast estimability
- [`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
  : Inspect inference methods available for this fit
- [`inference_table()`](https://bbuchsbaum.github.io/mixeff/reference/inference_table.md)
  : Fixed-effect inference table
- [`bootstrap_control()`](https://bbuchsbaum.github.io/mixeff/reference/bootstrap_control.md)
  : Fixed-effect bootstrap control
- [`parametric_bootstrap()`](https://bbuchsbaum.github.io/mixeff/reference/parametric_bootstrap.md)
  : Parametric bootstrap likelihood-ratio comparison

## Marginal effects & emmeans

- [`mm_grid()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
  [`mm_predictions()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
  [`mm_means()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
  [`mm_comparisons()`](https://bbuchsbaum.github.io/mixeff/reference/mm_grid.md)
  : Marginal grids, predictions, means, and comparisons
- [`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md)
  : Wald inference on a linear combination of fixed effects
- [`recover_data.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/emmeans-support.md)
  [`emm_basis.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/emmeans-support.md)
  [`recover_data.mm_glmm()`](https://bbuchsbaum.github.io/mixeff/reference/emmeans-support.md)
  [`emm_basis.mm_glmm()`](https://bbuchsbaum.github.io/mixeff/reference/emmeans-support.md)
  : Optional emmeans support for mixeff LMMs

## Model comparison & reporting

- [`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)
  : Compare fitted mixeff models
- [`compare_covariance()`](https://bbuchsbaum.github.io/mixeff/reference/compare_covariance.md)
  : Compare covariance parameterizations for current random terms
- [`model_report()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md)
  [`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md)
  : Produce reporting tables for a fitted mixeff model

## Audit-first verbs

Make the model’s design, reductions, and optimizer state explicit.

- [`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md) :
  Audit a compiled model spec or fitted model

- [`audit_design()`](https://bbuchsbaum.github.io/mixeff/reference/audit_design.md)
  :

  Deprecated alias for
  [`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md)

- [`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
  : Explain the random-effects structure of a compiled model

- [`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
  : Show requested, effective, and fitted model-state changes

- [`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
  [`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
  : Inspect mixeff diagnostics and fit status

- [`parameterization()`](https://bbuchsbaum.github.io/mixeff/reference/parameterization.md)
  : Inspect covariance parameterization

- [`random_options()`](https://bbuchsbaum.github.io/mixeff/reference/random_options.md)
  : Inspect nearby random-effect spellings for one grouping factor

- [`random_blocks()`](https://bbuchsbaum.github.io/mixeff/reference/random_blocks.md)
  : Inspect random-effect blocks

- [`roles()`](https://bbuchsbaum.github.io/mixeff/reference/roles.md) :
  Declare or inspect design roles

- [`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md)
  : Inspect the optimizer certificate

- [`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md)
  : Verify convergence of a fitted linear mixed model

- [`reproducibility()`](https://bbuchsbaum.github.io/mixeff/reference/reproducibility.md)
  : Inspect reproducibility metadata

## broom / broom.mixed

- [`mm_broom`](https://bbuchsbaum.github.io/mixeff/reference/mm_broom.md)
  : Tidy, glance, and augment methods for mixeff fits

## Persistence

- [`as_json()`](https://bbuchsbaum.github.io/mixeff/reference/as_json.md)
  : Serialize a mixeff spec or fit to JSON
- [`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
  : Revive a serialized mixeff object
- [`fit_handle_alive()`](https://bbuchsbaum.github.io/mixeff/reference/fit_handle_alive.md)
  : Test whether a mixeff fit has a live native handle

## Formula & schema internals

- [`mm_parse_formula()`](https://bbuchsbaum.github.io/mixeff/reference/mm_parse_formula.md)
  : Parse and canonicalize an lme4-style formula

- [`mm_formula_manifest()`](https://bbuchsbaum.github.io/mixeff/reference/mm_formula_manifest.md)
  : The wrapper's formula manifest

- [`mm_json_negotiate()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_negotiate.md)
  :

  Negotiate a JSON schema header against what `mixeff` supports

- [`mm_json_known_schemas()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_known_schemas.md)
  : Closed list of schema/version pairs the wrapper understands

## Package

- [`mixeff`](https://bbuchsbaum.github.io/mixeff/reference/mixeff-package.md)
  [`mixeff-package`](https://bbuchsbaum.github.io/mixeff/reference/mixeff-package.md)
  : mixeff: Audit-First Mixed-Effects Models via the 'mixedmodels' Rust
  Crate
