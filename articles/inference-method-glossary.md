# Inference Method Glossary

This page is the reference for the raw enum values that appear in
[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
and in the inference tables returned by
[`summary()`](https://rdrr.io/r/base/summary.html),
[`inference_table()`](https://bbuchsbaum.github.io/mixeff/reference/inference_table.md),
[`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md),
[`anova()`](https://rdrr.io/r/stats/anova.html),
[`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md),
and [`confint()`](https://rdrr.io/r/stats/confint.html).

Printed output may use reader-facing display columns such as
`display_status`, `display_reason`, and `what_to_do_next`. The values
below are the stable machine-readable contract for tests, reports, and
downstream tooling.

## `expected_status`

`expected_status` is used by
[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
to say whether a route is usable on the current fit.

| Value | Definition |
|----|----|
| `available` | The method can run on this fit and produce the requested quantity. |
| `not_assessed` | The method is known but refused for this fit, with the reason named in `expected_reliability_reason` or `reason_code`. |
| `not_yet_wired` | The method is part of the inference contract but is not exposed through the current R bridge. |

## `expected_reliability_reason` and `reliability_reason`

`expected_reliability_reason` is the reason
[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
attaches to a route on its map. `reliability_reason` is the reason the
engine attaches to an actual inference row when a verb runs. The
vocabularies overlap but are not identical, and where the engine has
recorded an observed reason for a route,
[`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
reports the observed value in preference to its own prediction — so on a
fitted LMM the Wald route usually shows the engine’s
`asymptotic_wald_z_fallback` rather than the two route-map entries
marked *(route map)* below.

| Value | Definition | Example |
|----|----|----|
| `not_available` | No additional reliability reason is attached to the row; inspect `status` and `reason_code` for the decision. | A refused cluster-bootstrap p-value row uses a specific `reason_code` while reliability itself is not available. |
| `asymptotic_wald_z_fallback` | A t reference distribution was the target, degrees of freedom could not be computed, and a standard normal was substituted; the row is graded `low`. | [`inference_table()`](https://bbuchsbaum.github.io/mixeff/reference/inference_table.md) and [`summary()`](https://rdrr.io/r/base/summary.html) on a boundary fit report their Wald rows with this reason. This is the reason code most users meet first. |
| `asymptotic_wald_z_at_boundary` *(route map)* | Wald z can run, but boundary or reduced-rank fits can make asymptotic uncertainty too optimistic. | Defined for the [`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md) route map; in practice the engine’s observed `asymptotic_wald_z_fallback` takes precedence. |
| `interior_converged_well_specified` *(route map)* | The fit is interior and the fast Wald route has no boundary caveat from the route map. | Defined for the route map; superseded by the observed reason in current output. |
| `satterthwaite_unavailable_at_boundary` | Satterthwaite degrees of freedom are not certified when variance parameters are on the boundary. | `test_effect(fit, term, method = "satterthwaite")` refuses on a singular fit. |
| `kenward_roger_unavailable_at_boundary` | Kenward-Roger degrees of freedom and covariance adjustment are not certified when variance parameters are on the boundary. | `test_effect(fit, term, method = "kenward_roger")` refuses on a singular fit. |
| `satterthwaite_finite_difference_approximation` | The route uses a finite-difference Satterthwaite degrees-of-freedom approximation. | The standard reason on Satterthwaite rows for interior fits. |
| `kenward_roger_approximation` | The route uses Kenward-Roger finite-sample degrees of freedom with covariance adjustment. | The standard reason on Kenward-Roger rows for interior fits. |
| `parametric_bootstrap_monte_carlo` | The p-value comes from parametric-bootstrap replicates; reliability is governed by replicate count and Monte-Carlo standard error. | `contrast(fit, L, method = "bootstrap")` rows carry this reason. |
| `glmm_joint_laplace_active_hessian_wald` | GLMM Wald inference from the joint-Laplace fit’s certified Hessian. | The `summary()$inference` rows of a `method = "joint_laplace"` GLMM fit (the printed coefficient block shows the statistics, not the reason). |
| `bootstrap_monte_carlo_replicates` *(route map)* | Reliability is governed by replicate count and Monte-Carlo standard error. | The [`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md) prediction for bootstrap routes; the observed row value is `parametric_bootstrap_monte_carlo`. |
| `bootstrap_insufficient_replicates` | The bootstrap route ran, but the requested replicate count is too small for a moderate reliability grade. | A demonstration call with `nsim = 10` uses this reason. |
| `bootstrap_not_certified` | A bootstrap payload was returned without a certified p-value. | A bootstrap LRT row with no engine-certified p-value uses this reason. |
| `inference_unavailable_by_policy` | The method is refused by an explicit policy for this fit and requested target. | A bootstrap LRT on a REML fit returns this reliability reason with a more specific `reason_code`. |
| `bootstrap_lrt_requires_ml` | A bootstrap likelihood-ratio test requires ML-fitted models, not REML-fitted models. | `test_effect(fit, term, method = "bootstrap_lrt")` refuses on a REML fit. |
| `bootstrap_cluster_resample_p_value_unavailable` | Cluster resampling describes estimator distributions and does not certify fixed-effect p-values in the current schema. | `test_effect(fit, term, method = "cluster_bootstrap", group = "subj")` refuses with this reason. |
| `profile_likelihood_ci` | Profile-likelihood confidence intervals are available for this requested parameter set. | [`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md) marks `profile_ci` available on an ML fit away from the boundary. |
| `profile_beta_unavailable_under_reml` | Fixed-effect beta profile intervals require an ML fit and are not certified from REML profile payloads. | `confint(fit, method = "profile")` on a REML fit records beta rows with this reason. |
| `profile_ci_unavailable_at_boundary` | Profile intervals are not certified for boundary or reduced-rank fits. | [`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md) refuses `profile_ci` on a singular fit. |

The `method` column has one value worth listing here as well:
`not_computed`, which appears on GLMM rows fitted by the default
`pirls_profiled` estimator, whose standard errors and p-values are
withheld (see
[`vignette("glmm")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md)).

## `reason_code`

`reason_code` is attached to actual unavailable or not-assessed rows and
to typed errors from inference verbs. It may match
`expected_reliability_reason`, but the two fields answer different
questions: `expected_reliability_reason` predicts or grades a route,
while `reason_code` names the concrete refusal that actually occurred.

| Value | Definition | Example |
|----|----|----|
| `satterthwaite_unavailable_at_boundary` | Satterthwaite inference is refused because the fit is boundary or reduced-rank. | `contrast(fit, L, method = "satterthwaite")` on a singular fit returns a not-assessed row. |
| `kenward_roger_unavailable_at_boundary` | Kenward-Roger inference is refused because the fit is boundary or reduced-rank. | `test_effect(fit, term, method = "kenward_roger")` on a singular fit returns a not-assessed row. |
| `bootstrap_lrt_requires_ml` | Bootstrap likelihood-ratio testing is refused because at least one fit was REML-fitted. | `compare(reduced, full, method = "bootstrap")` refuses unless the models are ML fits. |
| `bootstrap_lrt_reduced_formula_failed` | The reduced formula for a term bootstrap LRT could not be constructed. | A term name that cannot be dropped from the model formula produces this code. |
| `bootstrap_lrt_engine_refused` | R validation passed, but the Rust bootstrap-LRT bridge refused the payload. | A malformed or unsupported engine-side bootstrap LRT request is returned with this code. |
| `bootstrap_lrt_requires_nested_models` | Compared models are not ordered as a nested reduced and larger alternative pair. | A bootstrap LRT comparison with no degrees-of-freedom increase in the alternative fit refuses. |
| `bootstrap_lrt_requires_nested_model_frames` | The reduced model uses variables that are absent from the alternative model frame. | A comparison where the reduced formula references a variable not carried by the full fit refuses. |
| `bootstrap_lrt_requires_same_observations` | The compared fits do not share identical model-frame values for shared columns. | A bootstrap LRT comparison after different row filtering refuses. |
| `bootstrap_lrt_requires_same_weights` | The compared fits use different case weights. | A weighted reduced fit compared to an unweighted full fit refuses. |
| `bootstrap_cluster_resample_p_value_unavailable` | Cluster bootstrap is available as an estimator-distribution target but not as a fixed-effect p-value certificate. | `test_effect(..., method = "cluster_bootstrap", group = "subj")` returns a not-assessed p-value row. |
| `cluster_bootstrap_multifactor_ambiguous` | A cluster-bootstrap term test needs an explicit grouping factor when the model has multiple random-effect groups. | `test_effect(..., method = "cluster_bootstrap")` on a crossed-random-effect model refuses until `group` is supplied. |
| `profile_beta_unavailable_under_reml` | Fixed-effect beta profile intervals are absent by contract under REML. | `confint(fit, method = "profile")` on a REML fit adds beta rows with `NA` limits and this code. |
| `fixed_effect_inference_table_unavailable_legacy_object` | A legacy or revived object lacks the Rust fixed-effect inference artifact. | [`inference_table()`](https://bbuchsbaum.github.io/mixeff/reference/inference_table.md) on an older saved fit returns not-assessed rows with this code. |
