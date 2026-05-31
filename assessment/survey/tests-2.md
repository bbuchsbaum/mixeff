# Survey: tests-2 — test-inference.R, test-inference-options.R, test-inference-options-display.R

Generated: 2026-05-31

---

## File: test-inference.R

### Fixture

`mk_inference_fit()` — 9 subjects × 5 observations, `y ~ x + (1 | subject)`,
REML, seed 30. Interior (non-singular) fit. Used by every test in this file.

### Tests covered

| # | Test name | What it asserts |
|---|-----------|-----------------|
| 1 | `contrast() formats Rust fixed-effect contrast inference rows` | S3 class; estimate == fixef x; finite SE/statistic/p_value; method/requested_method/status/statistic_name/contrast_family payload |
| 2 | `contrast() preserves explicit Rust method outcomes` | All six method spellings (auto, kenward_roger, asymptotic, bootstrap-no-payload, none) round-trip through requested_method; statistic_name is "z" for asymptotic; bootstrap without nsim returns not_assessed+reason match; none returns not_computed+inference_not_requested |
| 3 | `contrast() can request Rust fixed-effect-null bootstrap rows` | nsim=30, seed=1; status available; finite p_value; replicate counts, failed_refit_policy, seed_rng, null_target.covariance_policy from detail payload |
| 4 | `contrast() preserves matrix rows, labels, and right-hand sides` | Named L matrix; rownames preserved; rhs subtracted from estimate; all-finite; method == "asymptotic_wald_z" for both rows |
| 5 | `estimability consumes the Rust fixed-contrast assessment` | S3 class; status "estimable"; estimable TRUE; rank/requested_rank 1; reason NA |
| 6 | `estimability reports a real engine status, never the unavailable placeholder` | Default L (identity); no "not_assessed" status; no "rust_estimability_certificate_unavailable" reason |
| 7 | `df_for_contrast pipes through the Rust inference-table df values` | satterthwaite: S3 class, finite, method attr; kenward_roger: finite, method attr; none: NA, method "not_requested", mm_unavailable_reason present |
| 8 | `test_effect() and single-model anova() consume Rust term rows` | KR test_effect: S3 class, term/method/requested_method/status/finite p_value/details; anova(bootstrap): NA p_value, not_assessed; anova(KR): available+finite; anova(fit,fit): mm_model_comparison; print not garbling list |
| 9 | `R inference surfaces preserve Rust detail payloads` | summary()==inference_table() details identity; bootstrap detail fields (target_kind, replicates, seed_rng, null_target, contrast_family); KR term details (contrast_family.family_id, restriction_rows, kenward_roger.restriction_rank); anova details == test_effect details |
| 10 | `Phase 3 prediction and covariance requests do not fabricate uncertainty` | vcov(theta): all NA + "theta_covariance_unavailable" attr; predict(se.fit=TRUE): NA se; predict(interval="confidence"): mm_inference_unavailable error |
| 11 | `confint(method = 'wald') is labelled as uncertified asymptotic output` | nrow == length(fixef); all finite; method attr "wald_asymptotic_from_stored_standard_errors"; status attr "not_certified_by_rust_inference_contract" |
| 12 | `confint(method = 'bootstrap') consumes full-model bootstrap intervals` | parm="x"; 1 row; finite; method attr "bootstrap_full_model_distribution"; interval "percentile"; status "available"; bootstrap payload fields (kind, replicates, seed_rng, replicate_statistics length); print contains "Bootstrap run:" and "Full bootstrap payload available"; no raw payload leakage in print |
| 13 | `inference_table() consumes Rust artifact rows when available` | S3 class; schema_name present; kind=="coefficient"; method in {satterthwaite,asymptotic_wald_z}; some available; finite p/stat for available; details column present |
| 14 | `summary() renders Rust coefficient inference rows` | S3 class; inference S3 class; requested_method "auto"; method in expected set; exactly one p-value column; some finite p |
| 15 | `summary() does not compute p-values missing from Rust rows` | Mutated artifact (NULL p_value, custom status/reason); NA in first p col; reason preserved |
| 16 | `saved fits preserve Rust artifact inference rows` | saveRDS/readRDS; inference_table identical |
| 17 | `saved fits preserve stored inference row details` | Detail injected into artifact; round-trip; contrast_family.family_id "c1" |
| 18 | `legacy fits without artifact inference table use unavailable fallback` | NULL inference table; all method "not_computed"; all status "not_assessed"; all reason "fixed_effect_inference_table_unavailable_legacy_object" |

### Tolerances / assertions used

- `expect_equal(ct$table$estimate, unname(fixef(fit)[["x"]]))` — default tolerance (`.Machine$double.eps^0.5`)
- `expect_equal(ct$table$estimate, as.numeric(L %*% fixef(fit)) - rhs)` — default tolerance
- All numeric checks via `is.finite()` / `is.na()` — no explicit absolute tolerance on inference quantities (SE, df, statistic, p)
- No lme4/lmerTest cross-validation tolerances in this file

### Skip directives

None. No `skip_on_cran`, no `skip_if_not_installed`.

---

## File: test-inference-options.R

### Fixture

No shared fixture function. Each test constructs its own data inline (subjects × days or subjects × conditions, seeded). Both REML and ML fits are used.

### Tests covered

| # | Test name | What it asserts |
|---|-----------|-----------------|
| 1 | `inference_options enumerates the six method routes` | S3 class; all six methods present (asymptotic_wald_z, satterthwaite, kenward_roger, bootstrap, bootstrap_lrt, cluster_bootstrap); expected_status / expected_reliability_reason / current columns present; no "recommended" column; exactly one current row |
| 2 | `inference_options marks satterthwaite/kenward_roger as not_assessed at boundary` | Singular fit; satt/KR expected_status "not_assessed"; satt reason "satterthwaite_unavailable_at_boundary"; wald expected_status "available" + current==TRUE |
| 3 | `summary() and inference_table() honor an explicit method on a singular fit` | auto → all asymptotic_wald_z available; explicit satterthwaite → all not_assessed + non-empty reason; inference_table(method="satterthwaite") same |
| 4 | `test_effect(method = 'bootstrap') works on a single-df term (singular fit)` | status available; reliability_reason in expected set; statistic_name "t"; finite p_value |
| 5 | `test_effect(method = 'bootstrap') produces a joint F test on a multi-df factor` | 3-level factor; status available; statistic_name "f"; num_df==2; finite p_value |
| 6 | `test_effect(method = 'bootstrap_lrt') refuses REML with a stable reason` | REML fit; status not_assessed; reason_code "bootstrap_lrt_requires_ml" |
| 7 | `test_effect(method = 'bootstrap_lrt') runs on ML fit and returns a chi-square row` | ML fit; status available; statistic_name "chi_square"; reliability "low"; reliability_reason "bootstrap_insufficient_replicates" at nsim=50; finite stat+p; bootstrap detail fields (successful_replicates, mcse, replicate_statistics length) |
| 8 | `inference_options mirrors bootstrap_lrt reliability threshold` | nsim=50 → "bootstrap_insufficient_replicates"; nsim=999 → "bootstrap_monte_carlo_replicates" |
| 9 | `test_effect(method = 'cluster_bootstrap') refuses p-values with a stable reason` | status not_assessed; NA p_value; reason_code "bootstrap_cluster_resample_p_value_unavailable"; detail fields (target_kind "cluster_resample", p_value_certified FALSE) |
| 10 | `test_effect(method = 'cluster_bootstrap') requires group for crossed models` | Two grouping factors, no group arg → not_assessed + "cluster_bootstrap_multifactor_ambiguous"; explicit bad group name → mm_arg_error |
| 11 | `inference_options rejects unknown terms` | term="definitely_not_a_term" → mm_arg_error |

### Tolerances / assertions used

- All numeric checks via `is.finite()` / `is.na()` — no absolute numeric tolerance
- No lme4 cross-validation in this file

### Skip directives

None.

---

## File: test-inference-options-display.R

### Fixture

`mm_inference_options_display_fit(REML=TRUE/FALSE)` — 10 subjects × 5 days,
`y ~ days + (1 | subj)`, seeded 31.

### Tests covered

| # | Test name | What it asserts |
|---|-----------|-----------------|
| 1 | `inference_options display columns are populated and readable` | display_status / display_reason / what_to_do_next columns present; all rows non-NA non-empty; no status/reliability_reason enum leakage; no snake_case enum pattern in display columns; profile_ci row present |
| 2 | `inference_options print uses display columns by default` | print output contains "runs now", "refused on this fit", "what_to_do_next", "raw enum columns" |
| 3 | `route-table refusal reasons match the verbs they advertise` | REML fit; cluster_bootstrap expected_status "not_assessed" + expected_reliability_reason == actual test_effect reason_code; bootstrap_lrt expected_status "not_assessed" + expected_reliability_reason == actual test_effect reason_code |
| 4 | `profile_ci route follows ML and REML profile contracts` | REML: expected_status "not_assessed" + reason "profile_beta_unavailable_under_reml"; ML: expected_status "available" + reason "profile_likelihood_ci" |

### Tolerances / assertions used

- No numeric assertions; purely structural/string checks
- `mm_raw_enum_like()` helper: `^[a-z][a-z0-9]*(_[a-z0-9]+)+$` used to detect leaked raw enum strings

### Skip directives

None.

---

## Cross-cutting observations

### What IS well-covered

1. **Rust wire contract integrity**: every method variant goes through; status/requested_method/reason round-trips are verified.
2. **Payload detail depth**: bootstrap detail fields (nsim, seed, null_target, target_kind, mcse, replicate_statistics), KR detail fields, estimability certificate — all probed.
3. **Boundary guard**: singular fit → satt/KR refused; auto resolves to wald; explicit satt on singular fit does not silently return the auto-cached row (regression guard).
4. **Inference options table shape**: six methods, current flag, no "recommended" column; display-column readability; route-table/execution consistency.
5. **Persistence**: saveRDS/readRDS roundtrip for inference rows and injected details.
6. **No-fabrication discipline**: theta_vcov returns NA, predict SE returns NA, predict interval throws.
7. **bootstrap_lrt ML vs REML gate and reliability grading** at the 999-replicate threshold.

### What is NOT covered (gaps)

See structured gaps below.

---

## Gaps (not tested but should be)

### G1 — `bootstrap_control()` argument validation
`bootstrap_control(nsim = 0)`, `bootstrap_control(nsim = -1)`, `bootstrap_control(seed = -1)`,
`bootstrap_control(nsim = NA)` should all throw `mm_arg_error`. Only the happy path is exercised.

### G2 — `confint()` argument validation
`confint(fit, level = 1.5)` and `confint(fit, level = 0)` should throw `mm_arg_error`. No test covers the level guard.

### G3 — `confint(method = "bootstrap", interval = "basic")` path
Only `interval = "percentile"` is exercised. The `"basic"` interval branch in `mm_select_bootstrap_interval` has no test.

### G4 — `confint()` with unknown `parm`
`confint(fit, parm = "nonexistent")` should throw `mm_arg_error`. Not covered.

### G5 — `confint()` with numeric `parm` index
The source accepts integer `parm` (translates to names). No test exercises that branch.

### G6 — `df_for_contrast()` on a boundary/singular fit
`df_for_contrast(singular_fit, L, method = "satterthwaite")` should return NA with a boundary reason. Only an interior fit is tested. The boundary path in `mm_boundary_df_method_unavailable` is exercised via `contrast()` and `test_effect()` but not via `df_for_contrast()` directly.

### G7 — `df_for_contrast()` with `method = "asymptotic"`
The source resolves asymptotic through the Rust bridge; the returned df value (Inf or a specific label) is untested.

### G8 — `test_effect()` with vector `term` (multiple terms in one call)
The function loops over `term` with `lapply`. No test passes a length-2 `term` vector to verify the loop produces a correctly combined table.

### G9 — `test_effect(method = "none")` path
`method = "none"` reaches `mm_unavailable_effect_table()`. No test covers the "none" method for `test_effect()` (only `contrast(..., method = "none")` is tested).

### G10 — `contrast(method = "boundary_lrt")` and `test_effect(method = "boundary_lrt")` produce correct refusal tables
These paths are exercised in `test-boundary-lrt.R` but not in `test-inference.R`. The inference.R file has no test that the `boundary_lrt` method returns `status = "unsupported"` with `reason_code = "boundary_lrt_not_applicable_to_fixed_effects"` from within this file's fixture.

### G11 — lme4 parity on inference quantities
None of the three files compare p-values, t/F statistics, or df against lme4/lmerTest reference values at the documented tolerances (fixef 1e-4, logLik 1e-3). The inference quantities are only checked for finiteness. A parity test pairing `summary(lmer(...), ddf="Satterthwaite")` against `summary(lmm(...), method="satterthwaite")` is absent.

### G12 — `inference_options()` for a boundary + ML fit (singular + REML=FALSE)
The boundary tests all use REML fits. A boundary ML fit would hit the path where `is_boundary=TRUE` and `is_reml=FALSE` simultaneously — e.g., profile_ci should be "not_assessed" with "profile_ci_unavailable_at_boundary" (not the REML reason). This combination is untested.

### G13 — `inference_options()` `approx_cost` column is present and non-empty
The table structure is checked but `approx_cost` content is never asserted. The `mm_inference_options_format_cost()` helper producing "~Xs @ nsim=N" vs "~N.Nmin @ nsim=N" output is untested.

### G14 — `inference_options()` `not_yet_wired` display path
`mm_inference_options_display_status("not_yet_wired")` returns "available upstream; R bridge pending". `mm_inference_options_next_step` has a branch for `not_yet_wired`. Neither branch is exercised by any test (no method currently returns that status in practice, but the code path is present).

### G15 — `estimability()` with a non-estimable contrast
Only estimable contrasts are tested. A contrast outside the column space should return `status = "not_estimable"` with a non-NA reason. Not covered.

### G16 — `summary(fit, tests = "coefficients", method = "bootstrap")` path
`summary()` with an explicit bootstrap method is not tested — only `method = "auto"` and mutated artifact (missing p_value) scenarios are covered.

### G17 — `inference_table(fit, method = "kenward_roger")` explicit KR path
`inference_table()` is only called with default method ("auto") or `method = "satterthwaite"` (in test-inference-options.R). An explicit `method = "kenward_roger"` call that returns available KR rows is not tested through `inference_table()` directly.

### G18 — `bootstrap_control(failed_refit_policy = "count_extreme")` and `"abort"`
Only `"exclude"` (the default) appears in tests. The `count_extreme` and `abort` policies are constructed in `bootstrap_control()` but never passed through a live inference call to verify the policy appears in the returned detail payload.
