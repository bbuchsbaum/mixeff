use extendr_api::prelude::*;
use mixeff_rs::compiler::{
    compile_formula_ir, CompiledModelArtifact, ContrastMatrix, ContrastRhs, FixedEffectHypothesis,
    FixedEffectInferenceRowKind, FixedEffectInferenceTable, FixedEffectTestMethod,
    FIXED_EFFECT_INFERENCE_TABLE_SCHEMA, FIXED_EFFECT_INFERENCE_TABLE_SCHEMA_VERSION,
};
use mixeff_rs::formula::parse_formula;
use mixeff_rs::model::{
    parametricbootstrap, BootstrapFailedRefitPolicy, BootstrapRefitOptions, BootstrapReplicate,
    BootstrapSeedRecord, BootstrapTarget, FixedEffectBootstrapOptions, LinearMixedModel,
    MixedModelBootstrap, MixedModelFit,
};
use nalgebra::{DMatrix, DVector};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};

mod data;

// ---------------------------------------------------------------------------
// Schema versions
//
// Per PRD §5.2 item 3, every artifact type that crosses the wire should
// expose a `pub const SCHEMA_VERSION` upstream so `mm_json_negotiate()` can
// fail fast on mismatch. Upstream now exposes these as `pub const`s in
// `mixeff_rs::compiler` (e.g. `COMPILED_ARTIFACT_SCHEMA_VERSION` is `1`).
// The wrapper-side schema-negotiator stores versions as strings, so each
// upstream `u32` is mirrored here as the matching string literal. Drift
// is caught by `mm_compile_audit_smoke` in the test suite, which reads the
// constant via the FFI manifest.
//
// Schema names are stable identifiers. Bumping a version is a breaking
// change negotiated by `mm_json_negotiate()` — never silently.

const SCHEMA_VERSION_FORMULA: &str = "v0";
const SCHEMA_VERSION_COMPILED_ARTIFACT: &str = "1";
const SCHEMA_VERSION_MODEL_AUDIT_REPORT: &str = "2";
const SCHEMA_VERSION_RANDOM_TERM_CARD: &str = "1";
const SCHEMA_VERSION_FIXED_EFFECT_INFERENCE_TABLE: &str =
    FIXED_EFFECT_INFERENCE_TABLE_SCHEMA_VERSION;
const SCHEMA_NAME_MARGINAL_QUANTITY_TABLE: &str = "mixedmodels.marginal_quantity_table";
const SCHEMA_VERSION_MARGINAL_QUANTITY_TABLE: &str = "1.0.0";
const MIXEDMODELS_CRATE_VERSION: &str = "0.1.0";

// Closed list of (schema_name, schema_version) the wrapper currently
// knows how to consume. Future phases extend this when artifacts gain
// schemas (compiled-model artifact, audit, theta_map, certificate,
// inference, reproducibility, prediction).
const KNOWN_SCHEMAS: &[(&str, &str)] = &[
    ("formula", SCHEMA_VERSION_FORMULA),
    (
        "mixedmodels.compiled_model_artifact",
        SCHEMA_VERSION_COMPILED_ARTIFACT,
    ),
    (
        "mixedmodels.model_audit_report",
        SCHEMA_VERSION_MODEL_AUDIT_REPORT,
    ),
    (
        "mixedmodels.random_term_card",
        SCHEMA_VERSION_RANDOM_TERM_CARD,
    ),
    (
        FIXED_EFFECT_INFERENCE_TABLE_SCHEMA,
        SCHEMA_VERSION_FIXED_EFFECT_INFERENCE_TABLE,
    ),
    (
        SCHEMA_NAME_MARGINAL_QUANTITY_TABLE,
        SCHEMA_VERSION_MARGINAL_QUANTITY_TABLE,
    ),
];

// ---------------------------------------------------------------------------
// Interrupt bridge
//
// extendr_api 0.9 does not wrap R_CheckUserInterrupt; declare the binding
// directly. Calling it inside a long-running Rust loop gives R a chance to
// terminate the call when the user presses Ctrl-C. R_CheckUserInterrupt
// longjmps out to R's error handler when a pending interrupt is observed,
// so callers do not need to handle a return value.
//
// Phase 0 ships only the bridge primitive plus a tiny demo (`mm_interrupt_demo`).
// Phase 1+ will use the same hook from inside `lmm()` / `glmm()` fit loops
// so that long PLS / PIRLS optimizations remain interruptible.

extern "C" {
    fn R_CheckUserInterrupt();
}

#[inline]
fn check_user_interrupt() {
    // SAFETY: R_CheckUserInterrupt is a stable, side-effect-only R API
    // function; the only thing it does on a pending interrupt is
    // longjmp out, which unwinds the Rust stack on the way to R's
    // error handler. No locals to drop because we're at a clean
    // suspension point.
    unsafe { R_CheckUserInterrupt() };
}

// ---------------------------------------------------------------------------
// Bridge functions exposed to R
// ---------------------------------------------------------------------------

/// Parse an lme4-style formula string and return its canonical form.
///
/// On success, returns the formula's `Display` rendering as a single
/// character string. On parse failure, throws an extendr error whose message
/// is prefixed with `mm_formula_error:` so the R wrapper can convert it to a
/// typed `mm_formula_error` condition.
///
/// @noRd
#[extendr]
fn mm_parse_formula(formula: &str) -> std::result::Result<String, String> {
    parse_formula(formula)
        .map(|f| format!("{}", f))
        .map_err(|e| format!("mm_formula_error: {}", e))
}

/// Return the package's formula manifest.
///
/// A small, versioned record of what the wrapper supports today. R callers
/// use this to gate behavior on capabilities rather than hard-coding
/// assumptions, and to record provenance on `mm_fit` artifacts.
///
/// @noRd
#[extendr]
fn mm_formula_manifest() -> List {
    let operators: Robj = vec!["+", "-", "*", ":", "/", "&", "|", "||"].into();
    let intercept_forms: Robj = vec!["1", "0", "-1"].into();
    let random_term_forms: Robj = vec![
        "(1 | g)",
        "(0 + x | g)",
        "(1 + x | g)",
        "(1 + x || g)",
        "(1 | g) + (0 + x | g)",
        "(1 | a/b)",
        "(1 | g:h)",
        "(1 | g) + (1 | h)",
    ]
    .into();
    let transformations: Robj = vec![
        "implicit_intercept",
        "nested_grouping_expansion",
        "interaction_grouping",
    ]
    .into();

    let formula_features = list!(
        operators = operators,
        intercept_forms = intercept_forms,
        random_term_forms = random_term_forms,
        transformations = transformations,
    );

    // Manifest keys mirror the negotiator's `KNOWN_SCHEMAS` entries
    // exactly; the agreement test in `test-schema-versioning.R` checks
    // every manifest key resolves to the same string in
    // `mm_json_known_schemas()`. The compiled-artifact schema name
    // contains a dot (`mixedmodels.compiled_model_artifact`) which
    // extendr's `list!` macro cannot accept as an identifier-style key,
    // so we build this inner list via `List::from_pairs`.
    let schema_versions = List::from_pairs([
        ("formula", Robj::from(SCHEMA_VERSION_FORMULA)),
        (
            "mixedmodels.compiled_model_artifact",
            Robj::from(SCHEMA_VERSION_COMPILED_ARTIFACT),
        ),
        (
            "mixedmodels.model_audit_report",
            Robj::from(SCHEMA_VERSION_MODEL_AUDIT_REPORT),
        ),
        (
            "mixedmodels.random_term_card",
            Robj::from(SCHEMA_VERSION_RANDOM_TERM_CARD),
        ),
        (
            FIXED_EFFECT_INFERENCE_TABLE_SCHEMA,
            Robj::from(SCHEMA_VERSION_FIXED_EFFECT_INFERENCE_TABLE),
        ),
        (
            SCHEMA_NAME_MARGINAL_QUANTITY_TABLE,
            Robj::from(SCHEMA_VERSION_MARGINAL_QUANTITY_TABLE),
        ),
    ]);

    let capabilities = list!(
        parse_formula = TRUE,
        compile_model = TRUE,
        audit_design = TRUE,
        explain_model = TRUE,
        random_options = TRUE,
        compare_covariance = TRUE,
        fit_lmm = TRUE,
        audit = TRUE,
        changes = TRUE,
        diagnostics = TRUE,
        fit_status = TRUE,
        parameterization = TRUE,
        roles = TRUE,
        as_json = TRUE,
        fit_glmm = FALSE,
        simulate = TRUE,
        inference = TRUE,
        fixed_effect_inference_table = TRUE,
        satterthwaite = TRUE,
        kenward_roger_explicit = TRUE,
        bootstrap_fixed_effect_payload = TRUE,
        marginal_quantity_table = TRUE,
        marginal_quantities = TRUE,
    );

    list!(
        mixeff_rust_version = env!("CARGO_PKG_VERSION"),
        crate_version = MIXEDMODELS_CRATE_VERSION,
        schema_versions = schema_versions,
        formula_features = formula_features,
        capabilities = capabilities,
    )
}

/// Negotiate a single schema declaration against the wrapper's known set.
///
/// Returns `TRUE` on a clean match. On mismatch or unknown schema, throws
/// an extendr error prefixed with `mm_schema_error:` for routing to the
/// typed condition on the R side. The full negotiation logic (multi-schema
/// headers, version-range comparisons) lives in R; this function is a
/// pin-point primitive over the closed `KNOWN_SCHEMAS` table.
///
/// @noRd
#[extendr]
fn mm_json_negotiate_one(name: &str, version: &str) -> std::result::Result<bool, String> {
    for (n, v) in KNOWN_SCHEMAS {
        if *n == name {
            if *v == version {
                return Ok(true);
            }
            return Err(format!(
                "mm_schema_error: schema '{}' version mismatch (wrapper expects '{}', got '{}')",
                name, v, version
            ));
        }
    }
    Err(format!("mm_schema_error: unknown schema '{}'", name))
}

/// Return the closed list of (schema_name, schema_version) the wrapper
/// currently knows how to consume.
///
/// @noRd
#[extendr]
fn mm_json_known_schemas() -> List {
    let names: Robj = KNOWN_SCHEMAS
        .iter()
        .map(|(n, _)| *n)
        .collect::<Vec<_>>()
        .into();
    let versions: Robj = KNOWN_SCHEMAS
        .iter()
        .map(|(_, v)| *v)
        .collect::<Vec<_>>()
        .into();
    list!(name = names, version = versions)
}

/// Compile a formula against a data manifest and return the artifact JSON.
///
/// Runs the upstream `parse_formula` -> `compile_formula_ir` ->
/// `CompiledModelArtifact::new` -> `attach_design_audit` pipeline and
/// returns the serialized `CompiledModelArtifact` (schema
/// `mixedmodels.compiled_model_artifact` v1). No fit happens here — that
/// lands in Phase 1.E.
///
/// Inputs follow the wire format described in `data.rs`. The R wrapper
/// (`compile_model()`) is the only supported caller; it validates the
/// shape before invoking this primitive, so most error paths surface as
/// `mm_data_error:` / `mm_formula_error:` / `mm_schema_error:` tagged
/// strings the R side routes to typed conditions.
///
/// @noRd
#[extendr]
fn mm_compile_model_json(
    formula: &str,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
) -> std::result::Result<String, String> {
    let parsed = parse_formula(formula).map_err(|e| format!("mm_formula_error: {}", e))?;
    let semantic = compile_formula_ir(&parsed);
    let df = data::build_dataframe(
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &column_order,
    )?;

    let mut artifact = CompiledModelArtifact::new(parsed.to_string(), semantic);
    artifact.attach_design_audit(&df);

    serde_json::to_string(&artifact)
        .map_err(|e| format!("mm_schema_error: failed to serialize artifact: {}", e))
}

/// Fit a linear mixed-effects model and return the fit payload JSON.
///
/// This is the Phase 1.E fit primitive: parse formula, build the upstream
/// `DataFrame`, construct `LinearMixedModel`, run `.fit(reml)`, then return
/// the post-fit `CompiledModelArtifact` plus flat extractor-friendly numeric
/// duplicates. R owns the user-facing S3 surface; Rust owns the numerical fit
/// and the compiler/audit artifact.
///
/// @noRd
#[extendr]
fn mm_fit_lmm_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    weights: Doubles,
    control_json: &str,
) -> std::result::Result<String, String> {
    let _control: Value = serde_json::from_str(control_json)
        .map_err(|e| format!("mm_fit_error: invalid control JSON: {}", e))?;
    let parsed = parse_formula(formula).map_err(|e| format!("mm_formula_error: {}", e))?;
    let df = data::build_dataframe(
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &column_order,
    )?;
    let weights = optional_case_weights(&weights, df.nrow())?;

    let mut model = LinearMixedModel::new(parsed, &df, weights.as_deref())
        .map_err(|e| format!("mm_fit_error: failed to construct LMM: {}", e))?;
    model
        .fit(reml)
        .map_err(|e| format!("mm_fit_error: failed to fit LMM: {}", e))?;

    let artifact_json = serde_json::to_string(model.compiler_artifact()).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize post-fit artifact: {}",
            e
        )
    })?;
    let artifact_value: Value = serde_json::from_str(&artifact_json).map_err(|e| {
        format!(
            "mm_schema_error: failed to inspect post-fit artifact JSON: {}",
            e
        )
    })?;
    let fit_status = artifact_value
        .get("optimizer_certificate")
        .and_then(|x| x.get("status"))
        .and_then(Value::as_str)
        .unwrap_or("not_assessed");

    let beta = model.coef();
    let beta_names = model.coef_names();
    let std_errors = model.stderror();
    let fixed_fitted = model.fixed_effect_fitted();
    let log_likelihood = model.loglikelihood();
    let dof = model.dof();
    let nobs = model.nobs();
    let df_residual = nobs.saturating_sub(dof);

    let opt = model.opt_summary();
    let payload = json!({
        "schema": {
            "schema_name": "mixeff.lmm_fit_result",
            "schema_version": 1
        },
        "artifact_json": artifact_json,
        "formula": model.formula.to_string(),
        "reml": reml,
        "beta": beta.iter().copied().collect::<Vec<_>>(),
        "beta_names": beta_names,
        "theta": model.theta(),
        "sigma": model.sigma(),
        "log_likelihood": log_likelihood,
        "deviance": -2.0 * log_likelihood,
        "aic": model.aic(),
        "bic": model.bic(),
        "nobs": nobs,
        "dof": dof,
        "df_residual": df_residual,
        "fit_status": fit_status,
        "std_errors": std_errors.iter().copied().collect::<Vec<_>>(),
        "fixed_fitted": fixed_fitted.iter().copied().collect::<Vec<_>>(),
        "fitted": model.fitted().iter().copied().collect::<Vec<_>>(),
        "residuals": model.residuals().iter().copied().collect::<Vec<_>>(),
        "ranef": random_effects_json(&model),
        "varcorr": varcorr_json(&model),
        "optimizer": {
            "backend": opt.backend_name(),
            "algorithm": opt.optimizer_name(),
            "return_value": opt.return_value.as_str(),
            "function_evaluations": opt.feval,
            "objective": opt.fmin,
            "reml": opt.reml
        }
    });

    serde_json::to_string(&payload)
        .map_err(|e| format!("mm_schema_error: failed to serialize fit result: {}", e))
}

/// Evaluate fixed-effect contrast rows through the Rust inference contract.
///
/// @noRd
#[extendr]
fn mm_fixed_effect_contrast_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    weights: Doubles,
    control_json: &str,
    l_values: Doubles,
    nrow: i32,
    ncol: i32,
    labels: Strings,
    rhs: Doubles,
    method: &str,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &weights,
        control_json,
    )?;
    let hypotheses = fixed_effect_hypotheses(l_values, nrow, ncol, labels, rhs)?;
    let method = fixed_effect_test_method(method)?;

    serde_json::to_string(&model.fixed_effect_contrast_inference_table(hypotheses, method))
        .map_err(|e| format!("mm_schema_error: failed to serialize contrast table: {}", e))
}

/// Evaluate fixed-effect-null bootstrap contrast rows through Rust.
///
/// @noRd
#[extendr]
fn mm_fixed_effect_bootstrap_contrast_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    weights: Doubles,
    control_json: &str,
    l_values: Doubles,
    nrow: i32,
    ncol: i32,
    labels: Strings,
    rhs: Doubles,
    bootstrap_options_json: &str,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &weights,
        control_json,
    )?;
    let hypotheses = fixed_effect_hypotheses(l_values, nrow, ncol, labels, rhs)?;
    let options = fixed_effect_bootstrap_options(bootstrap_options_json)?;

    serde_json::to_string(&model.fixed_effect_null_bootstrap_inference_table(hypotheses, options))
        .map_err(|e| {
            format!(
                "mm_schema_error: failed to serialize bootstrap contrast table: {}",
                e
            )
        })
}

/// Evaluate a full-model bootstrap contrast payload for fixed-effect
/// confidence intervals. This target simulates from the fitted model and
/// returns interval summaries plus replicate accounting; it does not certify
/// fixed-effect hypothesis-test p-values.
///
/// @noRd
#[extendr]
fn mm_full_model_bootstrap_contrast_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    weights: Doubles,
    control_json: &str,
    l_values: Doubles,
    nrow: i32,
    ncol: i32,
    labels: Strings,
    rhs: Doubles,
    bootstrap_options_json: &str,
    levels: Doubles,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &weights,
        control_json,
    )?;
    let mut hypotheses = fixed_effect_hypotheses(l_values, nrow, ncol, labels, rhs)?;
    if hypotheses.len() != 1 {
        return Err(
            "mm_inference_unavailable: full-model bootstrap intervals are currently certified only for scalar contrasts"
                .to_string(),
        );
    }
    let hypothesis = hypotheses.remove(0);
    let options = fixed_effect_bootstrap_options(bootstrap_options_json)?;
    let level_values: Vec<f64> = levels.iter().map(|v| v.0).collect();

    // Single-row contrast: extract L row and observed contrast estimate (L * beta).
    let l_row: Vec<f64> = hypothesis
        .l
        .values
        .row(0)
        .iter()
        .copied()
        .collect();
    let contrast_label = hypothesis.label.clone();
    let beta_observed = model.beta();
    let observed_estimate: f64 = l_row
        .iter()
        .zip(beta_observed.iter())
        .map(|(l, b)| l * b)
        .sum();

    let (mut rng, seed_record) = make_bootstrap_rng(options.seed);
    let bsamp = parametricbootstrap(&mut rng, options.requested_replicates, &model);

    // Contrast statistic per replicate: L * beta_b.
    let replicate_stats: Vec<f64> = bsamp
        .fits
        .iter()
        .map(|fit| {
            if fit.beta.iter().any(|x| !x.is_finite()) {
                f64::NAN
            } else {
                l_row.iter().zip(fit.beta.iter()).map(|(l, b)| l * b).sum()
            }
        })
        .collect();
    let mut finite_stats: Vec<f64> = replicate_stats
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    if finite_stats.is_empty() {
        return Err(
            "mm_inference_unavailable: bootstrap produced no finite contrast statistics".to_string(),
        );
    }
    finite_stats.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n_finite = finite_stats.len();

    let mut intervals = Vec::with_capacity(level_values.len() * 2);
    for &level in &level_values {
        if !(0.0 < level && level < 1.0) {
            return Err(format!(
                "mm_inference_unavailable: bootstrap level must be in (0, 1); got {level}"
            ));
        }
        let alpha = (1.0 - level) / 2.0;
        let lower_q = quantile_sorted(&finite_stats, alpha);
        let upper_q = quantile_sorted(&finite_stats, 1.0 - alpha);
        intervals.push(json!({
            "method": "percentile",
            "level": level,
            "lower": lower_q,
            "upper": upper_q,
            "n": n_finite,
        }));
        intervals.push(json!({
            "method": "basic",
            "level": level,
            "lower": 2.0 * observed_estimate - upper_q,
            "upper": 2.0 * observed_estimate - lower_q,
            "n": n_finite,
        }));
    }

    let refit_options = BootstrapRefitOptions::from_model(&model);
    let metadata = bsamp.run_metadata_for_model(
        &model,
        BootstrapTarget::full_model_distribution(format!("contrast: {}", contrast_label)),
        options.requested_replicates,
        options.failed_refit_policy,
        seed_record,
        refit_options,
        Some(format!("contrast: {}", contrast_label)),
        Some(&replicate_stats),
        None,
    );

    let payload = json!({
        "intervals": intervals,
        "metadata": metadata,
        "replicate_statistics": replicate_stats,
        "observed_estimate": observed_estimate,
        "contrast_label": contrast_label,
    });

    serde_json::to_string(&payload).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize full-model bootstrap contrast payload: {}",
            e
        )
    })
}

/// Evaluate a fixed-effect-null bootstrap *term* row (joint Wald/F over an
/// arbitrary L matrix) through Rust. Single-df hypotheses produce a t-form
/// row; multi-df hypotheses produce an F-form row with `numerator_df` set
/// to the effective restriction rank. Wraps
/// `LinearMixedModel::fixed_effect_null_bootstrap_inference_row` with
/// `kind = Term`.
///
/// @noRd
#[extendr]
fn mm_fixed_effect_bootstrap_term_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    weights: Doubles,
    control_json: &str,
    l_values: Doubles,
    nrow: i32,
    ncol: i32,
    label: &str,
    rhs: Doubles,
    bootstrap_options_json: &str,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &weights,
        control_json,
    )?;
    let nrow_us = usize::try_from(nrow)
        .map_err(|_| "mm_inference_unavailable: contrast row count must be non-negative")?;
    let ncol_us = usize::try_from(ncol)
        .map_err(|_| "mm_inference_unavailable: contrast column count must be non-negative")?;
    let l_vec = l_values.iter().map(|v| v.0).collect::<Vec<_>>();
    let rhs_vec = rhs.iter().map(|v| v.0).collect::<Vec<_>>();
    if l_vec.len() != nrow_us * ncol_us {
        return Err(format!(
            "mm_inference_unavailable: term contrast has {} value(s), expected {}",
            l_vec.len(),
            nrow_us * ncol_us
        ));
    }
    if rhs_vec.len() != nrow_us {
        return Err(format!(
            "mm_inference_unavailable: term contrast rhs has length {}, expected {}",
            rhs_vec.len(),
            nrow_us
        ));
    }
    let l_matrix = DMatrix::from_row_slice(nrow_us, ncol_us, &l_vec);
    let contrast = ContrastMatrix::new(l_matrix)
        .map_err(|e| format!("mm_inference_unavailable: invalid contrast matrix: {}", e))?;
    let rhs_struct = ContrastRhs::new(DVector::from_vec(rhs_vec))
        .map_err(|e| format!("mm_inference_unavailable: invalid rhs: {}", e))?;
    let hypothesis =
        FixedEffectHypothesis::new(label.to_string(), contrast, rhs_struct).map_err(|e| {
            format!(
                "mm_inference_unavailable: invalid term hypothesis '{}': {}",
                label, e
            )
        })?;
    let options = fixed_effect_bootstrap_options(bootstrap_options_json)?;

    let row = model.fixed_effect_null_bootstrap_inference_row(
        FixedEffectInferenceRowKind::Term,
        hypothesis,
        &options,
    );
    let table = FixedEffectInferenceTable::new(vec![row]);
    serde_json::to_string(&table).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize bootstrap term table: {}",
            e
        )
    })
}

/// Evaluate fixed-effect term rows through Rust-owned term hypotheses.
///
/// @noRd
#[extendr]
fn mm_fixed_effect_term_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    weights: Doubles,
    control_json: &str,
    method: &str,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &weights,
        control_json,
    )?;
    let method = fixed_effect_test_method(method)?;

    serde_json::to_string(&model.fixed_effect_term_inference_table(method))
        .map_err(|e| format!("mm_schema_error: failed to serialize term table: {}", e))
}

/// Run a parametric bootstrap likelihood-ratio test between two LMMs fitted
/// to the same data. Both models must be ML (the upstream engine refuses
/// REML). Returns a serialised `BootstrapLikelihoodRatioTest`.
///
/// @noRd
#[extendr]
fn mm_bootstrap_lrt_json(
    reduced_formula: &str,
    alternative_formula: &str,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    weights: Doubles,
    control_json: &str,
    bootstrap_options_json: &str,
) -> std::result::Result<String, String> {
    let reduced = fit_lmm_from_bridge_data(
        reduced_formula,
        false,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &weights,
        control_json,
    )?;
    let alternative = fit_lmm_from_bridge_data(
        alternative_formula,
        false,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &weights,
        control_json,
    )?;
    let options = fixed_effect_bootstrap_options(bootstrap_options_json)?;

    let observed_logl_red = reduced.loglikelihood();
    let observed_logl_alt = alternative.loglikelihood();
    let observed_lrt = 2.0 * (observed_logl_alt - observed_logl_red);
    if !observed_lrt.is_finite() {
        return Err(
            "mm_inference_unavailable: observed LRT statistic is not finite".to_string(),
        );
    }

    let (mut rng, seed_record) = make_bootstrap_rng(options.seed);

    // Custom bootstrap loop: simulate from the reduced (null) model, refit
    // both reduced and alternative on the simulated response, record the LRT.
    // We still build a MixedModelBootstrap recording the *reduced* refit
    // history so the existing run_metadata_for_model() telemetry path stays
    // intact (boundary counts, failed-refits, mcse).
    let mut fits: Vec<BootstrapReplicate> = Vec::with_capacity(options.requested_replicates);
    let mut replicate_stats: Vec<f64> = Vec::with_capacity(options.requested_replicates);

    for _ in 0..options.requested_replicates {
        let y_sim = reduced.simulate(&mut rng);
        let mut work_red = reduced.clone();
        let mut work_alt = alternative.clone();
        let refit_red = work_red.refit(y_sim.as_slice());
        let refit_alt = work_alt.refit(y_sim.as_slice());

        let stat = match (&refit_red, &refit_alt) {
            (Ok(()), Ok(())) => {
                let s = 2.0 * (work_alt.loglikelihood() - work_red.loglikelihood());
                if s.is_finite() {
                    s
                } else {
                    f64::NAN
                }
            }
            _ => f64::NAN,
        };
        replicate_stats.push(stat);

        let beta = work_red.beta();
        if refit_red.is_ok() {
            fits.push(BootstrapReplicate {
                objective: work_red.objective(),
                sigma: work_red.sigma(),
                beta,
                se: work_red.stderror(),
                theta: work_red.theta(),
            });
        } else {
            let n_beta = beta.len();
            fits.push(BootstrapReplicate {
                objective: f64::NAN,
                sigma: f64::NAN,
                se: DVector::from_element(n_beta, f64::NAN),
                beta,
                theta: work_red.theta(),
            });
        }
    }

    let bsamp = MixedModelBootstrap { fits };
    let successful: usize = replicate_stats.iter().filter(|v| v.is_finite()).count();
    let p_value = if successful > 0 {
        let count_ge = replicate_stats
            .iter()
            .filter(|v| v.is_finite() && **v >= observed_lrt)
            .count();
        Some((count_ge as f64) / (successful as f64))
    } else {
        None
    };
    let mcse = p_value.map(|p| (p * (1.0 - p) / successful as f64).sqrt());

    let refit_options = BootstrapRefitOptions::from_model(&reduced);
    let metadata = bsamp.run_metadata_for_model(
        &reduced,
        BootstrapTarget::fixed_effect_null(
            "bootstrap likelihood-ratio test",
            "alternative vs. reduced",
        ),
        options.requested_replicates,
        options.failed_refit_policy,
        seed_record,
        refit_options,
        Some("lrt_chi_square".to_string()),
        Some(&replicate_stats),
        p_value,
    );

    let payload = json!({
        "observed_statistic": observed_lrt,
        "p_value": p_value,
        "mcse": mcse,
        "notes": metadata.notes,
        "payload": {
            "metadata": metadata,
            "replicate_statistics": replicate_stats,
        },
    });

    serde_json::to_string(&payload)
        .map_err(|e| format!("mm_schema_error: failed to serialize bootstrap LRT: {}", e))
}

fn fit_lmm_from_bridge_data(
    formula: &str,
    reml: bool,
    column_order: &Strings,
    numeric_columns: &List,
    categorical_values: &List,
    categorical_levels: &List,
    weights: &Doubles,
    control_json: &str,
) -> std::result::Result<LinearMixedModel, String> {
    let _control: Value = serde_json::from_str(control_json)
        .map_err(|e| format!("mm_fit_error: invalid control JSON: {}", e))?;
    let parsed = parse_formula(formula).map_err(|e| format!("mm_formula_error: {}", e))?;
    let df = data::build_dataframe(
        numeric_columns,
        categorical_values,
        categorical_levels,
        column_order,
    )?;
    let weights = optional_case_weights(weights, df.nrow())?;
    let mut model = LinearMixedModel::new(parsed, &df, weights.as_deref())
        .map_err(|e| format!("mm_fit_error: failed to construct LMM: {}", e))?;
    model
        .fit(reml)
        .map_err(|e| format!("mm_fit_error: failed to fit LMM: {}", e))?;
    Ok(model)
}

fn optional_case_weights(
    weights: &Doubles,
    n_obs: usize,
) -> std::result::Result<Option<Vec<f64>>, String> {
    let values = weights.iter().map(|v| v.0).collect::<Vec<_>>();
    if values.is_empty() {
        return Ok(None);
    }
    if values.len() != n_obs {
        return Err(format!(
            "mm_data_error: weights length ({}) does not match number of observations ({n_obs})",
            values.len()
        ));
    }
    for (i, w) in values.iter().enumerate() {
        if !w.is_finite() || *w <= 0.0 {
            return Err(format!(
                "mm_data_error: weight at index {i} must be finite and positive (got {w})"
            ));
        }
    }
    Ok(Some(values))
}

fn fixed_effect_hypotheses(
    l_values: Doubles,
    nrow: i32,
    ncol: i32,
    labels: Strings,
    rhs: Doubles,
) -> std::result::Result<Vec<FixedEffectHypothesis>, String> {
    let nrow = usize::try_from(nrow)
        .map_err(|_| "mm_inference_unavailable: contrast row count must be non-negative")?;
    let ncol = usize::try_from(ncol)
        .map_err(|_| "mm_inference_unavailable: contrast column count must be non-negative")?;
    let l_values = l_values.iter().map(|value| value.0).collect::<Vec<_>>();
    let rhs = rhs.iter().map(|value| value.0).collect::<Vec<_>>();
    let labels = labels
        .iter()
        .map(|label| {
            let s: &str = label.as_ref();
            s.to_string()
        })
        .collect::<Vec<_>>();

    if l_values.len() != nrow * ncol {
        return Err(format!(
            "mm_inference_unavailable: contrast matrix has {} value(s), expected {}",
            l_values.len(),
            nrow * ncol
        ));
    }
    if rhs.len() != nrow {
        return Err(format!(
            "mm_inference_unavailable: contrast rhs has length {}, expected {nrow}",
            rhs.len()
        ));
    }
    if labels.len() != nrow {
        return Err(format!(
            "mm_inference_unavailable: contrast labels have length {}, expected {nrow}",
            labels.len()
        ));
    }

    (0..nrow)
        .map(|i| {
            let start = i * ncol;
            let l = ContrastMatrix::new(DMatrix::from_row_slice(
                1,
                ncol,
                &l_values[start..start + ncol],
            ))
            .map_err(|e| format!("mm_inference_unavailable: {e}"))?;
            let rhs = ContrastRhs::new(DVector::from_row_slice(&[rhs[i]]))
                .map_err(|e| format!("mm_inference_unavailable: {e}"))?;
            FixedEffectHypothesis::new(labels[i].clone(), l, rhs)
                .map_err(|e| format!("mm_inference_unavailable: {e}"))
        })
        .collect()
}

fn fixed_effect_test_method(method: &str) -> std::result::Result<FixedEffectTestMethod, String> {
    match method {
        "auto" => Ok(FixedEffectTestMethod::Auto),
        "asymptotic" | "asymptotic_wald_z" => Ok(FixedEffectTestMethod::AsymptoticWaldZ),
        "satterthwaite" => Ok(FixedEffectTestMethod::Satterthwaite),
        "kenward_roger" => Ok(FixedEffectTestMethod::KenwardRoger),
        "bootstrap" => Ok(FixedEffectTestMethod::ParametricBootstrap),
        other => Err(format!(
            "mm_inference_unavailable: unsupported fixed-effect inference method `{other}`"
        )),
    }
}

fn fixed_effect_bootstrap_options(
    options_json: &str,
) -> std::result::Result<FixedEffectBootstrapOptions, String> {
    let value: Value = serde_json::from_str(options_json).map_err(|e| {
        format!(
            "mm_inference_unavailable: invalid bootstrap options JSON: {}",
            e
        )
    })?;
    let requested_replicates = value
        .get("requested_replicates")
        .or_else(|| value.get("nsim"))
        .and_then(Value::as_u64)
        .unwrap_or(999) as usize;
    let seed = value.get("seed").and_then(Value::as_u64);
    let failed_refit_policy = match value
        .get("failed_refit_policy")
        .or_else(|| value.get("failed_refit"))
        .and_then(Value::as_str)
        .unwrap_or("exclude")
    {
        "exclude" => BootstrapFailedRefitPolicy::Exclude,
        "count_extreme" => BootstrapFailedRefitPolicy::CountExtreme,
        "abort" => BootstrapFailedRefitPolicy::Abort,
        other => {
            return Err(format!(
                "mm_inference_unavailable: unsupported bootstrap failed-refit policy `{other}`"
            ))
        }
    };
    Ok(FixedEffectBootstrapOptions {
        requested_replicates,
        failed_refit_policy,
        seed,
    })
}

fn random_effects_json(model: &LinearMixedModel) -> Value {
    let effects = model.ranef_b();
    let terms = model
        .reterms
        .iter()
        .zip(effects.iter())
        .map(|(rt, b)| {
            let values = (0..rt.levels.len())
                .map(|level_idx| {
                    (0..rt.vsize)
                        .map(|coef_idx| b[(coef_idx, level_idx)])
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            json!({
                "group": rt.grouping_name.as_str(),
                "levels": rt.levels.clone(),
                "names": rt.cnames.clone(),
                "values": values
            })
        })
        .collect::<Vec<_>>();
    json!(terms)
}

fn varcorr_json(model: &LinearMixedModel) -> Value {
    let vc = model.varcorr();
    let components = vc
        .components
        .iter()
        .map(|component| {
            json!({
                "group": component.group.as_str(),
                "names": component.names.clone(),
                "std_dev": component.std_dev.clone(),
                "correlations": component.correlations.clone()
            })
        })
        .collect::<Vec<_>>();
    json!({
        "components": components,
        "residual_sd": vc.residual_sd
    })
}

/// Render an `audit_design()` artifact as text.
///
/// Takes the JSON produced by `mm_compile_model_json` and returns the
/// upstream `ModelAuditReport::Display` rendering — the source of truth
/// for audit wording. R's `print.mm_audit()` calls this and emits the
/// result with `cat()`. Routing all English wording through this entry
/// point keeps the R9 "no advice creep" contract enforceable: reviewers
/// see drift in one place (the Rust crate), not scattered across R
/// formatters.
///
/// @noRd
#[extendr]
fn mm_audit_report_text(artifact_json: &str) -> std::result::Result<String, String> {
    let artifact: CompiledModelArtifact = serde_json::from_str(artifact_json)
        .map_err(|e| format!("mm_schema_error: failed to deserialize artifact: {}", e))?;
    Ok(artifact.audit_report().to_text())
}

/// Serialize the structured `ModelAuditReport` for an `audit_design()` artifact.
///
/// This is the structured counterpart to `mm_audit_report_text`: Rust still
/// authors the report and random-term-card wording, while R receives the
/// already-authored fields for later formatters such as `explain_model()`.
///
/// @noRd
#[extendr]
fn mm_audit_report_json(artifact_json: &str) -> std::result::Result<String, String> {
    let artifact: CompiledModelArtifact = serde_json::from_str(artifact_json)
        .map_err(|e| format!("mm_schema_error: failed to deserialize artifact: {}", e))?;
    serde_json::to_string(&artifact.audit_report())
        .map_err(|e| format!("mm_schema_error: failed to serialize audit report: {}", e))
}

/// Demo of the interrupt bridge — a no-op loop that yields to R between
/// iterations so Ctrl-C can terminate it. Returns `iters` on clean
/// completion.
///
/// Used in Phase 0 only as a smoke test that the interrupt FFI binding is
/// linked correctly. Phase 1+ replaces this with the hook called from
/// inside fit / inference loops.
///
/// @noRd
#[extendr]
fn mm_interrupt_demo(iters: i32) -> i32 {
    let n = iters.max(0);
    for _ in 0..n {
        check_user_interrupt();
    }
    n
}

fn make_bootstrap_rng(seed: Option<u64>) -> (StdRng, BootstrapSeedRecord) {
    match seed {
        Some(s) => (StdRng::seed_from_u64(s), BootstrapSeedRecord::std_rng(s)),
        None => (StdRng::from_entropy(), BootstrapSeedRecord::unspecified()),
    }
}

// Linear-interpolation quantile on a pre-sorted slice. Mirrors the upstream
// `quantile_sorted` helper (which is not pub). p must be in [0, 1].
fn quantile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let pos = p.clamp(0.0, 1.0) * (n as f64 - 1.0);
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let frac = pos - pos.floor();
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

extendr_module! {
    mod mixeff;
    fn mm_parse_formula;
    fn mm_formula_manifest;
    fn mm_json_negotiate_one;
    fn mm_json_known_schemas;
    fn mm_compile_model_json;
    fn mm_fit_lmm_json;
    fn mm_fixed_effect_contrast_json;
    fn mm_fixed_effect_bootstrap_contrast_json;
    fn mm_full_model_bootstrap_contrast_json;
    fn mm_fixed_effect_bootstrap_term_json;
    fn mm_fixed_effect_term_json;
    fn mm_bootstrap_lrt_json;
    fn mm_audit_report_text;
    fn mm_audit_report_json;
    fn mm_interrupt_demo;
}
