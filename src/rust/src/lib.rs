use extendr_api::prelude::*;
use mixeff_rs::compiler::{
    compile_formula_ir, CompiledModelArtifact, ContrastMatrix, ContrastRhs, FixedEffectHypothesis,
    FixedEffectInferenceRowKind, FixedEffectInferenceTable, FixedEffectTermTestType,
    FixedEffectTestMethod, FIXED_EFFECT_INFERENCE_TABLE_SCHEMA,
    FIXED_EFFECT_INFERENCE_TABLE_SCHEMA_VERSION,
};
use mixeff_rs::formula::parse_formula;
use mixeff_rs::model::linear::ConvergenceVerificationOptions;
use mixeff_rs::model::{
    parametricbootstrap, BootstrapFailedRefitPolicy, BootstrapRefitOptions, BootstrapReplicate,
    BootstrapSeedRecord, BootstrapTarget, Family, FitOptions, FitToleranceOverrides,
    FixedEffectBootstrapOptions, GeneralizedLinearMixedModel, GlmmFitOptions, GlmmPredictionScale,
    LinearMixedModel, LinkFunction, MixedModelBootstrap, MixedModelFit, NewReLevels,
    OptimizerControl,
};
use mixeff_rs::stats::{
    profile_confint_payload, BoundaryLikelihoodRatioTest, FitSummaryPayload, LinearModelFit,
    ModelComparisonMethod, ModelComparisonOptions, ModelComparisonRefitPolicy,
    ModelComparisonTable, BOUNDARY_LRT_SCHEMA, BOUNDARY_LRT_SCHEMA_VERSION, FIT_SUMMARY_SCHEMA,
    FIT_SUMMARY_SCHEMA_VERSION, PROFILE_LIKELIHOOD_CI_SCHEMA, PROFILE_LIKELIHOOD_CI_SCHEMA_VERSION,
};
use mixeff_rs::types::Optimizer;
use nalgebra::{DMatrix, DVector};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};
use std::collections::HashMap;

// ---- caller optimizer controls (mm_control optimizer/tolerance/start) -------
//
// mixeff-rs 368a3fa added a narrow, audit-recorded escape hatch over the
// "driver chooses the optimizer" default (OptimizerControl on FitOptions /
// GlmmFitOptions). These helpers translate the optional mm_control() fields
// from the control JSON into that engine surface. Absent fields keep the
// driver's automatic selection; the chosen optimizer/source is recorded in the
// optimizer certificate by the engine.

/// Map a caller optimizer name (mm_control(optimizer=)) to the engine enum.
fn parse_optimizer_name(name: &str) -> std::result::Result<Optimizer, String> {
    match name.trim().to_ascii_lowercase().as_str() {
        "bobyqa" | "nlopt_bobyqa" => Ok(Optimizer::NloptBobyqa),
        "newuoa" | "nlopt_newuoa" => Ok(Optimizer::NloptNewuoa),
        "cobyla" => Ok(Optimizer::Cobyla),
        "pattern_search" => Ok(Optimizer::PatternSearch),
        "trust_bq" => Ok(Optimizer::TrustBq),
        "prima_bobyqa" => Ok(Optimizer::PrimaBobyqa),
        "prima_cobyla" => Ok(Optimizer::PrimaCobyla),
        "prima_lincoa" => Ok(Optimizer::PrimaLincoa),
        "prima_newuoa" => Ok(Optimizer::PrimaNewuoa),
        other => Err(format!(
            "mm_arg_error: unknown optimizer '{}'; supported: auto, bobyqa, newuoa, cobyla, \
             pattern_search, trust_bq, prima_bobyqa, prima_cobyla, prima_lincoa, prima_newuoa",
            other
        )),
    }
}

/// Coerce a JSON value to a Vec<f64>, rejecting non-numeric entries. Accepts a
/// bare number too, since jsonlite(auto_unbox = TRUE) collapses a length-1 R
/// vector (e.g. a single-theta warm start) to a JSON scalar.
fn control_f64_vec(value: &Value) -> std::result::Result<Vec<f64>, String> {
    if let Some(scalar) = value.as_f64() {
        return Ok(vec![scalar]);
    }
    value
        .as_array()
        .ok_or_else(|| "mm_arg_error: expected a numeric array in control".to_string())?
        .iter()
        .map(|v| {
            v.as_f64()
                .ok_or_else(|| "mm_arg_error: control vector entries must be numeric".to_string())
        })
        .collect()
}

/// Build an engine `OptimizerControl` from the mm_control() JSON. Every field is
/// optional; when none are set this is `OptimizerControl::default()` (the
/// driver's automatic behaviour, unchanged).
fn parse_optimizer_control(control: &Value) -> std::result::Result<OptimizerControl, String> {
    let present = |key: &str| control.get(key).filter(|v| !v.is_null());
    let mut oc = OptimizerControl::default();
    if let Some(opt) = present("optimizer").and_then(Value::as_str) {
        if !opt.eq_ignore_ascii_case("auto") {
            oc = oc.with_optimizer(parse_optimizer_name(opt)?);
        }
    }
    let mut tol = FitToleranceOverrides::default();
    let mut tol_set = false;
    if let Some(v) = present("ftol_rel").and_then(Value::as_f64) {
        tol = tol.with_ftol_rel(v);
        tol_set = true;
    }
    if let Some(v) = present("ftol_abs").and_then(Value::as_f64) {
        tol = tol.with_ftol_abs(v);
        tol_set = true;
    }
    if let Some(v) = present("xtol_rel").and_then(Value::as_f64) {
        tol = tol.with_xtol_rel(v);
        tol_set = true;
    }
    if let Some(v) = present("xtol_abs") {
        tol = tol.with_xtol_abs(control_f64_vec(v)?);
        tol_set = true;
    }
    if let Some(v) = present("initial_step") {
        tol = tol.with_initial_step(control_f64_vec(v)?);
        tol_set = true;
    }
    if tol_set {
        oc = oc.with_tolerances(tol);
    }
    if let Some(v) = present("start") {
        oc = oc.with_start_theta(control_f64_vec(v)?);
    }
    if let Some(mf) = present("max_feval").and_then(Value::as_i64) {
        if mf > 0 {
            oc = oc.with_max_feval(mf as usize);
        }
    }
    Ok(oc)
}

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
const SCHEMA_NAME_MODEL_COMPARISON_TABLE: &str = "mixedmodels.model_comparison_table";
const SCHEMA_VERSION_MODEL_COMPARISON_TABLE: &str = "1.0.0";
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
    (
        SCHEMA_NAME_MODEL_COMPARISON_TABLE,
        SCHEMA_VERSION_MODEL_COMPARISON_TABLE,
    ),
    (BOUNDARY_LRT_SCHEMA, BOUNDARY_LRT_SCHEMA_VERSION),
    (FIT_SUMMARY_SCHEMA, FIT_SUMMARY_SCHEMA_VERSION),
    (
        PROFILE_LIKELIHOOD_CI_SCHEMA,
        PROFILE_LIKELIHOOD_CI_SCHEMA_VERSION,
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
        (
            SCHEMA_NAME_MODEL_COMPARISON_TABLE,
            Robj::from(SCHEMA_VERSION_MODEL_COMPARISON_TABLE),
        ),
        (FIT_SUMMARY_SCHEMA, Robj::from(FIT_SUMMARY_SCHEMA_VERSION)),
        (
            PROFILE_LIKELIHOOD_CI_SCHEMA,
            Robj::from(PROFILE_LIKELIHOOD_CI_SCHEMA_VERSION),
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
        fit_glmm = TRUE,
        simulate = TRUE,
        inference = TRUE,
        fixed_effect_inference_table = TRUE,
        satterthwaite = TRUE,
        kenward_roger_explicit = TRUE,
        bootstrap_fixed_effect_payload = TRUE,
        model_comparison_table = TRUE,
        fit_summary_payload = TRUE,
        marginal_quantity_table = TRUE,
        marginal_quantities = TRUE,
        verify_convergence = TRUE,
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
    categorical_ordered: Strings,
) -> std::result::Result<String, String> {
    let parsed = parse_formula(formula).map_err(|e| format!("mm_formula_error: {}", e))?;
    let semantic = compile_formula_ir(&parsed);
    let df = data::build_dataframe(
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
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
    categorical_ordered: Strings,
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
        &categorical_ordered,
        &column_order,
    )?;
    let weights = optional_case_weights(&weights, df.nrow())?;

    let mut model = LinearMixedModel::new(parsed, &df, weights.as_deref())
        .map_err(|e| format!("mm_fit_error: failed to construct LMM: {}", e))?;
    // Caller optimizer controls (mm_control optimizer/tolerance/start/max_feval);
    // default keeps the driver's automatic selection. Applied here and in the
    // recompute helper so refits (predict/inference) reproduce the same fit.
    let optimizer_control = parse_optimizer_control(&_control)?;
    let fit_options = if reml {
        FitOptions::reml()
    } else {
        FitOptions::ml()
    }
    .with_optimizer_control(optimizer_control);
    model
        .fit_with_options(fit_options)
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
    let fit_summary =
        serde_json::to_value(FitSummaryPayload::from_linear_model(&model)).map_err(|e| {
            format!(
                "mm_schema_error: failed to serialize fit-summary payload: {}",
                e
            )
        })?;
    let payload = json!({
        "schema": {
            "schema_name": "mixeff.lmm_fit_result",
            "schema_version": 1
        },
        "artifact_json": artifact_json,
        "formula": model.formula().to_string(),
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
        "fit_summary": fit_summary,
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

/// Fit a generalized linear mixed-effects model and return the fit payload JSON.
///
/// This is the Stage B.1 GLMM fit primitive. The R layer owns formula/data
/// validation and user-facing S3 shape; Rust owns construction of the upstream
/// `GeneralizedLinearMixedModel`, the PIRLS fit, and the compiler artifact.
/// `method = "pirls_profiled"` maps to upstream `fast = true`. `method =
/// "joint_laplace"` maps to the labelled upstream `fast = false, n_agq = 1`
/// joint route.
///
/// @noRd
#[extendr]
fn mm_fit_glmm_json(
    formula: &str,
    family: &str,
    link: &str,
    method: &str,
    n_agq: i32,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    categorical_ordered: Strings,
    weights: Doubles,
    offset: Doubles,
    control_json: &str,
) -> std::result::Result<String, String> {
    let _control: Value = serde_json::from_str(control_json)
        .map_err(|e| format!("mm_fit_error: invalid control JSON: {}", e))?;
    let parsed = parse_formula(formula).map_err(|e| format!("mm_formula_error: {}", e))?;
    let df = data::build_dataframe(
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
        &column_order,
    )?;
    let family = glmm_family(family)?;
    let link = glmm_link(link)?;
    let (fast, method_label) = glmm_method(method, n_agq)?;
    let n_agq = usize::try_from(n_agq)
        .map_err(|_| "mm_arg_error: nAGQ must be a positive integer".to_string())?;

    // Prior weights (e.g. binomial trial counts) and a fixed linear-predictor
    // offset are both supported by the engine; route to the matching
    // constructor when supplied (empty Doubles => None).
    let weights = optional_case_weights(&weights, df.nrow())?;
    let offset = optional_offset(&offset, df.nrow())?;
    let mut model = match (weights, offset) {
        (None, None) => GeneralizedLinearMixedModel::new(parsed, &df, family, Some(link)),
        (Some(w), None) => {
            GeneralizedLinearMixedModel::new_with_weights(parsed, &df, family, Some(link), w)
        }
        (None, Some(o)) => {
            GeneralizedLinearMixedModel::new_with_offset(parsed, &df, family, Some(link), o)
        }
        (Some(w), Some(o)) => GeneralizedLinearMixedModel::new_with_weights_and_offset(
            parsed,
            &df,
            family,
            Some(link),
            w,
            o,
        ),
    }
    .map_err(|e| format!("mm_fit_error: failed to construct GLMM: {}", e))?;

    // Caller optimizer controls (mm_control optimizer/tolerance/start/max_feval)
    // route through the GLMM optimizer-control surface; absent fields keep the
    // driver's automatic selection. max_feval still caps the joint path's
    // otherwise engine-chosen budget, now via OptimizerControl rather than a
    // direct optsum poke.
    let optimizer_control = parse_optimizer_control(&_control)?;
    let glmm_options = if fast {
        GlmmFitOptions::fast_laplace()
    } else {
        GlmmFitOptions::joint_laplace()
    }
    .with_n_agq(n_agq)
    .with_verbose(false)
    .with_optimizer_control(optimizer_control);
    model
        .fit_with_glmm_options(glmm_options)
        .map_err(|e| format!("mm_fit_error: failed to fit GLMM: {}", e))?;

    let artifact_json = serde_json::to_string(model.compiler_artifact()).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize GLMM post-fit artifact: {}",
            e
        )
    })?;
    let artifact_value: Value = serde_json::from_str(&artifact_json).map_err(|e| {
        format!(
            "mm_schema_error: failed to inspect GLMM post-fit artifact JSON: {}",
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
    let log_likelihood = model.loglikelihood();
    let dof = model.dof();
    let nobs = model.nobs();
    let df_residual = nobs.saturating_sub(dof);
    let opt = model.opt_summary();
    let fit_summary = serde_json::to_value(FitSummaryPayload::from_generalized_model(&model))
        .map_err(|e| {
            format!(
                "mm_schema_error: failed to serialize GLMM fit-summary payload: {}",
                e
            )
        })?;

    let payload = json!({
        "schema": {
            "schema_name": "mixeff.glmm_fit_result",
            "schema_version": 1
        },
        "artifact_json": artifact_json,
        "formula": model.formula_label().unwrap_or_else(|| formula.to_string()),
        "family": glmm_family_label(family),
        "link": glmm_link_label(link),
        "method": method_label,
        "n_agq": n_agq,
        "beta": beta.iter().copied().collect::<Vec<_>>(),
        "beta_names": beta_names,
        "theta": model.theta(),
        "dispersion": model.dispersion(false),
        "log_likelihood": log_likelihood,
        "deviance": -2.0 * log_likelihood,
        "aic": model.aic(),
        "bic": model.bic(),
        "nobs": nobs,
        "dof": dof,
        "df_residual": df_residual,
        "fit_status": fit_status,
        "std_errors": std_errors.iter().copied().collect::<Vec<_>>(),
        "fitted": model.fitted().iter().copied().collect::<Vec<_>>(),
        "residuals": model.residuals().iter().copied().collect::<Vec<_>>(),
        "ranef": random_effects_json_glmm(&model),
        "varcorr": varcorr_json_glmm(&model),
        "fit_summary": fit_summary,
        "optimizer": {
            "backend": opt.backend_name(),
            "algorithm": opt.optimizer_name(),
            "return_value": opt.return_value.as_str(),
            "function_evaluations": opt.feval,
            "objective": opt.fmin,
            "reml": opt.reml
        }
    });

    serde_json::to_string(&payload).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize GLMM fit result: {}",
            e
        )
    })
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
    categorical_ordered: Strings,
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
        &categorical_ordered,
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
    categorical_ordered: Strings,
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
        &categorical_ordered,
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
    categorical_ordered: Strings,
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
        &categorical_ordered,
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
    let l_row: Vec<f64> = hypothesis.l.values.row(0).iter().copied().collect();
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
            "mm_inference_unavailable: bootstrap produced no finite contrast statistics"
                .to_string(),
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
    categorical_ordered: Strings,
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
        &categorical_ordered,
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
    categorical_ordered: Strings,
    weights: Doubles,
    control_json: &str,
    method: &str,
    term_test_type: &str,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
        &weights,
        control_json,
    )?;
    let method = fixed_effect_test_method(method)?;
    let term_test_type = fixed_effect_term_test_type(term_test_type)?;

    serde_json::to_string(&model.fixed_effect_term_inference_table_for_type(method, term_test_type))
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
    categorical_ordered: Strings,
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
        &categorical_ordered,
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
        &categorical_ordered,
        &weights,
        control_json,
    )?;
    let options = fixed_effect_bootstrap_options(bootstrap_options_json)?;

    let observed_logl_red = reduced.loglikelihood();
    let observed_logl_alt = alternative.loglikelihood();
    let observed_lrt = 2.0 * (observed_logl_alt - observed_logl_red);
    if !observed_lrt.is_finite() {
        return Err("mm_inference_unavailable: observed LRT statistic is not finite".to_string());
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

/// Build an upstream model-comparison table for fitted LMM payloads.
///
/// Each list element is the R-side bridge payload for one fitted model
/// (`formula_string`, `REML`, `spec_data`, `weights`, and `control_json`).
/// Returning the upstream `ModelComparisonTable` keeps nestedness, ML-refit,
/// information-criteria, and reason-code rules owned by the Rust contract.
///
/// @noRd
#[extendr]
fn mm_compare_models_json(
    model_payloads: List,
    method: &str,
    refit_policy: &str,
) -> std::result::Result<String, String> {
    if model_payloads.len() < 2 {
        return Err(
            "mm_arg_error: model comparison requires at least two fitted models".to_string(),
        );
    }

    let method = model_comparison_method(method)?;
    let refit_policy = model_comparison_refit_policy(refit_policy)?;
    let mut models: Vec<LinearMixedModel> = Vec::with_capacity(model_payloads.len());
    for (idx, payload) in model_payloads.values().enumerate() {
        models.push(fit_lmm_from_bridge_payload_robj(&payload, idx + 1)?);
    }
    let model_refs = models
        .iter()
        .map(|model| model as &dyn MixedModelFit)
        .collect::<Vec<_>>();

    let table = ModelComparisonTable::compare_with_options(
        &model_refs,
        ModelComparisonOptions {
            method,
            refit_policy,
        },
    )
    .map_err(|e| format!("mm_inference_unavailable: model comparison failed: {}", e))?;

    let payload = json!({
        "schema": {
            "schema_name": SCHEMA_NAME_MODEL_COMPARISON_TABLE,
            "schema_version": SCHEMA_VERSION_MODEL_COMPARISON_TABLE,
        },
        "payload": table,
    });
    serde_json::to_string(&payload).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize model comparison: {}",
            e
        )
    })
}

/// Boundary-aware variance-component likelihood-ratio test.
///
/// `reduced_payload` may be NULL when the tested random term is the only
/// random-effect term. In that case the reduced model is the fixed-effect
/// linear model built from the full model's response and fixed-effect design.
/// Otherwise it is the bridge payload for the reduced LMM.
///
/// @noRd
#[extendr]
fn mm_boundary_lrt_json(
    reduced_payload: Robj,
    full_payload: Robj,
    reduced_formula: &str,
) -> std::result::Result<String, String> {
    let full = fit_lmm_from_bridge_payload_robj(&full_payload, 2)?;
    let result = if reduced_payload.is_null() {
        let reduced = LinearModelFit::fit(
            full.response().clone(),
            full.model_matrix().clone(),
            Some(reduced_formula.to_string()),
        )
        .map_err(|e| format!("mm_inference_unavailable: boundary LRT reduced LM failed: {e}"))?;
        BoundaryLikelihoodRatioTest::variance_component(&reduced, &full)
    } else {
        let reduced = fit_lmm_from_bridge_payload_robj(&reduced_payload, 1)?;
        BoundaryLikelihoodRatioTest::variance_component(&reduced, &full)
    };

    serde_json::to_string(&result)
        .map_err(|e| format!("mm_schema_error: failed to serialize boundary LRT: {e}"))
}

/// Bounded convergence verification for a fitted LMM.
///
/// Rebuilds and refits the model from the R-side bridge payload (the same
/// payload `mm_compare_models_json` / `mm_boundary_lrt_json` use), then runs
/// the engine's verification workflow: restart from the optimum, jittered
/// restarts, and (optionally) an alternate-optimizer consensus pass. The
/// verifier re-runs its own fits internally, so no optimizer state has to
/// cross the boundary. `options_json` carries R-side overrides; absent
/// fields keep the engine defaults.
///
/// @noRd
#[extendr]
fn mm_verify_convergence_json(
    fit_payload: Robj,
    options_json: &str,
) -> std::result::Result<String, String> {
    let overrides: Value = serde_json::from_str(options_json)
        .map_err(|e| format!("mm_arg_error: invalid verification options JSON: {}", e))?;

    let mut options = ConvergenceVerificationOptions::default();
    if let Some(v) = overrides
        .get("restart_from_optimum")
        .and_then(Value::as_bool)
    {
        options.restart_from_optimum = v;
    }
    if let Some(v) = overrides.get("jitter_starts").and_then(Value::as_u64) {
        options.jitter_starts = v as usize;
    }
    if let Some(v) = overrides.get("jitter_scale").and_then(Value::as_f64) {
        options.jitter_scale = v;
    }
    if let Some(v) = overrides
        .get("run_optimizer_consensus")
        .and_then(Value::as_bool)
    {
        options.run_optimizer_consensus = v;
    }
    if let Some(v) = overrides
        .get("max_function_evaluations")
        .and_then(Value::as_u64)
    {
        options.max_function_evaluations = v as usize;
    }
    if let Some(v) = overrides.get("objective_tolerance").and_then(Value::as_f64) {
        options.objective_tolerance = v;
    }
    if let Some(v) = overrides.get("theta_tolerance").and_then(Value::as_f64) {
        options.theta_tolerance = v;
    }
    if let Some(v) = overrides.get("beta_tolerance").and_then(Value::as_f64) {
        options.beta_tolerance = v;
    }

    let mut model = fit_lmm_from_bridge_payload_robj(&fit_payload, 1)?;
    let verification = model
        .verify_convergence_with_options(options)
        .map_err(|e| {
            format!(
                "mm_inference_unavailable: convergence verification failed: {}",
                e
            )
        })?;

    serde_json::to_string(&verification).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize convergence verification: {}",
            e
        )
    })
}

fn fit_lmm_from_bridge_data(
    formula: &str,
    reml: bool,
    column_order: &Strings,
    numeric_columns: &List,
    categorical_values: &List,
    categorical_levels: &List,
    categorical_ordered: &Strings,
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
        categorical_ordered,
        column_order,
    )?;
    let weights = optional_case_weights(weights, df.nrow())?;
    let mut model = LinearMixedModel::new(parsed, &df, weights.as_deref())
        .map_err(|e| format!("mm_fit_error: failed to construct LMM: {}", e))?;
    // Caller optimizer controls (mm_control optimizer/tolerance/start/max_feval);
    // default keeps the driver's automatic selection. Applied here and in the
    // recompute helper so refits (predict/inference) reproduce the same fit.
    let optimizer_control = parse_optimizer_control(&_control)?;
    let fit_options = if reml {
        FitOptions::reml()
    } else {
        FitOptions::ml()
    }
    .with_optimizer_control(optimizer_control);
    model
        .fit_with_options(fit_options)
        .map_err(|e| format!("mm_fit_error: failed to fit LMM: {}", e))?;
    Ok(model)
}

fn fit_lmm_from_bridge_payload_robj(
    payload: &Robj,
    index: usize,
) -> std::result::Result<LinearMixedModel, String> {
    let payload_list = List::try_from(payload).map_err(|e| {
        format!("mm_schema_error: comparison payload {index} must be a list: {e:?}")
    })?;
    let payload_map = list_to_map(&payload_list, &format!("comparison payload {index}"))?;
    let spec_data = required_list(&payload_map, "spec_data", index)?;
    let spec_map = list_to_map(&spec_data, &format!("comparison payload {index} spec_data"))?;

    let formula = required_string(&payload_map, "formula_string", index)?;
    let reml = required_bool(&payload_map, "REML", index)?;
    let column_order = required_strings(&spec_map, "column_order", index)?;
    let numeric_columns = required_list(&spec_map, "numeric_columns", index)?;
    let categorical_values = required_list(&spec_map, "categorical_values", index)?;
    let categorical_levels = required_list(&spec_map, "categorical_levels", index)?;
    let categorical_ordered = required_strings(&spec_map, "categorical_ordered", index)?;
    let weights = required_doubles(&payload_map, "weights", index)?;
    let control_json = required_string(&payload_map, "control_json", index)?;

    fit_lmm_from_bridge_data(
        &formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
        &weights,
        &control_json,
    )
}

fn glmm_family(family: &str) -> std::result::Result<Family, String> {
    match family {
        "bernoulli" => Ok(Family::Bernoulli),
        // Grouped binomial (proportion response + trial-count weights). The R
        // layer sends "bernoulli" for ungrouped 0/1 responses and "binomial"
        // only when trial weights are supplied.
        "binomial" => Ok(Family::Binomial),
        "poisson" => Ok(Family::Poisson),
        "gamma" => Ok(Family::Gamma),
        other => Err(format!("mm_fit_error: unsupported GLMM family `{other}`")),
    }
}

fn glmm_link(link: &str) -> std::result::Result<LinkFunction, String> {
    match link {
        "identity" => Ok(LinkFunction::Identity),
        "log" => Ok(LinkFunction::Log),
        "logit" => Ok(LinkFunction::Logit),
        "probit" => Ok(LinkFunction::Probit),
        "cloglog" => Ok(LinkFunction::Cloglog),
        "inverse" => Ok(LinkFunction::Inverse),
        "sqrt" => Ok(LinkFunction::Sqrt),
        other => Err(format!("mm_fit_error: unsupported GLMM link `{other}`")),
    }
}

fn glmm_method(method: &str, n_agq: i32) -> std::result::Result<(bool, &'static str), String> {
    if n_agq < 1 {
        return Err("mm_arg_error: nAGQ must be a positive integer".to_string());
    }
    match method {
        "pirls_profiled" => Ok((true, "pirls_profiled")),
        "joint_laplace" => {
            if n_agq > 1 {
                return Err(
                    "mm_arg_error: method='joint_laplace' currently requires nAGQ <= 1".to_string(),
                );
            }
            Ok((false, "joint_laplace"))
        }
        other => Err(format!("mm_arg_error: unsupported GLMM method `{other}`")),
    }
}

fn glmm_family_label(family: Family) -> &'static str {
    match family {
        Family::Normal => "normal",
        Family::Bernoulli => "bernoulli",
        Family::Binomial => "binomial",
        Family::Poisson => "poisson",
        Family::Gamma => "gamma",
        Family::InverseGaussian => "inverse_gaussian",
        _ => "unknown",
    }
}

fn glmm_link_label(link: LinkFunction) -> &'static str {
    match link {
        LinkFunction::Identity => "identity",
        LinkFunction::Log => "log",
        LinkFunction::Logit => "logit",
        LinkFunction::Probit => "probit",
        LinkFunction::Cloglog => "cloglog",
        LinkFunction::Inverse => "inverse",
        LinkFunction::Sqrt => "sqrt",
        _ => "unknown",
    }
}

fn list_to_map(list: &List, context: &str) -> std::result::Result<HashMap<String, Robj>, String> {
    list.try_into()
        .map_err(|e| format!("mm_schema_error: failed to read {context}: {e:?}"))
}

fn required_robj(
    map: &HashMap<String, Robj>,
    field: &str,
    index: usize,
) -> std::result::Result<Robj, String> {
    map.get(field)
        .cloned()
        .ok_or_else(|| format!("mm_schema_error: comparison payload {index} is missing `{field}`"))
}

fn required_string(
    map: &HashMap<String, Robj>,
    field: &str,
    index: usize,
) -> std::result::Result<String, String> {
    let robj = required_robj(map, field, index)?;
    String::try_from(&robj).map_err(|e| {
        format!(
            "mm_schema_error: comparison payload {index} field `{field}` must be a string: {e:?}"
        )
    })
}

fn required_bool(
    map: &HashMap<String, Robj>,
    field: &str,
    index: usize,
) -> std::result::Result<bool, String> {
    let robj = required_robj(map, field, index)?;
    bool::try_from(&robj).map_err(|e| {
        format!(
            "mm_schema_error: comparison payload {index} field `{field}` must be TRUE/FALSE: {e:?}"
        )
    })
}

fn required_list(
    map: &HashMap<String, Robj>,
    field: &str,
    index: usize,
) -> std::result::Result<List, String> {
    let robj = required_robj(map, field, index)?;
    List::try_from(&robj).map_err(|e| {
        format!("mm_schema_error: comparison payload {index} field `{field}` must be a list: {e:?}")
    })
}

fn required_strings(
    map: &HashMap<String, Robj>,
    field: &str,
    index: usize,
) -> std::result::Result<Strings, String> {
    let robj = required_robj(map, field, index)?;
    Strings::try_from(&robj).map_err(|e| {
        format!("mm_schema_error: comparison payload {index} field `{field}` must be a character vector: {e:?}")
    })
}

fn required_doubles(
    map: &HashMap<String, Robj>,
    field: &str,
    index: usize,
) -> std::result::Result<Doubles, String> {
    let robj = required_robj(map, field, index)?;
    Doubles::try_from(&robj).map_err(|e| {
        format!("mm_schema_error: comparison payload {index} field `{field}` must be a numeric vector: {e:?}")
    })
}

fn model_comparison_method(method: &str) -> std::result::Result<ModelComparisonMethod, String> {
    match method {
        "auto" => Ok(ModelComparisonMethod::Auto),
        "lrt" | "likelihood_ratio" => Ok(ModelComparisonMethod::LikelihoodRatio),
        "bootstrap" => Ok(ModelComparisonMethod::Auto),
        "aic" | "information_criteria" => Ok(ModelComparisonMethod::InformationCriteria),
        other => Err(format!(
            "mm_arg_error: unsupported model-comparison method `{other}`"
        )),
    }
}

fn model_comparison_refit_policy(
    refit_policy: &str,
) -> std::result::Result<ModelComparisonRefitPolicy, String> {
    match refit_policy {
        "error" => Ok(ModelComparisonRefitPolicy::Error),
        "auto" | "ml" => Ok(ModelComparisonRefitPolicy::Ml),
        "never" => Ok(ModelComparisonRefitPolicy::Never),
        other => Err(format!(
            "mm_arg_error: unsupported model-comparison refit policy `{other}`"
        )),
    }
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

fn optional_offset(
    offset: &Doubles,
    n_obs: usize,
) -> std::result::Result<Option<Vec<f64>>, String> {
    let values = offset.iter().map(|v| v.0).collect::<Vec<_>>();
    if values.is_empty() {
        return Ok(None);
    }
    if values.len() != n_obs {
        return Err(format!(
            "mm_data_error: offset length ({}) does not match number of observations ({n_obs})",
            values.len()
        ));
    }
    for (i, o) in values.iter().enumerate() {
        if !o.is_finite() {
            return Err(format!(
                "mm_data_error: offset at index {i} must be finite (got {o})"
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

fn fixed_effect_term_test_type(
    term_test_type: &str,
) -> std::result::Result<FixedEffectTermTestType, String> {
    match term_test_type {
        "I" | "1" | "type_i" | "type_1" => Ok(FixedEffectTermTestType::TypeI),
        "II" | "2" | "type_ii" | "type_2" => Ok(FixedEffectTermTestType::TypeII),
        "III" | "3" | "type_iii" | "type_3" => Ok(FixedEffectTermTestType::TypeIII),
        other => Err(format!(
            "mm_inference_unavailable: unsupported fixed-effect term test type `{other}`"
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
        .reterms()
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

fn random_effects_json_glmm(model: &GeneralizedLinearMixedModel) -> Value {
    let effects = model.ranef();
    let terms = model
        .lmm()
        .reterms()
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
    varcorr_value(&vc)
}

fn varcorr_json_glmm(model: &GeneralizedLinearMixedModel) -> Value {
    let vc = model.varcorr();
    varcorr_value(&vc)
}

fn varcorr_value(vc: &mixeff_rs::stats::VarCorr) -> Value {
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

/// Render the compact `audit_design()` summary as text.
///
/// The compact counterpart to `mm_audit_report_text`: it returns the
/// upstream `ModelAuditReport::render_summary` rendering (Audit Summary
/// plus the Requested Model section). R's default `print.mm_audit()`
/// emits this, so the compact view stays fully upstream-authored — the
/// R9 "no advice creep" contract holds without R slicing rendered text.
///
/// @noRd
#[extendr]
fn mm_audit_report_summary_text(artifact_json: &str) -> std::result::Result<String, String> {
    let artifact: CompiledModelArtifact = serde_json::from_str(artifact_json)
        .map_err(|e| format!("mm_schema_error: failed to deserialize artifact: {}", e))?;
    Ok(artifact.audit_report().render_summary())
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

/// Conditional variance matrices of the random effects, serialized for R.
///
/// Mirrors `lme4::condVar` / `ranef(condVar = TRUE)`. Each RE term yields a
/// `p × p × n` PSD array (one `p × p` block per grouping level), flattened
/// column-major so `array(payload$postvar, dim = payload$dim)` reconstructs the
/// 3-D array R expects. `names` indexes the leading two dimensions (slope
/// names per term) and `levels` indexes the trailing dimension.
///
/// @noRd
#[extendr]
fn mm_lmm_cond_var_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    categorical_ordered: Strings,
    weights: Doubles,
    control_json: &str,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
        &weights,
        control_json,
    )?;
    let condvar = model.cond_var();
    let terms: Vec<Value> = model
        .reterms()
        .iter()
        .zip(condvar.iter())
        .map(|(rt, blocks)| {
            let p = rt.vsize;
            let n = blocks.len();
            // Column-major flatten: R's array(x, dim = c(p,p,n)) reads
            // the slice row-by-row within each column, columns within each
            // slice. That maps to mat[(row, col)] in DMatrix iteration.
            let mut postvar: Vec<f64> = Vec::with_capacity(p * p * n);
            for level_idx in 0..n {
                let mat = &blocks[level_idx];
                for col in 0..p {
                    for row in 0..p {
                        postvar.push(mat[(row, col)]);
                    }
                }
            }
            json!({
                "group": rt.grouping_name.as_str(),
                "names": rt.cnames.clone(),
                "levels": rt.levels.clone(),
                "postvar": postvar,
                "dim": [p, p, n],
            })
        })
        .collect();
    let payload = json!({
        "schema": {
            "schema_name": "mixeff.lmm_cond_var",
            "schema_version": 1,
        },
        "terms": terms,
    });
    serde_json::to_string(&payload).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize cond_var payload: {}",
            e
        )
    })
}

/// New-data predictions through `LinearMixedModel::predict_new`.
///
/// `allow_new_levels_policy` is the string-form of `NewReLevels`:
///   "error"       -> NewReLevels::Error
///   "population"  -> NewReLevels::Population
///   "missing"     -> NewReLevels::Missing
///
/// Returns a JSON payload carrying one prediction per `newdata` row.
/// Predictions whose upstream value is `None` (unseen grouping level under
/// the "missing" policy) are encoded as JSON `null`, which the R side
/// translates to `NA_real_`.
///
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn mm_lmm_predict_new_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    categorical_ordered: Strings,
    weights: Doubles,
    control_json: &str,
    new_column_order: Strings,
    new_numeric_columns: List,
    new_categorical_values: List,
    new_categorical_levels: List,
    new_categorical_ordered: Strings,
    allow_new_levels_policy: &str,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
        &weights,
        control_json,
    )?;
    let newdf = data::build_dataframe(
        &new_numeric_columns,
        &new_categorical_values,
        &new_categorical_levels,
        &new_categorical_ordered,
        &new_column_order,
    )?;
    let policy = match allow_new_levels_policy {
        "error" => NewReLevels::Error,
        "population" => NewReLevels::Population,
        "missing" => NewReLevels::Missing,
        other => {
            return Err(format!(
                "mm_arg_error: unsupported allow_new_levels policy `{}`; expected one of error|population|missing",
                other
            ));
        }
    };
    let predictions = model
        .predict_new(&newdf, policy)
        .map_err(|e| format!("mm_inference_unavailable: predict_new failed: {}", e))?;
    let pred_array: Vec<Value> = predictions
        .iter()
        .map(|o| match o {
            Some(v) => json!(*v),
            None => Value::Null,
        })
        .collect();
    let payload = json!({
        "schema": {
            "schema_name": "mixeff.lmm_predict_new",
            "schema_version": 1,
        },
        "predictions": pred_array,
        "policy": allow_new_levels_policy,
        "n_new": predictions.len(),
    });
    serde_json::to_string(&payload).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize predict_new payload: {}",
            e
        )
    })
}

/// New-data prediction VARIANCE / intervals through
/// `LinearMixedModel::predict_new_variance_with_level`.
///
/// Returns the engine `PredictionVariancePayload` JSON (serde): one row per
/// `newdata` row carrying `se_fit`, the fixed / random / covariance variance
/// components, `confidence_*` / `prediction_*` bounds at `level`, and a
/// row-level `status` / `reason` (unavailable for new grouping levels). The R
/// side maps unavailable rows to `NA` with the reason, preserving the
/// no-fake-certainty contract.
///
/// `allow_new_levels_policy` is the string form of `NewReLevels`
/// ("error" | "population" | "missing").
///
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn mm_lmm_predict_new_variance_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    categorical_ordered: Strings,
    weights: Doubles,
    control_json: &str,
    new_column_order: Strings,
    new_numeric_columns: List,
    new_categorical_values: List,
    new_categorical_levels: List,
    new_categorical_ordered: Strings,
    allow_new_levels_policy: &str,
    level: f64,
) -> std::result::Result<String, String> {
    let model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
        &weights,
        control_json,
    )?;
    let newdf = data::build_dataframe(
        &new_numeric_columns,
        &new_categorical_values,
        &new_categorical_levels,
        &new_categorical_ordered,
        &new_column_order,
    )?;
    let policy = match allow_new_levels_policy {
        "error" => NewReLevels::Error,
        "population" => NewReLevels::Population,
        "missing" => NewReLevels::Missing,
        other => {
            return Err(format!(
                "mm_arg_error: unsupported allow_new_levels policy `{}`; expected one of error|population|missing",
                other
            ));
        }
    };
    let payload = model
        .predict_new_variance_with_level(&newdf, policy, level)
        .map_err(|e| {
            format!(
                "mm_inference_unavailable: predict_new_variance failed: {}",
                e
            )
        })?;
    serde_json::to_string(&payload).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize predict_new_variance payload: {}",
            e
        )
    })
}

/// Reconstruct and fit a `GeneralizedLinearMixedModel` from R-side bridge data.
/// Mirrors the construction/fit in `mm_fit_glmm_json` so post-fit queries
/// (e.g. prediction variance) reproduce the same fit. Kept separate from the
/// fit-result bridge so the latter's serialization path is left untouched.
#[allow(clippy::too_many_arguments)]
fn fit_glmm_from_bridge_data(
    formula: &str,
    family: &str,
    link: &str,
    method: &str,
    n_agq: i32,
    column_order: &Strings,
    numeric_columns: &List,
    categorical_values: &List,
    categorical_levels: &List,
    categorical_ordered: &Strings,
    weights: &Doubles,
    offset: &Doubles,
    control_json: &str,
) -> std::result::Result<GeneralizedLinearMixedModel, String> {
    let _control: Value = serde_json::from_str(control_json)
        .map_err(|e| format!("mm_fit_error: invalid control JSON: {}", e))?;
    let parsed = parse_formula(formula).map_err(|e| format!("mm_formula_error: {}", e))?;
    let df = data::build_dataframe(
        numeric_columns,
        categorical_values,
        categorical_levels,
        categorical_ordered,
        column_order,
    )?;
    let family = glmm_family(family)?;
    let link = glmm_link(link)?;
    let (fast, _method_label) = glmm_method(method, n_agq)?;
    let n_agq = usize::try_from(n_agq)
        .map_err(|_| "mm_arg_error: nAGQ must be a positive integer".to_string())?;
    let weights = optional_case_weights(weights, df.nrow())?;
    let offset = optional_offset(offset, df.nrow())?;
    let mut model = match (weights, offset) {
        (None, None) => GeneralizedLinearMixedModel::new(parsed, &df, family, Some(link)),
        (Some(w), None) => {
            GeneralizedLinearMixedModel::new_with_weights(parsed, &df, family, Some(link), w)
        }
        (None, Some(o)) => {
            GeneralizedLinearMixedModel::new_with_offset(parsed, &df, family, Some(link), o)
        }
        (Some(w), Some(o)) => GeneralizedLinearMixedModel::new_with_weights_and_offset(
            parsed,
            &df,
            family,
            Some(link),
            w,
            o,
        ),
    }
    .map_err(|e| format!("mm_fit_error: failed to construct GLMM: {}", e))?;
    let optimizer_control = parse_optimizer_control(&_control)?;
    let glmm_options = if fast {
        GlmmFitOptions::fast_laplace()
    } else {
        GlmmFitOptions::joint_laplace()
    }
    .with_n_agq(n_agq)
    .with_verbose(false)
    .with_optimizer_control(optimizer_control);
    model
        .fit_with_glmm_options(glmm_options)
        .map_err(|e| format!("mm_fit_error: failed to fit GLMM: {}", e))?;
    Ok(model)
}

/// New-data prediction VARIANCE / intervals for a GLMM through
/// `GeneralizedLinearMixedModel::predict_new_variance_with_level`.
///
/// `scale` is the string form of `GlmmPredictionScale` ("link" | "response");
/// `allow_new_levels_policy` is the string form of `NewReLevels`. Returns the
/// engine `PredictionVariancePayload` JSON. Certified (joint-Laplace) fits
/// return available conditional rows; fast-PIRLS / new-level rows are marked
/// degraded / unavailable with a row-level reason.
///
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn mm_glmm_predict_new_variance_json(
    formula: &str,
    family: &str,
    link: &str,
    method: &str,
    n_agq: i32,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    categorical_ordered: Strings,
    weights: Doubles,
    offset: Doubles,
    control_json: &str,
    new_column_order: Strings,
    new_numeric_columns: List,
    new_categorical_values: List,
    new_categorical_levels: List,
    new_categorical_ordered: Strings,
    scale: &str,
    allow_new_levels_policy: &str,
    level: f64,
) -> std::result::Result<String, String> {
    let model = fit_glmm_from_bridge_data(
        formula,
        family,
        link,
        method,
        n_agq,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
        &weights,
        &offset,
        control_json,
    )?;
    let newdf = data::build_dataframe(
        &new_numeric_columns,
        &new_categorical_values,
        &new_categorical_levels,
        &new_categorical_ordered,
        &new_column_order,
    )?;
    let scale = match scale {
        "link" => GlmmPredictionScale::Link,
        "response" => GlmmPredictionScale::Response,
        other => {
            return Err(format!(
                "mm_arg_error: unsupported prediction scale `{}`; expected one of link|response",
                other
            ));
        }
    };
    let policy = match allow_new_levels_policy {
        "error" => NewReLevels::Error,
        "population" => NewReLevels::Population,
        "missing" => NewReLevels::Missing,
        other => {
            return Err(format!(
                "mm_arg_error: unsupported allow_new_levels policy `{}`; expected one of error|population|missing",
                other
            ));
        }
    };
    let payload = model
        .predict_new_variance_with_level(&newdf, scale, policy, level)
        .map_err(|e| {
            format!(
                "mm_inference_unavailable: GLMM predict_new_variance failed: {}",
                e
            )
        })?;
    serde_json::to_string(&payload).map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize GLMM predict_new_variance payload: {}",
            e
        )
    })
}

/// Profile-likelihood confidence intervals through
/// `mixeff_rs::stats::profile_confint_payload`.
///
/// `level` is the confidence level in `(0, 1)`. The upstream payload carries
/// both the computed intervals and the raw profile rows. Under REML the
/// upstream contract omits beta from the profiled parameters; the R-side
/// wrapper turns this absence into a typed `profile_beta_unavailable_under_reml`
/// reason rather than fabricating beta CIs.
///
/// @noRd
#[extendr]
fn mm_lmm_profile_confint_json(
    formula: &str,
    reml: bool,
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    categorical_ordered: Strings,
    weights: Doubles,
    control_json: &str,
    level: f64,
) -> std::result::Result<String, String> {
    if !(level > 0.0 && level < 1.0) {
        return Err(format!(
            "mm_arg_error: profile confint level must be in (0, 1); got {}",
            level
        ));
    }
    let mut model = fit_lmm_from_bridge_data(
        formula,
        reml,
        &column_order,
        &numeric_columns,
        &categorical_values,
        &categorical_levels,
        &categorical_ordered,
        &weights,
        control_json,
    )?;
    let payload = profile_confint_payload(&mut model, level).map_err(|e| {
        format!(
            "mm_inference_unavailable: profile_confint_payload failed: {}",
            e
        )
    })?;
    payload.to_json().map_err(|e| {
        format!(
            "mm_schema_error: failed to serialize profile CI payload: {}",
            e
        )
    })
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
    fn mm_fit_glmm_json;
    fn mm_fixed_effect_contrast_json;
    fn mm_fixed_effect_bootstrap_contrast_json;
    fn mm_full_model_bootstrap_contrast_json;
    fn mm_fixed_effect_bootstrap_term_json;
    fn mm_fixed_effect_term_json;
    fn mm_bootstrap_lrt_json;
    fn mm_compare_models_json;
    fn mm_boundary_lrt_json;
    fn mm_verify_convergence_json;
    fn mm_audit_report_text;
    fn mm_audit_report_summary_text;
    fn mm_audit_report_json;
    fn mm_lmm_cond_var_json;
    fn mm_lmm_predict_new_json;
    fn mm_lmm_predict_new_variance_json;
    fn mm_glmm_predict_new_variance_json;
    fn mm_lmm_profile_confint_json;
    fn mm_interrupt_demo;
}
