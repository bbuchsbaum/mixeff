use std::sync::{Arc, atomic::{AtomicBool, Ordering}};
use mixeff_rs::error::MixedModelError;
use mixeff_rs::formula::parse_formula;
use mixeff_rs::model::{DataFrame, FitOptions, FitProgressCallback, LinearMixedModel,
    GeneralizedLinearMixedModel, GlmmFitOptions, Family, LinkFunction,
    try_parametricbootstrap};
use mixeff_rs::stats::bootstrap::parametricbootstrap_glmm;
use rand::{SeedableRng, rngs::StdRng};

fn callback(armed: Arc<AtomicBool>) -> FitProgressCallback {
    FitProgressCallback::new(move |_| {
        if armed.load(Ordering::SeqCst) {
            Err(MixedModelError::Interrupted("host interrupt".into()))
        } else { Ok(()) }
    })
}
fn main() {
    let mut data = DataFrame::new();
    let x: Vec<f64> = (0..48).map(|i| (i % 8) as f64 - 3.5).collect();
    let y: Vec<f64> = (0..48).map(|i| 2.0 + 0.4*x[i] + (i/8) as f64*0.3 + ((i*7)%11) as f64*0.13).collect();
    data.add_numeric("x", x).unwrap();
    data.add_numeric("y", y).unwrap();
    data.add_categorical("g", (0..48).map(|i| format!("g{}",i/8)).collect()).unwrap();
    let formula = parse_formula("y ~ 1 + x + (1 | g)").unwrap();
    let armed = Arc::new(AtomicBool::new(false));
    let mut live = LinearMixedModel::new(formula.clone(), &data, None).unwrap();
    live.fit_with_options(FitOptions::ml().with_progress_callback(callback(armed.clone()))).unwrap();
    let restored = LinearMixedModel::restore_json(&live.snapshot_json().unwrap()).unwrap();
    armed.store(true, Ordering::SeqCst);
    let live_result = try_parametricbootstrap(&mut StdRng::seed_from_u64(42), 2, &live);
    let restored_result = try_parametricbootstrap(&mut StdRng::seed_from_u64(42), 2, &restored);
    assert!(matches!(live_result, Err(MixedModelError::Interrupted(_))));
    assert!(restored_result.is_ok());
    println!("LMM: live=Interrupted; restored=Ok({} replicates)", restored_result.unwrap().fits.len());

    let mut counts = DataFrame::new();
    counts.add_numeric("x", (0..48).map(|i| (i%8) as f64 - 3.5).collect()).unwrap();
    counts.add_numeric("y", (0..48).map(|i| (1 + (i*7+i/8)%5) as f64).collect()).unwrap();
    counts.add_categorical("g", (0..48).map(|i| format!("g{}",i/8)).collect()).unwrap();
    armed.store(false, Ordering::SeqCst);
    let mut live = GeneralizedLinearMixedModel::new(formula, &counts, Family::Poisson, Some(LinkFunction::Log)).unwrap();
    live.fit_with_glmm_options(GlmmFitOptions::fast_laplace().with_progress_callback(callback(armed.clone()))).unwrap();
    let restored = GeneralizedLinearMixedModel::restore_json(&live.snapshot_json().unwrap()).unwrap();
    armed.store(true, Ordering::SeqCst);
    let live_result = parametricbootstrap_glmm(&mut StdRng::seed_from_u64(42), 2, &live);
    let restored_result = parametricbootstrap_glmm(&mut StdRng::seed_from_u64(42), 2, &restored);
    assert!(matches!(live_result, Err(MixedModelError::Interrupted(_))));
    assert!(restored_result.is_ok());
    println!("GLMM: live=Interrupted; restored=Ok({} replicates)", restored_result.unwrap().fits.len());
}
