//! Translate R column specs into the upstream `mixedmodels::DataFrame`.
//!
//! The R wire format is three parallel lists keyed by column name:
//!   - `numeric_columns`: list(name = REAL vector)
//!   - `categorical_values`: list(name = character vector of observed values)
//!   - `categorical_levels`: list(name = character vector of canonical levels)
//! plus a `column_order` character vector that names every column in the
//! original `data.frame` order and a `categorical_ordered` character vector
//! naming the subset of categorical columns that are ordered factors (coded
//! with `contr.poly` polynomial contrasts rather than treatment coding). Each
//! column appears in exactly one of the numeric/categorical lists.
//!
//! Rust trusts the R wrapper to have validated:
//!   - every value in `categorical_values[name]` is in `categorical_levels[name]`,
//!   - every name in `column_order` appears in exactly one of the column lists,
//!   - all columns have the same length (data.frame invariant).
//! Violations of these invariants produce typed `mm_data_error` strings the
//! R caller routes back to a typed condition.
//!
//! `DataFrame::add_numeric` and `add_categorical_with_levels` validate inputs
//! and return typed `MixedModelError`s (non-finite numeric, level mismatch,
//! column-length mismatch). The R wrapper's `mm_translate_data()` should
//! satisfy these preconditions, so Err is unexpected — but we propagate it
//! through `mm_data_error:` rather than ignore it, per the no-silent-surgery
//! contract.

use std::collections::{HashMap, HashSet};

use extendr_api::prelude::*;

use mixeff_rs::model::{CategoricalContrast, DataFrame};

pub(crate) fn build_dataframe(
    numeric_columns: &List,
    categorical_values: &List,
    categorical_levels: &List,
    categorical_ordered: &Strings,
    column_order: &Strings,
) -> std::result::Result<DataFrame, String> {
    let mut numeric_map: HashMap<String, Vec<f64>> = HashMap::with_capacity(numeric_columns.len());
    for (name, robj) in numeric_columns.iter() {
        let values = robj.as_real_slice().ok_or_else(|| {
            format!(
                "mm_data_error: numeric column '{name}' is not a REAL vector ({:?})",
                robj.rtype()
            )
        })?;
        numeric_map.insert(name.to_string(), values.to_vec());
    }

    let mut cat_values_map: HashMap<String, Vec<String>> =
        HashMap::with_capacity(categorical_values.len());
    for (name, robj) in categorical_values.iter() {
        let values = robj.as_string_vector().ok_or_else(|| {
            format!(
                "mm_data_error: categorical values for column '{name}' is not a character vector ({:?})",
                robj.rtype()
            )
        })?;
        cat_values_map.insert(name.to_string(), values);
    }

    let mut cat_levels_map: HashMap<String, Vec<String>> =
        HashMap::with_capacity(categorical_levels.len());
    for (name, robj) in categorical_levels.iter() {
        let levels = robj.as_string_vector().ok_or_else(|| {
            format!(
                "mm_data_error: levels for categorical column '{name}' is not a character vector ({:?})",
                robj.rtype()
            )
        })?;
        cat_levels_map.insert(name.to_string(), levels);
    }

    // Names of categorical columns the R wrapper flagged as ordered factors.
    // These are coded with orthonormal polynomial contrasts (R's `contr.poly`)
    // instead of the default treatment coding, matching lme4's ordered-factor
    // behaviour. The R side guarantees each name here also appears in
    // `categorical_values`/`categorical_levels`.
    // NB: `Strings::iter()` aborts (non-unwinding panic) on a zero-length
    // vector, so guard the common no-ordered-factors case before iterating.
    // `.len()` is safe on an empty `Strings`; only `.iter()` is not.
    let ordered_set: HashSet<&str> = if categorical_ordered.len() == 0 {
        HashSet::new()
    } else {
        categorical_ordered
            .iter()
            .map(|name_rstr| name_rstr.as_ref())
            .collect()
    };

    let mut df = DataFrame::new();
    for name_rstr in column_order.iter() {
        // `Rstr::as_str` is deprecated; per upstream guidance, take the
        // `&str` view through `AsRef`. Deref coercion would also work
        // but the explicit annotation makes the intent obvious.
        let name: &str = name_rstr.as_ref();

        if let Some(values) = numeric_map.remove(name) {
            df.add_numeric(name, values).map_err(|e| {
                format!("mm_data_error: failed to add numeric column '{name}': {e}")
            })?;
            continue;
        }

        if let Some(values) = cat_values_map.remove(name) {
            let levels = cat_levels_map.remove(name).ok_or_else(|| {
                format!("mm_data_error: categorical column '{name}' is missing its `levels` entry")
            })?;
            if ordered_set.contains(name) {
                // Ordered factor: polynomial (contr.poly) contrast basis. The
                // engine builds the QR-orthonormalized Vandermonde over the
                // level order the R wrapper supplied and labels columns .L/.Q/.C/^k.
                let contrast = CategoricalContrast::polynomial(levels.clone()).map_err(|e| {
                    format!(
                        "mm_data_error: failed to build polynomial contrast for ordered column '{name}': {e}"
                    )
                })?;
                df.add_categorical_with_contrast(name, values, levels, contrast)
                    .map_err(|e| {
                        format!(
                            "mm_data_error: failed to add ordered categorical column '{name}': {e}"
                        )
                    })?;
            } else {
                df.add_categorical_with_levels(name, values, levels)
                    .map_err(|e| {
                        format!("mm_data_error: failed to add categorical column '{name}': {e}")
                    })?;
            }
            continue;
        }

        return Err(format!(
            "mm_data_error: column '{name}' was named in `column_order` but not supplied as numeric or categorical"
        ));
    }

    Ok(df)
}
