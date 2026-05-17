//! Translate R column specs into the upstream `mixedmodels::DataFrame`.
//!
//! The R wire format is three parallel lists keyed by column name:
//!   - `numeric_columns`: list(name = REAL vector)
//!   - `categorical_values`: list(name = character vector of observed values)
//!   - `categorical_levels`: list(name = character vector of canonical levels)
//! plus a `column_order` character vector that names every column in the
//! original `data.frame` order. Each column appears in exactly one of the
//! numeric/categorical lists.
//!
//! Rust trusts the R wrapper to have validated:
//!   - every value in `categorical_values[name]` is in `categorical_levels[name]`,
//!   - every name in `column_order` appears in exactly one of the column lists,
//!   - all columns have the same length (data.frame invariant).
//! Violations of these invariants produce typed `mm_data_error` strings the
//! R caller routes back to a typed condition.
//!
//! Caution: `DataFrame::add_categorical_with_levels` panics if any value is
//! not in `levels`. The R wrapper's `mm_translate_data()` derives `levels`
//! either from `levels(factor)` (so values ⊆ levels by construction) or from
//! `unique(character)` (same), so the panic path is unreachable through the
//! supported entry point.

use std::collections::HashMap;

use extendr_api::prelude::*;

use mixeff_rs::model::DataFrame;

pub(crate) fn build_dataframe(
    numeric_columns: &List,
    categorical_values: &List,
    categorical_levels: &List,
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

    let mut df = DataFrame::new();
    for name_rstr in column_order.iter() {
        // `Rstr::as_str` is deprecated; per upstream guidance, take the
        // `&str` view through `AsRef`. Deref coercion would also work
        // but the explicit annotation makes the intent obvious.
        let name: &str = name_rstr.as_ref();

        if let Some(values) = numeric_map.remove(name) {
            df.add_numeric(name, values);
            continue;
        }

        if let Some(values) = cat_values_map.remove(name) {
            let levels = cat_levels_map.remove(name).ok_or_else(|| {
                format!("mm_data_error: categorical column '{name}' is missing its `levels` entry")
            })?;
            df.add_categorical_with_levels(name, values, levels);
            continue;
        }

        return Err(format!(
            "mm_data_error: column '{name}' was named in `column_order` but not supplied as numeric or categorical"
        ));
    }

    Ok(df)
}
