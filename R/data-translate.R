#' Translate a base R data.frame to the wrapper's typed column wire format
#'
#' Internal: `mm_translate_data()` decomposes a `data.frame` into the
#' three parallel lists that `mm_compile_model_json()` (and later
#' `mm_fit_lmm_json()`) expect on the FFI: `numeric_columns`,
#' `categorical_values`, and `categorical_levels`. The output preserves
#' the original column order via the returned `column_order` character
#' vector.
#'
#' Type rules:
#' * `numeric` / `integer` / `logical` → `"numeric"` (logicals coerce to 0/1).
#' * `factor` → `"categorical"` with `levels = levels(col)` (canonical
#'   factor order, not first-appearance).
#' * `character` → `"categorical"` with `levels = unique(col)`
#'   (first-appearance order — matches the upstream `CategoricalColumn::new`
#'   default).
#' * Anything else (`Date`, `POSIXct`, lists, etc.) raises an
#'   `mm_data_error`.
#'
#' Validation contract: every value in a categorical column is in its
#' returned `levels` by construction (factor: `as.character(col) ⊆ levels(col)`;
#' character: levels are derived from values themselves). The Rust side
#' relies on this and uses the panicking
#' `DataFrame::add_categorical_with_levels` constructor; that panic path
#' is unreachable through the supported R entry point.
#'
#' @param data A base R `data.frame` (or anything that satisfies
#'   `is.data.frame()`).
#'
#' @return A list with four named components:
#' \describe{
#'   \item{`column_order`}{character — names of the columns in original order}
#'   \item{`numeric_columns`}{named list of numeric vectors, one per
#'     numeric/integer/logical column}
#'   \item{`categorical_values`}{named list of character vectors of observed
#'     values, one per factor/character column}
#'   \item{`categorical_levels`}{named list of character vectors of canonical
#'     levels, parallel to `categorical_values`}
#' }
#'
#' @keywords internal
#' @noRd
mm_translate_data <- function(data) {
  if (!is.data.frame(data)) {
    mm_abort(
      message = "`data` must be a data.frame.",
      class = "mm_data_error",
      input = data
    )
  }
  if (ncol(data) == 0L) {
    mm_abort(
      message = "`data` must have at least one column.",
      class = "mm_data_error",
      input = data
    )
  }

  col_names <- names(data)
  if (is.null(col_names) || any(!nzchar(col_names)) || anyDuplicated(col_names) > 0L) {
    mm_abort(
      message = "`data` must have unique, non-empty column names.",
      class = "mm_data_error",
      input = col_names
    )
  }

  numeric_columns    <- list()
  categorical_values <- list()
  categorical_levels <- list()

  for (nm in col_names) {
    col <- data[[nm]]
    if (is.factor(col)) {
      categorical_values[[nm]] <- as.character(col)
      categorical_levels[[nm]] <- as.character(levels(col))
    } else if (is.character(col)) {
      vals <- as.character(col)
      categorical_values[[nm]] <- vals
      categorical_levels[[nm]] <- unique(vals)
    } else if (is.logical(col) || is.integer(col) || is.numeric(col)) {
      numeric_columns[[nm]] <- as.numeric(col)
    } else {
      mm_abort(
        message = sprintf(
          "Unsupported column type for `%s`: %s. mixeff v0.1 accepts numeric, integer, logical, factor, and character columns.",
          nm, paste(class(col), collapse = "/")
        ),
        class = "mm_data_error",
        column = nm,
        input = col
      )
    }
  }

  list(
    column_order       = col_names,
    numeric_columns    = numeric_columns,
    categorical_values = categorical_values,
    categorical_levels = categorical_levels
  )
}

#' Refuse compilation when a design variable contains NA
#'
#' mixeff is no-silent-surgery: it does not silently drop rows with NA
#' in design variables. Users must `na.omit()` (or more carefully select
#' complete cases) explicitly before handing data to `compile_model()` /
#' `lmm()`. This helper enumerates the offending columns and raises one
#' typed `mm_data_error` listing all of them rather than failing on the
#' first.
#'
#' @keywords internal
#' @noRd
mm_check_no_na <- function(data, vars, .call = rlang::caller_env()) {
  na_counts <- vapply(
    vars,
    function(v) sum(is.na(data[[v]])),
    integer(1)
  )
  affected <- vars[na_counts > 0L]
  if (length(affected) == 0L) {
    return(invisible(TRUE))
  }
  details <- paste(
    sprintf("`%s` (%d NA)", affected, na_counts[na_counts > 0L]),
    collapse = ", "
  )
  mm_abort(
    message = sprintf(
      "Missing values in design variable(s): %s. mixeff requires complete cases; pass na.omit(data) explicitly before fitting.",
      details
    ),
    class = "mm_data_error",
    columns = affected,
    na_counts = unname(na_counts[na_counts > 0L]),
    call = .call
  )
}
