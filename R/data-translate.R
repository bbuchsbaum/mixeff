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
#' * `ordered` factor → additionally recorded in `categorical_ordered` so the
#'   engine codes it with orthonormal polynomial contrasts (`contr.poly`),
#'   matching lme4's ordered-factor behaviour rather than treatment coding.
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
#' @return A list with five named components:
#' \describe{
#'   \item{`column_order`}{character — names of the columns in original order}
#'   \item{`numeric_columns`}{named list of numeric vectors, one per
#'     numeric/integer/logical column}
#'   \item{`categorical_values`}{named list of character vectors of observed
#'     values, one per factor/character column}
#'   \item{`categorical_levels`}{named list of character vectors of canonical
#'     levels, parallel to `categorical_values`}
#'   \item{`categorical_ordered`}{character — names of the categorical columns
#'     that are ordered factors (coded with `contr.poly`); a subset of the
#'     names in `categorical_values`, possibly empty}
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
  categorical_ordered <- character(0)

  for (nm in col_names) {
    col <- data[[nm]]
    if (is.factor(col)) {
      categorical_values[[nm]] <- as.character(col)
      categorical_levels[[nm]] <- as.character(levels(col))
      if (is.ordered(col)) {
        mm_assert_ordered_contrast_policy(nm, col)
        categorical_ordered <- c(categorical_ordered, nm)
      }
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
    column_order        = col_names,
    numeric_columns     = numeric_columns,
    categorical_values  = categorical_values,
    categorical_levels  = categorical_levels,
    categorical_ordered = categorical_ordered
  )
}

#' Refuse ordered factors whose contrast policy mixeff cannot honour
#'
#' Ordered factors are coded with `contr.poly` to match lme4/R's default. That
#' parity only holds when the ordered contrast is R's default: if the caller
#' has switched the global ordered-contrast option away from `contr.poly`, or
#' attached an explicit non-poly `contrasts` attribute to the factor, applying
#' `contr.poly` anyway would silently diverge from the coding the user (and
#' lme4) would expect. Per the no-silent-surgery contract we refuse with a
#' typed `mm_arg_error` rather than pick a coding behind the user's back.
#'
#' @keywords internal
#' @noRd
mm_assert_ordered_contrast_policy <- function(nm, col, .call = rlang::caller_env()) {
  # R resolves the ordered-factor contrast option POSITIONALLY (element 2 is
  # the ordered coding), honouring names only when present. The standard form
  # `options(contrasts = c("contr.treatment", "contr.poly"))` is unnamed, so
  # reading it by name (`[["ordered"]]`) would throw "subscript out of bounds"
  # and abort a correctly poly-coded fit. Mirror R: prefer the named entry,
  # else fall back to the second element.
  opt <- getOption("contrasts")
  global_ordered <- if (!is.null(names(opt)) && "ordered" %in% names(opt)) {
    opt[["ordered"]]
  } else if (length(opt) >= 2L) {
    as.character(opt)[[2L]]
  } else {
    NA_character_
  }
  if (!identical(as.character(global_ordered), "contr.poly")) {
    shown <- if (length(global_ordered) != 1L || is.na(global_ordered)) {
      "not contr.poly"
    } else {
      as.character(global_ordered)
    }
    mm_abort(
      message = sprintf(
        paste0(
          "Ordered factor `%s` is coded with `contr.poly` (to match lme4), but the ",
          "global ordered-factor contrast option resolves to `%s`. Reset it with ",
          "options(contrasts = c(\"contr.treatment\", \"contr.poly\")), or drop the ",
          "ordering with `%s <- factor(%s, ordered = FALSE)` for treatment coding."
        ),
        nm, shown, nm, nm
      ),
      class = "mm_arg_error",
      column = nm,
      call = .call
    )
  }
  # An ordered factor is honoured only when its explicit contrasts, if any, is
  # the string "contr.poly". An attached numeric contrast MATRIX is rejected
  # even when it numerically equals contr.poly, because we cannot cheaply
  # certify equivalence and will not silently substitute our own basis.
  col_contrasts <- attr(col, "contrasts")
  if (!is.null(col_contrasts) && !identical(col_contrasts, "contr.poly")) {
    mm_abort(
      message = sprintf(
        paste0(
          "Ordered factor `%s` carries an explicit `contrasts` attribute that is not ",
          "the string \"contr.poly\" (an attached contrast MATRIX is rejected even when ",
          "it numerically equals contr.poly). mixeff codes ordered factors with ",
          "contr.poly; set `contrasts(%s) <- \"contr.poly\"` (string form) or drop the ",
          "attribute with `attr(%s, \"contrasts\") <- NULL`."
        ),
        nm, nm, nm
      ),
      class = "mm_arg_error",
      column = nm,
      call = .call
    )
  }
  invisible(TRUE)
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
