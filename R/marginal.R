#' Marginal grids, predictions, means, and comparisons
#'
#' These helpers provide a small native marginal-quantities surface for
#' Gaussian LMM fits. They cover the common population-level workflow:
#' construct a reference grid, evaluate fixed-effect predictions, average them
#' into marginal means, and compare those means by simple differences.
#'
#' The returned tables use the `mixedmodels.marginal_quantity_table` row
#' contract. Inference is routed through [contrast()] so rows retain the same
#' method, status, reliability, estimability, and reason fields as fixed-effect
#' contrasts. Ordinary full-rank LMMs use the versioned
#' `mixedmodels.fixed_effect_covariance_matrix` payload for fixed-effect
#' uncertainty; rank-deficient or otherwise uncertified fits surface explicit
#' unavailable status and reasons instead of partial covariance numbers.
#'
#' @param fit A fitted `mm_lmm`.
#' @param specs Character vector, or a one-sided formula such as `~ trt` or
#'   `~ trt | group`, naming the displayed marginal dimensions.
#' @param by Optional character vector of grouping variables for marginal
#'   summaries or pairwise comparisons.
#' @param at Named list of fixed-predictor values to force in the grid.
#' @param cov.reduce Function used to reduce numeric fixed predictors that are
#'   not explicitly gridded.
#' @param grid Optional object returned by `mm_grid()`.
#' @param method Requested inference method, passed to [contrast()].
#' @param level Confidence level for intervals computed from contrast standard
#'   errors.
#' @param weights Averaging weights for `mm_means()` and `mm_comparisons()`.
#'   `"equal"` weights reference-grid cells equally; `"proportional"` weights
#'   cells by observed fixed-factor frequencies.
#' @param comparison Comparison scale. Only `"difference"` is implemented.
#' @param target Prediction target. Only `"population"` is implemented.
#' @param scale Prediction scale. Gaussian LMMs have identical `"link"` and
#'   `"response"` scales.
#' @param ... Reserved for future methods.
#'
#' @return `mm_grid()` returns an `mm_grid` object. The other helpers return an
#'   `mm_marginal_quantity` object with a contract-shaped `table`.
#'
#' @export
mm_grid <- function(fit, specs, by = NULL, at = list(), cov.reduce = mean, ...) {
  UseMethod("mm_grid")
}

#' @rdname mm_grid
#' @export
mm_grid.mm_lmm <- function(fit, specs, by = NULL, at = list(),
                           cov.reduce = mean, ...) {
  parsed <- mm_parse_marginal_specs(specs, by)
  vars <- mm_fixed_predictor_vars(fit)
  unknown <- setdiff(unique(c(parsed$specs, parsed$by, names(at))), vars)
  if (length(unknown)) {
    mm_abort(
      message = sprintf("Unknown fixed-effect predictor(s): %s.",
                        paste(unknown, collapse = ", ")),
      class = "mm_arg_error",
      input = unknown
    )
  }

  if (!is.list(at) ||
      (length(at) && (is.null(names(at)) || any(!nzchar(names(at)))))) {
    mm_abort(
      message = "`at` must be a named list of grid values.",
      class = "mm_arg_error",
      input = at
    )
  }

  frame <- fit$model_frame
  grid_vars <- mm_reference_grid_vars(fit, parsed$specs, parsed$by, at)
  values <- lapply(grid_vars, function(var) {
    if (var %in% names(at)) {
      return(mm_grid_values_like(frame[[var]], at[[var]], var))
    }
    mm_default_grid_values(frame[[var]], var, var %in% c(parsed$specs, parsed$by),
                           cov.reduce)
  })
  names(values) <- grid_vars
  grid <- do.call(expand.grid, c(values, KEEP.OUT.ATTRS = FALSE,
                                stringsAsFactors = FALSE))
  grid <- mm_restore_grid_classes(grid, frame)

  X <- mm_fixed_basis(fit, grid)
  grid_id <- paste0("grid", seq_len(nrow(grid)))
  out <- list(
    grid = grid,
    X = X,
    specs = parsed$specs,
    by = parsed$by,
    at = at,
    grid_id = grid_id,
    factor_vars = mm_fixed_factor_vars(fit),
    numeric_vars = mm_fixed_numeric_vars(fit)
  )
  class(out) <- "mm_grid"
  out
}

#' @rdname mm_grid
#' @export
mm_predictions <- function(fit, grid = NULL, specs = NULL, by = NULL,
                           at = list(),
                           method = c("auto", "satterthwaite",
                                      "kenward_roger", "bootstrap",
                                      "asymptotic", "none"),
                           level = 0.95,
                           target = c("population"),
                           scale = c("response", "link"), ...) {
  UseMethod("mm_predictions")
}

#' @rdname mm_grid
#' @export
mm_predictions.mm_lmm <- function(fit, grid = NULL, specs = NULL, by = NULL,
                                  at = list(),
                                  method = c("auto", "satterthwaite",
                                             "kenward_roger", "bootstrap",
                                             "asymptotic", "none"),
                                  level = 0.95,
                                  target = c("population"),
                                  scale = c("response", "link"), ...) {
  method <- match.arg(method)
  target <- match.arg(target)
  scale <- match.arg(scale)
  level <- mm_validate_marginal_level(level)
  grid <- mm_resolve_grid(fit, grid, specs, by, at)

  ct <- contrast(fit, grid$X, method = method)
  rows <- mm_marginal_rows_from_contrast(
    ct$table,
    quantity = "prediction",
    labels = mm_grid_labels(grid$grid, names(grid$grid)),
    rhs = rep(0, nrow(grid$X)),
    level = level,
    target = target,
    scale = scale,
    weights = "identity",
    comparison = "identity",
    by = grid$by,
    specs = grid$specs,
    grid_id = grid$grid_id,
    details_extra = lapply(seq_len(nrow(grid$grid)), function(i) {
      list(grid = as.list(grid$grid[i, , drop = FALSE]))
    })
  )
  mm_new_marginal_quantity(rows, grid = grid, L = grid$X,
                           requested_method = method)
}

#' @rdname mm_grid
#' @export
mm_means <- function(fit, specs, by = NULL, at = list(), grid = NULL,
                     method = c("auto", "satterthwaite",
                                "kenward_roger", "bootstrap",
                                "asymptotic", "none"),
                     level = 0.95,
                     weights = c("equal", "proportional"),
                     target = c("population"),
                     scale = c("response", "link"), ...) {
  UseMethod("mm_means")
}

#' @rdname mm_grid
#' @export
mm_means.mm_lmm <- function(fit, specs, by = NULL, at = list(), grid = NULL,
                            method = c("auto", "satterthwaite",
                                       "kenward_roger", "bootstrap",
                                       "asymptotic", "none"),
                            level = 0.95,
                            weights = c("equal", "proportional"),
                            target = c("population"),
                            scale = c("response", "link"), ...) {
  method <- match.arg(method)
  weights <- match.arg(weights)
  target <- match.arg(target)
  scale <- match.arg(scale)
  level <- mm_validate_marginal_level(level)
  grid <- mm_resolve_grid(fit, grid, specs, by, at)

  groups <- unique(grid$grid[, c(grid$specs, grid$by), drop = FALSE])
  groups <- groups[do.call(order, groups), , drop = FALSE]
  L <- mm_group_basis(fit, grid, groups, weights)
  rownames(L) <- mm_grid_labels(groups, c(grid$specs, grid$by))
  ct <- contrast(fit, L, method = method)
  rows <- mm_marginal_rows_from_contrast(
    ct$table,
    quantity = "mean",
    labels = rownames(L),
    rhs = rep(0, nrow(L)),
    level = level,
    target = target,
    scale = scale,
    weights = weights,
    comparison = "identity",
    by = grid$by,
    specs = grid$specs,
    grid_id = paste0("mean", seq_len(nrow(L))),
    details_extra = lapply(seq_len(nrow(groups)), function(i) {
      list(grid = as.list(groups[i, , drop = FALSE]))
    })
  )
  mm_new_marginal_quantity(rows, grid = grid, L = L,
                           requested_method = method)
}

#' @rdname mm_grid
#' @export
mm_comparisons <- function(fit, specs, by = NULL, at = list(), grid = NULL,
                           comparison = c("difference", "ratio", "odds_ratio"),
                           method = c("auto", "satterthwaite",
                                      "kenward_roger", "bootstrap",
                                      "asymptotic", "none"),
                           level = 0.95,
                           weights = c("equal", "proportional"),
                           target = c("population"),
                           scale = c("response", "link"), ...) {
  UseMethod("mm_comparisons")
}

#' @rdname mm_grid
#' @export
mm_comparisons.mm_lmm <- function(fit, specs, by = NULL, at = list(),
                                  grid = NULL,
                                  comparison = c("difference", "ratio", "odds_ratio"),
                                  method = c("auto", "satterthwaite",
                                             "kenward_roger", "bootstrap",
                                             "asymptotic", "none"),
                                  level = 0.95,
                                  weights = c("equal", "proportional"),
                                  target = c("population"),
                                  scale = c("response", "link"), ...) {
  comparison <- mm_match_marginal_comparison(comparison)
  method <- match.arg(method)
  weights <- match.arg(weights)
  target <- match.arg(target)
  scale <- match.arg(scale)
  level <- mm_validate_marginal_level(level)
  grid <- mm_resolve_grid(fit, grid, specs, by, at)

  means <- mm_means(fit, specs = grid$specs, by = grid$by, grid = grid,
                    method = "none", weights = weights, target = target,
                    scale = scale)
  L_means <- means$L
  group_frame <- mm_marginal_group_frame(means$table, grid$specs, grid$by)
  pairs <- mm_pairwise_rows(group_frame, grid$specs, grid$by)
  if (!nrow(pairs)) {
    mm_abort(
      message = "`mm_comparisons()` needs at least two marginal means per comparison group.",
      class = "mm_inference_unavailable",
      input = specs
    )
  }

  L <- L_means[pairs$left, , drop = FALSE] - L_means[pairs$right, , drop = FALSE]
  labels <- paste(rownames(L_means)[pairs$left], "-", rownames(L_means)[pairs$right])
  rownames(L) <- labels
  ct <- contrast(fit, L, method = method)
  rows <- mm_marginal_rows_from_contrast(
    ct$table,
    quantity = "comparison",
    labels = labels,
    rhs = rep(0, nrow(L)),
    level = level,
    target = target,
    scale = scale,
    weights = weights,
    comparison = comparison,
    by = grid$by,
    specs = grid$specs,
    grid_id = paste0("comparison", seq_len(nrow(L))),
    details_extra = lapply(seq_len(nrow(pairs)), function(i) {
      list(
        left = rownames(L_means)[pairs$left[[i]]],
        right = rownames(L_means)[pairs$right[[i]]]
      )
    })
  )
  mm_new_marginal_quantity(rows, grid = grid, L = L,
                           requested_method = method)
}

mm_match_marginal_comparison <- function(comparison) {
  comparison <- match.arg(comparison, c("difference", "ratio", "odds_ratio"))
  if (!identical(comparison, "difference")) {
    mm_abort(
      message = sprintf(
        "`comparison = \"%s\"` is not supported; only `\"difference\"` is implemented for marginal quantities.",
        comparison
      ),
      class = "mm_inference_unavailable",
      input = comparison
    )
  }
  comparison
}

#' @export
print.mm_grid <- function(x, ...) {
  cat("Marginal grid:\n")
  print(x$grid, row.names = FALSE)
  invisible(x)
}

#' @export
print.mm_marginal_quantity <- function(x, ...) {
  cat("Marginal quantities:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

mm_resolve_grid <- function(fit, grid, specs, by, at) {
  if (is.null(grid)) {
    if (is.null(specs)) {
      mm_abort(
        message = "`specs` is required when `grid` is not supplied.",
        class = "mm_arg_error",
        input = specs
      )
    }
    return(mm_grid(fit, specs = specs, by = by, at = at))
  }
  if (!inherits(grid, "mm_grid")) {
    mm_abort(
      message = "`grid` must be an object returned by `mm_grid()`.",
      class = "mm_inference_unavailable",
      input = grid
    )
  }
  grid
}

mm_parse_marginal_specs <- function(specs, by = NULL) {
  formula_by <- character()
  if (inherits(specs, "formula")) {
    rhs <- specs[[length(specs)]]
    if (is.call(rhs) && identical(as.character(rhs[[1L]]), "|")) {
      formula_by <- all.vars(rhs[[3L]])
      specs <- all.vars(rhs[[2L]])
    } else {
      specs <- all.vars(rhs)
    }
  }
  specs <- unique(as.character(specs %||% character()))
  by <- unique(c(formula_by, as.character(by %||% character())))
  if (!length(specs)) {
    mm_abort(
      message = "`specs` must name at least one fixed-effect predictor.",
      class = "mm_arg_error",
      input = specs
    )
  }
  list(specs = specs, by = by)
}

mm_fixed_predictor_vars <- function(fit) {
  all.vars(stats::delete.response(stats::terms(mm_fixed_formula(fit))))
}

mm_fixed_factor_vars <- function(fit) {
  vars <- mm_fixed_predictor_vars(fit)
  vars[vapply(fit$model_frame[vars], function(x) is.factor(x) || is.character(x),
              logical(1))]
}

mm_fixed_numeric_vars <- function(fit) {
  vars <- mm_fixed_predictor_vars(fit)
  vars[vapply(fit$model_frame[vars], is.numeric, logical(1))]
}

mm_reference_grid_vars <- function(fit, specs, by, at) {
  vars <- mm_fixed_predictor_vars(fit)
  unique(c(specs, by, names(at), mm_fixed_factor_vars(fit), vars))
}

mm_grid_values_like <- function(template, values, var) {
  if (!length(values)) {
    mm_abort(
      message = sprintf("`at$%s` must contain at least one value.", var),
      class = "mm_arg_error",
      input = values
    )
  }
  if (is.factor(template)) {
    unknown <- setdiff(as.character(values), levels(template))
    if (length(unknown)) {
      mm_abort(
        message = sprintf("Unknown level(s) for `%s`: %s.",
                          var, paste(unknown, collapse = ", ")),
        class = "mm_arg_error",
        input = values
      )
    }
    return(factor(as.character(values), levels = levels(template)))
  }
  values
}

mm_default_grid_values <- function(x, var, displayed, cov.reduce) {
  if (is.factor(x)) {
    return(factor(levels(x), levels = levels(x)))
  }
  if (is.character(x)) {
    return(sort(unique(x)))
  }
  if (is.logical(x)) {
    return(c(FALSE, TRUE))
  }
  if (is.numeric(x)) {
    if (displayed) {
      ux <- sort(unique(x))
      if (length(ux) <= 10L) {
        return(ux)
      }
    }
    value <- cov.reduce(x)
    if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
      mm_abort(
        message = sprintf("`cov.reduce` must return one non-missing number for `%s`.", var),
        class = "mm_arg_error",
        input = value
      )
    }
    return(as.numeric(value))
  }
  mm_abort(
    message = sprintf("Unsupported fixed predictor type for `%s`.", var),
    class = "mm_inference_unavailable",
    input = x
  )
}

mm_restore_grid_classes <- function(grid, frame) {
  for (nm in intersect(names(grid), names(frame))) {
    if (is.factor(frame[[nm]])) {
      grid[[nm]] <- factor(as.character(grid[[nm]]), levels = levels(frame[[nm]]))
    }
  }
  grid
}

mm_fixed_basis <- function(fit, grid) {
  # Reconstruct the reference-grid design in the engine's coefficient basis
  # (treatment contrasts for all factors, columns aligned to names(beta) by
  # name). Using R's default contrasts here silently mis-evaluated marginal
  # means for ordered-factor and interaction models. See
  # mm_engine_fixed_matrix() in predict.R.
  mm_engine_fixed_matrix(fit, grid)
}

mm_grid_labels <- function(grid, vars) {
  vars <- unique(vars)
  if (!length(vars)) {
    return(rep("(overall)", nrow(grid)))
  }
  apply(grid[, vars, drop = FALSE], 1L, function(row) {
    paste(paste(names(row), as.character(row), sep = "="), collapse = ", ")
  })
}

mm_group_basis <- function(fit, grid, groups, weights) {
  out <- matrix(NA_real_, nrow = nrow(groups), ncol = ncol(grid$X))
  colnames(out) <- colnames(grid$X)
  for (i in seq_len(nrow(groups))) {
    keep <- mm_rows_match(grid$grid, groups[i, , drop = FALSE])
    w <- mm_cell_weights(fit, grid, keep, weights)
    out[i, ] <- colSums(grid$X[keep, , drop = FALSE] * w)
  }
  out
}

mm_rows_match <- function(grid, values) {
  keep <- rep(TRUE, nrow(grid))
  for (nm in names(values)) {
    keep <- keep & as.character(grid[[nm]]) == as.character(values[[nm]][[1L]])
  }
  keep
}

mm_cell_weights <- function(fit, grid, keep, weights) {
  n <- sum(keep)
  if (!n) {
    return(numeric())
  }
  if (identical(weights, "equal")) {
    return(rep(1 / n, n))
  }
  factors <- intersect(mm_fixed_factor_vars(fit), names(grid$grid))
  if (!length(factors)) {
    return(rep(1 / n, n))
  }
  cells <- grid$grid[keep, factors, drop = FALSE]
  counts <- apply(cells, 1L, function(row) {
    obs_keep <- rep(TRUE, nrow(fit$model_frame))
    for (nm in factors) {
      obs_keep <- obs_keep & as.character(fit$model_frame[[nm]]) == as.character(row[[nm]])
    }
    sum(obs_keep)
  })
  if (!sum(counts)) {
    return(rep(1 / n, n))
  }
  counts / sum(counts)
}

mm_marginal_rows_from_contrast <- function(table, quantity, labels, rhs, level,
                                           target, scale, weights, comparison,
                                           by, specs, grid_id,
                                           details_extra = NULL) {
  n <- nrow(table)
  intervals <- mm_marginal_intervals(table, level)
  details_extra <- details_extra %||% rep(list(NULL), n)
  details <- Map(function(base, extra) {
    c(base %||% list(), extra %||% list())
  }, table$details %||% rep(list(NULL), n), details_extra)

  data.frame(
    quantity = rep(quantity, n),
    label = labels,
    estimate = table$estimate,
    rhs = rhs,
    std_error = table$std_error,
    df = table$df,
    statistic = table$statistic,
    statistic_name = table$statistic_name,
    p_value = table$p_value,
    conf_low = intervals[, 1L],
    conf_high = intervals[, 2L],
    method = table$method,
    requested_method = table$requested_method,
    status = table$status,
    reliability = table$reliability,
    estimability = table$estimability,
    reason = table$reason,
    target = rep(target, n),
    scale = rep(scale, n),
    weights = rep(weights, n),
    comparison = rep(comparison, n),
    by = I(rep(list(by), n)),
    specs = I(rep(list(specs), n)),
    grid_id = grid_id,
    details = I(details),
    notes = table$notes,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

mm_marginal_intervals <- function(table, level) {
  out <- matrix(NA_real_, nrow = nrow(table), ncol = 2L)
  alpha <- 1 - level
  for (i in seq_len(nrow(table))) {
    se <- table$std_error[[i]]
    est <- table$estimate[[i]]
    if (!is.finite(se) || !is.finite(est)) {
      next
    }
    df <- table$df[[i]]
    crit <- if (is.finite(df) && df > 0) {
      stats::qt(1 - alpha / 2, df = df)
    } else {
      stats::qnorm(1 - alpha / 2)
    }
    out[i, ] <- c(est - crit * se, est + crit * se)
  }
  out
}

mm_new_marginal_quantity <- function(table, grid, L, requested_method) {
  obj <- list(
    table = table,
    grid = grid,
    L = L,
    requested_method = requested_method,
    schema = list(
      schema_name = "mixedmodels.marginal_quantity_table",
      schema_version = "1.0.0"
    )
  )
  class(obj) <- "mm_marginal_quantity"
  obj
}

mm_validate_marginal_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    mm_abort(
      message = "`level` must be a single number between 0 and 1.",
      class = "mm_arg_error",
      input = level
    )
  }
  level
}

mm_marginal_group_frame <- function(table, specs, by) {
  vars <- unique(c(specs, by))
  pieces <- strsplit(table$label, ", ", fixed = TRUE)
  out <- lapply(pieces, function(piece) {
    values <- setNames(rep(NA_character_, length(vars)), vars)
    for (part in piece) {
      split <- strsplit(part, "=", fixed = TRUE)[[1L]]
      if (length(split) == 2L && split[[1L]] %in% vars) {
        values[[split[[1L]]]] <- split[[2L]]
      }
    }
    values
  })
  as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
}

mm_pairwise_rows <- function(groups, specs, by) {
  by <- unique(by)
  split_key <- if (length(by)) {
    apply(groups[, by, drop = FALSE], 1L, paste, collapse = "\r")
  } else {
    rep("all", nrow(groups))
  }
  idx <- split(seq_len(nrow(groups)), split_key)
  pairs <- do.call(rbind, lapply(idx, function(ii) {
    if (length(ii) < 2L) {
      return(NULL)
    }
    cmb <- utils::combn(ii, 2L)
    data.frame(left = cmb[2L, ], right = cmb[1L, ])
  }))
  if (is.null(pairs)) {
    return(data.frame(left = integer(), right = integer()))
  }
  rownames(pairs) <- NULL
  pairs
}
