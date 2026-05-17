#' Declare or inspect design roles
#'
#' `roles()` has two Phase 1.F uses. With named string arguments it constructs
#' a declared-role object, e.g. `roles(subject = "sampled_unit")`. With a
#' compiled spec or fit as its only unnamed argument, it returns the observed
#' role fallback inferred from the formula and model frame.
#'
#' @param ... Either named role strings or one unnamed `mm_spec` / `mm_fit`.
#'
#' @return An `mm_roles` object with a data-frame `table`.
#'
#' @export
roles <- function(...) {
  dots <- list(...)
  dot_names <- names(dots) %||% rep("", length(dots))

  if (length(dots) == 1L && !nzchar(dot_names[[1L]]) &&
      inherits(dots[[1L]], "mm_compiled")) {
    return(mm_roles_from_compiled(dots[[1L]]))
  }

  if (!length(dots)) {
    out <- list(table = mm_roles_empty_table(), source = "declared")
    class(out) <- "mm_roles"
    return(out)
  }

  if (is.null(names(dots)) || any(!nzchar(names(dots)))) {
    mm_abort(
      message = "`roles()` declarations must be named character strings.",
      class = "mm_arg_error",
      input = dots
    )
  }
  rows <- lapply(names(dots), function(nm) {
    role <- dots[[nm]]
    if (!is.character(role) || length(role) != 1L || is.na(role) || !nzchar(role)) {
      mm_abort(
        message = sprintf("Role for `%s` must be a single non-empty string.", nm),
        class = "mm_arg_error",
        input = role
      )
    }
    data.frame(
      variable = nm,
      role = role,
      origin = "declared_by_user",
      group = "",
      evidence = "",
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  out <- list(table = table, source = "declared")
  class(out) <- "mm_roles"
  out
}

#' @method print mm_roles
#' @export
print.mm_roles <- function(x, ...) {
  cat("Design roles:\n")
  if (!nrow(x$table)) {
    cat("  none declared\n")
    return(invisible(x))
  }
  print(x$table, row.names = FALSE)
  invisible(x)
}

mm_roles_from_compiled <- function(x) {
  artifact <- mm_compiled_artifact(x)
  frame <- x$model_frame %||% data.frame()
  random_terms <- artifact$semantic_model$random_terms %||% list()
  groups <- unique(vapply(random_terms, function(term) {
    mm_group_label(term$group)
  }, character(1)))
  groups <- groups[nzchar(groups) & groups != "(unknown)"]

  rows <- list()
  for (group in groups) {
    rows[[length(rows) + 1L]] <- data.frame(
      variable = group,
      role = "sampled_unit",
      origin = "observed_from_data",
      group = group,
      evidence = "appears as a grouping factor",
      stringsAsFactors = FALSE
    )
  }

  fixed <- unlist(artifact$semantic_model$fixed_terms %||% list(), use.names = FALSE)
  fixed <- unique(setdiff(fixed, "1"))
  fixed <- fixed[fixed %in% names(frame)]
  for (variable in fixed) {
    observed <- mm_observed_fixed_role(variable, groups, frame)
    rows[[length(rows) + 1L]] <- data.frame(
      variable = variable,
      role = observed$role,
      origin = "observed_from_data",
      group = observed$group,
      evidence = observed$evidence,
      stringsAsFactors = FALSE
    )
  }

  table <- if (length(rows)) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  } else {
    mm_roles_empty_table()
  }
  obj <- list(table = table, source = "observed_from_data")
  class(obj) <- "mm_roles"
  obj
}

mm_observed_fixed_role <- function(variable, groups, frame) {
  if (!length(groups) || !variable %in% names(frame)) {
    return(list(
      role = "observed_fixed_effect",
      group = "",
      evidence = "fixed-effect term observed in data"
    ))
  }
  usable_groups <- groups[groups %in% names(frame)]
  for (group in usable_groups) {
    split_values <- split(frame[[variable]], frame[[group]], drop = TRUE)
    unique_counts <- vapply(split_values, function(x) length(unique(x)), integer(1))
    if (any(unique_counts > 1L)) {
      return(list(
        role = "observed_within_group",
        group = group,
        evidence = sprintf("varies within `%s`", group)
      ))
    }
  }
  if (length(usable_groups)) {
    return(list(
      role = "observed_between_group",
      group = usable_groups[[1L]],
      evidence = sprintf("constant within `%s`", usable_groups[[1L]])
    ))
  }
  list(
    role = "observed_fixed_effect",
    group = "",
    evidence = "fixed-effect term observed in data"
  )
}

mm_roles_empty_table <- function() {
  data.frame(
    variable = character(),
    role = character(),
    origin = character(),
    group = character(),
    evidence = character(),
    stringsAsFactors = FALSE
  )
}
