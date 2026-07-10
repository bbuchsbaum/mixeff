#' Revive a serialized mixeff object
#'
#' `revive()` restores the process-local parts of a `mixeff` object after
#' `saveRDS()` / `readRDS()` or a worker restart. The fitted artifact and flat
#' extractor values are the durable source of truth; the Rust handle is only a
#' cache and may be absent. In the current bridge, revival recreates the lazy
#' R-side cache and explicitly leaves `rust_handle = NULL`.
#'
#' @param fit A fitted `mm_fit` object.
#' @param ... Reserved for future methods.
#'
#' @return A revived `mm_fit` object.
#'
#' @export
revive <- function(fit, ...) {
  UseMethod("revive")
}

#' @rdname revive
#' @export
revive.mm_fit <- function(fit, ...) {
  if (!is.list(fit)) {
    mm_abort(
      message = "`fit` must be a fitted mixeff object.",
      class = "mm_arg_error",
      input = fit
    )
  }

  if (!is.list(fit$artifact) && is.character(fit$fit$artifact_json)) {
    fit$artifact <- mm_json_parse_artifact(fit$fit$artifact_json)
  }
  if (!is.list(fit$artifact)) {
    mm_abort(
      message = "`fit` does not carry a parsed compiler artifact.",
      class = "mm_arg_error",
      input = fit
    )
  }
  if (is.null(attr(fit$artifact, "raw_json")) &&
      is.character(fit$fit$artifact_json) &&
      length(fit$fit$artifact_json) == 1L &&
      nzchar(fit$fit$artifact_json)) {
    attr(fit$artifact, "raw_json") <- fit$fit$artifact_json
  }

  fit$schema <- fit$schema %||% mm_object_schema(fit$artifact)
  fit$rust_handle <- NULL
  fit$lazy_cache <- mm_empty_lazy_cache()
  class(fit) <- unique(c(class(fit), "mm_fit", "mm_compiled"))
  fit
}

#' @export
revive.default <- function(fit, ...) {
  mm_abort(
    message = "`revive()` expects a fitted mixeff object.",
    class = "mm_arg_error",
    input = fit
  )
}

#' Test whether a mixeff fit has a live native handle
#'
#' The native handle is a process-local cache. A `FALSE` result does not mean
#' the fit is unusable: Phase 2 extractors read from the durable artifact and
#' flat R-side payload, and [revive()] recreates the lazy cache after
#' serialization.
#'
#' @param fit A fitted `mm_fit` object.
#' @param ... Reserved for future methods.
#'
#' @return A length-one logical value.
#'
#' @export
fit_handle_alive <- function(fit, ...) {
  UseMethod("fit_handle_alive")
}

#' @rdname fit_handle_alive
#' @export
fit_handle_alive.mm_fit <- function(fit, ...) {
  ptr <- fit$rust_handle
  !is.null(ptr) && identical(typeof(ptr), "externalptr")
}

#' @export
fit_handle_alive.default <- function(fit, ...) {
  FALSE
}

#' Extract low-level model components
#'
#' `getME()` provides a small, honest subset of the familiar lme4 extractor.
#' The fixed-effect design (`"X"`), random-effect design (`"Z"`), relative
#' covariance factor (`"Lambda"` / `"Lambdat"`), grouping factors (`"flist"`),
#' random coefficient names (`"cnms"`), response (`"y"`), fixed coefficients
#' (`"beta"` / `"fixef"`), and theta vector (`"theta"`) are rebuilt lazily
#' from the serialized R object.
#'
#' @param object A fitted `mm_lmm` object.
#' @param name Component name, or a character vector of names.
#' @param ... Reserved for future methods.
#'
#' @return The requested component, or a named list for multiple names.
#'
#' @export
getME <- function(object, name, ...) {
  UseMethod("getME")
}

#' @rdname getME
#' @export
getME.mm_lmm <- function(object, name, ...) {
  if (missing(name) || !is.character(name) || !length(name)) {
    mm_abort(
      message = "`name` must be a non-empty character vector.",
      class = "mm_arg_error",
      input = name
    )
  }

  values <- lapply(name, function(one) {
    switch(
      one,
      X = stats::model.matrix(object, type = "fixed"),
      Z = stats::model.matrix(object, type = "random"),
      Zt = Matrix::t(stats::model.matrix(object, type = "random")),
      Lambda = .mm_lazy(object, "Lambda", mm_lambda_matrix),
      Lambdat = Matrix::t(.mm_lazy(object, "Lambda", mm_lambda_matrix)),
      theta = object$theta,
      beta = object$beta,
      fixef = object$beta,
      y = mm_response_vector(object),
      mu = fitted(object),
      flist = .mm_lazy(object, "flist", mm_random_flist),
      cnms = .mm_lazy(object, "cnms", mm_random_cnms),
      mm_abort(
        message = sprintf("`getME()` component `%s` is not available.", one),
        class = "mm_arg_error",
        input = one
      )
    )
  })
  names(values) <- name
  if (length(values) == 1L) values[[1L]] else values
}

#' @export
getME.default <- function(object, name, ...) {
  mm_abort(
    message = "`getME()` expects a fitted mixeff LMM.",
    class = "mm_arg_error",
    input = object
  )
}

#' @rdname mm_lmm-methods
#' @param type For `model.matrix()`, `"fixed"` returns the fixed-effect design
#'   matrix and `"random"` returns the sparse random-effect design matrix. For
#'   `vcov()`, `"fixed"` returns the fixed-effect covariance surface and
#'   `"theta"` returns an unavailable theta-covariance matrix with a reason
#'   attribute.
#' @export
model.matrix.mm_lmm <- function(object, type = c("fixed", "random"), ...) {
  type <- match.arg(type)
  switch(
    type,
    fixed = .mm_lazy(object, "X", mm_fixed_model_matrix),
    random = .mm_lazy(object, "Z", mm_random_model_matrix)
  )
}

#' @rdname mm_lmm-methods
#' @export
model.matrix.mm_glmm <- model.matrix.mm_lmm

#' @rdname mm_lmm-methods
#' @export
vcov.mm_lmm <- function(object, type = c("fixed", "theta"),
                        correlation = FALSE, ...) {
  type <- match.arg(type)
  if (identical(type, "theta")) {
    theta_names <- parameterization(object)$table$theta_name
    if (length(theta_names) != length(object$theta)) {
      theta_names <- paste0("theta", seq_along(object$theta))
    }
    out <- matrix(NA_real_, nrow = length(object$theta), ncol = length(object$theta))
    dimnames(out) <- list(theta_names, theta_names)
    attr(out, "mm_unavailable_reason") <- "theta_covariance_unavailable"
    attr(out, "mm_method") <- "unavailable"
    return(out)
  }
  V <- object$fixed_effect_vcov %||%
    mm_fixed_effect_vcov_from_payload(
      object$artifact$fixed_effect_covariance_matrix,
      object$beta,
      object$std_errors,
      coef_map = object$coef_map
    )
  if (isTRUE(correlation)) {
    # Match lme4: attach the correlation matrix as a "correlation" attribute.
    attr(V, "correlation") <- stats::cov2cor(V)
  }
  V
}

#' @rdname mm_lmm-methods
#' @export
vcov.mm_glmm <- vcov.mm_lmm

mm_fixed_effect_vcov_from_payload <- function(payload, beta, std_errors,
                                              coef_map = NULL) {
  coef_names <- names(beta) %||% names(std_errors)
  if (!is.null(payload) &&
      identical(as.character(payload$schema_name %||% NA_character_),
                "mixedmodels.fixed_effect_covariance_matrix")) {
    status <- as.character(payload$status %||% "unavailable")
    if (!identical(status, "available") && !identical(status, "unavailable")) {
      mm_abort(
        message = sprintf("Fixed-effect covariance payload has unknown status `%s`.",
                          status),
        class = "mm_schema_error",
        input = payload
      )
    }
    matrix_payload <- payload$matrix %||% NULL
    if (identical(status, "available") && !is.null(matrix_payload)) {
      out <- mm_numeric_matrix_from_rows(matrix_payload)
      payload_names <- as.character(unlist(payload$coef_names %||% coef_names,
                                           use.names = FALSE))
      # Payloads produced by the engine label rows in the engine encoding;
      # when the caller's beta is already lme4-named (any normalized fit),
      # translate before validating/aligning.
      payload_names <- mm_coef_engine_to_lme4(payload_names, coef_map)
      mm_validate_fixed_effect_vcov_payload(
        out,
        payload_names,
        coef_names,
        payload
      )
      if (length(payload_names) == nrow(out) && length(payload_names) == ncol(out)) {
        dimnames(out) <- list(payload_names, payload_names)
      }
      out <- out[coef_names, coef_names, drop = FALSE]
      attr(out, "mm_method") <- as.character(payload$method %||% "model_based")
      attr(out, "mm_status") <- status
      attr(out, "mm_reliability") <- as.character(payload$reliability %||% NA_character_)
      attr(out, "mm_reason") <- payload$reason %||% NA_character_
      attr(out, "mm_details") <- payload$details %||% list()
      attr(out, "mm_notes") <- as.character(unlist(payload$notes %||% list(),
                                                    use.names = FALSE))
      attr(out, "mm_schema_name") <- as.character(payload$schema_name)
      attr(out, "mm_schema_version") <- as.character(payload$schema_version %||% NA_character_)
      return(out)
    }
    if (identical(status, "available")) {
      mm_abort(
        message = "Available fixed-effect covariance payload does not contain a matrix.",
        class = "mm_schema_error",
        input = payload
      )
    }
    if (!is.null(matrix_payload)) {
      mm_abort(
        message = "Unavailable fixed-effect covariance payload must not contain a matrix.",
        class = "mm_schema_error",
        input = payload
      )
    }
    reason <- as.character(payload$reason %||% NA_character_)
    if (length(reason) != 1L || is.na(reason) || !nzchar(reason)) {
      mm_abort(
        message = "Unavailable fixed-effect covariance payload must include a non-empty reason.",
        class = "mm_schema_error",
        input = payload
      )
    }
    out <- matrix(NA_real_, nrow = length(coef_names), ncol = length(coef_names),
                  dimnames = list(coef_names, coef_names))
    attr(out, "mm_method") <- as.character(payload$method %||% "model_based")
    attr(out, "mm_status") <- status
    attr(out, "mm_reliability") <- as.character(payload$reliability %||% "not_available")
    attr(out, "mm_unavailable_reason") <- reason
    attr(out, "mm_reason") <- attr(out, "mm_unavailable_reason")
    attr(out, "mm_details") <- payload$details %||% list()
    attr(out, "mm_notes") <- as.character(unlist(payload$notes %||% list(),
                                                  use.names = FALSE))
    attr(out, "mm_schema_name") <- as.character(payload$schema_name)
    attr(out, "mm_schema_version") <- as.character(payload$schema_version %||% NA_character_)
    return(out)
  }

  se <- as.numeric(std_errors)
  out <- diag(se^2, nrow = length(se), ncol = length(se))
  dimnames(out) <- list(names(std_errors), names(std_errors))
  attr(out, "mm_unavailable_reason") <- "fixed_effect_covariance_payload_unavailable"
  attr(out, "mm_method") <- "diagonal_from_stored_standard_errors"
  attr(out, "mm_status") <- "unavailable"
  attr(out, "mm_reliability") <- "not_available"
  out
}

mm_validate_fixed_effect_vcov_payload <- function(V, payload_names, coef_names, payload) {
  p <- length(coef_names)
  bad <- function(message) {
    mm_abort(message = message, class = "mm_schema_error", input = payload)
  }
  if (!is.matrix(V) || !is.numeric(V)) {
    bad("Available fixed-effect covariance payload must decode to a numeric matrix.")
  }
  if (!identical(dim(V), c(p, p))) {
    bad(sprintf(
      "Available fixed-effect covariance payload has shape %s, expected %s x %s.",
      paste(dim(V), collapse = " x "), p, p
    ))
  }
  if (length(payload_names) != p ||
      any(is.na(payload_names)) ||
      any(!nzchar(payload_names))) {
    bad("Available fixed-effect covariance payload must include one coefficient name per row and column.")
  }
  if (!identical(sort(payload_names), sort(coef_names))) {
    bad("Available fixed-effect covariance payload coefficient names do not match the fitted fixed effects.")
  }
  if (any(!is.finite(V))) {
    bad("Available fixed-effect covariance payload must be finite.")
  }
  if (!isSymmetric(unname(V), tol = sqrt(.Machine$double.eps))) {
    bad("Available fixed-effect covariance payload must be symmetric.")
  }
  invisible(TRUE)
}

mm_numeric_matrix_from_rows <- function(rows) {
  if (is.matrix(rows)) {
    return(matrix(as.numeric(rows), nrow = nrow(rows), ncol = ncol(rows)))
  }
  matrix(
    as.numeric(unlist(rows, use.names = FALSE)),
    nrow = length(rows),
    byrow = TRUE
  )
}

#' Inspect random-effect blocks
#'
#' `random_blocks()` summarizes the random-effect block structure recorded in
#' the compiler artifact: grouping factor, basis, covariance family, theta
#' parameter count, level counts, and design-support status.
#'
#' @param object A compiled `mm_spec` or fitted `mm_fit`.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_random_blocks` object with a data-frame `table`.
#'
#' @export
random_blocks <- function(object, ...) {
  UseMethod("random_blocks")
}

#' @rdname random_blocks
#' @export
random_blocks.mm_compiled <- function(object, ...) {
  artifact <- mm_compiled_artifact(object)
  rows <- lapply(artifact$design_audit$random_terms %||% list(), mm_random_block_row)
  table <- if (length(rows)) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  } else {
    mm_random_blocks_empty_table()
  }
  obj <- list(
    table = table,
    random_terms = artifact$design_audit$random_terms %||% list(),
    semantic_terms = artifact$semantic_model$random_terms %||% list()
  )
  class(obj) <- "mm_random_blocks"
  obj
}

#' @export
random_blocks.default <- function(object, ...) {
  mm_abort(
    message = "`random_blocks()` expects a compiled or fitted mixeff object.",
    class = "mm_arg_error",
    input = object
  )
}

#' @method print mm_random_blocks
#' @export
print.mm_random_blocks <- function(x, ...) {
  cat("Random-effect blocks:\n")
  if (!nrow(x$table)) {
    cat("  none\n")
    return(invisible(x))
  }
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Inspect the optimizer certificate
#'
#' @param object A compiled `mm_spec` or fitted `mm_fit`.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_optimizer_certificate` object containing the raw certificate
#'   and a compact table view.
#'
#' @export
optimizer_certificate <- function(object, ...) {
  UseMethod("optimizer_certificate")
}

#' @rdname optimizer_certificate
#' @export
optimizer_certificate.mm_compiled <- function(object, ...) {
  cert <- mm_compiled_artifact(object)$optimizer_certificate %||% list()
  obj <- list(raw = cert, table = mm_optimizer_certificate_table(cert))
  class(obj) <- "mm_optimizer_certificate"
  obj
}

#' @method print mm_optimizer_certificate
#' @export
print.mm_optimizer_certificate <- function(x, ...) {
  cat("Optimizer certificate:\n")
  if (!nrow(x$table)) {
    cat("  not assessed\n")
    return(invisible(x))
  }
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Fixed-effect inference table
#'
#' Fitted artifacts may carry Rust-owned fixed-effect inference rows. When
#' present, those rows are the source of truth for estimates, standard errors,
#' degrees of freedom, statistics, p-values, methods, status, reliability, and
#' unavailable reasons. Legacy objects without this artifact field fall back to
#' an unavailable table.
#'
#' @param fit A fitted `mm_lmm`.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_inference_table`.
#'
#' @export
inference_table <- function(fit, ...) {
  UseMethod("inference_table")
}

#' @rdname inference_table
#' @param method Inference method. `"auto"` (the default) returns the
#'   artifact-cached table that the engine resolved at fit time. Any other
#'   value (`"satterthwaite"`, `"kenward_roger"`, `"asymptotic"`, `"none"`)
#'   recomputes the table by dispatching one `contrast()` per fixed-effect
#'   term with the requested method, so refusals and reasons are surfaced
#'   honestly rather than silently swapped for the auto-resolved row.
#' @export
inference_table.mm_lmm <- function(fit,
                                   method = c("auto", "satterthwaite",
                                              "kenward_roger", "asymptotic",
                                              "none"),
                                   ...) {
  method <- match.arg(method)
  if (!identical(method, "auto")) {
    return(mm_inference_table_recompute(fit, method))
  }
  parsed <- mm_json_parse_fixed_effect_inference_table(
    fit$artifact$fixed_effect_inference_table %||% NULL
  )
  if (!is.null(parsed)) {
    tbl <- parsed$table
    if (!is.null(tbl$label)) {
      # Artifact rows carry engine-encoded labels in ENGINE column order;
      # translate to lme4 names and reorder coefficient rows to match
      # names(fit$beta), so positional pairing with fixef()/vcov() is safe.
      tbl$label <- mm_coef_engine_to_lme4(tbl$label, fit$coef_map)
      if (!is.null(tbl$kind) && !is.null(fit$coef_map)) {
        is_coef <- tbl$kind == "coefficient"
        coef_rows <- tbl[is_coef, , drop = FALSE]
        hit <- match(names(fit$beta), coef_rows$label)
        if (!anyNA(hit) && nrow(coef_rows) == length(hit)) {
          tbl <- rbind(coef_rows[hit, , drop = FALSE],
                       tbl[!is_coef, , drop = FALSE])
          rownames(tbl) <- NULL
        }
      }
    }
    obj <- list(table = tbl, raw = parsed$raw)
    class(obj) <- "mm_inference_table"
    return(obj)
  }

  table <- data.frame(
    term = names(fit$beta),
    label = names(fit$beta),
    kind = "coefficient",
    estimate = unname(fit$beta),
    std_error = unname(fit$std_errors),
    df = NA_real_,
    numerator_df = NA_real_,
    denominator_df = NA_real_,
    statistic = NA_real_,
    statistic_name = NA_character_,
    p_value = NA_real_,
    method = "not_computed",
    status = "not_assessed",
    reliability = "not_available",
    reliability_reason = NA_character_,
    estimability = "not_assessed",
    reason = "fixed_effect_inference_table_unavailable_legacy_object",
    notes = I(rep(list(character()), length(fit$beta))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  obj <- list(table = table, raw = NULL)
  class(obj) <- "mm_inference_table"
  obj
}

#' @method print mm_inference_table
#' @export
print.mm_inference_table <- function(x, ...) {
  cat("Inference table:\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Inspect reproducibility metadata
#'
#' @param object A compiled `mm_spec` or fitted `mm_fit`.
#' @param ... Reserved for future methods.
#'
#' @return An `mm_reproducibility` object.
#'
#' @export
reproducibility <- function(object, ...) {
  UseMethod("reproducibility")
}

#' @rdname reproducibility
#' @export
reproducibility.mm_compiled <- function(object, ...) {
  raw <- mm_compiled_artifact(object)$reproducibility %||% list()
  obj <- list(raw = raw, thresholds = mm_repro_threshold_table(raw$thresholds))
  class(obj) <- "mm_reproducibility"
  obj
}

#' @method print mm_reproducibility
#' @export
print.mm_reproducibility <- function(x, ...) {
  cat("Reproducibility:\n")
  cat(sprintf("  fit intent: %s\n", mm_scalar_text(x$raw$fit_intent, "not_recorded")))
  cat(sprintf("  random state used: %s\n",
              mm_scalar_text(x$raw$random_state_used, "not_recorded")))
  if (nrow(x$thresholds)) {
    cat("Thresholds:\n")
    print(x$thresholds, row.names = FALSE)
  }
  invisible(x)
}

#' Test whether a fit is singular or reduced-rank
#'
#' @param x A fitted `mm_lmm`.
#' @param tol Reserved for compatibility with lme4's `isSingular()`.
#' @param ... Reserved for future methods.
#'
#' @return A length-one logical value.
#'
#' @export
is_singular <- function(x, tol = 1e-4, ...) {
  UseMethod("is_singular")
}

#' @rdname is_singular
#' @export
is_singular.mm_lmm <- function(x, tol = 1e-4, ...) {
  status <- tolower(fit_status(x))
  cov_status <- vapply(
    x$artifact$effective_covariance %||% list(),
    function(one) tolower(mm_scalar_text(one$status)),
    character(1)
  )
  status %in% c("converged_boundary", "converged_reduced_rank",
                "boundary", "reduced_rank") ||
    any(cov_status %in% c("boundary", "reduced_rank", "singular"))
}

#' @export
is_singular.default <- function(x, tol = 1e-4, ...) {
  mm_abort(
    message = "`is_singular()` expects a fitted mixeff LMM.",
    class = "mm_arg_error",
    input = x
  )
}

# ---- lazy-cache helpers ----------------------------------------------------

mm_empty_lazy_cache <- function() {
  new.env(parent = emptyenv())
}

mm_object_schema <- function(artifact) {
  header <- artifact$schema %||% list()
  list(
    schema_name = as.character(header$schema_name %||% NA_character_),
    schema_version = as.character(header$schema_version %||% NA_character_),
    crate_version = as.character(header$crate_version %||% NA_character_),
    package_version = tryCatch(
      as.character(utils::packageVersion("mixeff")),
      error = function(cnd) NA_character_
    )
  )
}

.mm_lazy <- function(fit, key, producer) {
  if (!inherits(fit, "mm_fit")) {
    mm_abort(
      message = "Lazy extractors require a fitted mixeff object.",
      class = "mm_arg_error",
      input = fit
    )
  }
  if (!is.environment(fit$lazy_cache)) {
    fit <- revive(fit)
  }
  cache <- fit$lazy_cache
  if (!exists(key, envir = cache, inherits = FALSE)) {
    assign(key, producer(fit), envir = cache)
  }
  get(key, envir = cache, inherits = FALSE)
}

mm_fixed_model_matrix <- function(fit) {
  stats::model.matrix(mm_fixed_formula(fit), data = fit$model_frame)
}

mm_fixed_formula <- function(fit) {
  response <- mm_response_name(fit)
  fixed <- as.character(unlist(
    fit$artifact$semantic_model$fixed_terms %||% list("1"),
    use.names = FALSE
  ))
  has_intercept <- "1" %in% fixed
  fixed <- setdiff(fixed, "1")
  rhs <- if (length(fixed)) {
    paste(fixed, collapse = " + ")
  } else {
    "1"
  }
  if (!has_intercept && !identical(rhs, "1")) {
    rhs <- paste("0", rhs, sep = " + ")
  }
  stats::as.formula(paste(response, "~", rhs), env = environment(fit$formula))
}

mm_response_name <- function(fit) {
  response <- fit$artifact$semantic_model$response %||% NULL
  if (is.character(response) && length(response) == 1L && nzchar(response)) {
    return(response)
  }
  all.vars(fit$formula)[[1L]]
}

mm_response_vector <- function(fit) {
  out <- fit$model_frame[[mm_response_name(fit)]]
  names(out) <- rownames(fit$model_frame)
  out
}

mm_random_model_matrix <- function(fit) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  if (!length(terms)) {
    return(Matrix::Matrix(numeric(), nrow = nrow(fit$model_frame), ncol = 0L,
                          sparse = TRUE))
  }
  mats <- lapply(seq_along(terms), function(i) {
    mm_random_term_matrix(fit, terms[[i]], i)
  })
  out <- do.call(cbind, mats)
  Matrix::Matrix(out, sparse = TRUE)
}

mm_random_term_matrix <- function(fit, term, index) {
  group_label <- mm_random_term_group_label(fit, term, index)
  group <- mm_group_factor(fit$model_frame, group_label)
  levels <- levels(group)
  basis <- term$basis %||% list()
  basis_labels <- vapply(basis, mm_basis_label, character(1))
  basis_values <- lapply(basis, mm_basis_values, frame = fit$model_frame)
  if (!length(basis_values)) {
    basis_labels <- "(Intercept)"
    basis_values <- list(rep(1, nrow(fit$model_frame)))
  }

  cols <- vector("list", length(levels) * length(basis_values))
  col_names <- character(length(cols))
  k <- 0L
  for (level in levels) {
    group_mask <- as.numeric(group == level)
    for (j in seq_along(basis_values)) {
      k <- k + 1L
      cols[[k]] <- group_mask * basis_values[[j]]
      col_names[[k]] <- paste(group_label, level, basis_labels[[j]], sep = ".")
    }
  }
  out <- do.call(cbind, cols)
  colnames(out) <- col_names
  rownames(out) <- rownames(fit$model_frame)
  out
}

mm_random_term_group_label <- function(fit, term, index) {
  audit <- fit$artifact$design_audit$random_terms %||% list()
  audit_group <- audit[[index]]$group$name %||% NULL
  if (is.character(audit_group) && length(audit_group) == 1L && nzchar(audit_group)) {
    return(audit_group)
  }
  single <- term$group$single$name %||% NULL
  if (is.character(single) && length(single) == 1L && nzchar(single)) {
    return(single)
  }
  mm_scalar_text(term$group, sprintf("term_%d", index))
}

mm_group_factor <- function(frame, group_label) {
  if (group_label %in% names(frame)) {
    return(factor(frame[[group_label]]))
  }
  parts <- strsplit(group_label, ":", fixed = TRUE)[[1L]]
  if (length(parts) > 1L && all(parts %in% names(frame))) {
    return(interaction(frame[, parts, drop = FALSE], drop = TRUE, sep = ":"))
  }
  mm_abort(
    message = sprintf("Cannot rebuild grouping factor `%s` from the stored model frame.",
                      group_label),
    class = "mm_data_error",
    input = group_label
  )
}

mm_basis_label <- function(basis) {
  name <- mm_scalar_text(basis$name %||% basis$source)
  kind <- mm_scalar_text(basis$kind)
  if (identical(kind, "intercept") || identical(name, "intercept") ||
      identical(name, "1")) {
    "(Intercept)"
  } else {
    name
  }
}

mm_basis_values <- function(basis, frame) {
  label <- mm_basis_label(basis)
  if (identical(label, "(Intercept)")) {
    return(rep(1, nrow(frame)))
  }
  source <- mm_scalar_text(basis$source, label)
  candidates <- unique(c(label, source))
  hits <- candidates[candidates %in% names(frame)]
  hit <- if (length(hits)) hits[[1L]] else NA_character_
  if (is.na(hit)) {
    mm_abort(
      message = sprintf("Cannot rebuild random-effect basis `%s` from the stored model frame.",
                        label),
      class = "mm_data_error",
      input = label
    )
  }
  value <- frame[[hit]]
  if (!is.numeric(value) && !is.integer(value)) {
    mm_abort(
      message = sprintf("Random-effect basis `%s` is not numeric; lazy R-side rebuild is unavailable.",
                        label),
      class = "mm_inference_unavailable",
      input = label
    )
  }
  as.numeric(value)
}

mm_random_flist <- function(fit) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  out <- lapply(seq_along(terms), function(i) {
    label <- mm_random_term_group_label(fit, terms[[i]], i)
    mm_group_factor(fit$model_frame, label)
  })
  names(out) <- vapply(seq_along(terms), function(i) {
    mm_random_term_group_label(fit, terms[[i]], i)
  }, character(1))
  class(out) <- c("mm_flist", "list")
  out
}

mm_random_cnms <- function(fit) {
  re <- ranef(fit)
  out <- lapply(re, names)
  class(out) <- c("mm_cnms", "list")
  out
}

mm_lambda_matrix <- function(fit) {
  terms <- fit$artifact$semantic_model$random_terms %||% list()
  if (!length(terms)) {
    return(Matrix::Matrix(numeric(), nrow = 0L, ncol = 0L, sparse = TRUE))
  }
  theta <- as.numeric(fit$theta)
  offset <- 0L
  blocks <- vector("list", length(terms))
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    group <- mm_group_factor(fit$model_frame, mm_random_term_group_label(fit, term, i))
    p <- length(term$basis %||% list())
    if (!p) p <- 1L
    family <- mm_scalar_text(term$covariance, "full")
    n_theta <- switch(
      family,
      full = p * (p + 1L) / 2L,
      diagonal = p,
      diag = p,
      scalar = 1L,
      p * (p + 1L) / 2L
    )
    available <- max(0L, length(theta) - offset)
    piece <- theta[seq.int(offset + 1L, length.out = min(n_theta, available))]
    offset <- offset + n_theta
    if (!length(piece)) piece <- 1
    L <- matrix(0, nrow = p, ncol = p)
    if (identical(family, "scalar")) {
      diag(L) <- piece[[1L]]
    } else if (family %in% c("diagonal", "diag")) {
      diag(L) <- rep(piece, length.out = p)
    } else {
      L[lower.tri(L, diag = TRUE)] <- rep(piece, length.out = p * (p + 1L) / 2L)
    }
    blocks[[i]] <- kronecker(
      Matrix::Diagonal(n = length(levels(group))),
      Matrix::Matrix(L, sparse = TRUE)
    )
  }
  out <- do.call(Matrix::bdiag, blocks)
  attr(out, "mm_method") <- "rebuilt_from_stored_theta"
  out
}

mm_random_block_row <- function(term) {
  group <- term$group %||% list()
  budget <- term$information_budget %||% list()
  basis <- term$basis %||% list()
  data.frame(
    term_id = mm_scalar_text(term$term_id),
    group = mm_scalar_text(group$name),
    basis = paste(vapply(basis, mm_basis_label, character(1)), collapse = ", "),
    covariance = mm_scalar_text(budget$covariance_family),
    theta_parameters = as.integer(term$requested_covariance_parameters %||% NA_integer_),
    group_levels = as.integer(group$n_levels %||% NA_integer_),
    min_rows_per_group = as.integer(group$min_obs_per_level %||% NA_integer_),
    median_rows_per_group = as.numeric(group$median_obs_per_level %||% NA_real_),
    status = mm_scalar_text(budget$status, "not_assessed"),
    reason = mm_scalar_text(budget$reason),
    stringsAsFactors = FALSE
  )
}

mm_random_blocks_empty_table <- function() {
  data.frame(
    term_id = character(),
    group = character(),
    basis = character(),
    covariance = character(),
    theta_parameters = integer(),
    group_levels = integer(),
    min_rows_per_group = integer(),
    median_rows_per_group = numeric(),
    status = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

mm_optimizer_certificate_table <- function(cert) {
  if (!length(cert)) {
    return(data.frame(metric = character(), value = character(),
                      stringsAsFactors = FALSE))
  }
  fields <- c(
    status = "status",
    optimizer = "optimizer_name",
    objective = "objective_value",
    iterations = "iterations",
    free_gradient_norm = "free_gradient_norm",
    projected_gradient_norm = "projected_gradient_norm",
    hessian_eigen_min = "hessian_eigen_min",
    hessian_rank = "hessian_rank",
    information_rank = "information_rank"
  )
  rows <- lapply(names(fields), function(metric) {
    value <- cert[[fields[[metric]]]]
    data.frame(
      metric = metric,
      value = mm_scalar_text(value, "not_recorded"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mm_repro_threshold_table <- function(thresholds) {
  thresholds <- thresholds %||% list()
  if (!length(thresholds)) {
    return(data.frame(name = character(), value = character(),
                      stringsAsFactors = FALSE))
  }
  rows <- lapply(thresholds, function(item) {
    flat <- as.character(unlist(item, use.names = FALSE))
    data.frame(
      name = if (length(flat) >= 1L) flat[[1L]] else "",
      value = if (length(flat) >= 2L) flat[[2L]] else "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Recompute an inference table by dispatching one contrast() per fixed-effect
# term with an explicit method. Used by inference_table(method != "auto") and,
# transitively, by summary(method != "auto") so the user's method request is
# honored instead of silently replaced by the auto-resolved cached row.
mm_inference_table_recompute <- function(fit, method) {
  terms <- names(fit$beta)
  k <- length(terms)
  if (!k) {
    empty <- mm_fixed_effect_inference_empty_table()
    obj <- list(table = empty, raw = NULL,
                requested_method = method, recomputed = TRUE)
    class(obj) <- "mm_inference_table"
    return(obj)
  }
  rows <- lapply(seq_len(k), function(i) {
    L <- rep(0, k)
    L[i] <- 1
    names(L) <- terms
    ct <- contrast(fit, L, method = method)
    mm_inference_row_from_contrast(ct$table, term = terms[[i]])
  })
  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL
  obj <- list(table = tbl, raw = NULL,
              requested_method = method, recomputed = TRUE)
  class(obj) <- "mm_inference_table"
  obj
}

# Reshape a single-row contrast table into the inference-table row shape,
# preserving the engine-authored method, status, reliability_reason, reason,
# reason_code, reason_detail, estimability, and notes fields.
mm_inference_row_from_contrast <- function(ct, term) {
  pull <- function(col, default) {
    if (col %in% names(ct)) ct[[col]][[1L]] else default
  }
  pull_list <- function(col) {
    if (col %in% names(ct)) ct[[col]][1L] else list(NULL)
  }
  data.frame(
    term = term,
    label = term,
    kind = "coefficient",
    estimate = pull("estimate", NA_real_),
    std_error = pull("std_error", NA_real_),
    df = pull("df", NA_real_),
    numerator_df = NA_real_,
    denominator_df = pull("df", NA_real_),
    statistic = pull("statistic", NA_real_),
    statistic_name = pull("statistic_name", NA_character_),
    p_value = pull("p_value", NA_real_),
    method = pull("method", "not_computed"),
    status = pull("status", "not_assessed"),
    reliability = pull("reliability", "not_available"),
    reliability_reason = pull("reliability_reason", NA_character_),
    reason = pull("reason", NA_character_),
    reason_code = pull("reason_code", NA_character_),
    reason_detail = pull("reason_detail", NA_character_),
    estimability = I(pull_list("estimability")),
    details = I(pull_list("details")),
    notes = I(pull_list("notes")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
