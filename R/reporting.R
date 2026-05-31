#' Produce reporting tables for a fitted mixeff model
#'
#' `model_report()` assembles a structured, publication-oriented report from
#' the Rust artifact fields carried by a fitted model plus R-owned provenance
#' such as the call and session metadata. `reporting_table()` extracts one
#' section as a data-frame-compatible object.
#'
#' @param fit A fitted `mm_fit`, usually from [lmm()].
#' @param sections Character vector of report sections, or `"all"`.
#' @param section One section name, or `"all"`.
#' @param view `"compact"` for reader-facing columns, or `"audit"` for the
#'   full provenance table with `source`, `reason`, `details`, and related
#'   audit columns.
#' @param ... Reserved for future methods.
#'
#' @return `model_report()` returns an `mm_model_report`. `reporting_table()`
#'   returns a data frame for one section or a named list for `"all"`.
#'
#' @export
model_report <- function(fit, sections = "all", ...) {
  UseMethod("model_report")
}

#' @rdname model_report
#' @export
model_report.mm_fit <- function(fit, sections = "all", ...) {
  sections <- mm_report_sections_arg(sections)
  builders <- mm_report_builders()
  built <- lapply(builders, function(builder) builder(fit))
  built <- built[sections]

  unavailable <- do.call(rbind, lapply(built, function(section) {
    section$unavailable
  }))
  if (is.null(unavailable) || !nrow(unavailable)) {
    unavailable <- mm_report_unavailable_empty()
  }

  report <- list(
    metadata = mm_report_metadata(fit, names(built)),
    sections = lapply(built, `[[`, "table"),
    unavailable = unavailable,
    provenance = do.call(rbind, lapply(names(built), function(name) {
      data.frame(
        section = name,
        source = built[[name]]$source,
        stringsAsFactors = FALSE
      )
    }))
  )
  if ("unavailable" %in% names(report$sections)) {
    report$sections$unavailable <- unavailable
  }
  class(report) <- c("mm_model_report", "mm_report")
  report
}

#' @rdname model_report
#' @export
reporting_table <- function(fit, section = "all", view = c("compact", "audit"), ...) {
  UseMethod("reporting_table")
}

#' @rdname model_report
#' @export
reporting_table.mm_fit <- function(fit, section = "all",
                                   view = c("compact", "audit"), ...) {
  view <- match.arg(view)
  section <- mm_report_sections_arg(section, allow_many = FALSE)
  report <- model_report(fit, sections = if (identical(section, "all")) "all" else section)
  if (identical(section, "all")) {
    return(mm_report_view_sections(report$sections, view))
  }
  mm_report_view_table(report$sections[[section]], section, view)
}

#' @method reporting_table mm_model_report
#' @export
reporting_table.mm_model_report <- function(fit, section = "all",
                                            view = c("compact", "audit"), ...) {
  view <- match.arg(view)
  section <- mm_report_sections_arg(section, allow_many = FALSE)
  if (identical(section, "all")) {
    return(mm_report_view_sections(fit$sections, view))
  }
  mm_report_view_table(fit$sections[[section]], section, view)
}

#' @rdname model_report
#' @method reporting_table mm_model_comparison
#' @export
reporting_table.mm_model_comparison <- function(fit, section = "comparison_ledger",
                                                view = c("compact", "audit"),
                                                ...) {
  mm_report_comparison_object_table(fit, section = section, view = view)
}

#' @rdname model_report
#' @method reporting_table mm_drop1
#' @export
reporting_table.mm_drop1 <- function(fit, section = "comparison_ledger",
                                     view = c("compact", "audit"), ...) {
  mm_report_comparison_object_table(fit, section = section, view = view)
}

mm_report_comparison_object_table <- function(fit, section = "comparison_ledger",
                                              view = c("compact", "audit")) {
  view <- match.arg(view)
  section <- mm_report_sections_arg(section, allow_many = FALSE)
  table <- mm_report_section(fit$ledger %||% mm_comparison_ledger_empty(),
                             "mixeff.comparison")$table
  if (identical(section, "all")) {
    out <- list(comparison_ledger = mm_report_view_table(table, "comparison_ledger", view))
    return(out)
  }
  if (!identical(section, "comparison_ledger")) {
    mm_abort(
      message = "Comparison objects only expose the `comparison_ledger` reporting section.",
      class = "mm_schema_error",
      input = section
    )
  }
  mm_report_view_table(table, "comparison_ledger", view)
}

#' @rdname model_report
#' @method reporting_table mm_random_effect_test
#' @export
reporting_table.mm_random_effect_test <- function(fit, section = "all",
                                                  view = c("compact", "audit"),
                                                  ...) {
  view <- match.arg(view)
  table <- fit$table
  if (identical(view, "audit")) {
    return(table)
  }
  keep <- intersect(
    c("term", "group", "statistic", "statistic_name", "p_value",
      "reference_distribution", "method", "status", "reason_code"),
    names(table)
  )
  table[, keep, drop = FALSE]
}

#' @method print mm_model_report
#' @export
print.mm_model_report <- function(x, ...) {
  cat("mixeff model report\n")
  overview <- x$sections$overview
  if (!is.null(overview) && nrow(overview)) {
    show <- overview[overview$field %in% c("formula", "fit_method", "nobs",
                                           "fit_status", "inference"), ,
                     drop = FALSE]
    if (nrow(show)) {
      print(show[, c("field", "value", "status")], row.names = FALSE)
    }
  }
  cat("\nSections:\n")
  cat(paste0("  ", names(x$sections), collapse = "\n"), "\n", sep = "")
  if (nrow(x$unavailable)) {
    cat(sprintf("\nUnavailable/caveated fields: %d\n", nrow(x$unavailable)))
  }
  invisible(x)
}

mm_report_section_names <- function() {
  c(
    "overview",
    "model_specification",
    "data_design",
    "random_terms",
    "random_effects",
    "fixed_effects",
    "fit_statistics",
    "optimizer",
    "comparison_ledger",
    "reproducibility",
    "unavailable"
  )
}

mm_report_sections_arg <- function(sections, allow_many = TRUE) {
  known <- mm_report_section_names()
  if (!is.character(sections) || !length(sections) || anyNA(sections)) {
    mm_abort(
      message = "`sections` must be a non-empty character vector.",
      class = "mm_schema_error",
      input = sections
    )
  }
  if ("all" %in% sections) {
    if (!allow_many && length(sections) != 1L) {
      mm_abort(
        message = "`section = \"all\"` cannot be combined with other section names.",
        class = "mm_schema_error",
        input = sections
      )
    }
    return(if (allow_many) known else "all")
  }
  unknown <- setdiff(sections, known)
  if (length(unknown)) {
    mm_abort(
      message = sprintf("Unknown report section(s): %s.",
                        paste(unknown, collapse = ", ")),
      class = "mm_schema_error",
      input = sections
    )
  }
  if (!allow_many && length(sections) != 1L) {
    mm_abort(
      message = "`section` must be a single section name.",
      class = "mm_schema_error",
      input = sections
    )
  }
  sections
}

mm_report_view_sections <- function(sections, view) {
  out <- lapply(names(sections), function(section) {
    mm_report_view_table(sections[[section]], section, view)
  })
  names(out) <- names(sections)
  out
}

mm_report_view_table <- function(table, section, view = c("compact", "audit")) {
  view <- match.arg(view)
  if (identical(view, "audit")) {
    return(table)
  }
  cols <- switch(
    section,
    overview = c("field", "value"),
    fixed_effects = c("term", "estimate", "std_error", "df", "statistic",
                      "statistic_name", "p_value", "method", "status",
                      "reliability"),
    data_design = c("group", "role", "group_levels", "min_rows_per_group",
                    "median_rows_per_group", "max_rows_per_group", "status"),
    random_terms = c("term_id", "original_fragment", "group", "basis",
                     "covariance", "theta_parameters", "design_status",
                     "english"),
    random_effects = c("group", "basis_lhs", "kind", "variance", "std_dev",
                       "correlation", "status"),
    fit_statistics = c("field", "value"),
    optimizer = c("metric", "value", "status"),
    reproducibility = c("field", "value"),
    unavailable = c("section", "field", "status", "reason"),
    model_specification = c("field", "value", "status"),
    comparison_ledger = c("comparison_id", "formula", "comparison_method",
                          "statistic", "p_value", "status", "reason"),
    names(table)
  )
  cols <- intersect(cols, names(table))
  out <- table[, cols, drop = FALSE]
  out <- mm_drop_empty_report_columns(out)
  attr(out, "view") <- "compact"
  attr(out, "audit_columns") <- setdiff(names(table), names(out))
  out
}

mm_drop_empty_report_columns <- function(table) {
  keep <- vapply(table, function(x) {
    if (is.list(x)) {
      return(any(lengths(x) > 0L))
    }
    if (is.numeric(x) || is.integer(x) || is.logical(x)) {
      return(any(!is.na(x)))
    }
    values <- as.character(x)
    any(!is.na(values) & nzchar(values) & values != "reason_not_recorded")
  }, logical(1))
  table[, keep, drop = FALSE]
}

mm_report_builders <- function() {
  list(
    overview = mm_report_overview,
    model_specification = mm_report_model_specification,
    data_design = mm_report_data_design,
    random_terms = mm_report_random_terms,
    random_effects = mm_report_random_effects,
    fixed_effects = mm_report_fixed_effects,
    fit_statistics = mm_report_fit_statistics,
    optimizer = mm_report_optimizer,
    comparison_ledger = mm_report_comparison_ledger,
    reproducibility = mm_report_reproducibility,
    unavailable = mm_report_unavailable_section
  )
}

mm_report_metadata <- function(fit, sections) {
  schema <- fit$schema %||% mm_object_schema(fit$artifact)
  data.frame(
    field = c("created_at", "package_version", "crate_version",
              "artifact_schema", "artifact_schema_version", "fit_class",
              "sections"),
    value = c(
      format(Sys.time(), usetz = TRUE),
      mm_report_package_version(),
      mm_scalar_text(schema$crate_version, "not_recorded"),
      mm_scalar_text(schema$schema_name, "not_recorded"),
      mm_scalar_text(schema$schema_version, "not_recorded"),
      paste(class(fit), collapse = ", "),
      paste(sections, collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
}

mm_report_overview <- function(fit) {
  artifact <- mm_compiled_artifact(fit)
  inference <- inference_table(fit)$table
  available <- if (nrow(inference)) sum(inference$status == "available") else 0L
  total <- nrow(inference)
  schema <- fit$schema %||% mm_object_schema(artifact)
  table <- mm_report_kv(
    field = c("model_class", "formula", "effective_formula", "fit_method",
              "mode", "nobs", "fit_status", "inference",
              "artifact_schema", "crate_version", "package_version"),
    value = c(
      if (inherits(fit, "mm_lmm")) "LMM" else class(fit)[[1L]],
      deparse1(fit$formula),
      artifact$effective_formula %||% artifact$requested_formula %||% deparse1(fit$formula),
      if (isTRUE(fit$REML)) "REML" else "ML",
      mm_scalar_text(artifact$reproducibility$fit_intent, "not_recorded"),
      as.character(fit$nobs %||% nrow(fit$model_frame)),
      fit_status(fit),
      sprintf("%d/%d available fixed-effect rows", available, total),
      paste(mm_scalar_text(schema$schema_name, "not_recorded"),
            mm_scalar_text(schema$schema_version, "not_recorded")),
      mm_scalar_text(schema$crate_version, "not_recorded"),
      mm_report_package_version()
    ),
    source = c(rep("mixed", 2L), "CompiledModelArtifact", "mm_fit",
               "reproducibility", "mm_fit", "OptimizerCertificate",
               "fixed_effect_inference_table", "CompiledModelArtifact",
               "CompiledModelArtifact", "R package metadata")
  )
  mm_report_section(table, "mixed")
}

mm_report_model_specification <- function(fit) {
  artifact <- mm_compiled_artifact(fit)
  changes_table <- changes(fit)$table
  base <- mm_report_kv(
    field = c("call", "formula", "requested_formula", "effective_formula",
              "fixed_terms", "random_terms"),
    value = c(
      paste(deparse(fit$call), collapse = " "),
      deparse1(fit$formula),
      artifact$requested_formula %||% "",
      artifact$effective_formula %||% artifact$requested_formula %||% "",
      mm_list_text(artifact$semantic_model$fixed_terms),
      paste(vapply(artifact$semantic_model$random_terms %||% list(), function(term) {
        term$source_syntax$text %||% term$source_syntax$written %||% ""
      }, character(1)), collapse = " + ")
    ),
    source = c("R call", "R formula", "CompiledModelArtifact",
               "CompiledModelArtifact", "semantic_model", "semantic_model")
  )
  if (nrow(changes_table)) {
    changes_rows <- data.frame(
      field = paste0("change:", seq_len(nrow(changes_table))),
      value = paste(changes_table$stage, changes_table$status,
                    changes_table$detail, sep = " | "),
      source = "changes",
      status = changes_table$status,
      reason = changes_table$detail,
      stringsAsFactors = FALSE
    )
    base <- rbind(base, changes_rows)
  }
  mm_report_section(base, "mixed")
}

mm_report_data_design <- function(fit) {
  blocks <- random_blocks(fit)$table
  if (!nrow(blocks)) {
    table <- data.frame(
      group = character(),
      role = character(),
      group_levels = integer(),
      min_rows_per_group = integer(),
      median_rows_per_group = numeric(),
      max_rows_per_group = integer(),
      status = character(),
      reason = character(),
      source = character(),
      stringsAsFactors = FALSE
    )
    return(mm_report_section(table, "ModelAuditReport"))
  }

  blocks$role <- mm_report_group_roles(fit, blocks$group)
  blocks$max_rows_per_group <- vapply(blocks$group, function(group) {
    mm_report_group_max_rows(fit, group)
  }, integer(1))
  table <- blocks[, c("group", "role", "group_levels", "min_rows_per_group",
                      "median_rows_per_group", "max_rows_per_group", "status",
                      "reason"), drop = FALSE]
  table$source <- "ModelAuditReport/random_term_cards"
  mm_report_section(table, "ModelAuditReport/random_term_cards")
}

mm_report_random_terms <- function(fit) {
  audit <- audit_design(fit)
  cards <- audit$random_term_cards %||% list()
  rows <- unlist(lapply(cards, mm_report_random_term_card_rows), recursive = FALSE)
  table <- if (length(rows)) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(
      term_id = character(),
      original_fragment = character(),
      canonical_fragment = character(),
      group = character(),
      block = integer(),
      basis = character(),
      intercept = logical(),
      slopes = character(),
      covariance = character(),
      theta_parameters = integer(),
      english = character(),
      constraints = character(),
      design_status = character(),
      source = character(),
      stringsAsFactors = FALSE
    )
  }
  if (length(audit$cross_card_constraints %||% list())) {
    table <- rbind(table, mm_report_cross_card_rows(audit$cross_card_constraints))
  }
  mm_report_section(table, "ModelAuditReport.random_term_cards")
}

mm_report_random_effects <- function(fit) {
  has_fit_summary <- is.list(fit$fit_summary) && is.list(fit$fit_summary$varcorr)
  vc_obj <- if (has_fit_summary) {
    mm_varcorr_from_result(fit$fit_summary$varcorr)
  } else {
    VarCorr(fit)
  }
  varcorr_source <- if (has_fit_summary) {
    "mixedmodels.fit_summary.varcorr"
  } else {
    "fit$varcorr"
  }
  varcorr_status <- if (has_fit_summary) "available" else "available_from_varcorr"
  vc <- vc_obj$table
  table <- if (nrow(vc)) {
    data.frame(
      group = vc$group,
      term_id = NA_character_,
      basis_lhs = vc$name,
      basis_rhs = vc$name,
      kind = "variance",
      variance = vc$variance,
      std_dev = vc$std_dev,
      correlation = vc$correlation,
      covariance_family = NA_character_,
      status = varcorr_status,
      reason = NA_character_,
      source = varcorr_source,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      group = character(),
      term_id = character(),
      basis_lhs = character(),
      basis_rhs = character(),
      kind = character(),
      variance = numeric(),
      std_dev = numeric(),
      correlation = character(),
      covariance_family = character(),
      status = character(),
      reason = character(),
      source = character(),
      stringsAsFactors = FALSE
    )
  }
  residual <- vc_obj$residual_sd
  if (!is.na(residual)) {
    table <- rbind(table, data.frame(
      group = "Residual",
      term_id = NA_character_,
      basis_lhs = "Residual",
      basis_rhs = "Residual",
      kind = "residual_variance",
      variance = residual^2,
      std_dev = residual,
      correlation = "",
      covariance_family = "residual",
      status = varcorr_status,
      reason = NA_character_,
      source = varcorr_source,
      stringsAsFactors = FALSE
    ))
  }
  unavailable <- if (has_fit_summary) {
    mm_report_unavailable_empty()
  } else {
    mm_report_unavailable(
      section = "random_effects",
      field = "stable_random_effect_variance_covariance_payload",
      status = "schema_gap",
      reason = "using_fit_varcorr_until_rust_report_payload_is_available",
      source = "planning/reporting_artifact_requirements.md"
    )
  }
  mm_report_section(table, varcorr_source, unavailable = unavailable)
}

mm_report_fixed_effects <- function(fit) {
  table <- inference_table(fit)$table
  table$source <- "fixed_effect_inference_table"
  mm_report_section(table, "fixed_effect_inference_table")
}

mm_report_fit_statistics <- function(fit) {
  table <- mm_report_kv(
    field = c("logLik", "deviance", "AIC", "BIC", "nobs", "df_residual",
              "sigma"),
    value = c(fit$logLik, fit$deviance, fit$AIC, fit$BIC, fit$nobs,
              fit$df_residual, fit$sigma),
    source = "mm_fit"
  )
  mm_report_section(table, "mm_fit")
}

mm_report_optimizer <- function(fit) {
  table <- optimizer_certificate(fit)$table
  table$source <- "OptimizerCertificate"
  mm_report_section(table, "OptimizerCertificate")
}

mm_report_comparison_ledger <- function(fit) {
  table <- data.frame(
    comparison_id = character(),
    model_id = character(),
    formula = character(),
    fit_method = character(),
    refit = logical(),
    comparison_method = character(),
    statistic = numeric(),
    df = numeric(),
    p_value = numeric(),
    status = character(),
    reason = character(),
    source = character(),
    stringsAsFactors = FALSE
  )
  unavailable <- mm_report_unavailable(
    section = "comparison_ledger",
    field = "comparison_ledger",
    status = "not_applicable",
    reason = "no_model_comparison_recorded_on_this_fit",
    source = "mixeff"
  )
  mm_report_section(table, "mixeff", unavailable = unavailable)
}

mm_report_reproducibility <- function(fit) {
  repro <- reproducibility(fit)
  raw <- repro$raw %||% list()
  table <- mm_report_kv(
    field = c("fit_intent", "random_state_used", "mixeff_version",
              "r_version", "platform", "contrasts"),
    value = c(
      mm_scalar_text(raw$fit_intent, "not_recorded"),
      mm_scalar_text(raw$random_state_used, "not_recorded"),
      mm_report_package_version(),
      paste(R.version$major, R.version$minor, sep = "."),
      R.version$platform,
      paste(getOption("contrasts"), collapse = ", ")
    ),
    source = c("reproducibility", "reproducibility", "R package metadata",
               "R session", "R session", "R options")
  )
  if (nrow(repro$thresholds)) {
    threshold_rows <- data.frame(
      field = paste0("threshold:", repro$thresholds$name),
      value = repro$thresholds$value,
      source = "reproducibility.thresholds",
      status = "available",
      reason = NA_character_,
      stringsAsFactors = FALSE
    )
    table <- rbind(table, threshold_rows)
  }
  mm_report_section(table, "mixed")
}

mm_report_unavailable_section <- function(fit) {
  mm_report_section(mm_report_unavailable_empty(), "mixeff")
}

mm_report_section <- function(table, source, unavailable = mm_report_unavailable_empty()) {
  list(
    table = mm_report_normalize_table(table, source),
    source = source,
    unavailable = mm_report_normalize_unavailable(unavailable)
  )
}

mm_report_normalize_table <- function(table, source) {
  if (!is.data.frame(table)) {
    mm_abort(
      message = "Report sections must be data frames.",
      class = "mm_schema_error",
      input = table
    )
  }

  n <- nrow(table)
  if (!"source" %in% names(table)) {
    table$source <- rep(as.character(source), length.out = n)
  } else {
    table$source <- mm_report_fill_text(table$source, source, n)
  }
  if (!"status" %in% names(table)) {
    table$status <- rep("available", n)
  } else {
    table$status <- mm_report_fill_text(table$status, "available", n)
  }
  if (!"reason" %in% names(table)) {
    table$reason <- rep(NA_character_, n)
  } else {
    table$reason <- mm_report_optional_text(table$reason, n)
  }

  needs_reason <- !is.na(table$status) & nzchar(table$status) &
    !table$status %in% c("available", "available_from_varcorr")
  missing_reason <- is.na(table$reason) | !nzchar(table$reason)
  table$reason[needs_reason & missing_reason] <- "reason_not_recorded"
  table
}

mm_report_normalize_unavailable <- function(unavailable) {
  if (!is.data.frame(unavailable)) {
    mm_abort(
      message = "Report unavailable ledgers must be data frames.",
      class = "mm_schema_error",
      input = unavailable
    )
  }
  required <- c("section", "field", "status", "reason", "source", "action_taken")
  for (col in setdiff(required, names(unavailable))) {
    unavailable[[col]] <- rep(NA_character_, nrow(unavailable))
  }
  unavailable <- unavailable[, required, drop = FALSE]
  for (col in required) {
    unavailable[[col]] <- mm_report_optional_text(unavailable[[col]], nrow(unavailable))
  }
  missing_reason <- is.na(unavailable$reason) | !nzchar(unavailable$reason)
  unavailable$reason[missing_reason] <- "reason_not_recorded"
  missing_action <- is.na(unavailable$action_taken) | !nzchar(unavailable$action_taken)
  unavailable$action_taken[missing_action] <- "reported"
  unavailable
}

mm_report_fill_text <- function(x, fallback, n) {
  out <- mm_report_optional_text(x, n)
  missing <- is.na(out) | !nzchar(out)
  out[missing] <- as.character(fallback)
  out
}

mm_report_optional_text <- function(x, n) {
  out <- as.character(x)
  if (length(out) != n) {
    out <- rep(out, length.out = n)
  }
  out[is.na(out) | !nzchar(out)] <- NA_character_
  out
}

mm_report_kv <- function(field, value, source, status = "available",
                         reason = NA_character_) {
  n <- length(field)
  data.frame(
    field = as.character(field),
    value = as.character(value),
    source = rep(as.character(source), length.out = n),
    status = rep(as.character(status), length.out = n),
    reason = rep(as.character(reason), length.out = n),
    stringsAsFactors = FALSE
  )
}

mm_report_unavailable <- function(section, field, status, reason, source) {
  data.frame(
    section = section,
    field = field,
    status = status,
    reason = reason,
    source = source,
    action_taken = "reported",
    stringsAsFactors = FALSE
  )
}

mm_report_unavailable_empty <- function() {
  data.frame(
    section = character(),
    field = character(),
    status = character(),
    reason = character(),
    source = character(),
    action_taken = character(),
    stringsAsFactors = FALSE
  )
}

mm_report_random_term_card_rows <- function(card) {
  blocks <- card$blocks %||% list()
  constraints <- mm_report_constraint_text(card$implied_constraints)
  if (!length(blocks)) {
    blocks <- list(list())
  }
  lapply(seq_along(blocks), function(i) {
    block <- blocks[[i]]
    data.frame(
      term_id = mm_scalar_text(card$term_id),
      original_fragment = mm_scalar_text(card$original_fragment),
      canonical_fragment = mm_scalar_text(card$canonical_fragment),
      group = mm_group_ir_label(card$group),
      block = i,
      basis = mm_list_text(block$basis),
      intercept = isTRUE(block$intercept),
      slopes = mm_list_text(block$slopes),
      covariance = mm_scalar_text(block$covariance),
      theta_parameters = as.integer(block$theta_parameters %||% NA_integer_),
      english = mm_scalar_text(block$english),
      constraints = constraints,
      design_status = mm_scalar_text(card$design_support$status, "not_assessed"),
      source = "random_term_cards",
      stringsAsFactors = FALSE
    )
  })
}

mm_report_cross_card_rows <- function(constraints) {
  rows <- lapply(constraints %||% list(), function(constraint) {
    data.frame(
      term_id = mm_list_text(constraint$between_cards),
      original_fragment = "",
      canonical_fragment = "",
      group = "",
      block = NA_integer_,
      basis = mm_list_text(constraint$between_basis),
      intercept = NA,
      slopes = "",
      covariance = "",
      theta_parameters = NA_integer_,
      english = "",
      constraints = mm_scalar_text(constraint$reason),
      design_status = "constraint",
      source = "cross_card_constraints",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

mm_report_constraint_text <- function(constraints) {
  constraints <- constraints %||% list()
  if (!length(constraints)) {
    return("")
  }
  paste(vapply(constraints, function(constraint) {
    reason <- mm_scalar_text(constraint$reason)
    between <- mm_list_text(constraint$between)
    if (nzchar(between)) {
      paste0(reason, " [", between, "]")
    } else {
      reason
    }
  }, character(1)), collapse = "; ")
}

mm_group_ir_label <- function(group) {
  if (is.null(group)) return("")
  for (path in list(
    c("name"),
    c("single", "name"),
    c("interaction", "names"),
    c("cell", "names")
  )) {
    value <- group
    for (part in path) {
      value <- value[[part]]
      if (is.null(value)) break
    }
    if (!is.null(value)) {
      return(mm_list_text(value))
    }
  }
  mm_scalar_text(group)
}

mm_report_group_roles <- function(fit, groups) {
  audit <- audit_design(fit)
  card_roles <- lapply(audit$random_term_cards %||% list(), function(card) {
    data.frame(
      group = mm_group_ir_label(card$group),
      role = mm_scalar_text(card$role_origin$role, "not_recorded"),
      stringsAsFactors = FALSE
    )
  })
  roles <- if (length(card_roles)) {
    do.call(rbind, card_roles)
  } else {
    data.frame(group = character(), role = character(), stringsAsFactors = FALSE)
  }
  vapply(groups, function(group) {
    hit <- roles$role[roles$group == group]
    if (length(hit)) hit[[1L]] else "not_recorded"
  }, character(1))
}

mm_report_group_max_rows <- function(fit, group) {
  out <- tryCatch({
    factor <- mm_group_factor(fit$model_frame, group)
    max(tabulate(factor), na.rm = TRUE)
  }, error = function(cnd) NA_integer_)
  as.integer(out)
}

mm_report_package_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("mixeff")),
    error = function(cnd) NA_character_
  )
}
