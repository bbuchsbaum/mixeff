#' Verify convergence of a fitted mixed model
#'
#' `verify_convergence()` re-runs the fit under the engine's bounded
#' verification workflow and reports whether the extra runs agree with the
#' fitted optimum: a restart from the optimum, one or more jittered restarts,
#' and (LMM-only, opt-in) an alternate-optimizer consensus pass. Agreement is
#' judged by the engine against the objective/theta/beta tolerances below;
#' the verdict (`status`), the per-run deltas, and the wording are all owned
#' by the Rust contract — R only formats them.
#'
#' The verifier refits the model from the stored specification before it
#' starts, so a call costs roughly `2 + jitter_starts` fits (plus consensus
#' runs when enabled). GLMM verification uses wider default tolerances than
#' the LMM route (objective `1e-4`, beta `1e-3`): each objective evaluation
#' carries inner-PIRLS noise, and the joint path's derivative-free beta
#' search bounds beta reproducibility at its `ftol` stop. A GLMM whose
#' requested `joint_laplace` estimator was substituted by the labelled
#' fast-PIRLS fallback (see `optimizer_certificate()`) re-requests the joint
#' route; when the verification refit falls back the same way — the expected
#' case — every run verifies the profiled objective the fitted numbers came
#' from, with ordinary objective deltas. Runs are only reported as
#' substitution runs (no objective delta) in the asymmetric case where the
#' reference refit certifies the joint route but an individual verification
#' run falls back.
#'
#' @param fit A fitted `mm_lmm` from [lmm()] or `mm_glmm` from [glmm()].
#' @param ... Reserved for future methods.
#' @param restart Logical; re-optimize starting from the fitted optimum and
#'   compare against it.
#' @param jitter_starts Number of restarts from jittered copies of the fitted
#'   covariance parameters.
#' @param jitter_scale Relative scale of the jitter applied to theta.
#' @param consensus Logical; also refit with an engine-chosen alternate
#'   optimizer and compare. Default `FALSE`: this vendored build compiles
#'   without the optional `nlopt` backend, and for some models the engine's
#'   alternate choice is an nlopt optimizer — its absence would then be
#'   reported as a non-agreeing run (status `fragile`) that reflects the
#'   build, not the fit. Enable it when you want the consensus pass and will
#'   read the per-run diagnostics.
#' @param max_feval Positive integer cap on objective evaluations per
#'   verification run.
#' @param objective_tolerance,theta_tolerance,beta_tolerance Positive
#'   agreement tolerances on the objective value, the covariance parameters,
#'   and the fixed effects.
#'
#' @return An object of class `mm_convergence_verification` carrying:
#' \describe{
#'   \item{`status`}{the engine verdict: `not_run`, `restart_agrees`,
#'     `optimizer_consensus`, `fragile`, or `unstable`}
#'   \item{`message`}{the engine's one-line summary}
#'   \item{`table`}{a data frame with one row per verification run (label,
#'     optimizer, return code, objective/theta/beta deltas, agreement)}
#'   \item{`reference`}{the reference optimum the runs were compared to}
#'   \item{`tolerances`}{the agreement tolerances that were applied}
#'   \item{`raw`}{the parsed engine payload}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- lmm(y ~ t + (1 | s), df)
#' verify_convergence(fit)
#' }
#'
#' @seealso [optimizer_certificate()] for what the original fit ran;
#'   [mm_control()] to refit with a different optimizer or tolerances.
#'
#' @export
verify_convergence <- function(fit, ...) {
  UseMethod("verify_convergence")
}

#' @rdname verify_convergence
#' @export
verify_convergence.default <- function(fit, ...) {
  mm_abort(
    message = sprintf(
      "verify_convergence() supports mixed models fitted by lmm() or glmm(); got <%s>.",
      paste(class(fit), collapse = "/")
    ),
    class = "mm_schema_error",
    input = fit
  )
}

#' @rdname verify_convergence
#' @export
verify_convergence.mm_lmm <- function(fit, ...,
                                      restart = TRUE,
                                      jitter_starts = 1L,
                                      jitter_scale = 1e-4,
                                      consensus = FALSE,
                                      max_feval = 500L,
                                      objective_tolerance = 1e-5,
                                      theta_tolerance = 1e-3,
                                      beta_tolerance = 1e-4) {
  mm_verify_check_flag(restart, "restart")
  mm_verify_check_flag(consensus, "consensus")
  mm_verify_check_count(jitter_starts, "jitter_starts", min = 0)
  mm_verify_check_count(max_feval, "max_feval", min = 1)
  mm_verify_check_positive(jitter_scale, "jitter_scale")
  mm_verify_check_positive(objective_tolerance, "objective_tolerance")
  mm_verify_check_positive(theta_tolerance, "theta_tolerance")
  mm_verify_check_positive(beta_tolerance, "beta_tolerance")

  payload <- mm_rust_fit_bridge_payload(fit)
  payload$REML <- isTRUE(fit$REML)
  options_json <- jsonlite::toJSON(
    list(
      restart_from_optimum = isTRUE(restart),
      jitter_starts = as.integer(jitter_starts),
      jitter_scale = as.numeric(jitter_scale),
      run_optimizer_consensus = isTRUE(consensus),
      max_function_evaluations = as.integer(max_feval),
      objective_tolerance = as.numeric(objective_tolerance),
      theta_tolerance = as.numeric(theta_tolerance),
      beta_tolerance = as.numeric(beta_tolerance)
    ),
    auto_unbox = TRUE
  )

  json <- tryCatch(
    mm_verify_convergence_json(payload, as.character(options_json)),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json)
  }
  parsed <- mm_json_parse_convergence_verification(json)

  obj <- list(
    status = mm_scalar_text(parsed$status, "not_run"),
    message = mm_scalar_text(parsed$message),
    table = mm_verify_runs_table(parsed$runs %||% list()),
    reference = list(
      objective = parsed$reference_objective,
      theta = as.numeric(unlist(parsed$reference_theta %||% list())),
      beta = as.numeric(unlist(parsed$reference_beta %||% list())),
      effective_ranks = as.integer(unlist(parsed$reference_effective_ranks %||%
                                            list()))
    ),
    tolerances = list(
      objective = parsed$objective_tolerance,
      theta = parsed$theta_tolerance,
      beta = parsed$beta_tolerance
    ),
    raw = parsed
  )
  class(obj) <- "mm_convergence_verification"
  obj
}

#' @rdname verify_convergence
#' @export
verify_convergence.mm_glmm <- function(fit, ...,
                                       restart = TRUE,
                                       jitter_starts = 1L,
                                       jitter_scale = 1e-4,
                                       consensus = FALSE,
                                       max_feval = 1000L,
                                       objective_tolerance = 1e-4,
                                       theta_tolerance = 1e-3,
                                       beta_tolerance = 1e-3) {
  # GLMM defaults mirror the engine's glmm_defaults(): a larger evaluation
  # budget (each objective evaluation pays for inner PIRLS solves), an
  # objective tolerance sized to the inner-PIRLS noise floor, and a beta
  # tolerance sized to the joint path's derivative-free beta search. The
  # optimizer-consensus pass is LMM-only (the GLMM fit drivers own optimizer
  # selection), so consensus = TRUE is refused rather than ignored.
  mm_verify_check_flag(restart, "restart")
  mm_verify_check_flag(consensus, "consensus")
  if (isTRUE(consensus)) {
    mm_abort(
      message = paste0(
        "The optimizer-consensus verification pass is not available for ",
        "GLMM fits; call verify_convergence() with consensus = FALSE."
      ),
      class = "mm_arg_error",
      input = consensus
    )
  }
  mm_verify_check_count(jitter_starts, "jitter_starts", min = 0)
  mm_verify_check_count(max_feval, "max_feval", min = 1)
  mm_verify_check_positive(jitter_scale, "jitter_scale")
  mm_verify_check_positive(objective_tolerance, "objective_tolerance")
  mm_verify_check_positive(theta_tolerance, "theta_tolerance")
  mm_verify_check_positive(beta_tolerance, "beta_tolerance")

  payload <- mm_rust_glmm_verify_payload(fit)
  options_json <- jsonlite::toJSON(
    list(
      restart_from_optimum = isTRUE(restart),
      jitter_starts = as.integer(jitter_starts),
      jitter_scale = as.numeric(jitter_scale),
      run_optimizer_consensus = FALSE,
      max_function_evaluations = as.integer(max_feval),
      objective_tolerance = as.numeric(objective_tolerance),
      theta_tolerance = as.numeric(theta_tolerance),
      beta_tolerance = as.numeric(beta_tolerance)
    ),
    auto_unbox = TRUE
  )

  json <- tryCatch(
    mm_verify_convergence_glmm_json(payload, as.character(options_json)),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) {
    mm_abort_from_bridge(json)
  }
  parsed <- mm_json_parse_convergence_verification(json)

  obj <- list(
    status = mm_scalar_text(parsed$status, "not_run"),
    message = mm_scalar_text(parsed$message),
    table = mm_verify_runs_table(parsed$runs %||% list()),
    reference = list(
      objective = parsed$reference_objective,
      theta = as.numeric(unlist(parsed$reference_theta %||% list())),
      beta = as.numeric(unlist(parsed$reference_beta %||% list())),
      effective_ranks = as.integer(unlist(parsed$reference_effective_ranks %||%
                                            list()))
    ),
    tolerances = list(
      objective = parsed$objective_tolerance,
      theta = parsed$theta_tolerance,
      beta = parsed$beta_tolerance
    ),
    raw = parsed
  )
  class(obj) <- "mm_convergence_verification"
  obj
}

## Bridge payload for the GLMM verification refit. Carries everything the
## engine needs to rebuild the SAME model: the post-prep model frame and
## formula (grouped-binomial responses are already proportion + trial
## weights here), the exact engine family string recorded at fit time
## (bernoulli/binomial split, NB theta-mode spelling), and the REQUESTED
## estimator -- for a substituted (fallback) fit the verification re-requests
## the joint route; when the refit falls back the same way, all runs verify
## the profiled objective the fitted numbers came from.
mm_rust_glmm_verify_payload <- function(fit) {
  payload <- mm_rust_fit_bridge_payload(fit)
  payload$offset <- mm_bridge_weights(fit$offset)
  payload$family <- as.character(
    fit$engine_family %||% mm_glmm_engine_family_fallback(fit)
  )
  payload$link <- as.character(fit$family$link)
  method <- as.character(fit$method %||% "pirls_profiled")
  sub <- mm_estimator_substitution(fit)
  if (!is.null(sub) && !is.na(sub$requested_method)) {
    method <- sub$requested_method
  }
  if (identical(method, "fast_pirls_profiled")) {
    method <- "pirls_profiled"
  }
  payload$method <- method
  payload$n_agq <- as.numeric(fit$nAGQ %||% 1L)
  payload
}

## Engine-family reconstruction for fits saved before `engine_family` was
## stored on the object (pre-0.2.0 artifacts). Mirrors the fit-time prep
## rules: NB theta mode rides the family string; a binomial response without
## prior weights was fitted as bernoulli.
mm_glmm_engine_family_fallback <- function(fit) {
  info <- fit$family %||% list()
  fam <- as.character(info$family %||% NA_character_)
  if (identical(fam, "negative_binomial")) {
    return(mm_glmm_nb_engine_spec(info))
  }
  if (identical(fam, "binomial")) {
    return(if (is.null(fit$weights)) "bernoulli" else "binomial")
  }
  fam
}

#' @method print mm_convergence_verification
#' @export
print.mm_convergence_verification <- function(x, ...) {
  cat(sprintf("Convergence verification (status: %s)\n", x$status))
  if (nzchar(x$message)) {
    cat(sprintf("  %s\n", x$message))
  }
  if (nrow(x$table)) {
    cat("Runs:\n")
    out <- x$table
    for (col in c("objective_value", "objective_delta", "theta_delta",
                  "beta_delta")) {
      out[[col]] <- signif(out[[col]], 6)
    }
    if (!any(nzchar(out$diagnostics))) {
      out$diagnostics <- NULL
    }
    print(out, row.names = FALSE)
  }
  cat(sprintf(
    "Tolerances: objective %s; theta %s; beta %s\n",
    format(x$tolerances$objective %||% NA),
    format(x$tolerances$theta %||% NA),
    format(x$tolerances$beta %||% NA)
  ))
  invisible(x)
}

mm_json_parse_convergence_verification <- function(json) {
  if (!is.character(json) || length(json) != 1L || is.na(json) ||
      !nzchar(json)) {
    mm_abort(
      message = "`json` must be a single non-empty character string.",
      class = "mm_schema_error",
      input = json
    )
  }
  tryCatch(
    jsonlite::fromJSON(json, simplifyVector = FALSE),
    error = function(cnd) {
      mm_abort(
        message = sprintf(
          "Failed to parse convergence-verification JSON: %s",
          conditionMessage(cnd)
        ),
        class = "mm_schema_error",
        input = json,
        parent = cnd
      )
    }
  )
}

mm_verify_runs_table <- function(runs) {
  if (!length(runs)) {
    return(data.frame(
      label = character(),
      optimizer = character(),
      return_code = character(),
      objective_value = numeric(),
      objective_delta = numeric(),
      theta_delta = numeric(),
      beta_delta = numeric(),
      agrees = logical(),
      diagnostics = character(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(runs, function(run) {
    data.frame(
      label = mm_scalar_text(run$label),
      optimizer = mm_scalar_text(run$optimizer_name),
      return_code = mm_scalar_text(run$return_code),
      objective_value = mm_verify_number(run$objective_value),
      objective_delta = mm_verify_number(run$objective_delta),
      theta_delta = mm_verify_number(run$max_abs_theta_delta),
      beta_delta = mm_verify_number(run$max_abs_beta_delta),
      agrees = isTRUE(run$agrees),
      diagnostics = mm_list_text(run$diagnostics),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

mm_verify_number <- function(value) {
  if (is.numeric(value) && length(value) == 1L) as.numeric(value) else NA_real_
}

mm_verify_check_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    mm_abort(
      message = sprintf("`%s` must be TRUE or FALSE.", name),
      class = "mm_arg_error",
      input = value
    )
  }
}

mm_verify_check_count <- function(value, name, min) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value < min || value != as.integer(value)) {
    mm_abort(
      message = sprintf("`%s` must be a single integer >= %d.", name, min),
      class = "mm_arg_error",
      input = value
    )
  }
}

mm_verify_check_positive <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value <= 0) {
    mm_abort(
      message = sprintf("`%s` must be a single positive number.", name),
      class = "mm_arg_error",
      input = value
    )
  }
}
