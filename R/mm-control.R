#' Control mixeff fitting behavior
#'
#' `mm_control()` collects small R-side controls for [lmm()] and [glmm()].
#' `verbose = -1` suppresses the pre-fit [explain_model()] message;
#' non-negative values emit it once before optimization (it travels on the
#' message stream, so `suppressMessages()` and knitr's `message = FALSE`
#' also quiet it).
#'
#' By default the fit driver selects the optimizer and its tolerances
#' automatically (see [optimizer_certificate()] to inspect what ran). The
#' `optimizer`, `start`, and `ftol_*`/`xtol_rel` arguments are a narrow,
#' opt-in escape hatch — for recourse when the default fails to converge, for
#' warm starts, and for explicit tolerance overrides. Any override you supply is
#' recorded in the optimizer certificate, so the fit stays auditable.
#'
#' @param verbose Integer verbosity level. Use `-1` to suppress the automatic
#'   model explanation (and the GLMM estimator notice).
#' @param max_feval Optional positive integer capping the optimizer's objective
#'   evaluations. Most useful for [glmm()] with `method = "joint_laplace"`,
#'   whose native joint optimizer otherwise runs to an engine-chosen budget.
#'   `NULL` (default) leaves the engine default in place.
#' @param optimizer Optional optimizer name, overriding the driver's automatic
#'   choice. One of `"auto"` (default behavior), `"bobyqa"`, `"newuoa"`,
#'   `"cobyla"`, `"pattern_search"`, `"trust_bq"`, or the PRIMA variants
#'   (`"prima_bobyqa"`, `"prima_cobyla"`, `"prima_lincoa"`, `"prima_newuoa"`).
#'   An unsupported or not-compiled choice raises a typed error rather than
#'   silently falling back. `NULL`/`"auto"` keep automatic selection.
#' @param start Optional numeric warm-start vector for the covariance
#'   parameters (theta). Its length must match the model's theta dimension
#'   (the engine validates this). `NULL` (default) cold-starts.
#' @param ftol_rel,ftol_abs Optional positive relative/absolute convergence
#'   tolerances on the objective. `NULL` keeps the engine default.
#' @param xtol_rel Optional positive relative convergence tolerance on the
#'   optimizer parameters. `NULL` keeps the engine default.
#' @param optCtrl Optional named list in the style of lme4's
#'   `lmerControl(optCtrl = )`, translated to the engine's controls:
#'   `maxfun`/`maxeval` -> `max_feval`; `ftol_rel`/`ftol_abs`/`xtol_rel` ->
#'   the arguments of the same name; `xtol_abs` (scalar or one per theta)
#'   -> the absolute parameter tolerance; `rhobeg` (scalar or one per
#'   theta) -> the optimizer's initial step. Unknown names are refused with
#'   an `mm_arg_error` rather than ignored. lme4's `check.conv.*` options
#'   have no counterpart: mixeff does not emit post-hoc convergence
#'   warnings; the fit's convergence status is the engine's typed
#'   certificate ([optimizer_certificate()], [verify_convergence()]).
#'
#' @return A list of class `mm_control`.
#'
#' @seealso [optimizer_certificate()] to inspect which optimizer ran and whether
#'   a caller override was applied.
#'
#' @export
mm_control <- function(verbose = 0L, max_feval = NULL, optimizer = NULL,
                       start = NULL, ftol_rel = NULL, ftol_abs = NULL,
                       xtol_rel = NULL, optCtrl = NULL) {
  if (!is.null(optCtrl)) {
    oc <- mm_translate_optctrl(optCtrl)
    if (!is.null(oc$max_feval)) max_feval <- max_feval %||% oc$max_feval
    if (!is.null(oc$ftol_rel)) ftol_rel <- ftol_rel %||% oc$ftol_rel
    if (!is.null(oc$ftol_abs)) ftol_abs <- ftol_abs %||% oc$ftol_abs
    if (!is.null(oc$xtol_rel)) xtol_rel <- xtol_rel %||% oc$xtol_rel
  } else {
    oc <- list()
  }
  if (!is.numeric(verbose) || length(verbose) != 1L || is.na(verbose)) {
    mm_abort(
      message = "`verbose` must be a single numeric value.",
      class = "mm_arg_error",
      input = verbose
    )
  }
  out <- list(verbose = as.integer(verbose))

  if (!is.null(max_feval)) {
    if (!is.numeric(max_feval) || length(max_feval) != 1L ||
        is.na(max_feval) || max_feval < 1) {
      mm_abort(
        message = "`max_feval` must be `NULL` or a single positive integer.",
        class = "mm_arg_error",
        input = max_feval
      )
    }
    out$max_feval <- as.integer(max_feval)
  }

  if (!is.null(optimizer)) {
    if (!is.character(optimizer) || length(optimizer) != 1L ||
        is.na(optimizer) || !nzchar(optimizer)) {
      mm_abort(
        message = "`optimizer` must be `NULL` or a single optimizer name.",
        class = "mm_arg_error",
        input = optimizer
      )
    }
    out$optimizer <- optimizer
  }

  if (!is.null(start)) {
    if (!is.numeric(start) || !length(start) || anyNA(start) ||
        any(!is.finite(start))) {
      mm_abort(
        message = "`start` must be `NULL` or a finite numeric vector (warm-start theta).",
        class = "mm_arg_error",
        input = start
      )
    }
    out$start <- as.numeric(start)
  }

  for (nm in c("ftol_rel", "ftol_abs", "xtol_rel")) {
    val <- get(nm)
    if (!is.null(val)) {
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          !is.finite(val) || val <= 0) {
        mm_abort(
          message = sprintf("`%s` must be `NULL` or a single positive number.", nm),
          class = "mm_arg_error",
          input = val
        )
      }
      out[[nm]] <- as.numeric(val)
    }
  }

  # Engine-only vector controls reachable through optCtrl.
  if (!is.null(oc$xtol_abs)) out$xtol_abs <- oc$xtol_abs
  if (!is.null(oc$initial_step)) out$initial_step <- oc$initial_step

  class(out) <- "mm_control"
  out
}

# lme4/nloptr-style optCtrl -> engine control fields. Every name is mapped or
# refused; nothing is silently ignored.
mm_translate_optctrl <- function(optCtrl) {
  if (!is.list(optCtrl) || (length(optCtrl) && is.null(names(optCtrl))) ||
      any(!nzchar(names(optCtrl)))) {
    mm_abort(message = "`optCtrl` must be a named list.",
             class = "mm_arg_error", input = optCtrl)
  }
  map <- c(maxfun = "max_feval", maxeval = "max_feval",
           ftol_rel = "ftol_rel", ftol_abs = "ftol_abs",
           xtol_rel = "xtol_rel", xtol_abs = "xtol_abs",
           rhobeg = "initial_step")
  unknown <- setdiff(names(optCtrl), names(map))
  if (length(unknown)) {
    mm_abort(
      message = sprintf(
        paste0("`optCtrl` entries not supported by the engine: %s. ",
               "Supported: %s."),
        paste(sprintf("`%s`", unknown), collapse = ", "),
        paste(names(map), collapse = ", ")
      ),
      class = "mm_arg_error",
      input = optCtrl
    )
  }
  out <- list()
  for (nm in names(optCtrl)) {
    val <- optCtrl[[nm]]
    vector_ok <- map[[nm]] %in% c("xtol_abs", "initial_step")
    if (!is.numeric(val) || !length(val) || anyNA(val) || any(!is.finite(val)) ||
        any(val <= 0) || (!vector_ok && length(val) != 1L)) {
      mm_abort(
        message = sprintf("`optCtrl$%s` must be %s positive finite number%s.",
                          nm, if (vector_ok) "a vector of" else "a single",
                          if (vector_ok) "s" else ""),
        class = "mm_arg_error",
        input = val
      )
    }
    out[[map[[nm]]]] <- as.numeric(val)
  }
  out
}

mm_validate_control <- function(control) {
  if (missing(control) || is.null(control)) {
    return(mm_control())
  }
  if (!is.list(control)) {
    mm_abort(
      message = "`control` must be a list created by mm_control().",
      class = "mm_arg_error",
      input = control
    )
  }
  mm_control(
    verbose = control$verbose %||% 0L,
    max_feval = control$max_feval,
    optimizer = control$optimizer,
    start = control$start,
    ftol_rel = control$ftol_rel,
    ftol_abs = control$ftol_abs,
    xtol_rel = control$xtol_rel,
    optCtrl = c(
      if (!is.null(control$xtol_abs)) list(xtol_abs = control$xtol_abs),
      if (!is.null(control$initial_step)) list(rhobeg = control$initial_step)
    )
  )
}
