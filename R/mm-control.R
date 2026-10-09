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
#'   model explanation and the fit notices (the GLMM estimator notice, the
#'   `mm_rows_dropped` missing-value notice, grouping coercion and scaling
#'   advisories).
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
#' @param disp_method,disp_dof_correction,max_phi_iter lme4 2.1-0's
#'   `glmerControl()` dispersion controls for [glmm()] fits with a free
#'   dispersion parameter (Gamma, inverse Gaussian, Gaussian with a
#'   non-identity link). `disp_method = "moment"` (the default) profiles the
#'   dispersion `phi` in a damped fixed-point loop around PIRLS, with the
#'   working weights divided by `phi`; `"old/buggy"` restores lme4 < 2.1
#'   (working weights with `phi = 1`, theta relative to `sigma()`).
#'   `disp_dof_correction` (default `TRUE`) divides the deviance by
#'   `n - rank([X, Z])` rather than `n` in the moment estimator, and
#'   `max_phi_iter` (lme4's `maxPhiIter`, default 100) caps the `phi`
#'   iterations per objective evaluation. Both only matter for
#'   `"moment"`. `NULL` keeps the engine default (the lme4 2.1-0 default).
#'   Invalid values raise an `mm_arg_error`. For families without a free
#'   dispersion (binomial, Poisson, negative binomial) and for [lmm()] the
#'   settings have no effect: they are kept on the control, and the fit
#'   announces that they were ignored with an `mm_control_ignored_notice`
#'   message (silenced by `verbose = -1`). lme4 ignores them silently.
#'
#' @return A list of class `mm_control`.
#'
#' @seealso [optimizer_certificate()] to inspect which optimizer ran and whether
#'   a caller override was applied.
#'
#' @export
mm_control <- function(verbose = 0L, max_feval = NULL, optimizer = NULL,
                       start = NULL, ftol_rel = NULL, ftol_abs = NULL,
                       xtol_rel = NULL, optCtrl = NULL,
                       disp_method = NULL, disp_dof_correction = NULL,
                       max_phi_iter = NULL) {
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

  if (!is.null(disp_method)) {
    if (!is.character(disp_method) || length(disp_method) != 1L ||
        is.na(disp_method) || !disp_method %in% mm_disp_methods) {
      mm_abort(
        message = sprintf(
          "`disp_method` must be `NULL`, %s.",
          paste(sprintf("\"%s\"", mm_disp_methods), collapse = " or ")
        ),
        class = "mm_arg_error",
        input = disp_method
      )
    }
    out$disp_method <- disp_method
  }
  if (!is.null(disp_dof_correction)) {
    if (!is.logical(disp_dof_correction) || length(disp_dof_correction) != 1L ||
        is.na(disp_dof_correction)) {
      mm_abort(
        message = "`disp_dof_correction` must be `NULL`, `TRUE` or `FALSE`.",
        class = "mm_arg_error",
        input = disp_dof_correction
      )
    }
    out$disp_dof_correction <- disp_dof_correction
  }
  if (!is.null(max_phi_iter)) {
    if (!is.numeric(max_phi_iter) || length(max_phi_iter) != 1L ||
        is.na(max_phi_iter) || !is.finite(max_phi_iter) ||
        max_phi_iter < 1 || max_phi_iter != round(max_phi_iter) ||
        max_phi_iter > .Machine$integer.max) {
      mm_abort(
        message = "`max_phi_iter` must be `NULL` or a single positive integer.",
        class = "mm_arg_error",
        input = max_phi_iter
      )
    }
    out$max_phi_iter <- as.integer(max_phi_iter)
  }

  # Engine-only vector controls reachable through optCtrl.
  if (!is.null(oc$xtol_abs)) out$xtol_abs <- oc$xtol_abs
  if (!is.null(oc$initial_step)) out$initial_step <- oc$initial_step

  class(out) <- "mm_control"
  out
}

# lme4 2.1-0 glmerControl(disp_method = ) choices.
mm_disp_methods <- c("moment", "old/buggy")

mm_disp_control_names <- c("disp_method", "disp_dof_correction", "max_phi_iter")

# The dispersion controls the caller set, or character() if none.
mm_disp_controls_set <- function(control) {
  mm_disp_control_names[vapply(mm_disp_control_names, function(nm) {
    !is.null(control[[nm]])
  }, logical(1))]
}

# Dispersion controls only act on GLMMs with a free dispersion. Elsewhere
# they are kept but have no effect; say so (typed notice, verbose >= 0)
# rather than ignore them silently.
mm_disp_control_notice <- function(control, applies, what) {
  set <- mm_disp_controls_set(control)
  if (!length(set) || applies || !isTRUE(control$verbose >= 0L)) {
    return(invisible(NULL))
  }
  mm_inform(
    sprintf(
      paste0(
        "mm_control(%s) ignored: the dispersion controls only affect GLMMs ",
        "with a free dispersion parameter (Gamma, inverse Gaussian, Gaussian ",
        "with a non-identity link), not %s."
      ),
      paste(set, collapse = ", "), what
    ),
    class = "mm_control_ignored_notice"
  )
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
    disp_method = control$disp_method,
    disp_dof_correction = control$disp_dof_correction,
    max_phi_iter = control$max_phi_iter,
    optCtrl = c(
      if (!is.null(control$xtol_abs)) list(xtol_abs = control$xtol_abs),
      if (!is.null(control$initial_step)) list(rhobeg = control$initial_step)
    )
  )
}
