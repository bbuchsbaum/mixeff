# lme4-identical coefficient naming --------------------------------------
#
# The Rust engine labels fixed-effect coefficients in its own encoding
# ("recipe: B", "recipe: B:temperature: .L") and orders interaction columns
# with the LAST factor varying fastest, whereas R's model.matrix() (and
# therefore lme4) uses "recipeB", "recipeB:temperature.L" with the FIRST
# factor varying fastest. mixeff fits carry lme4-identical names and order on
# every programmatic surface (fixef, coef tables, vcov dimnames, lincombs,
# emmeans, broom) so downstream code written against lme4 is drop-in
# compatible; the engine encoding survives only inside engine-rendered
# explain/audit text and at the FFI boundary, where `fit$coef_map` translates.
#
# `fit$coef_map` fields (p = number of estimated fixed-effect coefficients):
#   engine_names  character(p), engine encoding, ENGINE order (beta_names as
#                 received from the fit result)
#   lme4_names    character(p), R/lme4 encoding, R model.matrix column order
#   perm          integer(p): lme4 position j corresponds to engine index
#                 perm[j], i.e. beta_lme4 = beta_engine[perm]
#   assign        integer(p), lme4 order: the terms() assign index of each
#                 coefficient (0 = intercept), for term -> coefficient
#                 selection without name parsing
#   term_labels   attr(terms, "term.labels") of the fixed formula

# Translate R-style model.matrix column names into the engine encoding.
# Factor components gain the engine's "var: level" separator; logical
# components drop model.matrix's "TRUE" suffix (the wire coerces logicals to
# numeric 0/1, so the engine emits the bare name for the same column R calls
# "varTRUE"); numeric components pass through. `factor_vars` must be sorted
# longest-first so a factor whose name prefixes another ("rec" vs "recipe")
# matches greedily.
mm_engine_encode_names <- function(r_names, factor_vars, logical_vars = character(0)) {
  translate_component <- function(comp) {
    for (v in factor_vars) {
      if (startsWith(comp, v)) {
        lev <- substring(comp, nchar(v) + 1L)
        if (nzchar(lev)) return(paste0(v, ": ", lev))
      }
    }
    for (v in logical_vars) {
      if (identical(comp, paste0(v, "TRUE"))) return(v)
    }
    comp
  }
  vapply(r_names, function(nm) {
    if (nm == "(Intercept)") return(nm)
    parts <- strsplit(nm, ":", fixed = TRUE)[[1L]]
    paste(vapply(parts, translate_component, character(1)), collapse = ":")
  }, character(1))
}

# Build the engine <-> lme4 coefficient map for a fit-like object that has
# `$artifact`, `$formula`, and `$model_frame` (mm_fixed_formula and the
# training-basis helpers need those) plus engine-named `engine_names`.
# Aborts (typed) if any engine coefficient cannot be reconstructed from R's
# design: silent misalignment corrupts every downstream number, so the only
# acceptable failure mode is a loud one.
mm_coef_name_map <- function(fit, engine_names = names(fit$beta)) {
  rhs <- stats::delete.response(stats::terms(mm_fixed_formula(fit)))
  # Forcing the fit-time per-factor coding over a factor that carries an
  # attached contrast matrix makes model.frame warn "contrasts dropped"; the
  # engine ignores attached matrices on unordered factors by documented
  # design (ordered ones are guarded at translate time), so the drop is
  # intentional here, not information loss.
  muffle_contrasts_dropped <- function(w) {
    if (grepl("contrasts dropped from factor", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  }
  X <- withCallingHandlers(
    {
      mf <- stats::model.frame(rhs, data = fit$model_frame,
                               na.action = stats::na.pass,
                               xlev = mm_training_xlevels(fit))
      stats::model.matrix(rhs, data = mf,
                          contrasts.arg = mm_training_contrasts(fit))
    },
    warning = muffle_contrasts_dropped
  )
  r_names <- colnames(X)
  x_assign <- attr(X, "assign")

  fe_vars <- all.vars(rhs)
  is_fac <- vapply(fit$model_frame, is.factor, logical(1))
  factor_vars <- intersect(names(is_fac)[is_fac], fe_vars)
  factor_vars <- factor_vars[order(nchar(factor_vars), decreasing = TRUE)]
  is_lgl <- vapply(fit$model_frame, is.logical, logical(1))
  logical_vars <- intersect(names(is_lgl)[is_lgl], fe_vars)

  encoded <- mm_engine_encode_names(r_names, factor_vars, logical_vars)

  # The engine pivots rank-deficient columns out of the fit (like lme4's
  # rank-deficiency drop), so engine_names may be a subset of R's columns.
  keep <- encoded %in% engine_names
  matched <- match(encoded[keep], engine_names)
  unmatched <- setdiff(engine_names, encoded)
  if (length(unmatched) || anyDuplicated(matched)) {
    mm_abort(
      message = paste0(
        "Could not build the lme4-compatible coefficient name map; engine ",
        "coefficient(s) with no unique R design column: ",
        paste(unmatched %||% engine_names[duplicated(matched)], collapse = ", "),
        ". The engine's fixed-effect design does not reconstruct from R's ",
        "model.matrix() for this formula; please report it with a reproducer."
      ),
      class = "mm_schema_error",
      expected = engine_names,
      observed = encoded
    )
  }

  list(
    engine_names = engine_names,
    lme4_names   = r_names[keep],
    perm         = matched,
    assign       = as.integer(x_assign[keep]),
    term_labels  = attr(stats::terms(rhs), "term.labels")
  )
}

# Fetch the map from a fit, loudly. Fits constructed by this version always
# carry one; an absent map means a foreign/legacy object we must not guess at.
mm_fit_coef_map <- function(fit) {
  map <- fit$coef_map
  if (is.null(map)) {
    mm_abort(
      message = paste0(
        "This fit object has no coefficient name map (`coef_map`); it was ",
        "likely saved by an older mixeff version. Re-fit the model to use ",
        "this function."
      ),
      class = "mm_schema_error"
    )
  }
  map
}

# Rename + reorder an engine-ordered vector into the lme4 basis.
mm_coef_apply_map <- function(x, map) {
  stats::setNames(as.numeric(x)[map$perm], map$lme4_names)
}

# Translate engine-encoded labels to lme4 names; labels that are not fixed
# effect coefficients (theta rows, "sigma", term names) pass through.
mm_coef_engine_to_lme4 <- function(labels, map) {
  if (is.null(map) || !length(labels)) return(labels)
  dict <- stats::setNames(map$lme4_names, map$engine_names[map$perm])
  hit <- labels %in% names(dict)
  labels[hit] <- dict[labels[hit]]
  labels
}

# Permute an L matrix whose columns align with the lme4-ordered coefficients
# into the engine's column order for the FFI boundary.
mm_coef_l_to_engine <- function(L, fit) {
  map <- fit$coef_map
  if (is.null(map)) return(L)
  out <- L
  out[, map$perm] <- L
  colnames(out) <- map$engine_names
  out
}

# Strip the engine's "var: level" separator for random-effect column labels
# ("modality: Audio-only" -> "modalityAudio-only", "temperature: .L" ->
# "temperature.L"), matching lme4's ranef()/VarCorr() naming. The separator
# never occurs inside engine-generated labels except as the var/level joint,
# so a fixed-string removal is exact for engine output.
mm_coef_strip_engine_sep <- function(labels) {
  gsub(": ", "", labels, fixed = TRUE)
}

# lme4-form random-effect column labels: strip the engine separator, then
# rename logical slopes to lme4's dummy form ("x" -> "xTRUE"; the wire
# coerces logicals to numeric 0/1 so the engine emits the bare name). Used
# for both ranef() column names and cond_var postvar dimnames so the two
# stay name-aligned.
mm_re_colnames_lme4 <- function(labels, model_frame) {
  labels <- mm_coef_strip_engine_sep(labels)
  is_lgl <- vapply(model_frame, is.logical, logical(1))
  lgl_vars <- names(is_lgl)[is_lgl]
  hit <- labels %in% lgl_vars
  labels[hit] <- paste0(labels[hit], "TRUE")
  labels
}

# Normalize a freshly-constructed fit (engine-named beta/std_errors/vcov,
# engine-named random-effect columns) into the lme4 naming contract. Runs
# once per fit constructor, immediately before the class is assigned.
mm_apply_lme4_coef_naming <- function(fit) {
  map <- mm_coef_name_map(fit)
  fit$coef_map <- map
  fit$beta <- mm_coef_apply_map(fit$beta, map)
  if (!is.null(fit$std_errors)) {
    fit$std_errors <- mm_coef_apply_map(fit$std_errors, map)
  }
  V <- fit$fixed_effect_vcov
  if (!is.null(V) && is.matrix(V)) {
    if (nrow(V) != length(map$perm) || ncol(V) != length(map$perm)) {
      # beta/std_errors were just reordered; a wrong-size vcov would silently
      # misalign against them, so this must be loud.
      mm_abort(
        message = sprintf(
          "Fixed-effect covariance is %d x %d but the fit has %d coefficients.",
          nrow(V), ncol(V), length(map$perm)
        ),
        class = "mm_schema_error"
      )
    }
    atts <- attributes(V)
    V <- V[map$perm, map$perm, drop = FALSE]
    dimnames(V) <- list(map$lme4_names, map$lme4_names)
    for (a in setdiff(names(atts), c("dim", "dimnames"))) {
      attr(V, a) <- atts[[a]]
    }
    fit$fixed_effect_vcov <- V
  }
  if (!is.null(fit$random_effects)) {
    # In-place per-element assignment: lapply() would strip the mm_ranef class.
    for (g in seq_along(fit$random_effects)) {
      cols <- colnames(fit$random_effects[[g]])
      if (!is.null(cols)) {
        colnames(fit$random_effects[[g]]) <-
          mm_re_colnames_lme4(cols, fit$model_frame)
      }
    }
  }
  fit
}

# Coefficients the engine marked aliased (rank-deficiency pivot): the fit
# stores them as hard zeros so predictions stay correct, but DISPLAYING a
# zero reads as "no effect" — the opposite of "not separately estimable".
# Display surfaces show these as NA (R's lm() convention for aliased terms).
mm_aliased_coefficients <- function(fit) {
  cols <- unlist(
    fit$artifact$design_audit$fixed_effects$aliased_columns %||% list(),
    use.names = FALSE
  )
  if (!length(cols)) return(character())
  mm_coef_engine_to_lme4(as.character(cols), fit$coef_map)
}
