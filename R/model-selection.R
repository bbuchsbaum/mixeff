# Fixed-effect F comparison of nested LMMs (pbkrtest::KRmodcomp /
# SATmodcomp) and lmerTest-style backward elimination (step()).

# Restriction matrix L (q x p, columns = big model's coefficients) whose null
# space is the small model's fixed-effect column space: X_big beta lies in
# span(X_small) iff L beta = 0 (pbkrtest's model2restrictionMatrix).
mm_restriction_matrix <- function(X_big, X_small, tol = 1e-8) {
  q0 <- qr(X_small)
  Q0 <- qr.Q(q0)[, seq_len(q0$rank), drop = FALSE]
  # Small model must be nested: its columns lie in span(X_big).
  q1 <- qr(X_big)
  resid <- X_small - qr.fitted(q1, X_small)
  if (max(abs(resid)) > tol * max(1, max(abs(X_small)))) {
    return(NULL)
  }
  X_perp <- X_big - Q0 %*% crossprod(Q0, X_big)
  sv <- svd(X_perp, nu = 0L)
  keep <- sv$d > tol * max(1, sv$d[1L])
  L <- t(sv$v[, keep, drop = FALSE])
  colnames(L) <- colnames(X_big)
  L
}

# compare(small, big, method = "kenward_roger" / "satterthwaite").
mm_compare_fixed_f <- function(small, big, method) {
  X1 <- stats::model.matrix(big, type = "fixed")
  X0 <- stats::model.matrix(small, type = "fixed")
  if (nrow(X1) != nrow(X0) ||
      !isTRUE(all.equal(unname(mm_response_vector(big)),
                        unname(mm_response_vector(small))))) {
    mm_abort(
      message = "F-test model comparison needs both models fitted to the same response and rows.",
      class = "mm_inference_unavailable",
      reason_code = "fixed_f_comparison_different_data",
      input = method
    )
  }
  re_text <- function(fit) {
    sort(vapply(fit$artifact$semantic_model$random_terms %||% list(),
                function(x) x$source_syntax$text %||% "", character(1)))
  }
  if (!identical(re_text(big), re_text(small))) {
    mm_abort(
      message = paste0(
        "`method = \"", method, "\"` compares fixed effects (pbkrtest's ",
        "KRmodcomp/SATmodcomp): the two models must have the same random ",
        "effects. Use method = \"lrt\" or \"bootstrap\" to compare random ",
        "structures."
      ),
      class = "mm_inference_unavailable",
      reason_code = "fixed_f_comparison_random_effects_differ",
      input = list(big = re_text(big), small = re_text(small))
    )
  }
  X1 <- X1[, names(big$beta), drop = FALSE]
  L <- mm_restriction_matrix(X1, X0)
  if (is.null(L) || !nrow(L)) {
    mm_abort(
      message = paste0(
        "The models are not nested in their fixed effects (the smaller ",
        "model's design is not contained in the larger one's, or they are ",
        "identical); an F comparison is not defined."
      ),
      class = "mm_inference_unavailable",
      reason_code = "fixed_f_comparison_not_nested",
      input = method
    )
  }
  refit_reml <- !isTRUE(big$REML)
  fit_kr <- if (refit_reml) {
    # Kenward-Roger / Satterthwaite are REML constructions; pbkrtest refits
    # an ML model by REML too. Recorded in the result (`refit_reml`).
    mm_internal_lmm(big$formula, big$model_frame, REML = TRUE,
                    weights = big$weights, offset = mm_fit_offset_arg(big),
                    control = mm_internal_control(big))
  } else {
    big
  }
  bridge <- mm_rust_fit_bridge_payload(fit_kr)
  rownames(L) <- paste0("restriction_", seq_len(nrow(L)))
  json <- tryCatch(
    mm_fixed_effect_joint_test_json(
      bridge$formula_string,
      TRUE,
      bridge$spec_data$column_order,
      bridge$spec_data$numeric_columns,
      bridge$spec_data$categorical_values,
      bridge$spec_data$categorical_levels,
      bridge$spec_data$categorical_ordered,
      bridge$weights,
      bridge$control_json,
      as.numeric(t(mm_coef_l_to_engine(L, fit_kr))),
      as.integer(nrow(L)),
      as.integer(ncol(L)),
      "small vs large model",
      rep(0, nrow(L)),
      method
    ),
    error = function(cnd) cnd
  )
  if (inherits(json, "condition")) mm_abort_from_bridge(json, method = method)
  row <- mm_json_parse_fixed_effect_inference_table(
    jsonlite::fromJSON(json, simplifyVector = FALSE)
  )$table
  row <- row[1L, , drop = FALSE]
  pull <- function(col, default = NA) {
    if (col %in% names(row)) row[[col]][[1L]] else default
  }
  kr <- pull("details", NULL)$kenward_roger %||% list()
  list(
    statistic = as.numeric(pull("statistic", NA_real_)),
    num_df = as.numeric(pull("numerator_df", NA_real_)),
    den_df = as.numeric(pull("denominator_df", pull("df", NA_real_))),
    f_scaling = as.numeric(kr$f_scaling %||% NA_real_),
    p_value = as.numeric(pull("p_value", NA_real_)),
    method = as.character(pull("method", method)),
    status = as.character(pull("status", "not_assessed")),
    reason = as.character(pull("reason", NA_character_)),
    L = L,
    refit_reml = refit_reml
  )
}

mm_apply_fixed_f_to_table <- function(table, f, method) {
  last <- nrow(table)
  table$den_df <- NA_real_
  table$statistic[last] <- f$statistic
  table$statistic_name[last] <- "F"
  table$df[last] <- f$num_df
  table$den_df[last] <- f$den_df
  table$p_value[last] <- f$p_value
  table$comparison_method[last] <- paste0(method, "_f")
  table$status[last] <- f$status
  table$reason[last] <- if (isTRUE(f$refit_reml)) {
    "F test on the REML refit of the larger model (pbkrtest semantics)"
  } else {
    f$reason
  }
  table
}

#' Stepwise model selection
#'
#' `stats::step()` is not a generic, so mixeff (like lmerTest) exports a
#' `step()` generic whose default method is `stats::step()`; the `mm_lmm`
#' method is the lmerTest-style backward elimination documented in
#' [step.mm_lmm()].
#'
#' @param object A fitted model.
#' @param ... Passed to the method (for the default, to `stats::step()`).
#' @return See the method.
#' @export
step <- function(object, ...) {
  UseMethod("step")
}

#' @rdname step
#' @export
step.default <- function(object, ...) {
  stats::step(object, ...)
}

#' Backward elimination of random and fixed effects (lmerTest's `step()`)
#'
#' `step()` for `mm_lmm` fits follows `lmerTest::step()`: it first removes
#' random-effect terms one at a time while the least significant has
#' `p > alpha.random`, then removes fixed-effect terms (respecting
#' marginality: a term is only dropped when no higher-order term contains
#' it) while the least significant has `p > alpha.fixed`, refitting after
#' every removal.
#'
#' Differences from lmerTest, by design: random terms are tested with
#' [test_random_effect()]'s boundary-corrected likelihood-ratio test (a
#' chi-bar-square mixture, so p-values are about half of lmerTest's
#' `ranova()` naive chi-square); only whole random terms are removed
#' (lmerTest's `ranova()` also tries reducing `(x | g)` to `(1 | g)` --
#' write the reduced term as its own term, e.g. `(1 | g) + (0 + x | g)`, to
#' make that reduction available); and the last random term is never
#' removed (a model without random effects is not an LMM). Fixed terms are
#' tested with [test_effect()] (`ddf` = `"Satterthwaite"` or
#' `"Kenward-Roger"`).
#'
#' @param object A fitted `mm_lmm`.
#' @param ddf Denominator degrees of freedom for the fixed-effect F tests.
#' @param alpha.random,alpha.fixed Elimination thresholds.
#' @param reduce.fixed,reduce.random Whether to eliminate fixed / random
#'   terms.
#' @param keep Character vector of fixed-effect terms never to remove.
#' @param ... Unused.
#'
#' @return An object of class `mm_step` with the elimination tables
#'   `random` and `fixed` and the final model in `$model` (also returned by
#'   `lmerTest::get_model()` when lmerTest is loaded).
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:10, each = 6)), x = rep(0:5, 10),
#'                 z = rnorm(60))
#' d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(60)
#' fit <- lmm(y ~ x + z + (1 | g) + (0 + x | g), d,
#'            control = mm_control(verbose = -1))
#' s <- step(fit)
#' s
#' s$model
#' @method step mm_lmm
#' @export
step.mm_lmm <- function(object, ddf = c("Satterthwaite", "Kenward-Roger"),
                        alpha.random = 0.1, alpha.fixed = 0.05,
                        reduce.fixed = TRUE, reduce.random = TRUE,
                        keep, ...) {
  ddf <- match.arg(ddf)
  method <- if (identical(ddf, "Satterthwaite")) "satterthwaite" else "kenward_roger"
  if (missing(keep)) keep <- character()
  fit <- object
  random_rows <- list()
  if (isTRUE(reduce.random)) {
    repeat {
      terms <- mm_random_effect_term_table(fit)
      if (nrow(terms) <= 1L) break
      tests <- lapply(terms$term_id, function(id) {
        tryCatch(test_random_effect(fit, id)$table, error = function(cnd) NULL)
      })
      ok <- !vapply(tests, is.null, logical(1))
      if (!any(ok)) break
      tab <- do.call(rbind, lapply(tests[ok], function(t) {
        t[, c("term", "statistic", "p_value", "method"), drop = FALSE]
      }))
      worst <- which.max(tab$p_value)
      if (!is.finite(tab$p_value[worst]) || tab$p_value[worst] <= alpha.random) {
        random_rows[[length(random_rows) + 1L]] <- cbind(tab, eliminated = 0L)
        break
      }
      idx <- which(terms$term_id == terms$term_id[ok][worst])
      random_rows[[length(random_rows) + 1L]] <- cbind(
        tab[worst, , drop = FALSE], eliminated = length(random_rows) + 1L
      )
      fit <- mm_internal_lmm(mm_drop_random_term_formula(fit, idx),
                             fit$model_frame, REML = isTRUE(fit$REML),
                             weights = fit$weights,
                             control = mm_internal_control(fit, keep_start = FALSE))
    }
  }
  fixed_rows <- list()
  if (isTRUE(reduce.fixed)) {
    repeat {
      cand <- setdiff(mm_droppable_terms(fit), keep)
      if (!length(cand)) break
      tab <- do.call(rbind, lapply(cand, function(term) {
        t <- test_effect(fit, term, method = method)$table
        # One-df terms come back as t tests; report F = t^2 on (1, ddf) as
        # lmerTest does.
        is_t <- identical(t$statistic_name[[1L]], "t")
        data.frame(term = term,
                   num_df = if (is_t) 1 else t$num_df %||% NA_real_,
                   den_df = t$den_df,
                   F_value = if (is_t) t$statistic^2 else t$statistic,
                   p_value = t$p_value, stringsAsFactors = FALSE)
      }))
      worst <- which.max(tab$p_value)
      if (!is.finite(tab$p_value[worst]) || tab$p_value[worst] <= alpha.fixed) {
        fixed_rows[[length(fixed_rows) + 1L]] <- cbind(tab, eliminated = 0L)
        break
      }
      fixed_rows[[length(fixed_rows) + 1L]] <- cbind(
        tab[worst, , drop = FALSE], eliminated = length(fixed_rows) + 1L
      )
      fit <- mm_internal_lmm(mm_drop_fixed_term_formula(fit, tab$term[worst]),
                             fit$model_frame, REML = isTRUE(fit$REML),
                             weights = fit$weights,
                             control = mm_internal_control(fit))
    }
  }
  bind <- function(rows) {
    if (!length(rows)) return(NULL)
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }
  structure(
    list(random = bind(random_rows), fixed = bind(fixed_rows), model = fit,
         ddf = ddf, alpha.random = alpha.random, alpha.fixed = alpha.fixed),
    class = "mm_step"
  )
}

# lmerTest::get_model(<mm_step>) support; registered into lmerTest's generic
# when lmerTest is loaded (see zzz.R), so mixeff does not mask it.
get_model.mm_step <- function(x, ...) {
  x$model
}

#' @method print mm_step
#' @export
print.mm_step <- function(x, ...) {
  cat("Backward reduced random-effect table:\n")
  if (is.null(x$random)) cat("  (none tested)\n") else print(x$random, row.names = FALSE)
  cat(sprintf("\nBackward reduced fixed-effect table (ddf: %s):\n", x$ddf))
  if (is.null(x$fixed)) cat("  (none tested)\n") else print(x$fixed, row.names = FALSE)
  cat("\nModel found:\n")
  print(x$model$formula, showEnv = FALSE)
  invisible(x)
}
