# lme4-style model-frame preparation shared by lmm() and glmm() ------------
#
# lme4 builds its model frame with stats::model.frame(): formula variables
# (including stateful transforms such as poly()/scale()/ns()) are evaluated
# on the full `data`, then `subset` and `na.action` select rows, and unused
# factor levels are dropped. mixeff mirrors that here, then hands the engine
# a plain data frame plus an "engine formula":
#
# * Fixed-part variables the engine's stateless formula language cannot
#   evaluate (factor(), poly(), scale(), ns(), cut(), relevel(), x > 0, ...)
#   are materialised as synthetic numeric/factor columns. Their names are
#   syntactic (".poly_x_2.1", ".factor_cyl") and recorded in `expansion` so
#   coefficient names shown to users match lme4's model.matrix() names
#   ("poly(x, 2)1", "factor(cyl)6"), and so prediction on new data re-uses
#   the training basis through the terms' `predvars` (as lm()/lme4 do).
# * Formula operators the engine does not parse (`^`, `%in%`, parenthesised
#   groups, `-`) are expanded by stats::terms() and sent term by term.
# * offset() terms are evaluated in R and summed with the `offset =`
#   argument. GLMMs pass the total offset to the engine; LMMs fit the
#   offset-adjusted response (y - offset) and add the offset back to fitted
#   values and predictions (lmer semantics).
#
# Formulas the engine accepts as written (no offset) are sent unchanged, so
# existing fits, explain()/audit() text and coefficient labels are identical.

# ---- formula helpers --------------------------------------------------------

mm_is_bar_call <- function(x) {
  is.call(x) && (identical(x[[1L]], quote(`|`)) || identical(x[[1L]], quote(`||`)))
}

# Remove random-effect bars from a formula RHS (a pared-down lme4::nobars()).
mm_nobars_term <- function(term) {
  if (!is.call(term)) return(term)
  if (mm_is_bar_call(term)) return(NULL)
  op <- term[[1L]]
  if (identical(op, quote(`(`))) {
    inner <- mm_nobars_term(term[[2L]])
    if (is.null(inner)) return(NULL)
    term[[2L]] <- inner
    return(term)
  }
  if (identical(op, quote(`+`))) {
    if (length(term) == 2L) {
      inner <- mm_nobars_term(term[[2L]])
      if (is.null(inner)) return(NULL)
      term[[2L]] <- inner
      return(term)
    }
    l <- mm_nobars_term(term[[2L]])
    r <- mm_nobars_term(term[[3L]])
    if (is.null(l)) return(r)
    if (is.null(r)) return(l)
    term[[2L]] <- l
    term[[3L]] <- r
    return(term)
  }
  term
}

# The fixed-effects-only version of an lme4 formula (response kept).
mm_fixed_only_formula <- function(formula) {
  rhs <- mm_nobars_term(formula[[length(formula)]])
  if (is.null(rhs)) rhs <- 1
  out <- formula
  out[[length(out)]] <- rhs
  out
}

# lme4's expandSlash(): `a/b` -> list(b:a, a); `a/b/c` -> list(c:(b:a), b:a, a).
mm_expand_slash_group <- function(g) {
  if (is.call(g) && identical(g[[1L]], quote(`/`))) {
    left <- mm_expand_slash_group(g[[2L]])
    return(c(list(call(":", g[[3L]], left[[1L]])), left))
  }
  list(g)
}

# lme4-style grouping-factor specs for every random term of `formula`:
# list(label = "b:a", vars = c("b", "a"), lhs = <bar lhs>, double = <is ||>).
# Labels and variable order follow lme4 (names(ranef(lmer_fit))).
mm_lme4_group_specs <- function(formula) {
  if (!inherits(formula, "formula")) return(list())
  bars <- mm_find_bars(formula[[length(formula)]])
  out <- list()
  for (b in bars) {
    for (g in mm_expand_slash_group(b[[3L]])) {
      out[[length(out) + 1L]] <- list(
        label = deparse1(g),
        vars = all.vars(g),
        lhs = b[[2L]],
        double = identical(b[[1L]], quote(`||`))
      )
    }
  }
  out
}

# ---- engine acceptance probe ----------------------------------------------

mm_engine_parses <- function(text) {
  ok <- tryCatch(
    {
      .Call(wrap__mm_parse_formula, text)
      TRUE
    },
    error = function(cnd) FALSE
  )
  isTRUE(ok)
}

# Synthetic, engine-safe column name for an expanded variable.
mm_synthetic_base <- function(label, taken) {
  base <- gsub("[^A-Za-z0-9]+", "_", label)
  base <- gsub("^_+|_+$", "", base)
  if (!nzchar(base)) base <- "term"
  base <- paste0(".", base)
  out <- base
  k <- 1L
  # A base must not equal or prefix-collide with any other name in play, so
  # the "var + level" coefficient names decode unambiguously.
  clash <- function(nm) any(startsWith(taken, nm) | startsWith(nm, taken))
  while (clash(out)) {
    k <- k + 1L
    out <- paste0(base, "_", k)
  }
  out
}

# Suffixes model.matrix() appends to a (numeric, possibly matrix-valued)
# variable's name: "" for a plain vector / unnamed 1-column matrix, the column
# names otherwise, or 1..k when a multi-column matrix has none.
mm_numeric_component_suffixes <- function(value) {
  if (is.null(dim(value))) return("")
  k <- ncol(value)
  cn <- colnames(value)
  if (is.null(cn)) {
    if (k == 1L) "" else as.character(seq_len(k))
  } else {
    cn
  }
}

# ---- expansion ----------------------------------------------------------------

# Analyse the fixed part of `formula`. Returns NULL when the engine accepts the
# formula as written and there is nothing to evaluate in R; otherwise a list
# describing every fixed-part variable and the term structure.
mm_plan_formula_expansion <- function(formula, data) {
  fixed <- mm_fixed_only_formula(formula)
  rhs_terms <- stats::delete.response(stats::terms(fixed, data = data))
  has_offset <- !is.null(attr(rhs_terms, "offset"))
  # The response is the engine's business (it evaluates log(y) etc. itself;
  # glmm() resolves cbind() responses separately), so probe the RHS only.
  probe <- paste(".mm_probe_response ~", deparse1(formula[[length(formula)]]))
  if (!has_offset && mm_engine_parses(probe)) {
    return(NULL)
  }
  vars <- as.list(attr(rhs_terms, "variables"))[-1L]
  labels <- vapply(vars, deparse1, character(1))
  offset_idx <- attr(rhs_terms, "offset") %||% integer()
  kinds <- vapply(seq_along(vars), function(i) {
    v <- vars[[i]]
    if (i %in% offset_idx) return("offset")
    if (is.name(v)) return("bare")
    if (mm_engine_parses(paste(".mm_probe_response ~", labels[[i]]))) {
      return("passthrough")
    }
    "expand"
  }, character(1))
  list(terms = rhs_terms, labels = labels, kinds = kinds)
}

# Evaluate the expanded / offset variables on the full data (model.frame()
# semantics: stateful transforms see every row; subset/NA come afterwards).
mm_evaluate_expansion <- function(plan, data) {
  # A fit's model frame carries the training basis of its stateful terms
  # (attr "mm_predvars"), so internal refits on it (REML->ML, drop1, update,
  # refit) reproduce exactly the same columns rather than re-deriving
  # data-dependent knots/centres from the retained rows.
  stored <- attr(data, "mm_predvars")
  if (is.list(stored) && length(stored)) {
    pv <- attr(plan$terms, "variables")
    for (i in seq_along(plan$labels)) {
      if (plan$labels[[i]] %in% names(stored)) {
        pv[[i + 1L]] <- stored[[plan$labels[[i]]]]
      }
    }
    attr(plan$terms, "predvars") <- pv
  }
  mf <- stats::model.frame(plan$terms, data = data, na.action = stats::na.pass)
  plan$terms <- stats::terms(mf)
  plan$values <- lapply(seq_along(plan$labels), function(i) {
    if (plan$kinds[[i]] %in% c("expand", "offset")) mf[[i]] else NULL
  })
  plan
}

# Build the synthetic columns (on the selected rows) and the engine formula.
mm_finish_expansion <- function(plan, formula, rows, data_names) {
  # Only the formula's own columns are kept next to the synthetic ones, so
  # those are the names a synthetic name must not collide with.
  taken <- data_names
  records <- list()
  spell <- list()
  for (i in seq_along(plan$labels)) {
    label <- plan$labels[[i]]
    kind <- plan$kinds[[i]]
    if (kind %in% c("bare", "passthrough")) {
      spell[[label]] <- label
      next
    }
    if (identical(kind, "offset")) next
    value <- plan$values[[i]]
    base <- mm_synthetic_base(label, taken)
    if (is.character(value)) value <- factor(value)
    if (is.factor(value) || is.logical(value)) {
      value <- if (is.logical(value)) {
        factor(value[rows], levels = c(FALSE, TRUE))
      } else {
        value[rows, drop = TRUE]
      }
      cols <- stats::setNames(list(value), base)
      rec <- list(label = label, kind = "factor", columns = base,
                  components = stats::setNames(label, base))
    } else if (is.numeric(value)) {
      sfx <- mm_numeric_component_suffixes(value)
      mat <- if (is.null(dim(value))) matrix(value, ncol = 1L) else value
      mat <- mat[rows, , drop = FALSE]
      names_k <- if (ncol(mat) == 1L) base else paste0(base, ".", seq_len(ncol(mat)))
      cols <- stats::setNames(
        lapply(seq_len(ncol(mat)), function(k) as.numeric(mat[, k])),
        names_k
      )
      rec <- list(label = label, kind = "numeric", columns = names_k,
                  components = stats::setNames(paste0(label, sfx), names_k))
    } else {
      mm_abort(
        message = sprintf(
          paste0("Formula term `%s` evaluates to an unsupported type (%s); ",
                 "fixed-effect terms must be numeric, logical, or factor."),
          label, paste(class(value), collapse = "/")
        ),
        class = "mm_formula_error",
        input = label
      )
    }
    rec$values <- cols
    records[[label]] <- rec
    spell[[label]] <- rec$columns
    taken <- c(taken, rec$columns)
  }

  # Engine fixed terms: each R term expanded into the product of its
  # components' engine spellings (the engine has no parenthesised groups).
  tt <- plan$terms
  fac <- attr(tt, "factors")
  term_labels <- attr(tt, "term.labels")
  term_map <- list()
  engine_terms <- character()
  for (j in seq_along(term_labels)) {
    comp <- rownames(fac)[fac[, j] > 0]
    parts <- lapply(comp, function(v) spell[[v]])
    combos <- Reduce(function(acc, p) {
      as.vector(outer(acc, p, paste, sep = ":"))
    }, parts[-1L], parts[[1L]])
    term_map[[term_labels[[j]]]] <- combos
    engine_terms <- c(engine_terms, combos)
  }
  intercept <- attr(tt, "intercept") == 1L
  fixed_rhs <- c(if (!intercept) "0", engine_terms)
  if (!length(fixed_rhs)) fixed_rhs <- "1"
  bars <- vapply(mm_find_bars(formula[[length(formula)]]),
                 function(b) paste0("(", deparse1(b), ")"), character(1))
  rhs <- paste(c(fixed_rhs, bars), collapse = " + ")
  lhs <- deparse1(formula[[2L]])

  offset_idx <- which(plan$kinds == "offset")
  offset <- NULL
  if (length(offset_idx)) {
    offset <- Reduce(`+`, lapply(offset_idx, function(i) {
      as.numeric(plan$values[[i]])[rows]
    }))
  }

  list(
    terms = tt,
    records = records,
    term_map = term_map,
    engine_rhs = rhs,
    engine_lhs = lhs,
    offset = offset,
    formula_env = environment(formula)
  )
}

# ---- the shared data preparation -----------------------------------------------

# Prepare `data` for a fit: evaluate formula transforms, apply `subset` and
# `na.action`, drop unused factor levels, coerce grouping columns, and build
# the engine formula. Returns
#   data         the model frame (raw formula variables + synthetic columns)
#   formula_engine  the formula the engine fits
#   expansion    NULL or the expansion record (see mm_finish_expansion())
#   weights, offset, offset_arg   aligned with the rows of `data`
#   na_action    the na.action record (class "omit"/"exclude") or NULL
mm_prepare_model_data <- function(formula, data, subset_expr, na.action,
                                  weights, offset_arg, enclos, verbose,
                                  lmm = TRUE) {
  if (!is.data.frame(data)) {
    mm_abort(message = "`data` must be a data.frame.", class = "mm_data_error",
             input = data)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    mm_abort(
      message = "`formula` must be a two-sided R formula (lhs ~ rhs).",
      class = "mm_formula_error",
      formula = formula
    )
  }
  n0 <- nrow(data)
  if (!is.null(offset_arg)) {
    if (!is.numeric(offset_arg) || length(offset_arg) != n0 ||
        any(!is.finite(offset_arg))) {
      mm_abort(
        message = sprintf(
          "`offset` must be a finite numeric vector with one value per row of `data` (%d).",
          n0
        ),
        class = "mm_arg_error",
        input = offset_arg
      )
    }
    offset_arg <- as.numeric(offset_arg)
  }

  plan <- mm_plan_formula_expansion(formula, data)
  if (!is.null(plan)) plan <- mm_evaluate_expansion(plan, data)

  # ---- subset ----
  rows <- seq_len(n0)
  if (!identical(subset_expr, quote(NULL)) && !is.null(subset_expr)) {
    keep <- eval(subset_expr, data, enclos)
    if (is.logical(keep)) {
      if (length(keep) != n0) {
        mm_abort(
          message = sprintf("`subset` must have one value per row of `data` (%d).", n0),
          class = "mm_arg_error", input = keep
        )
      }
      keep[is.na(keep)] <- FALSE
      rows <- which(keep)
    } else if (is.numeric(keep)) {
      rows <- as.integer(keep)
      rows <- rows[!is.na(rows)]
      if (any(rows < 0L)) rows <- seq_len(n0)[rows]
    } else {
      mm_abort(message = "`subset` must evaluate to a logical or integer index.",
               class = "mm_arg_error", input = keep)
    }
  }

  # ---- missing values (raw formula variables + evaluated transforms) ----
  raw_vars <- intersect(all.vars(formula), names(data))
  check <- data[rows, raw_vars, drop = FALSE]
  if (!is.null(plan)) {
    for (i in which(plan$kinds %in% c("expand", "offset"))) {
      v <- plan$values[[i]]
      check[[plan$labels[[i]]]] <- if (is.null(dim(v))) v[rows] else v[rows, , drop = FALSE]
    }
  }
  na_action <- NULL
  if (is.null(na.action)) {
    mm_check_no_na(check, names(check))
  } else {
    if (!is.function(na.action)) na.action <- match.fun(na.action)
    cleaned <- na.action(check)
    na_action <- attr(cleaned, "na.action")
    if (!is.null(na_action)) {
      rows <- rows[-as.integer(na_action)]
    }
  }

  out <- data[rows, raw_vars, drop = FALSE]
  if (!is.null(weights)) weights <- weights[rows]
  if (!is.null(offset_arg)) offset_arg <- offset_arg[rows]

  # ---- drop unused levels (lme4: model.frame(drop.unused.levels = TRUE)) ----
  response_vars <- all.vars(formula[[2L]])
  for (nm in setdiff(names(out), response_vars)) {
    x <- out[[nm]]
    if (is.factor(x) && nlevels(x) > length(unique(x[!is.na(x)]))) {
      ctr <- attr(x, "contrasts")
      x <- x[, drop = TRUE]
      if (is.character(ctr)) attr(x, "contrasts") <- ctr
      out[[nm]] <- x
    }
  }

  # ---- grouping factors: categorical, unordered, announced coercion ----
  out <- mm_apply_grouping_coercion(formula, out, verbose)
  fixed_vars <- all.vars(mm_fixed_only_formula(formula)[[3L]])
  for (g in setdiff(mm_formula_grouping_vars(formula), fixed_vars)) {
    # Grouping only needs levels; an ordered grouping column would otherwise
    # be shipped with a dense k x k contr.poly basis the engine never uses.
    if (is.ordered(out[[g]])) out[[g]] <- factor(out[[g]], ordered = FALSE)
  }

  # ---- logical predictors are factors (levels FALSE/TRUE) in model.matrix ----
  grouping <- mm_formula_grouping_vars(formula)
  bare <- mm_bare_term_vars(formula)
  for (nm in setdiff(intersect(names(out), bare), c(response_vars, grouping))) {
    if (is.logical(out[[nm]])) {
      out[[nm]] <- factor(out[[nm]], levels = c(FALSE, TRUE))
    }
  }
  out <- mm_factor_character_columns_except(out, response_vars)

  # ---- expansion columns + engine formula ----
  expansion <- NULL
  formula_engine <- formula
  offset <- offset_arg
  if (!is.null(plan)) {
    expansion <- mm_finish_expansion(plan, formula, rows, raw_vars)
    for (rec in expansion$records) {
      for (nm in names(rec$values)) out[[nm]] <- rec$values[[nm]]
      rec$values <- NULL
    }
    for (lab in names(expansion$records)) expansion$records[[lab]]$values <- NULL
    if (!is.null(expansion$offset)) {
      offset <- if (is.null(offset)) expansion$offset else offset + expansion$offset
    }
    lhs <- expansion$engine_lhs
    formula_engine <- stats::as.formula(
      paste(lhs, "~", expansion$engine_rhs),
      env = environment(formula)
    )
  }

  # LMM offsets: fit the offset-adjusted response; fitted values and
  # predictions add the offset back (lmer semantics).
  if (isTRUE(lmm) && !is.null(offset)) {
    y <- eval(formula[[2L]], out, environment(formula) %||% enclos)
    if (!is.numeric(y) || length(y) != nrow(out)) {
      mm_abort(
        message = "An LMM offset requires a numeric response with one value per row.",
        class = "mm_data_error"
      )
    }
    out[[".mm_offset_response"]] <- as.numeric(y) - offset
    formula_engine[[2L]] <- as.name(".mm_offset_response")
  }

  if (!is.null(na_action)) attr(out, "na.action") <- na_action
  if (!is.null(expansion)) {
    pv <- attr(expansion$terms, "predvars")
    if (!is.null(pv)) {
      idx <- which(plan$kinds %in% c("expand", "offset"))
      attr(out, "mm_predvars") <- stats::setNames(
        lapply(idx, function(i) pv[[i + 1L]]), plan$labels[idx]
      )
    }
  }
  list(
    data = out,
    formula_engine = formula_engine,
    expansion = expansion,
    weights = weights,
    offset = offset,
    offset_arg = offset_arg,
    na_action = na_action
  )
}

# Variables that enter the model as bare columns (fixed terms or random-term
# left-hand sides), as opposed to only inside a transform such as I(l * x).
mm_bare_term_vars <- function(formula) {
  bare_of <- function(f) {
    tt <- tryCatch(stats::terms(f), error = function(e) NULL)
    if (is.null(tt)) return(character())
    vars <- as.list(attr(tt, "variables"))[-1L]
    vapply(Filter(is.name, vars), as.character, character(1))
  }
  fixed <- mm_fixed_only_formula(formula)
  fixed[[2L]] <- NULL
  out <- bare_of(fixed)
  for (b in mm_find_bars(formula[[length(formula)]])) {
    out <- c(out, bare_of(stats::as.formula(call("~", b[[2L]]))))
  }
  unique(out)
}

mm_factor_character_columns_except <- function(data, except) {
  for (nm in setdiff(names(data), except)) {
    if (is.character(data[[nm]])) data[[nm]] <- factor(data[[nm]])
  }
  data
}

# ---- fit accessors ---------------------------------------------------------------

# Formula the engine fits (synthetic columns, expanded terms, offset-adjusted
# response). Equal to the user formula for fits that needed no expansion.
mm_engine_formula <- function(fit) {
  fit$engine_formula %||% fit$formula
}

# Columns of the model frame the engine sees (raw variables used only inside
# expanded transforms, e.g. a Date column under as.numeric(), stay R-side).
mm_engine_frame <- function(fit) {
  mf <- fit$model_frame
  vars <- intersect(all.vars(mm_engine_formula(fit)), names(mf))
  mf[, vars, drop = FALSE]
}

# User-facing fixed-effect terms (with `predvars` for stateful transforms).
mm_user_fixed_terms <- function(fit) {
  fit$expansion$terms %||%
    stats::delete.response(stats::terms(mm_fixed_formula(fit)))
}

# Translate model.matrix() column names built on the synthetic columns into
# the names model.matrix() gives the user's formula ("poly(x, 2)1:f2").
mm_expansion_user_names <- function(names, expansion) {
  if (is.null(expansion) || !length(expansion$records)) return(names)
  numeric_map <- character()
  factor_bases <- character()
  for (rec in expansion$records) {
    if (identical(rec$kind, "numeric")) {
      numeric_map <- c(numeric_map, rec$components)
    } else {
      factor_bases <- c(factor_bases, stats::setNames(rec$label, rec$columns))
    }
  }
  factor_order <- names(factor_bases)[order(nchar(names(factor_bases)),
                                            decreasing = TRUE)]
  translate <- function(comp) {
    if (comp %in% names(numeric_map)) return(numeric_map[[comp]])
    for (b in factor_order) {
      if (startsWith(comp, b)) {
        return(paste0(factor_bases[[b]], substring(comp, nchar(b) + 1L)))
      }
    }
    comp
  }
  vapply(names, function(nm) {
    if (identical(nm, "(Intercept)")) return(nm)
    parts <- strsplit(nm, ":", fixed = TRUE)[[1L]]
    paste(vapply(parts, translate, character(1)), collapse = ":")
  }, character(1), USE.NAMES = FALSE)
}

# Add the synthetic expansion columns to arbitrary data (newdata, a reference
# grid), evaluating the user's transforms with the TRAINING basis (`predvars`).
# Factor-valued transforms are re-levelled to the training levels.
mm_expand_newdata <- function(fit, data) {
  ex <- fit$expansion
  if (is.null(ex) || !length(ex$records)) return(data)
  needed <- unlist(lapply(ex$records, `[[`, "columns"), use.names = FALSE)
  if (all(needed %in% names(data))) return(data)
  mf <- stats::model.frame(ex$terms, data = data, na.action = stats::na.pass)
  for (rec in ex$records) {
    value <- mf[[rec$label]]
    if (identical(rec$kind, "factor")) {
      train <- fit$model_frame[[rec$columns]]
      lv <- levels(train)
      chr <- as.character(value)
      bad <- !is.na(chr) & !(chr %in% lv)
      if (any(bad)) {
        mm_abort(
          message = sprintf(
            "Term `%s` has value(s) in `newdata` not seen when fitting: %s.",
            rec$label, paste(unique(chr[bad]), collapse = ", ")
          ),
          class = "mm_data_error",
          input = unique(chr[bad])
        )
      }
      data[[rec$columns]] <- factor(chr, levels = lv, ordered = is.ordered(train))
    } else {
      mat <- if (is.null(dim(value))) matrix(value, ncol = 1L) else value
      for (k in seq_along(rec$columns)) {
        data[[rec$columns[[k]]]] <- as.numeric(mat[, k])
      }
    }
  }
  data
}

# Offset contribution of the formula's offset() terms for `data`.
mm_formula_offset <- function(fit, data) {
  ex <- fit$expansion
  if (is.null(ex)) return(NULL)
  tt <- ex$terms
  idx <- attr(tt, "offset")
  if (is.null(idx)) return(NULL)
  vars <- attr(tt, "variables")
  pv <- attr(tt, "predvars") %||% vars
  env <- environment(tt) %||% ex$formula_env %||% globalenv()
  Reduce(`+`, lapply(idx, function(i) {
    as.numeric(eval(pv[[i + 1L]], data, env))
  }))
}

# Pad a per-observation vector/matrix for na.exclude fits (lm/lme4 semantics).
mm_napredict <- function(fit, x) {
  na <- fit$na.action
  if (is.null(na)) return(x)
  stats::napredict(na, x)
}

mm_naresid <- function(fit, x) {
  na <- fit$na.action
  if (is.null(na)) return(x)
  stats::naresid(na, x)
}

# ---- lme4 grouping-factor labels -------------------------------------------------

# The engine labels an interaction / nested grouping factor "a & b" with
# levels "x_y"; lme4 names it after the grouping expression ("a:b", or "b:a"
# for the nested term generated by `a/b`) with levels "x:y" in that variable
# order. Build the map once per fit from the user formula and the model frame.
mm_lme4_group_map <- function(engine_groups, formula, model_frame) {
  specs <- mm_lme4_group_specs(formula)
  out <- list()
  for (eg in unique(engine_groups)) {
    vars <- strsplit(eg, " & ", fixed = TRUE)[[1L]]
    if (length(vars) < 2L || !all(vars %in% names(model_frame))) next
    hit <- NULL
    for (s in specs) {
      if (setequal(s$vars, vars) && length(s$vars) == length(vars)) {
        hit <- s
        break
      }
    }
    order_vars <- if (is.null(hit)) vars else hit$vars
    label <- if (is.null(hit)) paste(vars, collapse = ":") else hit$label
    cols <- model_frame[vars]
    engine_lab <- do.call(paste, c(lapply(cols, as.character), sep = "_"))
    lme4_lab <- do.call(paste, c(lapply(model_frame[order_vars], as.character),
                                 sep = ":"))
    keep <- !duplicated(engine_lab)
    ord_keys <- lapply(model_frame[order_vars][keep, , drop = FALSE],
                       function(x) as.integer(factor(x)))
    rank <- do.call(order, unname(ord_keys))
    out[[eg]] <- list(
      label = label,
      levels = stats::setNames(lme4_lab[keep], engine_lab[keep]),
      order = lme4_lab[keep][rank]
    )
  }
  out
}

mm_rename_ranef_groups <- function(ranef_list, group_map) {
  if (!length(group_map)) return(ranef_list)
  cls <- class(ranef_list)
  nms <- names(ranef_list)
  for (i in seq_along(ranef_list)) {
    m <- group_map[[nms[[i]]]]
    if (is.null(m)) next
    df <- ranef_list[[i]]
    rn <- rownames(df)
    mapped <- unname(m$levels[rn])
    mapped[is.na(mapped)] <- rn[is.na(mapped)]
    rownames(df) <- mapped
    ord <- intersect(m$order, mapped)
    if (length(ord) == nrow(df)) df <- df[ord, , drop = FALSE]
    ranef_list[[i]] <- df
    nms[[i]] <- m$label
  }
  names(ranef_list) <- nms
  class(ranef_list) <- cls
  ranef_list
}

mm_rename_postvar_groups <- function(postvars, group_map) {
  if (!length(group_map)) return(postvars)
  nms <- names(postvars)
  for (i in seq_along(postvars)) {
    m <- group_map[[nms[[i]]]]
    if (is.null(m)) next
    arr <- postvars[[i]]
    lv <- dimnames(arr)[[3L]]
    mapped <- unname(m$levels[lv])
    mapped[is.na(mapped)] <- lv[is.na(mapped)]
    dimnames(arr)[[3L]] <- mapped
    postvars[[i]] <- arr
    nms[[i]] <- m$label
  }
  names(postvars) <- nms
  postvars
}

# Apply lme4 grouping-factor labels to a freshly constructed fit.
mm_apply_lme4_group_labels <- function(fit) {
  engine_groups <- names(fit$random_effects %||% list())
  for (comp in fit$varcorr$components_raw %||% list()) {
    engine_groups <- c(engine_groups, as.character(comp$group))
  }
  engine_groups <- unique(gsub(":", " & ", engine_groups, fixed = TRUE))
  map <- mm_lme4_group_map(engine_groups, fit$formula, fit$model_frame)
  fit$group_map <- map
  if (!length(map)) return(fit)
  fit$random_effects <- mm_rename_ranef_groups(fit$random_effects, map)
  relabel <- function(g) {
    key <- gsub(":", " & ", g, fixed = TRUE)
    hit <- key %in% names(map)
    g[hit] <- vapply(map[key[hit]], `[[`, character(1), "label")
    g
  }
  if (!is.null(fit$varcorr)) {
    vc <- fit$varcorr
    vc$table$group <- relabel(vc$table$group)
    vc$components_raw <- lapply(vc$components_raw, function(comp) {
      comp$group <- relabel(comp$group)
      comp
    })
    fit$varcorr <- vc
  }
  fit
}

# ---- random-effect contribution from stored BLUPs -----------------------------

# Random-effect part of the linear predictor for `data`, for the random terms
# of `re_formula` (a formula whose bars name a subset of the model's terms),
# reconstructed from the stored conditional modes like lme4's
# mkNewReTrms(): each term's model matrix (its bar's LHS evaluated on
# `data`) times the BLUP rows of each observation's group level. Unseen
# levels contribute zero when allow_new_levels = TRUE, else refuse.
mm_re_eta <- function(fit, data, re_formula, allow_new_levels,
                      missing_class = "mm_arg_error") {
  specs <- mm_lme4_group_specs(re_formula)
  eta <- numeric(nrow(data))
  ranefs <- fit$random_effects %||% list()
  xlev <- lapply(Filter(is.factor, fit$model_frame), levels)
  for (s in specs) {
    blup <- ranefs[[s$label]]
    if (is.null(blup)) {
      mm_abort(
        message = sprintf(
          "Grouping factor `%s` in `re.form` is not a grouping factor of the fitted model (available: %s).",
          s$label, paste(names(ranefs), collapse = ", ")
        ),
        class = missing_class,
        input = s$label
      )
    }
    missing_vars <- setdiff(c(s$vars, all.vars(s$lhs)), names(data))
    if (length(missing_vars)) {
      mm_abort(
        message = sprintf(
          "`newdata` is missing variable(s) required by the random-effects formula: %s.",
          paste(missing_vars, collapse = ", ")
        ),
        class = "mm_data_error",
        missing = missing_vars
      )
    }
    lhs_formula <- stats::as.formula(call("~", s$lhs), env = environment(fit$formula))
    lhs_vars <- intersect(all.vars(s$lhs), names(xlev))
    mf <- stats::model.frame(lhs_formula, data = data, na.action = stats::na.pass,
                             xlev = xlev[lhs_vars])
    Z <- stats::model.matrix(lhs_formula, mf)
    colnames(Z) <- mm_re_colnames_lme4(colnames(Z), fit$model_frame)
    absent <- setdiff(colnames(Z), names(blup))
    if (length(absent)) {
      mm_abort(
        message = sprintf(
          "Random effect(s) %s for grouping factor `%s` in `re.form` are not in the fitted model.",
          paste(sprintf("`%s`", absent), collapse = ", "), s$label
        ),
        class = "mm_arg_error",
        input = absent
      )
    }
    group <- do.call(paste, c(lapply(data[s$vars], as.character), sep = ":"))
    idx <- match(group, rownames(blup))
    unseen <- is.na(idx)
    if (any(unseen) && !isTRUE(allow_new_levels)) {
      mm_abort(
        message = sprintf(
          paste0("`newdata` has %d row(s) with grouping level(s) of `%s` not ",
                 "seen during fitting. Set allow.new.levels = TRUE to predict ",
                 "them at the population mean (zero random effect)."),
          sum(unseen), s$label
        ),
        class = "mm_inference_unavailable",
        group = s$label
      )
    }
    for (col in colnames(Z)) {
      modes <- blup[[col]][idx]
      modes[is.na(modes)] <- 0
      eta <- eta + Z[, col] * modes
    }
  }
  eta
}

# ---- term labels for expanded fits ---------------------------------------------

# Order-insensitive key for an interaction label ("b:a" == "a:b").
mm_term_key <- function(x) {
  vapply(strsplit(x, ":", fixed = TRUE), function(p) {
    paste(sort(trimws(p)), collapse = ":")
  }, character(1))
}

# Map engine fixed-term labels (synthetic columns) back to the user's terms:
# a term that is the whole of one user term gets that term's label
# ("factor(f)"); a piece of a multi-column term gets its column name
# ("poly(x, 2)1"). Labels of fits without expansion pass through.
mm_user_term_label <- function(fit, terms) {
  ex <- fit$expansion
  if (is.null(ex) || !length(terms)) return(terms)
  keys <- mm_term_key(terms)
  out <- mm_expansion_user_names(terms, ex)
  for (lab in names(ex$term_map)) {
    combos <- ex$term_map[[lab]]
    if (length(combos) == 1L) {
      out[keys == mm_term_key(combos)] <- lab
    }
  }
  out
}

# Fixed terms drop1() offers: the user's terms for expanded fits.
mm_drop1_terms <- function(fit) {
  if (!is.null(fit$expansion)) {
    return(attr(fit$expansion$terms, "term.labels"))
  }
  setdiff(mm_fixed_effect_terms(fit), "1")
}

# Display formula for a drop1() row.
mm_drop_formula_label <- function(fit, term, reduced_formula) {
  if (is.null(fit$expansion)) return(deparse1(reduced_formula))
  deparse1(stats::update.formula(
    fit$formula, stats::as.formula(paste(". ~ . -", term))
  ))
}
