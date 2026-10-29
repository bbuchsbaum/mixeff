# Ordinary chi-square LRTs here compare fixed effects on the same likelihood.
mm_assert_glmm_lrt_pair <- function(null, alternative) {
  refuse <- function(reason, code) {
    mm_abort(paste0("GLMM likelihood-ratio comparison requires ", reason, "."),
             class = "mm_inference_unavailable", reason_code = code)
  }
  same_values <- function(x, y) {
    isTRUE(all.equal(unname(x), unname(y), tolerance = 0, check.attributes = FALSE))
  }
  n <- null$nobs
  if (length(n) != 1L || is.na(n) || !identical(n, alternative$nobs) ||
      !identical(rownames(null$model_frame), rownames(alternative$model_frame))) {
    refuse("the same observed rows", "glmm_lrt_rows_differ")
  }
  y0 <- mm_response_vector(null)
  y1 <- mm_response_vector(alternative)
  if (is.null(y0) || is.null(y1) || length(y0) != n || length(y1) != n ||
      !same_values(y0, y1)) {
    refuse("the same response values", "glmm_lrt_response_differs")
  }
  same_default <- function(x, y, default) {
    x <- x %||% rep(default, n)
    y <- y %||% rep(default, n)
    length(x) == n && length(y) == n && same_values(x, y)
  }
  if (!same_default(null$weights, alternative$weights, 1)) {
    refuse("the same prior weights and binomial trial counts", "glmm_lrt_weights_differ")
  }
  if (!same_default(null$offset, alternative$offset, 0)) {
    refuse("the same offsets", "glmm_lrt_offsets_differ")
  }
  f0 <- null$family
  f1 <- alternative$family
  if (is.null(f0$family) || is.null(f0$link) ||
      !identical(f0$family, f1$family) || !identical(f0$link, f1$link)) {
    refuse("the same family and link", "glmm_lrt_family_differs")
  }
  if (identical(f0$family, "negative_binomial")) {
    if (is.null(f0$nb_theta_estimated) || is.null(f1$nb_theta_estimated) ||
        !identical(f0$nb_theta_estimated, f1$nb_theta_estimated) ||
        (!isTRUE(f0$nb_theta_estimated) &&
         (is.null(f0$nb_theta) || is.null(f1$nb_theta) ||
          !same_values(f0$nb_theta, f1$nb_theta)))) {
      refuse("the same negative-binomial theta mode and fixed theta", "glmm_lrt_nb_theta_differs")
    }
  }
  if (mm_glmm_free_dispersion(null)) {
    p0 <- null$dispersion_policy
    p1 <- alternative$dispersion_policy
    if (is.null(p0) || is.null(p1) || !identical(p0, p1)) {
      refuse("the same known fit-time dispersion policy", "glmm_lrt_dispersion_differs")
    }
  }
  meta <- function(fit) {
    m <- fit$artifact$glmm_fit_metadata
    list(method = m$effective_method %||% m$estimation_method,
         objective = m$objective_definition, constants = m$response_constants,
         nagq = m$n_agq)
  }
  m0 <- meta(null)
  m1 <- meta(alternative)
  if (any(vapply(m0, is.null, logical(1))) ||
      any(vapply(m1, is.null, logical(1))) || !identical(m0, m1)) {
    refuse("the same known effective estimator, objective, and nAGQ", "glmm_lrt_estimator_differs")
  }
  random_text <- function(fit) {
    terms <- fit$artifact$semantic_model$random_terms
    if (!length(terms)) return(NULL)
    text <- vapply(terms, function(term) term$source_syntax$text %||% "", character(1))
    if (any(!nzchar(text))) return(NULL)
    sort(text)
  }
  r0 <- random_text(null)
  r1 <- random_text(alternative)
  if (is.null(r0) || !identical(r0, r1)) {
    refuse("the same random-effects structure", "glmm_lrt_random_effects_differ")
  }
  vars <- tryCatch(unique(unlist(lapply(r0, function(x) all.vars(stats::as.formula(paste("~", x)))))),
                   error = function(e) NULL)
  if (is.null(vars) || !all(vars %in% names(null$model_frame)) ||
      !all(vars %in% names(alternative$model_frame)) ||
      !all(vapply(vars, function(v) same_values(null$model_frame[[v]], alternative$model_frame[[v]]), logical(1)))) {
    refuse("the same random-effects design values", "glmm_lrt_random_design_differs")
  }
  X0 <- tryCatch(stats::model.matrix(null, type = "fixed"), error = function(e) NULL)
  X1 <- tryCatch(stats::model.matrix(alternative, type = "fixed"), error = function(e) NULL)
  if (is.null(X0) || is.null(X1) || nrow(X0) != n || nrow(X1) != n ||
      any(!is.finite(X0)) || any(!is.finite(X1))) {
    refuse("available finite fixed-effect design matrices", "glmm_lrt_design_unavailable")
  }
  # Normalize columns before SVD so the rank test is relative to each column's
  # scale. Project the smaller space onto the larger, independent of basis order.
  basis <- function(X, tol = 1e-8) {
    norms <- sqrt(colSums(X * X))
    X <- X[, norms > 0, drop = FALSE]
    norms <- norms[norms > 0]
    if (!ncol(X)) return(matrix(numeric(), nrow(X), 0L))
    s <- svd(sweep(X, 2L, norms, "/"), nu = min(dim(X)), nv = 0L)
    rank <- sum(s$d > tol * s$d[1L])
    s$u[, seq_len(rank), drop = FALSE]
  }
  U0 <- basis(X0)
  U1 <- basis(X1)
  if (ncol(U1) <= ncol(U0) ||
      (ncol(U0) && max(abs(U0 - U1 %*% crossprod(U1, U0))) > 1e-8)) {
    refuse("strictly nested fixed-effect column spaces with increasing rank", "glmm_lrt_not_nested")
  }
  invisible(TRUE)
}
