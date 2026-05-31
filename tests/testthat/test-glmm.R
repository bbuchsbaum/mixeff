mm_glmm_cases <- function(ids = NULL) {
  cases <- mm_lme4_parity_cases(ids = ids, model = "glmm")
  Filter(function(case) identical(case$model %||% "lmm", "glmm"), cases)
}

mm_expand_binomial_trials <- function(data) {
  required <- c("incidence", "size")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(sprintf(
      "Cannot expand binomial trials; missing column(s): %s",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }

  rows <- lapply(seq_len(nrow(data)), function(i) {
    successes <- as.integer(data$incidence[[i]])
    failures <- as.integer(data$size[[i]]) - successes
    n <- successes + failures
    if (n < 1L) {
      return(NULL)
    }
    out <- data[rep(i, n), setdiff(names(data), required), drop = FALSE]
    out$y <- c(rep.int(1L, successes), rep.int(0L, failures))
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  droplevels(out)
}

mm_glmm_case_data <- function(case) {
  data <- mm_lme4_case_data(case)
  transform <- case$data_transform %||% "none"
  switch(
    transform,
    none = data,
    expand_binomial_trials = mm_expand_binomial_trials(data),
    stop(sprintf("Unsupported GLMM parity data transform `%s`.", transform),
         call. = FALSE)
  )
}

mm_glmm_case_family <- function(case) {
  family <- case$family %||% stop("GLMM parity case is missing `family`.",
                                  call. = FALSE)
  link <- case$link %||% switch(
    family,
    binomial = "logit",
    poisson = "log",
    gamma = "log",
    stop(sprintf("GLMM parity case has unsupported family `%s`.", family),
         call. = FALSE)
  )

  switch(
    family,
    binomial = stats::binomial(link = link),
    poisson = stats::poisson(link = link),
    gamma = stats::Gamma(link = link),
    stop(sprintf("GLMM parity case has unsupported family `%s`.", family),
         call. = FALSE)
  )
}

mm_glmm_fit_pair <- function(case) {
  mm_skip_if_no_lme4()
  data <- mm_glmm_case_data(case)
  formula <- mm_lme4_case_formula(case)
  family <- mm_glmm_case_family(case)
  method <- case$method %||% "pirls_profiled"

  list(
    case = case,
    data = data,
    formula = formula,
    family = family,
    mixeff = glmm(formula, data, family = family, method = method,
                  control = mm_control(verbose = -1)),
    lme4 = suppressMessages(suppressWarnings(
      lme4::glmer(formula, data = data, family = family, nAGQ = 1L)
    ))
  )
}

mm_expect_glmm_lme4_parity <- function(case) {
  pair <- mm_glmm_fit_pair(case)
  tol <- mm_lme4_case_tolerances(case)
  fit <- pair$mixeff
  ref <- pair$lme4

  expect_s3_class(fit, "mm_glmm")
  expect_identical(nobs(fit), stats::nobs(ref),
                   info = sprintf("nobs parity failed for case `%s`", case$id))
  expect_identical(fit$family$family, case$family)
  expect_identical(fit$family$link, case$link)
  expect_identical(fit$method, case$method)

  mm_assert_parity(fixef(fit), lme4::fixef(ref),
                   case$id, "fixef", tol$fixef, "GLMM fixef")
  mm_assert_parity(fit$theta, lme4::getME(ref, "theta"),
                   case$id, "theta", tol$theta, "GLMM theta")
  mm_assert_parity(as.numeric(logLik(fit)), as.numeric(stats::logLik(ref)),
                   case$id, "logLik", tol$logLik, "GLMM logLik")
  mm_assert_parity(deviance(fit), -2 * as.numeric(stats::logLik(ref)),
                   case$id, "deviance", tol$deviance, "GLMM deviance")

  invisible(pair)
}

test_that("GLMM lme4 parity cases are declared in the manifest", {
  ids <- vapply(mm_glmm_cases(), `[[`, character(1), "id")
  expect_setequal(
    ids,
    c("cbpp_binomial_logit_profiled_pirls",
      "grouseticks_poisson_log_profiled_pirls")
  )

  for (case in mm_glmm_cases()) {
    expect_identical(case$model, "glmm")
    expect_true(case$family %in% c("binomial", "poisson", "gamma"))
    expect_identical(case$expected_status, "expected_mismatch")
    for (field in c("fixef", "theta", "logLik", "deviance")) {
      entry <- mm_parity_lookup(case$id, field)
      expect_false(is.null(entry),
                   info = sprintf("missing GLMM parity ledger entry for %s/%s",
                                  case$id, field))
      expect_identical(entry$status, "expected_mismatch")
    }
  }
})

test_that("GLMM cases match lme4 within classified ledger bounds", {
  mm_skip_if_no_lme4()

  for (case in mm_glmm_cases()) {
    mm_expect_glmm_lme4_parity(case)
  }
})
