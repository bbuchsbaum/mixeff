#!/usr/bin/env Rscript

# Fit mixeff LMMs on classic lme4 datasets and compare them with lme4::lmer().
# Run from the package root with:
#   Rscript examples/lme4-standard-datasets.R

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("This example requires the lme4 package.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(mixeff)
})

benchmark_reps <- as.integer(Sys.getenv("MIXEFF_EXAMPLE_BENCHMARK_REPS", "3"))
if (is.na(benchmark_reps) || benchmark_reps < 1L) {
  benchmark_reps <- 3L
}

cases <- list(
  list(
    id = "sleepstudy_random_intercept",
    dataset = "sleepstudy",
    formula = Reaction ~ Days + (1 | Subject),
    reml = TRUE
  ),
  list(
    id = "sleepstudy_random_intercept_slope",
    dataset = "sleepstudy",
    formula = Reaction ~ Days + (1 + Days | Subject),
    reml = TRUE
  ),
  list(
    id = "dyestuff_random_intercept",
    dataset = "Dyestuff",
    formula = Yield ~ 1 + (1 | Batch),
    reml = TRUE
  ),
  list(
    id = "dyestuff2_singular_random_intercept",
    dataset = "Dyestuff2",
    formula = Yield ~ 1 + (1 | Batch),
    reml = TRUE
  ),
  list(
    id = "penicillin_crossed_intercepts",
    dataset = "Penicillin",
    formula = diameter ~ 1 + (1 | plate) + (1 | sample),
    reml = TRUE
  ),
  list(
    id = "pastes_two_intercepts",
    dataset = "Pastes",
    formula = strength ~ 1 + (1 | batch) + (1 | cask),
    reml = TRUE
  ),
  list(
    id = "cake_recipe_temp",
    dataset = "cake",
    formula = angle ~ recipe * temp + (1 | recipe:replicate),
    reml = TRUE
  )
)

load_lme4_data <- function(name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "lme4", envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) {
    stop(sprintf("Could not load lme4 dataset `%s`.", name), call. = FALSE)
  }
  get(name, envir = env, inherits = FALSE)
}

comparison_key <- function(x) {
  x <- gsub(": ", "", x, fixed = TRUE)
  gsub(" & ", ":", x, fixed = TRUE)
}

level_key <- function(x) {
  gsub("_", ":", x, fixed = TRUE)
}

fixed_effect_comparison <- function(mm_fit, lme4_fit) {
  mm_beta <- fixef(mm_fit)
  lme4_beta <- lme4::fixef(lme4_fit)
  names(mm_beta) <- comparison_key(names(mm_beta))
  names(lme4_beta) <- comparison_key(names(lme4_beta))
  terms <- union(names(mm_beta), names(lme4_beta))
  data.frame(
    term = terms,
    mixeff = unname(mm_beta[terms]),
    lme4 = unname(lme4_beta[terms]),
    difference = unname(mm_beta[terms] - lme4_beta[terms]),
    row.names = NULL,
    check.names = FALSE
  )
}

fit_stat_comparison <- function(mm_fit, lme4_fit) {
  data.frame(
    statistic = c("sigma", "logLik", "AIC", "BIC", "nobs"),
    mixeff = c(
      sigma(mm_fit),
      as.numeric(logLik(mm_fit)),
      AIC(mm_fit),
      BIC(mm_fit),
      nobs(mm_fit)
    ),
    lme4 = c(
      sigma(lme4_fit),
      as.numeric(logLik(lme4_fit)),
      AIC(lme4_fit),
      BIC(lme4_fit),
      nobs(lme4_fit)
    ),
    row.names = NULL,
    check.names = FALSE
  )
}

variance_components <- function(mm_fit, lme4_fit) {
  lme4_vc <- as.data.frame(lme4::VarCorr(lme4_fit))
  lme4_vc <- lme4_vc[, c("grp", "var1", "vcov", "sdcor"), drop = FALSE]
  names(lme4_vc) <- c("group", "name", "variance", "std_dev")

  mm_vc <- VarCorr(mm_fit)$table
  mm_vc <- mm_vc[, c("group", "name", "variance", "std_dev"), drop = FALSE]
  mm_vc <- rbind(
    mm_vc,
    data.frame(
      group = "Residual",
      name = NA_character_,
      variance = sigma(mm_fit)^2,
      std_dev = sigma(mm_fit),
      check.names = FALSE
    )
  )

  list(mixeff = mm_vc, lme4 = lme4_vc)
}

random_effect_summary <- function(mm_fit, lme4_fit) {
  mm_re <- ranef(mm_fit)
  lme4_re <- lme4::ranef(lme4_fit)
  names(mm_re) <- comparison_key(names(mm_re))
  names(lme4_re) <- comparison_key(names(lme4_re))
  common_groups <- intersect(names(mm_re), names(lme4_re))
  if (!length(common_groups)) {
    return(data.frame())
  }

  rows <- lapply(common_groups, function(group) {
    mm_df <- as.data.frame(mm_re[[group]])
    lme4_df <- as.data.frame(lme4_re[[group]])
    rownames(mm_df) <- level_key(rownames(mm_df))
    rownames(lme4_df) <- level_key(rownames(lme4_df))
    common_cols <- intersect(names(mm_df), names(lme4_df))
    common_levels <- intersect(rownames(mm_df), rownames(lme4_df))
    if (!length(common_cols) || !length(common_levels)) {
      return(NULL)
    }

    diff <- as.matrix(mm_df[common_levels, common_cols, drop = FALSE]) -
      as.matrix(lme4_df[common_levels, common_cols, drop = FALSE])
    data.frame(
      group = group,
      terms = paste(common_cols, collapse = ", "),
      levels = length(common_levels),
      max_abs_difference = max(abs(diff)),
      row.names = NULL,
      check.names = FALSE
    )
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

round_numeric_columns <- function(x, digits = 7) {
  numeric_cols <- vapply(x, is.numeric, logical(1))
  x[numeric_cols] <- lapply(x[numeric_cols], signif, digits = digits)
  x
}

timed_fit <- function(expr, reps = benchmark_reps) {
  expr <- substitute(expr)
  env <- parent.frame()
  timings <- numeric(reps)
  fit <- NULL
  for (i in seq_len(reps)) {
    timing <- system.time({
      fit <- eval(expr, envir = env)
    })
    timings[[i]] <- unname(timing[["elapsed"]])
  }
  list(fit = fit, elapsed = timings)
}

benchmark_comparison <- function(mm_elapsed, lme4_elapsed) {
  data.frame(
    engine = c("mixeff", "lme4"),
    reps = c(length(mm_elapsed), length(lme4_elapsed)),
    med_sec = c(stats::median(mm_elapsed), stats::median(lme4_elapsed)),
    min_sec = c(min(mm_elapsed), min(lme4_elapsed)),
    max_sec = c(max(mm_elapsed), max(lme4_elapsed)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

print_case <- function(case) {
  data <- load_lme4_data(case$dataset)
  mm_timed <- timed_fit(
    lmm(case$formula, data, REML = case$reml,
        control = mm_control(verbose = -1))
  )
  lme4_timed <- timed_fit(
    suppressMessages(suppressWarnings(
      lme4::lmer(case$formula, data = data, REML = case$reml)
    ))
  )
  mm_fit <- mm_timed$fit
  lme4_fit <- lme4_timed$fit

  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(case$id, "\n", sep = "")
  cat("Formula: ", deparse1(case$formula), "\n", sep = "")
  cat("Dataset: ", case$dataset, "; REML: ", case$reml, "\n", sep = "")
  cat("Benchmark repetitions per engine: ", benchmark_reps, "\n", sep = "")

  cat("\nFit timing\n")
  timing <- benchmark_comparison(mm_timed$elapsed, lme4_timed$elapsed)
  timing$ratio <- timing$med_sec / timing$med_sec[timing$engine == "lme4"]
  print(round_numeric_columns(timing), row.names = FALSE)

  cat("\nFixed effects\n")
  print(round_numeric_columns(fixed_effect_comparison(mm_fit, lme4_fit)),
        row.names = FALSE)

  cat("\nFit statistics\n")
  stats <- fit_stat_comparison(mm_fit, lme4_fit)
  stats$difference <- stats$mixeff - stats$lme4
  print(round_numeric_columns(stats), row.names = FALSE)

  cat("\nVariance components: mixeff\n")
  vc <- variance_components(mm_fit, lme4_fit)
  print(round_numeric_columns(vc$mixeff), row.names = FALSE)

  cat("\nVariance components: lme4\n")
  print(round_numeric_columns(vc$lme4), row.names = FALSE)

  re_summary <- random_effect_summary(mm_fit, lme4_fit)
  if (is.data.frame(re_summary) && nrow(re_summary)) {
    cat("\nRandom-effect modes\n")
    print(round_numeric_columns(re_summary), row.names = FALSE)
  }

  invisible(list(
    mixeff = mm_fit,
    lme4 = lme4_fit,
    timing = timing,
    mixeff_elapsed = mm_timed$elapsed,
    lme4_elapsed = lme4_timed$elapsed
  ))
}

results <- lapply(cases, print_case)
names(results) <- vapply(cases, `[[`, character(1), "id")

cat("\nCompleted ", length(results), " lme4 parity examples.\n", sep = "")
