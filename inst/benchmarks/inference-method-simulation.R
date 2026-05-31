#!/usr/bin/env Rscript

# Small calibration scaffold for inference route behavior.
# Fast mode is meant for package checks and vignette data. Slow mode adds
# bootstrap routes and more replications for local evidence gathering.

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--args")) next
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    if (!grepl("^--[^=]+=", arg)) {
      stop(sprintf("Unknown argument `%s`. Use --help for usage.", arg), call. = FALSE)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    out[[gsub("-", "_", key, fixed = TRUE)]] <- value
  }
  out
}

usage <- function() {
  cat(
    "Usage: Rscript inst/benchmarks/inference-method-simulation.R [options]\n",
    "\n",
    "Options:\n",
    "  --mode=fast|slow   fast excludes bootstrap routes [fast]\n",
    "  --reps=N           replications per fixture/effect condition [fast: 3, slow: 100]\n",
    "  --nsim=N           bootstrap replicates in slow mode [199]\n",
    "  --seed=N           base random seed [2026]\n",
    "  --out=FILE         output CSV [inst/extdata/inference-method-simulation-summary.csv]\n",
    "  --help             show this message\n",
    sep = ""
  )
}

scalar_arg <- function(args, name, default) {
  value <- args[[name]] %||% ""
  if (!nzchar(value)) default else value
}

int_arg <- function(args, name, default, min = 1L) {
  value <- suppressWarnings(as.integer(scalar_arg(args, name, default)))
  if (length(value) != 1L || is.na(value) || value < min) {
    stop(sprintf("`%s` must be an integer >= %d.", name, min), call. = FALSE)
  }
  value
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("This simulation requires the %s package.", package), call. = FALSE)
  }
}

make_fixture <- function(fixture, beta, seed) {
  set.seed(seed)
  n_subject <- if (identical(fixture, "small_group")) 5L else 14L
  n_per <- if (identical(fixture, "small_group")) 4L else 6L
  subject <- factor(rep(seq_len(n_subject), each = n_per))
  x_base <- seq(-0.75, 0.75, length.out = n_per)
  x <- rep(x_base, n_subject)
  b0 <- stats::rnorm(n_subject, sd = 0.45)
  b1 <- switch(
    fixture,
    interior = stats::rnorm(n_subject, sd = 0.25),
    boundary = rep(0, n_subject),
    reduced_rank = 0.55 * b0,
    small_group = rep(0, n_subject),
    stop(sprintf("unknown fixture `%s`", fixture), call. = FALSE)
  )
  y <- 1 + beta * x + b0[as.integer(subject)] + b1[as.integer(subject)] * x +
    stats::rnorm(length(x), sd = 0.35)
  data.frame(y = y, x = x, subject = subject)
}

fixture_formula <- function(fixture) {
  if (identical(fixture, "small_group")) {
    y ~ x + (1 | subject)
  } else {
    y ~ x + (1 + x | subject)
  }
}

fit_fixture <- function(fixture, beta, seed) {
  dat <- make_fixture(fixture, beta = beta, seed = seed)
  mixeff::lmm(
    fixture_formula(fixture),
    dat,
    REML = FALSE,
    control = mixeff::mm_control(verbose = -1)
  )
}

ci_from_contrast <- function(row, true_beta) {
  if (!identical(row$status[[1L]], "available") ||
      !is.finite(row$estimate[[1L]]) ||
      !is.finite(row$std_error[[1L]])) {
    return(NA)
  }
  crit <- if (is.finite(row$df[[1L]]) && row$df[[1L]] > 0) {
    stats::qt(0.975, df = row$df[[1L]])
  } else {
    stats::qnorm(0.975)
  }
  lower <- row$estimate[[1L]] - crit * row$std_error[[1L]]
  upper <- row$estimate[[1L]] + crit * row$std_error[[1L]]
  isTRUE(lower <= true_beta && true_beta <= upper)
}

contrast_result <- function(fit, method, true_beta, nsim, seed) {
  L <- c(0, 1)
  names(L) <- names(mixeff::fixef(fit))
  method_arg <- switch(
    method,
    asymptotic_wald_z = "asymptotic",
    satterthwaite = "satterthwaite",
    kenward_roger = "kenward_roger",
    bootstrap = "bootstrap",
    stop(sprintf("unsupported contrast method `%s`", method), call. = FALSE)
  )
  out <- if (identical(method, "bootstrap")) {
    mixeff::contrast(
      fit,
      L,
      method = method_arg,
      bootstrap = mixeff::bootstrap_control(nsim = nsim, seed = seed)
    )
  } else {
    mixeff::contrast(fit, L, method = method_arg)
  }
  row <- out$table[1L, , drop = FALSE]
  list(
    p_value = row$p_value[[1L]],
    coverage = ci_from_contrast(row, true_beta)
  )
}

profile_result <- function(fit, true_beta) {
  ci <- stats::confint(fit, parm = "x", method = "profile")
  if (!("x" %in% rownames(ci)) || any(!is.finite(ci["x", ]))) {
    return(list(p_value = NA_real_, coverage = NA))
  }
  list(
    p_value = NA_real_,
    coverage = isTRUE(ci["x", 1L] <= true_beta && true_beta <= ci["x", 2L])
  )
}

bootstrap_lrt_result <- function(fit, nsim, seed) {
  out <- mixeff::test_effect(
    fit,
    "x",
    method = "bootstrap_lrt",
    bootstrap = mixeff::bootstrap_control(nsim = nsim, seed = seed)
  )
  list(
    p_value = out$table$p_value[[1L]],
    coverage = NA
  )
}

evaluate_method <- function(fit, method, true_beta, nsim, seed) {
  tryCatch(
    switch(
      method,
      asymptotic_wald_z = contrast_result(fit, method, true_beta, nsim, seed),
      satterthwaite = contrast_result(fit, method, true_beta, nsim, seed),
      kenward_roger = contrast_result(fit, method, true_beta, nsim, seed),
      bootstrap = contrast_result(fit, method, true_beta, nsim, seed),
      profile = profile_result(fit, true_beta),
      bootstrap_lrt = bootstrap_lrt_result(fit, nsim, seed),
      boundary_lrt = list(p_value = NA_real_, coverage = NA),
      stop(sprintf("unknown method `%s`", method), call. = FALSE)
    ),
    error = function(cnd) list(p_value = NA_real_, coverage = NA)
  )
}

metric_mean <- function(x) {
  keep <- !is.na(x)
  if (!any(keep)) return(NA_real_)
  mean(x[keep])
}

simulate_cell <- function(fixture, method, mode, reps, nsim, seed, alt_beta) {
  slow_only <- c("bootstrap", "bootstrap_lrt")
  not_wired <- "boundary_lrt"
  if ((identical(mode, "fast") && method %in% slow_only) ||
      method %in% not_wired) {
    return(data.frame(
      method = method,
      fixture = fixture,
      mode = mode,
      n_reps = 0L,
      type_I_error = NA_real_,
      power_at_alt = NA_real_,
      coverage_at_alt = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  null_p <- alt_p <- alt_coverage <- numeric(reps)
  for (i in seq_len(reps)) {
    null_fit <- tryCatch(fit_fixture(fixture, beta = 0, seed = seed + i),
                         error = function(cnd) NULL)
    alt_fit <- tryCatch(fit_fixture(fixture, beta = alt_beta, seed = seed + 10000L + i),
                        error = function(cnd) NULL)
    null_res <- if (is.null(null_fit)) {
      list(p_value = NA_real_, coverage = NA)
    } else {
      evaluate_method(null_fit, method, true_beta = 0, nsim = nsim,
                      seed = seed + 20000L + i)
    }
    alt_res <- if (is.null(alt_fit)) {
      list(p_value = NA_real_, coverage = NA)
    } else {
      evaluate_method(alt_fit, method, true_beta = alt_beta, nsim = nsim,
                      seed = seed + 30000L + i)
    }
    null_p[[i]] <- null_res$p_value
    alt_p[[i]] <- alt_res$p_value
    alt_coverage[[i]] <- alt_res$coverage
  }

  data.frame(
    method = method,
    fixture = fixture,
    mode = mode,
    n_reps = reps,
    type_I_error = metric_mean(null_p < 0.05),
    power_at_alt = metric_mean(alt_p < 0.05),
    coverage_at_alt = metric_mean(alt_coverage),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(TRUE))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0L)
}

require_package("mixeff")

mode <- match.arg(scalar_arg(args, "mode", "fast"), c("fast", "slow"))
reps <- int_arg(args, "reps", if (identical(mode, "fast")) "3" else "100")
nsim <- int_arg(args, "nsim", if (identical(mode, "fast")) "19" else "199")
seed <- int_arg(args, "seed", "2026", min = 0L)
out <- scalar_arg(args, "out", "inst/extdata/inference-method-simulation-summary.csv")
alt_beta <- as.numeric(scalar_arg(args, "alt_beta", "0.35"))

fixtures <- c("interior", "boundary", "reduced_rank", "small_group")
methods <- c(
  "asymptotic_wald_z",
  "satterthwaite",
  "kenward_roger",
  "bootstrap",
  "profile",
  "bootstrap_lrt",
  "boundary_lrt"
)

rows <- vector("list", length(fixtures) * length(methods))
k <- 1L
for (fixture in fixtures) {
  for (method in methods) {
    rows[[k]] <- simulate_cell(
      fixture = fixture,
      method = method,
      mode = mode,
      reps = reps,
      nsim = nsim,
      seed = seed + k * 1000L,
      alt_beta = alt_beta
    )
    k <- k + 1L
  }
}
summary <- do.call(rbind, rows)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, out, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", out, nrow(summary)))
print(summary)
