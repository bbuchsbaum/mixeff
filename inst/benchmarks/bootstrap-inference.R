#!/usr/bin/env Rscript

# Benchmark mixeff bootstrap inference routes against common R baselines.
# Run from the package root with:
#   Rscript inst/benchmarks/bootstrap-inference.R --nsim=200 --reps=3

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--args")) {
      next
    }
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
    "Usage: Rscript inst/benchmarks/bootstrap-inference.R [options]\n",
    "\n",
    "Options:\n",
    "  --nsim=N       Bootstrap replicates per timing run [200]\n",
    "  --reps=N       Timing repetitions per route [3]\n",
    "  --seed=N       Base random seed [1]\n",
    "  --out=DIR      Output directory [benchmarks/bootstrap-inference]\n",
    "  --help         Show this message\n",
    sep = ""
  )
}

scalar_arg <- function(args, name, default) {
  value <- args[[name]] %||% Sys.getenv(paste0("MIXEFF_BOOT_", toupper(name)), unset = "")
  if (!nzchar(value)) value <- default
  value
}

int_arg <- function(args, name, default, min = 1L) {
  value <- suppressWarnings(as.integer(scalar_arg(args, name, default)))
  if (length(value) != 1L || is.na(value) || value < min) {
    stop(sprintf("`%s` must be an integer >= %d.", name, min), call. = FALSE)
  }
  value
}

`%||%` <- function(x, y) if (is.null(x)) y else x

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("This benchmark requires the %s package.", package), call. = FALSE)
  }
}

make_data <- function(seed = 1L) {
  set.seed(seed)
  n_subj <- 18L
  days <- 0:9
  b0 <- stats::rnorm(n_subj, sd = 30)
  b1 <- 0.5 * b0 / 30 * 10 + stats::rnorm(n_subj, sd = 0.5)
  do.call(rbind, lapply(seq_len(n_subj), function(i) {
    data.frame(
      subj = factor(i),
      days = days,
      rt = 250 + b0[i] + (10 + b1[i]) * days +
        stats::rnorm(length(days), sd = 20)
    )
  }))
}

elapsed_once <- function(expr) {
  gc()
  unname(system.time(force(expr))[["elapsed"]])
}

time_route <- function(route, reps, expr_factory) {
  rows <- vector("list", reps)
  for (i in seq_len(reps)) {
    expr <- expr_factory(i)
    rows[[i]] <- data.frame(
      route = route,
      rep = i,
      elapsed_sec = elapsed_once(eval(expr, parent.frame())),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarise <- function(raw) {
  split_rows <- split(raw, raw$route)
  out <- lapply(split_rows, function(x) {
    data.frame(
      route = x$route[[1L]],
      reps = nrow(x),
      median_sec = stats::median(x$elapsed_sec),
      min_sec = min(x$elapsed_sec),
      max_sec = max(x$elapsed_sec),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

args <- parse_args(commandArgs(TRUE))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0L)
}

require_package("mixeff")
require_package("lme4")

nsim <- int_arg(args, "nsim", "200")
reps <- int_arg(args, "reps", "3")
seed <- int_arg(args, "seed", "1", min = 0L)
out_dir <- scalar_arg(args, "out", "benchmarks/bootstrap-inference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- make_data(seed)
formula <- rt ~ days + (1 + days | subj)
reduced_formula <- rt ~ 1 + (1 + days | subj)

mixeff_fit <- mixeff::lmm(formula, dat, REML = FALSE,
                          control = mixeff::mm_control(verbose = -1))
mixeff_reduced <- mixeff::lmm(reduced_formula, dat, REML = FALSE,
                              control = mixeff::mm_control(verbose = -1))
lme4_fit <- suppressMessages(suppressWarnings(
  lme4::lmer(formula, data = dat, REML = FALSE)
))
lme4_reduced <- suppressMessages(suppressWarnings(
  lme4::lmer(reduced_formula, data = dat, REML = FALSE)
))

routes <- list(
  time_route("mixeff_test_effect_bootstrap", reps, function(i) substitute(
    mixeff::test_effect(
      mixeff_fit,
      "days",
      method = "bootstrap",
      bootstrap = mixeff::bootstrap_control(nsim = NSIM, seed = SEED)
    ),
    list(NSIM = nsim, SEED = seed + i)
  )),
  time_route("mixeff_confint_bootstrap", reps, function(i) substitute(
    stats::confint(
      mixeff_fit,
      parm = "days",
      method = "bootstrap",
      bootstrap = mixeff::bootstrap_control(nsim = NSIM, seed = SEED)
    ),
    list(NSIM = nsim, SEED = seed + 1000L + i)
  )),
  time_route("mixeff_compare_bootstrap_lrt", reps, function(i) substitute(
    mixeff::compare(mixeff_reduced, mixeff_fit, method = "bootstrap",
                    nsim = NSIM, seed = SEED),
    list(NSIM = nsim, SEED = seed + 2000L + i)
  )),
  time_route("lme4_bootMer_fixef_distribution", reps, function(i) substitute(
    lme4::bootMer(lme4_fit, FUN = lme4::fixef, nsim = NSIM,
                  seed = SEED, use.u = FALSE, type = "parametric"),
    list(NSIM = nsim, SEED = seed + 3000L + i)
  ))
)

if (requireNamespace("pbkrtest", quietly = TRUE)) {
  routes[[length(routes) + 1L]] <- time_route("pbkrtest_PBmodcomp", reps, function(i) substitute(
    pbkrtest::PBmodcomp(lme4_fit, lme4_reduced, nsim = NSIM, seed = SEED),
    list(NSIM = nsim, SEED = seed + 4000L + i)
  ))
}

raw <- do.call(rbind, routes)
raw$nsim <- nsim
raw$nobs <- nrow(dat)
summary <- summarise(raw)
summary$nsim <- nsim
summary$nobs <- nrow(dat)

utils::write.csv(raw, file.path(out_dir, "bootstrap-inference-raw.csv"),
                 row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "bootstrap-inference-summary.csv"),
                 row.names = FALSE)

cat("Wrote:\n")
cat(sprintf("  %s\n", file.path(out_dir, "bootstrap-inference-raw.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "bootstrap-inference-summary.csv")))
print(summary)
