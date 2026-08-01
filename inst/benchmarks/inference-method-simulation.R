#!/usr/bin/env Rscript

# Calibration study for mixeff inference routes (Audit1 WI-3.1, DEC-2).
#
# Two modes with different epistemic status:
#   * smoke   -- 3 replications; proves the pipeline executes. NOT evidence.
#                The vignette refuses to present smoke output as calibration.
#   * release -- 2,500 replications for analytic routes, 500 for bootstrap
#                routes (Monte Carlo SE at a true 0.05 rate: ~0.0044 and
#                ~0.0097). This is the mode that produces the shipped
#                artifact, alongside a machine-readable manifest recording
#                package/engine provenance and the exact invocation.
#
# Design principles (audit section 3):
#   * requested / successful / failed / refused / boundary replicate counts
#     are all recorded; nothing is silently dropped;
#   * every rate is reported with a Wilson 95% interval, both conditional on
#     a successful evaluation and unconditional (refusals and failures count
#     as non-rejections / non-coverage);
#   * deterministic per-cell seeds; the artifact is reproducible from the
#     documented invocation printed into the manifest;
#   * routes that do not run carry an explicit `reason`, so "not wired" is
#     distinguishable from "excluded in this mode".

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
    "  --mode=smoke|release  smoke = pipeline check (NOT evidence) [smoke]\n",
    "  --reps=N              analytic-route reps [smoke: 3, release: 2500]\n",
    "  --reps-bootstrap=N    bootstrap-route reps [smoke: 0, release: 500]\n",
    "  --reps-joint=N        certified joint-GLMM reps [smoke: 3, release: 1000]\n",
    "  --nsim=N              bootstrap replicates per bootstrap rep [199]\n",
    "  --seed=N              base random seed [2026]\n",
    "  --out=FILE            output CSV [inst/extdata/inference-method-simulation-summary.csv]\n",
    "  --manifest=FILE       manifest JSON [inst/extdata/inference-method-simulation-manifest.json]\n",
    "  --help                show this message\n",
    sep = ""
  )
}

scalar_arg <- function(args, name, default) {
  value <- args[[name]] %||% ""
  if (!nzchar(value)) default else value
}

int_arg <- function(args, name, default, min = 0L) {
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

wilson <- function(k, n, level = 0.95) {
  if (!is.finite(k) || !is.finite(n) || n <= 0) {
    return(c(NA_real_, NA_real_))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- k / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, center - half), min(1, center + half))
}

# ---- fixtures ---------------------------------------------------------------

make_lmm_fixture <- function(fixture, beta, seed) {
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

lmm_fixture_formula <- function(fixture) {
  if (identical(fixture, "small_group")) {
    y ~ x + (1 | subject)
  } else {
    y ~ x + (1 + x | subject)
  }
}

fit_lmm_fixture <- function(fixture, beta, seed, reml = FALSE) {
  dat <- make_lmm_fixture(fixture, beta = beta, seed = seed)
  mixeff::lmm(
    lmm_fixture_formula(fixture),
    dat,
    REML = reml,
    control = mixeff::mm_control(verbose = -1)
  )
}

make_glmm_fixture <- function(beta, seed) {
  set.seed(seed)
  n_group <- 20L
  n_per <- 8L
  g <- factor(rep(seq_len(n_group), each = n_per))
  x <- rep(seq(-1, 1, length.out = n_per), n_group)
  b0 <- stats::rnorm(n_group, sd = 0.5)
  eta <- -0.2 + beta * x + b0[as.integer(g)]
  data.frame(y = stats::rbinom(length(eta), 1L, stats::plogis(eta)), x = x, g = g)
}

fit_glmm_fixture <- function(beta, seed, method = "pirls_profiled",
                             inference = "auto") {
  dat <- make_glmm_fixture(beta = beta, seed = seed)
  mixeff::glmm(
    y ~ x + (1 | g),
    dat,
    family = stats::binomial(),
    method = method,
    inference = inference,
    control = mixeff::mm_control(verbose = -1)
  )
}

# ---- per-replicate evaluation ----------------------------------------------
#
# Every evaluator returns list(outcome, p_value, reject, coverage):
#   outcome  "evaluated" | "refused" | "failed"
#   reject   TRUE/FALSE under the cell's test (p < .05, or CI excludes truth
#            for interval-only routes); NA when not evaluated
#   coverage CI contains the true beta; NA when the route has no interval

ci_covers <- function(lower, upper, truth) {
  isTRUE(is.finite(lower) && is.finite(upper) && lower <= truth && truth <= upper)
}

eval_outcome <- function(expr) {
  tryCatch(
    expr,
    mm_inference_unavailable = function(cnd) {
      list(outcome = "refused", p_value = NA_real_, reject = NA, coverage = NA)
    },
    error = function(cnd) {
      list(outcome = "failed", p_value = NA_real_, reject = NA, coverage = NA)
    }
  )
}

lmm_contrast_eval <- function(fit, method_arg, truth, nsim, seed) {
  eval_outcome({
    L <- c(0, 1)
    names(L) <- names(mixeff::fixef(fit))
    out <- if (identical(method_arg, "bootstrap")) {
      mixeff::contrast(fit, L, method = method_arg,
                       bootstrap = mixeff::bootstrap_control(nsim = nsim, seed = seed))
    } else {
      mixeff::contrast(fit, L, method = method_arg)
    }
    row <- out$table[1L, , drop = FALSE]
    if (!identical(row$status[[1L]], "available") ||
        !is.finite(row$p_value[[1L]])) {
      list(outcome = "refused", p_value = NA_real_, reject = NA, coverage = NA)
    } else {
      crit <- if (is.finite(row$df[[1L]]) && row$df[[1L]] > 0) {
        stats::qt(0.975, df = row$df[[1L]])
      } else {
        stats::qnorm(0.975)
      }
      lower <- row$estimate[[1L]] - crit * row$std_error[[1L]]
      upper <- row$estimate[[1L]] + crit * row$std_error[[1L]]
      list(outcome = "evaluated", p_value = row$p_value[[1L]],
           reject = isTRUE(row$p_value[[1L]] < 0.05),
           coverage = ci_covers(lower, upper, truth))
    }
  })
}

lmm_profile_eval <- function(fit, truth) {
  eval_outcome({
    ci <- stats::confint(fit, parm = "x", method = "profile")
    if (!("x" %in% rownames(ci)) || any(!is.finite(ci["x", ]))) {
      list(outcome = "refused", p_value = NA_real_, reject = NA, coverage = NA)
    } else {
      list(outcome = "evaluated", p_value = NA_real_,
           reject = !ci_covers(ci["x", 1L], ci["x", 2L], 0),
           coverage = ci_covers(ci["x", 1L], ci["x", 2L], truth))
    }
  })
}

lmm_bootstrap_lrt_eval <- function(fit, nsim, seed) {
  eval_outcome({
    out <- mixeff::test_effect(
      fit, "x", method = "bootstrap_lrt",
      bootstrap = mixeff::bootstrap_control(nsim = nsim, seed = seed)
    )
    p <- out$table$p_value[[1L]]
    if (!is.finite(p)) {
      list(outcome = "refused", p_value = NA_real_, reject = NA, coverage = NA)
    } else {
      list(outcome = "evaluated", p_value = p, reject = isTRUE(p < 0.05),
           coverage = NA)
    }
  })
}

glmm_wald_eval <- function(fit, truth) {
  eval_outcome({
    L <- c(0, 1)
    names(L) <- names(mixeff::fixef(fit))
    out <- mixeff::contrast(fit, L)
    row <- out$table[1L, , drop = FALSE]
    if (!is.finite(row$p_value[[1L]])) {
      list(outcome = "refused", p_value = NA_real_, reject = NA, coverage = NA)
    } else {
      crit <- stats::qnorm(0.975)
      lower <- row$estimate[[1L]] - crit * row$std_error[[1L]]
      upper <- row$estimate[[1L]] + crit * row$std_error[[1L]]
      list(outcome = "evaluated", p_value = row$p_value[[1L]],
           reject = isTRUE(row$p_value[[1L]] < 0.05),
           coverage = ci_covers(lower, upper, truth))
    }
  })
}

glmm_bootstrap_eval <- function(fit, truth, nsim, seed) {
  eval_outcome({
    ci <- stats::confint(fit, parm = "x", method = "bootstrap",
                         nsim = nsim, seed = seed)
    list(outcome = "evaluated", p_value = NA_real_,
         reject = !ci_covers(ci["x", 1L], ci["x", 2L], 0),
         coverage = ci_covers(ci["x", 1L], ci["x", 2L], truth))
  })
}

# ---- cell registry ----------------------------------------------------------
#
# Each cell: fixture label, method label, model family ("lmm"/"glmm"), a
# fitter(beta, seed) and an evaluator(fit, truth, nsim, seed), plus which rep
# budget it draws from ("analytic", "bootstrap", "joint").

lmm_cells <- function() {
  fixtures <- c("interior", "boundary", "reduced_rank", "small_group")
  analytic <- list(
    asymptotic_wald_z = function(fit, truth, nsim, seed)
      lmm_contrast_eval(fit, "asymptotic", truth, nsim, seed),
    satterthwaite = function(fit, truth, nsim, seed)
      lmm_contrast_eval(fit, "satterthwaite", truth, nsim, seed),
    kenward_roger = function(fit, truth, nsim, seed)
      lmm_contrast_eval(fit, "kenward_roger", truth, nsim, seed),
    profile = function(fit, truth, nsim, seed) lmm_profile_eval(fit, truth)
  )
  boot <- list(
    bootstrap = function(fit, truth, nsim, seed)
      lmm_contrast_eval(fit, "bootstrap", truth, nsim, seed),
    bootstrap_lrt = function(fit, truth, nsim, seed)
      lmm_bootstrap_lrt_eval(fit, nsim, seed)
  )
  cells <- list()
  for (fx in fixtures) {
    fitter <- local({
      fx_local <- fx
      function(beta, seed) fit_lmm_fixture(fx_local, beta, seed)
    })
    # Kenward-Roger is certified for REML fits only (the study's other
    # routes need ML: profile-beta and bootstrap-LRT are ML-only). Fitting
    # the KR cell under ML would make the whole cell a wall of policy
    # refusals carrying zero calibration information.
    fitter_reml <- local({
      fx_local <- fx
      function(beta, seed) fit_lmm_fixture(fx_local, beta, seed, reml = TRUE)
    })
    for (m in names(analytic)) {
      cells[[length(cells) + 1L]] <- list(
        fixture = fx, method = m, model = "lmm", budget = "analytic",
        fitter = if (identical(m, "kenward_roger")) fitter_reml else fitter,
        evaluator = analytic[[m]]
      )
    }
    for (m in names(boot)) {
      cells[[length(cells) + 1L]] <- list(
        fixture = fx, method = m, model = "lmm", budget = "bootstrap",
        fitter = fitter, evaluator = boot[[m]]
      )
    }
  }
  cells
}

glmm_cells <- function() {
  list(
    list(
      fixture = "glmm_binomial", method = "glmm_wald_joint_laplace",
      model = "glmm", budget = "joint",
      fitter = function(beta, seed)
        fit_glmm_fixture(beta, seed, method = "joint_laplace"),
      evaluator = function(fit, truth, nsim, seed) glmm_wald_eval(fit, truth)
    ),
    list(
      fixture = "glmm_binomial", method = "wald_z_working_hessian",
      model = "glmm", budget = "analytic",
      fitter = function(beta, seed)
        fit_glmm_fixture(beta, seed, inference = "working_hessian"),
      evaluator = function(fit, truth, nsim, seed) glmm_wald_eval(fit, truth)
    ),
    list(
      fixture = "glmm_binomial", method = "glmm_parametric_bootstrap",
      model = "glmm", budget = "bootstrap",
      fitter = function(beta, seed) fit_glmm_fixture(beta, seed),
      evaluator = function(fit, truth, nsim, seed)
        glmm_bootstrap_eval(fit, truth, nsim, seed)
    )
  )
}

# ---- cell simulation --------------------------------------------------------

boundary_status <- function(fit) {
  status <- as.character(fit$fit_status %||% "")
  status %in% c("converged_boundary", "converged_reduced_rank")
}

simulate_cell <- function(cell, reps, nsim, seed, alt_beta, mode) {
  method <- cell$method
  fixture <- cell$fixture
  if (reps < 1L) {
    return(cell_row(
      cell, mode, requested = 0L,
      reason = sprintf("excluded_in_%s_mode", mode)
    ))
  }

  counts <- list(
    null = c(evaluated = 0L, refused = 0L, failed = 0L, fit_failed = 0L,
             boundary = 0L, reject = 0L),
    alt = c(evaluated = 0L, refused = 0L, failed = 0L, fit_failed = 0L,
            boundary = 0L, reject = 0L, cover = 0L, has_ci = 0L)
  )

  for (i in seq_len(reps)) {
    for (arm in c("null", "alt")) {
      truth <- if (identical(arm, "null")) 0 else alt_beta
      arm_offset <- if (identical(arm, "null")) 0L else 10000L
      fit <- tryCatch(
        cell$fitter(truth, seed + arm_offset + i),
        error = function(cnd) NULL
      )
      if (is.null(fit)) {
        counts[[arm]][["fit_failed"]] <- counts[[arm]][["fit_failed"]] + 1L
        next
      }
      if (boundary_status(fit)) {
        counts[[arm]][["boundary"]] <- counts[[arm]][["boundary"]] + 1L
      }
      res <- cell$evaluator(fit, truth, nsim, seed + 20000L + arm_offset + i)
      counts[[arm]][[res$outcome]] <- counts[[arm]][[res$outcome]] + 1L
      if (identical(res$outcome, "evaluated")) {
        if (isTRUE(res$reject)) {
          counts[[arm]][["reject"]] <- counts[[arm]][["reject"]] + 1L
        }
        if (identical(arm, "alt") && !is.na(res$coverage)) {
          counts[[arm]][["has_ci"]] <- counts[[arm]][["has_ci"]] + 1L
          if (isTRUE(res$coverage)) {
            counts[[arm]][["cover"]] <- counts[[arm]][["cover"]] + 1L
          }
        }
      }
    }
  }

  cell_row(cell, mode, requested = reps, counts = counts)
}

cell_row <- function(cell, mode, requested, counts = NULL, reason = NA_character_) {
  empty <- c(evaluated = NA_integer_, refused = NA_integer_,
             failed = NA_integer_, fit_failed = NA_integer_,
             boundary = NA_integer_, reject = NA_integer_)
  null_c <- counts$null %||% empty
  alt_c <- counts$alt %||% c(empty, cover = NA_integer_, has_ci = NA_integer_)

  rate <- function(k, n) if (is.na(k) || is.na(n) || n < 1) NA_real_ else k / n
  wl <- function(k, n) if (is.na(k) || is.na(n) || n < 1) c(NA_real_, NA_real_) else wilson(k, n)

  n_eval_null <- null_c[["evaluated"]]
  n_eval_alt <- alt_c[["evaluated"]]
  t1_w <- wl(null_c[["reject"]], n_eval_null)
  pw_w <- wl(alt_c[["reject"]], n_eval_alt)
  cv_w <- wl(alt_c[["cover"]], alt_c[["has_ci"]])

  data.frame(
    method = cell$method,
    fixture = cell$fixture,
    model = cell$model,
    mode = mode,
    n_requested = requested,
    n_evaluated_null = n_eval_null,
    n_refused_null = null_c[["refused"]],
    n_failed_null = null_c[["failed"]] ,
    n_fit_failed_null = null_c[["fit_failed"]],
    n_boundary_null = null_c[["boundary"]],
    n_evaluated_alt = n_eval_alt,
    n_refused_alt = alt_c[["refused"]],
    n_failed_alt = alt_c[["failed"]],
    n_fit_failed_alt = alt_c[["fit_failed"]],
    n_boundary_alt = alt_c[["boundary"]],
    type_I_error = rate(null_c[["reject"]], n_eval_null),
    type_I_lower = t1_w[[1L]],
    type_I_upper = t1_w[[2L]],
    type_I_unconditional = rate(null_c[["reject"]], requested),
    power_at_alt = rate(alt_c[["reject"]], n_eval_alt),
    power_lower = pw_w[[1L]],
    power_upper = pw_w[[2L]],
    power_unconditional = rate(alt_c[["reject"]], requested),
    coverage_at_alt = rate(alt_c[["cover"]], alt_c[["has_ci"]]),
    coverage_lower = cv_w[[1L]],
    coverage_upper = cv_w[[2L]],
    reason = reason,
    stringsAsFactors = FALSE
  )
}

# ---- provenance manifest ----------------------------------------------------

file_md5 <- function(path) {
  if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_
}

git_head <- function() {
  sha <- suppressWarnings(
    tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
             error = function(cnd) character())
  )
  if (length(sha) == 1L && nzchar(sha)) sha else NA_character_
}

build_manifest <- function(mode, reps, reps_bootstrap, reps_joint, nsim, seed,
                           alt_beta, out, elapsed_s, n_rows) {
  lock <- "src/rust/upstream/mixeff-rs.lock"
  lock_lines <- if (file.exists(lock)) readLines(lock, warn = FALSE) else character()
  rust_pin <- gsub('"', "", fixed = TRUE,
                   sub("^resolved_sha[ =:]+", "",
                       grep("^resolved_sha", lock_lines, value = TRUE)[1L] %||%
                         NA_character_))
  list(
    artifact = basename(out),
    mode = mode,
    invocation = sprintf(
      paste0("Rscript inst/benchmarks/inference-method-simulation.R ",
             "--mode=%s --reps=%d --reps-bootstrap=%d --reps-joint=%d ",
             "--nsim=%d --seed=%d"),
      mode, reps, reps_bootstrap, reps_joint, nsim, seed
    ),
    reps_analytic = reps,
    reps_bootstrap = reps_bootstrap,
    reps_joint = reps_joint,
    nsim = nsim,
    seed = seed,
    alt_beta = alt_beta,
    n_rows = n_rows,
    elapsed_seconds = round(elapsed_s, 1),
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = as.character(utils::packageVersion("mixeff")),
    package_commit = git_head(),
    rust_pin = rust_pin,
    cargo_lock_md5 = file_md5("src/rust/upstream/mixeff-rs/Cargo.lock"),
    vendor_archive_md5 = file_md5("src/rust/vendor.tar.xz"),
    r_version = R.version.string,
    platform = R.version$platform
  )
}

# ---- main --------------------------------------------------------------------

args <- parse_args(commandArgs(TRUE))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0L)
}

require_package("mixeff")
require_package("jsonlite")

mode <- match.arg(scalar_arg(args, "mode", "smoke"), c("smoke", "release"))
reps <- int_arg(args, "reps", if (identical(mode, "smoke")) "3" else "2500", min = 1L)
reps_bootstrap <- int_arg(args, "reps_bootstrap",
                          if (identical(mode, "smoke")) "0" else "500")
reps_joint <- int_arg(args, "reps_joint",
                      if (identical(mode, "smoke")) "3" else "1000", min = 1L)
nsim <- int_arg(args, "nsim", "199", min = 19L)
seed <- int_arg(args, "seed", "2026", min = 0L)
out <- scalar_arg(args, "out", "inst/extdata/inference-method-simulation-summary.csv")
manifest_path <- scalar_arg(args, "manifest",
                            "inst/extdata/inference-method-simulation-manifest.json")
alt_beta <- as.numeric(scalar_arg(args, "alt_beta", "0.35"))

cells <- c(lmm_cells(), glmm_cells())
budget <- c(analytic = reps, bootstrap = reps_bootstrap, joint = reps_joint)

t0 <- Sys.time()
rows <- vector("list", length(cells))
for (k in seq_along(cells)) {
  cell <- cells[[k]]
  cell_reps <- budget[[cell$budget]]
  message(sprintf("[%d/%d] %s / %s (%d reps)",
                  k, length(cells), cell$fixture, cell$method, cell_reps))
  rows[[k]] <- simulate_cell(
    cell = cell,
    reps = cell_reps,
    nsim = nsim,
    seed = seed + k * 100000L,
    alt_beta = alt_beta,
    mode = mode
  )
}
summary <- do.call(rbind, rows)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, out, row.names = FALSE)
manifest <- build_manifest(mode, reps, reps_bootstrap, reps_joint, nsim, seed,
                           alt_beta, out, elapsed, nrow(summary))
jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("wrote %s (%d rows) and %s in %.1fs\n",
            out, nrow(summary), manifest_path, elapsed))
print(summary[, c("method", "fixture", "mode", "n_requested",
                  "n_evaluated_null", "type_I_error", "type_I_lower",
                  "type_I_upper", "coverage_at_alt")])
