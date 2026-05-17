#!/usr/bin/env Rscript

# Benchmark mixeff::lmm() against lme4::lmer() across synthetic LMM designs.
# Run from the package root with:
#   Rscript inst/benchmarks/lme4-scaling.R

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    if (identical(arg, "--no-plots")) {
      out$plots <- FALSE
      next
    }
    if (!grepl("^--[^=]+=", arg)) {
      stop(sprintf("Unknown argument `%s`. Use --help for usage.", arg), call. = FALSE)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- value
  }
  out
}

usage <- function() {
  cat(
    "Usage: Rscript inst/benchmarks/lme4-scaling.R [options]\n",
    "\n",
    "Options:\n",
    "  --reps=N                    Timing repetitions per engine/design [3]\n",
    "  --warmup=N                  Untimed warmup fits per engine/design [1]\n",
    "  --rows=CSV                  Row counts for row-scaling cases [500,1000,2500,5000]\n",
    "  --groups=CSV                Group counts for grouped cases [25,50,100,200]\n",
    "  --crossed-levels=CSV        Subject/item levels for crossed cases [10,15,20,30]\n",
    "  --scenarios=CSV             Any of rows,groups,slopes,crossed,crossed_slope [all]\n",
    "  --out=DIR                   Output directory [benchmarks/lme4-scaling]\n",
    "  --no-plots                  Write CSV files only; otherwise write ggplot2 PDFs\n",
    "  --help                      Show this message\n",
    "\n",
    "Environment variables with matching names are also honored:\n",
    "  MIXEFF_BENCH_REPS, MIXEFF_BENCH_WARMUP, MIXEFF_BENCH_ROWS,\n",
    "  MIXEFF_BENCH_GROUPS, MIXEFF_BENCH_CROSSED_LEVELS, MIXEFF_BENCH_OUT.\n",
    sep = ""
  )
}

value_or_env <- function(args, name, env, default) {
  if (!is.null(args[[name]])) {
    return(args[[name]])
  }
  value <- Sys.getenv(env, unset = NA_character_)
  if (!is.na(value) && nzchar(value)) {
    return(value)
  }
  default
}

parse_int_list <- function(x, name) {
  values <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1L]]))
  if (!length(values) || anyNA(values) || any(values <= 0L)) {
    stop(sprintf("`%s` must be a comma-separated list of positive integers.", name),
         call. = FALSE)
  }
  unique(values)
}

parse_int_scalar <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < 0L) {
    stop(sprintf("`%s` must be a non-negative integer.", name), call. = FALSE)
  }
  value
}

parse_scenarios <- function(x) {
  scenarios <- strsplit(x, ",", fixed = TRUE)[[1L]]
  scenarios <- trimws(scenarios)
  if (identical(scenarios, "all")) {
    scenarios <- c("rows", "groups", "slopes", "crossed", "crossed_slope")
  }
  valid <- c("rows", "groups", "slopes", "crossed", "crossed_slope")
  unknown <- setdiff(scenarios, valid)
  if (length(unknown)) {
    stop(sprintf("Unknown scenario(s): %s", paste(unknown, collapse = ", ")),
         call. = FALSE)
  }
  unique(scenarios)
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("This benchmark requires the %s package.", package), call. = FALSE)
  }
}

scenario_display <- function(x) {
  labels <- c(
    rows = "Rows",
    groups = "Random-intercept groups",
    slopes = "Random slopes",
    crossed = "Crossed intercepts",
    crossed_slope = "Crossed intercepts + subject slope"
  )
  unname(labels[x])
}

balanced_factor <- function(n, levels, prefix) {
  factor(sprintf("%s%04d", prefix, rep(seq_len(levels), length.out = n)))
}

seed_for <- function(scenario, scale_value) {
  sum(utf8ToInt(scenario)) * 1000L + as.integer(scale_value)
}

simulate_rows <- function(n_rows) {
  set.seed(seed_for("rows", n_rows))
  n_subject <- 50L
  subject <- balanced_factor(n_rows, n_subject, "s")
  x <- stats::rnorm(n_rows)
  b_subject <- stats::rnorm(n_subject, sd = 0.8)
  y <- 0.4 + 1.1 * x + b_subject[subject] + stats::rnorm(n_rows, sd = 1.0)
  data.frame(y = y, x = x, subject = subject)
}

simulate_groups <- function(n_subject, obs_per_subject = 8L) {
  n_rows <- n_subject * obs_per_subject
  set.seed(seed_for("groups", n_subject))
  subject <- factor(rep(sprintf("s%04d", seq_len(n_subject)), each = obs_per_subject))
  x <- stats::rnorm(n_rows)
  b_subject <- stats::rnorm(n_subject, sd = 0.8)
  y <- 0.4 + 1.1 * x + b_subject[subject] + stats::rnorm(n_rows, sd = 1.0)
  data.frame(y = y, x = x, subject = subject)
}

simulate_slopes <- function(n_subject, obs_per_subject = 10L) {
  n_rows <- n_subject * obs_per_subject
  set.seed(seed_for("slopes", n_subject))
  subject <- factor(rep(sprintf("s%04d", seq_len(n_subject)), each = obs_per_subject))
  x <- stats::rnorm(n_rows)
  b0 <- stats::rnorm(n_subject, sd = 0.7)
  b1 <- stats::rnorm(n_subject, sd = 0.25)
  y <- 0.4 + 1.1 * x + b0[subject] + b1[subject] * x +
    stats::rnorm(n_rows, sd = 1.0)
  data.frame(y = y, x = x, subject = subject)
}

simulate_crossed <- function(n_level, slope = FALSE, reps = 2L) {
  set.seed(seed_for(if (slope) "crossed_slope" else "crossed", n_level))
  grid <- expand.grid(
    subject = factor(sprintf("s%04d", seq_len(n_level))),
    item = factor(sprintf("i%04d", seq_len(n_level))),
    rep = seq_len(reps),
    KEEP.OUT.ATTRS = FALSE
  )
  n_rows <- nrow(grid)
  x <- stats::rnorm(n_rows)
  b_subject <- stats::rnorm(n_level, sd = 0.7)
  b_item <- stats::rnorm(n_level, sd = 0.45)
  eta <- 0.4 + 1.1 * x + b_subject[grid$subject] + b_item[grid$item]
  if (isTRUE(slope)) {
    b_slope <- stats::rnorm(n_level, sd = 0.25)
    eta <- eta + b_slope[grid$subject] * x
  }
  grid$y <- eta + stats::rnorm(n_rows, sd = 1.0)
  grid$x <- x
  grid[, c("y", "x", "subject", "item")]
}

make_specs <- function(rows, groups, crossed_levels, scenarios) {
  specs <- list()
  add <- function(spec) {
    specs[[length(specs) + 1L]] <<- spec
  }

  if ("rows" %in% scenarios) {
    for (n in rows) {
      add(list(
        scenario = "rows",
        scale_var = "rows",
        scale_value = n,
        formula = y ~ x + (1 | subject),
        data = simulate_rows(n),
        n_random_effect_coefficients = 50L,
        random_terms = "(1 | subject)"
      ))
    }
  }

  if ("groups" %in% scenarios) {
    for (n_group in groups) {
      add(list(
        scenario = "groups",
        scale_var = "subject levels",
        scale_value = n_group,
        formula = y ~ x + (1 | subject),
        data = simulate_groups(n_group),
        n_random_effect_coefficients = n_group,
        random_terms = "(1 | subject)"
      ))
    }
  }

  if ("slopes" %in% scenarios) {
    for (n_group in groups) {
      add(list(
        scenario = "slopes",
        scale_var = "subject levels",
        scale_value = n_group,
        formula = y ~ x + (1 + x | subject),
        data = simulate_slopes(n_group),
        n_random_effect_coefficients = 2L * n_group,
        random_terms = "(1 + x | subject)"
      ))
    }
  }

  if ("crossed" %in% scenarios) {
    for (n_level in crossed_levels) {
      add(list(
        scenario = "crossed",
        scale_var = "subject/item levels",
        scale_value = n_level,
        formula = y ~ x + (1 | subject) + (1 | item),
        data = simulate_crossed(n_level, slope = FALSE),
        n_random_effect_coefficients = 2L * n_level,
        random_terms = "(1 | subject) + (1 | item)"
      ))
    }
  }

  if ("crossed_slope" %in% scenarios) {
    for (n_level in crossed_levels) {
      add(list(
        scenario = "crossed_slope",
        scale_var = "subject/item levels",
        scale_value = n_level,
        formula = y ~ x + (1 + x | subject) + (1 | item),
        data = simulate_crossed(n_level, slope = TRUE),
        n_random_effect_coefficients = 3L * n_level,
        random_terms = "(1 + x | subject) + (1 | item)"
      ))
    }
  }

  specs
}

fit_once <- function(engine, formula, data) {
  if (identical(engine, "mixeff")) {
    mixeff::lmm(formula, data, REML = FALSE,
                control = mixeff::mm_control(verbose = -1))
  } else {
    suppressMessages(suppressWarnings(
      lme4::lmer(formula, data = data, REML = FALSE)
    ))
  }
}

time_fit <- function(engine, spec, rep_id) {
  fit <- NULL
  error <- NA_character_
  timing <- system.time({
    fit <- tryCatch(
      fit_once(engine, spec$formula, spec$data),
      error = function(cnd) {
        error <<- conditionMessage(cnd)
        NULL
      }
    )
  })
  ok <- !is.null(fit)
  rm(fit)
  invisible(gc(verbose = FALSE))

  data.frame(
    scenario = spec$scenario,
    scale_var = spec$scale_var,
    scale_value = spec$scale_value,
    n_rows = nrow(spec$data),
    n_subject = if ("subject" %in% names(spec$data)) nlevels(spec$data$subject) else NA_integer_,
    n_item = if ("item" %in% names(spec$data)) nlevels(spec$data$item) else NA_integer_,
    n_random_effect_coefficients = spec$n_random_effect_coefficients,
    random_terms = spec$random_terms,
    formula = deparse1(spec$formula),
    engine = engine,
    rep = rep_id,
    elapsed_sec = unname(timing[["elapsed"]]),
    ok = ok,
    error = error,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

warmup_fits <- function(spec, engines, warmup) {
  if (warmup <= 0L) {
    return(invisible(NULL))
  }
  for (engine in engines) {
    for (i in seq_len(warmup)) {
      try(fit_once(engine, spec$formula, spec$data), silent = TRUE)
      invisible(gc(verbose = FALSE))
    }
  }
  invisible(NULL)
}

summarise_timings <- function(raw) {
  key_cols <- c(
    "scenario", "scale_var", "scale_value", "n_rows", "n_subject", "n_item",
    "n_random_effect_coefficients", "random_terms", "formula", "engine"
  )
  split_keys <- raw[key_cols]
  split_keys[] <- lapply(split_keys, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    x
  })
  parts <- split(raw, do.call(paste, c(split_keys, sep = "\r")), drop = TRUE)
  rows <- lapply(parts, function(x) {
    ok_elapsed <- x$elapsed_sec[isTRUE(any(x$ok)) & x$ok & !is.na(x$elapsed_sec)]
    first <- x[1L, key_cols, drop = FALSE]
    first$n_reps <- nrow(x)
    first$n_ok <- length(ok_elapsed)
    first$n_fail <- nrow(x) - length(ok_elapsed)
    first$median_sec <- if (length(ok_elapsed)) stats::median(ok_elapsed) else NA_real_
    first$min_sec <- if (length(ok_elapsed)) min(ok_elapsed) else NA_real_
    first$max_sec <- if (length(ok_elapsed)) max(ok_elapsed) else NA_real_
    first$fits_per_sec <- if (!is.na(first$median_sec) && first$median_sec > 0) {
      1 / first$median_sec
    } else {
      NA_real_
    }
    first
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary$speedup_vs_lme4 <- NA_real_

  for (i in seq_len(nrow(summary))) {
    ref <- summary[
      summary$engine == "lme4" &
        summary$scenario == summary$scenario[[i]] &
        summary$scale_value == summary$scale_value[[i]],
      ,
      drop = FALSE
    ]
    if (nrow(ref) == 1L && !is.na(ref$median_sec) &&
        !is.na(summary$median_sec[[i]]) && summary$median_sec[[i]] > 0) {
      summary$speedup_vs_lme4[[i]] <- ref$median_sec / summary$median_sec[[i]]
    }
  }

  summary[order(summary$scenario, summary$scale_value, summary$engine), , drop = FALSE]
}

plot_metric <- function(summary, metric, file, ylab, subtitle, log_y = FALSE,
                        hline = NULL) {
  dat <- summary[!is.na(summary[[metric]]), , drop = FALSE]
  if (!nrow(dat)) {
    return(FALSE)
  }

  dat$plot_y <- dat[[metric]]
  dat$engine <- factor(dat$engine, levels = c("mixeff", "lme4"))
  dat$scenario_label <- factor(
    paste0(scenario_display(dat$scenario), "\n", dat$random_terms),
    levels = unique(paste0(scenario_display(dat$scenario), "\n", dat$random_terms))
  )
  line_dat <- dat[
    stats::ave(dat$scale_value, dat$scenario_label, dat$engine, FUN = length) > 1L,
    ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = scale_value,
      y = plot_y,
      color = engine,
      shape = engine,
      group = engine
    )
  ) +
    ggplot2::geom_point(size = 2.4, stroke = 0.7, na.rm = TRUE) +
    ggplot2::facet_wrap(
      stats::as.formula("~ scenario_label"),
      scales = "free_x",
      ncol = 2
    ) +
    ggplot2::scale_color_manual(
      values = c(mixeff = "#0072B2", lme4 = "#D55E00"),
      breaks = c("mixeff", "lme4")
    ) +
    ggplot2::scale_shape_manual(
      values = c(mixeff = 16, lme4 = 17),
      breaks = c("mixeff", "lme4")
    ) +
    ggplot2::labs(
      title = "mixeff vs lme4 scaling benchmark",
      subtitle = subtitle,
      x = "Scale value (rows, grouping levels, or crossed subject/item levels)",
      y = ylab,
      color = NULL,
      shape = NULL,
      caption = "Synthetic Gaussian LMMs fit with REML = FALSE. Points are median elapsed time over completed repetitions."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(color = "grey25", margin = ggplot2::margin(b = 8)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey88", linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_line(color = "grey88", linewidth = 0.3),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      legend.position = "top",
      plot.caption = ggplot2::element_text(color = "grey35", hjust = 0)
    )

  if (nrow(line_dat)) {
    p <- p + ggplot2::geom_line(
      data = line_dat,
      linewidth = 0.8,
      alpha = 0.9,
      na.rm = TRUE
    )
  }

  if (isTRUE(log_y)) {
    p <- p + ggplot2::scale_y_log10()
  }
  if (!is.null(hline)) {
    p <- p + ggplot2::geom_hline(
      yintercept = hline,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.45
    )
  }

  grDevices::pdf(file, width = 11, height = 8.5, onefile = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  TRUE
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(0L))
  }

  require_package("mixeff")
  require_package("lme4")

  reps <- parse_int_scalar(value_or_env(args, "reps", "MIXEFF_BENCH_REPS", "3"), "reps")
  warmup <- parse_int_scalar(value_or_env(args, "warmup", "MIXEFF_BENCH_WARMUP", "1"),
                             "warmup")
  rows <- parse_int_list(value_or_env(args, "rows", "MIXEFF_BENCH_ROWS",
                                      "500,1000,2500,5000"), "rows")
  groups <- parse_int_list(value_or_env(args, "groups", "MIXEFF_BENCH_GROUPS",
                                        "25,50,100,200"), "groups")
  crossed_levels <- parse_int_list(
    value_or_env(args, "crossed_levels", "MIXEFF_BENCH_CROSSED_LEVELS",
                 "10,15,20,30"),
    "crossed-levels"
  )
  scenarios <- parse_scenarios(value_or_env(args, "scenarios", "MIXEFF_BENCH_SCENARIOS",
                                            "all"))
  out_dir <- value_or_env(args, "out", "MIXEFF_BENCH_OUT",
                          file.path("benchmarks", "lme4-scaling"))
  plots <- !identical(args$plots, FALSE)
  if (plots) {
    require_package("ggplot2")
  }
  engines <- c("mixeff", "lme4")

  if (reps < 1L) {
    stop("`reps` must be at least 1.", call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  specs <- make_specs(rows, groups, crossed_levels, scenarios)
  raw_rows <- vector("list", length(specs) * length(engines) * reps)
  idx <- 0L

  cat("Benchmarking ", length(specs), " designs; ", reps,
      " repetitions per engine/design.\n", sep = "")
  cat("Output directory: ", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")

  for (spec in specs) {
    cat("\n", spec$scenario, " scale=", spec$scale_value,
        " rows=", nrow(spec$data),
        " random_effect_coefficients=", spec$n_random_effect_coefficients,
        "\n", sep = "")
    warmup_fits(spec, engines, warmup)
    for (engine in engines) {
      for (rep_id in seq_len(reps)) {
        idx <- idx + 1L
        cat("  ", engine, " rep ", rep_id, "... ", sep = "")
        row <- time_fit(engine, spec, rep_id)
        raw_rows[[idx]] <- row
        if (isTRUE(row$ok)) {
          cat(sprintf("%.4f sec\n", row$elapsed_sec))
        } else {
          cat("failed: ", row$error, "\n", sep = "")
        }
      }
    }
  }

  raw <- do.call(rbind, raw_rows)
  summary <- summarise_timings(raw)

  raw_path <- file.path(out_dir, "lme4-scaling-raw.csv")
  summary_path <- file.path(out_dir, "lme4-scaling-summary.csv")
  utils::write.csv(raw, raw_path, row.names = FALSE)
  utils::write.csv(summary, summary_path, row.names = FALSE)

  cat("\nSummary\n")
  print(summary[, c("scenario", "scale_value", "n_rows", "engine", "median_sec",
                    "fits_per_sec", "speedup_vs_lme4", "n_ok", "n_fail")],
        row.names = FALSE)
  cat("\nWrote:\n  ", raw_path, "\n  ", summary_path, "\n", sep = "")

  if (plots) {
    speed_path <- file.path(out_dir, "lme4-scaling-fits-per-sec.pdf")
    time_path <- file.path(out_dir, "lme4-scaling-median-sec.pdf")
    speedup_path <- file.path(out_dir, "lme4-scaling-speedup-vs-lme4.pdf")
    plot_metric(
      summary,
      "fits_per_sec",
      speed_path,
      "Fits per second (higher is faster)",
      "Throughput curves by design complexity."
    )
    plot_metric(summary, "median_sec", time_path,
                "Median elapsed seconds (lower is faster)",
                "Elapsed wall-clock time on a log scale.",
                log_y = TRUE)
    plot_metric(summary, "speedup_vs_lme4", speedup_path,
                "Speed relative to lme4 (>1 is faster)",
                "Relative speed uses lme4 median time as the baseline.",
                hline = 1)
    cat("  ", speed_path, "\n  ", time_path, "\n  ", speedup_path, "\n", sep = "")
  }

  invisible(summary)
}

main()
