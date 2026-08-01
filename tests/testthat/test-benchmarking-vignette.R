test_that("benchmarking vignette stays light and points at benchmark scripts", {
  skip_on_cran() # duplicates R CMD check's own vignette build
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is unavailable")

  path <- test_path("../../vignettes/benchmarking.Rmd")
  # Source-tree-only test: the vignette source is not on the installed
  # tests path during R CMD check, only during devtools::test().
  skip_if_not(file.exists(path), "vignette source path is not available")
  text <- readLines(path, warn = FALSE)
  expect_true(any(grepl("inst/benchmarks/lme4-scaling.R", text, fixed = TRUE)))
  expect_true(any(grepl("inst/benchmarks/bootstrap-inference.R", text, fixed = TRUE)))
  expect_true(any(grepl("eval = FALSE", text, fixed = TRUE)))

  out <- rmarkdown::render(path, output_file = tempfile(fileext = ".html"),
                           quiet = TRUE)
  expect_true(file.exists(out))
})

test_that("the shipped calibration artifact is release-grade with honest accounting", {
  # Audit1 WI-3.1 (DEC-2): the artifact behind the vignette's calibration
  # section must be release-mode with the agreed replication floor, carry
  # Wilson intervals, close its replicate accounting, and ship with a
  # provenance manifest whose invocation reproduces it. A smoke-mode or
  # under-replicated CSV must fail here, not render as evidence.
  csv <- system.file("extdata", "inference-method-simulation-summary.csv",
                     package = "mixeff")
  man <- system.file("extdata", "inference-method-simulation-manifest.json",
                     package = "mixeff")
  expect_true(nzchar(csv))
  expect_true(nzchar(man))

  d <- read.csv(csv, stringsAsFactors = FALSE)
  expect_true(all(d$mode == "release"))
  expect_true(all(c("n_requested", "n_evaluated_null", "n_refused_null",
                    "n_failed_null", "n_fit_failed_null", "n_boundary_null",
                    "type_I_lower", "type_I_upper", "power_lower",
                    "power_upper", "coverage_lower", "coverage_upper",
                    "type_I_unconditional", "reason") %in% names(d)))

  ran <- d[d$n_requested > 0, ]
  expect_true(nrow(ran) > 0L)
  # Replication floors: analytic >= 2500, bootstrap >= 500 (DEC-2).
  boot_methods <- c("bootstrap", "bootstrap_lrt", "glmm_parametric_bootstrap")
  joint_methods <- "glmm_wald_joint_laplace"
  analytic <- ran[!(ran$method %in% c(boot_methods, joint_methods)), ]
  bootstrapped <- ran[ran$method %in% boot_methods, ]
  expect_true(all(analytic$n_requested >= 2500L))
  expect_true(all(bootstrapped$n_requested >= 500L))
  expect_true(all(ran[ran$method %in% joint_methods, "n_requested"] >= 1000L))
  # Accounting closes exactly; nothing dropped.
  expect_true(all(ran$n_evaluated_null + ran$n_refused_null +
                    ran$n_failed_null + ran$n_fit_failed_null ==
                    ran$n_requested))
  # Wilson bounds bracket their point estimates.
  fin <- is.finite(ran$type_I_error)
  expect_true(all(ran$type_I_lower[fin] <= ran$type_I_error[fin] &
                    ran$type_I_error[fin] <= ran$type_I_upper[fin]))
  # Rows that did not run say why.
  skipped <- d[d$n_requested == 0, ]
  expect_true(all(nzchar(skipped$reason) & !is.na(skipped$reason)))

  m <- jsonlite::read_json(man)
  expect_identical(m$mode, "release")
  expect_true(m$reps_analytic >= 2500L)
  expect_true(m$reps_bootstrap >= 500L)
  expect_true(grepl("--mode=release", m$invocation, fixed = TRUE))
  expect_true(nzchar(m$package_commit %||% ""))
  expect_true(nzchar(m$rust_pin %||% ""))
})
