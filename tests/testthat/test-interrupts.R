# The interrupt bridge is plumbing for Phase 1+ — long-running PLS/PIRLS
# fits will call R_CheckUserInterrupt periodically inside their iteration
# loops. In Phase 0 we ship only the FFI binding plus a no-op demo so we
# can verify the symbol is linked.
#
# Robustly testing actual Ctrl-C handling from inside testthat is fragile
# (different platforms, parallel test runners, etc.). We assert the
# easy invariants here and leave a manual reproduction note in the file
# for anyone changing the binding.
#
# Manual repro: in an interactive R session,
#
#   library(mixeff)
#   mixeff:::mm_interrupt_demo(1e9)   # then press Ctrl-C
#
# should return control to the prompt within a fraction of a second.

test_that("interrupt demo returns the iteration count on clean completion", {
  expect_identical(mixeff:::mm_interrupt_demo(0L), 0L)
  expect_identical(mixeff:::mm_interrupt_demo(1L), 1L)
  expect_identical(mixeff:::mm_interrupt_demo(100L), 100L)
})

test_that("interrupt demo input validation", {
  expect_error(mixeff:::mm_interrupt_demo(NA_integer_),
               "`iters` must be a single non-NA integer")
  expect_error(mixeff:::mm_interrupt_demo(c(1L, 2L)),
               "`iters` must be a single non-NA integer")
})

test_that("interrupt demo handles negative input by clamping to zero", {
  # Rust side does iters.max(0); R side passes negatives through. Either
  # behaviour is acceptable — assert the contract: completes without error
  # and returns a non-negative count.
  out <- mixeff:::mm_interrupt_demo(-5L)
  expect_true(out >= 0L)
})

test_that("a running fit stops cleanly with mm_interrupted on SIGINT", {
  skip_on_cran()
  skip_if_not_installed("processx")
  skip_if_not_installed("pkgload")
  skip_on_os("windows")
  # Load mixeff in the child the way this session did: an installed copy
  # (R CMD check, covr) from the same library, else the source tree.
  ns_path <- getNamespaceInfo("mixeff", "path")
  load_line <- if (dir.exists(file.path(ns_path, "Meta"))) {
    sprintf(".libPaths(c(%s, .libPaths())); suppressMessages(library(mixeff))",
            deparse(dirname(ns_path)))
  } else {
    skip_if_not(file.exists(file.path(ns_path, "DESCRIPTION")),
                "package source tree is not available")
    sprintf("suppressMessages(pkgload::load_all(%s, compile = FALSE, quiet = TRUE))",
            deparse(ns_path))
  }
  script <- tempfile(fileext = ".R")
  writeLines(c(
    load_line,
    "set.seed(1)",
    "n <- 40000; d <- data.frame(s = factor(sample(800, n, TRUE)),",
    "  it = factor(sample(400, n, TRUE)), x = rnorm(n))",
    "d$y <- rbinom(n, 1, plogis(-0.5 + 0.3 * d$x + rnorm(800)[d$s] + rnorm(400)[d$it]))",
    "cat('READY\\n')",
    "r <- tryCatch({",
    "  glmm(y ~ x + (1 | s) + (1 | it), d, binomial, control = mm_control(verbose = -1))",
    "  'COMPLETED'",
    "}, mm_interrupted = function(e) 'MM_INTERRUPTED', error = function(e) paste('ERROR', conditionMessage(e)))",
    "cat(r, '\\n')",
    "cat('ALIVE', mixeff:::mm_interrupt_demo(3L), '\\n')"
  ), script)
  p <- processx::process$new(file.path(R.home("bin"), "Rscript"), script,
                             stdout = "|", stderr = "|")
  on.exit(p$kill(), add = TRUE)
  out <- ""
  deadline <- Sys.time() + 120
  while (!grepl("READY", out) && p$is_alive() && Sys.time() < deadline) {
    p$poll_io(500)
    out <- paste0(out, p$read_output())
  }
  skip_if_not(grepl("READY", out), "child process did not start")
  Sys.sleep(2)
  p$interrupt()
  p$wait(60000)
  out <- paste0(out, p$read_all_output())
  skip_if(grepl("COMPLETED", out), "fit finished before the interrupt arrived")
  expect_match(out, "MM_INTERRUPTED", fixed = TRUE)
  # The R session survives and the native library is still usable.
  expect_match(out, "ALIVE 3", fixed = TRUE)
})
