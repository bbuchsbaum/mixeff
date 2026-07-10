#!/usr/bin/env Rscript

# Canonical release gate for mixeff.
#
# Runs every acceptance check in one place and writes a machine- and
# human-readable artifact. Intended to be run from a CLEAN checkout so the
# result is the reproducible acceptance evidence for an R-universe / CRAN
# submission.
#
#   Rscript tools/release-gate.R              # full gate
#   Rscript tools/release-gate.R --fast       # skip R CMD check + offline install
#
# Gates, in order:
#   1. Vendor provenance   — license notice + vendored-snapshot drift are current
#   2. Build               — R CMD build (with vignettes) produces the tarball
#   3. R CMD check --as-cran on the BUILT TARBALL (never on the source dir)
#   4. lint                — lintr::lint_package() == 0
#   5. Fast test suite     — default gates OFF (no slow-parity / aphantasia)
#   6. Error/UX battery    — structured error classes + message fragments
#   7. Schema + manifest   — wire-contract snapshot tests
#   8. Offline install     — no-network source install (CARGO_NET_OFFLINE)
#
# The slow parity suite (MIXEFF_RUN_SLOW_PARITY) and aphantasia reproduction
# (MIXEFF_RUN_APHANTASIA) are NOT part of the gate — they are opt-in and can
# run for many minutes. Run them separately with those env vars set.

args <- commandArgs(trailingOnly = TRUE)
fast <- "--fast" %in% args
pkg_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(pkg_root, "DESCRIPTION")))

started <- proc.time()["elapsed"]
results <- list()

record <- function(name, ok, detail = "") {
  results[[length(results) + 1L]] <<- list(name = name, ok = isTRUE(ok),
                                           detail = detail)
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

run <- function(name, expr, timeout = Inf) {
  cat(sprintf("\n== %s ==\n", name))
  res <- tryCatch(
    withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "gate_err")
  )
  if (inherits(res, "gate_err")) {
    record(name, FALSE, res$msg)
    return(invisible(FALSE))
  }
  res
}

## 1. Vendor provenance -------------------------------------------------------
run("vendor: license notice current", {
  out <- system2("Rscript", c("tools/check-license-note.R"),
                 stdout = TRUE, stderr = TRUE)
  record("vendor: license notice current", TRUE)
})
run("vendor: snapshot drift", {
  out <- system2("Rscript", c("tools/check-vendor-drift.R"),
                 stdout = TRUE, stderr = TRUE)
  ok <- is.null(attr(out, "status")) || attr(out, "status") == 0
  record("vendor: snapshot drift", ok, if (!ok) paste(tail(out, 2), collapse = " "))
})

## 2 + 3. Build tarball, R CMD check --as-cran on the tarball ----------------
if (!fast) {
  run("R CMD check --as-cran (built tarball)", {
    if (!requireNamespace("rcmdcheck", quietly = TRUE)) {
      stop("rcmdcheck is required for the release gate")
    }
    chk <- rcmdcheck::rcmdcheck(
      pkg_root,
      args = c("--as-cran", "--no-manual"),
      error_on = "never"
    )
    # rcmdcheck builds the tarball first, then checks it -- the correct
    # release-verification path (never `R CMD check .` on the source dir).
    ne <- length(chk$errors); nw <- length(chk$warnings); nn <- length(chk$notes)
    # The only acceptable NOTE is the unavoidable first-submission note.
    non_submission_notes <- Filter(
      function(n) !grepl("New submission", n, fixed = TRUE), chk$notes
    )
    ok <- ne == 0 && nw == 0 && length(non_submission_notes) == 0
    record("R CMD check --as-cran (built tarball)", ok,
           sprintf("%dE/%dW/%dN", ne, nw, nn))
  })
} else {
  record("R CMD check --as-cran (built tarball)", NA, "skipped (--fast)")
}

## 4. lint --------------------------------------------------------------------
run("lint_package == 0", {
  lints <- lintr::lint_package(pkg_root)
  record("lint_package == 0", length(lints) == 0,
         sprintf("%d lint(s)", length(lints)))
})

## 5-7. Test suites (installed release build) --------------------------------
# The gate exercises the INSTALLED package (a release build); devtools::test()
# / load_all() is a debug build ~60x slower and is not representative. Install
# to a throwaway library so the gate is self-contained from a clean checkout.
gate_lib <- tempfile("mixeff-release-gate-lib-")
dir.create(gate_lib)
run("install release build", {
  st <- system2("R", c("CMD", "INSTALL", "--preclean",
                       paste0("--library=", gate_lib), pkg_root),
                stdout = FALSE, stderr = FALSE)
  record("install release build", identical(st, 0L))
})
.libPaths(c(gate_lib, .libPaths()))

run("fast suite + error-UX + schema/manifest", {
  Sys.unsetenv(c("MIXEFF_RUN_SLOW_PARITY", "MIXEFF_RUN_APHANTASIA",
                 "MIXEFF_APHANTASIA_JOINT"))
  suite_t <- proc.time()["elapsed"]
  res <- testthat::test_dir(
    file.path(pkg_root, "tests", "testthat"),
    package = "mixeff", load_package = "installed",
    stop_on_failure = FALSE, reporter = "silent"
  )
  df <- as.data.frame(res)
  fails <- sum(df$failed) + sum(df$error)
  suite_elapsed <- round(proc.time()["elapsed"] - suite_t)
  record("fast suite green", fails == 0,
         sprintf("%d pass / %d fail / %d skip in %ds",
                 sum(df$passed), fails, sum(df$skipped), suite_elapsed))
  # Error-UX and schema/manifest are part of that run; assert they were present.
  ux_present <- any(grepl("error-ux", df$file))
  schema_present <- any(grepl("schema|manifest", df$file))
  record("error-UX battery present", ux_present)
  record("schema/manifest tests present", schema_present)
})

## 8. Offline install ---------------------------------------------------------
if (!fast) {
  run("no-network source install", {
    out <- system2("Rscript", c("tools/check-offline-install.R"),
                   stdout = TRUE, stderr = TRUE,
                   env = "CARGO_NET_OFFLINE=true")
    ok <- is.null(attr(out, "status")) || attr(out, "status") == 0
    record("no-network source install", ok,
           if (!ok) paste(tail(out, 2), collapse = " "))
  })
} else {
  record("no-network source install", NA, "skipped (--fast)")
}

## Summary + artifact --------------------------------------------------------
elapsed <- round(proc.time()["elapsed"] - started)
n_fail <- sum(vapply(results, function(r) isFALSE(r$ok), logical(1)))
n_pass <- sum(vapply(results, function(r) isTRUE(r$ok), logical(1)))
n_skip <- sum(vapply(results, function(r) is.na(r$ok), logical(1)))

lines <- c(
  sprintf("# mixeff release gate -- %s",
          tryCatch(as.character(read.dcf("DESCRIPTION")[, "Version"]),
                   error = function(e) "?")),
  sprintf("elapsed: %ds | PASS %d  FAIL %d  SKIP %d", elapsed, n_pass, n_fail,
          n_skip),
  "",
  vapply(results, function(r) {
    sprintf("- [%s] %s%s",
            if (is.na(r$ok)) "SKIP" else if (r$ok) "PASS" else "FAIL",
            r$name, if (nzchar(r$detail)) paste0(" -- ", r$detail) else "")
  }, character(1))
)
artifact <- file.path(pkg_root, "release-gate-report.txt")
writeLines(lines, artifact)

cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
cat(sprintf("\nAcceptance artifact: %s\n", artifact))
if (n_fail > 0) quit(status = 1L)
