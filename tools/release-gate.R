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
#   5. Fast test suite     — slow-parity / aphantasia gates OFF, but
#                            NOT_CRAN=true so skip_on_cran blocks DO run here
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
  detail <- paste(as.character(detail), collapse = " ")
  if (is.na(detail)) detail <- ""
  normalized_ok <- if (length(ok) == 1L && is.na(ok)) NA else isTRUE(ok)
  results[[length(results) + 1L]] <<- list(name = name, ok = normalized_ok,
                                           detail = detail)
  label <- if (is.na(normalized_ok)) {
    "SKIP"
  } else if (normalized_ok) {
    "PASS"
  } else {
    "FAIL"
  }
  cat(sprintf("[%s] %s%s\n", label, name,
              if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

run <- function(name, expr) {
  cat(sprintf("\n== %s ==\n", name))
  # Belt and braces: expr evaluates at top level, so a step that changes the
  # working directory and then errors could leak its cwd into later steps.
  wd <- getwd()
  on.exit(setwd(wd), add = TRUE)
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
  out <- system2("Rscript", file.path(pkg_root, "tools", "check-license-note.R"),
                 stdout = TRUE, stderr = TRUE)
  ok <- is.null(attr(out, "status")) || attr(out, "status") == 0
  record("vendor: license notice current", ok,
         if (!ok) paste(tail(out, 2), collapse = " ") else "")
})
run("vendor: snapshot drift", {
  out <- system2("Rscript", c("tools/check-vendor-drift.R"),
                 stdout = TRUE, stderr = TRUE)
  ok <- is.null(attr(out, "status")) || attr(out, "status") == 0
  record("vendor: snapshot drift", ok,
         if (!ok) paste(tail(out, 2), collapse = " ") else "")
})

## 2 + 3. Build tarball, R CMD check --as-cran on the tarball ----------------
if (!fast) {
  run("R CMD check --as-cran (built tarball)", {
    # Use base R only. The acceptance gate must work on a clean machine with
    # the package's declared dependencies; requiring the optional rcmdcheck
    # wrapper would make the gate itself less portable than the package.
    desc <- read.dcf(file.path(pkg_root, "DESCRIPTION"))
    pkg <- unname(desc[1L, "Package"])
    version <- unname(desc[1L, "Version"])
    # Keep the check dir under the package root (not tempdir): on a FAIL the
    # .Rcheck directory IS the acceptance evidence, and the session tempdir
    # is wiped on exit. .Rbuildignore/.gitignore exclude it.
    check_root <- file.path(pkg_root, ".release-gate")
    unlink(check_root, recursive = TRUE)
    dir.create(check_root)
    old_wd <- setwd(check_root)
    tryCatch({
      r_bin <- file.path(R.home("bin"), "R")
      # Hermetic CRAN-like env: NOT_CRAN= (empty) is the only value that
      # reproduces CRAN here, because tools/config.R tests nzchar(), so even
      # NOT_CRAN=false would count as "not CRAN" and drop cargo --offline.
      check_env <- c("LANG=C", "LC_ALL=C", "LC_CTYPE=C", "NOT_CRAN=",
                     "_R_CHECK_SYSTEM_CLOCK_=0")

      build_out <- system2(
        r_bin,
        c("CMD", "build", shQuote(pkg_root)),
        stdout = TRUE, stderr = TRUE, env = check_env
      )
      build_status <- attr(build_out, "status")
      if (!is.null(build_status) && build_status != 0L) {
        stop(paste(tail(build_out, 20L), collapse = "\n"))
      }

      tarball <- file.path(check_root, sprintf("%s_%s.tar.gz", pkg, version))
      if (!file.exists(tarball)) {
        stop("R CMD build succeeded but did not produce ", basename(tarball))
      }

      check_out <- system2(
        r_bin,
        c("CMD", "check", "--as-cran", "--no-manual", shQuote(tarball)),
        stdout = TRUE, stderr = TRUE, env = check_env
      )
      check_status <- attr(check_out, "status")
      if (is.null(check_status)) check_status <- 0L
      log_path <- file.path(check_root, paste0(pkg, ".Rcheck"), "00check.log")
      if (!file.exists(log_path)) {
        stop(paste(
          "R CMD check did not produce 00check.log:",
          paste(tail(check_out, 20L), collapse = "\n"),
          sep = "\n"
        ))
      }

      log <- readLines(log_path, warn = FALSE)
      status_line <- tail(grep("^Status:", log, value = TRUE), 1L)
      if (!length(status_line)) {
        stop(paste(
          "00check.log has no Status line (check killed or truncated):",
          paste(tail(check_out, 20L), collapse = "\n"),
          sep = "\n"
        ))
      }
      status_count <- function(kind) {
        hit <- regmatches(
          status_line,
          regexec(sprintf("([0-9]+) %s", kind), status_line)
        )[[1L]]
        if (length(hit) >= 2L) as.integer(hit[[2L]]) else 0L
      }
      ne <- status_count("ERROR")
      nw <- status_count("WARNING")
      nn <- status_count("NOTE")

      # The only acceptable NOTE is the CRAN incoming-feasibility note whose
      # body contains nothing beyond "New submission" and the tarball-size
      # line (both documented in cran-comments.md). Matching "New submission"
      # anywhere in the log would also whitelist bundled URL/spelling/DOI
      # complaints, which CRAN bounces.
      submission_only <- FALSE
      if (nn == 1L) {
        idx <- grep("^\\* checking .*NOTE$", log)
        if (length(idx) == 1L &&
            grepl("CRAN incoming feasibility", log[idx])) {
          nxt <- grep("^\\* ", log)
          end <- suppressWarnings(min(nxt[nxt > idx]))
          if (!is.finite(end)) end <- length(log) + 1L
          body <- trimws(log[seq(idx + 1L, end - 1L)])
          body <- body[nzchar(body)]
          body <- body[!grepl("^Maintainer:", body)]
          body <- body[body != "New submission"]
          body <- body[!grepl("^Size of tarball:", body)]
          submission_only <- length(body) == 0L
        }
      }
      ok <- check_status == 0L && ne == 0L && nw == 0L &&
        (nn == 0L || submission_only)
      if (!ok) cat(paste(check_out, collapse = "\n"), "\n")
      record("R CMD check --as-cran (built tarball)", ok,
             sprintf("%dE/%dW/%dN%s", ne, nw, nn,
                     if (nn == 1L && !submission_only)
                       " (NOTE not submission-only)" else ""))
    }, finally = setwd(old_wd))
  })
} else {
  record("R CMD check --as-cran (built tarball)", NA, "skipped (--fast)")
}

## 4. Install release build ---------------------------------------------------
# The gate exercises the INSTALLED package (a release build); devtools::test()
# / load_all() is a debug build ~60x slower and is not representative. Install
# to a throwaway library so the gate is self-contained from a clean checkout.
# Install BEFORE linting: lintr::object_usage_linter resolves internal
# functions from the installed namespace, so linting an uninstalled package
# yields false "no visible global function" positives.
gate_lib <- tempfile("mixeff-release-gate-lib-")
dir.create(gate_lib)
install_ok <- FALSE
run("install release build", {
  install_out <- system2(file.path(R.home("bin"), "R"),
                         c("CMD", "INSTALL", "--preclean",
                           shQuote(paste0("--library=", gate_lib)),
                           shQuote(pkg_root)),
                         stdout = TRUE, stderr = TRUE)
  st <- attr(install_out, "status")
  install_ok <<- is.null(st) || st == 0L
  record("install release build", install_ok,
         if (!install_ok) paste(tail(install_out, 3), collapse = " ") else "")
})
if (install_ok) {
  # Only prepend the gate library when it actually holds the fresh build;
  # otherwise test_dir would silently resolve a stale mixeff from the user
  # library and the artifact would report on the wrong binary.
  .libPaths(c(gate_lib, .libPaths()))
} else {
  record("lint_package == 0", FALSE, "skipped: install failed")
  record("fast suite green", FALSE, "skipped: install failed")
}

## 5. lint --------------------------------------------------------------------
if (install_ok) run("lint_package == 0", {
  lints <- lintr::lint_package(pkg_root)
  record("lint_package == 0", length(lints) == 0,
         sprintf("%d lint(s)", length(lints)))
})

## 6-8. Test suites (installed release build) --------------------------------
if (install_ok) run("fast suite + error-UX + schema/manifest", {
  Sys.unsetenv(c("MIXEFF_RUN_SLOW_PARITY", "MIXEFF_RUN_APHANTASIA",
                 "MIXEFF_APHANTASIA_JOINT"))
  # The gate is the local acceptance run: skip_on_cran blocks (profile
  # likelihood, parity scoreboard) must execute here even though they are
  # skipped on CRAN's machines. Scope the override to this step: a leaked
  # NOT_CRAN=true reaches tools/config.R in the offline-install subprocess
  # below and silently drops cargo --offline.
  # No on.exit here: this block is a promise evaluated at top level, where
  # on.exit handlers do not fire. Restore explicitly below; the offline
  # install step also pins NOT_CRAN= in its own subprocess env in case this
  # step errors before the restore line.
  old_not_cran <- Sys.getenv("NOT_CRAN", unset = NA_character_)
  Sys.setenv(NOT_CRAN = "true")
  suite_t <- proc.time()["elapsed"]
  res <- testthat::test_dir(
    file.path(pkg_root, "tests", "testthat"),
    package = "mixeff", load_package = "installed",
    stop_on_failure = FALSE, reporter = "silent"
  )
  if (is.na(old_not_cran)) Sys.unsetenv("NOT_CRAN") else
    Sys.setenv(NOT_CRAN = old_not_cran)
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
    offline_check <- file.path(pkg_root, "tools", "check-offline-install.R")
    out <- system2("Rscript", offline_check,
                   stdout = TRUE, stderr = TRUE,
                   env = c("CARGO_NET_OFFLINE=true", "NOT_CRAN="))
    ok <- is.null(attr(out, "status")) || attr(out, "status") == 0
    record("no-network source install", ok,
           if (!ok) paste(tail(out, 2), collapse = " ") else "")
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
          tryCatch(as.character(
            read.dcf(file.path(pkg_root, "DESCRIPTION"))[, "Version"]
          ),
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
