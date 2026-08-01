#!/usr/bin/env Rscript

# Canonical release gate for mixeff.
#
# Runs every acceptance check in one place and writes a machine- and
# human-readable artifact. Intended to be run from a CLEAN checkout so the
# result is the reproducible acceptance evidence for an R-universe / CRAN
# submission.
#
#   Rscript tools/release-gate.R                      # full gate
#   Rscript tools/release-gate.R --fast               # skip as-cran + offline
#   Rscript tools/release-gate.R --release-candidate  # full gate + slow tiers
#                                                     # + skip audit + live
#                                                     # aphantasia vignette
#
# Gates, in order:
#   1. Vendor provenance   — license notice + vendored-snapshot drift are current
#   2. Build               — R CMD build (with vignettes) produces THE tarball;
#                            every later step (as-cran check, install, suites)
#                            runs against that one artifact, whose sha256 and
#                            build identity (commit, pin, toolchain) go in the
#                            report header
#   3. R CMD check --as-cran on the BUILT TARBALL (never on the source dir)
#   4. Install the tarball into a throwaway gate library
#   5. lint                — lintr::lint_package() == 0
#   6. Fast test suite     — all six MIXEFF_* tier flags cleared, but
#                            NOT_CRAN=true so skip_on_cran blocks DO run here
#   7. Error/UX battery    — structured error classes + message fragments
#   8. Schema + manifest   — wire-contract snapshot tests
#   9. Offline install     — no-network source install (CARGO_NET_OFFLINE)
#
# --release-candidate additionally sets MIXEFF_RUN_SLOW_PARITY=true and
# MIXEFF_RUN_APHANTASIA=true for the suite, renders the aphantasia vignette
# live (MIXEFF_RUN_APHANTASIA_VIGNETTE=true), and FAILS on any suite skip
# outside a short allowlist (missing Suggests / pandoc, and the deeper opt-in
# stress / joint-proof tiers, which stay off even here).
#
# --fast is for iteration only; its report is stamped
# "FAST RUN -- NOT RELEASE EVIDENCE".

args <- commandArgs(trailingOnly = TRUE)
fast <- "--fast" %in% args
release_candidate <- "--release-candidate" %in% args
if (fast && release_candidate) {
  stop("--fast and --release-candidate are mutually exclusive")
}
gate_mode <- if (release_candidate) "release-candidate" else
  if (fast) "fast" else "full"
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

## Identity helpers ------------------------------------------------------------
sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  # tools::sha256sum arrived in R 4.5; fall back to the platform CLI on older R.
  sha_fun <- asNamespace("tools")$sha256sum
  if (is.function(sha_fun)) return(unname(sha_fun(path)))
  for (cmd in c("shasum", "sha256sum")) {
    if (!nzchar(Sys.which(cmd))) next
    cmd_args <- if (cmd == "shasum") c("-a", "256", shQuote(path)) else shQuote(path)
    out <- suppressWarnings(system2(cmd, cmd_args, stdout = TRUE, stderr = TRUE))
    if ((is.null(attr(out, "status")) || attr(out, "status") == 0) &&
        length(out)) {
      return(strsplit(trimws(out[[1L]]), "\\s+")[[1L]][[1L]])
    }
  }
  NA_character_
}

cmd_line1 <- function(cmd, cmd_args = "--version") {
  out <- tryCatch(
    suppressWarnings(system2(cmd, cmd_args, stdout = TRUE, stderr = TRUE)),
    error = function(e) character(0)
  )
  ok <- length(out) && (is.null(attr(out, "status")) || attr(out, "status") == 0)
  if (ok) trimws(out[[1L]]) else "unavailable"
}

or_unknown <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) "unknown" else x
}

## 1. Vendor provenance --------------------------------------------------------
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

## 2. Build the release tarball ------------------------------------------------
# ONE immutable artifact: the tarball built here is what gets checked with
# --as-cran, installed into the gate library, and exercised by every suite
# below — the checked artifact and the tested artifact are the same build.
desc <- read.dcf(file.path(pkg_root, "DESCRIPTION"))
pkg <- unname(desc[1L, "Package"])
version <- unname(desc[1L, "Version"])
r_bin <- file.path(R.home("bin"), "R")
# Hermetic CRAN-like env: NOT_CRAN= (empty) is the only value that reproduces
# CRAN here, because tools/config.R tests nzchar(), so even NOT_CRAN=false
# would count as "not CRAN" and drop cargo --offline.
check_env <- c("LANG=C", "LC_ALL=C", "LC_CTYPE=C", "NOT_CRAN=",
               "_R_CHECK_SYSTEM_CLOCK_=0")
# Keep the build/check dir under the package root (not tempdir): on a FAIL the
# .Rcheck directory IS the acceptance evidence, and the session tempdir is
# wiped on exit. .Rbuildignore/.gitignore exclude it.
check_root <- file.path(pkg_root, ".release-gate")
unlink(check_root, recursive = TRUE)
dir.create(check_root)
tarball <- file.path(check_root, sprintf("%s_%s.tar.gz", pkg, version))
build_ok <- FALSE
run("R CMD build (release tarball)", {
  old_wd <- setwd(check_root)
  tryCatch({
    build_out <- system2(
      r_bin,
      c("CMD", "build", shQuote(pkg_root)),
      stdout = TRUE, stderr = TRUE, env = check_env
    )
    build_status <- attr(build_out, "status")
    if (!is.null(build_status) && build_status != 0L) {
      stop(paste(tail(build_out, 20L), collapse = "\n"))
    }
    if (!file.exists(tarball)) {
      stop("R CMD build succeeded but did not produce ", basename(tarball))
    }
    build_ok <<- TRUE
    record("R CMD build (release tarball)", TRUE, basename(tarball))
  }, finally = setwd(old_wd))
})

## Build identity (recorded in the report header) ------------------------------
git_commit <- or_unknown(cmd_line1("git", c("-C", pkg_root, "rev-parse", "HEAD")))
git_dirty <- length(suppressWarnings(
  system2("git", c("-C", pkg_root, "status", "--porcelain"), stdout = TRUE)
)) > 0L
rust_pin <- local({
  lines <- readLines(file.path(pkg_root, "tools", "vendor-rust.R"), warn = FALSE)
  hit <- grep("^PINNED_REV\\s*<-", lines, value = TRUE)
  if (length(hit)) {
    return(sub('^PINNED_REV\\s*<-\\s*"([^"]+)".*$', "\\1", hit[[1L]]))
  }
  lock <- file.path(pkg_root, "src", "rust", "upstream", "mixeff-rs.lock")
  if (file.exists(lock)) {
    hit <- grep("^resolved_sha\\s*=", readLines(lock, warn = FALSE), value = TRUE)
    if (length(hit)) {
      return(trimws(sub('^resolved_sha\\s*=\\s*"?([^"#]+)"?.*$', "\\1", hit[[1L]])))
    }
  }
  NA_character_
})
tarball_sha256 <- if (build_ok) sha256_file(tarball) else NA_character_
tarball_bytes <- if (build_ok) file.size(tarball) else NA_real_
vendor_tar_sha256 <- sha256_file(file.path(pkg_root, "src", "rust", "vendor.tar.xz"))
cargo_lock_sha256 <- sha256_file(
  file.path(pkg_root, "src", "rust", "upstream", "mixeff-rs", "Cargo.lock")
)
identity_lines <- c(
  sprintf("mode: %s", gate_mode),
  sprintf("commit: %s%s", git_commit,
          if (git_dirty) "  [DIRTY WORKING TREE]" else ""),
  sprintf("mixeff-rs pin: %s", or_unknown(rust_pin)),
  sprintf("tarball: %s (%s bytes)",
          if (build_ok) basename(tarball) else "not built",
          if (is.na(tarball_bytes)) "?" else
            format(tarball_bytes, scientific = FALSE, trim = TRUE)),
  sprintf("tarball sha256: %s", or_unknown(tarball_sha256)),
  sprintf("src/rust/vendor.tar.xz sha256: %s", or_unknown(vendor_tar_sha256)),
  sprintf("upstream Cargo.lock sha256: %s", or_unknown(cargo_lock_sha256)),
  sprintf("R: %s (%s)", R.version.string, R.version$platform),
  sprintf("rustc: %s | cargo: %s", cmd_line1("rustc"), cmd_line1("cargo"))
)

## cran-comments.md tarball size (WI-6.4) --------------------------------------
# cran-comments.md quotes the exact "Size of tarball" line from the incoming
# NOTE; a stale figure means the recorded evidence describes some other
# tarball. The gate owns that number: refresh it from the tarball just built.
if (build_ok) run("cran-comments.md tarball size current", {
  cc_path <- file.path(pkg_root, "cran-comments.md")
  size_str <- format(tarball_bytes, scientific = FALSE, trim = TRUE)
  cc <- if (file.exists(cc_path)) readLines(cc_path, warn = FALSE) else character(0)
  idx <- grep("Size of tarball: [0-9]+ bytes", cc)
  if (!length(idx)) {
    record("cran-comments.md tarball size current", FALSE,
           "no 'Size of tarball: NNN bytes' line found in cran-comments.md")
  } else {
    old <- sub(".*Size of tarball: ([0-9]+) bytes.*", "\\1", cc[idx[[1L]]])
    if (identical(old, size_str)) {
      record("cran-comments.md tarball size current", TRUE,
             sprintf("%s bytes", size_str))
    } else {
      cc[idx] <- gsub("Size of tarball: [0-9]+ bytes",
                      sprintf("Size of tarball: %s bytes", size_str), cc[idx])
      writeLines(cc, cc_path)
      record("cran-comments.md tarball size current", TRUE,
             sprintf("updated %s -> %s bytes", old, size_str))
    }
  }
})

## 3. R CMD check --as-cran on the built tarball -------------------------------
if (fast) {
  record("R CMD check --as-cran (built tarball)", NA, "skipped (--fast)")
} else if (!build_ok) {
  record("R CMD check --as-cran (built tarball)", FALSE, "skipped: build failed")
} else {
  run("R CMD check --as-cran (built tarball)", {
    old_wd <- setwd(check_root)
    tryCatch({
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
}

## 4. Install the built tarball ------------------------------------------------
# The gate exercises the INSTALLED package (a release build; devtools::test()
# / load_all() is a debug build ~60x slower and not representative) and
# installs the TARBALL, not the source dir, so the artifact checked above is
# the artifact the suites run against. Install BEFORE linting:
# lintr::object_usage_linter resolves internal functions from the installed
# namespace, so linting an uninstalled package yields false "no visible
# global function" positives.
suite_gate <- if (release_candidate) {
  "suite green (slow-parity + aphantasia tiers ON)"
} else {
  "fast suite green"
}
gate_lib <- tempfile("mixeff-release-gate-lib-")
dir.create(gate_lib)
install_ok <- FALSE
if (!build_ok) {
  record("install built tarball", FALSE, "skipped: build failed")
} else {
  run("install built tarball", {
    install_out <- system2(r_bin,
                           c("CMD", "INSTALL", "--preclean",
                             shQuote(paste0("--library=", gate_lib)),
                             shQuote(tarball)),
                           stdout = TRUE, stderr = TRUE)
    st <- attr(install_out, "status")
    install_ok <<- is.null(st) || st == 0L
    record("install built tarball", install_ok,
           if (!install_ok) paste(tail(install_out, 3), collapse = " ") else "")
  })
}
if (install_ok) {
  # Only prepend the gate library when it actually holds the fresh build;
  # otherwise test_dir would silently resolve a stale mixeff from the user
  # library and the artifact would report on the wrong binary.
  .libPaths(c(gate_lib, .libPaths()))
} else {
  record("lint_package == 0", FALSE, "skipped: install failed")
  record(suite_gate, FALSE, "skipped: install failed")
}

## 5. lint ---------------------------------------------------------------------
if (install_ok) run("lint_package == 0", {
  lints <- lintr::lint_package(pkg_root)
  record("lint_package == 0", length(lints) == 0,
         sprintf("%d lint(s)", length(lints)))
})

## 6-8. Test suites (installed release build) ----------------------------------
# All six tier flags are cleared first so an operator's exported shell vars
# cannot silently change what the gate runs; --release-candidate then turns
# the two release tiers back on explicitly.
mixeff_tier_flags <- c(
  "MIXEFF_RUN_SLOW_PARITY", "MIXEFF_RUN_APHANTASIA", "MIXEFF_APHANTASIA_JOINT",
  "MIXEFF_RUN_APHANTASIA_STRESS", "MIXEFF_APHANTASIA_JOINT_PROOF",
  "MIXEFF_APHANTASIA_JOINT_BUDGET"
)
# Skips acceptable in a release-candidate run. testthat prefixes reasons with
# "Reason: "; matching is fixed-substring. Anything else fails the gate.
rc_skip_allowlist <- c(
  "is not installed",                       # skip_if_not_installed() Suggests
  "pandoc is unavailable",                  # rmarkdown::pandoc_available()
  "Set MIXEFF_RUN_APHANTASIA_STRESS=true",  # deeper opt-in tier, off here
  "Set MIXEFF_APHANTASIA_JOINT_PROOF=true", # deeper opt-in tier, off here
  "empty test"
)
if (install_ok) run("suite + error-UX + schema/manifest", {
  Sys.unsetenv(mixeff_tier_flags)
  if (release_candidate) {
    Sys.setenv(MIXEFF_RUN_SLOW_PARITY = "true", MIXEFF_RUN_APHANTASIA = "true")
  }
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
  if (release_candidate) {
    Sys.unsetenv(c("MIXEFF_RUN_SLOW_PARITY", "MIXEFF_RUN_APHANTASIA"))
  }
  df <- as.data.frame(res)
  fails <- sum(df$failed) + sum(df$error)
  suite_elapsed <- round(proc.time()["elapsed"] - suite_t)
  record(suite_gate, fails == 0,
         sprintf("%d pass / %d fail / %d skip in %ds",
                 sum(df$passed), fails, sum(df$skipped), suite_elapsed))
  # Error-UX and schema/manifest are part of that run; assert they were present.
  ux_present <- any(grepl("error-ux", df$file))
  schema_present <- any(grepl("schema|manifest", df$file))
  record("error-UX battery present", ux_present)
  record("schema/manifest tests present", schema_present)
  if (release_candidate) {
    skip_msgs <- unlist(lapply(res, function(t) {
      vapply(Filter(function(x) inherits(x, "expectation_skip"), t$results),
             conditionMessage, character(1))
    }))
    if (is.null(skip_msgs)) skip_msgs <- character(0)
    allowed <- vapply(skip_msgs, function(m) {
      any(vapply(rc_skip_allowlist, grepl, logical(1), x = m, fixed = TRUE))
    }, logical(1))
    unexpected <- skip_msgs[!allowed]
    if (length(unexpected)) {
      cat("Skips outside the release-candidate allowlist:\n")
      cat(paste0("  - ", unique(unexpected), "\n"), sep = "")
    }
    record("release-candidate: skips all allowlisted",
           length(unexpected) == 0L,
           sprintf("%d skip(s), %d outside allowlist",
                   length(skip_msgs), length(unexpected)))
  }
})

## release-candidate: aphantasia vignette rendered live ------------------------
if (release_candidate) {
  if (!install_ok) {
    record("aphantasia vignette (live)", FALSE, "skipped: install failed")
  } else run("aphantasia vignette (live)", {
    vig <- file.path(pkg_root, "vignettes", "reproducing-aphantasia.Rmd")
    if (!requireNamespace("rmarkdown", quietly = TRUE) ||
        !rmarkdown::pandoc_available()) {
      record("aphantasia vignette (live)", FALSE, "rmarkdown/pandoc unavailable")
    } else {
      out_dir <- file.path(check_root, "vignette-live")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      old_flag <- Sys.getenv("MIXEFF_RUN_APHANTASIA_VIGNETTE",
                             unset = NA_character_)
      Sys.setenv(MIXEFF_RUN_APHANTASIA_VIGNETTE = "true")
      # Restore the env var before re-raising: this block is a top-level
      # promise, so on.exit would not fire (same caveat as the suite step).
      rendered <- tryCatch(
        rmarkdown::render(vig, output_dir = out_dir,
                          envir = new.env(parent = globalenv()), quiet = TRUE),
        error = function(e) e
      )
      if (is.na(old_flag)) Sys.unsetenv("MIXEFF_RUN_APHANTASIA_VIGNETTE") else
        Sys.setenv(MIXEFF_RUN_APHANTASIA_VIGNETTE = old_flag)
      if (inherits(rendered, "error")) stop(conditionMessage(rendered))
      record("aphantasia vignette (live)", file.exists(rendered),
             basename(rendered))
    }
  })
}

## 9. Offline install ----------------------------------------------------------
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

## Summary + artifact -----------------------------------------------------------
elapsed <- round(proc.time()["elapsed"] - started)
n_fail <- sum(vapply(results, function(r) isFALSE(r$ok), logical(1)))
n_pass <- sum(vapply(results, function(r) isTRUE(r$ok), logical(1)))
n_skip <- sum(vapply(results, function(r) is.na(r$ok), logical(1)))

fast_banner <-
  "FAST RUN -- NOT RELEASE EVIDENCE (as-cran and offline-install gates skipped)"
lines <- c(
  sprintf("# mixeff release gate -- %s", version),
  identity_lines,
  if (fast) fast_banner,
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
if (fast) cat("\n", fast_banner, "\n", sep = "")
cat(sprintf("\nAcceptance artifact: %s\n", artifact))
if (n_fail > 0) quit(status = 1L)
