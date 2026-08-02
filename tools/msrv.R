# Rust toolchain preflight for configure/configure.win (Audit1 WI-6.3).
#
# CRAN's Rust policy asks packages to check BOTH cargo and rustc, including
# their declared minimum versions, and to look under ~/.cargo/bin for
# personal installations. This script:
#   * parses the declared minima from DESCRIPTION SystemRequirements
#     (robust to commas inside parenthesized clauses -- the naive ", "
#     split fragmented "Cargo (Rust's package manager, >= 1.78.0)");
#   * extends PATH with the user's cargo bin dir using the PLATFORM path
#     separator (the previous literal ":" corrupted PATH on Windows);
#   * detects missing binaries via Sys.which()/exit status (system(intern=)
#     signals a warning, not an error, for status 127 -- the old
#     tryCatch(error=) handlers could never fire);
#   * compares BOTH detected versions against BOTH declared minima and
#     fails with a plain installation instruction (never installs Rust);
#   * prints both detected versions into the installation log.
#
# Sourced by tools/config.R at configure time. Everything is defined as a
# function and executed by the mm_msrv_main() call at the bottom, so tests
# can source this file with MIXEFF_MSRV_NO_RUN=true and unit-test the
# parsing/comparison logic directly.

mm_msrv_declared_minimum <- function(sysreqs, tool) {
  # Matches e.g. "Cargo (Rust's package manager, >= 1.78.0)" and
  # "rustc (>= 1.78.0)"; returns NA_character_ when no minimum is declared.
  pattern <- paste0(tool, "[^,()]*\\([^)]*>=\\s*([0-9][0-9.]*)\\)")
  m <- regexec(pattern, sysreqs, ignore.case = TRUE)
  hit <- regmatches(sysreqs, m)[[1L]]
  if (length(hit) >= 2L) hit[[2L]] else NA_character_
}

mm_msrv_extract_semver <- function(text) {
  if (length(text) != 1L || !nzchar(text) ||
      !grepl("\\d+\\.\\d+(\\.\\d+)?", text)) {
    return(NA_character_)
  }
  sub(".*?(\\d+\\.\\d+(\\.\\d+)?).*", "\\1", text)
}

mm_msrv_extend_path <- function() {
  sep <- .Platform$path.sep
  home <- Sys.getenv("HOME")
  candidates <- file.path(home, ".cargo", "bin")
  if (identical(.Platform$OS.type, "windows")) {
    profile <- Sys.getenv("USERPROFILE")
    if (nzchar(profile)) {
      candidates <- c(candidates, file.path(profile, ".cargo", "bin"))
    }
  }
  candidates <- unique(candidates[nzchar(candidates) & dir.exists(candidates)])
  if (length(candidates)) {
    Sys.setenv(PATH = paste(c(Sys.getenv("PATH"), candidates),
                            collapse = sep))
  }
  invisible(candidates)
}

mm_msrv_tool_version <- function(tool, install_msg) {
  if (!nzchar(Sys.which(tool))) {
    stop(paste(install_msg, collapse = "\n"), call. = FALSE)
  }
  out <- suppressWarnings(
    tryCatch(system2(tool, "--version", stdout = TRUE, stderr = TRUE),
             error = function(e) character())
  )
  status <- attr(out, "status") %||% 0L
  if (!length(out) || !identical(as.integer(status), 0L)) {
    stop(paste(install_msg, collapse = "\n"), call. = FALSE)
  }
  out[[1L]]
}

mm_msrv_check_minimum <- function(tool, detected_line, minimum) {
  if (is.na(minimum)) {
    return(invisible(TRUE))
  }
  detected <- mm_msrv_extract_semver(detected_line)
  if (is.na(detected)) {
    stop(sprintf(
      "Could not parse a version from `%s --version` output: %s",
      tool, detected_line
    ), call. = FALSE)
  }
  if (utils::compareVersion(minimum, detected) == 1L) {
    stop(sprintf(paste0(
      "\n------------------ [UNSUPPORTED %s VERSION] ------------------\n",
      "- Minimum supported %s version is %s.\n",
      "- Installed %s version is %s.\n",
      "Please update your Rust toolchain (https://rustup.rs) and retry.\n",
      "----------------------------------------------------------------"
    ), toupper(tool), tool, minimum, tool, detected), call. = FALSE)
  }
  invisible(TRUE)
}

mm_msrv_install_msg <- function(tool, minimum) {
  version_clause <- if (is.na(minimum)) "" else sprintf(" (>= %s)", minimum)
  c(
    sprintf("--------------------- [%s NOT FOUND] ---------------------",
            toupper(tool)),
    sprintf("The '%s' command%s was not found on the PATH or under", tool,
            version_clause),
    "~/.cargo/bin. Please install a Rust toolchain, e.g. from:",
    "  https://rustup.rs",
    "or from your OS package manager:",
    " - Debian/Ubuntu: apt-get install cargo rustc",
    " - Fedora/CentOS: dnf install cargo rustc",
    " - macOS: brew install rust",
    "and re-run the installation. This package never installs Rust itself.",
    "------------------------------------------------------------------"
  )
}

mm_msrv_main <- function() {
  desc <- read.dcf("DESCRIPTION")
  if (!"SystemRequirements" %in% colnames(desc)) {
    stop(paste(
      "`SystemRequirements` not found in `DESCRIPTION`.",
      "Please specify `SystemRequirements: Cargo (Rust's package manager), rustc`",
      sep = "\n"
    ), call. = FALSE)
  }
  sysreqs <- desc[, "SystemRequirements"]
  if (!grepl("cargo", sysreqs, ignore.case = TRUE) ||
      !grepl("rustc", sysreqs, ignore.case = TRUE)) {
    stop("`SystemRequirements` must declare both Cargo and rustc.",
         call. = FALSE)
  }

  cargo_min <- mm_msrv_declared_minimum(sysreqs, "cargo")
  rustc_min <- mm_msrv_declared_minimum(sysreqs, "rustc")

  mm_msrv_extend_path()

  cargo_line <- mm_msrv_tool_version("cargo",
                                     mm_msrv_install_msg("cargo", cargo_min))
  rustc_line <- mm_msrv_tool_version("rustc",
                                     mm_msrv_install_msg("rustc", rustc_min))

  mm_msrv_check_minimum("cargo", cargo_line, cargo_min)
  mm_msrv_check_minimum("rustc", rustc_line, rustc_min)

  message(sprintf("Using %s\nUsing %s", cargo_line, rustc_line))
  invisible(list(cargo = cargo_line, rustc = rustc_line))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!identical(tolower(Sys.getenv("MIXEFF_MSRV_NO_RUN")), "true")) {
  mm_msrv_main()
}
