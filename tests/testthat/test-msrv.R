# Configure-time Rust toolchain preflight (tools/msrv.R, Audit1 WI-6.3).
# Source-tree-only: the script ships in tools/, not in the installed package.

mm_msrv_env <- function() {
  path <- test_path("../../tools/msrv.R")
  skip_if_not(file.exists(path), "tools/msrv.R is not available")
  env <- new.env(parent = baseenv())
  old <- Sys.getenv("MIXEFF_MSRV_NO_RUN", unset = NA_character_)
  Sys.setenv(MIXEFF_MSRV_NO_RUN = "true")
  on.exit({
    if (is.na(old)) Sys.unsetenv("MIXEFF_MSRV_NO_RUN")
    else Sys.setenv(MIXEFF_MSRV_NO_RUN = old)
  })
  sys.source(path, envir = env)
  env
}

test_that("declared minima parse despite commas inside parentheses", {
  env <- mm_msrv_env()
  sysreqs <- "Cargo (Rust's package manager, >= 1.78.0), rustc (>= 1.78.0), GNU make"
  expect_identical(env$mm_msrv_declared_minimum(sysreqs, "cargo"), "1.78.0")
  expect_identical(env$mm_msrv_declared_minimum(sysreqs, "rustc"), "1.78.0")
  # No declared minimum -> NA, and the check is then a no-op.
  expect_true(is.na(env$mm_msrv_declared_minimum("Cargo, rustc", "cargo")))

  # The shipped DESCRIPTION must itself parse to BOTH minima -- the old
  # ", " split fragmented the cargo clause so its minimum was never checked.
  desc <- read.dcf(test_path("../../DESCRIPTION"))
  shipped <- desc[, "SystemRequirements"]
  expect_false(is.na(env$mm_msrv_declared_minimum(shipped, "cargo")))
  expect_false(is.na(env$mm_msrv_declared_minimum(shipped, "rustc")))
})

test_that("version comparison enforces minima with clear failures", {
  env <- mm_msrv_env()
  expect_true(env$mm_msrv_check_minimum("rustc", "rustc 1.82.0 (f6e511eec)",
                                        "1.78.0"))
  expect_error(
    env$mm_msrv_check_minimum("cargo", "cargo 1.70.0 (ec8a8a0ca)", "1.78.0"),
    "Minimum supported cargo version"
  )
  expect_error(
    env$mm_msrv_check_minimum("cargo", "garbage output", "1.78.0"),
    "Could not parse"
  )
})

test_that("PATH extension appends cargo bin dirs as clean entries", {
  env <- mm_msrv_env()
  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path))
  appended <- env$mm_msrv_extend_path()
  parts <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1L]]
  expect_true(all(nzchar(parts)))
  # Every existing candidate dir must land as its OWN entry in the split
  # PATH — the literal-":" bug glued the dir onto the previous entry, so it
  # would never appear as an exact element after a platform-separator split.
  for (dir in appended) {
    expect_true(dir %in% parts,
                info = sprintf("%s not a standalone PATH entry", dir))
  }
})

test_that("missing tools fail with the plain installation instruction", {
  env <- mm_msrv_env()
  msg <- env$mm_msrv_install_msg("cargo", "1.78.0")
  expect_true(any(grepl("rustup.rs", msg, fixed = TRUE)))
  expect_true(any(grepl("never installs Rust", msg, fixed = TRUE)))
  expect_error(
    env$mm_msrv_tool_version("definitely-not-a-real-tool-xyz",
                             env$mm_msrv_install_msg("cargo", "1.78.0")),
    "NOT FOUND"
  )
})
