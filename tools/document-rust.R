#!/usr/bin/env Rscript

# Regenerate R/extendr-wrappers.R from the #[extendr] surface.
#
# This used to run inside src/Makevars during installation. It must not:
# installing a package should never rewrite its own R sources, the source
# tree may be read-only, and `cargo run` has to EXECUTE the binary it just
# built -- which fails when the host and --target triples differ (observed
# on win-builder R-devel: "%1 is not a valid Win32 application", os error
# 193, after the library itself had built cleanly).
#
# So wrapper generation is a DEVELOPMENT step. Run it after changing the
# #[extendr] surface, then commit the regenerated wrappers:
#
#   Rscript tools/document-rust.R
#
# tests/testthat/test-extendr-wrappers.R guards against forgetting.
#
# Note the deliberate absence of --target: the document binary runs on this
# machine, so it must be built for this machine.

pkg_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(pkg_root, "DESCRIPTION")))

manifest <- file.path(pkg_root, "src", "rust", "Cargo.toml")
if (!file.exists(manifest)) {
  stop("src/rust/Cargo.toml not found; run tools/vendor-rust.R first.")
}

if (!nzchar(Sys.which("cargo"))) {
  stop("cargo not found on PATH; install a Rust toolchain (https://rustup.rs).")
}

target_dir <- file.path(pkg_root, "src", "rust", "target")
args <- c("run", "--bin", "document",
          "--manifest-path", shQuote(manifest),
          "--target-dir", shQuote(target_dir))

message("[document-rust.R] regenerating R/extendr-wrappers.R ...")
status <- system2("cargo", args)
if (!identical(status, 0L)) {
  stop("[document-rust.R] cargo run --bin document failed (status ", status, ")")
}

wrappers <- file.path(pkg_root, "R", "extendr-wrappers.R")
if (!file.exists(wrappers)) {
  stop("[document-rust.R] the document binary did not write ", wrappers)
}
message(sprintf("[document-rust.R] done: %s (%d lines)",
                wrappers, length(readLines(wrappers, warn = FALSE))))
message("[document-rust.R] commit the result if it changed.")
