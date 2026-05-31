# cran-comments.md — mixeff 0.1.0 (R-Universe pre-release)

## Submission target

This document tracks notes for the next CRAN submission. mixeff 0.1.0
is still in R-Universe pre-release hardening; no CRAN submission has
been made from this checkout. The upstream Rust engine is bundled as a
pinned, vendored snapshot and its optional `nlopt` C-library path is
feature-gated off for this build, so the source package builds with the
declared Rust toolchain requirements only.

## Test environments

- macOS arm64, R 4.5.1 — `R CMD build --no-manual` with vignette
  rebuilding enabled: clean.
- macOS arm64, R 4.5.1 — `R CMD check --as-cran --no-manual` with
  `CARGO_NET_OFFLINE=true`: 0 errors, 0 warnings, 2 NOTEs.
- R-Universe API (`bbuchsbaum.r-universe.dev/api/packages/mixeff`) —
  404 at dry-run time; no current job matrix is available yet.
- Windows UCRT / mac-builder remote dry-runs — pending external
  submission.

## R CMD check results

There were 0 errors, 0 warnings, and 2 NOTEs:

1. **New submission / tarball size** — informational for a first
   submission; the dry-run tarball was about 5.8 MB.
2. **Unable to verify current time** — environmental; would not fire on
   the CRAN clock-synced runner.

## Downstream dependencies

None.

## CRAN policy notes

- **System requirements**: `Cargo` (Rust's package manager, >= 1.78.0),
  `rustc` (>= 1.78.0), `GNU make`. Documented in `SystemRequirements:`.
- **Vendoring**: The upstream `mixedmodels` Rust crate and its
  transitive Cargo registry dependencies are vendored under
  `src/rust/upstream/` and `src/vendor/` (via `tools/vendor-rust.R`).
  The full vendored set is reconstituted from `src/rust/vendor.tar.xz`
  at `R CMD INSTALL` time so `R CMD build`'s clean step does not break
  the offline build.
- **Tarball size**: ~5.8 MB source tarball; installed size is about
  8.1 MB, dominated by the Rust-compiled library. `.Rbuildignore`
  aggressively trims tests, datasets, and examples from vendored
  crates; LICENSE / NOTICE files are preserved in `inst/LICENSE.note`.
- **Cross-platform**: macOS arm64 + x86_64, Ubuntu LTS, Windows UCRT
  (Rtools43+ MinGW, `x86_64-pc-windows-gnu`). Windows i386 dropped.
- **Optimizer**: ships with the upstream's pure-Rust optimizer (cobyla
  / pattern_search). The `nlopt` C library (which lme4 uses for BOBYQA)
  is feature-gated upstream and disabled in this build; CRAN does not
  need to install nlopt or CMake.
- **Reproducibility**: every printed claim traces to a versioned JSON
  artifact emitted by the Rust compiler. The R object survives
  `saveRDS()` / `readRDS()` without a live Rust handle.

## Reverse dependencies

None at this version; mixeff is a new package.
