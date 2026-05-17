# cran-comments.md — mixeff 0.1.0 (R-Universe pre-release)

## Submission target

This document tracks notes for the next CRAN submission. mixeff 0.1.0
is distributed through R-Universe (bbuchsbaum.r-universe.dev) only;
CRAN submission is held until the upstream `mixedmodels` Rust crate's
`nlopt` feature-gate PR lands so the package can build on CRAN's
toolchain without requiring a C `nlopt` library at build time (CRAN
policy on system dependencies).

## Test environments

- macOS arm64, R 4.5.1 — `R CMD check --as-cran` clean (3 NOTEs, all
  environmental: new submission, future-file-timestamp, HTML Tidy
  version).
- Linux (R-Universe build farm) — TBD on first publish.
- Windows UCRT (R-Universe build farm) — TBD on first publish.

## R CMD check results

There were 0 errors, 0 warnings, and 3 NOTEs:

1. **New submission** — informational; appears once.
2. **Unable to verify current time** — environmental; would not fire on
   the CRAN clock-synced runner.
3. **HTML Tidy version** — local toolchain note; would not fire on
   the CRAN runner.

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
- **Tarball size**: ~5.8 MB installed, dominated by the Rust-compiled
  library. `.Rbuildignore` aggressively trims tests, datasets, and
  examples from vendored crates; LICENSE / NOTICE files are preserved
  in `inst/LICENSE.note`.
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
