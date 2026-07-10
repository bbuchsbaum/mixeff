# cran-comments.md — mixeff 0.2.0

## Submission target

First CRAN submission of mixeff (0.2.0), following an R-Universe
pre-release cycle. The upstream Rust engine is bundled as a pinned,
vendored snapshot (`tools/vendor-rust.R`; provenance in
`src/rust/upstream/mixeff-rs.lock`) and its optional `nlopt` C-library
path is feature-gated off for this build, so the source package builds
with the declared Rust toolchain requirements only.

## Test environments

- macOS arm64, R 4.5.x — `R CMD check --as-cran --no-manual` with
  vignette rebuilding: results recorded below (this checkout).
- GitHub Actions (ubuntu, macOS, windows UCRT) — R CMD check green on
  the pre-release branch head; R-hub v2 (linux/windows/macos R-devel)
  green at the previous pin.
- win-builder / mac-builder dry-runs — to be run at submission time.

## R CMD check results

`R CMD check --as-cran --no-manual` on mixeff 0.2.0 (macOS arm64,
vignettes rebuilt): 0 errors, 0 warnings, 1 NOTE.

1. **New submission / tarball size** — informational for a first
   submission; the tarball is about 6.4 MB (Rust sources vendored for
   fully offline builds).

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
