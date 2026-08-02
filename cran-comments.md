# cran-comments.md — mixeff 0.2.0

## Submission target

First CRAN submission of mixeff (0.2.0), following an R-Universe
pre-release cycle. The upstream Rust engine (`mixeff-rs`) is bundled as
a pinned, vendored snapshot (`tools/vendor-rust.R`; provenance in
`src/rust/upstream/mixeff-rs.lock`), and its optional `nlopt` C-library
feature is disabled for this build, so the source package builds with
the declared Rust toolchain requirements only.

## Test environments

- macOS arm64, R 4.5.1 — `R CMD check --as-cran` on the built tarball,
  with vignette rebuilding: results recorded below (current release-prep
  checkout).
- macOS arm64, R 4.5.1 — the same built tarball with
  `_R_CHECK_DEPENDS_ONLY_=true` and Suggested packages not forced:
  results recorded below.
- GitHub Actions (ubuntu R-devel + release, macOS, windows UCRT) —
  R CMD check green on the release-prep head.
- R-hub v2 extra-check court (sanitizers, valgrind, noSuggests, MKL/ATLAS,
  LTO, rchk; the platform list in `.github/workflows/rhub.yaml`) and
  win-builder / mac-builder — to be run at submission time; results will
  be recorded here.

## R CMD check results

`R CMD check --as-cran` on the built mixeff 0.2.0 tarball (macOS arm64,
vignettes rebuilt): 0 errors, 0 warnings, 1 NOTE. The stricter
depends-only check produced the same 0 errors, 0 warnings, 1 NOTE result.

1. **CRAN incoming feasibility** — the single NOTE carries exactly two
   lines beyond the maintainer address: "New submission", and
   "Size of tarball: 6240378 bytes" (the Rust sources are vendored for
   fully offline builds).

## Downstream dependencies

None.

## CRAN policy notes

- **System requirements**: `Cargo` (Rust's package manager, >= 1.78.0),
  `rustc` (>= 1.78.0), `GNU make`. Documented in `SystemRequirements:`.
- **Vendoring**: The upstream `mixeff-rs` Rust crate is vendored under
  `src/rust/upstream/`, and its transitive Cargo registry dependencies
  ship as `src/rust/vendor.tar.xz` (about 4 MB, the largest single item
  in the tarball). `src/Makevars` unpacks that archive to `src/vendor/`
  at `R CMD INSTALL` time — reconstituting it at install rather than
  shipping it unpacked keeps `R CMD build`'s cleanup from breaking the
  offline build. Compilation runs `cargo build --offline` against the
  vendored registry (`vendor-config.toml`); no network access is needed
  or attempted.
- **Vendored-crate trimming**: `tools/vendor-rust.R` strips tests,
  benchmarks, examples, datasets, and CI configuration from the
  vendored crates before they are committed; upstream LICENSE / NOTICE
  files are preserved in `inst/LICENSE.note`.
- **Installed size**: about 8.6 MB, dominated by the Rust-compiled
  library.
- **Cross-platform**: macOS arm64 + x86_64, Ubuntu LTS, Windows UCRT
  (Rtools43+ MinGW, `x86_64-pc-windows-gnu`). Windows i386 is not
  supported.
- **Optimizer**: the engine ships pure-Rust optimizers only (LMM fits
  route to `trust_bq` or `pattern_search` depending on the design;
  `cobyla` is also bundled). The `nlopt` C library is feature-gated off
  upstream and disabled in this build; CRAN machines do not need nlopt
  or CMake.
- **Serialization**: fitted model objects survive `saveRDS()` /
  `readRDS()`; all reported quantities are stored in the R object and
  no live external handle is required after reloading.

## Reverse dependencies

None at this version; mixeff is a new package.
