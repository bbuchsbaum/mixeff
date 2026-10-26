# cran-comments.md — mixeff 0.2.0 (draft; not ready to submit)

## Submission status

First CRAN submission, following an R-universe pre-release cycle. **Submission
is on hold:** an MKL reference-parity failure remains unexplained, Windows
builder results and merged-head R-universe qualification are outstanding.
PR 6 also requires an approving reviewer with write access before merge.

## Candidate identity

Source `806528106473594aa7342c6585c2b56be5ed2649`, engine pin
`004a5d9b12a5d71f483c6c5fc15dbf93ffdb8b7b`. Later release-document updates
are excluded from the package. The source tarball SHA-256 is
`a742c9973c1ed64a04d7d3aad3f1088cf53ff4c37cd14f1758c5539995ec0881`.

Size of tarball: 6807157 bytes

The vendored Rust sources allow offline installation. The expected local
incoming NOTE concerns a new submission and tarball size.

## Current executed checks (26 October 2026)

- macOS arm64, R 4.5.1: tarball `R CMD check --as-cran --no-manual`,
  0 errors, 0 warnings, 1 expected incoming NOTE.
- Canonical release-candidate gate: 14 passing steps, no failed/skipped
  steps; installed extended suite 6,881 passes, no failures/warnings,
  three documented allowlisted skips. Live aphantasia vignette and
  offline source installation pass. External xz and vendor-tar access
  were deliberately blocked during this run.
- GitHub Actions at 8065281: Ubuntu R-release/R-devel, macOS R-release,
  Windows R-release, PDF manual, release gates, vendor regeneration,
  pkgdown and coverage pass. R-code coverage: 94.2%.
- mac-builder job 1793045109-362c7021c79baa3a: Status OK including PDF
  manual, actual R 4.6.1 Patched. API flavor r-devel does not change
  the reported runtime version.
- win-builder: final R-release upload accepted; result pending. Final
  R-devel upload rejected with FTP550 while an older same-name upload
  remained queued. The bounded queue wait expired without retry.

## Current R-hub court

Run 40389108548 checked bf612d6; packaging-only 8065281 has identical
estimation runtime and tests. Ubuntu GCC/Clang, gcc16, Windows, both macOS
architectures, m1-san, donttest, lto and atlas pass. The four platforms
that exposed an external xz assumption (clang-asan, clang-ubsan, nold,
rchk) all pass repaired retry 40390625134. Valgrind remains pending.

gcc-asan and nosuggests failed before checking mixeff: lmerTest could not
load lme4 because libRlapack.so was missing. These are environmental gaps,
not passing package checks.

MKL job 121805238156 executed one failed Brown variance-component SD
comparison (mixeff about14.7, lme4 about8.2). A standalone MKL probe and
instrumented full-check run40398326331 pass unchanged assertions and
runtime; saved check fits agree at14.707846 versus14.719799. Independent
lme4 REML evaluation corroborates the latter fit. The original failure
saved no fitted objects and remains unexplained. No tolerance was relaxed.
This remains a submission blocker, not an environment-only exclusion.

## Remaining release work

Obtain the required PR review and merge; qualify the merged source through
R-universe; obtain final Windows builder receipts and the remaining court
results; resolve or explicitly adjudicate the MKL uncertainty. The source
currently on R-universe is older (6e9a62a), so its native successes do not
qualify this candidate. Its pre-existing wasm failure is tracked separately.

Exact receipts and the full readiness decision are in
`assessment/2026-10-26-cancellation-release.md` and
`release-evidence/cancellation-2026-10-26/`. Historical July/earlier October
results are retained in prior Git revisions and assessments, not presented
as current qualification here. No CRAN submission has occurred.

## CRAN policy notes

- **System requirements**: `Cargo` (Rust's package manager, >= 1.85.0),
  `rustc` (>= 1.85.0), `GNU make`. Documented in `SystemRequirements:`.
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
- **Installed size (historical measurements, refresh before submission)**: about 11 MB on Windows (libs 8.7 MB) and up to
  about 18 MB on Linux when the shared library is not stripped (libs
  15.4 MB; 12.2 MB after `strip --strip-debug`). Nearly all of it is the
  statically linked Rust engine, built with `lto = true` and
  `codegen-units = 1`; the R code is about 1 MB.
- **Cross-platform**: macOS arm64 + x86_64, Ubuntu LTS, Windows UCRT
  (Rtools43+ MinGW, `x86_64-pc-windows-gnu`). Windows i386 is not
  supported.
- **Optimizer**: the engine ships pure-Rust optimizers only (LMM fits
  route to `trust_bq` or `pattern_search` depending on the design;
  `cobyla` is also bundled). The `nlopt` C library is feature-gated off
  upstream and disabled in this build; CRAN machines do not need nlopt
  or CMake.
- **Serialization**: fitted model objects survive `saveRDS()` /
  `readRDS()`; all reported quantities are stored in the R object. Within
  a session a fit also holds an external pointer to the live engine model
  (`fit$rust_handle`) so follow-on computations (contrasts, prediction
  intervals, bootstraps) skip a refit. The pointer does not survive
  serialization; after reloading, queries restore the recorded native
  fitted-state snapshot without optimization. An explicitly different
  estimand (for example ML instead of REML) still requires a refit.
  `options(mixeff.keep_handle = FALSE)` disables the live cache.
  The bridge reattaches the R interrupt callback after restoration. Serial
  and threaded LMM/GLMM bootstrap cancellation is covered by R SIGINT tests.
- **Threads**: single-threaded by default. The parametric bootstrap
  uses more threads only when asked (`bootstrap_control(threads =)`);
  results are bit-identical for every thread count
  (including on Windows, where worker threads inherit R's x87 and MXCSR
  floating-point state).

## Reverse dependencies

None at this version; mixeff is a new package.
