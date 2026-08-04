# cran-comments.md — mixeff 0.2.0

## Submission target

First CRAN submission of mixeff (0.2.0), following an R-Universe
pre-release cycle. The upstream Rust engine (`mixeff-rs`) is bundled as
a pinned, vendored snapshot (`tools/vendor-rust.R`; provenance in
`src/rust/upstream/mixeff-rs.lock`), and its optional `nlopt` C-library
feature is disabled for this build, so the source package builds with
the declared Rust toolchain requirements only.

## Test environments

Checked as the built source tarball throughout (never the source
directory); the acceptance run is `tools/release-gate.R
--release-candidate`, which builds one tarball, records its SHA-256,
checks that artifact with `--as-cran`, installs that same artifact, and
runs the full test suite plus the slow parity and real-data
reproduction tiers against it.

- macOS arm64, R 4.5.1 — `--as-cran` on the built tarball with
  vignettes rebuilt: 0 errors, 0 warnings, 1 NOTE (below).
- macOS arm64, R 4.5.1 — the same tarball with
  `_R_CHECK_DEPENDS_ONLY_=true`: Status OK.
- mac-builder, macOS arm64, R 4.6.1 — **Status OK** (0 errors,
  0 warnings, 0 notes), vignettes rebuilt, PDF manual built.
- win-builder, R-devel (2026-07-30 r90327 ucrt) and R-release — 0
  errors, 0 warnings, 2 NOTEs (both below).
- GitHub Actions — Ubuntu R-devel and R-release, macOS, Windows UCRT:
  all green, including the PDF manual leg.
- R-hub v2 extra-check court, 18 platforms
  (`.github/workflows/rhub.yaml`): primary R-devel compilers
  (`ubuntu-gcc12`, `ubuntu-clang`, `gcc16`), Windows and macOS R-devel,
  sanitizers (`clang-ubsan`, `m1-san`, `clang-asan`, `gcc-asan`),
  `valgrind`, `rchk`, `nosuggests`, `nold`, `lto`, `mkl`, `atlas`,
  `donttest`. Results and the four non-green platforms are itemized
  below.

## R CMD check results

1. **CRAN incoming feasibility** (all platforms) — "New submission", plus
   the tarball size on some machines. The Rust sources are vendored so
   the package builds fully offline, which accounts for the size.

2. **`checking compiled code`, win-builder only** — this check does not
   run there:

   ```
   Error in ccE(lines, flags = new_flags, include = include) :
     'cc' is not on the path
   ```

   `tools:::ccE` shells out to `cc -E` to preprocess R's own headers
   while assembling the list of public R API entry points; win-builder's
   Rtools provides `gcc` but no `cc`, so the check aborts before any
   object of ours is examined and reports no finding about the package.
   The same check returns **OK** wherever `cc` is present: macOS
   (local and mac-builder), Linux R-devel (R-hub `clang-asan` and
   `nosuggests` containers), and the Windows UCRT CI leg. The adjacent
   `checking Rust compilation` is OK on win-builder.

## R-hub court results

The court was run twice: once on the pre-fix tree, which surfaced four
non-green platforms, and once after the fixes below, which reproduces
**16 of 18 platforms green**. Green includes R-devel on three compilers
(`ubuntu-gcc12`, `ubuntu-clang`, `gcc16`), Windows and both macOS
R-devel builds, the sanitizers `clang-asan`, `clang-ubsan` and `m1-san`,
plus `rchk`, `nosuggests`, `lto`, `mkl`, `atlas`, `nold` and `donttest`.

Three of the original four were real and are fixed:

- **`rchk`** and the **win-builder installation ERROR** shared one root
  cause: an install-time `cargo run --bin document` step that
  regenerated the extendr wrappers. It failed two different ways — on
  Windows `cargo run` with a `--target` must execute a cross-compiled
  binary ("%1 is not a valid Win32 application", os error 193), and in
  the rchk container linking an *executable* against a **static**
  `libR.a` needs GNU readline symbols that are not on the link line
  (every undefined symbol referenced by R's own `sys-std.c`, none by
  mixeff). Installing a package must not rewrite its own R sources in
  any case; the wrappers are committed and ship in the tarball, and the
  step is gone. Both platforms now pass, and a guard test keeps the
  committed wrappers in step with the `#[extendr]` surface.

- **`clang-asan`** — the R process was aborted by the V8 JavaScript
  engine, reached through the `jsonvalidate` Suggests dependency used by
  one schema-validation test; no mixeff frame appeared in the trace.
  Because a fatal abort inside a dependency cannot be trapped, that test
  is now `skip_on_cran()`; the schema contract is still verified in the
  release gate, which runs the suite with `NOT_CRAN=true`. `clang-asan`
  now passes.

- **`nosuggests`** — vignette re-building needed `rmarkdown`, which that
  container omits. `rmarkdown` genuinely is required to build these
  vignettes (`rmarkdown::html_vignette`), so it is declared in
  `VignetteBuilder:` alongside `knitr`; both are in `Suggests`. It now
  passes, and the same condition reproduced locally with
  `_R_CHECK_DEPENDS_ONLY_=true` gives Status OK.

The two platforms that are still not green never executed a check
against the package:

- **`valgrind`** — failed in `r-hub/actions/setup-deps` on both runs, at
  the same step, before the check started. No check log is produced.

- **`gcc-asan`** — hit the GitHub Actions six-hour job limit on both
  runs (22:05:14 to 04:05:20 UTC, `run-check` still in progress when the
  runner terminated it). AddressSanitizer plus a from-source Rust build
  and the full test suite does not fit the runner budget. Sanitizer
  coverage is provided by `clang-asan` and `m1-san`, which pass.

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
