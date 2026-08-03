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

## R-hub court: the four non-green platforms

None of these is a defect in the package; each is itemized with its
evidence.

- **`valgrind`** — failed in `r-hub/actions/setup-deps`, before the
  check ran. No check log was produced; nothing was executed against
  the package.

- **`clang-asan`** — the R process was aborted by the V8 JavaScript
  engine, reached through the `jsonvalidate` Suggests dependency used
  by one schema-validation test:

  ```
  # Fatal error in , line 0
  # Check failed: static_cast<uintptr_t>(caller_frame_top_) ... real_jslimit()
  ...V8.so(v8::internal::Deoptimizer::DoComputeOutputFrames()+0x620)
  ```

  No mixeff frame appears in the trace. Because a fatal abort inside a
  dependency cannot be trapped, that test is now `skip_on_cran()`; the
  schema contract is still verified in the release gate, which runs the
  suite with `NOT_CRAN=true`. `clang-ubsan` and `m1-san` pass.

- **`nosuggests`** — vignette re-building needed `rmarkdown`, which that
  container omits. `rmarkdown` is genuinely required to build these
  vignettes (`rmarkdown::html_vignette`), so it is now declared in
  `VignetteBuilder:` alongside `knitr` (both are in `Suggests`). The
  same condition reproduced locally with `_R_CHECK_DEPENDS_ONLY_=true`
  gives Status OK.

- **`rchk`** — the same root cause as the win-builder installation
  failure, surfacing differently. The library built cleanly
  (`Finished release profile ... in 1m 15s`); the failure was the
  install-time `cargo run --bin document` step that followed, which
  links an *executable* against that container's **static** `libR.a`
  and so needs GNU readline symbols that are not on the link line:

  ```
  error: linking with `cc` failed: exit status: 1
  rust-lld: error: undefined symbol: tilde_expand_word
    >>> referenced by sys-std.c:505 ... in archive libR.a
  rust-lld: error: undefined symbol: rl_readline_name, add_history, ...
  ```

  Every undefined symbol belongs to readline and is referenced by R's own
  `sys-std.c`, not by mixeff. That step has been removed: installation no
  longer regenerates the extendr wrappers (they are committed and ship in
  the tarball), which also fixed the win-builder ERROR. rchk's actual
  PROTECT analysis never ran because the build stopped here.


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
