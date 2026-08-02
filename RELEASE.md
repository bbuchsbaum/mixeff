# Release process — mixeff

The one canonical acceptance command, run from a **clean checkout**:

``` sh
Rscript tools/release-gate.R
```

It builds ONE tarball and runs everything against that artifact:
`R CMD check --as-cran` on the tarball (never on the source directory),
`R CMD INSTALL <tarball>` into a throwaway gate library, then lint, the
fast test suite + error-UX battery + schema/manifest tests against that
installed build, and the offline-install gate. It writes
`release-gate-report.txt` — that file is the acceptance artifact to
attach to an R-universe / CRAN submission. Its header records the
evidence identity: mode, git commit (with a dirty-tree flag), the
mixeff-rs pin, tarball size + sha256, sha256 of `src/rust/vendor.tar.xz`
and the upstream `Cargo.lock`, and the R / rustc / cargo versions. The
gate also refreshes the “Size of tarball” figure in `cran-comments.md`
from the tarball it just built. Exit status is non-zero if any gate
fails.

Modes:

- default — the full gate above.
- `--fast` — skips the two slowest gates (as-cran check, offline
  install) during iteration; the report is stamped
  `FAST RUN -- NOT RELEASE EVIDENCE`.
- `--release-candidate` — the pre-submission run: the full gate, plus
  `MIXEFF_RUN_SLOW_PARITY=true` and `MIXEFF_RUN_APHANTASIA=true` in the
  suite, a live render of the aphantasia vignette
  (`MIXEFF_RUN_APHANTASIA_VIGNETTE=true`), and a skip audit — any suite
  skip outside a short allowlist (missing Suggests / pandoc, the deeper
  stress / joint-proof tiers) fails the gate.

## Fast suite vs. slow suites

The default test run is the **fast release suite** — a few minutes on an
installed build. Three opt-in tiers are gated behind environment
variables and are NOT part of the default gate (they can run for many
minutes); `--release-candidate` turns the first two on. The gate clears
all six `MIXEFF_*` tier flags before the suite, so exported shell
variables cannot change what a gate run tests:

| Tier | Enable with | What it covers |
|----|----|----|
| Slow parity | `MIXEFF_RUN_SLOW_PARITY=true` | large crossed / maximal random-slope LMM parity vs lme4 (Brown, sdamr, iamciera, OSF) |
| Aphantasia core | `MIXEFF_RUN_APHANTASIA=true` | Loo aphantasia reproduction |
| Aphantasia joint | `+ MIXEFF_APHANTASIA_JOINT=true` | the expensive joint-Laplace route |

Always run the suite against the **installed** package
(`load_package = "installed"`), not `devtools::test()` / `load_all()` —
the latter is a debug build ~60x slower and not representative of
release timing.

## Release checklist

Local gates (automated by `tools/release-gate.R`):

Vendor provenance current (`inst/LICENSE.note`, vendored-snapshot drift)

Tarball built; identity (commit, pin, sha256s, toolchain) recorded in
`release-gate-report.txt`; `cran-comments.md` tarball size current

`R CMD check --as-cran` on the built tarball: 0E / 0W / only the
new-submission NOTE

`lint_package()` == 0

Fast suite + error-UX battery + schema/manifest tests green (against the
installed tarball)

Offline source install (no network, `CARGO_NET_OFFLINE=true`)

Pre-submission run (`--release-candidate`, record the report):

Slow-parity suite green (`MIXEFF_RUN_SLOW_PARITY=true`) — any misses are
documented Decision-D near-ties, not regressions

Aphantasia reproduction green (`MIXEFF_RUN_APHANTASIA=true`) and the
vignette renders live

No suite skips outside the gate’s allowlist

CRAN extra-check court (run the manual **R-hub** workflow with its
default platform list before submission):

Primary R-devel compilers: `ubuntu-gcc16`, `ubuntu-clang`, `windows`,
`macos`, `macos-arm64`

Memory/undefined behavior: `clang-asan`, `clang-ubsan`, `gcc-asan`,
`m1-san`, `valgrind`

Dependency and numeric variants: `nosuggests`, `nold`, `lto`, `mkl`,
`atlas`

Expanded examples and native-code analysis: `donttest`, `rchk`

Every job is green, or every non-package failure is linked and
adjudicated in `cran-comments.md`; a green ordinary three-OS check is
not a substitute for this court

The default list intentionally mirrors the applicable public extra-check
classes maintained alongside Brian Ripley’s CRAN checks. Some historical
extra checks are not applicable here (for example C++-language-mode
checks: the package bridge is C plus Rust), and some are now defaults in
R-devel (`STRICT_R_HEADERS`). Keep the explicit list synchronized with
<https://r-hub.github.io/containers/> and
<https://www.stats.ox.ac.uk/pub/bdr/>.

External state — **each depends on something outside this checkout and
must be confirmed explicitly, not implied by local green tests**:

> **Self-contained repo:** the vendored engine (`src/rust/upstream/`,
> `vendor.tar.xz`, `vendor-config.toml`) is committed so R-universe
> (which builds from the repo, not a pre-vendored tarball) can build.
> Re-run `tools/vendor-rust.R` **and re-commit those paths** on every
> `PINNED_REV` bump, or the R-universe build breaks. CI never
> regenerates the snapshot before testing: every job drift-checks and
> then tests the committed tree, and a separate `vendor-drift` job
> (R-CMD-check.yaml) re-runs `tools/vendor-rust.R` in isolation and
> fails on any content difference (the lock file’s `vendored_utc`
> timestamp is compared field-wise without that line, and
> `vendor.tar.xz` by extracted content, since neither is
> byte-reproducible across runs/tar implementations).

`mixeff-rs` pin (`PINNED_REV` in `tools/vendor-rust.R`) is published on
the mixeff-rs GitHub default branch (CI vendoring needs it on the
remote)

`release/0.2.0-prep` merged to `main` and pushed; GitHub Actions (R CMD
check on 3 OS incl. the PDF-manual leg, vendor-drift, Release gates with
the full fast suite, pkgdown) green on the pushed head

R-hub v2 default CRAN extra-check court green

win-builder dry-run clean

mac-builder R-devel dry-run clean (closest public reproduction of the
Oxford arm64 M1 check; that machine has Cargo under `~/.cargo`, not on
its initial `PATH`)

R-universe: package registered in the bbuchsbaum registry; remote build
matrix green

CRAN: `cran-comments.md` current; submit; address reviewer findings

## Notes

- `NEWS.md` heads the current version’s section.
- `DESCRIPTION` `Version:` bumped; `Authors@R` present
  (Author/Maintainer are generated at build time — this is why the check
  must run on the tarball).
