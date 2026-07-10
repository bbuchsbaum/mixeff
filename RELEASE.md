# Release process — mixeff

The one canonical acceptance command, run from a **clean checkout**:

``` sh
Rscript tools/release-gate.R
```

It builds the tarball, runs `R CMD check --as-cran` on that tarball
(never on the source directory), lints, runs the fast test suite +
error-UX battery + schema/manifest tests against the installed release
build, and runs the offline-install gate. It writes
`release-gate-report.txt` — that file is the acceptance artifact to
attach to an R-universe / CRAN submission. Exit status is non-zero if
any gate fails. Use `--fast` to skip the two slowest gates (R CMD check,
offline install) during iteration.

## Fast suite vs. slow suites

The default test run is the **fast release suite** — a few minutes on an
installed build. Three opt-in tiers are gated behind environment
variables and are NOT part of the release gate (they can run for many
minutes):

| Tier | Enable with | What it covers |
|----|----|----|
| Slow parity | `MIXEFF_RUN_SLOW_PARITY=true` | large crossed / maximal random-slope LMM parity vs lme4 (Brown, sdamr, iamciera, OSF) |
| Aphantasia core | `MIXEFF_RUN_APHANTASIA=true` | Loo aphantasia reproduction |
| Aphantasia joint | `+ MIXEFF_APHANTASIA_JOINT=true` | the expensive joint-Laplace route |

Always run the suite against the **installed** package
(`load_package = "installed"`), not
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) /
`load_all()` — the latter is a debug build ~60x slower and not
representative of release timing.

## Release checklist

Local gates (automated by `tools/release-gate.R`):

Vendor provenance current (`inst/LICENSE.note`, vendored-snapshot drift)

`R CMD check --as-cran` on the built tarball: 0E / 0W / only the
new-submission NOTE

`lint_package()` == 0

Fast suite + error-UX battery + schema/manifest tests green

Offline source install (no network, `CARGO_NET_OFFLINE=true`)

Opt-in verification (run separately, record outcomes):

Slow-parity suite reviewed (`MIXEFF_RUN_SLOW_PARITY=true`) — any misses
are documented Decision-D near-ties, not regressions

Aphantasia reproduction reviewed if the fixture changed

External state — **each depends on something outside this checkout and
must be confirmed explicitly, not implied by local green tests**:

> **Self-contained repo:** the vendored engine (`src/rust/upstream/`,
> `vendor.tar.xz`, `vendor-config.toml`) is committed so R-universe
> (which builds from the repo, not a pre-vendored tarball) can build.
> Re-run `tools/vendor-rust.R` **and re-commit those paths** on every
> `PINNED_REV` bump, or the R-universe build breaks.

`mixeff-rs` pin (`PINNED_REV` in `tools/vendor-rust.R`) is published on
the mixeff-rs GitHub default branch (CI vendoring needs it on the
remote)

`release/0.2.0-prep` merged to `main` and pushed; GitHub Actions (R CMD
check on 3 OS, Release gates, pkgdown) green on the pushed head

R-hub v2 (linux / windows / macos, R-devel) green

win-builder + mac-builder dry-runs clean

R-universe: package registered in the bbuchsbaum registry; remote build
matrix green

CRAN: `cran-comments.md` current; submit; address reviewer findings

## Notes

- `NEWS.md` heads the current version’s section.
- `DESCRIPTION` `Version:` bumped; `Authors@R` present
  (Author/Maintainer are generated at build time — this is why the check
  must run on the tarball).
