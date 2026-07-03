# mixeff: Audit-First Mixed-Effects Models via the 'mixedmodels' Rust Crate

An R wrapper for the `mixedmodels` Rust crate. The package is
audit-first: every printed claim traces back to a versioned JSON
artifact produced by the Rust compiler, and the package refuses to
fabricate inference results the engine cannot certify. See
[`vignette("intro", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/intro.md)
for an overview and the demystification surface for random-effects
syntax (Phase 1+).

## Author

**Maintainer**: Brad Buchsbaum <bbuchsbaum@research.baycrest.org>
