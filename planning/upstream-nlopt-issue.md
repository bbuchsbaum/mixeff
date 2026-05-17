# Upstream issue draft — feature-gate `nlopt`

Filed against `/Users/bbuchsbaum/code/rust/mixedmodels`. The upstream repo
tracks work in `mote` (per its `AGENTS.md`); there is no `.github/` dir, so
the canonical form is a mote bead. A `mote new` invocation is provided at the
end; the body text is the issue.

---

## Title

`Feature-gate nlopt behind a Cargo feature for CRAN compatibility`

## Tags

`build, cargo, optimizer, downstream-cran`

## Priority

`p1` (downstream R wrapper is blocked from CRAN until this lands; not
blocking R-Universe distribution)

## Body

### Problem

`nlopt = "0.8.1"` is currently a non-optional dependency in `Cargo.toml`. The
`nlopt` crate wraps a C library and requires `CMake` at install time when no
system `libnlopt` is present. The CRAN build farm does not provide `CMake`
and rejects packages that need it, which means a downstream R wrapper
(`mixeff`, at `/Users/bbuchsbaum/code/mixeff`) cannot pass
`R CMD check --as-cran` while `mixedmodels` has `nlopt` as a hard dep.

This blocks the `mixeff` v0.2 CRAN submission. v0.1 ships to R-Universe and
is unaffected.

### Scope is narrow — abstraction already exists

The optimizer abstraction is already in place (`src/types/opt_summary.rs`),
with four variants:

```rust
enum Optimizer {
    Cobyla,           // pure-Rust, current default for small theta fallback
    PatternSearch,    // pure-Rust
    NloptNewuoa,      // large theta
    NloptBobyqa,      // small theta
}
```

`backend_name()` already groups these as `"nlopt"` vs `"native"`. The test
`test_large_theta_nlopt_matches_or_beats_cobyla_baseline` confirms the
pure-Rust path is competitive. **No new optimizer trait is needed.** The work
is to put the existing `nlopt` code paths behind `#[cfg(feature = "nlopt")]`
and route selection to the native optimizers when the feature is off.

### Required changes

#### 1. `Cargo.toml`

```toml
[features]
default = []
nlopt = ["dep:nlopt"]

[dependencies]
nlopt = { version = "0.8.1", optional = true }
# all other deps unchanged
```

#### 2. `src/model/linear.rs`

`#[cfg(feature = "nlopt")]` on:

- the `use nlopt::{...}` import (line 13)
- `fit_nlopt_small_theta_with_maxeval` (line 904 region, def at 2428)
- `fit_nlopt_large_theta_with_maxeval` (def at 2403)
- `fit_nlopt_with_algorithm` (def at 2458)
- `fit_nlopt_small_theta` (def at 2442)
- `fit_nlopt_large_theta` (def at 2630)
- `nlopt_ok` (def at 1801)
- `nlopt_status_label` (def at 1810)

For the two selectors:

- `use_nlopt_bobyqa_small_theta_optimizer` (line 1536) — `cfg`-gate the body
  to return `false` when the feature is off.
- `use_large_theta_nlopt_optimizer` (line 1556) — same.

The dispatch in `fit()` (around line 2658) then naturally routes to the
existing `cobyla` / `pattern_search` paths when the nlopt selectors return
`false`.

#### 3. `src/types/opt_summary.rs`

The cleanest approach: keep the `NloptNewuoa` and `NloptBobyqa` enum variants
unconditionally (so serialized opt summaries remain stable across builds)
and `cfg`-gate only the construction sites. Alternative: gate the variants
themselves and add a build-time-checked panic in deserializers — heavier and
serializer-fragile, recommend against.

Update:

- `backend_name()` (line 179) — unchanged.
- `optimizer_name()` (line 187) — unchanged; the strings are pure metadata.
- The `assert_eq!(opt.backend_name(), "nlopt")` test at line 430 and related
  assertions (lines 445, 466) → wrap with `#[cfg(feature = "nlopt")]`.

#### 4. Tests

- `test_large_theta_fit_uses_nlopt_newuoa` (linear.rs:8034) →
  `#[cfg(feature = "nlopt")]`
- `test_large_theta_nlopt_matches_or_beats_cobyla_baseline` (linear.rs:8052)
  → `#[cfg(feature = "nlopt")]`
- Any other test that constructs `Optimizer::Nlopt*` variants → same gate.

#### 5. CI

Add two build matrix entries:

| Job | Command | Purpose |
| --- | --- | --- |
| default | `cargo build && cargo test` | confirms the no-nlopt path is the default and works |
| nlopt | `cargo build --features nlopt && cargo test --features nlopt` | confirms parity path still works |

If CI currently passes a feature flag implicitly, drop it.

### Acceptance criteria

- `cargo build --no-default-features` succeeds without `CMake` available
  on the host, and the produced binary has no link-time dependency on
  `libnlopt`.
- `cargo test --no-default-features` runs the full test suite minus the
  `nlopt`-gated tests; all remaining tests pass.
- `cargo build --features nlopt` produces a binary behaviorally identical
  to today's default build.
- `cargo test --features nlopt` passes the full existing suite including
  `test_large_theta_nlopt_matches_or_beats_cobyla_baseline`.
- `cargo doc --no-default-features` succeeds.

### Numerical consequence (for downstream documentation only)

With default features off, fit dispatch routes to `cobyla` (small θ) and
`pattern_search` (large θ) instead of `BOBYQA` / `NEWUOA`. The existing
parity test confirms `cobyla` matches or beats the `nlopt` baseline on the
covered datasets, but minima are not bit-identical. The `mixeff` R wrapper
documents this in its non-goals (`mixeff/planning/PRD.md` §3, Decision D
resolved 2026-04-27) and uses loose tolerances on parity tests (`fixef`
1e-4, `theta` 1e-3, `logLik` 1e-3, `sigma` 1e-4 against `lme4` on
`sleepstudy`, `Pastes`, `cake`).

### Out of scope

- New optimizers beyond what's already implemented. (`argmin` 0.10 was
  mentioned in early downstream planning but the existing `cobyla` +
  `pattern_search` are sufficient — adding `argmin` is a separate decision.)
- Changing the default `nlopt` algorithm choice.
- Renaming `Optimizer` enum variants.

### Cross-references

- Downstream tracking bead:
  `bd-01KQ906S43Q5T2WD7GRZDAK7VZ` in `/Users/bbuchsbaum/code/mixeff/.mote`
- Downstream PRD:
  `/Users/bbuchsbaum/code/mixeff/planning/PRD.md` §5.2 item 1, §3 Decision D,
  §11 parity tolerances, §13 Decision A resolved, Risk R1

---

## Filing it

From `/Users/bbuchsbaum/code/rust/mixedmodels`:

```sh
mote new "Feature-gate nlopt behind a Cargo feature for CRAN compatibility" \
  -p 1 \
  --tag build --tag cargo --tag optimizer --tag downstream-cran \
  --body "$(cat /Users/bbuchsbaum/code/mixeff/planning/upstream-nlopt-issue.md | sed -n '/^## Body/,/^---$/p' | sed '1d;$d')"
```

Or paste the body manually:

```sh
mote new "Feature-gate nlopt behind a Cargo feature for CRAN compatibility" \
  -p 1 \
  --tag build --tag cargo --tag optimizer --tag downstream-cran
mote note <returned-id> --kind decision "<paste body>"
```

(`mote new --body` accepts the text inline; the `sed` extraction is just to
avoid hand-copying.)
