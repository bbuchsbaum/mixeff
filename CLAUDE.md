# CLAUDE.md — `mixeff`

Project-specific guidance for Claude Code working in this repository.

## What this is

`mixeff` is an R package wrapping the `mixeff-rs` Rust crate
(GitHub: `bbuchsbaum/mixeff-rs`), with a local peer checkout at
`/Users/bbuchsbaum/code/rust/mixeff-rs`. The crate is bundled into the
package as a **pinned snapshot** by `tools/vendor-rust.R`; the pinned commit
SHA is the `PINNED_REV` constant in that script (the committed source of
truth), and `src/rust/upstream/mixeff-rs.lock` records the resolved
provenance. The PRD, vision, and mission live in `planning/`. Read those
before making non-trivial changes:

- `planning/vision.md`
- `planning/mission.md`
- `planning/PRD.md` — the canonical specification

The Rust compiler and inference contracts are at
`/Users/bbuchsbaum/code/rust/mixeff-rs/docs/`.

## Issue tracking: mote

This project uses [`mote`](https://github.com/bbuchsbaum/mote) for local issue
tracking, path reservations, and inter-agent coordination. **Always** start a
working session with:

```bash
mote doctor
mote actor show
mote board
mote ready
```

Reserve paths before editing:

```bash
mote preflight --issue <bd-id> --paths <path> ...
mote begin    <bd-id> --paths <path> ... --note "starting"
```

Record progress, decisions, and blockers as `mote note` entries on the active
bead. Close with `mote done` or hand off with `mote handoff`. Full protocol
(commands, exit codes, repository policy) is in `AGENTS.md`.

`.mote/` is local state and is not committed to git.

### Upstream `mixeff-rs` bugs

If work in `mixeff` uncovers a concrete upstream Rust-engine bug, you may create
a mote issue in `/Users/bbuchsbaum/code/rust/mixeff-rs`. Use that path only
for actionable bugs with evidence: a minimal reproducer or fixture case,
reference behavior, observed Rust behavior, commits, and tolerance/status
context.

For feature requests, design/API requests, prioritization questions, or
speculative improvements, discuss with the user before creating an upstream
mote issue.

## Working principles for this repo

- **Phases come from `planning/PRD.md` §10.** Do not invent new phase scope
  without updating the PRD first.
- **Decisions A–D (`planning/PRD.md` §13) are unresolved.** Surface them, do
  not assume answers.
- **The upstream Rust crate is a peer, not a submodule.** Treat
  `/Users/bbuchsbaum/code/rust/mixeff-rs` as a separate repository; if you
  need to change it, reserve there explicitly and coordinate via `mote msg`.
  Bumping the bundled snapshot means updating `PINNED_REV` in
  `tools/vendor-rust.R` and re-running it — never hand-edit
  `src/rust/upstream/`.
- **No silent surgery.** The package's contract with users is that every model
  reduction or refusal crosses the boundary as a structured diagnostic. Any
  code path that hides a transformation from the user is a bug, not a
  convenience.
- **JSON artifacts are the source of truth.** External pointers are caches.
  Code that assumes a live Rust handle without the handle-revival path is a
  bug.

## When in doubt

- Architecture or scope question → `planning/PRD.md`.
- Why we're doing this at all → `planning/vision.md`, `planning/mission.md`.
- mote workflow → `AGENTS.md`.
- The upstream Rust contract → `/Users/bbuchsbaum/code/rust/mixeff-rs/docs/`.
