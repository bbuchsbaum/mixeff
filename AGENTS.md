# AGENTS.md — `mixeff`

This file is the protocol for any coding agent (Claude Code, Codex, Cursor, etc.)
working in this repository.

## Project orientation

`mixeff` is an R wrapper for the `mixedmodels` Rust crate at
`/Users/bbuchsbaum/code/rust/mixedmodels`. Background reading, in order:

- `planning/vision.md` — the long-arc world the package is building toward.
- `planning/mission.md` — what we do, who we serve, five operating principles.
- `planning/PRD.md` — full product requirements, API surface, phasing, risks.
- `planning/r_layer_proposal.md` — the upstream R-layer proposal (source of truth
  for design intent).
- `/Users/bbuchsbaum/code/rust/mixedmodels/docs/` — the Rust compiler and
  inference contracts the wrapper speaks to.

Phases 0–5 in `planning/PRD.md` §10 are the canonical work breakdown.

## Issue tracking and planning: mote

This project uses [`mote`](https://github.com/bbuchsbaum/mote) as a local,
daemonless issue and coordination tracker. Treat `.mote/ops/` as the append-only
source of truth. Reservations are advisory; agents must voluntarily check and
respect them.

### One-time setup

```bash
cargo install --git https://github.com/bbuchsbaum/mote --locked
mote --version
mote init
mote actor set <actor>      # e.g. claude-fixer, claude-docs, bbuchsbaum
mote doctor
```

### Starting a session

Run from the repo root or a subdirectory:

```bash
mote doctor
mote actor show
mote board
mote ready
```

Use a stable actor name. Do not invent a new actor each turn unless the work is
intentionally separate.

### Starting work

Identify the bead and the likely path scope. If no suitable bead exists, create
one:

```bash
mote new "Short task title" -p <0|1|2|3> --tag <area>
```

Reserve before editing:

```bash
mote preflight --issue <bd-id> --paths <path> [<path> ...]
mote begin    <bd-id> --paths <path> [<path> ...] --note "starting"
```

Keep reservations narrow. Prefer exact files; reserve directories only for broad
changes. If `preflight` or `begin` exits `2`, inspect the conflict before
overriding:

```bash
mote who-has <path>
```

For long or broad work touching shared git state, prefer `git worktree` over
relying solely on mote reservations.

### During work

Record material progress, blockers, and decisions:

```bash
mote note <bd-id> --kind progress "what changed"
mote note <bd-id> --kind decision "decision and why"
mote note <bd-id> --kind blocker "what is blocked"
```

Coordinate with other actors via messages, not by guessing intent:

```bash
mote inbox
mote msg send --to <actor> --issue <bd-id> --kind request "short request"
mote msg ack <msg-id>
```

Before expanding scope, re-run `preflight` and reserve the added paths.

### Upstream `mixedmodels` issues

When mixeff work exposes a concrete upstream Rust-engine bug, agents may lodge
an issue in the peer repository's mote store at
`/Users/bbuchsbaum/code/rust/mixedmodels`. Reserve this for bugs with actionable
evidence: a minimal reproducer or fixture case, expected/reference behavior,
observed Rust behavior, relevant commits, and any tolerance/status context.

Do **not** create upstream mote issues for feature requests, design/API
requests, prioritization questions, or speculative improvements without first
discussing them with the user. Those should usually remain mixeff planning notes
until the user agrees they should become upstream work.

### Finishing or handoff

```bash
mote done    <bd-id> --note "finished"
mote handoff <bd-id> --to <actor> --note "state and next step" --release
```

If stopping without completing or handing off, release leases explicitly:

```bash
mote release   <bd-id>
mote unreserve <rv-id>
```

### Exit codes

- `0` success
- `2` reducer rejected (stale scalar clock, path overlap, duplicate ack) — do
  *not* retry blindly; inspect with `mote show <bd-id>`,
  `mote history <bd-id> --include-rejected`, `mote who-has <path>`.
- `3` invalid command — check actor identity and arguments.
- `4` repository or storage problem — `mote doctor`, `mote fsck`
  (use `mote fsck --clean-tmp` only when stale tmp files should be cleaned).

### Repository policy

`.mote/` is **not** committed to git by default. The op log is local
coordination state, not a deliverable. Never hand-edit `.mote/ops/*.json`;
publish changes through the CLI.

## Editing conventions

- The R package source lives under `R/`; Rust glue under `src/rust/src/`. Both
  are detailed in `planning/PRD.md` §6.
- The upstream Rust crate at `/Users/bbuchsbaum/code/rust/mixedmodels` is a
  *peer* repository, not a submodule. Do not edit it from this repo unless you
  reserve paths there explicitly and coordinate via mote message.
- Decisions A–D in `planning/PRD.md` §13 must be resolved before Phase 0
  begins. Surface them, do not assume them.
