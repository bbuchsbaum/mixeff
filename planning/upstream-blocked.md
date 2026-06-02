# Upstream-blocked features — revisit on every `mixeff-rs` re-pin

This is the durable registry of mixeff (R) work that is **blocked on
`mixeff-rs` (Rust) engine changes**. Each row links a mixeff mote bead to the
upstream bead that must land first, and says what to do here once it ships.

**Trigger to consult this file:** any time you bump `PINNED_REV` in
`tools/vendor-rust.R` (i.e. pull in new engine code). `tools/vendor-rust.R`
prints a pointer to this file at the end of a re-vendor run, and `AGENTS.md`
lists "consult `planning/upstream-blocked.md`" as a re-pin step. The mixeff
beads below are also tagged `upstream-blocked` (`mote ls --tag upstream-blocked`)
and held at `status=blocked` so they never surface as "ready" prematurely.

**On re-pin:** for each row, check whether the upstream bead is closed in the
peer store (`mote --store /Users/bbuchsbaum/code/rust/mixeff-rs/.mote show <id>`).
If it shipped: unblock the mixeff bead (`mote set <id> status=open`), do the
wrapper work, add a routine test, and remove the row here.

| mixeff bead | feature | upstream `mixeff-rs` bead | wrapper work once it lands |
|---|---|---|---|
| `bd-01KT2N1GGD61RHZ68KZFM9KD49` (T11) | negative-binomial GLMM family | `bd-01KT40XMMR7QC1TXKNX7DXZ6K5` | add `negative.binomial`/`MASS::negative.binomial(theta)` to the `glmm()` family allow-list (`R/glmm.R:204`), map to the engine NB family, pass `theta`, add a `glmer.nb` parity test |

### Landed and wired (kept for history)

| mixeff bead | feature | upstream bead | resolution |
|---|---|---|---|
| `bd-01KT2N1GKAQCZ1K8FYWHDC98KF` (T12) | caller optimizer / tolerance / `start=` warm-start control | `bd-01KT40AHAY9Z7SYM2B0FTVF53Z` ✅ shipped (`368a3fa`) | re-pinned to `368a3fa`; `mm_control(optimizer=, start=, ftol_rel=/ftol_abs=/xtol_rel=)` wired through the bridge; override recorded in `optimizer_certificate()`; PRD §8.3 reconciled. T12 closed. |

## Related (R-side already works; upstream is an enhancement, not a blocker)

These are NOT blocking anything — the R wrapper already delivers the capability;
the upstream bead would let us drop an R-side reconstruction for an
engine-certified path. Revisit opportunistically, not urgently.

| capability | upstream bead | why it's only an enhancement |
|---|---|---|
| GLMM `predict_new` (engine-certified new-data prediction) | `bd-01KT3XJD6SE54G7GT3Z2T6HMEX` | `predict.mm_glmm()` already computes population & conditional predictions R-side, validated against engine `fitted()`; an engine path would restore provenance + drop the R reconstruction |
