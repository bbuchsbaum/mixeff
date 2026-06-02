# Plan — closing the lme4 *functional-equivalence* gap

**Status:** proposal (2026-06-01) · **Branch context:** `fix/lme4-parity-wrapper-beads`
**Goal:** make `mixeff` *functionally equivalent, or close,* to `lme4` — i.e. anything an
lme4 user can compute, `mixeff` can compute with equivalent results — **without**
pursuing a literal verbatim-script drop-in or bit-exact numerics (both explicit PRD §3
non-goals).

Evidence base: 10-agent audit (2026-06-01). Where a claim cites `file:line`, it comes
from that audit and should be re-verified against current source before implementation.

---

## 1. Re-ranked blockers (functional-equivalence lens)

The literal-drop-in friction (verb names `lmm`/`glmm`, no masking, coefficient-label
format) is **demoted** — it doesn't affect what you can compute. The capability and
correctness gaps are **promoted**.

| # | Blocker | Why it breaks equivalence | Effort | Owner | Evidence |
|---|---------|---------------------------|--------|-------|----------|
| 1 | **GLMM default ≠ glmer** (`pirls_profiled` vs joint Laplace) | Same model, *different coefficients* — deepest non-equivalence | M | upstream + R | `R/glmm.R:48`; `inst/extdata/expected-mismatches.json` (cbpp fixef ~3% rel, grouseticks ~170% rel) |
| 2 | **GLMM second-class** — no `predict`/`confint`/`contrast`/`anova`/`drop1`/`simulate`/`refit` | Whole workflows can't run on a GLMM | M–L | mostly R | `predict.R:152` (refusal); NAMESPACE (inference S3 = `mm_lmm` only) |
| 3 | **`predict()` SE / intervals** (LMM + GLMM) | Blocks plotting/reporting lme4 users expect | M | both | `predict.R:61-101` (NA SE, interval refused) |
| 4 | **Estimation control collapsed** (`mm_control` = `verbose`+`max_feval`) | No recourse when default optimizer struggles; regresses PRD §8.3 | M | upstream-led | `R/mm-control.R:20` |
| 5 | **`update()` missing** | Core model-building/comparison idiom errors | S | R-only | no `update` method in NAMESPACE |
| 6 | **No `broom.mixed` `tidy()`/`glance()`/`augment()`** | Breaks tidy-table / marginaleffects / ggeffects ecosystem | S | R-only | `methods-extract.R:517` (shaped, not registered) |
| 7 | **`na.action`/`subset`/`contrasts`/LMM `offset`** | Routine data handling refused | M | mostly R | `R/glmm.R:60-67`; `lmm()` lacks them |
| 8 | **Few GLMM families** (no negative-binomial) | Overdispersed counts can't be fit | M | upstream | `R/glmm.R:204-208` (closed set: binom/pois/Gamma-log) |

**Deprioritized (not equivalence blockers):** verb aliases, masking, coef-name format,
bit-exact numerics, AR(1)/spatial residual covariance, model-selection. A thin
`lmer <- lmm` / `glmer <- glmm` alias is a cheap optional ergonomics nicety, not a
milestone item.

---

## 2. Ownership split (drives sequencing)

- **R-wrapper-only (this repo, no engine change):** #5, #6, GLMM Wald `confint` (part of #2),
  fixed-effect Wald prediction CI (part of #3), `na.action`/`subset`/`contrasts` + LMM
  `offset` via R-side `model.frame`/response-shift (#7), the default-estimator *decision*
  and call-site labeling (part of #1), all docs/release hygiene.
- **Upstream `mixeff-rs` (peer repo — reserve + coordinate via `mote msg`):** full
  estimation control surface (#4), RE-inclusive prediction variance (#3), routine
  `joint_laplace` certification payload (#1), negative-binomial family (#8), GLMM
  newdata design construction (part of #2).

---

## 3. Milestones

### M1 — R-only quick wins + honesty/hygiene (no upstream; highest equivalence-per-effort)
- [ ] `update.mm_lmm` / `update.mm_glmm`: reconstruct call from stored formula/data/control → refit. [#5]
- [ ] Register `tidy()` / `glance()` / `augment()` for `broom.mixed` (accessors already shaped). [#6]
- [ ] `confint.mm_glmm` (Wald): `est ± z·SE` from the existing GLMM summary covariance payload. [#2]
- [ ] Reconcile the **"drop-in" framing**: README / `docs/index.md` / vignette intros →
      "audit-first; lme4-*equivalent* formulas + extractors, faster and honest," matching PRD §3 / mission.md.
- [ ] Clear **R CMD check** to 0 WARN / 0 NOTE: duplicated `\argument{x}` in
      `man/mm_lmm-methods.Rd`; `.Rbuildignore` the `assessment/` dir, `vignettes/*.css|*.js`,
      and committed `tests/testthat/testthat-problems.rds`.
- [ ] Fix the **parity-harness footgun**: bare unexported `mm_parity_lookup`
      (`R/parity-ledger.R`) under a swallowing `tryCatch` (`helper-parity-scoreboard.R:164`)
      can silently no-op the parity safety net — qualify it / make lookup failure loud.

### M2 — GLMM functional equivalence  *(runs in parallel with M3)*
- [ ] **Certify `joint_laplace` vs glmer** as a *routine* (non-env-gated) parity test; native
      joint-Laplace landed in commit `5db2612`. [#1]
- [ ] **Default decision:** flip `glmm()` to `joint_laplace` once certified, *or* emit a loud
      call-site notice that `pirls_profiled` is a different estimator. Update PRD §glmm either way. [#1]
- [ ] `contrast` / `anova` / `drop1` for GLMM, reusing the Wald-z machinery the emmeans
      bridge already drives. [#2]
- [ ] `predict.mm_glmm` (`type="link"`/`"response"`, in-sample + newdata). [#2, #3]
- [ ] Negative-binomial family — upstream `mixeff-rs` work + R validation. [#8]

### M3 — Estimation control + reporting completeness  *(runs in parallel with M2; upstream-led)*
- [ ] Expand `mm_control()` toward `lmerControl` parity: optimizer choice, maxit/tolerances,
      `start=`. Needs an FFI surface → **mote coordination with `mixeff-rs`**. [#4]
- [ ] Prediction **SE/intervals**: fixed-effect Wald CI is R-side (`X·vcov(β)·Xᵀ`); full
      intervals incl. RE need an engine variance payload (upstream). [#3]
- [ ] `na.action` / `subset` / `contrasts` at the fit front door (R-side `model.frame`);
      LMM `offset` (R-side response shift for identity link). [#7]

### M4 — Adoption surface + release
- [ ] Write the planned-but-missing `lme4-migration.Rmd` (verb/arg map; the
      `groupaphant` vs `group: aphant` coef-name note + `gsub` recipe; "what's
      NA-with-reason and why").
- [ ] Runnable `\examples{}` on key man pages (`glmm.Rd`, `mm_lmm-methods.Rd` — ~80%
      currently have none), `lme4::sleepstudy`/`cbpp`-flavored; fix `mm_lmm-methods.Rd` `\value`.
- [ ] Flesh out `_pkgdown.yml` (reference groups + articles ordering) and rebuild the stale site.
- [ ] Cut R-universe → CRAN once check is clean and the upstream `nlopt` gate resolves.

---

## 4. Sequencing (decision: GLMM + LMM-control in parallel)

```
now ──► M1 (R-only)  ───────────────► release? (R-universe)
        │
        ├─ M2 GLMM equivalence ──┐
        │   (R + upstream cert)  ├──► M4 docs/site ──► CRAN
        └─ M3 control + predict ─┘
            (upstream-led + R)
```

- **M1** and **M4's framing/check items** start immediately and concurrently.
- **M2 and M3 run in parallel** (per scope decision). They are mostly independent: M2 is
  GLMM-centric, M3 is optimizer/data-handling-centric. Open the **upstream `mixeff-rs`
  coordination early** (longest lead time: estimation-control FFI, prediction variance,
  joint_laplace cert, neg-binomial).

---

## 5. Open decisions (defaults assumed; change anytime)

- **Scope:** ✅ DECIDED — M2 and M3 in parallel.
- **Upstream appetite:** *default* = R-first, queue upstream via mote beads in
  `/Users/bbuchsbaum/code/rust/mixeff-rs` and coordinate async. (Alternatives: R-only this
  cycle / full-court upstream now.)
- **Release timing:** *default* = cut an R-universe release after M1 quick wins for early
  visible progress; CRAN after M3 + nlopt gate. (Alternatives: hold for GLMM parity / hold
  for full parity / no release pressure.)

---

## 6. PRD / governance notes

- M1 framing item resolves a real contradiction: README/index say "drop-in replacement"
  while PRD §3 + mission.md say "**not** a drop-in lme4 replacement." Align to the PRD's
  honest framing.
- #4 (estimation control) is **drift against the project's own PRD §8.3**, which specified
  an `lmerControl`-mirroring `mm_control()`. Restoring it is bringing the impl back to spec,
  not new scope.
- `planning/HANDOFF.md` is stale (claims "Phase 1.A closed" while Phases 2–5 largely
  shipped) and PRD §13 Decisions A–D are already resolved — refresh/retire as part of M4.
- Per CLAUDE.md: any `mixeff-rs` change requires explicit path reservation in the peer repo
  and `mote msg` coordination; bump the bundled snapshot only via `PINNED_REV` in
  `tools/vendor-rust.R` + re-run, never hand-edit `src/rust/upstream/`.
