# OSF Brannick & LaCroix 2025 — blinking / music GLMM+LMM corpus

An in-the-wild check that `mixeff` reproduces the published mixed models from
Brannick & LaCroix (2025), *Blinking indexes dynamic attending during and
after music listening*, Scientific Reports 15:27262
<https://doi.org/10.1038/s41598-025-12200-6>. Tracked by mote
**bd-01M0842MM295RFXXN4C1ZDE9CV**.

## Source

OSF project **[gf3km](https://osf.io/gf3km/)** (*Music and Attention*),
view-only token `0c0e8abeb85c4f97868a922eff278fa6`. Analysis component
**[7ryjh](https://osf.io/7ryjh/)**:

| OSF file | guid | role |
| -------- | ---- | ---- |
| `8-25_Database.csv` (131 MB) | [9x5rz](https://osf.io/9x5rz/) | cleaned ANT blink / RT time series |
| `MUSICxANT1-2_Script_10-9-24.R` | [35dfe](https://osf.io/35dfe/) | Oct 2024 analysis (WOI search + model selection) |
| `SR_R1_SupplementalCode.R` | [685eda7ea310effdd5d1adee](https://osf.io/7ryjh/) | June 2025 revision (CCF + paper-like EC / RT) |
| `MOA_ParticipantRawData.xlsx` | [tdsjr](https://osf.io/tdsjr/) | participant raw workbook (parent node) |

Local copies of the two R scripts live in this directory as
`upstream-MUSICxANT1-2_Script_10-9-24.R` and
`upstream-SR_R1_SupplementalCode.R`.

## What is (and is not) reproducible from the public node

The paper fits binomial GLMMs on blink probability and Gaussian LMMs on
reaction time / blink onset. Random-effect structures were chosen by a
forward “best path” search (Barr et al. 2013); test statistics were
generated only after the winning model was locked.

**In the committed corpus (Cue-Onset slices only):**

| paper table | model | RE | nobs |
| ----------- | ----- | -- | ---- |
| Table 3 | `Blink ~ Music * ANT * Cue + (1 + ANT \| P)` alerting | random intercept + ANT slope | 34,257 |
| Table 4 | same formula, orienting (Center vs Valid) | `(1 + ANT \| P)` | 33,449 |

Factor coding uses High-dynamic / Happy as the music reference, ANT T1 as
the session reference, and Double (alerting) or Center (orienting) as the
cue reference. Executive-control and RT slices are **not** committed: the
public event-level file does not recover Table 5/6 nobs.

**Not on the public node, so not in this corpus:**

- `MusicSaccadeReport.csv` / `MusicTrialBinnedData(100ms)_8-26-24.csv` —
  listening-period blink series used for Table 2
  (`Peak_CCF ~ poly(Window_Start, 3) * Music + (1 | P)`).
- `Mendelson10.wav` / `Shostavich10.wav` — audio needed to rebuild spectral
  novelty and the CCF predictor in Table 6.

Table 6’s published formula also includes `CCFmean`, `∆Blink`, and
`Blinktime`. Those columns cannot be rebuilt without the missing listening
file, so Table 6 is out of the corpus.

## Script vs paper (do not collapse these)

The October 2024 script and the published tables are not the same analysis:

- Cue-Onset search windows in the script are 2400–2700 ms and 2750–3050 ms;
  the paper locks 2401–2902 ms.
- Alerting model selection in the script prefers
  `(1 + ANT_Cx + CUE_TYPE | P)`; Table 3 reports `(1 + ANT | P)`.
- Music reference in the script starts as Silence; the paper tables use High
  / Happy, then relevel to produce the extra Low-ref / Silence-ref rows.

This corpus reconstructs the **Table 3–4 Cue-Onset designs**. The script is
kept as provenance for the forward-search path.

The paper’s own Table 3 three-way High × T2 × No-cue term is internally
inconsistent: the table is β = −0.877, z = −2.389, p = .017, while the
Results prose quotes z = −2.657, p = .008.

## Replication findings (2026-08-17)

`8-25_Database.csv` is an **event-level** file (662,898 rows: fixations,
saccades, blinks, RT), not the 7M-sample stream described in the methods.

| check | result |
| ----- | ------ |
| Table 3 nobs | **exact** 34,257 at Cue-Onset 2401–2902 + `ALERTING` |
| Table 4 nobs | **exact** 33,449 at the same WOI + `Valid_Center` |
| Table 5 nobs | **not recovered** (64,323 Neutral/Incongruent events from 2903 ms to trial end vs paper 116,994) |
| Table 6 nobs | **not recovered** (5,683 RT/BlinkStart rows vs paper 13,176) |
| Table 3 coefficients | **not recovered.** Same formula, High-ref coding, and nobs give intercept −4.498 (paper −4.487) but Cue-No 0.618 vs 0.781 and several Music terms sign-flipped relative to the table. Adding the script’s extra Cue slope does not restore the table. |
| Table 4 Cue-Valid | close (−0.418 vs −0.437); ANT T2 is not (−0.554 vs −0.051). |
| mixeff vs glmer | **holds** on the Table 3 slice: `joint_laplace` max \|Δfixef\| ≈ 8e-4, ΔlogLik ≈ 2e-3. |

The public node therefore supports an in-the-wild **engine-parity** case on
the published alerting/orienting designs. It does not support a digit-level
reproduction of the printed Table 3–6 estimates. Table 2 remains blocked on
the missing listening-period files.

## Files

| file | purpose |
| ---- | ------- |
| `reconstruct.R` | download `8-25_Database.csv`, slice the paper WOIs, write slim fixtures. Needs network. |
| `published-models.R` | thin wrapper that sources the testthat helper. |
| `../../tests/testthat/helper-osf-brannick-lacroix.R` | shared formulas, coding, and table anchors (ships with tests). |
| `parity-harness.R` | offline `glmer`/`lmer` vs `mixeff` comparison on the committed fixtures. |
| `upstream-*.R` | local copies of the OSF analysis scripts. |
| `../../tests/fixtures/osf_bl_alerting.csv.gz` | Cue-Onset alerting slice (~10 KB). |
| `../../tests/fixtures/osf_bl_orienting.csv.gz` | Cue-Onset orienting slice (~10 KB). |
| `../../tests/testthat/test-osf-brannick-lacroix-2025-parity.R` | gated regression test. |

Regenerate fixtures from the package root:

```sh
Rscript data-raw/osf-brannick-lacroix-2025/reconstruct.R
# or, if the 131 MB CSV is already local:
Rscript data-raw/osf-brannick-lacroix-2025/reconstruct.R /path/to/8-25_Database.csv
# OSF_BL_DATABASE=/path/to/8-25_Database.csv Rscript data-raw/osf-brannick-lacroix-2025/reconstruct.R
```

The 131 MB source CSV is cached under `tmp/osf-brannick-lacroix-2025/`
(gitignored) and is not part of the corpus. Only the two Cue-Onset slices
(~20 KB gzipped) are committed.
