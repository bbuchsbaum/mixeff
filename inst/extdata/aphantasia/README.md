# Aphantasia Revision 3 Fixture

This directory contains an anonymized, trial-level fixture generated from
the aphantasia manuscript analysis (revision 3; local manuscript archive,
not distributed).

Files:

- `trials.rds`: trial-level analysis data with stable hashed participant IDs.
- `metadata.rds`: participant metadata keyed by the same hashed IDs.
- `reference.json`: frozen lme4 reference fits and cached manuscript inference
  summaries used by the mixeff reproduction tests.

Regenerate with:

```sh
Rscript tools/build-aphantasia-fixture.R
```

Set `MIXEFF_APHANTASIA_ROOT` to point at a different local manuscript checkout.

Participant anonymization uses a fixed salt label
`mixeff-aphantasia-revision3-v1`, an MD5 hash of `salt:id`, the prefix `p_`,
and the first 16 hexadecimal characters. The raw participant identifiers are
not stored in the fixture.

The fixture keeps only analysis columns needed for the reproduction tests:
participant, bubbled, back_masked, SOA, block_num, trial_image, category,
correct, rt, aphantasia, age, vviq_standard, source, and source_folder.
`source_folder` is retained for the S9 age-matched folder analysis.

The test file keeps ordinary checks fast. Set `MIXEFF_RUN_APHANTASIA=true` to
run the core model refits (primary, sensitivity, intact, combined, RT, S7, and
S9). Primary, intact, and combined use glmm()'s default joint-Laplace
estimator and reach strict lme4 parity with no parity-ledger exemption
(about 2, 2, and 7.5 minutes per fit); the sensitivity and S7/S9 variants
use the profiled fast-PIRLS path, which also clears the strict tolerances
there at ~10-20 seconds per fit. Set `MIXEFF_APHANTASIA_JOINT=false` as a
debug-build escape hatch to route every case through profiled fast-PIRLS —
expect strict-tolerance parity failures for intact and combined. Set
`MIXEFF_RUN_APHANTASIA_STRESS=true` for the S1 random-effects stability
variants (profiled, ~100 seconds in total).
