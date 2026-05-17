# Upstream mixedmodels parity mismatch report

Use this only for concrete Rust-engine bugs found while testing mixeff parity.
Feature requests, design/API requests, and prioritization questions should stay
in mixeff planning until discussed with the user.

## Summary

- Title:
- Source: mixeff parity test / manual reduction / generated design
- Severity: blocker / regression / numeric discrepancy / metadata discrepancy
- One-line diagnosis:

## Case

- Dataset or generator:
- Generator seed, if simulated:
- Formula:
- REML: TRUE / FALSE
- Family/link, if GLMM:
- Reference packages:
  - lme4:
  - lmerTest:
  - pbkrtest:
- mixeff commit:
- mixedmodels commit:
- Platform/R/Rust:

## Compared Fields

| field | tolerance | reference value | observed Rust value | absolute/relative diff | status |
| --- | ---: | ---: | ---: | ---: | --- |
|  |  |  |  |  |  |

## Expected Behavior

State the reference behavior and why it is the right comparison target. Include
whether the comparison is against lme4, lmerTest, pbkrtest, a closed-form
identity, or a documented Rust contract.

## Observed Behavior

State the Rust-backed mixeff behavior. Include status labels, reasons,
reliability labels, and any structured details that disagree with the contract.

## Minimal Reproducer

```r
library(mixeff)

# Include data construction or a classic dataset load.
# Include the exact formula, REML flag, method, compared field, and tolerance.
```

## Reduction Notes

- Smallest dataset/design that reproduces:
- Does the mismatch persist under ML and REML:
- Does the mismatch persist with random slopes removed:
- Does the mismatch persist with correlations removed:
- Boundary/singularity status:

## Filing Notes

- Upstream mote id:
- Link or path to failing mixeff test/output:
- Related mixeff mote id:
