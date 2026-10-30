# MKL Brown reference repair — 26 October 2026

**The numerical discrepancy is reproduced and diagnosed.** Default lme4 BOBYQA reached an inferior nearly singular reference endpoint. The targeted reference repair is published at `e7697417017e738eb371098a8332ad475740c9d8`; final canonical MKL acceptance passes at `d5fdc54`. Production mixeff estimation code and all existing parity tolerances are unchanged.

## Decisive evidence

[Run 40407299066](https://github.com/bbuchsbaum/mixeff/actions/runs/40407299066) reproduced the original SD mismatch and saved the original fitted pair. It also reproduced the same reference endpoint in a separate R process, so the trigger is not exclusive to the full-check process. The default reference reported optimizer code 0 and normal BOBYQA exit, together with a singular-boundary message. Its stimulus slope SD was 8.175123 and stimulus correlation 0.999999988.

| Covariance endpoint | Stimulus slope SD | Common lme4 REML criterion |
|---|---:|---:|
| Captured default reference | 8.175123 | 302385.8003502963 |
| mixeff | 14.707846 | 302385.6692150966 |
| Tight refit from the exact captured covariance | 14.719220 | 302385.6692145155 |

Restarting the captured reference improves its criterion by **0.1311357808**, much larger than measured evaluator variation around 1e-9. This demonstrates an inferior attainable reference endpoint; it is not a claim about global optimality or flat-profile equivalence. No production engine defect is supported by this comparison.

The saved fixture records the exact relative-Cholesky coordinates, group/coefficient order, input-file MD5 and GitHub artifact provenance. It is replayed on the current fixture rather than requiring a default optimizer run to fail intermittently.

## Repair and regression

The `rt_modality_full` lme4 reference alone now uses BOBYQA `rhobeg=.01`, `rhoend=1e-9`, `maxfun=20000`. The other cases retain their optimizer controls. The new regression guards fixture identity and coordinate order, independently evaluates the captured endpoint with absolute tolerance 1e-6, requires a fresh objective improvement greater than 0.1, checks the repaired objective within absolute 1e-5, and applies the existing Brown parity checks against mixeff. No parity tolerance was loosened.

Reference validation also addresses two demonstrated blind spots. The old helper suppressed all warnings; matching lme4 2.1-0 defaults to skipping derivative calculations above 10,000 observations, while this fixture contains 21,679. Even explicitly calculating derivatives leaves `checkConv`'s separate count guards active. The helper now enables derivatives and removes both count limits for these reference fits, then checks optimizer status, finite likelihood deviance and derivatives and convergence diagnostics. Legitimate singular-boundary notices remain permitted. Validation-only run 40408493648 reproduced the mismatch and saved its pair, confirming why validation alone was insufficient: a singular-boundary label does not establish that an endpoint is the best attainable solution. The scoped tight controls supply the numerical repair.

Invalid references and parity failures save both available fitted objects and diagnostics to RDS. Reference-fitting errors preserve the engine fit and error context. Receipts include session information, controls and messages. The LRT references use the same validated path.

The final focused installed Brown suite passes **128 assertions, 0 failures, 0 test warnings, 0 skips**. Real witnesses cover a one-evaluation optimizer failure, an accepted singular boundary, and an actual failed testthat expectation whose fitted objects reload. A package-startup warning says testthat was built under R 4.5.2; it is distinct from test/optimizer warnings. Independent review verified the boundary/capture behavior, exact-start objective evidence and final two-file repair.

## Bounded multistart evidence

The full Brown REML model was fit from the ordinary lme4 endpoint, an independent unit-diagonal/zero-correlation start whose stimulus slope is 8.2 divided by the residual SD of an ordinary fixed-effects model, and the engine covariance converted to lme4's relative-Cholesky convention. The independent start is a nominal SD 8.2 stress case, not a reconstruction of the failing covariance; sigma is profiled during fitting. The captured endpoint is separately replayed exactly.

Local matched-lme4 fits converge near 14.719, with common-objective spread 3.32e-9. The engine's objective excess is 5.81e-7. Saved derivatives pass full convergence checking with count guards disabled.

[Run 40408493648](https://github.com/bbuchsbaum/mixeff/actions/runs/40408493648) repeats the bounded experiment on MKL with all derivative/convergence checks active, in fresh default-thread and one-thread processes. Both settings produce SD 14.719220,14.719487,14.719005; all three optimizer codes are 0 with no warnings/convergence messages. Their common objective spread is 1.75e-9. This is multistart validation, not a fixed-SD likelihood profile.

## MKL dependency-cache repair

The first canonical repaired-source run 40410870077 could not load lme4: its
cached nloptr binary required `libflexiblas.so.3`, absent from the unchanged MKL
image. The cache key contained the Fedora/R version but no MKL identity. That
run skipped 196 lme4-dependent tests and is not numerical qualification. It also
exposed a missing optional-dependency guard in the new receipt test; that guard
is now fixed. A private-copy witness with an unloadable lme4 verifies 12 data
assertions and 7 expected dependency skips, with no failures or test warnings.
The installed user libraries were unchanged.

The [official r-hub wrapper](https://github.com/r-hub/actions/blob/v1/setup-deps/action.yaml)
uses a shared cache version when no RHUB repository is configured and does not
expose a cache-version override. For MKL only, the workflow now invokes the same
pinned underlying r-lib action with `cache-version: mkl-native-v1` and source
requirements for nloptr/mvtnorm. It requires nloptr, lme4 and emmeans to load
before starting the package check. Other R-hub platforms retain the wrapper.
There is no BLAS substitution or numerical-test bypass. An intermediate run
40413434042 was cancelled during initialization when this durable cache repair
superseded the initial per-package rebuild approach.

## Final qualification and release boundary

Candidate `d5fdc546637848e9fe23314b50532478ba5916fb` includes the numerical repair,
optional-dependency guard and MKL cache isolation. Canonical R-hub
[run 40413879004](https://github.com/bbuchsbaum/mixeff/actions/runs/40413879004)
checked this actual release branch and **passed: Status OK, 0 errors, 0 warnings,
0 notes; 6,194 test assertions pass, 0 fail, 0 warn**. Its 62 skips comprise
44 CRAN-gated tests, 16 source-availability checks and two opt-in slow sweeps.
No Brown test was skipped. The archived Brown test file and endpoint fixture
are byte-identical to the candidate. The preflight confirms lme4 2.1-0 and
Intel MKL 2025.3; the check uses R-devel r90763 on Fedora 44. This run excludes
the PDF manual and vignette-output rebuilding; vignette R code ran successfully.
The source tarball is 6,806,217 bytes, SHA-256
`3312b93e45be02f57a3fa79bb5f554843499bd668a855cde6c29013c330c19ca`.
This is the MKL-built artifact, not a new canonical cross-builder submission receipt.

The core numerical repair at e769741 also passed all five primary R CMD jobs
(Windows, macOS, Ubuntu release/devel and vendor regeneration). Primary CI for
d5fdc54 was still running when this MKL receipt was recorded. The earlier instrumented
failures and dependency failure are retained as evidence. Packaged tests and a
fixture changed; the prior 8065281 full-release tarball and builder receipts are
not byte-identical qualification of this candidate.

Exact raw logs, local test results, artifact identities and summary CSVs are
under `release-evidence/mkl-reference-2026-10-26/`. Large RDS/ZIP artifacts are
retained locally and git-ignored. The endpoint fixture is committed under
`tests/fixtures/brown-mkl-reference-endpoint.json`. GitHub PR 6 still requires an
approving review; no remote-main bypass or CRAN submission is part of this repair.
