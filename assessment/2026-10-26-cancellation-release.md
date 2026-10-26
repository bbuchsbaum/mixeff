# Restored cancellation and release qualification — 26 October 2026

**Decision: hold CRAN submission.** The repaired candidate is qualified by local and primary hosted checks, but one MKL reference-parity failure remains unexplained. Publication is independently blocked by GitHub's required approving review on [PR 6](https://github.com/bbuchsbaum/mixeff/pull/6). A merge attempt returned HTTP 405; branch protection was not bypassed.

## Candidate and repairs

Runtime/build source: `806528106473594aa7342c6585c2b56be5ed2649`. Engine pin: `004a5d9b12a5d71f483c6c5fc15dbf93ffdb8b7b`. Final tarball: 6,807,157 bytes, SHA-256 `a742c9973c1ed64a04d7d3aad3f1088cf53ff4c37cd14f1758c5539995ec0881`. Subsequent assessment, receipt and submission-draft updates are excluded from the package build.

The October remediation and current interface are integrated in PR 6. The Rust engine exposes callback setters that attach or clear a host callback without refitting or changing stored fitted state. The R bridge reattaches its interrupt callback after restoring LMM/GLMM snapshots. [Rust PR 13](https://github.com/bbuchsbaum/mixeff-rs/pull/13) merged as `33c8639d4e4eaf055dfa09ac36794c5b11912332`; the pinned repair is reachable from upstream main.

Native tests verify serial/two-worker interruption, caller-thread callbacks, byte-identical snapshots and successful execution after clearing the callback. Actual macOS R subprocess SIGINT tests cover RDS-restored LMM/GLMM with one/two threads, unchanged baselines and a subsequent successful native call. Windows native callback tests pass; R-console SIGINT delivery on Windows is not certified by these POSIX subprocess tests.

Warning cleanup fixes literal documentation matching, deprecated matrix assertions, undersized bootstrap examples and incomplete benchmark HTML. Contrast validation rejects nonnumeric/complex input before coercion. Release gates now reject unhandled warnings. Archive inspection and installation use R's internal support, removing an undeclared external xz dependency found on four R-hub images.

## Executed qualification

- Final canonical release-candidate gate: **14 PASS, 0 FAIL, 0 SKIP** in 2,198 seconds, under a PATH that rejects external xz and tar access to the vendor archive. Vendor inspection validates 4,401 entries.
- Extended installed suite: **6,881 passed, 0 failed, 0 warnings, 3 allowlisted skips** in 997 seconds. Slow parity/core aphantasia enabled; live vignette and separate offline installation passed.
- Local macOS arm64/R 4.5.1 as-CRAN: **0 errors, 0 warnings, 1 expected incoming note**. The separate offline gate rebuilds a different tarball from the same source; it is not byte identity evidence for the release artifact.
- Primary GitHub workflows at 8065281 all pass: [R CMD check](https://github.com/bbuchsbaum/mixeff/actions/runs/40390580451), [release gates](https://github.com/bbuchsbaum/mixeff/actions/runs/40390580478), [pkgdown](https://github.com/bbuchsbaum/mixeff/actions/runs/40390580762), [coverage](https://github.com/bbuchsbaum/mixeff/actions/runs/40390580456). R-code coverage: 94.2%.
- Final [mac-builder](https://mac.R-project.org/macbuilder/results/1793045109-362c7021c79baa3a/): **Status OK**, including PDF manual, R 4.6.1 Patched. The API flavor was r-devel; the actual runtime version governs this claim.
- Upstream Rust [CI](https://github.com/bbuchsbaum/mixeff-rs/actions/runs/40386105800): 16 executed jobs pass; two scheduled external-parity jobs were skipped. Local callback/default/R-feature, snapshot and bootstrap tests pass.

## R-hub findings and remaining uncertainty

[Initial executed court](https://github.com/bbuchsbaum/mixeff/actions/runs/40389108548) tested bf612d6; 8065281 changes archive handling only, with identical estimation runtime and tests. Ubuntu GCC/Clang, gcc16, Windows, both macOS architectures, m1-san, donttest, lto and atlas pass. All four archive-affected platforms (clang-asan, clang-ubsan, nold, rchk) pass [the repaired retry](https://github.com/bbuchsbaum/mixeff/actions/runs/40390625134). gcc-asan/nosuggests fail before checking mixeff because lmerTest cannot load its lme4 dependency's libRlapack.so. Valgrind remains pending at handoff. Setup failures are not package passes.

MKL job 121805238156 genuinely failed one Brown `rt_modality_full` SD comparison: mixeff approximately 14.7 versus lme4 8.2, alongside 6,156 passes and zero warnings. R 4.7.0, lme4 2.1-0, Matrix 1.7-6, MKL 2025.3. No test tolerance or skip was changed.

Local diagnostics with matched lme4, a [standalone MKL probe](https://github.com/bbuchsbaum/mixeff/actions/runs/40395281983), and an [instrumented full MKL check](https://github.com/bbuchsbaum/mixeff/actions/runs/40398326331) all agree near SD14.72. The instrumented check captured its original test-process fits before assertions and finished Status OK. Reference SD14.719799 versus mixeff14.707846; their common lme4 REML objectives are 302385.6692145177 versus 302385.6692150966. Tight reference refits corroborate this agreement. The original failure artifact retains the assertion/trace but no original fitted objects, so its 8.2 endpoint cannot be evaluated retrospectively.

These are non-reproductions, not an explanation of the failure. Independent numerical review supports merging the verified repairs while retaining the CRAN hold. Do not attribute the discrepancy to a flat ridge or weaken tolerances. One bounded next diagnostic is reference fits from predetermined default, engine and SD8.2 starts under MKL, comparing common objectives and statuses. Hold resolution requires a diagnosed recurrence/correction or an explicit release decision accepting the residual uncertainty; accumulating passing reruns is insufficient.

## Publication and external gates

PR 6 is open. GitHub requires one approving reviewer with write access before merging. Root main retains the earlier local remediation integration; the full repaired candidate is published on `fix/restored-cancellation`. No safeguards were changed and no CRAN submission occurred.

Final win-builder R-release upload received FTP226; no result email was found at handoff. Final R-devel upload received FTP550 because an older same-name job remained queued. A bounded 30-minute queue watcher timed out without retrying; no watcher remains active. The earlier tarball's submissions do not qualify the final packaging change.

R-universe still reports source commit 6e9a62a. Its native jobs passed that older revision; its wasm job failed in cargo-version wrapper parsing. These do not qualify this candidate. Recheck native R-universe jobs after the approved merge.

Tracking: repair/publication `bd-01M5WCP9PDSE6MKK9QB14KDEBR`; MKL investigation `bd-01M5WXZQV8HV23ZEERAGRYQEGB`; release qualification `bd-01KRCKFC72WP98JKSAVKXXB668`. Raw logs, artifact hashes, rejected attempts and builder receipts are under `release-evidence/cancellation-2026-10-26/`. Large tarballs/RDS/ZIP files are retained locally and git-ignored. Diagnostic source/captures are on branch `audit/mkl-parity-20261026` at 3f8051b; that branch is not part of the release change.
