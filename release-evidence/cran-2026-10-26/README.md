# Local CRAN readiness receipts, 26 October 2026

`artifact.json` binds the checked runtime source, vendored engine and tarball digest. `check.log` is the final 0E/0W/1N as-CRAN result; `check-tests.log` contains its CRAN-mode suite. `extended-suite.log` records the full slow-parity/core-aphantasia suite against the installed tarball and lme4 2.1-0, plus lint and live HTML checks.

The auxiliary script copied one skip-policy string incorrectly, so its final policy assertion exited nonzero after all tests and rendering had completed. `final-result-audit.R` reads the unchanged canonical allowlist directly and validates the saved results: 6851 pass, zero failures, three allowed skips, complete HTML. `final-result-audit.log` records this correction and all warning categories. The original script/log are retained rather than overwritten.

Intermediate failures are intentionally retained and are superseded by the final receipts. The old gate's live-render PASS checked only file existence and hid a truncated albersdown C-locale document; the final portable render is complete and the gate now checks the closing HTML tag.

The tarball and full testthat result RDS are retained locally here but excluded from git. Their SHA-256 values are in `artifact.json`. The scripts retain the actual temporary paths used during execution; they are run receipts, not general installation tooling. No remote publication or CRAN submission was performed.
