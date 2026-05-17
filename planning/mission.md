# Mission — `mixeff`

`mixeff` is an R package that makes the `mixedmodels` Rust compiler usable from R
without compromising what the compiler is for.

## What we do

We give applied analysts in psycholinguistics, ecology, neuroscience, education,
clinical trials, and the social sciences a way to fit linear and generalized linear
mixed-effects models that is:

- **Familiar** — lme4-style formulas (`y ~ x + (x | g)`), S3 methods (`fixef`,
  `ranef`, `VarCorr`, `predict`, `simulate`, `anova`), and the extractors users
  already know.
- **Honest** — every reported number is backed by a named method, a documented
  status, and a versioned artifact. Where inference is unavailable, we return `NA`
  with the reason, not a plausible substitute.
- **Auditable** — `audit()`, `changes()`, `diagnostics()`, `parameterization()`,
  `optimizer_certificate()`, `inference_table()`, and `reproducibility()` are
  first-class verbs, not afterthoughts.
- **Durable** — the fitted object round-trips through `saveRDS()` and survives
  cluster restarts, forks, and process death without depending on a live Rust
  handle to remain interpretable.

## Who we serve

The user we keep in mind is a working scientist who needs to *defend* a mixed-model
result — to a co-author, a reviewer, a regulator, or themselves six months later.
They are fluent in lme4 syntax, comfortable with R's S3 ecosystem, and tired of
warning messages that vanish without explanation.

We do not optimize for users who want a black-box prediction engine, nor for
methodologists prototyping novel covariance structures — those audiences have
better-suited tools.

## How we work

Five operating principles govern every design decision:

1. **Rust owns model semantics; R captures intent.** Decisions about convergence,
   identifiability, covariance reductions, and inference availability live inside
   the compiler. R packages user input, presents results, and never overrules the
   engine.
2. **The JSON artifact is the source of truth.** The external pointer is a cache.
   If the cache and the artifact disagree, the artifact wins. Every printed claim
   must trace back to an artifact field.
3. **No fake certainty.** When a method is unavailable for a contrast, a boundary
   fit, or a missing degree of freedom, we return `NA` with the Rust reason code.
   We do not fabricate.
4. **No silent surgery.** Every model reduction, refusal, or canonicalization
   crosses the boundary as a structured diagnostic. Users see what the engine did,
   why, and what the requested → effective transition looked like.
5. **Audit before fit, and after.** `compile_model()`, `explain_model()`, and
   `audit_design()` work without optimization. The first-class workflow is
   *explain, fit, audit, contrast* — not *fit, debug, regret*.

## What we are not

- **Not a drop-in `lme4` replacement.** We share the formula language and many
  generic names; we do not share the philosophy or the numerical conventions.
- **Not bit-exact with `lme4`.** Statistical equivalence within documented
  tolerances on parity datasets is the bar; matching the sixth decimal of a
  profile-deviance evaluation is not.
- **Not a research playground.** GAM smooths, residual covariance structures,
  multivariate cross-outcome residuals, and regularized model search are out of
  scope for v1. The compiler may grow them later; this package will follow, not
  lead.
- **Not a masking package.** We do not displace `lme4::lmer` or `lme4::glmer` on
  attach. Both packages can coexist in a session.

## Success criteria

We will know we have succeeded when:

- A working scientist can fit a model, save it, hand the `.rds` to a collaborator,
  and have that collaborator inspect the audit trail without the original session,
  the original optimizer, or the original Rust handle.
- A reviewer can ask "what method produced this p-value, and was it reliable?" and
  the answer comes from `inference_table(fit)`, not from the analyst's memory.
- A boundary or singular fit is a *fact about the model*, not an embarrassment to
  be re-fit until it disappears.
- The package is on CRAN, building reproducibly from vendored Rust sources, on
  every platform R itself supports.
