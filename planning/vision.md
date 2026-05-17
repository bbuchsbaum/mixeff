# Vision — `mixeff`

A future in which fitting a mixed model in R is an act of *evidence*, not an act of
faith.

Today, the standard tools (`lme4`, `glmmTMB`, `nlme`) ship a culture of quiet
compromise: convergence warnings that scroll off-screen, singular fits printed without
ceremony, p-values produced by methods the software never names, optimizer choices
buried six function calls deep, and a fitted object that contains numbers but not
provenance. Practitioners learn to read between the lines, copy-paste reassuring
diagnostics from blog posts, and hope that "boundary (singular)" is not the same as
"wrong."

`mixeff` exists for the next generation of that work.

We see a world in which:

- **Every printed claim has a provenance.** A coefficient, a standard error, a
  variance component, a predicted value — each can be traced to a versioned artifact
  produced by a named compiler at a known crate version against a known schema.
- **The model object is a record, not a black box.** A fitted model survives
  `saveRDS()`, a cluster restart, a knitr cache, or a colleague's laptop without
  losing its audit trail. Reproducibility is a property of the artifact, not a
  discipline imposed on the user.
- **Refusal is a feature.** When the design is non-identifiable, when inference is
  not defensible, when the optimizer landed on a boundary the requested method
  cannot certify, the software says so — by name, with a stable diagnostic code, and
  without manufacturing a number to fill the gap.
- **The audit happens *before* the fit.** Users learn whether their data can support
  the model they wrote down before they wait for an optimizer. "Explain before fit"
  is the first-class workflow, not a debugging step after a warning.
- **The lme4 formula on-ramp is preserved, but the ceiling is higher.** Anyone
  fluent in `(x | g)` can use the package on day one; anyone who needs to defend a
  mixed-model analysis to a methods reviewer, a regulator, or their future self can
  do so with structured artifacts in hand.

We are not building a faster `lme4`. We are building the package that practitioners
will reach for when they need to *show their work* — to a co-author, a referee, an
auditor, or a model-comparison committee — and know that what they show is the same
thing the software actually did.

The long arc: applied mixed-model practice should look, ten years from now, less
like "I ran lmer and it converged" and more like "I ran the model in confirmatory
mode against schema v1.4; the audit refused term X for insufficient grouping levels;
the certificate reports interior convergence; Satterthwaite degrees of freedom were
available for these contrasts and unavailable for those." `mixeff` is one of the
packages that gets us there.
