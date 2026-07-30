# mixeff: Mixed-Effects Models via the 'mixeff-rs' Rust Crate

An R wrapper for the `mixeff-rs` Rust crate.
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fit
linear and generalized linear mixed-effects models from lme4-style
formulas. Fitted models store the compiled model, optimizer result, and
inference metadata as data, so a result can be inspected, saved, and
reloaded; where an inference method is unavailable, the result is either
withheld with a stable reason code or labeled with the method actually
used, never silently swapped. See
[`vignette("intro", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/intro.md)
for an overview and
[`vignette("demystifying-formulas", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/demystifying-formulas.md)
for what random-effects formulas mean.

## See also

Useful links:

- <https://bbuchsbaum.github.io/mixeff/>

- <https://github.com/bbuchsbaum/mixeff>

- Report bugs at <https://github.com/bbuchsbaum/mixeff/issues>

## Author

**Maintainer**: Brad Buchsbaum <bbuchsbaum@research.baycrest.org>
