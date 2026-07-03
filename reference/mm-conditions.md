# Typed conditions raised by mixeff

All structured errors and informational signals raised by the package
inherit from the `mm_condition` class so callers can catch them
generically (`tryCatch(..., mm_condition = handler)`), or by a more
specific subclass for finer-grained handling.

## Details

The package ships `mm_formula_error`, `mm_schema_error`, and the
`mm_bridge_error` fallback for untagged Rust errors; `mm_data_error`
(raised by
[`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
on data shape / type problems); `mm_fit_error` for fit construction /
optimization failures; and `mm_inference_unavailable` for inference,
extractor, or prediction requests the engine cannot certify on this fit.
Further classes are `mm_not_identifiable`, `mm_design_refusal`, and
`mm_fit_not_optimized` (see PRD §8.1).

`mm_arg_error` is distinct from all of the above: it signals that the
*caller passed an invalid or malformed argument* (wrong type, value out
of range, unknown option name, missing required field), independent of
any model, engine, schema, or data semantics. Catching
`mm_inference_unavailable` must mean "the engine refuses to certify
this", not "you typed it wrong" — so plain argument validation raises
`mm_arg_error`, never the domain-refusal classes.
