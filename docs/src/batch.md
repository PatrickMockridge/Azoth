# The batch API

<!-- Hand-written. Listed in SUMMARY.md by tools/gen_docs.py. -->

**This page was written as a decision record before any batch code existed**, because the
decisions below were cheap to make then and expensive to change once `BatchResult` was
public and someone had written `batch.warnings[i]` into their own code. It is kept in that
voice: each section says what was decided and why, and where the code now lives.

The implementation is `azoth.batch`, and `python/tests/test_batch.py` is what holds it to
the decisions.

```python
import azoth.batch

r = azoth.batch.hydraulics.darcy_weisbach(
    f=[0.02, 0.02], L=[10.0, 25.0], D=[0.05, 0.08], rho=[998.0, 1000.0], v=[2.0, 3.0],
    mu=[1e-3, 1e-3],
)
r.dp                              # array('d', [...]) of SI base magnitudes, pascals
r.units["dp"]                     # 'Pa'
r.warnings[1]                     # element 1's warnings
r.regime                          # ('turbulent', 'turbulent')
r.has_warning(WarningCode.TRANSITIONAL_FLOW)
r.warning_indices(WarningCode.TRANSITIONAL_FLOW)
r.quantities("dp")                # the same values as pint quantities; allocates
```

## What a batch call is

N independent evaluations of one existing calculation, over arrays of inputs, with one
call crossing the language boundary instead of N.

**Inputs are plain numbers in the spec's canonical unit; outputs are SI base magnitudes
with a unit map.** This is the one place the batch API is *not* the scalar API, and it is
a real loss: the scalar API's `UnitMismatchError` has no counterpart here. Passing metres
where the spec declares millimetres is a thousand-fold error that nothing catches.
`orifice_flow` is the live example - its `d` is in **millimetres**, so `d=[50.0]` means
50 mm. It is stated there, in the batch module docstrings, and here.

Converting a sequence of `pint` quantities would cost one `pint` operation per element,
the same order as the per-element boundary crossing this API exists to remove, so the
conversion is **one factor applied to a whole array**, taken from `pint` rather than
written down (`batch/_core.py::_si_per_canonical`).

## Decision 1: the batch API is a loop over the scalar kernels, not a new kernel

**Decided.** Both implementations loop over their existing scalar function.

This is the decision that matters most, and the reason is not performance. `README.md`
and `CONTRIBUTING.md` both stake this project on there being exactly **two**
implementations that check each other: Python is the reference and Rust is the core,
and every spec case runs through both. A vectorised Rust kernel would be a *third*
implementation — a different shape, a different arithmetic order, and one nobody
cross-checks — so the claim that two independent implementations agree would quietly
stop being true while every test still passed.

That is a large price for an optimisation. And it is not even the optimisation that
matters: the README's own argument is that the PyO3 call overhead exceeds the cost of
the arithmetic, so removing N−1 boundary crossings is where the win is. Looping in Rust
over the same scalar function captures it.

If a genuinely vectorised kernel is wanted later, it ships as its own separate,
separately validated step whose agreement test is against *both* existing scalar
implementations. It does not arrive as a batch feature.

## Decision 2: struct of arrays, with warnings positional

**Decided.** One array per output field, and a parallel structure for warnings:

```python
r = azoth.batch.hydraulics.darcy_weisbach(f=[...], L=[...], D=[...], rho=[...], v=[...])
r.dp                 # array of SI magnitudes, Pa
r.warnings[7]        # the tuple of warnings for element 7
r.has_warning(WarningCode.TRANSITIONAL_FLOW)      # any element
r.warning_indices(WarningCode.TRANSITIONAL_FLOW)  # which ones
```

**Rejected: an array of result objects.** It allocates N dataclasses *and* N pint
quantities, which is most of what the boundary crossing cost. The whole point is to
cross once.

**Rejected: bare arrays of floats with no result object.** It loses the field-name
contract every caller already knows, and leaves nowhere to put per-element warnings.

Per-element warnings are why the shape is struct-of-arrays rather than a table:
`warnings[i]` is element `i`'s warnings, empty for a clean element. `Warning` itself
does **not** gain an index field — `azoth/core/warnings.py` is a cross-language
contract asserted structurally by the test helpers, and widening it for a batch feature
would leak into the scalar API.

**Settled: an array of a dimensioned quantity is `array('d')` of SI magnitudes plus a
unit map.** The three candidates were:

| Option | What `r.dp` is | Cost |
|---|---|---|
| **SI magnitudes plus a declared unit** | `array('d')` in Pa, with `r.units["dp"] == "Pa"` | Not a `pint` quantity, so it must not be mistaken for one — but it is the reason the API is fast |
| An array of `pint` quantities | `list[Q]` | N objects, which is what the batch shape exists to avoid |
| `pint` with an array backend | `Q[np.ndarray, "Pa"]` | Needs numpy as a dependency, which decision 4 refuses |

The first, with `BatchResult.quantities(field)` as the explicit way back to `pint` for a
caller who wants it. `quantities` allocates one quantity per element and says so in its
docstring, so choosing it is choosing to spend what the batch shape saved.

**Also settled, and not in the original three: an enum output is a label column.** A flow
regime is `laminar`, `transitional` or `turbulent`, and there is no number it should be.
Encoding it as an index into a table would make a caller look the mapping up to read a
value and would need a sentinel for "not available" that is not one of the three real
answers; `None` in a label column says that unambiguously.

**A column's kind comes from the field's declared type, never from its values.** This was
not in the original note and had to be decided during implementation, because the obvious
implementation — look at the first element — is wrong. `darcy_weisbach`'s `re` is
`float | None`: sniffing values would make it a *numeric* column when viscosity was
supplied and a *label* column when it was not, so the field's kind would depend on which
argument the caller passed, and the two backends would disagree. A numeric field is always
numeric; an absent value is `NaN`, which is what the extension sends.

## Decision 3: errors are fail-fast, and both backends agree

**Decided.** If element 7 of 1000 violates an `error`-severity bound, the whole call
raises. No partial result, no sentinel, no zeros in the gap.

A partial result with silent gaps is the failure this project is organised against: a
caller who receives 999 numbers and does not notice which one is missing has a wrong
answer that looks complete. The alternative — an `error_mask` alongside the arrays — is
a reasonable design for a later version if callers need it, and is noted here so that
choosing it later is a decision rather than an oversight.

The reason this needs stating, rather than being obvious: **both backends must make the
same choice.** The cross-language agreement test compares warning lists structurally, so
a Python implementation that returned a partial result where Rust raised would surface
as a confusing mismatch rather than as the design question it is.

Range *warnings* are unaffected: an element outside a validated range still returns its
value and carries a warning, in batch exactly as in scalar.

## Decision 4: no numpy dependency

**Decided.** Batch functions return `array.array('d')`, which is
buffer-protocol-compatible, so `numpy.asarray` reads it without a copy. "numpy interop"
and "a numpy dependency" are separable, and this project takes the first and refuses the
second. `CONTRIBUTING.md` requires new dependencies to be justified, and the argument for
this one fails: the `numpy` crate's C-API access is version-sensitive and awkward next to
the `abi3-py312` feature that lets this project ship one wheel covering every supported
CPython, and that is a large packaging cost for a marshalling optimisation that does not
need it.

**What was built differs from the plan here, and the difference is worth recording.** The
plan was `PyBuffer<f64>`, accepting any buffer-protocol object for free. What shipped is
plain iteration: `sequence()` in `batch/_core.py` accepts any iterable of numbers, so a
numpy array is accepted — numpy arrays are iterable — but element by element rather than
as a memcpy. The buffer path was not taken because it would not have paid for itself yet:
the Python side materialises a list before the call regardless (it has to, to apply the
unit factor), so the bytes are already being moved one float at a time, and switching the
Rust side to a buffer would leave that untouched.

It remains the right optimisation if profiling ever says extraction dominates. It is
recorded as a deviation rather than quietly dropped, because "we accept buffers" and "we
accept iterables" are different claims and only one of them is true.

What *is* tested, without numpy installed: that the inputs accept any iterable, and that
the output column really is buffer-protocol-readable — see
`test_inputs_take_any_iterable_and_outputs_are_buffers`.

## Decision 5: no broadcasting

**Decided.** Every input is a sequence of the same length. A scalar where an array is
expected is an error, not a value to broadcast.

Broadcasting is convenient and it is also how a batch call silently does something
other than what the caller pictured — a length-one list where a length-N list was meant
produces a plausible array of answers rather than an error. A caller who wants one
element repeated can say so.

This is the same instinct as the scalar API's refusal to accept a bare float where a
quantity is expected, which is documented as the mistake the library exists to make
impossible.

## One calculation has no batch form

`hydraulics.crane_k_factors` takes `fittings`: a list of registry ids, declared in the
spec with no unit, because a fitting is a category and not a magnitude. There is no
column shape for it — N elements each carrying their own list of names is a sequence of
sequences, and one shared list would compute the same answer N times, which is not what a
batch API is for.

So `azoth.batch` refuses it by name rather than half-supporting it, and
`test_the_excluded_set_is_exactly_crane_k_factors` asserts the excluded set is exactly
that one. A future calc with a categorical input lands in that assertion and has to be
argued for instead of quietly joining the exclusion.

## What this is not

- **Not new mathematics.** No batch-specific calc, no batched spec, no batched range
  check. Each element goes through the same checks the scalar call applies.
- **Not a performance promise about the arithmetic.** The arithmetic cost is identical.
  What disappears is the per-element boundary crossing.
- **Not a flowsheet or a solver.** `docs/src/index.md` is explicit that composition is
  the caller's job; a batch call is N independent evaluations of one calculation, not a
  connected system.
- **Not automatic for a new calc, in the sense of being generated.** The plan expected
  the per-calc batch wrappers to fall out of the registry. They did not: a wrapper's
  signature, result class and docstring carry judgement the spec does not, so each is
  hand-written, and a coverage test
  (`test_every_batchable_calc_has_a_batch_arm`) is what makes omitting one fail the build
  naming the calc and the file. That is the same mechanism
  `python/tests/test_registration_completeness.py` uses for the rest of the wiring, and
  the same honest position: a check that fails loudly, rather than a generator that would
  have to invent the API.
