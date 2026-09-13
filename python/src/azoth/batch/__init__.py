"""The same calculations, over arrays.

One call evaluates a calculation N times, crossing the language boundary once instead of
N times. The README has named this as the missing performance story since the beginning;
`docs/src/batch.md` records the decisions behind it, and most of that page is about what
this deliberately does *not* do.

* :mod:`azoth.batch.hydraulics` - the hydraulics calculations, batched
* :mod:`azoth.batch.thermal` - the thermal calculations, batched
* :mod:`azoth.batch.eos` - the equations-of-state calculations, batched

# A batch call is a loop over the scalar kernels

Not a vectorised kernel. `README.md` and `CONTRIBUTING.md` both stake this project on
there being exactly two implementations that check each other, and a vectorised Rust
kernel would be a third - a different shape and a different arithmetic order, checked by
nothing. The optimisation does not need one either: the boundary crossing is the cost,
not the arithmetic, so looping captures the win.

# One calculation has no batch form, and it is not an omission

`hydraulics.crane_k_factors` takes `fittings`: a list of registry ids, declared in the
spec with no unit because a fitting is a category and not a magnitude. There is no column
shape for that - N elements each carrying their own list of names is a sequence of
sequences, and one shared list would compute the same answer N times. So the batch API
refuses it by name rather than half-supporting it, and `tests/test_batch.py` asserts the
excluded set is exactly that one, so the exclusion stays a decision.

``from azoth.batch._core import batchable`` lists what is covered.
"""

from __future__ import annotations

from azoth.batch import eos, hydraulics, thermal

__all__ = ["eos", "hydraulics", "thermal"]
