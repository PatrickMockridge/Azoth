"""The pure-Python reference implementation for the process namespace.

Same contract as :mod:`azoth.thermal.reference`: always importable, with or without
the compiled extension, and written to mirror the Rust line for line so a reviewer can
read the two side by side against the port source.

**This namespace had no reference until P11.** The unit-operation tier that existed before
was Rust-only behind a binding, so the library's two-kernel rule — every model written
twice, `test_cross_impl` comparing them — was not applied to it. The modules here are that
rule arriving.
"""

from __future__ import annotations

from azoth.process.reference.mixer import mixer
from azoth.process.reference.pump import pump
from azoth.process.reference.splitter import splitter

__all__ = [
    "mixer",
    "pump",
    "splitter",
]
