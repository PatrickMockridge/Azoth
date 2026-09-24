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

from azoth.process.reference.component_splitter import component_splitter
from azoth.process.reference.compressor import compressor
from azoth.process.reference.distillation_column import distillation_column
from azoth.process.reference.cooler import cooler
from azoth.process.reference.expander import expander
from azoth.process.reference.filter import filter
from azoth.process.reference.gas_scrubber import gas_scrubber
from azoth.process.reference.heat_exchanger import heat_exchanger
from azoth.process.reference.heater import heater
from azoth.process.reference.manifold import manifold
from azoth.process.reference.mixer import mixer
from azoth.process.reference.pipe import pipe
from azoth.process.reference.pump import pump
from azoth.process.reference.separator import separator
from azoth.process.reference.shortcut_distillation_column import shortcut_distillation_column
from azoth.process.reference.splitter import splitter
from azoth.process.reference.tank import tank
from azoth.process.reference.throttling_valve import throttling_valve

__all__ = [
    "component_splitter",
    "compressor",
    "distillation_column",
    "cooler",
    "expander",
    "filter",
    "gas_scrubber",
    "heat_exchanger",
    "heater",
    "manifold",
    "mixer",
    "pipe",
    "pump",
    "separator",
    "shortcut_distillation_column",
    "splitter",
    "tank",
    "throttling_valve",
]
