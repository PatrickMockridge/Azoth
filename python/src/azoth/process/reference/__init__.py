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

from azoth.process.reference.absorption_column import absorption_column
from azoth.process.reference.component_splitter import component_splitter
from azoth.process.reference.compressor import compressor
from azoth.process.reference.cooler import cooler
from azoth.process.reference.distillation_column import distillation_column
from azoth.process.reference.ejector import ejector
from azoth.process.reference.expander import expander
from azoth.process.reference.filter import filter
from azoth.process.reference.flare import flare
from azoth.process.reference.gas_scrubber import gas_scrubber
from azoth.process.reference.heat_exchanger import heat_exchanger
from azoth.process.reference.heater import heater
from azoth.process.reference.manifold import manifold
from azoth.process.reference.mixer import mixer
from azoth.process.reference.packed_column import packed_column
from azoth.process.reference.pipe import pipe
from azoth.process.reference.plug_flow_reactor import plug_flow_reactor
from azoth.process.reference.pump import pump
from azoth.process.reference.rate_based_packed_column import rate_based_packed_column
from azoth.process.reference.separator import separator
from azoth.process.reference.shortcut_distillation_column import shortcut_distillation_column
from azoth.process.reference.splitter import splitter
from azoth.process.reference.stirred_tank_reactor import stirred_tank_reactor
from azoth.process.reference.stripping_column import stripping_column
from azoth.process.reference.tank import tank
from azoth.process.reference.three_phase_separator import three_phase_separator
from azoth.process.reference.throttling_valve import throttling_valve

__all__ = [
    "absorption_column",
    "component_splitter",
    "compressor",
    "cooler",
    "distillation_column",
    "ejector",
    "expander",
    "filter",
    "flare",
    "gas_scrubber",
    "heat_exchanger",
    "heater",
    "manifold",
    "mixer",
    "packed_column",
    "pipe",
    "plug_flow_reactor",
    "pump",
    "rate_based_packed_column",
    "separator",
    "shortcut_distillation_column",
    "splitter",
    "stirred_tank_reactor",
    "stripping_column",
    "tank",
    "three_phase_separator",
    "throttling_valve",
]
