"""The pure-Python reference implementation for the reactions namespace.

Same contract as :mod:`azoth.thermal.reference`: always importable, with or without
the compiled extension, and written to mirror the Rust line for line so a reviewer can
read the two side by side against the published equation.
"""

from __future__ import annotations

from azoth.reactions.reference.chemical_equilibrium import chemical_equilibrium
from azoth.reactions.reference.equilibrium_constant import equilibrium_constant
from azoth.reactions.reference.reference_potentials import reference_potentials

__all__ = [
    "chemical_equilibrium",
    "equilibrium_constant",
    "reference_potentials",
]
