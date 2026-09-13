"""The pure-Python reference implementation for the thermal namespace.

Same contract as :mod:`azoth.hydraulics.reference`: always importable, with or
without the compiled extension, and written to mirror the Rust line for line so a
reviewer can read the two side by side against the published equation.
"""

from __future__ import annotations

from azoth.thermal.reference.conduction_plane_wall import conduction_plane_wall

__all__ = [
    "conduction_plane_wall",
]
