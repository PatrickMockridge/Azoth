"""The pure-Python reference implementation for the thermal namespace.

Same contract as :mod:`azoth.hydraulics.reference`: always importable, with or
without the compiled extension, and written to mirror the Rust line for line so a
reviewer can read the two side by side against the published equation.
"""

from __future__ import annotations

from azoth.thermal.reference.boiling_heat_transfer import boiling_heat_transfer
from azoth.thermal.reference.conduction_plane_wall import conduction_plane_wall
from azoth.thermal.reference.fouling_model import fouling_model
from azoth.thermal.reference.lmtd_correction_factor import lmtd_correction_factor
from azoth.thermal.reference.shah_condensation import shah_condensation
from azoth.thermal.reference.shell_and_tube_rating import shell_and_tube_rating
from azoth.thermal.reference.shell_and_tube_vibration import shell_and_tube_vibration
from azoth.thermal.reference.tema_standard import tema_standard
from azoth.thermal.reference.tube_insert_model import tube_insert_model
from azoth.thermal.reference.two_phase_pressure_drop import two_phase_pressure_drop
from azoth.thermal.reference.zone_rating import zone_rating

__all__ = [
    "boiling_heat_transfer",
    "conduction_plane_wall",
    "fouling_model",
    "lmtd_correction_factor",
    "shah_condensation",
    "shell_and_tube_rating",
    "shell_and_tube_vibration",
    "tema_standard",
    "tube_insert_model",
    "two_phase_pressure_drop",
    "zone_rating",
]
