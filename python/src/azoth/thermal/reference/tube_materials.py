"""The shell-and-tube material table, vendored from NeqSim's design database.

``ShellAndTubeDesignCalculator.loadMaterialFromDB`` reads the ``HeatExchangerTubeMaterials``
table by exact grade; this is that table, compiled to
``data/heat_exchanger/tube_materials.csv``. **``nace_compliant`` stays the file's text**: the
class tests it against ``Yes``, which no row holds, and the port reproduces that defect.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from functools import cache

from azoth._data import find

TABLE_PATH = "data/heat_exchanger/tube_materials.csv"


@dataclass(frozen=True, slots=True)
class TubeMaterial:
    """One material grade, as the class reads it."""

    material_grade: str
    specification: str
    material_type: str
    smys_mpa: float
    allowable_stress_mpa: float
    thermal_conductivity_wpm_k: float
    nace_compliant: str

    def is_nace_compliant(self) -> bool:
        """The class's own test, which no row satisfies."""
        return self.nace_compliant.lower() == "yes"


@cache
def table() -> dict[str, TubeMaterial]:
    """The compiled table by exact grade."""
    text = find(TABLE_PATH).read_text(encoding="utf-8")
    rows: dict[str, TubeMaterial] = {}
    for row in csv.DictReader(io.StringIO(text)):
        rows[row["material_grade"]] = TubeMaterial(
            material_grade=row["material_grade"],
            specification=row["specification"],
            material_type=row["material_type"],
            smys_mpa=float(row["smys_mpa"]),
            allowable_stress_mpa=float(row["allowable_stress_mpa"]),
            thermal_conductivity_wpm_k=float(row["thermal_conductivity_wpm_k"]),
            nace_compliant=row["nace_compliant"],
        )
    return rows


def lookup(grade: str) -> TubeMaterial | None:
    """A grade by exact equality, or ``None``."""
    return table().get(grade)


__all__ = ["TABLE_PATH", "TubeMaterial", "lookup", "table"]
