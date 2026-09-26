"""NeqSim's packing specifications: the built-in table, the compiled CSV and the lookup.

The Python twin of ``crates/azoth-hydraulics/src/packing.rs``, mirroring
``PackingSpecificationLibrary``.

# The file wins a collision

The library registers its 22 built-ins and **then** loads ``designdata/Packing.csv``, so a row
whose normalized name matches an earlier registration replaces it. Measured, ``Pall-Ring-50``
resolves to the file's *plastic* row - 111.1 m**2/m**3, void 0.919, factor 180 - and not to the
built-in's metal one (120.0, 0.96, 66), because both normalize to ``pallring50``.

The file is ``data/packing/packing.csv``, compiled from NeqSim's ``designdata/Packing.csv``;
``databank/manifest.toml`` records all twelve of its columns and which three were dropped.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from functools import cache

from azoth._data import find

#: The compiled table, which is the manifest's `neqsim/Packing.csv` entry.
PACKING_CSV = "data/packing/packing.csv"


@dataclass(frozen=True, slots=True, eq=False)
class PackingSpecification:
    """A packing's geometry and material constants.

    Every one of these is positive by the class's own validation: a non-positive area, void
    fraction, packing factor, critical surface tension or Billet constant is refused rather than
    carried.
    """

    #: The display name, which for a random packing the library derives from the raw name and the
    #: nominal size - `Pall-Ring-50`, not `pallring`.
    name: str
    #: `random` or `structured`.
    category: str
    #: The material, which the critical surface tension is read from.
    material: str
    #: The nominal size in millimetres; zero for a structured packing.
    nominal_size_mm: float
    #: The specific surface area, in m**2/m**3.
    specific_surface_area: float
    #: The void fraction.
    void_fraction: float
    #: The packing factor, in 1/m.
    packing_factor: float
    #: The critical surface tension of the material, in N/m.
    critical_surface_tension: float
    #: The Billet-Schultes liquid constant.
    billet_liquid_constant: float
    #: The Billet-Schultes gas constant.
    billet_gas_constant: float


def normalize(name: str) -> str:
    """NeqSim's own alias normalisation: lowercased with every non-alphanumeric removed.

    Example:
        >>> normalize("PALL RING 50") == normalize("Pall-Ring-50")
        True
    """
    return "".join(character for character in name.lower() if character.isalnum())


def critical_surface_tension(material: str) -> float:
    """The critical surface tension a material's packing takes, in N/m.

    A plastic packing 0.033, a ceramic one 0.061, and anything else 0.075.
    """
    normalized = material.lower()
    if "plastic" in normalized:
        return 0.033
    if "ceramic" in normalized:
        return 0.061
    return 0.075


def display_name(raw_name: str, category: str, size_mm: float) -> str:
    """The display name a CSV row's raw name, category and size make.

    A structured packing keeps the raw name; a random one is named after its family, with the
    rounded size appended.
    """
    if category.lower() == "structured":
        return raw_name
    normalized = normalize(raw_name)
    rounded = round(size_mm)
    if "pallring" in normalized:
        return f"Pall-Ring-{rounded}"
    if any(family in normalized for family in ("rashig", "raschig", "rachig")):
        return f"Raschig-Ring-{rounded}"
    return f"{raw_name}-{rounded}"


#: `registerBuiltIns`' thirteen random packings: name, material, size in mm, area, void fraction,
#: packing factor, and the two Billet constants.
_RANDOM = (
    ("Pall-Ring-25", "metal", 25.0, 210.0, 0.94, 157.0, 1.0, 0.40),
    ("Pall-Ring-38", "metal", 38.0, 164.0, 0.95, 92.0, 1.0, 0.40),
    ("Pall-Ring-50", "metal", 50.0, 120.0, 0.96, 66.0, 1.0, 0.40),
    ("Raschig-Ring-25", "ceramic", 25.0, 190.0, 0.68, 580.0, 1.0, 0.40),
    ("Raschig-Ring-50", "ceramic", 50.0, 95.0, 0.74, 155.0, 1.0, 0.40),
    ("IMTP-25", "metal", 25.0, 226.0, 0.97, 134.0, 1.05, 0.42),
    ("IMTP-40", "metal", 40.0, 151.0, 0.97, 79.0, 1.05, 0.42),
    ("IMTP-50", "metal", 50.0, 102.0, 0.98, 56.0, 1.05, 0.42),
    ("IMTP-70", "metal", 70.0, 72.0, 0.98, 36.0, 1.05, 0.42),
    ("Berl-Saddle-25", "ceramic", 25.0, 260.0, 0.68, 360.0, 0.95, 0.38),
    ("Berl-Saddle-38", "ceramic", 38.0, 165.0, 0.70, 220.0, 0.95, 0.38),
    ("Berl-Saddle-50", "ceramic", 50.0, 105.0, 0.72, 150.0, 0.95, 0.38),
    ("Intalox-Saddle-25", "ceramic", 25.0, 255.0, 0.78, 200.0, 1.0, 0.40),
)

#: `registerBuiltIns`' nine structured packings, whose nominal size is zero.
_STRUCTURED = (
    ("Mellapak-125Y", "metal", 125.0, 0.99, 33.0, 1.2, 0.45),
    ("Mellapak-250Y", "metal", 250.0, 0.98, 66.0, 1.2, 0.45),
    ("Mellapak-350Y", "metal", 350.0, 0.97, 105.0, 1.2, 0.45),
    ("Mellapak-500Y", "metal", 500.0, 0.96, 180.0, 1.2, 0.45),
    ("Flexipac-1Y", "metal", 135.0, 0.99, 36.0, 1.15, 0.44),
    ("Flexipac-2Y", "metal", 220.0, 0.98, 60.0, 1.15, 0.44),
    ("Flexipac-3Y", "metal", 340.0, 0.97, 100.0, 1.15, 0.44),
    ("Sulzer-BX", "metal", 500.0, 0.90, 140.0, 1.3, 0.50),
    ("Sulzer-CY", "metal", 750.0, 0.88, 220.0, 1.3, 0.50),
)


@cache
def _table() -> dict[str, PackingSpecification]:
    """The registry: the built-ins first, then the file, which is what lets a file row replace."""
    table: dict[str, PackingSpecification] = {}
    for name, material, size, area, void, factor, liquid, gas in _RANDOM:
        table[normalize(name)] = PackingSpecification(
            name=name,
            category="random",
            material=material,
            nominal_size_mm=size,
            specific_surface_area=area,
            void_fraction=void,
            packing_factor=factor,
            critical_surface_tension=critical_surface_tension(material),
            billet_liquid_constant=liquid,
            billet_gas_constant=gas,
        )
    for name, material, area, void, factor, liquid, gas in _STRUCTURED:
        table[normalize(name)] = PackingSpecification(
            name=name,
            category="structured",
            material=material,
            nominal_size_mm=0.0,
            specific_surface_area=area,
            void_fraction=void,
            packing_factor=factor,
            critical_surface_tension=critical_surface_tension(material),
            billet_liquid_constant=liquid,
            billet_gas_constant=gas,
        )
    for row in csv.DictReader(io.StringIO(find(PACKING_CSV).read_text(encoding="utf-8"))):
        size = float(row["size_mm"])
        area = float(row["surface_area_pr_volume"])
        void = float(row["void_fraction"])
        factor = float(row["packing_factor"])
        # A row whose geometry is not positive is skipped, which is `parseCsvLine`'s own check.
        if area <= 0.0 or void <= 0.0 or factor <= 0.0:
            continue
        category = row["type"]
        material = row["material"]
        cp = float(row["cp"])
        ch = float(row["ch"])
        name = display_name(row["name"], category, size)
        table[normalize(name)] = PackingSpecification(
            name=name,
            category=category,
            material=material,
            nominal_size_mm=size,
            specific_surface_area=area,
            void_fraction=void,
            packing_factor=factor,
            critical_surface_tension=critical_surface_tension(material),
            billet_liquid_constant=cp if cp > 0.0 else 1.0,
            billet_gas_constant=ch / 6.0 if ch > 0.0 else 0.4,
        )
    return table


def packing(name: str) -> PackingSpecification | None:
    """A packing by name or alias, or ``None``."""
    return _table().get(normalize(name))


def packing_or_default(name: str) -> PackingSpecification:
    """A packing by name, falling back to ``Pall-Ring-50``, which is the class's own default.

    Example:
        >>> packing_or_default("nothing-like-this").name
        'Pall-Ring-50'
    """
    found = packing(name) or packing("Pall-Ring-50")
    assert found is not None, "Pall-Ring-50 is a built-in packing"
    return found


def packing_names() -> tuple[str, ...]:
    """Every registered packing's name, sorted."""
    return tuple(sorted(specification.name for specification in _table().values()))
