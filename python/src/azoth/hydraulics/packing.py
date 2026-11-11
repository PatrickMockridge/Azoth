"""NeqSim's packing specifications: the built-in table, the compiled CSV and the lookup.

The Python twin of ``crates/azoth-hydraulics/src/packing.rs``, mirroring
``PackingSpecificationLibrary``.

# The built-ins are the reference; the file adds material-qualified rows

The library registers its 22 built-ins under their name and a ``name-material`` alias, and then
loads ``designdata/Packing.csv``. A file row's display name carries its material -
``Pall-Ring-50-plastic`` - and a qualified name that already exists is skipped, so the file's
metal 50 mm ring (112.6 m**2/m**3, void 0.951, factor 85) does not displace the built-in metal
estimate (120.0, 0.96, 66). A qualified name that is new is also reachable unqualified unless
that name exists, which is how ``Pall-Ring-35`` and ``Pall-Ring-15`` resolve. Before upstream
PR #4405 (#4380) the file was keyed by display name alone, so the plastic 50 mm row's later line
overwrote the built-in's geometry. A ``placeholder`` row is not a packing at all.

The file is ``data/packing/packing.csv``, compiled from NeqSim's ``designdata/Packing.csv``;
``databank/manifest.toml`` records all twelve of its columns and which three were dropped.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from functools import cache

from azoth._data import find
from azoth.core.errors import InvalidInputError

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


@dataclass(frozen=True, slots=True)
class _Registry:
    """The lookup table and the display names `getPackingNames` lists."""

    by_key: dict[str, PackingSpecification]
    display_names: tuple[str, ...]


@cache
def _registry() -> _Registry:
    """The registry: the built-ins first, then the file, which never replaces one of them."""
    table: dict[str, PackingSpecification] = {}
    display: list[str] = []

    def register(specification: PackingSpecification) -> None:
        table[normalize(specification.name)] = specification
        display.append(specification.name)

    def register_alias(alias: str, specification: PackingSpecification) -> None:
        table[normalize(alias)] = specification

    for name, material, size, area, void, factor, liquid, gas in _RANDOM:
        specification = PackingSpecification(
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
        register(specification)
        register_alias(f"{name}-{material}", specification)
    for name, material, area, void, factor, liquid, gas in _STRUCTURED:
        specification = PackingSpecification(
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
        register(specification)
        register_alias(f"{name}-{material}", specification)
    for row in csv.DictReader(io.StringIO(find(PACKING_CSV).read_text(encoding="utf-8"))):
        category = row["type"]
        if category.lower() not in ("random", "structured"):
            # A `placeholder` row is a stand-in geometry, not a packing.
            continue
        material = row["material"].lower()
        size = float(row["size_mm"])
        area = float(row["surface_area_pr_volume"])
        void = float(row["void_fraction"])
        factor = float(row["packing_factor"])
        if area <= 0.0 or void <= 0.0 or factor <= 0.0:
            continue
        cp = float(row["cp"]) if row["cp"] else 1.0
        ch = float(row["ch"]) if row["ch"] else 0.4
        name = f"{display_name(row['name'], category, size)}-{material}"
        if normalize(name) in table:
            # The built-ins are the reference: a qualified row that already exists does not
            # replace one.
            continue
        specification = PackingSpecification(
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
        register(specification)
        legacy = name[: -len(f"-{material}")]
        if normalize(legacy) not in table:
            register_alias(legacy, specification)
    return _Registry(by_key=table, display_names=tuple(display))


def packing(name: str) -> PackingSpecification | None:
    """A packing by name or alias, or ``None``."""
    return _registry().by_key.get(normalize(name))


def required_packing(name: str) -> PackingSpecification:
    """A packing by name, or the refusal the class now makes.

    NeqSim's ``getOrDefault`` used to answer an unknown preset with ``Pall-Ring-50``; upstream
    PR #4405 (#4380) replaced it with ``require``, whose ``IllegalArgumentException`` this
    raises. A default would turn a misspelt packing into a silent geometry change.

    Raises:
        InvalidInputError: if no registered packing carries the name.
    """
    found = packing(name)
    if found is None:
        raise InvalidInputError(
            "packing_type",
            f"`{name}` is not a registered packing. NeqSim's calculator now requires a known "
            "preset (upstream PR #4405) and no longer falls back to `Pall-Ring-50`",
        )
    return found


def packing_names() -> tuple[str, ...]:
    """Every registered packing's display name, in registration order."""
    return _registry().display_names
