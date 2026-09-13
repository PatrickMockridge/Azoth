"""The component databank, and looking a substance up by name.

`data/components/components.csv` ships with the library, generated from NeqSim
(Equinor/NTNU, Apache-2.0) by `tools/gen_databank.py`; `NOTICE` carries the
attribution. Before this existed every calculation in the `eos` namespace took
`Tc`, `Pc` and `omega` as caller arguments, because there was no licence to ship a
databank. There is now.

    >>> from azoth.eos import component, from_names
    >>> methane = component("methane")
    >>> methane.Tc
    <Quantity(190.56, 'kelvin')>
    >>> fluid = from_names(["methane", "n-butane"])
    >>> len(fluid)
    2

# What is in it, and what is not

173 substances: NeqSim's `COMP.csv` filtered to the types a cubic equation of state
can describe. The ions are excluded because a cubic has no notion of one and NeqSim
fills their critical properties with a shared default - 29 rows carrying the same
`Pc`, `omega` and `Vc` is not a coincidence, and shipping them would ship
plausible-looking wrong numbers. `COMP_EXT.csv` is not vendored at all: 86 MB of
heavy fluids this library cannot characterise.

# Looking up is a call the caller makes

:class:`~azoth.eos.mixture.Component` still carries **no name**, and no calculation
gained a name parameter. A caller who wants a substance by name calls
:func:`component` and passes the result, which keeps the lookup visible in their own
code rather than hidden inside a flash. That was the property `mixture.py` was
protecting when it said a name would tempt something to look a value up; the answer
is not to refuse names, it is to make the lookup an explicit step.

# Provenance

Every row cites the same source. It is deliberately not a per-value citation: the
provenance here is institutional - a project, a version, a file - and a per-row URL
for 173 rows would be a citation-shaped thing that is not a citation.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from functools import cache

from azoth._data import find
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.units import Q, ureg
from azoth.eos.mixture import Component, Mixture

COMPONENTS_CSV = "data/components/components.csv"
KIJ_CSV = "data/components/kij.csv"

#: Columns the loader reads, in order. Named rather than positional because this
#: file's shape is the generator's contract, and a column inserted in the middle
#: should fail loudly rather than shift every value one field left.
COLUMNS = (
    "name",
    "cas",
    "formula",
    "molar_mass_kg_per_mol",
    "tc_k",
    "pc_pa",
    "acentric_factor",
    "critical_volume_m3_per_mol",
    "liquid_density_kg_per_m3",
    "citation",
)


@dataclass(frozen=True, slots=True)
class DatabankEntry:
    """One substance, as the file records it.

    Wider than :class:`~azoth.eos.mixture.Component` on purpose: the file carries
    the CAS number, the formula and the liquid density, which are worth reporting
    and are not part of what a cubic takes.
    """

    name: str
    cas: str
    formula: str
    Tc: Q
    Pc: Q
    omega: float
    molar_mass: Q
    critical_volume: Q
    liquid_density: Q
    citation: str

    def component(self) -> Component:
        """This entry as the thing a calculation takes."""
        return Component(Tc=self.Tc, Pc=self.Pc, omega=self.omega)

    def __repr__(self) -> str:
        return f"DatabankEntry({self.name!r}, Tc={self.Tc}, Pc={self.Pc}, omega={self.omega})"


def _rows(text: str) -> list[dict[str, str]]:
    body = [
        line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")
    ]
    return [dict(row) for row in csv.DictReader(io.StringIO("\n".join(body)))]


@cache
def _table() -> dict[str, DatabankEntry]:
    """Every component, keyed by lower-case name.

    Cached because the file is immutable within a process and every lookup would
    otherwise re-read and re-parse it.
    """
    entries: dict[str, DatabankEntry] = {}
    for row in _rows(find(COMPONENTS_CSV).read_text(encoding="utf-8")):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise InvalidInputError(
                COMPONENTS_CSV, f"a row is missing {missing}; the file is malformed"
            )
        name = row["name"]
        entries[name] = DatabankEntry(
            name=name,
            cas=row["cas"],
            formula=row["formula"],
            Tc=ureg.Quantity(float(row["tc_k"]), "K"),
            Pc=ureg.Quantity(float(row["pc_pa"]), "Pa"),
            omega=float(row["acentric_factor"]),
            molar_mass=ureg.Quantity(float(row["molar_mass_kg_per_mol"]), "kg/mol"),
            critical_volume=ureg.Quantity(float(row["critical_volume_m3_per_mol"]), "m**3/mol"),
            liquid_density=ureg.Quantity(float(row["liquid_density_kg_per_m3"]), "kg/m**3"),
            citation=row["citation"],
        )
    return entries


@cache
def _kij() -> dict[tuple[str, str], float]:
    """Interaction parameters, keyed by the ordered pair of names."""
    pairs: dict[tuple[str, str], float] = {}
    for row in _rows(find(KIJ_CSV).read_text(encoding="utf-8")):
        pair = (row["component_a"], row["component_b"])
        value = float(row["kij_pr"])
        pairs[pair] = value
        pairs[(pair[1], pair[0])] = value
    return pairs


def available() -> tuple[str, ...]:
    """Every name the databank has, sorted."""
    return tuple(sorted(_table()))


def entry(name: str) -> DatabankEntry:
    """One substance's full record.

    Raises:
        PropertyUnavailableError: if the name is not in the databank. An error rather
            than a default: substituting a similar substance would produce a
            plausible answer for the wrong fluid, and the caller would have no way to
            see it.
    """
    key = name.strip().lower()
    try:
        return _table()[key]
    except KeyError:
        raise PropertyUnavailableError(
            name,
            "critical constants",
            f"not in the component databank ({len(_table())} substances). "
            f"`available()` lists them; `component()` takes a Tc, Pc and omega "
            f"directly for anything else.",
        ) from None


def component(name: str) -> Component:
    """One substance as a calculation takes it - `Tc`, `Pc` and `omega`.

    Raises:
        PropertyUnavailableError: as :func:`entry`.
    """
    return entry(name).component()


def kij_for(names: tuple[str, ...]) -> dict[tuple[int, int], float]:
    """The interaction pairs the databank knows, for a list of components.

    Only pairs where both names are present are returned, and only where the
    databank has a value - an unlisted pair is zero, which is the ideal-mixture
    assumption and is what `mixture()` already does with an omitted pair.

    Keyed by index into `names`, which is the form `mixture()` takes.
    """
    pairs: dict[tuple[int, int], float] = {}
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            value = _kij().get((a.strip().lower(), names[j].strip().lower()))
            if value is not None and value != 0.0:
                pairs[(i, j)] = value
    return pairs


def from_names(names: list[str]) -> Mixture:
    """A :class:`~azoth.eos.mixture.Mixture` from a list of databank names.

    The interaction parameters come from the databank too, so a caller writing
    `from_names(["methane", "n-butane"])` gets the published `kij` rather than a
    silent zero. It is a convenience over the databank and never a requirement:
    `Component` and `mixture()` are what every spec case uses, and the library works
    with no data file at all.

    Raises:
        PropertyUnavailableError: if any name is not in the databank.
        InvalidInputError: if the list is empty, or a pair is malformed.
    """
    resolved = [name.strip().lower() for name in names]
    components = tuple(component(name) for name in resolved)
    return mixture(components, kij=kij_for(tuple(resolved)))


# Imported at the bottom because `mixture` lives with the types this module builds
# on, and importing it at the top would make the cycle explicit for no gain.
from azoth.eos.mixture import mixture  # noqa: E402

__all__ = ["DatabankEntry", "available", "component", "entry", "from_names", "kij_for"]
