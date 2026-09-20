"""The PHREEQC Pitzer catalogue, and the rule that chooses between it and the legacy CSV.

The Python mirror of ``crates/azoth-eos/src/pitzer_catalog.rs``, held to the same
constants by the tests on both sides.

``PhasePitzer.loadParametersFromDatabase`` tries
``PitzerParameterDatasets.tryApplyCompletePhreeqcPitzerCatalog`` **and returns if it
succeeds**, reading ``PitzerParameters.csv`` - the 30 rows this library has always had -
only when it fails. Measured, the choice is per-mixture: ``water + Na+ + Cl-`` takes the
catalogue and ``water + Na+ + HCO3-`` takes the CSV, because the catalogue has no ``C0``
row for that pair. **And the two disagree on pairs they share**: Na+/Cl- is
``0.07534 / 0.2769 / 0.00148`` under PHREEQC against the CSV's
``0.0765 / 0.2664 / 0.00127``.

The coverage rule is pairwise:

============================  ==================================
topology                      required
============================  ==================================
an opposite-sign ion pair     ``B0``, ``B1``, ``C0``
a same-sign ion pair          ``THETA``
three ions, one or two +      ``PSI``
a neutral with a neutral      ``LAMBDA``
a neutral with a cation+anion ``ZETA``
============================  ==================================

``B2`` is the one that is *optional*. And a phase with **no ions** returns false before
the catalogue is consulted at all, so an ion-free mixture always takes the CSV.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Sequence
from functools import cache
from typing import NamedTuple

from azoth._data import find as _find_file

#: The compiled catalogue.
PITZER_PHREEQC_CSV = "data/components/PitzerPhreeqc.csv"

#: How many species each family's rows name, from `PhreeqcPitzerParameterCatalog.Family`.
FAMILY_SPECIES: dict[str, int] = {
    "B0": 2,
    "B1": 2,
    "B2": 2,
    "C0": 2,
    "THETA": 2,
    "PSI": 3,
    "LAMBDA": 2,
    "ZETA": 3,
    "MU": 3,
    "ETA": 3,
    "ALPHAS": 2,
}

#: Moles below which a component is not in the topology.
ACTIVE_MOLES = 1.0e-20

#: Charge below which a component is a neutral rather than an ion.
ACTIVE_CHARGE = 0.5


def canonical_species(name: str) -> str:
    """One species name in NeqSim's spelling: `canonicalSpeciesName`.

    A numeric charge suffix becomes the repeated-sign form and nothing else changes. The
    order of the tests matters: `+3` before `+2`, or `Fe+3` would come back `Fe+` with a
    stray `3`.
    """
    for suffix, sign in (("+3", "+++"), ("-3", "---"), ("+2", "++"), ("-2", "--")):
        if name.endswith(suffix):
            return name[: -len(suffix)] + sign
    return name


def species_key(species: Sequence[str]) -> str:
    """The catalogue's lookup key: the canonical species **sorted** and joined with ``|``.

    Sorted, so a row is found whichever order its species are written in - the catalogue
    writes ``Ba+2 Cl-`` and ``Cl- H+``, and both must be found by a caller holding the
    names in either order.
    """
    return "|".join(sorted(canonical_species(name) for name in species))


@cache
def _rows() -> dict[tuple[str, str], tuple[float, ...]]:
    """The compiled catalogue, by ``(family, species_key)``.

    Read by the file's own header rather than by position, so a column inserted in the
    middle resolves to the wrong *name* rather than shifting every value one field left.
    """
    text = _find_file(PITZER_PHREEQC_CSV).read_text(encoding="utf-8")
    body = [line for line in text.splitlines() if line.strip()]
    out: dict[tuple[str, str], tuple[float, ...]] = {}
    for row in csv.DictReader(io.StringIO("\n".join(body))):
        out[(row["family"], row["species_key"])] = tuple(float(row[f"a{i}"]) for i in range(6))
    return out


def family_rows(family: str) -> tuple[tuple[str, tuple[float, ...]], ...]:
    """Every row of one family, sorted by key, so this is a function of the data."""
    rows = [(key, a) for (row_family, key), a in _rows().items() if row_family == family]
    return tuple(sorted(rows))


def find(family: str, species: Sequence[str]) -> tuple[float, ...] | None:
    """One row's six coefficients, or ``None`` where the catalogue carries no such row.

    ``find`` and not ``require``: the caller decides whether an absent row is a fallback,
    an optional term or an error, which is exactly the distinction ``B2`` and the coverage
    rule are built on.
    """
    return _rows().get((family, species_key(species)))


class Species(NamedTuple):
    """One component as the selection rule sees it."""

    name: str
    moles: float
    charge: float
    formula: str
    #: NeqSim's `componentType == "HC"` clause of `isHydrocarbon`.
    hydrocarbon: bool


def is_hydrocarbon(species: Species) -> bool:
    """Whether a neutral component is a hydrocarbon the automatic catalogue excludes.

    `isHydrocarbonForAutomaticCatalog`: NeqSim's own flag, or **a formula of carbon,
    hydrogen and digits with both present**. The formula clause is what the method exists
    for - a database component such as methane keeps the type ``normal`` in a GE phase, so
    the type alone would let it into a topology the catalogue has no rows for.

    **Two of NeqSim's three flag clauses have no counterpart here.** ``isIsTBPfraction``
    and ``isPlusFraction`` are set through the API rather than read from the databank, and
    the components they mark are typed ``HC`` in the table anyway.
    """
    if species.hydrocarbon:
        return True
    if not species.formula:
        return False
    carbon = False
    hydrogen = False
    for character in species.formula:
        if character == "C":
            carbon = True
        elif character == "H":
            hydrogen = True
        elif not character.isdigit():
            return False
    return carbon and hydrogen


def select_dataset(species: Sequence[Species]) -> tuple[str, str | None]:
    """Which dataset covers a phase's topology, and why not when the answer is the CSV.

    Returns ``("phreeqc", None)`` or ``("legacy", reason)``. A pure function of the
    composition, which is the point: NeqSim decides this inside
    ``loadParametersFromDatabase`` against a live ``PhasePitzer``, and the rule is about
    the *topology* rather than the phase's state, so it can be answered - and tested -
    without one.
    """
    active = [s for s in species if s.moles > ACTIVE_MOLES]
    ions = [s for s in active if abs(s.charge) >= ACTIVE_CHARGE]
    if not ions:
        return ("legacy", "NoIons")

    neutrals = [
        s
        for s in active
        if abs(s.charge) < ACTIVE_CHARGE and s.name.lower() != "water" and not is_hydrocarbon(s)
    ]

    def missing(family: str, names: Sequence[str]) -> str | None:
        if find(family, names) is None:
            return f"{family} for {', '.join(names)}"
        return None

    for index, one in enumerate(ions):
        for other in ions[index + 1 :]:
            names = [one.name, other.name]
            if one.charge * other.charge < 0.0:
                for family in ("B0", "B1", "C0"):
                    if (reason := missing(family, names)) is not None:
                        return ("legacy", reason)
            elif (reason := missing("THETA", names)) is not None:
                return ("legacy", reason)

    # A triple with one or two positives, which is the shape PSI describes.
    for first, one in enumerate(ions):
        for second in range(first + 1, len(ions)):
            for third in ions[second + 1 :]:
                positives = sum(1 for s in (one, ions[second], third) if s.charge > 0.0)
                if positives in (1, 2):
                    names = [one.name, ions[second].name, third.name]
                    if (reason := missing("PSI", names)) is not None:
                        return ("legacy", reason)

    for position, neutral in enumerate(neutrals):
        for other in neutrals[position:]:
            if (reason := missing("LAMBDA", [neutral.name, other.name])) is not None:
                return ("legacy", reason)
        for ion in ions:
            if (reason := missing("LAMBDA", [neutral.name, ion.name])) is not None:
                return ("legacy", reason)
        for cation in (s for s in ions if s.charge > 0.0):
            for anion in (s for s in ions if s.charge < 0.0):
                names = [neutral.name, cation.name, anion.name]
                if (reason := missing("ZETA", names)) is not None:
                    return ("legacy", reason)

    return ("phreeqc", None)
