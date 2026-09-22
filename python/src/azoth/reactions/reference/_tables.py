"""NeqSim's reaction tables, compiled into ``data/reactions/``.

The Python twin of ``crates/azoth-reactions/src/databank.rs``: five tables under
``data/reactions/``, written by ``tools/gen_reaction_data.py`` from the vendored copies
of NeqSim's own ``element``, ``STOCCOEFDATA`` and ``REACTIONDATA`` tables. The two sides
read the same bytes, and ``python/tests/test_data_agreement.py`` checks that.

**The element table is what makes a component reactive.** ``Component.getElements()``
builds an ``Element`` from a row of it, keyed by name, and the element matrix comes from
those rows. It is not parsed from the name: a substance the table does not carry has no
composition, so it cannot enter the balance. The table's 80 components bound what a fluid
can react at all, and MEG, DEG and TEG have no row between them.

**The three sources are different standard states, not a fallback chain.**
``REACTIONDATA``, ``REACTIONDATAPITZER`` and ``REACTIONDATAKENTEISENBERG`` carry the same
reaction names on constants that part company away from 298 K, and the Pitzer table alone
has a ``ValidationStatus`` column.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from functools import cache
from typing import Final, TypedDict

from azoth._data import find
from azoth.core.errors import InvalidInputError, PropertyUnavailableError

ELEMENTS_CSV: Final[str] = "data/reactions/elements.csv"
STOICHIOMETRY_CSV: Final[str] = "data/reactions/stoichiometry.csv"

#: The three sources, by the identifier ``ChemicalReactionDataSource`` gives them.
SOURCES: Final[dict[str, str]] = {
    "standard": "data/reactions/REACTIONDATA.csv",
    "pitzer": "data/reactions/REACTIONDATAPITZER.csv",
    "kent-eisenberg": "data/reactions/REACTIONDATAKENTEISENBERG.csv",
}


class _Row(TypedDict):
    id: str
    index: str
    name: str
    k1: str
    k2: str
    k3: str
    k4: str
    tref: str
    r: str
    actenergy: str
    reference: str
    usereaction: str


@dataclass(frozen=True, slots=True)
class ReactionRow:
    """One reaction's fitted constants, in the units the table holds them."""

    #: The reaction's name, which is what selects it.
    name: str
    #: The four coefficients of ``ln K = K1 + K2/T + K3*ln T + K4*T``.
    coefficients: tuple[float, float, float, float]
    #: The reference temperature the fit is anchored at, in K. **The ``ln K`` expression
    #: does not read it.**
    reference_temperature: float
    #: The pre-exponential rate factor, read by the reference-Arrhenius rate law.
    rate_factor: float
    #: The activation energy, in the table's own units, read by the same law.
    activation_energy: float
    #: The literature citation for the fit, a per-row provenance string.
    reference: str
    #: Whether NeqSim loads the row at all. A zero is why the one combustion row is
    #: dormant.
    use_reaction: bool
    #: The Pitzer source's evidence column. ``None`` for the other two, which do not
    #: carry it rather than carrying an empty one.
    validation_status: str | None


@cache
def _rows(relative: str) -> tuple[dict[str, str], ...]:
    """One compiled table, as dicts keyed by its header."""
    text = find(relative).read_text(encoding="utf-8")
    return tuple(csv.DictReader(io.StringIO(text)))


@cache
def elements() -> tuple[tuple[str, str, float], ...]:
    """Every row of the element table, as ``(component, element, count)``."""
    return tuple(
        (row["componentname"], row["atomelement"], float(row["number"]))
        for row in _rows(ELEMENTS_CSV)
    )


@cache
def element_composition(component: str) -> tuple[tuple[str, float], ...] | None:
    """One component's formula, or ``None`` where the table has no row for it.

    ``None`` is the state that keeps a component out of the element matrix, and it is
    the answer for every glycol: MEG, DEG and TEG have no row.
    """
    found = tuple((element, count) for name, element, count in elements() if name == component)
    return found or None


def source_path(source: str) -> str:
    """The compiled file one source identifier names.

    Raises:
        InvalidInputError: if ``source`` is not one of the three. A refusal rather than
            a fallback to ``standard``, because the three are different standard states
            and a fallback would answer a question the caller did not ask.
    """
    try:
        return SOURCES[source]
    except KeyError:
        raise InvalidInputError(
            "source",
            f"{source!r} is not one of {sorted(SOURCES)}; the three reaction tables are "
            f"different standard states and are not interchangeable",
        ) from None


@cache
def reactions(source: str) -> tuple[ReactionRow, ...]:
    """Every row of one source's table."""
    return tuple(
        ReactionRow(
            name=row["name"],
            coefficients=(
                float(row["k1"]),
                float(row["k2"]),
                float(row["k3"]),
                float(row["k4"]),
            ),
            reference_temperature=float(row["tref"]),
            rate_factor=float(row["r"]),
            activation_energy=float(row["actenergy"]),
            reference=row["reference"],
            use_reaction=float(row["usereaction"]) != 0.0,
            # The Pitzer table's thirteenth column; the other two stop at twelve.
            validation_status=row.get("validationstatus") or None,
        )
        for row in _rows(source_path(source))
    )


def reaction(source: str, name: str) -> ReactionRow:
    """One reaction by name.

    **The first row wins where the name is duplicated**, matching
    ``ChemicalReactionFactory.getChemicalReaction``, whose ``dataSet.next()`` takes the
    first match. ``MDEAprot`` is duplicated in the standard table - an ``Austgen1989``
    row with ``usereaction`` 0 and a ``Huttenhuis2005`` row with 1 - so the disabled row
    is the one that answers.

    Raises:
        PropertyUnavailableError: if the source carries no row by that name. A name the
            data could not answer, which is a different failure from a bad argument.
    """
    for row in reactions(source):
        if row.name == name:
            return row
    raise PropertyUnavailableError(
        fluid=name,
        property_name="reaction",
        reason=f"the {source!r} source carries no row by that name",
    )
