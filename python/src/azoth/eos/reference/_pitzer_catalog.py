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
from azoth.core.errors import InvalidInputError
from azoth.eos.components import pitzer_pair

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
    """The catalogue's lookup key: the canonical species **sorted**, **folded**, joined.

    Sorted, so a row is found whichever order its species are written in - the catalogue
    writes ``Ba+2 Cl-`` and ``Cl- H+``, and both must be found by a caller holding the
    names in either order.

    **Folded, because the catalogue's namespace is NeqSim's and this library's is not.**
    The catalogue writes ``Na+`` and ``HCO3-``, which are ``getComponentName()``'s
    spellings, while the component table writes ``na+`` and ``hco3-``.
    """
    return "|".join(sorted(canonical_species(name).lower() for name in species))


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
    return species.hydrocarbon or is_hydrocarbon_formula(species.formula)


def is_hydrocarbon_formula(formula: str) -> bool:
    """``ComponentGePitzer.hasHydrocarbonFormula``: carbon, hydrogen and digits, both present.

    Public because ``ComponentGePitzer.isHydrocarbon`` is ``super.isHydrocarbon() ||
    hasHydrocarbonFormula()``, and that method decides whether a Henry coefficient is
    **capped to the insoluble limit whatever its row says**: a hydrocarbon whose row carries
    no correlation is not the same state as a solvent whose row carries none.
    """
    if not formula:
        return False
    carbon = False
    hydrogen = False
    for character in formula:
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


#: The identity NeqSim gives the legacy dataset: `PhasePitzer.DEFAULT_PARAMETER_DATASET_ID`.
LEGACY_DATASET_ID = "neqsim-legacy-pitzer-parameters-v1"

#: The identity NeqSim gives the catalogue, carrying the commit the file came from.
PHREEQC_DATASET_ID = "usgs-phreeqc-pitzer-b0b3be767158ccc3322d2c816625cf470045e67e-catalog-v1"

#: Molality above which an ion is in the audited topology.
#:
#: **A different threshold from the selection rule's** ``ACTIVE_MOLES``, and deliberately:
#: `tryApplyCompletePhreeqcPitzerCatalog` tests an absolute mole count while
#: `activeIonIndexes` tests a *molality*, so the audit's topology is the narrower one. A
#: trace ion can be absent from the audit and still have forced the fallback.
ACTIVE_ION_MOLALITY = 1.0e-8

#: The four species the primary-salt audit excludes: `isPrimarySaltCoverageSpecies`.
REACTION_SPECIES = ("H3O+", "OH-", "HCO3-", "CO3--")

_REACTION_SPECIES_LOWER = frozenset(species.lower() for species in REACTION_SPECIES)


def is_primary_salt_species(name: str) -> bool:
    """Whether a component is covered by the primary-salt audit rather than the reaction one."""
    return name.lower() not in _REACTION_SPECIES_LOWER


def dataset_id(selection: tuple[str, str | None]) -> str:
    """The identity of the dataset a selection loads."""
    return PHREEQC_DATASET_ID if selection[0] == "phreeqc" else LEGACY_DATASET_ID


class Coverage(NamedTuple):
    """What the loaded dataset fails to cover, in the shape `PitzerParameterCoverage` reports."""

    dataset_id: str
    active_cations: tuple[str, ...]
    active_anions: tuple[str, ...]
    missing_binary: tuple[str, ...]
    missing_theta: tuple[str, ...]
    missing_psi: tuple[str, ...]

    def is_complete(self) -> bool:
        """Whether every interaction the active topology needs is defined."""
        return not (self.missing_binary or self.missing_theta or self.missing_psi)

    def diagnostic(self) -> str:
        """The diagnostic, byte for byte as `formatDiagnostic` writes it."""
        return (
            f"Pitzer parameter coverage incomplete for dataset '{self.dataset_id}': "
            f"activeCations={_java_list(self.active_cations)}, "
            f"activeAnions={_java_list(self.active_anions)}, "
            f"missingBinary={_java_list(self.missing_binary)}, "
            f"missingTheta={_java_list(self.missing_theta)}, "
            f"missingPsi={_java_list(self.missing_psi)}"
        )


def _java_list(values: Sequence[str]) -> str:
    """`PitzerParameterCoverage.immutableSortedCopy`'s rendering: Java's ``List.toString``."""
    return "[" + ", ".join(values) + "]"


def require_complete(audit: Coverage) -> None:
    """`requireCompletePitzerParameterCoverage`: the audit's refusal.

    **A refusal rather than a zero**: an absent same-sign or ternary parameter is not a
    fitted ideal solution, so evaluating it as one would return a number with no symptom.

    **NeqSim does not always enforce this.** `validateParameterCoverageOncePerState` calls
    it only when :func:`has_mixed_primary_salt_topology` holds, so a single-salt brine with
    an absent binary pair initializes and evaluates the pair at zero. The audit still
    reports it incomplete.

    :raises InvalidInputError: where an interaction the topology needs is absent.
    """
    if not audit.is_complete():
        raise InvalidInputError("composition", audit.diagnostic())


def has_mixed_primary_salt_topology(species: Sequence[Species], solvent_mass: float) -> bool:
    """Whether more than one active primary-salt cation or anion is present.

    The condition under which NeqSim enforces the audit at all: a single cation and a
    single anion keep the "established binary behavior", and a mixed topology does not.
    """
    active_moles = ACTIVE_ION_MOLALITY * solvent_mass
    cations = 0
    anions = 0
    for one in species:
        if (
            one.charge == 0.0
            or not is_primary_salt_species(one.name)
            or not one.moles > active_moles
        ):
            continue
        if one.charge > 0.0:
            cations += 1
        else:
            anions += 1
        if cations > 1 or anions > 1:
            return True
    return False


def coverage(
    species: Sequence[Species], selection: tuple[str, str | None], solvent_mass: float
) -> Coverage:
    """The coverage audit of a phase's primary-salt ions: `getPitzerParameterCoverage`."""
    return _audit(species, selection, solvent_mass, False)


def reaction_coverage(
    species: Sequence[Species], selection: tuple[str, str | None], solvent_mass: float
) -> Coverage:
    """The same audit including the reaction solver's acid-base species."""
    return _audit(species, selection, solvent_mass, True)


def _audit(
    species: Sequence[Species],
    selection: tuple[str, str | None],
    solvent_mass: float,
    include_reaction_species: bool,
) -> Coverage:
    active_moles = ACTIVE_ION_MOLALITY * solvent_mass
    cations = [s for s in species if _is_active(s, True, include_reaction_species, active_moles)]
    anions = [s for s in species if _is_active(s, False, include_reaction_species, active_moles)]

    missing_binary = [
        _report_pair_key(cation.name, anion.name)
        for cation in cations
        for anion in anions
        if not _binary_is_defined(selection, cation.name, anion.name)
    ]

    missing_theta: list[str] = []
    missing_psi: list[str] = []
    _mixed_interactions(cations, anions, selection, missing_theta, missing_psi)
    _mixed_interactions(anions, cations, selection, missing_theta, missing_psi)

    return Coverage(
        dataset_id=dataset_id(selection),
        active_cations=tuple(sorted(s.name for s in cations)),
        active_anions=tuple(sorted(s.name for s in anions)),
        missing_binary=tuple(sorted(missing_binary)),
        missing_theta=tuple(sorted(missing_theta)),
        missing_psi=tuple(sorted(missing_psi)),
    )


def _is_active(
    one: Species, positive: bool, include_reaction_species: bool, active_moles: float
) -> bool:
    wanted = one.charge > 0.0 if positive else one.charge < 0.0
    return (
        wanted
        and (include_reaction_species or is_primary_salt_species(one.name))
        and one.moles > active_moles
    )


def _mixed_interactions(
    same_sign: Sequence[Species],
    opposite: Sequence[Species],
    selection: tuple[str, str | None],
    missing_theta: list[str],
    missing_psi: list[str],
) -> None:
    """The absent theta/psi definitions of every same-sign pair against every opposite ion."""
    for index, one in enumerate(same_sign):
        for other in same_sign[index + 1 :]:
            if not _mixed_is_defined(selection, "THETA", [one.name, other.name]):
                missing_theta.append(_report_pair_key(one.name, other.name))
            for third in opposite:
                if not _mixed_is_defined(selection, "PSI", [one.name, other.name, third.name]):
                    missing_psi.append(_report_psi_key(one.name, other.name, third.name))


def _binary_is_defined(selection: tuple[str, str | None], first: str, second: str) -> bool:
    """Whether the loaded dataset defines a cation-anion pair.

    **Per pair and not per family**, which is NeqSim's own model: the loader calls
    ``setBinaryParameters(i, j, b0, b1, c)`` once for a pair and records one key, so the
    three families are defined together or not at all.
    """
    if selection[0] != "phreeqc":
        return pitzer_pair(first, second) is not None
    return all(find(family, [first, second]) is not None for family in ("B0", "B1", "C0"))


def _mixed_is_defined(selection: tuple[str, str | None], family: str, names: Sequence[str]) -> bool:
    """Whether the loaded dataset defines a same-sign or ternary interaction.

    **False for the legacy dataset however the row reads.** Its loader calls no theta or
    psi setter at all, so ``definedThetaPairs`` and ``definedPsiTuples`` stay empty and
    every mixed topology is incomplete.
    """
    return selection[0] == "phreeqc" and find(family, names) is not None


def _report_pair_key(first: str, second: str) -> str:
    """`PhasePitzer.pairKey`: the two names sorted and joined with ``|``.

    **Not** :func:`species_key`, although the two agree on two names. This one is a *report*
    key and never a lookup - the catalogue's own key sorts all three of a psi row's species
    and this one must not, or the diagnostic would name a triple NeqSim never writes.
    """
    return f"{first}|{second}" if first <= second else f"{second}|{first}"


def _report_psi_key(first: str, second: str, opposite: str) -> str:
    """`PhasePitzer.psiKey`: the sorted same-sign pair, then the opposite ion unsorted."""
    return f"{_report_pair_key(first, second)}|{opposite}"
