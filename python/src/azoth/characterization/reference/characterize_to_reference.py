"""``characterization.characterize_to_reference`` - one fluid re-cut onto another's slate.

Spec: ``specs/models/characterization/characterize_to_reference.toml``. Oracle:
``validation/neqsim/captures/pseudo_component_combiner_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/characterize_to_reference.rs`` line for line.

# Read out of the middle of ``PseudoComponentCombiner``

The class is 1841 lines and every public method takes or returns a ``SystemInterface``. This is
the part of it that is a function of two tables: ``extractComponents`` lifts each fluid's
pseudo-components into rows keyed by a boiling point that falls back to the molar mass,
``determineReferenceBoundaries`` puts a cut at the midpoint of each adjacent pair of the
*reference's* keys, and ``distributeToProfiles`` walks the source binning each row into the group
its key falls in.

A group's molar mass is its accumulated mass over its accumulated moles; its density is that mass
over the volume the rows would occupy, ``sum(m) / sum(m/rho)``.

# The empty groups are dropped

A group whose accumulated mass or moles is at or below ``1e-12`` produces no profile and the
caller skips it, so a source entirely below the reference's first cut resolves to one component
against a reference of many. ``reference_index`` is what says which cut that one is.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CharacterizeToReferenceResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.characterize_to_reference"

#: ``PseudoComponentCombiner.MASS_TOLERANCE``.
MASS_TOLERANCE = 1.0e-12


def sorting_key(boiling_point: float, molar_mass: float) -> float:
    """``PseudoComponentContribution.sortingKey``: the boiling point, or the molar mass."""
    return boiling_point if boiling_point > 0.0 and boiling_point != float("inf") else molar_mass


def _rows(
    moles: list[float],
    molar_mass: list[float],
    density: list[float],
    boiling_point: list[float],
    field: str,
) -> list[tuple[float, float, float, float, float]]:
    """A fluid's rows, ascending by sorting key: `(key, moles, M, rho, mass)`."""
    count = len(moles)
    if len(molar_mass) != count or len(density) != count or len(boiling_point) != count:
        raise InvalidInputError(
            field,
            f"a fluid's table is four aligned vectors: {count} mole amounts against "
            f"{len(molar_mass)} molar masses, {len(density)} densities and "
            f"{len(boiling_point)} boiling points",
        )
    rows = [
        (sorting_key(tb, mass), amount, mass, rho, amount * mass)
        for amount, mass, rho, tb in zip(moles, molar_mass, density, boiling_point, strict=True)
        if amount > 0.0
    ]
    rows.sort(key=lambda row: row[0])
    return rows


def characterize_to_reference(
    source_moles: list[float],
    source_molar_mass: list[Q],
    source_density: list[Q],
    source_boiling_point: list[float],
    reference_molar_mass: list[Q],
    reference_boiling_point: list[float],
) -> CharacterizeToReferenceResult:
    """A fluid's pseudo-components re-cut onto another fluid's slate.

    ``reference_*`` decides where the cuts are and nothing else - its densities and amounts are
    not read. **The groups are not the reference rows**: a group that comes out empty is dropped,
    so ``reference_index`` is what re-aligns the answer with the caller's own slate.

    Raises:
        InvalidInputError: if any table is not internally aligned.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    apply_checks(checks.on_input, lambda _name: None, warnings)

    source = _rows(
        [input_to_si(spec, "source_moles", value) for value in source_moles],
        [input_to_si(spec, "source_molar_mass", value) for value in source_molar_mass],
        [input_to_si(spec, "source_density", value) for value in source_density],
        [input_to_si(spec, "source_boiling_point", value) for value in source_boiling_point],
        "source_moles",
    )
    if len(reference_molar_mass) != len(reference_boiling_point):
        raise InvalidInputError(
            "reference_molar_mass",
            f"the reference is two aligned vectors: {len(reference_molar_mass)} molar masses "
            f"against {len(reference_boiling_point)} boiling points",
        )
    reference = sorted(
        sorting_key(input_to_si(spec, "reference_boiling_point", tb), mass)
        for mass, tb in zip(
            [input_to_si(spec, "reference_molar_mass", value) for value in reference_molar_mass],
            reference_boiling_point,
            strict=True,
        )
    )

    # The midpoint of each adjacent pair, or the larger of the two when either is not finite.
    boundaries = [
        0.5 * (low + high) if low != float("inf") and high != float("inf") else max(low, high)
        for low, high in zip(reference, reference[1:], strict=False)
    ]

    # A forward walk, advancing while the row's key is **greater than** the boundary.
    groups = [{"mass": 0.0, "moles": 0.0, "volume": 0.0, "density_mass": 0.0} for _ in reference]
    index = 0
    boundary = boundaries[0] if boundaries else float("inf")
    for key, moles, _mass, rho, row_mass in source:
        while index + 1 < len(reference) and key > boundary:
            index += 1
            boundary = boundaries[index] if index < len(boundaries) else float("inf")
        group = groups[index]
        group["mass"] += row_mass
        group["moles"] += moles
        if rho > 0.0:
            group["volume"] += row_mass / rho
            group["density_mass"] += row_mass * rho

    reference_index: list[float] = []
    group_moles: list[float] = []
    group_molar_mass: list[float] = []
    group_density: list[float] = []
    for position, group in enumerate(groups):
        if not (group["mass"] > MASS_TOLERANCE) or not (group["moles"] > MASS_TOLERANCE):
            continue
        reference_index.append(float(position))
        group_moles.append(group["moles"])
        group_molar_mass.append(group["mass"] / group["moles"])
        group_density.append(
            group["mass"] / group["volume"]
            if group["volume"] > MASS_TOLERANCE
            else group["density_mass"] / group["mass"]
        )

    apply_checks(checks.derived, lambda _name: None, warnings)

    return CharacterizeToReferenceResult(
        reference_index=tuple(reference_index),
        group_moles=tuple(from_si(value, "mol") for value in group_moles),
        group_molar_mass=tuple(from_si(value, "kg/mol") for value in group_molar_mass),
        group_density=tuple(from_si(value, "kg/m**3") for value in group_density),
        warnings=tuple(warnings),
    )
