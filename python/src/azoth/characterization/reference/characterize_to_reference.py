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

# The rows are records

``PseudoComponentContribution`` carries twenty-one fields and this id reads four of them. Two are
optional *per row* - a boiling point, without which the row sorts by its molar mass, and a
density, without which it contributes to no volume - and a parallel vector can say "absent" per
vector and not per entry.

# The empty groups are dropped

A group whose accumulated mass or moles is at or below ``1e-12`` produces no profile and the
caller skips it, so a source entirely below the reference's first cut resolves to one component
against a reference of many. ``reference_index`` is what says which cut that one is.

# The class's own options, less the two that need a fluid

``characterize_to_reference`` carries four of ``CharacterizationOptions``' switches:
``inheritReferenceProperties``, ``delumpBeforeRecharacterization`` with its ``delumpResolution``,
and ``sharedImaginaryBoundaries``. The first reads the reference's densities, the last its mole
amounts, and the middle is :mod:`azoth.characterization.reference.delumping`'s arithmetic run over
the source first.

``combineReservoirFluids`` and ``characterizeToCommonSlate`` weight several fluids into one
imaginary pool *before* cutting, and every one of the rest takes or returns a
``SystemInterface`` - so what is left out is what needs a fluid. The delumping was on that list
too and is not, because it is a function of a table.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from itertools import pairwise
from typing import Any

from azoth import _models_gen
from azoth.characterization.reference.delumping import (
    equal_mass_boundaries,
    sorting_key,
    split,
)
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CharacterizeToReferenceResult
from azoth.core.units import from_si, input_records_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.characterize_to_reference"

#: ``PseudoComponentCombiner.MASS_TOLERANCE``.
MASS_TOLERANCE = 1.0e-12

#: ``CharacterizationOptions.Builder``'s own ``delumpResolution``, which an absent input takes.
DELUMP_RESOLUTION = 12


def _rows(source: Sequence[Mapping[str, Any]]) -> list[tuple[float, float, float, float, float]]:
    """The source's rows, ascending by sorting key: `(key, moles, M, rho, mass)`.

    A row with a non-positive mole amount is skipped rather than carried: ``extractComponents``
    drops it, so a zero here is a row that is not on the grid at all. A row that states no density
    is carried at zero, which is what makes it contribute to the group's mass and to neither of
    its volume sums.
    """
    rows = [
        (
            sorting_key(row["boiling_point"], float(row["molar_mass"])),
            float(row["moles"]),
            float(row["molar_mass"]),
            0.0 if row["density"] is None else float(row["density"]),
            float(row["moles"]) * float(row["molar_mass"]),
        )
        for row in source
        if float(row["moles"]) > 0.0
    ]
    rows.sort(key=lambda row: row[0])
    return rows


def characterize_to_reference(
    source: Sequence[Mapping[str, Any]],
    reference: Sequence[Mapping[str, Any]],
    inherit_reference_properties: bool | None = None,
    delump_before_recharacterization: bool | None = None,
    delump_resolution: int | None = None,
    shared_imaginary_boundaries: bool | None = None,
) -> CharacterizeToReferenceResult:
    """A fluid's pseudo-components re-cut onto another fluid's slate.

    ``source`` and ``reference`` are lists of records, one per pseudo-component, keyed by the field
    names the spec declares - each in the unit that field's declaration fixes.

    ``inherit_reference_properties`` is the class's own variant of the same walk: set, each group
    keeps its own moles and takes the **reference row's** molar mass and density rather than the
    ones its own rows average to. Absent or false is the class's default entry point, where the
    reference decides only where the cuts are.

    ``delump_before_recharacterization`` splits each source row into a fine grid before the walk
    bins it, by :func:`azoth.characterization.reference.delumping.delumping`'s arithmetic;
    ``delump_resolution`` sets how fine, and an absent one takes the class's own twelve.
    ``shared_imaginary_boundaries`` places the cuts equally *by mass* on the reference's own fine
    composition instead of at the boiling-point midpoints - the one variant that reads the
    reference's mole amounts.

    **The groups are not the reference rows**: a group that comes out empty is dropped, so
    ``reference_index`` is what re-aligns the answer with the caller's own slate.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    # The one conversion site, shared with the Rust bridge - a record field's unit is in the spec's
    # `fields` block, so the reader is `input_records_to_si` rather than `input_to_si`.
    source_rows: list[dict[str, Any]] = input_records_to_si(spec, "source", source)
    reference_rows: list[dict[str, Any]] = input_records_to_si(spec, "reference", reference)

    apply_checks(checks.on_input, lambda _name: None, warnings)

    inherit = bool(inherit_reference_properties)
    delump = bool(delump_before_recharacterization)
    shared = bool(shared_imaginary_boundaries)
    resolution = DELUMP_RESOLUTION if delump_resolution is None else int(delump_resolution)

    # **The class's own delumping stage**: the source table expanded before the walk bins it. The
    # rows carry no name, so the sub-fractions are named nothing - see the spec's assumption.
    if delump:
        source_rows = [
            {
                "moles": fraction["moles"],
                "molar_mass": fraction["molar_mass"],
                "density": fraction["density"],
                "boiling_point": fraction["boiling_point"],
            }
            for fraction in split(source_rows, resolution)
        ]
    rows = _rows(source_rows)
    # **The rows are sorted, and not only their keys.** The walk reads the keys; the inheriting
    # variant also reads the row that sits at a group's own index - and the shared-boundaries
    # variant reads the same rows for their mole amounts.
    ordered = sorted(
        reference_rows, key=lambda row: sorting_key(row["boiling_point"], float(row["molar_mass"]))
    )
    keys = [sorting_key(row["boiling_point"], float(row["molar_mass"])) for row in ordered]

    if shared:
        # The equal-mass cuts are placed on the masses the reference carries, so this is the one
        # variant that reads its mole amounts - and it refuses a slate that states none rather
        # than placing the cuts on nothing.
        for index, row in enumerate(ordered):
            if row["moles"] is None:
                raise InvalidInputError(
                    f"reference[{index}].moles",
                    "the shared-boundaries variant places the cuts on the reference's own masses, "
                    "and this row states no mole amount",
                )
        boundaries = equal_mass_boundaries(ordered, resolution)
    else:
        # The midpoint of each adjacent pair, or the larger of the two when either is not finite.
        boundaries = [
            0.5 * (low + high) if math.isfinite(low) and math.isfinite(high) else max(low, high)
            for low, high in pairwise(keys)
        ]

    # A forward walk, advancing while the row's key is **greater than** the boundary.
    groups = [{"mass": 0.0, "moles": 0.0, "volume": 0.0, "density_mass": 0.0} for _ in keys]
    index = 0
    boundary = boundaries[0] if boundaries else float("inf")
    for key, moles, _mass, rho, row_mass in rows:
        while index + 1 < len(keys) and key > boundary:
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
        if inherit:
            density = ordered[position]["density"]
            if density is None:
                raise InvalidInputError(
                    f"reference[{position}].density",
                    "the inheriting variant takes each group's density from the reference row, and "
                    "this one states none",
                )
            group_molar_mass.append(float(ordered[position]["molar_mass"]))
            group_density.append(float(density))
        else:
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
