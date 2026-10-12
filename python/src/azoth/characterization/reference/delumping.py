"""``characterization.delumping`` - a lumped pseudo-component table split into a finer grid.

Spec: ``specs/models/characterization/delumping.toml``. Oracle:
``validation/neqsim/captures/pseudo_component_delumping_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/delumping.rs`` line for line.

# Two conservations, and they are different ones

``delumpContributions`` splits each lump into ``resolution`` sub-fractions whose mole *weights*
sum to one, so the parent's moles are conserved by the weights themselves - and then rescales the
molar-mass grid by a single factor so that ``sum(n_k * M_k)`` is the parent's mass exactly. The
moles are what the light-end bias was stated in, so it is the molar masses that move.

The sub-fractions' normal boiling points are spread with those molar masses,
``Tb_k = Tb * M_k / M``, which is what lets them cross a reference slate's cut boundaries - and
the expanded list is re-sorted by sorting key so it can be fed straight to a cut.

# One read of five fields, and sixteen copies

The class's row carries twenty-one numbers. This stage *reads* five of them and *copies* the rest
to every sub-fraction unchanged, which is why the row here is five fields wide: the helper is
shared with the walk that redistributes the result, and that walk is where the copies would be
read.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from itertools import pairwise
from typing import Any

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import DelumpingResult, SubFraction
from azoth.core.units import from_si, input_records_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.delumping"

#: ``PseudoComponentCombiner.MASS_TOLERANCE``, on the rescale's denominator.
MASS_TOLERANCE = 1.0e-12

#: The class's own ``window`` and ``decay``: the grid's half-width, and how fast the bias falls.
WINDOW = 0.4
DECAY = 1.5


def sorting_key(boiling_point: float | None, molar_mass: float) -> float:
    """``PseudoComponentContribution.sortingKey``: the boiling point, or the molar mass.

    **The absence has two spellings and they are the same thing.** A row that states no boiling
    point, and a row that states one which is not finite and positive, both fall back - the first
    because the record has no field for it and the second because the class's own test is on the
    value. A sentinel a caller wrote into a vector is now a field the caller leaves out.

    Shared with ``characterize_to_reference``, which sorts the same table by the same key.
    """
    if boiling_point is None or not (math.isfinite(boiling_point) and boiling_point > 0.0):
        return molar_mass
    return boiling_point


def _midpoint(low: float, high: float) -> float:
    """The midpoint of two keys, or the larger of them where either is not finite."""
    if math.isfinite(low) and math.isfinite(high):
        return 0.5 * (low + high)
    return max(low, high)


def _quantile_boundaries(rows: Sequence[Mapping[str, Any]], count: int) -> list[float]:
    """``determineQuantileBoundaries``: the sorting key of the row each mass target falls in.

    The targets are equal **masses**, ``total * (i + 1) / count``, and the boundary is the key of
    the sub-fraction whose running mass first reaches one - not an interpolation, so a boundary is
    always a key some row carries. A target nothing reaches is padded with an infinity, which the
    clamp below then discards because it is not finite.
    """
    if count <= 1 or not rows:
        return []
    total = sum(float(row["moles"]) * float(row["molar_mass"]) for row in rows)
    if not (total > MASS_TOLERANCE):
        return []
    targets = [total * (i + 1) / count for i in range(count - 1)]

    boundaries: list[float] = []
    cumulative = 0.0
    index = 0
    for row in rows:
        following = cumulative + float(row["moles"]) * float(row["molar_mass"])
        while index < len(targets) and following >= targets[index] - MASS_TOLERANCE:
            boundaries.append(sorting_key(row["boiling_point"], float(row["molar_mass"])))
            index += 1
        cumulative = following
    while len(boundaries) < len(targets):
        boundaries.append(float("inf"))
    return boundaries


def equal_mass_boundaries(reference: Sequence[Mapping[str, Any]], resolution: int) -> list[float]:
    """``determineReferenceEqualMassBoundaries``: the cuts placed equally by mass.

    The slate is delumped into its imaginary fine composition, equal-mass quantiles are placed on
    it, and each cut is then **clamped strictly inside the gap between the two rows it separates** -
    falling back to the midpoint where no quantile falls in the gap. The clamp is what keeps every
    reference row inside its own cut, which the property inheritance needs, so the two rules can
    only ever differ for a cut that lands between them.

    A ``resolution`` of one or less has nothing to delump and gives back the midpoints, which is the
    class's own fall-through.
    """
    keys = [sorting_key(row["boiling_point"], float(row["molar_mass"])) for row in reference]
    midpoints = [_midpoint(low, high) for low, high in pairwise(keys)]
    if resolution <= 1:
        return midpoints
    equal_mass = _quantile_boundaries(split(reference, resolution), len(reference))

    out: list[float] = []
    for index, midpoint in enumerate(midpoints):
        lower, upper = keys[index], keys[index + 1]
        if not (math.isfinite(lower) and math.isfinite(upper) and upper > lower):
            out.append(midpoint)
            continue
        candidate = equal_mass[index] if index < len(equal_mass) else float("nan")
        out.append(
            candidate if math.isfinite(candidate) and lower < candidate < upper else midpoint
        )
    return out


def _spread(boiling_point: float | None, ratio: float) -> float | None:
    """``parentHasTb ? Tb * ratio : Tb``: the sub-fraction's boiling point, spread with its mass.

    A parent that states no boiling point, or one that is not finite and positive, passes its own
    value through - which is the class's own test, and the same two spellings ``sorting_key``
    treats alike.
    """
    if boiling_point is None or not (math.isfinite(boiling_point) and boiling_point > 0.0):
        return boiling_point
    return boiling_point * ratio


def _unsplit(lump: Mapping[str, Any]) -> dict[str, Any]:
    """One lump, as the sub-fraction that is the lump itself - what an unsplit one becomes.

    **A row may carry no name**, which is what the walk's own source table is: it takes the split's
    arithmetic and drops the names, so a row without one is a sub-fraction that is not named.
    """
    return {
        "name": str(lump.get("name", "")),
        "moles": float(lump["moles"]),
        "molar_mass": float(lump["molar_mass"]),
        "density": None if lump["density"] is None else float(lump["density"]),
        "boiling_point": None if lump["boiling_point"] is None else float(lump["boiling_point"]),
    }


def split(lumps: Sequence[Mapping[str, Any]], resolution: int) -> list[dict[str, Any]]:
    """``delumpContributions``, without the sort: the expanded table in the input's own order.

    **A ``resolution`` of one or less returns the input unchanged, order included** - which is the
    class's early return, and the reason this is a separate function: the sort at the bottom is
    part of the splitting and not of the identity.
    """
    if resolution <= 1:
        return [_unsplit(lump) for lump in lumps]

    expanded: list[dict[str, Any]] = []
    for lump in lumps:
        moles = float(lump["moles"])
        molar_mass = float(lump["molar_mass"])
        # The class's own guard: a lump it cannot split is passed through rather than refused,
        # because the method defines what it does with one.
        if not (moles > 0.0) or not (molar_mass > 0.0):
            expanded.append(_unsplit(lump))
            continue

        position = [k / (resolution - 1) for k in range(resolution)]
        weights = [math.exp(-DECAY * t) for t in position]
        weight_sum = sum(weights)
        grid = [molar_mass * (1.0 - WINDOW + 2.0 * WINDOW * t) for t in position]

        raw_mass = sum(moles * (weights[k] / weight_sum) * grid[k] for k in range(resolution))
        parent_mass = moles * molar_mass
        scale = parent_mass / raw_mass if raw_mass > MASS_TOLERANCE else 1.0

        boiling_point = None if lump["boiling_point"] is None else float(lump["boiling_point"])
        for k in range(resolution):
            sub_moles = moles * (weights[k] / weight_sum)
            sub_mass = grid[k] * scale
            ratio = sub_mass / molar_mass
            expanded.append(
                {
                    "name": f"{lump.get('name', '')}_d{k}",
                    "moles": sub_moles,
                    "molar_mass": sub_mass,
                    "density": None if lump["density"] is None else float(lump["density"]),
                    "boiling_point": _spread(boiling_point, ratio),
                }
            )
    expanded.sort(key=lambda row: sorting_key(row["boiling_point"], row["molar_mass"]))
    return expanded


def delumping(
    lumps: Sequence[Mapping[str, Any]],
    resolution: int,
) -> DelumpingResult:
    """A lumped pseudo-component table split into a finer grid.

    ``lumps`` is a list of records, one per pseudo-component, keyed by the field names the spec
    declares - each in the unit that field's declaration fixes. ``resolution`` is how many
    sub-fractions each lump becomes; one or less is the class's identity and returns the input
    unchanged.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    # The one conversion site, shared with the Rust bridge - a record field's unit is in the spec's
    # `fields` block, so the reader is `input_records_to_si` rather than `input_to_si`.
    rows: list[dict[str, Any]] = input_records_to_si(spec, "lumps", lumps)

    apply_checks(checks.on_input, lambda _name: None, warnings)

    fractions = [
        SubFraction(
            name=row["name"],
            moles=from_si(row["moles"], "mol"),
            molar_mass=from_si(row["molar_mass"], "kg/mol"),
            density=None if row["density"] is None else from_si(row["density"], "kg/m**3"),
            boiling_point=(
                None if row["boiling_point"] is None else from_si(row["boiling_point"], "K")
            ),
        )
        for row in split(rows, int(resolution))
    ]

    apply_checks(checks.derived, lambda _name: None, warnings)

    return DelumpingResult(fractions=tuple(fractions), warnings=tuple(warnings))
