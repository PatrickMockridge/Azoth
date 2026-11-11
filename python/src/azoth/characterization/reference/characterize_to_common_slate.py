"""``characterization.characterize_to_common_slate`` - several fluids onto one shared slate.

Spec: ``specs/models/characterization/characterize_to_common_slate.toml``. Oracle:
``validation/neqsim/captures/characterization_fluid_rebuild_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/characterize_to_common_slate.rs``; the arithmetic it
calls is :mod:`azoth.characterization.reference._combiner`'s.

# What the slate is

Unlike ``combine_reservoir_fluids``, the fluids are not blended: each comes back separate, with
its **own** mole amounts but an identical set of pseudo-component properties. That is Pedersen's
"common EoS" slate - several fluids, one equation of state - and it is what makes two fluids cut
differently comparable. A fluid's physical mass can move, because it keeps its moles and takes
the slate's molar mass; ``fluid_moles`` is reported so that is visible.

# The cut names are positions in the slate

A cut invalid for any fluid is dropped from **every** fluid, and the names keep their slate
positions - ``PC1_PC``, ``PC3_PC`` - so a cut can be matched across the fluids by name.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import Any

from azoth import _models_gen
from azoth.characterization.reference import _combiner
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CharacterizeToCommonSlateResult
from azoth.core.units import from_si, input_records_to_si, to_count
from azoth.core.warnings import Warning

MODEL_ID = "characterization.characterize_to_common_slate"


def _resolve_weights(weights: Sequence[float] | None, count: int) -> list[float]:
    """``resolveWeights``: one weight per fluid, all ones when none is given."""
    if weights is None:
        return [1.0] * count
    resolved = [float(weight) for weight in weights]
    if len(resolved) != count:
        raise InvalidInputError(
            "weights",
            f"weights length ({len(resolved)}) must match number of fluids ({count})",
        )
    total = 0.0
    for weight in resolved:
        if weight < 0.0:
            raise InvalidInputError("weights", "weights must be non-negative")
        total += weight
    if not total > 0.0:
        raise InvalidInputError("weights", "at least one weight must be positive")
    return resolved


def _resolve_target(
    target_pseudo_components: int | None, extractions: Sequence[_combiner.Extraction]
) -> int:
    """The slate size: supplied, or the largest pseudo-component count, at least one."""
    if target_pseudo_components is None:
        return max(max((len(e.pseudo) for e in extractions), default=0), 1)
    target = to_count(target_pseudo_components, "target_pseudo_components")
    if target == 0:
        raise OutOfRangeError(
            "target_pseudo_components",
            0.0,
            "a slate needs at least one pseudo-component",
        )
    return target


def characterize_to_common_slate(
    components: Sequence[Mapping[str, Any]],
    weights: Sequence[float] | None = None,
    target_pseudo_components: int | None = None,
) -> CharacterizeToCommonSlateResult:
    """Several reservoir fluids re-characterized onto one shared pseudo-component slate.

    ``components`` is a list of records, one per component of every fluid, keyed by the field
    names the spec declares. A row states its fluid's index, whether it is a pseudo-component,
    its name and its mole amount; a pseudo-component states its molar mass, density and the
    properties the slate averages.

    ``weights`` is one non-negative weight per fluid, in first-appearance order; absent is equal
    weighting, and at least one must be positive. ``target_pseudo_components`` is the slate's
    size; absent takes the largest pseudo-component count among the fluids, and at least one.

    The answer is one row list with a parallel ``fluid_index``; ``fluid_moles`` is each fluid's
    total after re-characterization, and ``resolved_pseudo_components`` says which slate size was
    used.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    rows: list[dict[str, Any]] = input_records_to_si(spec, "components", components)

    extractions = _combiner.extract(rows)
    if not extractions:
        raise InvalidInputError("components", "a common slate needs at least one fluid")
    resolved_weights = _resolve_weights(weights, len(extractions))
    target = _resolve_target(target_pseudo_components, extractions)

    fluid_mole_totals = [
        sum(extraction.base.values()) + sum(row.moles for row in extraction.pseudo)
        for extraction in extractions
    ]

    result_rows = []
    fluid_index: list[float] = []
    fluid_moles: list[float] = []
    any_pseudo = any(extraction.pseudo for extraction in extractions)

    if any_pseudo:
        pool: list[_combiner.Contribution] = []
        for fluid, extraction in enumerate(extractions):
            if not extraction.pseudo or not resolved_weights[fluid] > 0.0:
                continue
            pseudo_mass = sum(row.mass for row in extraction.pseudo)
            if not pseudo_mass > _combiner.MASS_TOLERANCE:
                continue
            scale = resolved_weights[fluid] / pseudo_mass
            pool.extend(row.scaled_by_mass(scale) for row in extraction.pseudo)
        pool.sort(key=_combiner.Contribution.sorting_key)

        boundaries = _combiner.quantile_boundaries(pool, target)
        per_fluid = [
            _combiner.distribute(extraction.pseudo, boundaries, target)
            for extraction in extractions
        ]
        common = _combiner.compute_common_slate(
            per_fluid, fluid_mole_totals, resolved_weights, target
        )

        for fluid, extraction in enumerate(extractions):
            total = 0.0
            for name, moles in extraction.base.items():
                result_rows.append(_combiner.base_row(name, moles))
                fluid_index.append(float(fluid))
                total += moles
            for cut, (shared, own) in enumerate(zip(common, per_fluid[fluid], strict=True)):
                if not shared.is_valid() or not own.is_valid():
                    continue
                result_rows.append(
                    _combiner.profile_row(f"PC{cut + 1}_PC", replace(shared, moles=own.moles))
                )
                fluid_index.append(float(fluid))
                total += own.moles
            fluid_moles.append(total)
    else:
        # The class's own early return: the fluids come back as themselves.
        for fluid, extraction in enumerate(extractions):
            total = 0.0
            for name, moles in extraction.base.items():
                result_rows.append(_combiner.base_row(name, moles))
                fluid_index.append(float(fluid))
                total += moles
            fluid_moles.append(total)

    apply_checks(checks.derived, lambda _name: None, warnings)

    return CharacterizeToCommonSlateResult(
        components=tuple(result_rows),
        fluid_index=tuple(fluid_index),
        fluid_moles=tuple(from_si(moles, "mol") for moles in fluid_moles),
        resolved_pseudo_components=target,
        warnings=tuple(warnings),
    )
