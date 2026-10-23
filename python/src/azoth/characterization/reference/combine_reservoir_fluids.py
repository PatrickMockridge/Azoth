"""``characterization.combine_reservoir_fluids`` - several fluids blended onto one slate.

Spec: ``specs/models/characterization/combine_reservoir_fluids.toml``. Oracle:
``validation/neqsim/captures/characterization_fluid_rebuild_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/combine_reservoir_fluids.rs`` line for line; the
arithmetic it calls is :mod:`azoth.characterization.reference._combiner`'s, which is the class's
own extraction, boundary, distribution and blend helpers.

# What the blend is

Each fluid's base components merge by name, its pseudo-components are pooled, and the pooled mass
is cut into ``target_pseudo_components`` equal-mass groups. Each fluid's own list is then
distributed onto **those** groups, and the per-fluid groups blended: so a cut's properties are a
mixture of the fluids' own cuts on either side of it, weighted by how much of the blend each fluid
is. The cut points are keys rather than interpolations, and a target nothing reaches is an
infinity, so a slate larger than the pool leaves groups empty rather than being refused.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from azoth import _models_gen
from azoth.characterization.reference import _combiner
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CombineReservoirFluidsResult
from azoth.core.units import from_si, input_records_to_si, to_count
from azoth.core.warnings import Warning

MODEL_ID = "characterization.combine_reservoir_fluids"


def combine_reservoir_fluids(
    target_pseudo_components: int,
    components: Sequence[Mapping[str, Any]],
) -> CombineReservoirFluidsResult:
    """Several reservoir fluids blended into one fluid with a common pseudo-component slate.

    ``components`` is a list of records, one per component of each fluid, keyed by the field names
    the spec declares - each in the unit that field's declaration fixes. A row states which fluid
    it belongs to, whether it is a pseudo-component, its name and its mole amount; a
    pseudo-component states its molar mass, its density and the properties the blend averages,
    and a base component states none of them.

    The answer merges the base components by name and renames the combined pseudo-components
    ``PC1_PC``, ``PC2_PC`` and so on. A pool with no pseudo-component is the class's early return:
    the merged base components and nothing else.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    rows: list[dict[str, Any]] = input_records_to_si(spec, "components", components)

    apply_checks(
        checks.on_input,
        {"target_pseudo_components": float(target_pseudo_components)}.get,
        warnings,
    )
    target = to_count(target_pseudo_components, "target_pseudo_components")

    extractions = _combiner.extract(rows)
    base = _combiner.merge_bases(extractions)
    total_moles = sum(base.values())
    result_rows = [_combiner.base_row(name, moles) for name, moles in base.items()]

    # The class's own early return: the pool is empty, so the requested slate never appears.
    pooled = [contribution for extraction in extractions for contribution in extraction.pseudo]
    if pooled:
        boundaries = _combiner.quantile_boundaries(pooled, target)
        per_fluid = [
            _combiner.distribute(extraction.pseudo, boundaries, target)
            for extraction in extractions
        ]
        fluid_mass_totals = [sum(profile.mass() for profile in profiles) for profiles in per_fluid]
        fluid_mole_totals = [sum(profile.moles for profile in profiles) for profiles in per_fluid]

        counter = 1
        for profile in _combiner.combine_profiles(per_fluid, fluid_mass_totals, fluid_mole_totals):
            if not profile.is_valid():
                continue
            total_moles += profile.moles
            result_rows.append(_combiner.profile_row(f"PC{counter}_PC", profile))
            counter += 1

    apply_checks(checks.derived, lambda _name: None, warnings)

    return CombineReservoirFluidsResult(
        components=tuple(result_rows),
        total_moles=from_si(total_moles, "mol"),
        warnings=tuple(warnings),
    )
