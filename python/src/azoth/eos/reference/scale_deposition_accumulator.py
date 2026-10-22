"""``eos.scale_deposition_accumulator`` - NeqSim's ``ScaleDepositionAccumulator``.

A segmented pipeline's scale deposition from a temperature and pressure profile: the CaCO3
index per segment, its deposit mass and wall thickness, the total, and the time to half-block
the worst segment.

Spec: ``specs/models/eos/scale_deposition_accumulator.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import DepositionSegmentRow
from azoth.core.result import ScaleDepositionAccumulatorResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.scale_prediction import _state

CALC_ID = "eos.scale_deposition_accumulator"


def scale_deposition_accumulator(
    calcium: Q,
    bicarbonate: Q,
    sulphate: Q,
    barium: Q,
    sodium: Q,
    total_dissolved_solids: Q,
    ph: float,
    diameter: Q,
    length: Q,
    segment_temperatures: list[Q],
    segment_pressures: list[Q],
    service_time: Q,
    scale_density: Q,
    base_rate: Q,
    inhibitor_efficiency: float,
) -> ScaleDepositionAccumulatorResult:
    """The segmented deposition of a pipeline profile.

    Args:
        calcium: the calcium ion's concentration.
        bicarbonate: the bicarbonate ion's concentration.
        sulphate: the sulphate ion's concentration.
        barium: the barium ion's concentration.
        sodium: the sodium ion's concentration.
        total_dissolved_solids: TDS.
        ph: the aqueous phase's pH.
        diameter: the pipe's internal diameter.
        length: the pipe's length.
        segment_temperatures: the profile's segment temperatures.
        segment_pressures: the profile's segment pressures; the shorter list governs.
        service_time: the exposure time the masses are scaled by.
        scale_density: the assumed deposit density.
        base_rate: the deposition rate at SI = 1.
        inhibitor_efficiency: the uniform fractional efficiency, 0..1.

    Returns:
        One row per segment, the total mass, the worst thickness and the blockage time.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = scale_deposition_accumulator(
        ...     q(1.0, "kg/m**3"), q(0.5, "kg/m**3"), q(0.1, "kg/m**3"),
        ...     q(0.01, "kg/m**3"), q(20.0, "kg/m**3"), q(50.0, "kg/m**3"), 6.5,
        ...     q(0.1, "m"), q(1000.0, "m"), [q(333.15, "K")], [q(1.0e7, "Pa")],
        ...     q(31_557_600.0, "s"), q(2700.0, "kg/m**3"),
        ...     q(1.5844043908606944e-8, "kg/(m**2*s)"), 0.0)
        >>> abs(r.total_deposition_mass.magnitude - 5.112073003689538) < 1e-9
        True
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "calcium": calcium,
            "bicarbonate": bicarbonate,
            "sulphate": sulphate,
            "barium": barium,
            "sodium": sodium,
            "total_dissolved_solids": total_dissolved_solids,
            "diameter": diameter,
            "length": length,
            "service_time": service_time,
            "scale_density": scale_density,
            "base_rate": base_rate,
        }.items()
    }
    values["ph"] = ph
    values["inhibitor_efficiency"] = inhibitor_efficiency
    apply_checks(checks.on_input, values.get, warnings)

    temperatures = [
        input_to_si(spec, "segment_temperatures", value) for value in segment_temperatures
    ]
    pressures = [input_to_si(spec, "segment_pressures", value) for value in segment_pressures]

    service = max(0.0, values["service_time"])
    efficiency = min(1.0, max(0.0, inhibitor_efficiency))
    count = min(len(temperatures), len(pressures))

    if values["diameter"] <= 0.0 or values["length"] <= 0.0 or count == 0:
        return ScaleDepositionAccumulatorResult(
            segments=(),
            total_deposition_mass=from_si(0.0, "kg"),
            max_thickness=from_si(0.0, "m"),
            time_to_blockage=from_si(math.inf, "s"),
            warnings=tuple(warnings),
        )

    segment_length = values["length"] / count
    segment_area = math.pi * values["diameter"] * segment_length

    rows: list[DepositionSegmentRow] = []
    total_mass = 0.0
    max_thickness = 0.0
    for index in range(count):
        point = _state(
            {
                "calcium": values["calcium"],
                "barium": values["barium"],
                "strontium": 0.0,
                "iron": 0.0,
                "magnesium": 0.0,
                "sodium": values["sodium"],
                "bicarbonate": values["bicarbonate"],
                "sulphate": values["sulphate"],
                "total_dissolved_solids": values["total_dissolved_solids"],
            },
            temperatures[index],
            pressures[index],
            0.0,
            ph,
            False,
            False,
        )
        saturation_index = point["si"][0]
        driving = max(0.0, math.pow(10.0, saturation_index) - 1.0)
        rate = values["base_rate"] * driving * (1.0 - efficiency)
        mass = rate * segment_area * service
        thickness = mass / (values["scale_density"] * segment_area)

        total_mass += mass
        max_thickness = max(max_thickness, thickness)
        rows.append(
            DepositionSegmentRow(
                saturation_index=saturation_index,
                deposition_mass=from_si(mass, "kg"),
                thickness=from_si(thickness, "m"),
            )
        )

    blockage = (
        values["diameter"] / 4.0 * service / max_thickness
        if max_thickness > 0.0 and service > 0.0
        else math.inf
    )

    return ScaleDepositionAccumulatorResult(
        segments=tuple(rows),
        total_deposition_mass=from_si(total_mass, "kg"),
        max_thickness=from_si(max_thickness, "m"),
        time_to_blockage=from_si(blockage, "s"),
        warnings=tuple(warnings),
    )
