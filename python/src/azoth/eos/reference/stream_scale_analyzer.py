"""``eos.stream_scale_analyzer`` - NeqSim's ``StreamScaleAnalyzer``.

The coupled multi-mineral equilibrium at one litre of water, with each mineral's
precipitated mass converted to a rate by the aqueous throughput, and the dominant
mineral reported.

Spec: ``specs/models/eos/stream_scale_analyzer.toml``.
"""

from __future__ import annotations

from azoth.core.record_outputs_gen import StreamScaleMineralRow
from azoth.core.result import StreamScaleAnalyzerResult
from azoth.core.units import Q, from_si
from azoth.eos.reference.multi_mineral_scale_equilibrium import (
    multi_mineral_scale_equilibrium,
)

CALC_ID = "eos.stream_scale_analyzer"


def stream_scale_analyzer(
    calcium: Q,
    barium: Q,
    strontium: Q,
    iron: Q,
    magnesium: Q,
    sodium: Q,
    bicarbonate: Q,
    sulphate: Q,
    total_dissolved_solids: Q,
    temperature: Q,
    pressure: Q,
    co2_partial_pressure: Q,
    ph: float,
    auto_ph: bool,
    second_order_pressure_correction: bool,
    activity_model: str,
    water_flow: Q,
    bdot_ion_size: float | None = None,
    bdot_b: float | None = None,
    pitzer_ionic_strength_molal: float | None = None,
) -> StreamScaleAnalyzerResult:
    """The scale mass rates of an aqueous stream.

    Args:
        calcium: the calcium ion's concentration.
        barium: the barium ion's concentration.
        strontium: the strontium ion's concentration.
        iron: the ferrous iron's concentration.
        magnesium: the magnesium ion's concentration.
        sodium: the sodium ion's concentration.
        bicarbonate: the bicarbonate ion's concentration.
        sulphate: the sulphate ion's concentration.
        total_dissolved_solids: TDS.
        temperature: absolute temperature.
        pressure: absolute pressure.
        co2_partial_pressure: CO2 partial pressure, read by the auto-pH path.
        ph: the water's pH.
        auto_ph: whether to estimate the pH with the class's own Henry constant.
        second_order_pressure_correction: whether the Ksp correction carries its term.
        activity_model: ``davies``, ``bdot`` or ``pitzer_binary``.
        water_flow: the aqueous volumetric throughput.
        bdot_ion_size: the B-dot ion size; absent takes the class's 5 Angstrom.
        bdot_b: the B-dot Helgeson coefficient; absent takes 0.041.
        pitzer_ionic_strength_molal: the Pitzer strength; absent uses the class's estimate.

    Returns:
        One row per registered mineral and the total mass, rate and dominant mineral.

    Raises:
        OutOfRangeError: if an input is outside a declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = stream_scale_analyzer(
        ...     q(1.0, "kg/m**3"), q(1.0, "kg/m**3"), q(3.0, "kg/m**3"),
        ...     q(0.0, "kg/m**3"), q(0.0, "kg/m**3"), q(0.0, "kg/m**3"),
        ...     q(0.2, "kg/m**3"), q(1.5, "kg/m**3"), q(60.0, "kg/m**3"),
        ...     q(343.15, "K"), q(1.0e7, "Pa"), q(1.0e5, "Pa"), 0.0,
        ...     False, False, "davies",
        ...     water_flow=q(0.0011574074074074073, "m**3/s"))
        >>> r.dominant_scale
        'BaSO4'
    """
    volume = from_si(1.0e-3, "m**3")  # one litre
    equilibrium = multi_mineral_scale_equilibrium(
        calcium,
        barium,
        strontium,
        iron,
        magnesium,
        sodium,
        bicarbonate,
        sulphate,
        total_dissolved_solids,
        temperature,
        pressure,
        co2_partial_pressure,
        ph,
        auto_ph,
        second_order_pressure_correction,
        activity_model,
        bdot_ion_size,
        bdot_b,
        pitzer_ionic_strength_molal,
        water_volume=volume,
    )
    flow = float(water_flow.to("m**3/s").magnitude)

    rows: list[StreamScaleMineralRow] = []
    dominant = "none"
    best = 0.0
    for mineral in equilibrium.minerals:
        mass = float(mineral.precipitated_mass.to("kg/m**3").magnitude)
        rows.append(
            StreamScaleMineralRow(
                mineral=mineral.name,
                initial_si=mineral.initial_saturation_index,
                final_si=mineral.final_saturation_index,
                precipitated_mass=mineral.precipitated_mass,
                scale_rate=from_si(mass * flow, "kg/s"),
            )
        )
        if mass > best:
            best = mass
            dominant = mineral.name

    total_mass = float(equilibrium.total_mass.to("kg/m**3").magnitude)
    return StreamScaleAnalyzerResult(
        minerals=tuple(rows),
        total_scale_mass=equilibrium.total_mass,
        total_scale_rate=from_si(total_mass * flow, "kg/s"),
        dominant_scale=dominant,
        has_scaling=total_mass > 0.0,
        warnings=tuple(equilibrium.warnings),
    )
