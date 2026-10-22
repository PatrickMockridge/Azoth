"""``eos.pipe_beggs_and_brills`` - NeqSim's ``PipeBeggsAndBrills``, isothermal march.

A pipeline is divided into increments; each is flashed at the inlet temperature and its own
pressure, its phase properties feed the Beggs-and-Brill segment correlation, and the resulting
drop sets the next segment's pressure.

Spec: ``specs/models/eos/pipe_beggs_and_brills.toml``.
"""

from __future__ import annotations

from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import BeggsBrillSegmentResult, PipeBeggsAndBrillsResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import mixture_of
from azoth.eos.reference.aqueous_viscosity import aqueous_viscosity
from azoth.eos.reference.beggs_brill_segment import beggs_brill_segment
from azoth.eos.reference.phase_transport import phase_transport
from azoth.eos.reference.tp_multiflash import tp_multiflash
from azoth.eos.reference.viscosity import viscosity

CALC_ID = "eos.pipe_beggs_and_brills"

MOLAR_GAS_CONSTANT = 8.314_462_618_153_24
SRK_OMEGA_B = 0.086_64
WATER_NAMES = {"water", "water_PC", "H2O", "waterpc"}


def _phase_kind(
    composition: list[float],
    components: list[Any],
    names: list[str],
    z_factor: float,
    t: float,
    p: float,
) -> str:
    volume = z_factor * MOLAR_GAS_CONSTANT * t / p
    b = sum(
        fraction
        * SRK_OMEGA_B
        * MOLAR_GAS_CONSTANT
        * component.Tc.to("K").magnitude
        / component.Pc.to("Pa").magnitude
        for fraction, component in zip(composition, components, strict=True)
    )
    if b > 0.0 and volume / b > 1.75:
        return "gas"
    water = sum(
        fraction for fraction, name in zip(composition, names, strict=True) if name in WATER_NAMES
    )
    return "oil" if 1.0 - water > water else "aqueous"


def _viscosity(
    mixture: Any, names: list[str], kind: str, t: Q, p: Q, composition: list[float]
) -> float:
    if len(names) > 1:
        return phase_transport(names, kind, t, p, composition).mu.to("Pa*s").magnitude
    if kind == "aqueous":
        return aqueous_viscosity(mixture, t, p, composition).viscosity.to("Pa*s").magnitude
    return viscosity(mixture, t, p, composition).mu.to("Pa*s").magnitude


def _phase_properties(
    mixture: Any, names: list[str], t: Q, p: Q, n: float, overall: list[float]
) -> dict[str, float]:
    flash = tp_multiflash(mixture, t, p, overall)
    components = mixture.components
    t_k = t.to("K").magnitude
    p_pa = p.to("Pa").magnitude
    gas: dict[str, float] | None = None
    liquid_flow = 0.0
    liquid_mass = 0.0
    liquid_weight = 0.0
    liquid_mu = 0.0
    for index, fraction in enumerate(flash.phase_fractions):
        if fraction <= 0.0:
            continue
        composition = list(flash.x[index])
        z_factor = flash.z_factor[index]
        shift = sum(f * c.volume_shift_srk for f, c in zip(composition, components, strict=True))
        volume = z_factor * MOLAR_GAS_CONSTANT * t_k / p_pa - shift
        mass = sum(
            f * c.molar_mass.to("kg/mol").magnitude
            for f, c in zip(composition, components, strict=True)
        )
        density = mass / volume
        phase_flow = n * fraction * volume
        kind = _phase_kind(composition, components, names, z_factor, t_k, p_pa)
        mu = _viscosity(mixture, names, kind, t, p, composition)
        if kind == "gas" and gas is None:
            gas = {"flow": phase_flow, "density": density, "mu": mu}
        else:
            liquid_flow += phase_flow
            liquid_mass += phase_flow * density
            liquid_weight += phase_flow
            liquid_mu += phase_flow * mu
    if liquid_weight > 0.0:
        liquid_density = liquid_mass / liquid_weight
        liquid_viscosity = liquid_mu / liquid_weight
        return {
            "gas_flow": 0.0 if gas is None else gas["flow"],
            "liquid_flow": liquid_flow,
            "gas_density": liquid_density if gas is None else gas["density"],
            "liquid_density": liquid_density,
            "gas_viscosity": liquid_viscosity if gas is None else gas["mu"],
            "liquid_viscosity": liquid_viscosity,
        }
    if gas is None:
        raise InvalidInputError("inlet_n", "the flash left no phase to flow")
    return {
        "gas_flow": gas["flow"],
        "liquid_flow": 0.0,
        "gas_density": gas["density"],
        "liquid_density": gas["density"],
        "gas_viscosity": gas["mu"],
        "liquid_viscosity": gas["mu"],
    }


def pipe_beggs_and_brills(
    components: list[str],
    inlet_n: float,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    length: Q,
    diameter: Q,
    elevation: Q,
    number_of_increments: int,
    roughness: Q,
    eos: str | None = None,
) -> PipeBeggsAndBrillsResult:
    """The isothermal pressure profile of a Beggs-and-Brills pipeline.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature, which the isothermal march holds.
        length: the pipe's total length.
        diameter: the pipe's internal diameter.
        elevation: the total elevation change, positive uphill.
        number_of_increments: the number of segments.
        roughness: the pipe wall's absolute roughness.
        eos: the cubic the flash runs; absent takes ``srk``.

    Returns:
        The pressure, drop, holdup, velocity and regime profiles and the outlet state.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.
        InvalidInputError: if a segment's outlet pressure is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pipe_beggs_and_brills(
        ...     ["methane", "CO2", "water"], 80.0, [0.85, 0.05, 0.10],
        ...     q(1.0e7, "Pa"), q(333.15, "K"), q(1000.0, "m"), q(0.1, "m"),
        ...     q(0.0, "m"), 10, q(1.0e-5, "m"), "srk")
        >>> len(r.pressure_profile)
        11
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "inlet_n": inlet_n,
            "inlet_p": inlet_p,
            "inlet_t": inlet_t,
            "length": length,
            "diameter": diameter,
            "elevation": elevation,
            "roughness": roughness,
        }.items()
    }
    values["number_of_increments"] = number_of_increments
    apply_checks(checks.on_input, values.get, warnings)

    cubic = eos or "srk"
    if cubic not in ("pr", "srk"):
        raise InvalidInputError("eos", f"`{cubic}` is not one of `pr` or `srk`")
    mixture, _ideal = mixture_of(components, eos=cubic)

    increments = int(number_of_increments)
    segment_length = values["length"] / increments
    segment_elevation = values["elevation"] / increments

    pressure = values["inlet_p"]
    pressure_profile = [from_si(pressure, "Pa")]
    drop_profile = [from_si(0.0, "Pa")]
    holdup_profile: list[float] = []
    velocity_profile: list[Q] = []
    regime_profile: list[str] = []

    def drop_at(state_pressure: float) -> BeggsBrillSegmentResult:
        properties = _phase_properties(
            mixture,
            components,
            inlet_t,
            from_si(state_pressure, "Pa"),
            values["inlet_n"],
            list(inlet_z),
        )
        return beggs_brill_segment(
            from_si(properties["gas_flow"], "m**3/s"),
            from_si(properties["liquid_flow"], "m**3/s"),
            from_si(properties["gas_density"], "kg/m**3"),
            from_si(properties["liquid_density"], "kg/m**3"),
            from_si(properties["gas_viscosity"], "Pa*s"),
            from_si(properties["liquid_viscosity"], "Pa*s"),
            inlet_t,
            from_si(state_pressure, "Pa"),
            diameter,
            from_si(segment_length, "m"),
            from_si(segment_elevation, "m"),
            roughness,
        )

    for _ in range(increments):
        drop = drop_at(pressure)
        holdup_profile.append(drop.liquid_holdup)
        velocity_profile.append(drop.superficial_mixture_velocity)
        regime_profile.append(drop.regime.value)
        outlet = pressure - drop.total_drop.to("Pa").magnitude
        if outlet <= 0.0:
            raise InvalidInputError(
                "inlet_p",
                f"a segment's outlet pressure is {outlet} Pa: "
                "the line does not deliver at this inlet",
            )
        pressure = outlet
        pressure_profile.append(from_si(pressure, "Pa"))
        drop_profile.append(drop.total_drop)

    final_drop = drop_at(pressure)
    holdup_profile.append(final_drop.liquid_holdup)
    velocity_profile.append(final_drop.superficial_mixture_velocity)
    regime_profile.append(final_drop.regime.value)

    return PipeBeggsAndBrillsResult(
        pressure_profile=tuple(pressure_profile),
        pressure_drop_profile=tuple(drop_profile),
        temperature_profile=(from_si(values["inlet_t"], "K"),),
        liquid_holdup_profile=tuple(holdup_profile),
        mixture_superficial_velocity_profile=tuple(velocity_profile),
        flow_regime_profile=tuple(regime_profile),
        outlet_pressure=from_si(pressure, "Pa"),
        total_pressure_drop=from_si(values["inlet_p"] - pressure, "Pa"),
        warnings=tuple(warnings),
    )
