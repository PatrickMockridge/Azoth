"""``eos.surf_cooldown`` - NeqSim's ``SurfCooldownAnalyzer``.

The hydrate equilibrium temperature at the operating pressure, the fluid's density and
specific heat at the operating point, and the lumped cooldown that follows: the no-touch
time is when the fluid reaches the hydrate temperature plus the margin.

Spec: ``specs/models/eos/surf_cooldown.toml``.
"""

from __future__ import annotations

from typing import Any

from azoth.core.errors import AzothError, InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CooldownVerdict, SurfCooldownResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.hydrate_formation_temperature import hydrate_formation_temperature
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.pipeline_cooldown import pipeline_cooldown
from azoth.eos.reference.tp_multiflash import tp_multiflash

CALC_ID = "eos.surf_cooldown"

MOLAR_GAS_CONSTANT = 8.314_462_618_153_24
DEFAULT_TOTAL_TIME = 259_200.0
DEFAULT_TIME_STEP = 300.0
DEFAULT_HYDRATE_MARGIN = 3.0
DEFAULT_STEEL_DENSITY = 7850.0
DEFAULT_STEEL_SPECIFIC_HEAT = 500.0
WATER_PRESENCE_TOLERANCE = 1.0e-12


def _fluid_properties(
    mixture: Any, ideal_gas: Any, z: list[float], t: Q, p: Q, cubic: str
) -> tuple[float, float]:
    """The flashed phases' total density and mass specific heat."""
    flash = tp_multiflash(mixture, t, p, z)
    components = mixture.components
    t_k = t.to("K").magnitude
    p_pa = p.to("Pa").magnitude
    total_mass = 0.0
    total_volume = 0.0
    total_cp = 0.0
    for index, fraction in enumerate(flash.phase_fractions):
        if fraction <= 0.0:
            continue
        composition = list(flash.x[index])
        z_factor = flash.z_factor[index]
        shift = sum(
            xi * (c.volume_shift_srk if cubic == "srk" else c.volume_shift)
            for xi, c in zip(composition, components, strict=True)
        )
        volume = z_factor * MOLAR_GAS_CONSTANT * t_k / p_pa - shift
        mass = sum(
            xi * c.molar_mass.to("kg/mol").magnitude
            for xi, c in zip(composition, components, strict=True)
        )
        cp = (
            molar_enthalpy_entropy(mixture, ideal_gas, t, p, composition, z_factor)
            .cp.to("J/(mol*K)")
            .magnitude
        )
        total_mass += fraction * mass
        total_volume += fraction * volume
        total_cp += fraction * cp
    return total_mass / total_volume, total_cp / total_mass


def surf_cooldown(
    components: list[str],
    z: list[float],
    operating_pressure: Q,
    flash_temperature: Q,
    internal_diameter: Q,
    wall_thickness: Q,
    insulation_thickness: Q,
    coating_thickness: Q,
    insulation_conductivity: Q,
    seabed_temperature: Q,
    operating_temperature: Q | None = None,
    overall_u_value: Q | None = None,
    external_htc: Q | None = None,
    hydrate_margin: Q | None = None,
    required_no_touch_time: Q | None = None,
    total_time: Q | None = None,
    time_step: Q | None = None,
    require_water: bool | None = None,
    eos: str | None = None,
    hydrate_model: str | None = None,
) -> SurfCooldownResult:
    """The surface-cooldown no-touch-time analysis.

    Args:
        components: the fluid's substances, by name.
        z: the overall mole fractions.
        operating_pressure: the pressure the fluid and the hydrate equilibrium are read at.
        flash_temperature: the temperature the density and specific heat are flashed at.
        internal_diameter: the pipe's internal diameter.
        wall_thickness: the steel wall thickness.
        insulation_thickness: the insulation thickness.
        coating_thickness: the external coating thickness.
        insulation_conductivity: the insulation's thermal conductivity.
        seabed_temperature: the ambient temperature.
        operating_temperature: the cooldown's initial temperature; absent takes the flash one.
        overall_u_value: a supplied U value; absent computes the layers.
        external_htc: the external coefficient; absent takes 300.
        hydrate_margin: the margin above the hydrate temperature; absent takes 3 K.
        required_no_touch_time: the requirement; absent uses the absolute bands.
        total_time: the horizon; absent takes 72 hours.
        time_step: the march step; absent takes five minutes.
        require_water: whether an absent water is refused; absent is false.
        eos: the cubic; absent takes srk.
        hydrate_model: the cage model; absent takes pvtsim.

    Returns:
        The verdict, the water state, the no-touch time and the properties the march used.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.
        InvalidInputError: if the fluid carries no water and `require_water` is set.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = surf_cooldown(
        ...     ["methane", "ethane", "propane", "water"], [0.82, 0.06, 0.03, 0.09],
        ...     q(1.2e7, "Pa"), q(338.15, "K"), q(0.254, "m"), q(0.0159, "m"),
        ...     q(0.060, "m"), q(0.006, "m"), q(0.22, "W/(m*K)"), q(277.15, "K"),
        ...     overall_u_value=q(2.5, "W/(m**2*K)"))
        >>> r.water_present
        True
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "operating_pressure": operating_pressure,
            "flash_temperature": flash_temperature,
            "internal_diameter": internal_diameter,
            "wall_thickness": wall_thickness,
            "insulation_thickness": insulation_thickness,
            "coating_thickness": coating_thickness,
            "insulation_conductivity": insulation_conductivity,
            "seabed_temperature": seabed_temperature,
        }.items()
    }
    optional = {
        "operating_temperature": operating_temperature,
        "overall_u_value": overall_u_value,
        "external_htc": external_htc,
        "hydrate_margin": hydrate_margin,
        "required_no_touch_time": required_no_touch_time,
        "total_time": total_time,
        "time_step": time_step,
    }
    for name, value in optional.items():
        if value is not None:
            values[name] = input_to_si(spec, name, value)
    apply_checks(checks.on_input, values.get, warnings)

    if len(z) != len(components):
        raise InvalidInputError(
            "z",
            f"the composition has {len(z)} entries for {len(components)} components",
        )
    water_index = next(
        (index for index, name in enumerate(components) if name.lower() == "water"), None
    )
    water_mole_fraction = 0.0 if water_index is None else z[water_index]
    water_present = water_mole_fraction > WATER_PRESENCE_TOLERANCE
    if not water_present and require_water:
        raise InvalidInputError(
            "components",
            "the fluid carries no water, so no hydrate equilibrium temperature can be "
            "computed and the no-touch time is not meaningful; add water or leave "
            "`require_water` unset",
        )

    from azoth.eos.components import mixture_of

    cubic = eos or "srk"
    mixture, ideal_gas = mixture_of(components, eos=cubic)
    fluid_density, fluid_specific_heat = _fluid_properties(
        mixture,
        ideal_gas,
        z,
        from_si(values["flash_temperature"], "K"),
        from_si(values["operating_pressure"], "Pa"),
        cubic,
    )

    try:
        hydrate_temperature = (
            hydrate_formation_temperature(
                components,
                from_si(values["operating_pressure"], "Pa"),
                z,
                cubic,
                hydrate_model or "pvtsim",
            )
            .temperature.to("K")
            .magnitude
        )
    except AzothError:
        hydrate_temperature = None

    margin = DEFAULT_HYDRATE_MARGIN if hydrate_margin is None else values["hydrate_margin"]
    initial_temperature = values.get("operating_temperature", values["flash_temperature"])
    target = None if hydrate_temperature is None else from_si(hydrate_temperature + margin, "K")
    cooldown = pipeline_cooldown(
        from_si(values["internal_diameter"], "m"),
        from_si(values["wall_thickness"], "m"),
        from_si(values["insulation_thickness"], "m"),
        from_si(values["coating_thickness"], "m"),
        from_si(DEFAULT_STEEL_DENSITY, "kg/m**3"),
        from_si(DEFAULT_STEEL_SPECIFIC_HEAT, "J/(kg*K)"),
        from_si(values["insulation_conductivity"], "W/(m*K)"),
        from_si(fluid_density, "kg/m**3"),
        from_si(fluid_specific_heat, "J/(kg*K)"),
        from_si(initial_temperature, "K"),
        from_si(values["seabed_temperature"], "K"),
        from_si(values.get("time_step", DEFAULT_TIME_STEP), "s"),
        from_si(values.get("total_time", DEFAULT_TOTAL_TIME), "s"),
        overall_u_value=(
            None if overall_u_value is None else from_si(values["overall_u_value"], "W/(m**2*K)")
        ),
        use_layer_calculation=overall_u_value is None,
        external_htc=(
            None if external_htc is None else from_si(values["external_htc"], "W/(m**2*K)")
        ),
        target_temperature=target,
    )

    required = values.get("required_no_touch_time")

    def classify(duration: float) -> CooldownVerdict:
        if required is not None and required > 0.0:
            if duration >= required:
                return CooldownVerdict.OK
            if duration >= 0.75 * required:
                return CooldownVerdict.MARGINAL
            return CooldownVerdict.CRITICAL
        if duration >= 12.0 * 3600.0:
            return CooldownVerdict.OK
        if duration >= 6.0 * 3600.0:
            return CooldownVerdict.MARGINAL
        return CooldownVerdict.CRITICAL

    no_touch: float | None
    verdict: CooldownVerdict
    if hydrate_temperature is None or hydrate_temperature <= values["seabed_temperature"]:
        no_touch, verdict = None, CooldownVerdict.NO_HYDRATE_RISK
    elif initial_temperature <= hydrate_temperature + margin:
        no_touch, verdict = 0.0, classify(0.0)
    elif cooldown.time_to_reach_target is None:
        no_touch, verdict = None, CooldownVerdict.OK
    else:
        no_touch = cooldown.time_to_reach_target.to("s").magnitude
        verdict = classify(no_touch)

    return SurfCooldownResult(
        verdict=verdict,
        water_present=water_present,
        water_mole_fraction=water_mole_fraction,
        no_touch_time=None if no_touch is None else from_si(no_touch, "s"),
        initial_fluid_temperature=from_si(initial_temperature, "K"),
        hydrate_equilibrium_temperature=(
            None if hydrate_temperature is None else from_si(hydrate_temperature, "K")
        ),
        hydrate_margin=from_si(margin, "K"),
        seabed_temperature=from_si(values["seabed_temperature"], "K"),
        overall_u_value=cooldown.effective_u_value,
        time_constant=cooldown.time_constant,
        fluid_density=from_si(fluid_density, "kg/m**3"),
        fluid_specific_heat=from_si(fluid_specific_heat, "J/(kg*K)"),
        internal_diameter=from_si(values["internal_diameter"], "m"),
        insulation_thickness=from_si(values["insulation_thickness"], "m"),
        warnings=tuple(warnings),
    )
