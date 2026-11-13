"""``process.gas_dispersion`` - the ambient-flash screening, as ``GasDispersionAnalyzer``.

Spec: ``specs/models/process/gas_dispersion.toml``

The Python twin of ``crates/azoth-process/src/models/gas_dispersion.rs``: the ambient TP
flash, the gas phase's fuel data by Le Chatelier's rule, the density-ratio model selection,
and the flammable and toxic endpoints over :mod:`azoth.process.reference.gaussian_plume` or
the class's Britter-McQuaid fit.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GasDispersionResult, Phase
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._entropy_state import molar_mass_of
from azoth.eos.reference.pr_mass_density import pr_mass_density
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.pt_flash import check_composition, pt_flash
from azoth.process.reference.gaussian_plume import (
    _sigma_y,
    _sigma_z,
    _stable,
    _terrain,
    concentration,
    distance_to_threshold,
)

MODEL_ID = "process.gas_dispersion"

#: ``GasDispersionAnalyzer.GAS_CONSTANT``, J/(mol K).
_GAS_CONSTANT = 8.314462618
#: ``ToxicLibrary``'s gas constant, J/(mol K), which the class does not share with the above.
_TOXIC_GAS_CONSTANT = 8.314
#: ``GasDispersionAnalyzer.AIR_GAS_CONSTANT``, J/(kg K).
_AIR_GAS_CONSTANT = 287.05
#: ``GasDispersionAnalyzer.DENSE_GAS_RATIO``.
_DENSE_GAS_RATIO = 1.2
#: ``GasDispersionAnalyzer.DEFAULT_RELEASE_HEIGHT_M``.
_DEFAULT_RELEASE_HEIGHT = 0.0
#: The BoundaryConditions defaults the calculation reads.
_DEFAULT_AMBIENT_TEMPERATURE = 288.15
_DEFAULT_ATMOSPHERIC_PRESSURE = 101325.0
_DEFAULT_WIND_SPEED = 5.0
_DEFAULT_SURFACE_ROUGHNESS = 0.1
_DEFAULT_STABILITY = "D"
#: The class's ``flammableEndpointFractionOfLfl`` default.
_DEFAULT_ENDPOINT_FRACTION = 0.5
#: The Gaussian cloud-volume integral's slab count.
_CLOUD_STEPS = 200

#: The class's lower-flammable-limit table, by normalized name, volume fraction.
_LOWER_FLAMMABLE_LIMITS = {
    "methane": 0.044,
    "ethane": 0.030,
    "propane": 0.021,
    "nbutane": 0.018,
    "ibutane": 0.018,
    "npentane": 0.014,
    "ipentane": 0.014,
    "nhexane": 0.012,
    "hydrogen": 0.040,
    "h2s": 0.040,
    "co": 0.125,
    "methanol": 0.060,
    "benzene": 0.012,
}

_SELECTIONS = ("auto", "gaussian_plume", "heavy_gas")


def _normalize(name: str) -> str:
    """The class's name normalization: lower case, with ``-``, ``_`` and spaces stripped."""
    return name.lower().replace("-", "").replace("_", "").replace(" ", "")


def _lfl_of(name: str) -> float | None:
    """The LFL of a component, by the class's table."""
    return _LOWER_FLAMMABLE_LIMITS.get(_normalize(name))


def _fuel_data(
    names: list[str],
    phase_z: list[float],
    molar_masses: list[float | None],
    phase_molar_mass: float,
    manual: float | None,
) -> tuple[float, float, float, float | None]:
    """The class's ``calculateFuelData``: mole fraction, mass fraction, molar mass, LFL."""
    if manual is not None:
        return 1.0, 1.0, phase_molar_mass, manual
    total_mass = 0.0
    fuel_mole = 0.0
    fuel_mass = 0.0
    reciprocal = 0.0
    for index, (name, fraction) in enumerate(zip(names, phase_z, strict=True)):
        mass = molar_masses[index]
        if mass is None:
            raise InvalidInputError("molar_mass", f"component {index} carries no molar mass")
        fraction = max(0.0, float(fraction))
        total_mass += fraction * mass
        lfl = _lfl_of(name)
        if lfl is not None and fraction > 0.0:
            fuel_mole += fraction
            fuel_mass += fraction * mass
            reciprocal += fraction / lfl
    if fuel_mole <= 0.0 or reciprocal <= 0.0 or total_mass <= 0.0:
        return 0.0, 0.0, math.nan, None
    return fuel_mole, fuel_mass / total_mass, fuel_mass / fuel_mole, fuel_mole / reciprocal


def _component_data(
    names: list[str],
    phase_z: list[float],
    molar_masses: list[float | None],
    target: str,
) -> tuple[float, float, float]:
    """The class's ``componentData``: mole fraction, mass fraction, molar mass."""
    wanted = _normalize(target)
    total_mass = 0.0
    component_mole = 0.0
    component_mass = 0.0
    component_molar_mass = math.nan
    for index, (name, fraction) in enumerate(zip(names, phase_z, strict=True)):
        mass = molar_masses[index]
        if mass is None:
            raise InvalidInputError("molar_mass", f"component {index} carries no molar mass")
        fraction = max(0.0, float(fraction))
        total_mass += fraction * mass
        if _normalize(name) == wanted:
            component_mole = fraction
            component_mass = fraction * mass
            component_molar_mass = mass
    if component_mole <= 0.0 or total_mass <= 0.0:
        return 0.0, 0.0, math.nan
    return component_mole, component_mass / total_mass, component_molar_mass


def _dense_distance(
    volumetric_rate: float,
    source_density: float,
    air_density: float,
    wind_speed: float,
    ratio: float,
) -> float | None:
    """The dense-gas distance for a concentration ratio; ``None`` outside ``(0, 1)``."""
    if math.isnan(ratio) or ratio <= 0.0 or ratio >= 1.0:
        return None
    g0 = 9.81 * (source_density - air_density) / air_density
    alpha = math.pow(g0 * g0 * volumetric_rate / math.pow(wind_speed, 5), 0.2)
    scale = math.sqrt(volumetric_rate / wind_speed)
    b = -0.25 * math.log10(max(alpha, 0.01))
    log_d_over_l = -0.85 * math.log10(ratio) + b
    return scale * math.pow(10.0, log_d_over_l)


def _gaussian_cloud_volume(
    plume_rate: float,
    stack_height: float,
    wind_speed: float,
    stability: str,
    terrain: str,
    threshold: float,
    distance: float,
) -> float:
    """The class's Gaussian cloud-volume integral above the LFL."""
    if not math.isfinite(distance) or distance <= 1.0 or threshold <= 0.0:
        return 0.0
    steps = _CLOUD_STEPS
    start = 1.0
    dx = (distance - start) / steps
    volume = 0.0
    for index in range(steps):
        x = start + (index + 0.5) * dx
        centreline = concentration(
            plume_rate, stack_height, wind_speed, stability, terrain, x, 0.0, 0.0
        )
        if centreline <= threshold:
            continue
        factor = math.sqrt(2.0 * math.log(centreline / threshold))
        semi_width = _sigma_y(x, stability, terrain) * factor
        semi_height = _sigma_z(x, stability, terrain) * factor
        volume += 0.5 * math.pi * semi_width * semi_height * dx
    return volume


def _heavy_cloud_volume(distance: float) -> float:
    """The class's dense-gas cloud-volume estimate."""
    if not math.isfinite(distance) or distance <= 0.0:
        return 0.0
    width = max(1.0, 0.25 * distance)
    height = min(max(1.0, 0.05 * distance), 5.0)
    return 0.5 * distance * width * height


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)


def gas_dispersion(
    components: list[str],
    composition: list[float],
    mass_release_rate: Q,
    ambient_temperature: Q | None = None,
    atmospheric_pressure: Q | None = None,
    wind_speed: Q | None = None,
    pasquill_stability_class: str | None = None,
    surface_roughness: Q | None = None,
    offshore: bool | None = None,
    release_height: Q | None = None,
    model_selection: str | None = None,
    lower_flammable_limit: float | None = None,
    flammable_endpoint_fraction_of_lfl: float | None = None,
    toxic_component: str | None = None,
    toxic_threshold: float | None = None,
) -> GasDispersionResult:
    """The flammable and toxic endpoints of a continuous gas release at ambient conditions.

    Args:
        components: the released fluid's substances, by name.
        composition: the released fluid's overall mole fractions.
        mass_release_rate: the continuous release rate.
        ambient_temperature: the ambient temperature the clone is flashed at.
        atmospheric_pressure: the ambient pressure, for the flash and the air density.
        wind_speed: the wind speed at the 10 m reference.
        pasquill_stability_class: the Pasquill class, ``A``-``F``.
        surface_roughness: the boundary conditions' surface roughness.
        offshore: whether the release is offshore, which selects rural terrain.
        release_height: the plume's effective release height.
        model_selection: ``auto``, ``gaussian_plume`` or ``heavy_gas``.
        lower_flammable_limit: a manual LFL override by volume fraction.
        flammable_endpoint_fraction_of_lfl: the endpoint fraction, in ``(0, 1]``.
        toxic_component: the component whose toxic endpoint is screened.
        toxic_threshold: the toxic threshold as a **volume fraction**, so 100 ppm is
            ``1.0e-4``.

    Returns:
        The densities, fuel data, flammable distances and cloud volume, and the toxic
        distance where one was asked and found.

    Raises:
        InvalidInputError: for a non-positive release rate, wind speed or stated manual LFL,
            an unknown stability or model word, a composition that is not one entry per
            component and sums to one, an endpoint fraction outside ``(0, 1]``, or whatever
            the databank or the ambient flash refuses.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    rate = input_to_si(spec, "mass_release_rate", mass_release_rate)
    ambient_t = (
        _DEFAULT_AMBIENT_TEMPERATURE
        if ambient_temperature is None
        else input_to_si(spec, "ambient_temperature", ambient_temperature)
    )
    ambient_p = (
        _DEFAULT_ATMOSPHERIC_PRESSURE
        if atmospheric_pressure is None
        else input_to_si(spec, "atmospheric_pressure", atmospheric_pressure)
    )
    wind = (
        _DEFAULT_WIND_SPEED if wind_speed is None else input_to_si(spec, "wind_speed", wind_speed)
    )
    roughness = (
        _DEFAULT_SURFACE_ROUGHNESS
        if surface_roughness is None
        else input_to_si(spec, "surface_roughness", surface_roughness)
    )
    apply_checks(
        checks.on_input,
        {
            "mass_release_rate": rate,
            "ambient_temperature": ambient_t,
            "atmospheric_pressure": ambient_p,
            "wind_speed": wind,
            "lower_flammable_limit": lower_flammable_limit,
        }.get,
        warnings,
    )

    stability = _stable(pasquill_stability_class or _DEFAULT_STABILITY)
    selection = model_selection or "auto"
    if selection not in _SELECTIONS:
        raise InvalidInputError(
            "model_selection",
            f"the model is `auto`, `gaussian_plume` or `heavy_gas`, and {selection!r} is not",
        )
    endpoint_fraction = (
        _DEFAULT_ENDPOINT_FRACTION
        if flammable_endpoint_fraction_of_lfl is None
        else float(flammable_endpoint_fraction_of_lfl)
    )
    if not 0.0 < endpoint_fraction <= 1.0:
        raise InvalidInputError(
            "flammable_endpoint_fraction_of_lfl",
            f"the endpoint fraction must be in (0, 1], and {endpoint_fraction} is not",
        )
    stack = (
        _DEFAULT_RELEASE_HEIGHT
        if release_height is None
        else max(0.0, input_to_si(spec, "release_height", release_height))
    )
    terrain = "rural" if bool(offshore) or roughness < 0.5 else "urban"
    _terrain(terrain)

    mixture, _ideal_gas = _components.mixture_of(components, eos="pr")
    if len(composition) != len(mixture.components):
        raise InvalidInputError(
            "composition",
            f"a mixture of {len(mixture.components)} components needs that many mole "
            f"fractions, and {composition!r} has {len(composition)}",
        )
    check_composition(list(composition), len(mixture.components), "composition")

    flash = pt_flash(
        mixture, quantity(ambient_t, "K"), quantity(ambient_p, "Pa"), list(composition)
    )
    if flash.phase == Phase.ALL_LIQUID:
        phase_z = list(flash.x)
        phase_root = flash.z_liquid
    else:
        phase_z = list(flash.y)
        phase_root = flash.z_vapour
    molar_masses = [
        None if component.molar_mass is None else component.molar_mass.to("kg/mol").magnitude
        for component in mixture.components
    ]
    phase_molar_mass = molar_mass_of(mixture, phase_z)
    molar_volume = (
        pr_molar_volume(phase_root, quantity(ambient_t, "K"), quantity(ambient_p, "Pa"))
        .v.to("m**3/mol")
        .magnitude
    )
    shift = mixture.volume_shift(phase_z)
    source_density = (
        pr_mass_density(
            quantity(phase_molar_mass, "kg/mol"),
            quantity(molar_volume - shift, "m**3/mol"),
        )
        .rho.to("kg/m**3")
        .magnitude
    )
    air_density = ambient_p / (_AIR_GAS_CONSTANT * ambient_t)

    fuel_mole, fuel_mass, fuel_molar_mass, fuel_lfl = _fuel_data(
        components, phase_z, molar_masses, phase_molar_mass, lower_flammable_limit
    )
    selected = selection
    if selected == "auto":
        selected = (
            "heavy_gas" if source_density / air_density > _DENSE_GAS_RATIO else "gaussian_plume"
        )
    has_fuel = fuel_mole > 0.0 and fuel_mass > 0.0 and fuel_lfl is not None

    distance_endpoint: float | None = None
    distance_half: float | None = None
    distance_lfl: float | None = None
    cloud_volume = 0.0
    if has_fuel:
        lfl = fuel_lfl if fuel_lfl is not None else math.nan
        if selected == "heavy_gas":
            volumetric_rate = rate / max(source_density, 1.0e-9)
            lfl_ratio = lfl / fuel_mole
            distance_lfl = _dense_distance(
                volumetric_rate, source_density, air_density, wind, lfl_ratio
            )
            distance_half = _dense_distance(
                volumetric_rate, source_density, air_density, wind, 0.5 * lfl_ratio
            )
            distance_endpoint = _dense_distance(
                volumetric_rate, source_density, air_density, wind, endpoint_fraction * lfl_ratio
            )
            cloud_volume = 0.0 if distance_lfl is None else _heavy_cloud_volume(distance_lfl)
        else:
            fuel_rate = rate * fuel_mass
            lfl_kg = lfl * fuel_molar_mass * ambient_p / (_GAS_CONSTANT * ambient_t)
            distance_lfl = distance_to_threshold(fuel_rate, stack, wind, stability, terrain, lfl_kg)
            distance_half = distance_to_threshold(
                fuel_rate, stack, wind, stability, terrain, 0.5 * lfl_kg
            )
            distance_endpoint = distance_to_threshold(
                fuel_rate, stack, wind, stability, terrain, endpoint_fraction * lfl_kg
            )
            cloud_volume = (
                0.0
                if distance_lfl is None
                else _gaussian_cloud_volume(
                    fuel_rate, stack, wind, stability, terrain, lfl_kg, distance_lfl
                )
            )

    stated_toxic = "" if toxic_component is None else toxic_component
    asked_toxic = stated_toxic.strip() != "" and (
        toxic_threshold is not None and toxic_threshold > 0.0
    )
    toxic_distance: float | None = None
    reported_threshold: float | None = None
    reported_name = stated_toxic if asked_toxic else ""
    if asked_toxic:
        reported_threshold = float(toxic_threshold)  # type: ignore[arg-type]
        component_mole, component_mass, component_molar_mass = _component_data(
            components, phase_z, molar_masses, stated_toxic
        )
        if component_mole > 0.0 and component_mass > 0.0 and math.isfinite(component_molar_mass):
            if selected == "heavy_gas":
                volumetric_rate = rate / max(source_density, 1.0e-9)
                toxic_distance = _dense_distance(
                    volumetric_rate,
                    source_density,
                    air_density,
                    wind,
                    reported_threshold / component_mole,
                )
            else:
                toxic_rate = rate * component_mass
                threshold_kg = (
                    reported_threshold
                    * component_molar_mass
                    * ambient_p
                    / (_TOXIC_GAS_CONSTANT * ambient_t)
                )
                toxic_distance = distance_to_threshold(
                    toxic_rate, stack, wind, stability, terrain, threshold_kg
                )

    basis = "heavy_gas_britter_mcquaid" if selected == "heavy_gas" else "gaussian_plume_briggs"
    return GasDispersionResult(
        selected_model=selected,
        mass_release_rate=from_si(rate, "kg/s"),
        flammable_mass_release_rate=from_si(rate * fuel_mass, "kg/s"),
        source_density=from_si(source_density, "kg/m**3"),
        air_density=from_si(air_density, "kg/m**3"),
        fuel_mole_fraction=fuel_mole,
        fuel_mass_fraction=fuel_mass,
        lower_flammable_limit=None if fuel_lfl is None else fuel_lfl,
        flammable_endpoint_fraction=endpoint_fraction,
        distance_to_flammable_endpoint=(
            None if distance_endpoint is None else from_si(distance_endpoint, "m")
        ),
        distance_to_half_lfl=None if distance_half is None else from_si(distance_half, "m"),
        distance_to_lfl=None if distance_lfl is None else from_si(distance_lfl, "m"),
        flammable_cloud_volume=from_si(cloud_volume, "m**3"),
        toxic_component=reported_name,
        toxic_threshold=reported_threshold,
        toxic_distance=None if toxic_distance is None else from_si(toxic_distance, "m"),
        wind_speed=from_si(wind, "m/s"),
        stability_class=stability,
        screening_basis=basis,
        has_flammable_cloud=distance_lfl is not None and distance_lfl > 0.0,
        has_toxic_endpoint=toxic_distance is not None and toxic_distance > 0.0,
        warnings=tuple(warnings),
    )
