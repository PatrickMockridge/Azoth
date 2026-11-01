"""``process.flare`` - the flare's kernel.

Spec: ``specs/models/process/flare.toml``

The Python twin of ``crates/azoth-process/src/models/flare.rs`` over
``crates/azoth-process/src/kernels/flare.rs``, written to mirror them.

# The machine is a pass-through with a report

``Flare.run`` clones the inlet's thermo system into the outlet and changes nothing in it, so
the record's five fields are the outlet's exactly. What the machine computes is its own two
numbers:

``heatDuty    = inStream.LCV() * inStream.getFlowRate("Sm3/sec")``
``co2Emission = (sum_i z_i * n * nC_i) * 44.01e-3``

**The report has a surface beside them, and `run` reads none of it.** The radiation model is a
point source at the flame height, its inverse a horizontal distance; the dispersion surrogate
is the tip's exit velocity and momentum flux over the inlet's system density; and the capacity
check divides the operating point by up to three design limits. The setters feed the first and
the third, and they clamp rather than refuse.

# The duty multiplies two different volumetric references

``LCV()`` is ``new Standard_ISO6976(fluid, 0, 15.55, "volume").getValue("InferiorCalorificValue")
* 1.0e3`` - **joules per normal cubic metre at 0 °C**. ``Sm3/sec`` is
``n * R * 288.15 / atm`` - **a volume at 15 °C**, ideal. The product is therefore an energy
density at one reference times a volumetric flow at another, and the duty it reports is
``288.15 / 273.15`` - 5.5 per cent - above the same gas measured consistently. That is what
the class computes and what this reproduces, and the spec's assumptions say so.

Two molar gas constants are used, and the class really uses both: the standard's density
reaches ``8.314510`` and the Sm3 conversion divides by ``8.3144621``.

# The carbon comes from the element table

``getElements().getNumberOfElements("C")`` reads the component's formula, not the standard's
table - so a gas ISO 6976 has no row for still gets its CO2 number, and only its duty is
refused.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FlareResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.process.reference import capacity
from azoth.reactions.reference import _tables
from azoth.standards.reference.iso6976 import _route as quality_route

#: `Standard_ISO6976`'s own molar gas constant, J/(mol*K).
R_STANDARD = 8.314510

#: `ThermodynamicConstantsInterface.R`, which `getFlowRate("Sm3/sec")` uses.
R_THERMO = 8.3144621

#: The standard's reference pressure, Pa.
REFERENCE_PRESSURE = 101_325.0

#: The standard's normal (0 °C) reference, K - the basis the calorific value is per.
NORMAL_TEMPERATURE = 273.15

#: NeqSim's `standardStateTemperature`, K - 15 °C, the basis `Sm3` is measured at.
STANDARD_STATE_TEMPERATURE = 288.15

#: The 60 °F reference `Stream.LCV()` states the energy at, K.
SIXTY_FAHRENHEIT = 288.7

#: Molar mass of CO2, kg/mol, as `Flare.run` hard-codes it.
CO2_MOLAR_MASS = 44.01e-3

#: `Flare.flameHeight`'s initialiser, m - the point source's height.
DEFAULT_FLAME_HEIGHT = 30.0

#: `Flare.radiantFraction`'s initialiser: the share of the duty radiated.
DEFAULT_RADIANT_FRACTION = 0.18

#: `Flare.tipDiameter`'s initialiser, m.
DEFAULT_TIP_DIAMETER = 0.3

#: `setTipDiameter`'s floor, m.
MIN_TIP_DIAMETER = 1.0e-4

#: The surrogate's density floor, kg/m3.
MIN_DENSITY = 1.0e-3

#: The surrogate's tip-area floor, m2.
MIN_TIP_AREA = 1.0e-6

#: `CapacityCheckResult.exceeds`: a utilisation above it is an overload.
OVERLOAD_TOLERANCE = 1.0e-6


def flare(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    flame_height: Q | None = None,
    radiant_fraction: float | None = None,
    tip_diameter: Q | None = None,
    ground_distance: Q | None = None,
    radiation_flux_threshold: Q | None = None,
    design_heat_duty_capacity: Q | None = None,
    design_mass_flow_capacity: Q | None = None,
    design_molar_flow_capacity: Q | None = None,
) -> FlareResult:
    """A flare's steady state: the record through, and the numbers beside it.

    Args:
        components: the gas's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        flame_height: the point source's height; absent is the class's 30 m, clamped at zero.
        radiant_fraction: the share of the duty radiated; absent is 0.18, clamped to ``[0, 1]``.
        tip_diameter: the tip's outlet diameter; absent is 0.3 m, clamped at 1e-4 m.
        ground_distance: the distance the radiation flux is read at; absent leaves it absent.
        radiation_flux_threshold: the flux the distance is solved for; absent leaves it absent.
        design_heat_duty_capacity: the design heat release; absent or non-positive leaves the
            heat utilisation absent.
        design_mass_flow_capacity: the design mass flow, as the heat's is.
        design_molar_flow_capacity: the design molar flow, as the heat's is.

    Returns:
        The product's record - the inlet's, unchanged - the flare's duty and emission, and the
        radiation, dispersion and capacity numbers the class reports beside them.

    Raises:
        InvalidInputError: where a component has no row in the standard's table or none in the
            element table.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = flare(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.9, 0.1],
        ...     q(1.01325, "bar"),
        ...     q(288.15, "K"),
        ... )
        >>> round(r.co2_emission.to("kg/s").magnitude, 6)
        0.057213
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    p = input_to_si(spec, "inlet_p", inlet_p)
    t = input_to_si(spec, "inlet_t", inlet_t)
    distance = (
        None if ground_distance is None else input_to_si(spec, "ground_distance", ground_distance)
    )
    threshold = (
        None
        if radiation_flux_threshold is None
        else input_to_si(spec, "radiation_flux_threshold", radiation_flux_threshold)
    )
    design_heat = (
        None
        if design_heat_duty_capacity is None
        else input_to_si(spec, "design_heat_duty_capacity", design_heat_duty_capacity)
    )
    design_mass = (
        None
        if design_mass_flow_capacity is None
        else input_to_si(spec, "design_mass_flow_capacity", design_mass_flow_capacity)
    )
    design_molar = (
        None
        if design_molar_flow_capacity is None
        else input_to_si(spec, "design_molar_flow_capacity", design_molar_flow_capacity)
    )
    apply_checks(
        checks.on_input,
        {
            "inlet_n": n,
            "ground_distance": distance,
            "radiation_flux_threshold": threshold,
        }.get,
        warnings,
    )

    states = _route(components, n, inlet_z, p, t)

    # **The setters clamp**, so the port clamps: `setFlameHeight` at zero, `setRadiantFraction`
    # to `[0, 1]` with a NaN leaving the field standing, and `setTipDiameter` at 1e-4 m.
    height = (
        DEFAULT_FLAME_HEIGHT
        if flame_height is None
        else max(0.0, input_to_si(spec, "flame_height", flame_height))
    )
    fraction = DEFAULT_RADIANT_FRACTION if radiant_fraction is None else float(radiant_fraction)
    if not math.isnan(fraction):
        fraction = min(max(fraction, 0.0), 1.0)
    else:
        fraction = DEFAULT_RADIANT_FRACTION
    tip = (
        DEFAULT_TIP_DIAMETER
        if tip_diameter is None
        else max(MIN_TIP_DIAMETER, input_to_si(spec, "tip_diameter", tip_diameter))
    )
    flux = (
        None
        if distance is None
        else _radiation_heat_flux(states.heat_duty, distance, height, fraction)
    )
    radius = (
        None
        if threshold is None
        else _radiation_distance(states.heat_duty, threshold, height, fraction)
    )
    velocity, momentum = _dispersion(states, tip)
    heat, mass, molar, overloaded = _capacity(states, design_heat, design_mass, design_molar)

    return FlareResult(
        product_n=from_si(states.product_n, "mol/s"),
        product_z=states.product_z,
        product_p=from_si(p, "Pa"),
        product_t=inlet_t,
        product_h=from_si(states.product_h, "J/mol"),
        heat_duty=from_si(states.heat_duty, "W"),
        co2_emission=from_si(states.co2_emission, "kg/s"),
        radiation_flux=None if flux is None else from_si(flux, "W/m**2"),
        radiation_distance=None if radius is None else from_si(radius, "m"),
        exit_velocity=from_si(velocity, "m/s"),
        momentum_flux=from_si(momentum, "Pa"),
        standard_volume_rate=from_si(states.standard_volumetric_flow, "m**3/s"),
        heat_utilization=heat,
        mass_utilization=mass,
        molar_utilization=molar,
        overloaded=overloaded,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class FlareStates:
    """The outlet's record, the flare's two numbers and the operating point, in SI."""

    product_n: float
    product_z: tuple[float, ...]
    product_h: float
    heat_duty: float
    co2_emission: float
    standard_volumetric_flow: float
    density: float
    mass_rate: float
    molar_rate: float


def _route(
    components: list[str], inlet_n: float, inlet_z: list[float], inlet_p: float, inlet_t: float
) -> FlareStates:
    """The flare's arithmetic, in SI, in the order it computes it.

    **One arithmetic, two consumers**: :func:`flare` builds the result from this, and it is
    the twin of ``crates/azoth-process/src/kernels/flare.rs``.

    The duty's two pieces come from :func:`azoth.standards.reference.iso6976._route` - the
    *route* rather than the model, because a model's warnings are statements about its own
    input surface and this one has no reference temperatures of its own to warn about.
    """
    if inlet_n < 0.0:
        raise ValueError(f"a molar flow cannot be negative, and this one is {inlet_n}")

    quality = quality_route(components, inlet_z, NORMAL_TEMPERATURE, SIXTY_FAHRENHEIT)
    moles_per_normal_cubic_metre = REFERENCE_PRESSURE / (
        R_STANDARD * NORMAL_TEMPERATURE * quality.compression_factor
    )
    volumetric_calorific_value = quality.inferior_calorific_value * moles_per_normal_cubic_metre
    standard_volumetric_flow = inlet_n * R_THERMO * STANDARD_STATE_TEMPERATURE / REFERENCE_PRESSURE
    heat_duty = volumetric_calorific_value * standard_volumetric_flow

    carbon = 0.0
    for name, fraction in zip(components, inlet_z, strict=True):
        composition = _tables.element_composition(name)
        if composition is None:
            raise InvalidInputError(
                "components",
                f"{name} is not in the element table, so the carbon it carries is unknown "
                f"and the CO2 emission cannot be computed",
            )
        atoms = sum(count for element, count in composition if element == "C")
        carbon += fraction * inlet_n * atoms

    # The inlet's own enthalpy, which the outlet carries: `Stream::from_pt`'s state.
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    product_h, _ = enthalpy_at(mixture, ideal_gas, inlet_t, inlet_p, inlet_z)

    # **The dispersion surrogate's density is the inlet's own system density**,
    # `getThermoSystem().getDensity`'s reading: the phases' volume-weighted mean, which answers on
    # a two-phase feed where a single phase's cubic volume would not.
    molar_mass = capacity._molar_mass(mixture, list(inlet_z))
    density = capacity.system_mass_density(components, inlet_t, inlet_p, list(inlet_z))

    return FlareStates(
        product_n=inlet_n,
        product_z=tuple(inlet_z),
        product_h=float(product_h),
        heat_duty=heat_duty,
        co2_emission=carbon * CO2_MOLAR_MASS,
        standard_volumetric_flow=standard_volumetric_flow,
        density=density,
        mass_rate=inlet_n * molar_mass,
        molar_rate=inlet_n,
    )


def _radiation_heat_flux(
    heat_duty: float, ground_distance: float, flame_height: float, radiant_fraction: float
) -> float:
    """`estimateRadiationHeatFlux`: the point source at the flame height over a sphere."""
    radial = math.sqrt(max(0.0, ground_distance * ground_distance + flame_height * flame_height))
    if radial < 1.0e-6:
        return 0.0
    return radiant_fraction * heat_duty / (4.0 * math.pi * radial * radial)


def _radiation_distance(
    heat_duty: float, threshold: float, flame_height: float, radiant_fraction: float
) -> float:
    """`radiationDistanceForFlux`: the inverse, zero where the contour is inside the sphere."""
    if threshold <= 0.0:
        return 0.0
    numerator = radiant_fraction * heat_duty
    if numerator <= 0.0:
        return 0.0
    radial_squared = numerator / (4.0 * math.pi * threshold)
    horizontal_squared = radial_squared - flame_height * flame_height
    if horizontal_squared <= 0.0:
        return 0.0
    return math.sqrt(horizontal_squared)


def _dispersion(states: FlareStates, tip_diameter: float) -> tuple[float, float]:
    """`getDispersionSurrogate`: the tip's exit velocity and momentum flux, in SI.

    **The operating point is the reference**, so the surrogate's ``massRate/referenceMassRate``
    scaling is one and its standard volume rate is the inlet's own ``Sm3/sec``.
    """
    density = max(MIN_DENSITY, states.density)
    area = max(MIN_TIP_AREA, math.pi * tip_diameter * tip_diameter / 4.0)
    velocity = max(0.0, states.mass_rate / (density * area))
    return velocity, density * velocity * velocity


def _capacity(
    states: FlareStates,
    design_heat: float | None,
    design_mass: float | None,
    design_molar: float | None,
) -> tuple[float | None, float | None, float | None, bool | None]:
    """`evaluateCapacity`: the three utilisations and the overload flag.

    **A utilisation exists only where its design capacity is finite and positive**, which is the
    class's own test; and the flag exists only where a utilisation does, because the class's
    `false` with no limits is a default rather than a check.
    """

    def utilization(actual: float, design: float | None) -> float | None:
        if design is None or not math.isfinite(design) or design <= 0.0:
            return None
        return actual / design

    heat = utilization(states.heat_duty, design_heat)
    mass = utilization(states.mass_rate, design_mass)
    molar = utilization(states.molar_rate, design_molar)
    present = [value for value in (heat, mass, molar) if value is not None]
    if not present:
        return heat, mass, molar, None
    return heat, mass, molar, any(value > 1.0 + OVERLOAD_TOLERANCE for value in present)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.flare")
