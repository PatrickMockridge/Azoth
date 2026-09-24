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

from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FlareResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
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


def flare(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
) -> FlareResult:
    """A flare's steady state: the record through, and the two numbers beside it.

    Args:
        components: the gas's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.

    Returns:
        The product's record - the inlet's, unchanged - and the flare's duty and emission.

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
    apply_checks(checks.on_input, {"inlet_n": n}.get, warnings)

    states = _route(components, n, inlet_z, p, t)

    return FlareResult(
        product_n=from_si(states.product_n, "mol/s"),
        product_z=states.product_z,
        product_p=from_si(p, "Pa"),
        product_t=inlet_t,
        product_h=from_si(states.product_h, "J/mol"),
        heat_duty=from_si(states.heat_duty, "W"),
        co2_emission=from_si(states.co2_emission, "kg/s"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class FlareStates:
    """The outlet's record and the flare's two numbers, in SI."""

    product_n: float
    product_z: tuple[float, ...]
    product_h: float
    heat_duty: float
    co2_emission: float


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

    quality = quality_route(
        components, inlet_z, NORMAL_TEMPERATURE, SIXTY_FAHRENHEIT
    )
    moles_per_normal_cubic_metre = (
        REFERENCE_PRESSURE / (R_STANDARD * NORMAL_TEMPERATURE * quality.compression_factor)
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

    return FlareStates(
        product_n=inlet_n,
        product_z=tuple(inlet_z),
        product_h=float(product_h),
        heat_duty=heat_duty,
        co2_emission=carbon * CO2_MOLAR_MASS,
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.flare")
