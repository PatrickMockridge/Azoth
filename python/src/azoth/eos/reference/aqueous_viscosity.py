"""``eos.aqueous_viscosity`` - the liquid viscosity NeqSim gives an aqueous phase.

Spec: ``specs/models/eos/aqueous_viscosity.toml``

The Python twin of ``crates/azoth-eos/src/aqueous_viscosity.rs``, written to mirror it.

# The phase type is the whole point

``getPhysicalProperties()`` dispatches on the phase's type: ``PhaseType.AQUEOUS`` takes
``WaterPhysicalProperties``, whose viscosity model is the liquid ``Viscosity`` class - the
``polynom`` one - while a gas or a hydrocarbon liquid takes ``PFCTViscosityMethodHeavyOil``,
which ``eos.viscosity`` ports. So this is not "the liquid viscosity"; it is the one an
aqueous phase reaches.

# The mixing rule's interaction term is zero

``PhysicalPropertyMixingRule.initMixingRules`` fills ``Gij`` with an inner loop that starts
at ``k = l`` and breaks on ``k == l``, so its body never runs and every phase reports a zero
matrix - the probe prints the whole matrix for a water/methanol mixture whose ``gijvisc``
column is not empty. Grunberg-Nissan therefore reduces to ``exp(sum_i w_i ln mu_i)`` over
**mass** fractions.

# An ion is refused

``na+`` and ``cl-`` carry the same four ``LIQVISC`` numbers as methanol in NeqSim's own
table, and a brine computed from them comes out *less* viscous than pure water - the wrong
direction. The refusal names the ion and whose numbers its row carries.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import AqueousViscosityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Component, Mixture

#: The viscosity a component takes above its critical temperature, in cP.
ABOVE_CRITICAL_CP = 0.5
#: And the one it takes when its row names no model at all: NeqSim's `else` branch.
NO_MODEL_CP = 0.7


def aqueous_viscosity(mixture: Mixture, T: Q, P: Q, z: list[float]) -> AqueousViscosityResult:
    """The liquid viscosity of an aqueous phase.

    Args:
        mixture: the phase's components.
        T: absolute temperature of the state.
        P: absolute pressure of the state, which enters only through the correction.
        z: the phase's composition in mole fractions. The class reads it as *mass*
            fractions, which is where the mixing rule's weights come from.

    Returns:
        The phase's dynamic viscosity.

    Raises:
        InvalidInputError: if ``z`` is the wrong length, or if any component is an ion.
        PropertyUnavailableError: if a component carries no molar mass.
        OutOfRangeError: for a temperature or a pressure the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> from azoth.eos import mixture_of
        >>> r = aqueous_viscosity(
        ...     mixture_of(["water"])[0], q(300.0, "K"), q(1.0, "bar"), [1.0]
        ... )
        >>> round(r.viscosity.to("Pa*s").magnitude, 10)
        0.0008547585
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t = input_to_si(spec, "T", T)
    p = input_to_si(spec, "P", P)

    apply_checks(checks.on_input, {"T": t, "P": p}.get, warnings)

    components = mixture.components
    if len(z) != len(components):
        raise InvalidInputError(
            "z",
            f"a mixture of {len(components)} components needs {len(components)} mole "
            f"fractions, got {len(z)}",
        )

    # Imported here rather than at module scope: `components` imports this package's own
    # `__init__`, so a top-level import of it is a cycle.
    from azoth.eos.components import ION

    ions = [
        component_name(mixture, index)
        for index, component in enumerate(components)
        if component.component_class == ION
    ]
    if ions:
        raise InvalidInputError(
            "components",
            f"{', '.join(ions)} is an ion, and NeqSim's `LIQVISC` row for an ion carries the "
            f"same four numbers as methanol's - so a brine computed from it comes out less "
            f"viscous than pure water, which is the wrong direction. The class that would "
            f"close this is a salt-aware viscosity, not this correlation",
        )

    masses: list[float] = []
    for component in components:
        if component.molar_mass is None:
            raise PropertyUnavailableError(
                "component",
                "molar mass",
                "a card-added component needs its own molar mass: the mixing rule weights by "
                "mass fraction",
            )
        masses.append(component.molar_mass.to("kg/mol").magnitude)

    weights = [zi * mass for zi, mass in zip(z, masses, strict=True)]
    total_mass = sum(weights)
    if total_mass <= 0.0:
        raise InvalidInputError(
            "z", "a phase whose mass is zero has no mass fractions to weight by"
        )

    weighted = 0.0
    for index, component in enumerate(components):
        weight = weights[index] / total_mass
        if weight <= 0.0:
            # A zero weight contributes nothing, and `ln(0)` would poison the sum.
            continue
        weighted += weight * math.log(_pure_viscosity(component, t, p))

    return AqueousViscosityResult(
        viscosity=from_si(math.exp(weighted) * 1.0e-3, "Pa*s"),
        warnings=tuple(warnings),
    )


def component_name(mixture: Mixture, index: int) -> str:
    """The component's name, or its index where the mixture carries none."""
    if mixture.names is not None and index < len(mixture.names):
        return mixture.names[index]
    return f"the component at index {index}"


def _pure_viscosity(component: Component, temperature: float, pressure: float) -> float:
    """One component's pure-liquid viscosity in cP, with NeqSim's correction applied."""
    tc = component.Tc.to("K").magnitude
    pc = component.Pc.to("Pa").magnitude
    l1, l2, l3, l4 = component.liqvisc

    if temperature > tc:
        uncorrected = ABOVE_CRITICAL_CP
    elif component.liqvisc_model == 1:
        uncorrected = l1 * temperature**l2
    elif component.liqvisc_model == 2:
        uncorrected = math.exp(l1 + l2 / temperature)
    elif component.liqvisc_model == 3:
        uncorrected = math.exp(l1 + l2 / temperature + l3 * temperature + l4 * temperature**2)
    elif component.liqvisc_model == 4:
        uncorrected = 10.0 ** (l1 * (1.0 / temperature - 1.0 / l2))
    else:
        uncorrected = NO_MODEL_CP

    corrected = _pressure_correction(temperature, pressure, tc, pc, component.omega)
    return uncorrected * (corrected + 1.0) / 2.0


def _pressure_correction(
    temperature: float, pressure: float, tc: float, pc: float, omega: float
) -> float:
    """``getViscosityPressureCorrection``: NeqSim's own four-coefficient form."""
    reduced_t = temperature / tc
    if reduced_t > 1.0:
        return 1.0
    delta_pr = pressure / pc
    # `math.pow` rather than `**`: a float raised to a float exponent is `Any` to a type
    # checker, because Python's `float.__pow__` can return a complex.
    a = 0.9991 - (4.674e-4 / (1.0523 * math.pow(reduced_t, -0.03877) - 1.0513))
    d = (0.3257 / math.pow(1.0039 - reduced_t**2.573, 0.2906)) - 0.2086
    c = (
        -0.07921
        + 2.1616 * reduced_t
        - 13.4040 * reduced_t**2
        + 44.1706 * reduced_t**3
        - 84.8291 * reduced_t**4
        + 96.1209 * reduced_t**5
        - 59.8127 * reduced_t**6
        + 15.6719 * reduced_t**7
    )
    numerator = 1.0 + d * math.pow(delta_pr / 2.118, a)
    denominator = 1.0 + c * omega * delta_pr
    return numerator / denominator


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("eos.aqueous_viscosity")
