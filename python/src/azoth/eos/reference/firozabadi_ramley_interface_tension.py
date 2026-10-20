"""``eos.firozabadi_ramley_interface_tension`` - the interface tension from the
Firoozabadi-Ramey corresponding-states correlation.

Spec: ``specs/calcs/eos/firozabadi_ramley_interface_tension.toml``, which records
NeqSim's ``FirozabadiRamleyInterfaceTension.calcSurfaceTension``.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FirozabadiRamleyInterfaceTensionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.firozabadi_ramley_interface_tension"


def firozabadi_ramley_interface_tension(
    T: Q, pseudo_critical_temperature: Q, rho_1: Q, rho_2: Q
) -> FirozabadiRamleyInterfaceTensionResult:
    """The interface tension between two phases from their mass densities.

    ``pseudo_critical_temperature`` is ``sum(x_i Tc_i)`` over the **first**
    phase's composition; the order of the two densities is part of the input.

    Args:
        T: absolute temperature at the interface.
        pseudo_critical_temperature: the first phase's pseudo-critical temperature.
        rho_1: mass density of the first phase.
        rho_2: mass density of the second phase.

    Returns:
        The interfacial tension, in ``N/m``.

    Raises:
        OutOfRangeError: if ``T``, the pseudo-critical temperature, or either
            density is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = firozabadi_ramley_interface_tension(
        ...     q(300.0, "K"), q(540.2432308699056, "K"),
        ...     q(690.6031246582897, "kg/m**3"),
        ...     q(1000.0024520207028, "kg/m**3"))
        >>> round(r.sigma.magnitude, 16)
        0.0508336185936405
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "T": input_to_si(spec, "T", T),
        "pseudo_critical_temperature": input_to_si(
            spec, "pseudo_critical_temperature", pseudo_critical_temperature
        ),
        "rho_1": input_to_si(spec, "rho_1", rho_1),
        "rho_2": input_to_si(spec, "rho_2", rho_2),
    }
    apply_checks(checks.on_input, values.get, warnings)

    delta_over_1000 = abs(values["rho_2"] - values["rho_1"]) / 1000.0
    if delta_over_1000 < 0.2:
        a1, b1 = 2.2062, -0.94716
    elif delta_over_1000 < 0.5:
        a1, b1 = 2.915, -0.76852
    else:
        a1, b1 = 3.3858, -0.62590
    reduced_temperature = values["T"] / values["pseudo_critical_temperature"]
    temp1 = a1 * delta_over_1000 ** (b1 + 1.0) / reduced_temperature**0.3125
    sigma = temp1**4 / 1000.0

    apply_checks(checks.derived, lambda name: sigma if name == "sigma" else None, warnings)

    return FirozabadiRamleyInterfaceTensionResult(
        sigma=from_si(sigma, "N/m"), warnings=tuple(warnings)
    )
