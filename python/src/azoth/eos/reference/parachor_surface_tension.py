"""``eos.parachor_surface_tension`` - the surface tension from the parachor
(Macleod-Sugden) correlation.

Spec: ``specs/calcs/eos/parachor_surface_tension.toml``, which records the pure
component form of NeqSim's ``ParachorSurfaceTension.calcPureComponentSurfaceTension``.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ParachorSurfaceTensionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.parachor_surface_tension"


def parachor_surface_tension(
    parachor: float, rho_l: Q, rho_v: Q, M: Q
) -> ParachorSurfaceTensionResult:
    """The surface tension of a pure component, from the Macleod-Sugden correlation.

    ``parachor`` is in the mixed unit ``(mN/m)**(1/4) * cm**3/mol`` NeqSim stores;
    the ``1e-6`` in the formula converts its ``cm**3/mol`` to ``m**3/mol``, and the
    ``1e-3`` recovers ``N/m`` from ``mN/m``.

    Args:
        parachor: the parachor parameter, in ``(mN/m)**(1/4) * cm**3/mol``.
        rho_l: the liquid mass density at the interface state.
        rho_v: the vapour mass density at the interface state.
        M: the molar mass.

    Returns:
        The surface tension, in ``N/m``.

    Raises:
        OutOfRangeError: if ``M`` or ``rho_l`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = parachor_surface_tension(77.3, q(422.0, "kg/m**3"),
        ...     q(1.82, "kg/m**3"), q(0.016043, "kg/mol"))
        >>> round(r.sigma.magnitude, 16)
        0.01680030432034039
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "parachor": parachor,
        "rho_l": input_to_si(spec, "rho_l", rho_l),
        "rho_v": input_to_si(spec, "rho_v", rho_v),
        "M": input_to_si(spec, "M", M),
    }
    apply_checks(checks.on_input, values.get, warnings)

    molar_density_l = values["rho_l"] / values["M"]
    molar_density_v = values["rho_v"] / values["M"]
    sigma = 1e-3 * (values["parachor"] * 1e-6 * (molar_density_l - molar_density_v)) ** 4

    apply_checks(checks.derived, lambda name: sigma if name == "sigma" else None, warnings)

    return ParachorSurfaceTensionResult(sigma=from_si(sigma, "N/m"), warnings=tuple(warnings))
