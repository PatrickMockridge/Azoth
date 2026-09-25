"""``eos.parachor_mixture_surface_tension`` - the interface surface tension from the parachor.

Spec: ``specs/models/eos/parachor_mixture_surface_tension.toml``

The Python twin of ``crates/azoth-eos/src/parachor_mixture_surface_tension.rs``, written to
mirror it line for line.

# The pair's form, not the pure component's

``eos.parachor_surface_tension`` ports ``calcPureComponentSurfaceTension``, where the mole
fractions divide out. The interface form - ``ParachorSurfaceTension.calcSurfaceTension``,
reached through ``InterfaceProperties.getSurfaceTension`` - sums a per-component
**molar density** difference:

```text
sigma = 1e-3 * [ sum_i P_i * 1e-6 * (rho_l/M_l * x_l_i - rho_g/M_g * x_g_i) ]**4
```

so the gas's contribution is the negative one, which is the order
``InterfaceProperties.getSurfaceTension`` passes its two phases in.

# What the class answers where this refuses

The class catches the arithmetic exception and returns ``0.0``, and returns ``0.0`` for a
system with fewer than two phases besides. A silent zero reads as a tension of zero rather
than as a failure, so a non-positive density or molar mass is refused by name here.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ParachorMixtureSurfaceTensionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning


def parachor_mixture_surface_tension(
    parachors: Sequence[float],
    rho_gas: Q,
    M_gas: Q,
    x_gas: Sequence[float],
    rho_liquid: Q,
    M_liquid: Q,
    x_liquid: Sequence[float],
) -> ParachorMixtureSurfaceTensionResult:
    """The surface tension of the interface between a gas and a liquid.

    Args:
        parachors: each component's ``PARACHOR``, in NeqSim's mixed unit
            ``(mN/m)**(1/4) * cm**3/mol``.
        rho_gas: the gas phase's mass density.
        M_gas: the gas phase's molar mass.
        x_gas: the gas phase's mole fractions, in ``parachors``' order.
        rho_liquid: the liquid phase's mass density.
        M_liquid: the liquid phase's molar mass.
        x_liquid: the liquid phase's mole fractions.

    Returns:
        The interface's surface tension, in ``N/m``.

    Raises:
        InvalidInputError: if the three vectors differ in length.
        OutOfRangeError: if a density or a molar mass is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = parachor_mixture_surface_tension(
        ...     [77.3, 191.7],
        ...     q(19.938328330315997, "kg/m**3"), q(0.022959902730870334, "kg/mol"),
        ...     [0.8356249351028912, 0.1643750648971088],
        ...     q(549.3851668858096, "kg/m**3"), q(0.054107691623917306, "kg/mol"),
        ...     [0.09542082642782022, 0.9045791735721797])
        >>> round(r.sigma.magnitude, 18)
        0.009424889197268286
    """
    from azoth._models_gen import model

    spec = model("eos.parachor_mixture_surface_tension")
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(parachors)
    if len(x_gas) != count or len(x_liquid) != count:
        raise InvalidInputError(
            "parachors",
            f"{count} parachor(s), {len(x_gas)} gas fraction(s) and {len(x_liquid)} "
            f"liquid fraction(s): the sum is over the components, so the three vectors "
            f"are one entry each",
        )

    rho_g = input_to_si(spec, "rho_gas", rho_gas)
    m_g = input_to_si(spec, "M_gas", M_gas)
    rho_l = input_to_si(spec, "rho_liquid", rho_liquid)
    m_l = input_to_si(spec, "M_liquid", M_liquid)

    apply_checks(
        checks.on_input,
        {
            "rho_gas": rho_g,
            "M_gas": m_g,
            "rho_liquid": rho_l,
            "M_liquid": m_l,
        }.get,
        warnings,
    )

    # `rho/M` per phase, which is the molar density each component's fraction weights.
    molar_density_gas = rho_g / m_g
    molar_density_liquid = rho_l / m_l

    total = 0.0
    for i, parachor in enumerate(parachors):
        total += (
            float(parachor)
            * 1.0e-6
            * (molar_density_liquid * float(x_liquid[i]) - molar_density_gas * float(x_gas[i]))
        )
    sigma = 1.0e-3 * total**4

    apply_checks(checks.derived, lambda name: sigma if name == "sigma" else None, warnings)

    return ParachorMixtureSurfaceTensionResult(
        sigma=from_si(sigma, "N/m"), warnings=tuple(warnings)
    )
