"""``eos.chapman_enskog_diffusivity`` - the gas binary diffusivity from the Chapman-Enskog
theory.

Spec: ``specs/calcs/eos/chapman_enskog_diffusivity.toml``, which records the correlation, the
Neufeld collision integral and the two routes over the Lennard-Jones parameters.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ChapmanEnskogDiffusivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.chapman_enskog_diffusivity"

#: Neufeld, Janzen & Aziz (1972)'s eight coefficients of the collision integral.
_A, _B, _C, _D = 1.06036, 0.15610, 0.19300, 0.47635
_E, _F, _G, _H = 1.03587, 1.52996, 1.76474, 3.89411


def collision_integral(t: float) -> float:
    """``Omega_D(T/eps)``, Neufeld's fit - the pair's collision integral."""
    return (
        _A / math.pow(t, _B) + _C / math.exp(_D * t) + _E / math.exp(_F * t) + _G / math.exp(_H * t)
    )


def chapman_enskog_diffusivity(
    MA: Q, MB: Q, sigma: Q, eps: Q, T: Q, P: Q
) -> ChapmanEnskogDiffusivityResult:
    """The binary diffusivity of a gas pair, from the Chapman-Enskog theory.

    ``sigma`` and ``eps`` are the **pair's** Lennard-Jones parameters, combined with
    :func:`azoth.eos.components.lennard_jones_pair`.

    **This is a gas phase's default diffusivity.** NeqSim assigns the base ``Diffusivity``
    class to a gas, and the Fuller correlation is reached only by selecting a model; the
    rate-based column selects nothing, so its segment fluxes take this - over the database's
    Lennard-Jones parameters, which answer ``3.555e-5`` m**2/s for methane/nitrogen at 298 K
    where Poling's textbook parameters answer ``2.185e-5`` and the measurement is ``2.2e-5``.

    Args:
        MA: molar mass of the first component.
        MB: molar mass of the second.
        sigma: the pair's collision diameter, in ``m``.
        eps: the pair's energy parameter over Boltzmann's constant, in ``K``.
        T: absolute temperature.
        P: absolute pressure, which the correlation takes in bar.

    Returns:
        The pair's binary diffusion coefficient, in ``m**2/s``.

    Raises:
        OutOfRangeError: if ``T``, ``sigma`` or ``eps`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = chapman_enskog_diffusivity(
        ...     q(0.016043, "kg/mol"), q(0.0280135, "kg/mol"),
        ...     q(2.826253374e-10, "m"), q(140.43522570075095, "K"),
        ...     q(298.15, "K"), q(101325.0, "Pa"))
        >>> round(r.d.magnitude, 18)
        3.555070057344887e-05
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "MA": input_to_si(spec, "MA", MA),
        "MB": input_to_si(spec, "MB", MB),
        "sigma": input_to_si(spec, "sigma", sigma),
        "eps": input_to_si(spec, "eps", eps),
        "T": input_to_si(spec, "T", T),
        "P": input_to_si(spec, "P", P),
    }
    apply_checks(checks.on_input, values.get, warnings)

    # The published constants are tuned to g/mol, angstrom, bar and cm**2/s.
    pair_mass_g = 2.0 / (1.0 / (values["MA"] * 1000.0) + 1.0 / (values["MB"] * 1000.0))
    sigma_angstrom = values["sigma"] * 1.0e10
    p_bar = values["P"] * 1.0e-5

    omega = collision_integral(values["T"] / values["eps"])
    d_cm2s = (
        0.00266
        * values["T"] ** 1.5
        / (p_bar * pair_mass_g**0.5 * sigma_angstrom * sigma_angstrom * omega)
    )
    d = d_cm2s * 1.0e-4

    apply_checks(checks.derived, lambda name: d if name == "d" else None, warnings)

    return ChapmanEnskogDiffusivityResult(d=from_si(d, "m**2/s"), warnings=tuple(warnings))
