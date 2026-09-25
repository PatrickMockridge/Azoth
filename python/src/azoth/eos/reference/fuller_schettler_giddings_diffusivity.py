"""``eos.fuller_schettler_giddings_diffusivity`` - the gas binary diffusivity from the
Fuller-Schettler-Giddings correlation.

Spec: ``specs/calcs/eos/fuller_schettler_giddings_diffusivity.toml``, which records the
correlation, the diffusion-volume ladder and the two ways the class states its pressure unit.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FullerSchettlerGiddingsDiffusivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.fuller_schettler_giddings_diffusivity"


def fuller_schettler_giddings_diffusivity(
    MA: Q, MB: Q, VA: Q, VB: Q, T: Q, P: Q
) -> FullerSchettlerGiddingsDiffusivityResult:
    """The binary diffusion coefficient of a gas pair, from Fuller-Schettler-Giddings.

    ``VA`` and ``VB`` are *diffusion volumes*, not molar volumes: the caller resolves
    each from the component's name with ``azoth.eos.components.fuller_diffusion_volume``,
    which is the class's own ladder (the special-molecule table, then ``0.285*Vc``, then
    ``max(10, 0.95*M)``).

    **The pressure the correlation takes is in bar** and the constant is ``1.013e-3``, which
    is what ``calcBinaryDiffusionCoefficient`` uses; the class's javadoc writes ``1.013e-2``
    with ``P`` in atm, ten times the code's constant and inconsistent with its own unit.

    Args:
        MA: molar mass of the first component.
        MB: molar mass of the second.
        VA: diffusion volume of the first, in ``m**3/mol``.
        VB: diffusion volume of the second.
        T: absolute temperature.
        P: absolute pressure, which the correlation takes in bar.

    Returns:
        The pair's binary diffusion coefficient, in ``m**2/s``.

    Raises:
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = fuller_schettler_giddings_diffusivity(
        ...     q(0.016043, "kg/mol"), q(0.0280135, "kg/mol"),
        ...     q(2.514e-5, "m**3/mol"), q(1.85e-5, "m**3/mol"),
        ...     q(298.15, "K"), q(101325.0, "Pa"))
        >>> round(r.d.magnitude, 18)
        2.155058583907974e-05
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "MA": input_to_si(spec, "MA", MA),
        "MB": input_to_si(spec, "MB", MB),
        "VA": input_to_si(spec, "VA", VA),
        "VB": input_to_si(spec, "VB", VB),
        "T": input_to_si(spec, "T", T),
        "P": input_to_si(spec, "P", P),
    }
    apply_checks(checks.on_input, values.get, warnings)

    # The published constants are tuned to g/mol, cm**3/mol, bar and cm**2/s.
    ma_g = values["MA"] * 1000.0
    mb_g = values["MB"] * 1000.0
    va_cm3 = values["VA"] * 1.0e6
    vb_cm3 = values["VB"] * 1.0e6
    p_bar = values["P"] * 1.0e-5

    sigma_v = va_cm3 ** (1.0 / 3.0) + vb_cm3 ** (1.0 / 3.0)
    d_cm2s = (
        1.013e-3
        * values["T"] ** 1.75
        * (1.0 / ma_g + 1.0 / mb_g) ** 0.5
        / (p_bar * sigma_v * sigma_v)
    )
    d = d_cm2s * 1.0e-4

    apply_checks(checks.derived, lambda name: d if name == "d" else None, warnings)

    return FullerSchettlerGiddingsDiffusivityResult(
        d=from_si(d, "m**2/s"), warnings=tuple(warnings)
    )
