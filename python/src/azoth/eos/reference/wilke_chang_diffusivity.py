"""``eos.wilke_chang_diffusivity`` - the liquid binary diffusivity from the
Wilke-Chang correlation.

Spec: ``specs/calcs/eos/wilke_chang_diffusivity.toml``, which records the
correlation, the association parameter and the clamps.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WilkeChangDiffusivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.wilke_chang_diffusivity"


def wilke_chang_diffusivity(phi: float, M: Q, T: Q, eta: Q, VA: Q) -> WilkeChangDiffusivityResult:
    """The binary diffusion coefficient at infinite dilution, from Wilke-Chang.

    ``phi`` is the solvent association parameter; ``VA`` is the solute molar volume
    at the normal boiling point and ``eta`` the solvent viscosity. ``VA`` and
    ``eta`` are clamped to NeqSim's `[20, 600]` cm**3/mol and `[0.01, 500]` cP
    before the correlation.

    Args:
        phi: solvent association parameter; see ``components.wilke_chang_phi``.
        M: solvent molar mass.
        T: absolute temperature.
        eta: solvent dynamic viscosity at ``T``.
        VA: solute molar volume at its normal boiling point.

    Returns:
        The binary diffusion coefficient, in ``m**2/s``.

    Raises:
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> m = q(0.018015, "kg/mol")
        >>> va = q(4.0203262233375156e-5, "m**3/mol")
        >>> eta = q(8.915447896200597e-4, "Pa*s")
        >>> r = wilke_chang_diffusivity(2.26, m, q(298.15, "K"), eta, va)
        >>> round(r.d.magnitude, 16)
        1.7212261801806e-09
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "phi": phi,
        "M": input_to_si(spec, "M", M),
        "T": input_to_si(spec, "T", T),
        "eta": input_to_si(spec, "eta", eta),
        "VA": input_to_si(spec, "VA", VA),
    }
    apply_checks(checks.on_input, values.get, warnings)

    m_g = values["M"] * 1000.0
    va_cm3 = min(600.0, max(20.0, values["VA"] * 1.0e6))
    eta_cp = min(500.0, max(0.01, values["eta"] * 1000.0))

    d_cm2s = 7.4e-8 * (phi * m_g) ** 0.5 * values["T"] / (eta_cp * va_cm3**0.6)
    d = d_cm2s * 1.0e-4

    apply_checks(checks.derived, lambda name: d if name == "d" else None, warnings)

    return WilkeChangDiffusivityResult(d=from_si(d, "m**2/s"), warnings=tuple(warnings))
