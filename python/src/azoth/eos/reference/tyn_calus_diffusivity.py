"""``eos.tyn_calus_diffusivity`` - the liquid binary diffusivity from the Tyn-Calus
correlation.

Spec: ``specs/calcs/eos/tyn_calus_diffusivity.toml``, which records the correlation
and the molar-volume and viscosity clamps.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TynCalusDiffusivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.tyn_calus_diffusivity"


def tyn_calus_diffusivity(VA: Q, VB: Q, T: Q, eta: Q) -> TynCalusDiffusivityResult:
    """The binary diffusion coefficient at infinite dilution, from Tyn-Calus.

    ``VA`` and ``VB`` are the solute and solvent liquid molar volumes at the normal
    boiling point; ``eta`` is the solvent viscosity. All three are clamped to
    NeqSim's `[20, 600]` cm**3/mol and `[0.01, 500]` cP before the correlation.

    Args:
        VA: solute molar volume at its normal boiling point.
        VB: solvent molar volume at its normal boiling point.
        T: absolute temperature.
        eta: solvent dynamic viscosity at ``T``.

    Returns:
        The binary diffusion coefficient, in ``m**2/s``.

    Raises:
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> va = q(8.816478555304741e-5, "m**3/mol")
        >>> vb = q(1.0578760045924226e-4, "m**3/mol")
        >>> eta = q(9.163501315189954e-4, "Pa*s")
        >>> r = tyn_calus_diffusivity(va, vb, q(298.15, "K"), eta)
        >>> round(r.d.magnitude, 16)
        1.4502274002446e-09
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "VA": input_to_si(spec, "VA", VA),
        "VB": input_to_si(spec, "VB", VB),
        "T": input_to_si(spec, "T", T),
        "eta": input_to_si(spec, "eta", eta),
    }
    apply_checks(checks.on_input, values.get, warnings)

    va_cm3 = min(600.0, max(20.0, values["VA"] * 1.0e6))
    vb_cm3 = min(600.0, max(20.0, values["VB"] * 1.0e6))
    eta_cp = min(500.0, max(0.01, values["eta"] * 1000.0))

    d_cm2s = 8.93e-8 * vb_cm3**0.267 * values["T"] / (eta_cp * va_cm3**0.433)
    d = d_cm2s * 1.0e-4

    apply_checks(checks.derived, lambda name: d if name == "d" else None, warnings)

    return TynCalusDiffusivityResult(d=from_si(d, "m**2/s"), warnings=tuple(warnings))
