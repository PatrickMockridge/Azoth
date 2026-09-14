"""``eos.rk_departure`` - the Redlich-Kwong fugacity coefficient and departures.

```text
psi      = -1/2
I        = ln((z + B) / z)
ln_phi   = z - 1 - ln(z - B) - C*I
h_dep_rt = (z - 1) + C*(psi - 1)*I
s_dep_r  = ln(z - B) + C*psi*I          where C = A / B
```

Spec: ``specs/calcs/eos/rk_departure.toml``. The geometry is SRK's ``delta = (1, 0)``,
and the alpha derivative is the constant ``-1/2``.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import RkDepartureResult
from azoth.core.warnings import Warning
from azoth.eos.alpha_term import RkAlpha
from azoth.eos.cubic import RK

CALC_ID = "eos.rk_departure"


def rk_departure(a_reduced: float, b_reduced: float, z: float) -> RkDepartureResult:
    """The Redlich-Kwong fugacity coefficient and departure functions, for one state.

    Args:
        a_reduced: the cubic's attraction parameter ``A``.
        b_reduced: the cubic's repulsion parameter ``B``.
        z: the compressibility factor.

    Raises:
        OutOfRangeError: if ``b_reduced <= 0`` or ``z <= b_reduced``.

    Example:
        >>> d = rk_departure(0.1866943088347048, 0.02707510936404929, 0.8119920001727409)
        >>> round(d.ln_phi, 12)
        -0.172001805597
        >>> abs(d.h_dep_rt - d.s_dep_r - d.ln_phi) < 1e-15
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "a_reduced": a_reduced,
        "b_reduced": b_reduced,
        "z": z,
    }

    apply_checks(checks.on_input, values.get, warnings)

    # RK's logarithmic derivative is the constant -1/2 and its temperature
    # derivative is zero.
    term = RkAlpha()
    psi = term.psi(1.0)
    t_dpsi = term.psi_t(1.0)

    i_term = RK.i_term(z, b_reduced)
    coefficient = RK.coefficient(a_reduced, b_reduced)
    ln_z_minus_b = math.log(z - b_reduced)

    ln_phi = z - 1.0 - ln_z_minus_b - coefficient * i_term
    h_dep_rt = (z - 1.0) + coefficient * (psi - 1.0) * i_term
    s_dep_r = ln_z_minus_b + coefficient * psi * i_term

    t_da = a_reduced * (psi - 2.0)
    t_db = -b_reduced
    t_dc = coefficient * (psi - 1.0)

    d_f_dz = RK.df_dz(z, a_reduced, b_reduced)
    t_dfdt = RK.t_dfdt(z, a_reduced, b_reduced, t_da, t_db)
    t_dz = -t_dfdt / d_f_dz

    n_plus = z + RK.delta1 * b_reduced
    n_minus = z + RK.delta2 * b_reduced
    t_di = (t_dz + RK.delta1 * t_db) / n_plus - (t_dz + RK.delta2 * t_db) / n_minus

    cp_dep_r = (
        h_dep_rt
        + t_dz
        + t_dc * (psi - 1.0) * i_term
        + coefficient * t_dpsi * i_term
        + coefficient * (psi - 1.0) * t_di
    )

    def derived(quantity: str) -> float | None:
        if quantity == "z_minus_b_reduced":
            return z - b_reduced
        return None

    apply_checks(checks.derived, derived, warnings)

    return RkDepartureResult(
        ln_phi=ln_phi,
        h_dep_rt=h_dep_rt,
        s_dep_r=s_dep_r,
        cp_dep_r=cp_dep_r,
        warnings=tuple(warnings),
    )
