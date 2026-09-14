"""``eos.pr_departure`` - the Peng-Robinson fugacity coefficient and departures.

```text
psi      = -kappa*sqrt(Tr) / (1 + kappa*(1 - sqrt(Tr)))
I        = ln((z + (1 + sqrt(2))*B) / (z + (1 - sqrt(2))*B))
ln_phi   = z - 1 - ln(z - B) - C*I
h_dep_rt = (z - 1) + C*(psi - 1)*I
s_dep_r  = ln(z - B) + C*psi*I          where C = A / (2*sqrt(2)*B)
```

Spec: ``specs/calcs/eos/pr_departure.toml``, which carries the provenance, the
derivation of the Gibbs identity ``h_dep_rt - s_dep_r = ln_phi``, and the assumptions
about ``z``, ``kappa`` and ``Tr`` this calc cannot check.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrDepartureResult
from azoth.core.warnings import Warning
from azoth.eos.reference.alpha_term import Soave
from azoth.eos.reference.cubic import PR

CALC_ID = "eos.pr_departure"


def pr_departure(
    a_reduced: float,
    b_reduced: float,
    z: float,
    kappa: float,
    Tr: float,
) -> PrDepartureResult:
    """The Peng-Robinson fugacity coefficient and departure functions, for one state.

    Args:
        a_reduced: the cubic's attraction parameter ``A``, from
            :func:`azoth.eos.pr_alpha_ab`.
        b_reduced: the cubic's repulsion parameter ``B``, likewise.
        z: the compressibility factor, from :func:`azoth.eos.pr_z_factor`. Not
            checked against the cubic - see the module documentation.
        kappa: the alpha-function coefficient, from :func:`azoth.eos.pr_kappa` or
            :func:`azoth.eos.prsv_kappa`. Either serves.
        Tr: reduced temperature.

    Returns:
        ``ln_phi``, and the departure enthalpy and entropy made dimensionless as
        ``h_dep_rt`` and ``s_dep_r``. The multiplication by ``R`` and ``T`` happens
        where those live - the model layer - not here, so this namespace stays
        unit-free throughout.

    Raises:
        OutOfRangeError: if ``b_reduced <= 0`` (it is a divisor) or if
            ``z <= b_reduced`` (which makes ``ln(z - B)`` the logarithm of a
            negative number).

    Example:
        >>> d = pr_departure(
        ...     0.20206500174625697, 0.02431127309496514, 0.7907789662973796, 0.60282728832, 0.8
        ... )
        >>> round(d.ln_phi, 12)
        -0.191310556843
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
        "kappa": kappa,
        "Tr": Tr,
    }

    apply_checks(checks.on_input, values.get, warnings)

    # The logarithmic derivative of the alpha function, from the same Soave term the
    # mixture layer uses, so the two languages cannot evaluate it twice in different
    # orders or different places.
    term = Soave(kappa=kappa)
    psi = term.psi(Tr)

    i_term = PR.i_term(z, b_reduced)
    coefficient = PR.coefficient(a_reduced, b_reduced)
    ln_z_minus_b = math.log(z - b_reduced)

    ln_phi = z - 1.0 - ln_z_minus_b - coefficient * i_term
    h_dep_rt = (z - 1.0) + coefficient * (psi - 1.0) * i_term
    s_dep_r = ln_z_minus_b + coefficient * psi * i_term

    # The heat-capacity departure, from ``Cp^R/R = y + T*(dy/dT)_P`` with
    # ``y = h_dep_rt``. Every derivative below is already multiplied by ``T`` - that
    # is what the ``t_`` prefix means - so the ``1/T`` each of them carries cancels
    # and no term needs the absolute temperature, only ``Tr``.
    t_da = a_reduced * (psi - 2.0)
    t_db = -b_reduced
    t_dpsi = term.psi_t(Tr)
    t_dc = coefficient * (psi - 1.0)

    # ``z`` is a root of ``F(z, T) = 0``, so ``dz/dT = -(dF/dT)/(dF/dz)``.
    d_f_dz = PR.df_dz(z, a_reduced, b_reduced)
    t_dfdt = PR.t_dfdt(z, a_reduced, b_reduced, t_da, t_db)
    t_dz = -t_dfdt / d_f_dz

    n_plus = z + PR.delta1 * b_reduced
    n_minus = z + PR.delta2 * b_reduced
    t_di = (t_dz + PR.delta1 * t_db) / n_plus - (t_dz + PR.delta2 * t_db) / n_minus

    cp_dep_r = (
        h_dep_rt
        + t_dz
        + t_dc * (psi - 1.0) * i_term
        + coefficient * t_dpsi * i_term
        + coefficient * (psi - 1.0) * t_di
    )

    def derived(quantity: str) -> float | None:
        # The bound that matters is on the difference, not on `z`: `z = B` is the
        # zero-volume limit and nothing about `z` alone says where it is.
        if quantity == "z_minus_b_reduced":
            return z - b_reduced
        return None

    apply_checks(checks.derived, derived, warnings)

    return PrDepartureResult(
        ln_phi=ln_phi,
        h_dep_rt=h_dep_rt,
        s_dep_r=s_dep_r,
        cp_dep_r=cp_dep_r,
        warnings=tuple(warnings),
    )
