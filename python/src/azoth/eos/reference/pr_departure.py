"""``eos.pr_departure`` - the Peng-Robinson fugacity coefficient and departures.

```text
psi      = -kappa*sqrt(Tr) / (1 + kappa*(1 - sqrt(Tr)))
I        = ln((z + (1 + sqrt(2))*B) / (z + (1 - sqrt(2))*B))
ln_phi   = z - 1 - ln(z - B) - C*I
h_dep_rt = (z - 1) + C*(psi - 1)*I
s_dep_r  = ln(z - B) + C*psi*I          where C = A / (2*sqrt(2)*B)
```

Spec: ``specs/calcs/eos/pr_departure.yaml``, which carries the provenance, the
derivation of the Gibbs identity ``h_dep_rt - s_dep_r = ln_phi``, and the assumptions
about ``z``, ``kappa`` and ``Tr`` this calc cannot check.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrDepartureResult
from azoth.core.warnings import Warning

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

    # The logarithmic derivative of the alpha function. Hoisted rather than
    # recomputed for each use so the two languages cannot evaluate it twice in
    # different orders.
    sqrt_tr = Tr**0.5
    psi = -kappa * sqrt_tr / (1.0 + kappa * (1.0 - sqrt_tr))

    sqrt_2 = math.sqrt(2.0)
    i_term = math.log((z + (1.0 + sqrt_2) * b_reduced) / (z + (1.0 - sqrt_2) * b_reduced))
    coefficient = a_reduced / (2.0 * sqrt_2 * b_reduced)
    ln_z_minus_b = math.log(z - b_reduced)

    ln_phi = z - 1.0 - ln_z_minus_b - coefficient * i_term
    h_dep_rt = (z - 1.0) + coefficient * (psi - 1.0) * i_term
    s_dep_r = ln_z_minus_b + coefficient * psi * i_term

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
        warnings=tuple(warnings),
    )
