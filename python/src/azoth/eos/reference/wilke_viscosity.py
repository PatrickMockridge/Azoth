"""``eos.wilke_viscosity`` - the gas mixture viscosity from Wilke's rule over
Chung's pure-component viscosities.

Spec: ``specs/models/eos/wilke_viscosity.toml``, which records the mixing rule and
points the pure-component half at :func:`azoth.eos.chung_viscosity`. It is a
*direct* model: no iteration, so no algorithm block.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WilkeViscosityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.chung_viscosity import chung_viscosity

MODEL_ID = "eos.wilke_viscosity"


def wilke_viscosity(
    Tc: Sequence[Q],
    Vc: Sequence[Q],
    M: Sequence[Q],
    omega: Sequence[float],
    dipole: Sequence[float],
    kappa: Sequence[float],
    T: Q,
    V: Q,
    z: Sequence[float],
) -> WilkeViscosityResult:
    """The gas dynamic viscosity of a mixture, by Wilke's rule over Chung viscosities.

    Each pure-component viscosity is :func:`azoth.eos.chung_viscosity`'s at the
    *mixture* molar volume ``V``, then mixed: ``mu = sum_i z_i mu_i / sum_j z_j
    phi_ij`` with ``phi_ij = (1 + sqrt(mu_i/mu_j) (M_j/M_i)**0.25)**2 / sqrt(8 (1 +
    M_i/M_j))``.

    Args:
        Tc: critical temperature of each component.
        Vc: critical volume of each component.
        M: molar mass of each component.
        omega: acentric factor of each component.
        dipole: dipole moment of each component, in debye.
        kappa: viscosity correction factor of each component.
        T: absolute temperature.
        V: molar volume of the gas mixture at ``T``.
        z: mole fractions; non-negative and summing to one.

    Returns:
        The gas mixture dynamic viscosity, in ``Pa*s``.

    Raises:
        InvalidInputError: if the vectors disagree in length, or ``z`` is not a
            composition.
        OutOfRangeError: if ``T`` or ``V`` is not positive, or a component's
            ``Tc``/``Vc`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> tc = [q(190.56, "K"), q(369.83, "K")]
        >>> vc = [q(9.9e-5, "m**3/mol"), q(0.000203, "m**3/mol")]
        >>> m = [q(0.016043, "kg/mol"), q(0.044097, "kg/mol")]
        >>> omega = [0.0115, 0.1523]
        >>> dipole = [0.0, 0.0]
        >>> kappa = [0.0, 0.0]
        >>> t = q(300.0, "K")
        >>> v = q(2.2987856717900145e-3, "m**3/mol")
        >>> z = [0.5, 0.5]
        >>> r = wilke_viscosity(tc, vc, m, omega, dipole, kappa, t, v, z)
        >>> round(r.mu.magnitude, 16)
        9.5430743692361e-06
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "T": input_to_si(spec, "T", T),
        "V": input_to_si(spec, "V", V),
    }
    apply_checks(checks.on_input, values.get, warnings)

    tc = [input_to_si(spec, "Tc", value) for value in Tc]
    vc = [input_to_si(spec, "Vc", value) for value in Vc]
    m = [input_to_si(spec, "M", value) for value in M]
    omega = list(omega)
    dipole = list(dipole)
    kappa = list(kappa)
    z = list(z)

    n = len(tc)
    if any(len(v) != n for v in (vc, m, omega, dipole, kappa, z)):
        raise InvalidInputError(
            "Tc",
            f"the per-component vectors disagree in length: Tc has {n} entries, "
            f"Vc {len(vc)}, M {len(m)}, omega {len(omega)}, dipole {len(dipole)}, "
            f"kappa {len(kappa)}, z {len(z)}",
        )
    if n == 0:
        raise InvalidInputError("Tc", "a mixture of zero components has no viscosity")
    if any(value < 0.0 for value in z):
        bad = next(i for i, value in enumerate(z) if value < 0.0)
        raise InvalidInputError("z", f"z[{bad}] is {z[bad]} but a mole fraction cannot be negative")
    if abs(sum(z) - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "z",
            f"the mole fractions sum to {sum(z)}, not to one. Renormalising them here "
            f"would make a composition error invisible in every number downstream, so "
            f"it is refused instead",
        )

    pure_mu = []
    for i in range(n):
        result = chung_viscosity(omega[i], Tc[i], Vc[i], M[i], dipole[i], kappa[i], T, V)
        warnings.extend(result.warnings)
        pure_mu.append(result.mu.magnitude)

    mu = 0.0
    for i in range(n):
        denominator = 0.0
        for j in range(n):
            ratio = math.sqrt(pure_mu[i] / pure_mu[j]) * math.pow(m[j] / m[i], 0.25)
            phi = (1.0 + ratio) ** 2 / math.sqrt(8.0 * (1.0 + m[i] / m[j]))
            denominator += z[j] * phi
        mu += z[i] * pure_mu[i] / denominator

    return WilkeViscosityResult(mu=from_si(mu, "Pa*s"), warnings=tuple(warnings))
