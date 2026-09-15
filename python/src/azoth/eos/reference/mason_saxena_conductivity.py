"""``eos.mason_saxena_conductivity`` - the gas mixture conductivity from
Mason-Saxena mixing over Chung's pure-component conductivities.

Spec: ``specs/models/eos/mason_saxena_conductivity.toml``, which records the mixing
rule and points the pure-component half at :func:`azoth.eos.chung_conductivity`. It
is a *direct* model: no iteration, so no algorithm block.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MasonSaxenaConductivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.chung_conductivity import chung_conductivity

MODEL_ID = "eos.mason_saxena_conductivity"


def mason_saxena_conductivity(
    Cv0: Sequence[Q],
    M: Sequence[Q],
    omega: Sequence[float],
    Tc: Sequence[Q],
    Vc: Sequence[Q],
    dipole: Sequence[float],
    kappa: Sequence[float],
    T: Q,
    z: Sequence[float],
) -> MasonSaxenaConductivityResult:
    """The gas thermal conductivity of a mixture, by Mason-Saxena mixing over Chung.

    Each pure-component conductivity is :func:`azoth.eos.chung_conductivity`'s at
    ``T``, then mixed: ``k = sum_i z_i k_i / sum_j z_j A_ij`` with ``A_ij = (1 +
    sqrt(k_i/k_j) (M_i/M_j)**0.25)**2 / sqrt(8 (1 + M_i/M_j))``.

    Args:
        Cv0: ideal-gas heat capacity at constant volume of each component at ``T``.
        M: molar mass of each component.
        omega: acentric factor of each component.
        Tc: critical temperature of each component.
        Vc: critical volume of each component.
        dipole: dipole moment of each component, in debye.
        kappa: viscosity correction factor of each component.
        T: absolute temperature.
        z: mole fractions; non-negative and summing to one.

    Returns:
        The gas mixture thermal conductivity, in ``W/(m*K)``.

    Raises:
        InvalidInputError: if the vectors disagree in length, or ``z`` is not a
            composition.
        OutOfRangeError: if ``T`` is not positive, or a component's ``Tc``/``Vc`` is
            not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> cv = [q(27.544151394, "J/(mol*K)"), q(65.809299386, "J/(mol*K)")]
        >>> m = [q(0.016043, "kg/mol"), q(0.044097, "kg/mol")]
        >>> tc = [q(190.56, "K"), q(369.83, "K")]
        >>> vc = [q(9.9e-5, "m**3/mol"), q(0.000203, "m**3/mol")]
        >>> omega = [0.0115, 0.1523]
        >>> t = q(300.0, "K")
        >>> z = [0.5, 0.5]
        >>> r = mason_saxena_conductivity(cv, m, omega, tc, vc, [0.0, 0.0], [0.0, 0.0], t, z)
        >>> round(r.k.magnitude, 16)
        0.02506205948004617
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "T": input_to_si(spec, "T", T),
    }
    apply_checks(checks.on_input, values.get, warnings)

    cv0 = [input_to_si(spec, "Cv0", value) for value in Cv0]
    m = [input_to_si(spec, "M", value) for value in M]
    omega = list(omega)
    tc = [input_to_si(spec, "Tc", value) for value in Tc]
    vc = [input_to_si(spec, "Vc", value) for value in Vc]
    dipole = list(dipole)
    kappa = list(kappa)
    z = list(z)

    n = len(cv0)
    if any(len(v) != n for v in (m, omega, tc, vc, dipole, kappa, z)):
        raise InvalidInputError(
            "Cv0",
            f"the per-component vectors disagree in length: Cv0 has {n} entries, "
            f"M {len(m)}, omega {len(omega)}, Tc {len(tc)}, Vc {len(vc)}, "
            f"dipole {len(dipole)}, kappa {len(kappa)}, z {len(z)}",
        )
    if n == 0:
        raise InvalidInputError("Cv0", "a mixture of zero components has no conductivity")
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

    pure_k = []
    for i in range(n):
        result = chung_conductivity(Cv0[i], M[i], omega[i], Tc[i], Vc[i], dipole[i], kappa[i], T)
        warnings.extend(result.warnings)
        pure_k.append(result.k.magnitude)

    k = 0.0
    for i in range(n):
        denominator = 0.0
        for j in range(n):
            ratio = math.sqrt(pure_k[i] / pure_k[j]) * math.pow(m[i] / m[j], 0.25)
            aij = (1.0 + ratio) ** 2 / math.sqrt(8.0 * (1.0 + m[i] / m[j]))
            denominator += z[j] * aij
        k += z[i] * pure_k[i] / denominator

    return MasonSaxenaConductivityResult(k=from_si(k, "W/(m*K)"), warnings=tuple(warnings))
