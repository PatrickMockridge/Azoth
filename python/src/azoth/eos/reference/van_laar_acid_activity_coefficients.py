"""``eos.van_laar_acid_activity_coefficients`` - the activity coefficients from the
Taleb-Ponche-Mirabel Van Laar model for the water/nitric/sulfuric acid system.

Spec: ``specs/models/eos/van_laar_acid_activity_coefficients.toml``. A *direct* model:
no iteration, so no algorithm block.

This is the pure-Python reference: a second, independent expression of the same
physics as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import VanLaarAcidActivityCoefficientsResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import VanLaarAcidParameters

MODEL_ID = "eos.van_laar_acid_activity_coefficients"

#: The activity a species the model does not cover is given, so a liquid phase rejects
#: it. NeqSim's ``ComponentGEVanLaarAcid.NON_MODELED_COMPONENT_PENALTY``.
NON_MODELED_COMPONENT_PENALTY = 1.0e12

#: System I (H2O-HNO3) Van Laar exponent ``B`` for water. ``B_I_2 = 1/B_I_1``.
B_I_1 = 0.5695
#: System II (H2O-H2SO4) Van Laar exponent ``B`` for water. ``B_II_3 = 1/B_II_1``.
B_II_1 = 0.527
#: System III (HNO3-H2SO4) Van Laar exponent ``B`` for nitric acid. ``B_III_3 = 1/B_III_2``.
B_III_2 = 0.4
#: System III (HNO3-H2SO4) Van Laar parameter ``A`` for nitric acid, fitted at 273 K.
A_III_2 = -250.52
#: System III (HNO3-H2SO4) Van Laar parameter ``A`` for sulfuric acid, fitted at 273 K.
A_III_3 = -100.21


def _a_i1(t: float) -> float:
    return -391.43 - 7.44e4 / t


def _a_i2(t: float) -> float:
    return -627.739 - 1.406e5 / t


def _a_ii1(t: float) -> float:
    return 2.989e3 - 2.147e6 / t + 2.33e8 / (t * t)


def _a_ii3(t: float) -> float:
    return 5.672e3 - 4.074e6 / t + 4.421e8 / (t * t)


def _t_log10_gamma_water(x1: float, x2: float, x3: float, t: float) -> float:
    """``T log10(gamma_1)``, equation (10a) of Taleb et al. (1996)."""
    num1 = _a_i1(t) * (x2 * x2 + B_III_2 * x2 * x3) - 0.5 * A_III_2 * B_I_1 * x2 * x3
    den1 = B_I_1 * x1 + x2 + B_III_2 * x3
    num2 = _a_ii1(t) * (x3 * x3 + B_III_2 * x2 * x3) - 0.5 * A_III_3 * B_II_1 * x2 * x3
    den2 = B_II_1 * x1 + B_III_2 * x2 + x3
    return num1 / (den1 * den1) + num2 / (den2 * den2)


def _t_log10_gamma_nitric(x1: float, x2: float, x3: float, t: float) -> float:
    """``T log10(gamma_2)``, equation (10b)."""
    b_i_2 = 1.0 / B_I_1
    b_ii_3 = 1.0 / B_II_1
    num1 = _a_i2(t) * (x1 * x1 + b_ii_3 * x1 * x3) - 0.5 * _a_ii3(t) * b_i_2 * x1 * x3
    den1 = x1 + x2 * b_i_2 + b_ii_3 * x3
    num2 = A_III_3 * (x3 * x3 + B_II_1 * x1 * x3) - 0.5 * _a_ii1(t) * B_III_2 * x1 * x3
    den2 = B_II_1 * x1 + B_III_2 * x2 + x3
    return num1 / (den1 * den1) + num2 / (den2 * den2)


def _t_log10_gamma_sulfuric(x1: float, x2: float, x3: float, t: float) -> float:
    """``T log10(gamma_3)``, equation (10c)."""
    b_i_2 = 1.0 / B_I_1
    b_ii_3 = 1.0 / B_II_1
    b_iii_3 = 1.0 / B_III_2
    num1 = _a_ii3(t) * (x1 * x1 + b_i_2 * x1 * x2) - 0.5 * _a_i2(t) * b_ii_3 * x1 * x2
    den1 = x1 + x2 * b_i_2 + b_ii_3 * x3
    num2 = A_III_2 * (x2 * x2 + B_I_1 * x1 * x2) - 0.5 * _a_i1(t) * b_iii_3 * x1 * x2
    den2 = B_I_1 * x1 + x2 + b_iii_3 * x3
    return num1 / (den1 * den1) + num2 / (den2 * den2)


def van_laar_acid_activity_coefficients(
    params: VanLaarAcidParameters,
    T: Q,
    x: Sequence[float],
) -> VanLaarAcidActivityCoefficientsResult:
    """The activity coefficients of a mixture containing water, nitric acid and
    sulfuric acid, from the Van Laar model of Taleb, Ponche and Mirabel (1996).

    The three modelled species are evaluated on their own mole-fraction basis: the
    fractions of the acids in ``x`` are renormalised to sum to one, so a dissolved
    carrier gas does not enter the ternary expression. A component the model does not
    cover is given :data:`NON_MODELED_COMPONENT_PENALTY` rather than being refused,
    which is NeqSim's behaviour and pushes it out of the liquid phase.

    All the model's logarithms are base-10, and ``gamma = 10**(T log10 gamma / T)``.
    ``params`` is resolved by name through
    :func:`azoth.eos.components.van_laar_acid_parameters`; ``x`` is checked rather than
    renormalised.

    Args:
        params: the resolved acid identities, by component.
        T: absolute temperature.
        x: mole fractions; non-negative and summing to one.

    Returns:
        The natural logarithm and the value of each activity coefficient.

    Raises:
        InvalidInputError: if ``params`` and ``x`` disagree in length, or ``x`` is not a
            composition.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import van_laar_acid_parameters
        >>> q = azoth.ureg.Quantity
        >>> r = van_laar_acid_activity_coefficients(
        ...     van_laar_acid_parameters(["water", "hno3", "h2so4"]),
        ...     q(250.0, "K"), [0.5, 0.3, 0.2],
        ... )
        >>> f"{r.gamma[0]:.18f}"
        '0.008700227753081819'
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"T": input_to_si(spec, "T", T)}
    apply_checks(checks.on_input, values.get, warnings)

    t = input_to_si(spec, "T", T)
    x = list(x)
    index = list(params.acid_index)

    n = len(x)
    if n == 0:
        raise InvalidInputError("x", "a mixture of zero components has no activity coefficient")
    if len(index) != n:
        raise InvalidInputError(
            "components",
            f"the mixture has {len(index)} components but `x` has {n} entries",
        )
    if any(value < 0.0 for value in x):
        bad = next(i for i, value in enumerate(x) if value < 0.0)
        raise InvalidInputError("x", f"x[{bad}] is {x[bad]} but a mole fraction cannot be negative")
    if abs(sum(x) - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "x",
            f"the mole fractions sum to {sum(x)}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so "
            "it is refused instead",
        )

    # The ternary expression is defined on the three acids' own basis. Dividing through
    # by the total is NeqSim's form rather than a step that moves the answer: each of the
    # three expressions is homogeneous of degree zero, so it takes the same value on any
    # scaling of the basis. Kept because it is what the port source does, and because it
    # makes the empty basis below a case rather than a division by zero.
    acids = [0.0, 0.0, 0.0]
    for kind, fraction in zip(index, x, strict=True):
        if 1 <= kind <= 3:
            acids[kind - 1] += fraction
    total = sum(acids)
    # No acid at all leaves the basis undefined; NeqSim reads that as pure water, which
    # makes every acid's gamma its own infinite-dilution value in water.
    if total <= 0.0:
        x1, x2, x3 = 1.0, 0.0, 0.0
    else:
        x1, x2, x3 = (acids[0] / total, acids[1] / total, acids[2] / total)

    ternary = [
        10.0 ** (_t_log10_gamma_water(x1, x2, x3, t) / t),
        10.0 ** (_t_log10_gamma_nitric(x1, x2, x3, t) / t),
        10.0 ** (_t_log10_gamma_sulfuric(x1, x2, x3, t) / t),
    ]

    ln_gamma: list[float] = []
    gamma: list[float] = []
    for kind in index:
        value = ternary[kind - 1] if 1 <= kind <= 3 else NON_MODELED_COMPONENT_PENALTY
        gamma.append(value)
        ln_gamma.append(math.log(value))

    return VanLaarAcidActivityCoefficientsResult(
        ln_gamma=tuple(ln_gamma),
        gamma=tuple(gamma),
        warnings=tuple(warnings),
    )
