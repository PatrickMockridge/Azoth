"""``eos.rachford_rice`` - the vapour fraction that solves the Rachford-Rice equation.

```text
g(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))  =  0
```

A *model* rather than a calculation: what the spec pins down is the procedure, and this
module reads the procedure from ``azoth._models_gen`` rather than choosing it.

Spec: ``specs/models/eos/rachford_rice.toml``, which carries the scheme, the stopping
rule, the tolerance and the cap, since two implementations running different procedures
reach different roots.

NeqSim's ``RachfordRice`` carries two apportioning schemes behind a static ``method``
field and this is ``calcBetaNielsen2023``, its default - see the spec's assumptions for
why the alternative is named rather than ported, and for the one difference that matters:
NeqSim answers a root outside ``[0, 1]`` with a clamp, and this returns the root.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import RachfordRiceResult
from azoth.core.warnings import Warning

MODEL_ID = "eos.rachford_rice"

#: A K-value below this is an ion and does not take part in the split.
_ION_THRESHOLD = 1.0e-30

#: NeqSim's ``ThermodynamicModelSettings.phaseFractionMinimumLimit``. The answer for a
#: feed that cannot split, at whichever end the feed lies.
_PHASE_FRACTION_MINIMUM_LIMIT = 1.0e-12


def _sides(K: list[float]) -> tuple[bool, bool]:
    """Whether any K-value lies above one, and whether any lies below it."""
    above = any(value > 1.0 for value in K if value >= _ION_THRESHOLD)
    below = any(value < 1.0 for value in K if value >= _ION_THRESHOLD)
    return above, below


def _residual(z: list[float], K: list[float], beta: float) -> float:
    """``g(beta)``, ions excluded."""
    return sum(
        z[i] * (K[i] - 1.0) / (1.0 + beta * (K[i] - 1.0))
        for i in range(len(K))
        if K[i] >= _ION_THRESHOLD
    )


def _nielsen_2023(z: list[float], K: list[float], tolerance: float, max_iterations: int) -> float:
    """Nielsen & Lia's reformulation, transcribed from NeqSim's ``calcBetaNielsen2023``.

    The final clamp and the two single-phase returns are left out, so this is called only
    where a root exists and always returns it.
    """
    # `h` is `g` at the starting point, and its sign selects which unknown pair is used.
    # Solving in `1 / K` when the root is above one half turns the search back around.
    h = _residual(z, K, 0.5)
    if h > 0.0:
        work_K = [value if value < _ION_THRESHOLD else 1.0 / value for value in K]
    else:
        work_K = list(K)

    Kmax = 0.0
    Kmin = math.inf
    found = False
    for value in work_K:
        if value < _ION_THRESHOLD:
            continue
        if not found:
            Kmax = value
            Kmin = value
            found = True
        elif value < Kmin:
            Kmin = value
        elif value > Kmax:
            Kmax = value
    if not found:
        return _PHASE_FRACTION_MINIMUM_LIMIT

    alpha_min = 1.0 / (1.0 - Kmax)
    alpha_max = 1.0 / (1.0 - Kmin)

    alpha = 0.5
    a = (alpha - alpha_min) / (alpha_max - alpha)
    b = 1.0 / (alpha - alpha_min)

    # The per-component constants of the rescaled form. A component with `K` exactly one
    # has no distance to rescale, and the denominator floor sends its term to zero in
    # `funk` and in `hb` alike - which is what the equation says it contributes.
    c = [0.0] * len(work_K)
    d = [0.0] * len(work_K)
    for i, raw in enumerate(work_K):
        if raw < _ION_THRESHOLD:
            continue
        value = min(max(raw, 1.0e-25), 1.0e25)
        denom = 1.0 - value
        if abs(denom) < 1.0e-25:
            denom = -1.0e-25 if denom < 0.0 else 1.0e-25
        c[i] = 1.0 / denom
        d[i] = (alpha_min - c[i]) / (alpha_max - alpha_min)

    a_max = alpha
    b_max = 1.0e20
    a_min = 0.0
    b_min = 1.0 / (alpha_max - alpha_min)

    for _ in range(max_iterations):
        funk = 0.0
        funk_deriv = 0.0
        hb = 0.0
        hb_deriv = 0.0
        for i, value in enumerate(work_K):
            if value < _ION_THRESHOLD:
                continue
            funk -= z[i] * a * (1.0 + a) / (d[i] + a * (1.0 + d[i]))
            funk_deriv -= (
                z[i] * (a * a + (1.0 + a) * (1.0 + a) * d[i]) / (d[i] + a * (1.0 + d[i])) ** 2
            )
            hb += z[i] * b / (1.0 + b * (alpha_min - c[i]))
            hb_deriv += z[i] / (1.0 + b * (alpha_min - c[i])) ** 2

        if abs(funk) < tolerance and abs(hb) < tolerance:
            break

        if funk > 0.0:
            a_max = a
        else:
            a_min = a
        if hb > 0.0:
            b_max = b
        else:
            b_min = b

        a -= funk / funk_deriv
        if a > a_max or a < a_min:
            a = 0.5 * (a_max + a_min)
        b -= hb / hb_deriv
        if b > b_max or b < b_min:
            b = 0.5 * (b_max + b_min)

    beta = -((1.0 / b) / a - alpha_max)
    if h > 0.0:
        beta = 1.0 - beta
    return beta


def rachford_rice(z: Sequence[float], K: Sequence[float]) -> RachfordRiceResult:
    """The root of the Rachford-Rice equation for a feed and a set of K-values.

    Args:
        z: the overall mole fractions of the feed.
        K: the K-values ``K_i = y_i / x_i`` at the same temperature and pressure. A
            value below ``1e-30`` marks an ion, which is skipped as NeqSim skips it.

    Returns:
        ``beta``, the vapour fraction. Inside ``[0, 1]`` it is the split; outside it is
        the negative flash, returned as the equation gives it. A feed whose K-values do
        not straddle one has no root, and comes back as NeqSim's ``1e-12`` clamp at the
        end the feed lies towards.

    Raises:
        InvalidInputError: if ``z`` and ``K`` differ in length, or ``z`` is empty.
        OutOfRangeError: if any K-value is not positive.
        SolverNotConvergedError: if the iteration hits its cap.

    Example:
        >>> r = rachford_rice([0.6, 0.4], [7.304244305324782, 0.33749596785762953])
        >>> round(r.beta, 12)
        0.84220554758
        >>> r.is_clean
        True
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    feed = list(z)
    ratios = list(K)
    if len(feed) != len(ratios):
        raise InvalidInputError(
            field="K",
            reason=(
                f"the feed has {len(feed)} component(s) and there are {len(ratios)} K-value(s)"
            ),
        )
    if not feed:
        raise InvalidInputError(
            field="z", reason="a feed of zero components has no vapour fraction"
        )

    apply_checks(checks.on_input, {"K": min(ratios)}.get, warnings)

    algorithm = spec["algorithm"]
    above, below = _sides(ratios)
    if above and below:
        beta = _nielsen_2023(feed, ratios, algorithm["tolerance"], algorithm["max_iterations"])
    elif above:
        beta = 1.0 - _PHASE_FRACTION_MINIMUM_LIMIT
    else:
        beta = _PHASE_FRACTION_MINIMUM_LIMIT

    return RachfordRiceResult(beta=beta, warnings=tuple(warnings))
