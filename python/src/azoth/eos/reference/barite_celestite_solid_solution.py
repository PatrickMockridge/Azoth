"""``eos.barite_celestite_solid_solution`` - NeqSim's ``BariteCelestiteSolidSolution``.

A Margules one-parameter solid between BaSO4 and SrSO4, with the composition found by the
class's bisection and the reported index taken from the Ba arm alone.

Spec: ``specs/models/eos/barite_celestite_solid_solution.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import BariteCelestiteSolidSolutionResult
from azoth.core.warnings import Warning

CALC_ID = "eos.barite_celestite_solid_solution"


def _gamma(margules: float, x: float) -> float:
    other = 1.0 - x
    return math.exp(margules * other * other)


def barite_celestite_solid_solution(
    margules: float,
    a_ba: float,
    a_sr: float,
    a_so4: float,
    ksp_baso4: float,
    ksp_srso4: float,
) -> BariteCelestiteSolidSolutionResult:
    """The solid's composition and the class's saturation index.

    Args:
        margules: the Margules interaction parameter ``W/(RT)``.
        a_ba: the aqueous Ba2+ activity.
        a_sr: the aqueous Sr2+ activity.
        a_so4: the aqueous SO4 2- activity.
        ksp_baso4: the pure BaSO4 solubility product.
        ksp_srso4: the pure SrSO4 solubility product.

    Returns:
        The mole fractions, the Ba arm's index and whether it is positive.

    Raises:
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> r = barite_celestite_solid_solution(2.3, 0.01, 0.0, 0.05, 1.08e-10, 3.44e-7)
        >>> r.x_baso4
        1.0
        >>> abs(r.total_saturation_index - 6.665546248849069) < 1e-12
        True
    """
    from azoth._models_gen import model

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    apply_checks(
        checks.on_input,
        {
            "margules": margules,
            "a_ba": a_ba,
            "a_sr": a_sr,
            "a_so4": a_so4,
            "ksp_baso4": ksp_baso4,
            "ksp_srso4": ksp_srso4,
        }.get,
        warnings,
    )

    def pure_index(iap: float, ksp: float) -> float | None:
        return math.log10(iap / ksp) if iap > 0.0 and ksp > 0.0 else None

    if a_ba <= 0.0 and a_sr <= 0.0:
        x_ba: float = 0.5
        total: float | None = None
    elif a_sr <= 1.0e-30:
        x_ba = 1.0
        total = pure_index(a_ba * a_so4, ksp_baso4)
    elif a_ba <= 1.0e-30:
        x_ba = 0.0
        total = pure_index(a_sr * a_so4, ksp_srso4)
    else:

        def objective(x: float) -> float:
            x_sr = 1.0 - x
            r_ba = a_ba * a_so4 / (ksp_baso4 * x * _gamma(margules, x))
            r_sr = a_sr * a_so4 / (ksp_srso4 * x_sr * _gamma(margules, x_sr))
            return r_ba - r_sr

        x_lo, x_hi = 0.001, 0.999
        x_ba = 0.5
        for _ in range(100):
            x_mid = 0.5 * (x_lo + x_hi)
            f_mid = objective(x_mid)
            if abs(f_mid) < 1.0e-12 or (x_hi - x_lo) < 1.0e-10:
                x_ba = x_mid
                break
            if f_mid * objective(x_lo) < 0.0:
                x_hi = x_mid
            else:
                x_lo = x_mid
            x_ba = x_mid
        si_ba = a_ba * a_so4 / (ksp_baso4 * x_ba * _gamma(margules, x_ba))
        total = math.log10(max(si_ba, 1.0e-30))

    return BariteCelestiteSolidSolutionResult(
        x_baso4=x_ba,
        x_srso4=1.0 - x_ba,
        total_saturation_index=total,
        supersaturated=total is not None and total > 0.0,
        warnings=tuple(warnings),
    )
