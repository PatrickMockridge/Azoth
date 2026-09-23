"""``eos.iapws_henry_law`` - the Henry constant of a gas in water.

Spec: ``specs/calcs/eos/iapws_henry_law.toml``. A *direct* calculation: no iteration,
so no algorithm block.

This is the pure-Python reference: a second, independent expression of the same
physics as the Rust kernel, held to it by ``test_cross_impl.py``. The table is written
here as data and the equation once; the Rust kernel spells the rows out as a match.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HenryStatus, IapwsHenryLawResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.iapws_henry_law"

#: Water's critical temperature, the reduced temperature's denominator.
WATER_CRITICAL_TEMPERATURE = 647.096

#: Water's critical pressure in MPa, the guideline's reference pressure.
WATER_CRITICAL_PRESSURE_MPA = 22.064

#: The guideline is stated for liquid water, and this is where the liquid ends.
CORRELATION_MINIMUM_TEMPERATURE = 273.15

#: The six coefficients of the Wagner saturation-pressure series and their exponents.
VAPOR_PRESSURE_SERIES = (
    (-7.85951783, 1.0),
    (1.84408259, 1.5),
    (-11.7866497, 3.0),
    (22.6807411, 3.5),
    (-15.9618719, 4.0),
    (1.80122502, 7.5),
)

#: The guideline's gas table: ``(a, b, c, fitted_min, fitted_max, rms log residual)``.
#:
#: The fit windows differ per gas and are **not** the domain: a row evaluated at a
#: temperature outside its window is returned with ``GUIDELINE_EXTRAPOLATION``.
TABLE: dict[str, tuple[float, float, float, float, float, float]] = {
    "he": (-3.52839, 7.12983, 4.47770, 273.21, 553.18, 0.0341),
    "ne": (-3.18301, 5.31448, 5.43774, 273.20, 543.36, 0.0577),
    "ar": (-8.40954, 4.29587, 10.52779, 273.19, 568.36, 0.0443),
    "kr": (-8.97358, 3.61508, 11.29963, 273.19, 525.56, 0.0434),
    "xe": (-14.21635, 4.00041, 15.60999, 273.22, 574.85, 0.0363),
    "h2": (-4.73284, 6.08954, 6.06066, 273.15, 636.09, 0.0517),
    "n2": (-9.67578, 4.72162, 11.70585, 278.12, 636.46, 0.0372),
    "o2": (-9.44833, 4.43822, 11.42005, 274.15, 616.52, 0.0377),
    "co": (-10.52862, 5.13259, 12.01421, 278.15, 588.67, 0.0039),
    "co2": (-8.55445, 4.01195, 9.52345, 274.19, 642.66, 0.0528),
    "h2s": (-4.51499, 5.23538, 4.42126, 273.15, 533.09, 0.0408),
    "ch4": (-10.44708, 4.66491, 12.12986, 275.46, 633.11, 0.0386),
    "c2h6": (-19.67563, 4.51222, 20.62567, 275.44, 473.46, 0.0259),
    "sf6": (-16.56118, 2.15289, 20.35440, 283.14, 505.55, 0.0505),
}

#: The name a row answers to, which is NeqSim's ``findGas``.
ALIASES: dict[str, str] = {
    "he": "he",
    "helium": "he",
    "ne": "ne",
    "neon": "ne",
    "ar": "ar",
    "argon": "ar",
    "kr": "kr",
    "krypton": "kr",
    "xe": "xe",
    "xenon": "xe",
    "h2": "h2",
    "hydrogen": "h2",
    "n2": "n2",
    "nitrogen": "n2",
    "o2": "o2",
    "oxygen": "o2",
    "co": "co",
    "carbon monoxide": "co",
    "co2": "co2",
    "carbon dioxide": "co2",
    "h2s": "h2s",
    "hydrogen sulfide": "h2s",
    "hydrogen sulphide": "h2s",
    "ch4": "ch4",
    "methane": "ch4",
    "c2h6": "c2h6",
    "ethane": "c2h6",
    "sf6": "sf6",
    "sulfur hexafluoride": "sf6",
    "sulphur hexafluoride": "sf6",
}

#: Megapascals per bar, and pascals per bar: the guideline evaluates in MPa, NeqSim
#: returns bar, and this model returns pascals.
MPA_TO_BAR = 10.0
BAR_TO_PA = 1.0e5


def gas_from_name(name: str) -> str | None:
    """The table row a component name means, or ``None``.

    The table is its own vocabulary rather than this library's component list:
    ``krypton``, ``xenon`` and ``carbon monoxide`` are rows here and are not components
    ``azoth.eos.components`` carries data for, so the resolution is by string.
    """
    return ALIASES.get(name.strip().lower())


def iapws_henry_law(gas: str, T: Q) -> IapwsHenryLawResult:
    """The Henry constant of a gas in water, and the range it was fitted over.

    The standard state is limiting ``f/x`` at pure-water saturation, so the constant is
    a pressure per mole fraction and is large for a sparingly soluble gas. This is the
    **second** Henry arm this library carries: ``azoth.eos.reference._henry`` is the
    database correlation out of the compiled component table, and which one a phase
    takes is the phase's decision.

    Args:
        gas: one of the guideline's 14 rows, by formula or by name.
        T: absolute temperature, inside liquid water.

    Returns:
        ``kH`` in pascals, its logarithm, ``d(ln kH)/dT``, whether ``T`` is inside the
        row's fitted window, and the row's reported fit residual.

    Raises:
        OutOfRangeError: if ``T`` is outside ``[273.15, 647.096)`` K, where the
            guideline is not liquid water and has no reference state.
        InvalidInputError: if ``gas`` names no row.

    Example:
        >>> import azoth
        >>> r = iapws_henry_law("methane", azoth.ureg.Quantity(298.15, "K"))
        >>> round(r.henry.magnitude, 3)
        3947965646.001
        >>> r.status
        <HenryStatus.WITHIN_FITTED_RANGE: 'within_fitted_range'>
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    apply_checks(checks.on_input, {"T": t_si}.get, warnings)

    row = gas_from_name(gas)
    if row is None:
        raise InvalidInputError(
            "gas",
            f"unknown gas `{gas}`; the guideline's table carries `he`, `ne`, `ar`, `kr`, "
            f"`xe`, `h2`, `n2`, `o2`, `co`, `co2`, `h2s`, `ch4`, `c2h6` and `sf6`",
        )
    a, b, c, fitted_min, fitted_max, rms_log_residual = TABLE[row]

    # The domain is liquid water. It is refused here rather than left to the range
    # table, so that the two bounds are stated once, next to the equation that needs
    # them.
    t = t_si
    if (
        not math.isfinite(t)
        or t < CORRELATION_MINIMUM_TEMPERATURE
        or t >= WATER_CRITICAL_TEMPERATURE
    ):
        raise OutOfRangeError(
            "T",
            t,
            f"the guideline is stated for liquid water, "
            f"[{CORRELATION_MINIMUM_TEMPERATURE}, {WATER_CRITICAL_TEMPERATURE}) K, and the "
            f"reference state ceases to exist above it",
        )

    tr = t / WATER_CRITICAL_TEMPERATURE
    tau = 1.0 - tr
    series = sum(coefficient * tau**exponent for coefficient, exponent in VAPOR_PRESSURE_SERIES)
    log_pressure_mpa = math.log(WATER_CRITICAL_PRESSURE_MPA) + series / tr
    log_henry_mpa = log_pressure_mpa + a / tr + b * tau**0.355 / tr + c * tr**-0.41 * math.exp(tau)
    henry_bar = math.exp(log_henry_mpa) * MPA_TO_BAR
    henry_pa = henry_bar * BAR_TO_PA

    # The logarithmic derivative, in 1/K. Written against `Tr` and divided once, which
    # is how the guideline states it and how NeqSim evaluates it.
    series_derivative = sum(
        -coefficient * exponent * tau ** (exponent - 1.0)
        for coefficient, exponent in VAPOR_PRESSURE_SERIES
    )
    derivative_by_tr = (
        series_derivative / tr
        - series / tr**2
        - a / tr**2
        + b * (-0.355 * tau**-0.645 / tr - tau**0.355 / tr**2)
    )
    derivative_by_tr += c * tr**-0.41 * math.exp(tau) * (-0.41 / tr - 1.0)

    status = (
        HenryStatus.WITHIN_FITTED_RANGE
        if fitted_min <= t <= fitted_max
        else HenryStatus.GUIDELINE_EXTRAPOLATION
    )

    return IapwsHenryLawResult(
        # The logarithm is taken of the pascals this model returns rather than of the
        # bar figure, so that no two parts of one result are in different scales.
        henry=from_si(henry_pa, "Pa"),
        ln_henry=math.log(henry_pa),
        d_ln_henry_d_t=derivative_by_tr / WATER_CRITICAL_TEMPERATURE,
        status=status,
        rms_log_residual=rms_log_residual,
        warnings=tuple(warnings),
    )
