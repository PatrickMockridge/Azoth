"""``thermal.lmtd_correction_factor`` - the Bowman-Mueller-Nagle correction.

Spec: ``specs/calcs/thermal/lmtd_correction_factor.toml``

The Python twin of ``crates/azoth-thermal/src/lmtd_correction_factor.rs``, and both mirror
``LMTDcorrectionFactor``. The four temperatures enter as differences; ``R``, ``P`` and
``F_t`` are the class's own intermediates.

**Two behaviours not to tidy.** The ``R = 1`` branch is not clamped to one, unlike the
general branch. And an undefined per-shell ``P1`` in the multi-pass conversion answers 1.0
rather than zero, so a temperature cross reads 0.0 at one pass and 1.0 at two.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import LmtdCorrectionFactorResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.lmtd_correction_factor"
MIN_ACCEPTABLE_FT = 0.75


def calc_r(hot_in: float, hot_out: float, cold_in: float, cold_out: float) -> float:
    """`R = (Th_in - Th_out)/(Tc_out - Tc_in)`, a zero cold change answering 1."""
    cold_change = cold_out - cold_in
    if abs(cold_change) < 1e-10:
        return 1.0
    return (hot_in - hot_out) / cold_change


def calc_p(hot_in: float, _hot_out: float, cold_in: float, cold_out: float) -> float:
    """`P = (Tc_out - Tc_in)/(Th_in - Tc_in)`, a zero inlet difference answering 0."""
    max_change = hot_in - cold_in
    if abs(max_change) < 1e-10:
        return 0.0
    return (cold_out - cold_in) / max_change


def calc_ft_single_shell(r: float, p: float) -> float:
    """The single-shell (1-2N) correction factor."""
    if p <= 0.0 or p >= 1.0:
        return 1.0
    w = math.sqrt(r * r + 1.0)
    if abs(r - 1.0) < 1e-6:
        sqrt_two = math.sqrt(2.0)
        numerator = p * sqrt_two
        arg1 = 2.0 - p * (2.0 - sqrt_two)
        arg2 = 2.0 - p * (2.0 + sqrt_two)
        if arg1 <= 0.0 or arg2 <= 0.0 or arg1 / arg2 <= 0.0:
            return 0.0
        denominator = (1.0 - p) * math.log(arg1 / arg2)
        if abs(denominator) < 1e-15:
            return 1.0
        return numerator / denominator
    rp = r * p
    if rp >= 1.0:
        return 0.0
    log_arg1 = (1.0 - p) / (1.0 - rp)
    if log_arg1 <= 0.0:
        return 0.0
    arg2_num = 2.0 - p * (r + 1.0 - w)
    arg2_den = 2.0 - p * (r + 1.0 + w)
    if arg2_num <= 0.0 or arg2_den <= 0.0 or arg2_num / arg2_den <= 0.0:
        return 0.0
    numerator = w * math.log(log_arg1)
    denominator = (r - 1.0) * math.log(arg2_num / arg2_den)
    if abs(denominator) < 1e-15:
        return 1.0
    return min(max(numerator / denominator, 0.0), 1.0)


def convert_p_to_per_shell(p: float, r: float, n: int) -> float:
    """The overall `P` converted to a per-shell `P1` for `n` shells in series."""
    if n <= 1:
        return p
    if abs(r - 1.0) < 1e-6:
        return p / (n - (n - 1.0) * p)
    x = (1.0 - r * p) / (1.0 - p)
    if x <= 0.0:
        return math.nan
    x_to_one_over_n = math.pow(x, 1.0 / n)
    denominator = x_to_one_over_n - r
    if abs(denominator) < 1e-15:
        return math.nan
    return (x_to_one_over_n - 1.0) / denominator


def calc_ft_from_rp(r: float, p: float, shell_passes: int) -> float:
    """The correction factor for `shell_passes` shell passes."""
    if p <= 0.0 or p >= 1.0:
        return 1.0
    if r <= 0.0:
        return 1.0
    passes = max(1, shell_passes)
    per_shell = p
    if passes > 1:
        per_shell = convert_p_to_per_shell(p, r, passes)
        if per_shell <= 0.0 or per_shell >= 1.0 or math.isnan(per_shell):
            return 1.0
    return calc_ft_single_shell(r, per_shell)


def required_shell_passes(hot_in: float, hot_out: float, cold_in: float, cold_out: float) -> float:
    """The first pass count from one to six reaching 0.75, or -1."""
    r = calc_r(hot_in, hot_out, cold_in, cold_out)
    p = calc_p(hot_in, hot_out, cold_in, cold_out)
    for passes in range(1, 7):
        if calc_ft_from_rp(r, p, passes) >= MIN_ACCEPTABLE_FT:
            return float(passes)
    return -1.0


def lmtd_correction_factor(
    hot_inlet_temperature: Q,
    hot_outlet_temperature: Q,
    cold_inlet_temperature: Q,
    cold_outlet_temperature: Q,
    shell_passes: float | None = None,
) -> LmtdCorrectionFactorResult:
    """The LMTD correction factor for a shell-and-tube exchanger.

    Raises:
        InvalidInputError: for a non-integral shell-pass count.
        OutOfRangeError: for a temperature outside the spec's bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    values = {
        "hot_inlet_temperature": input_to_si(spec, "hot_inlet_temperature", hot_inlet_temperature),
        "hot_outlet_temperature": input_to_si(
            spec, "hot_outlet_temperature", hot_outlet_temperature
        ),
        "cold_inlet_temperature": input_to_si(
            spec, "cold_inlet_temperature", cold_inlet_temperature
        ),
        "cold_outlet_temperature": input_to_si(
            spec, "cold_outlet_temperature", cold_outlet_temperature
        ),
    }
    apply_checks(checks.on_input, values.get, warnings)

    if shell_passes is None:
        passes = 1
    else:
        if not math.isfinite(shell_passes) or shell_passes != math.floor(shell_passes):
            raise InvalidInputError(
                "shell_passes", f"{shell_passes} is not a whole shell-pass count"
            )
        passes = max(1, int(shell_passes))

    hot_in = values["hot_inlet_temperature"]
    hot_out = values["hot_outlet_temperature"]
    cold_in = values["cold_inlet_temperature"]
    cold_out = values["cold_outlet_temperature"]
    r = calc_r(hot_in, hot_out, cold_in, cold_out)
    p = calc_p(hot_in, hot_out, cold_in, cold_out)
    return LmtdCorrectionFactorResult(
        r=r,
        p=p,
        ft=calc_ft_from_rp(r, p, passes),
        ft_1_shell_pass=calc_ft_from_rp(r, p, 1),
        ft_2_shell_passes=calc_ft_from_rp(r, p, 2),
        required_shell_passes=required_shell_passes(hot_in, hot_out, cold_in, cold_out),
        warnings=tuple(warnings),
    )
