"""``eos.brine_mixing_scale_evaluator`` - NeqSim's ``BrineMixingScaleEvaluator``.

Two incompatible brines blended linearly at a set of mixing fractions, each blend
evaluated by the activity-corrected electrolyte calculator, and the worst saturation
index reported.

Spec: ``specs/models/eos/brine_mixing_scale_evaluator.toml``.
"""

from __future__ import annotations

import math
from typing import Any, cast

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import BrineMixingPointRow
from azoth.core.result import BrineMixingScaleEvaluatorResult
from azoth.core.units import Q, input_records_to_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.electrolyte_scale_calculator import _state

CALC_ID = "eos.brine_mixing_scale_evaluator"

#: The nine ions the class's `BrineComposition` carries, in its constructor's order.
IONS = (
    "calcium",
    "barium",
    "strontium",
    "magnesium",
    "sodium",
    "potassium",
    "chloride",
    "sulphate",
    "bicarbonate",
)

#: The minerals the worst case is chosen over, in the class's order.
MINERALS = ("CaCO3", "BaSO4", "CaSO4", "SrSO4")


def brine_mixing_scale_evaluator(
    brine_a: list[dict[str, Q]],
    brine_b: list[dict[str, Q]],
    temperature: Q,
    pressure: Q,
    ph: float,
    steps: int,
) -> BrineMixingScaleEvaluatorResult:
    """The mineral saturation of two brines across their mixing range.

    Args:
        brine_a: the first brine end member, exactly one row.
        brine_b: the second brine end member, exactly one row.
        temperature: the temperature every blend is evaluated at.
        pressure: the pressure every blend is evaluated at.
        ph: the pH every blend is evaluated at.
        steps: the number of mixing points, floored at two as the class floors it.

    Returns:
        One row per mixing fraction and the worst case over them.

    Raises:
        InvalidInputError: for a brine list that is not exactly one row.
        OutOfRangeError: for a negative ion, or a non-positive temperature or pressure.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> formation = [{
        ...     "calcium": q(5.0, "kg/m**3"), "barium": q(0.3, "kg/m**3"),
        ...     "strontium": q(0.8, "kg/m**3"), "magnesium": q(0.5, "kg/m**3"),
        ...     "sodium": q(30.0, "kg/m**3"), "potassium": q(0.4, "kg/m**3"),
        ...     "chloride": q(55.0, "kg/m**3"), "sulphate": q(0.0, "kg/m**3"),
        ...     "bicarbonate": q(0.4, "kg/m**3")}]
        >>> injection = [{
        ...     "calcium": q(0.4, "kg/m**3"), "barium": q(0.0, "kg/m**3"),
        ...     "strontium": q(0.008, "kg/m**3"), "magnesium": q(1.3, "kg/m**3"),
        ...     "sodium": q(11.0, "kg/m**3"), "potassium": q(0.4, "kg/m**3"),
        ...     "chloride": q(19.0, "kg/m**3"), "sulphate": q(2.7, "kg/m**3"),
        ...     "bicarbonate": q(0.14, "kg/m**3")}]
        >>> r = brine_mixing_scale_evaluator(
        ...     formation, injection, q(333.15, "K"), q(1.0e7, "Pa"), 6.0, 11)
        >>> r.worst_mineral
        'BaSO4'
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    def single(rows: list[dict[str, Q]], role: str) -> dict[str, float]:
        if len(rows) != 1:
            raise InvalidInputError(role, f"`{role}` is one brine and this is {len(rows)} row(s)")
        return cast("dict[str, float]", input_records_to_si(spec, role, rows)[0])

    first = single(brine_a, "brine_a")
    second = single(brine_b, "brine_b")
    temperature_si = input_to_si(spec, "temperature", temperature)
    pressure_si = input_to_si(spec, "pressure", pressure)
    apply_checks(
        checks.on_input,
        {"temperature": temperature_si, "pressure": pressure_si, "ph": ph}.get,
        warnings,
    )

    for label, brine in (("brine_a", first), ("brine_b", second)):
        for name in IONS:
            if brine[name] < 0.0:
                raise OutOfRangeError(f"{label}.{name}", brine[name], "must be non-negative")

    steps = max(2, int(steps))
    rows: list[BrineMixingPointRow] = []
    worst_si = -math.inf
    worst_mineral = "NONE"
    worst_fraction_a = math.nan

    def blend(name: str, fraction_a: float, fraction_b: float) -> float:
        return first[name] * fraction_a + second[name] * fraction_b

    for index in range(steps):
        fraction_a = index / (steps - 1)
        fraction_b = 1.0 - fraction_a

        point = _state(
            {
                **{name: blend(name, fraction_a, fraction_b) for name in IONS},
                "iron": 0.0,
                "carbonate": 0.0,
                "temperature": temperature_si,
                "pressure": pressure_si,
                "ph": ph,
            }
        )
        row = BrineMixingPointRow(
            fraction_a=fraction_a,
            si_caco3=point["si_caco3"],
            si_baso4=point["si_baso4"],
            si_caso4=point["si_caso4"],
            si_srso4=point["si_srso4"],
        )
        for mineral, index_value in zip(
            MINERALS,
            (row.si_caco3, row.si_baso4, row.si_caso4, row.si_srso4),
            strict=True,
        ):
            if index_value > worst_si:
                worst_si = index_value
                worst_mineral = mineral
                worst_fraction_a = fraction_a
        rows.append(row)

    return BrineMixingScaleEvaluatorResult(
        points=tuple(rows),
        worst_fraction_a=worst_fraction_a,
        worst_mineral=worst_mineral,
        worst_saturation_index=worst_si,
        warnings=tuple(warnings),
    )
