"""``characterization.tbp_cut_boundaries`` - a cumulative TBP curve cut into an assay.

Spec: ``specs/models/characterization/tbp_cut_boundaries.toml``. Oracle:
``validation/neqsim/captures/tbp_cut_boundaries_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/tbp_cut_boundaries.rs`` line for line.

# The curve is stated as boundaries

``addTBPCutBoundariesKelvin`` takes N cumulative yields, N boiling points and N-1 gravities, and
makes N-1 cuts of them. Each interval's share is ``(cumulative[i + 1] - cumulative[i]) / 100``, its
gravity is the one stated beside it, and ``withBoilingRangeKelvin(lower, upper)`` keeps both
boundaries and the midpoint between them - so a cut carries a range and an average rather than only
an average.

# The Celsius entry point is the same arithmetic

``addTBPCutBoundariesCelsius`` adds 273.15 to its temperatures and calls this one, which is the
crossing ``to_si`` makes at the boundary. So there is one id here and a caller states its
temperatures in whatever unit it holds them in.
"""

from __future__ import annotations

from collections.abc import Sequence
from itertools import pairwise

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TbpCutBoundariesResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.tbp_cut_boundaries"

#: ``OilAssayCharacterisation.PERCENT_TOLERANCE``, on the curve's two ends and on each yield.
PERCENT_TOLERANCE = 1.0e-8


def _dimensionless(value: float | Q) -> float:
    """A declared dimensionless input as its magnitude, quantity or not.

    A cumulative yield and a specific gravity carry no unit, so a bare number is one - there is no
    dimension to be wrong about. The temperatures stay strict: a bare number where a boiling point
    is expected is the mistake this library exists to make impossible.
    """
    if isinstance(value, int | float):
        return float(value)
    return float(value.to("dimensionless").magnitude)


def tbp_cut_boundaries(
    name_prefix: str,
    cumulative_volume_percent: Sequence[float],
    boiling_point: Sequence[Q],
    specific_gravity: Sequence[float],
) -> TbpCutBoundariesResult:
    """A cumulative TBP curve cut into an assay, one row per interval.

    ``cumulative_volume_percent`` and ``boiling_point`` are the curve's boundaries - N of each -
    and ``specific_gravity`` is one per interval, so one shorter than either.

    Raises:
        InvalidInputError: if the prefix is empty, if the three vectors are not N, N and N-1 long,
            if the yields do not span zero to a hundred, if either the yields or the temperatures
            fail to increase strictly, or if a gravity is not finite and positive.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    apply_checks(checks.on_input, lambda _name: None, warnings)

    yields = [_dimensionless(value) for value in cumulative_volume_percent]
    temperatures = [input_to_si(spec, "boiling_point", value) for value in boiling_point]
    gravities = [_dimensionless(value) for value in specific_gravity]

    if not name_prefix.strip():
        raise InvalidInputError(
            "name_prefix",
            "a cut's name is its prefix and its position, so an empty prefix names nothing",
        )
    count = len(yields)
    if count < 2 or len(temperatures) != count or len(gravities) != count - 1:
        raise InvalidInputError(
            "cumulative_volume_percent",
            f"a TBP curve is N cumulative yields, N temperatures and N-1 gravities: {count} "
            f"yields, {len(temperatures)} temperatures and {len(gravities)} gravities",
        )
    if abs(yields[0]) > PERCENT_TOLERANCE or abs(yields[count - 1] - 100.0) > PERCENT_TOLERANCE:
        raise InvalidInputError(
            "cumulative_volume_percent",
            f"a cumulative volume yield spans zero to a hundred percent, and this curve runs from "
            f"{yields[0]} to {yields[count - 1]}",
        )
    for index, yield_percent in enumerate(yields):
        if (
            yield_percent != yield_percent
            or yield_percent < -PERCENT_TOLERANCE
            or yield_percent > 100.0 + PERCENT_TOLERANCE
        ):
            raise InvalidInputError(
                f"cumulative_volume_percent[{index}]",
                f"a cumulative yield is finite and between 0 and 100, and this is {yield_percent}",
            )
        if index > 0 and not yield_percent > yields[index - 1]:
            raise InvalidInputError(
                f"cumulative_volume_percent[{index}]",
                f"the cumulative yields must strictly increase, and {yield_percent} follows "
                f"{yields[index - 1]}",
            )
    for index, temperature in enumerate(temperatures):
        if temperature != temperature or not temperature > 0.0:
            raise InvalidInputError(
                f"boiling_point[{index}]",
                f"a boiling-point boundary is finite and positive, and this is {temperature}",
            )
        if index > 0 and not temperature > temperatures[index - 1]:
            raise InvalidInputError(
                f"boiling_point[{index}]",
                f"the boiling-point boundaries must strictly increase, and {temperature} follows "
                f"{temperatures[index - 1]}",
            )
    for index, gravity in enumerate(gravities):
        if gravity != gravity or not gravity > 0.0:
            raise InvalidInputError(
                f"specific_gravity[{index}]",
                f"a specific gravity is finite and positive, and this is {gravity}",
            )

    apply_checks(checks.derived, lambda _name: None, warnings)

    return TbpCutBoundariesResult(
        cut_name=tuple(f"{name_prefix}{index + 1}" for index in range(len(gravities))),
        volume_fraction=tuple(
            (yields[index + 1] - yields[index]) / 100.0 for index in range(count - 1)
        ),
        specific_gravity=tuple(gravities),
        lower_boiling_point=tuple(from_si(value, "K") for value in temperatures[:-1]),
        upper_boiling_point=tuple(from_si(value, "K") for value in temperatures[1:]),
        average_boiling_point=tuple(
            from_si(0.5 * (lower + upper), "K") for lower, upper in pairwise(temperatures)
        ),
        warnings=tuple(warnings),
    )
