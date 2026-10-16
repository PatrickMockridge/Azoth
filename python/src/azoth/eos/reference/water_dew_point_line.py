"""``eos.water_dew_point_line``: the Python reference.

The mirror of ``crates/azoth-eos/src/water_dew_point_line.rs``. NeqSim's
``WaterDewPointEquilibriumLine``: step the pressure in ``10`` equal increments from
``minimum_pressure`` to ``maximum_pressure``, call the multiphase water dew point at each, and
record the temperature the fluid came back at.

**It holds no arithmetic of its own** - it is ``eos.water_dew_point_multiphase`` ten times - and
**a point that does not solve keeps the previous point's temperature**, which is upstream's own
``try``-and-log behaviour. ``failed_points`` is the marker upstream does not give.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import AzothError, InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WaterDewPointLineResult
from azoth.core.units import Q, from_si, input_to_si, ureg
from azoth.core.warnings import Warning
from azoth.eos.reference.water_dew_point_multiphase import water_dew_point_multiphase

MODEL_ID = "eos.water_dew_point_line"

#: How many points the class sweeps. Its own field, with no setter.
NUMBER_OF_POINTS = 10


def water_dew_point_line(
    components: list[str],
    z: list[float],
    minimum_pressure: Q,
    maximum_pressure: Q,
    eos: str = "srk",
) -> WaterDewPointLineResult:
    """The water dew point at each of ten pressures from ``minimum_pressure`` to
    ``maximum_pressure``.

    Raises:
        InvalidInputError: if ``minimum_pressure`` is not positive, if it is above
            ``maximum_pressure``, or if the fluid names no ``water`` - each of which makes every
            point of the line the same failure.
        OutOfRangeError: from the range checks.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    low_si = float(input_to_si(spec, "minimum_pressure", minimum_pressure))
    high_si = float(input_to_si(spec, "maximum_pressure", maximum_pressure))
    apply_checks(
        checks.on_input,
        {"minimum_pressure": low_si, "maximum_pressure": high_si}.get,
        warnings,
    )

    if low_si <= 0.0 or high_si <= 0.0:
        raise InvalidInputError(
            "minimum_pressure",
            "a pressure is a divisor in the reduced variables and both ends of the sweep are "
            "finite and positive by definition",
        )
    if low_si > high_si:
        raise InvalidInputError(
            "minimum_pressure",
            f"the sweep starts at {low_si} Pa and ends at {high_si} Pa, so it runs backwards",
        )

    step = (high_si - low_si) / (NUMBER_OF_POINTS - 1)
    temperatures: list[Q] = []
    pressures: list[Q] = []
    failed = 0
    # The class's own last-resort behaviour: a point that does not solve keeps the temperature
    # the fluid was already at, which is the previous point's.
    last = float("nan")
    for index in range(NUMBER_OF_POINTS):
        p_si = low_si + step * index
        try:
            point = water_dew_point_multiphase(components, z, ureg.Quantity(p_si, "Pa"), eos)
        except AzothError:
            failed += 1
        else:
            last = float(point.temperature.to("K").magnitude)
            warnings.extend(point.warnings)
        temperatures.append(from_si(last, "K"))
        pressures.append(from_si(p_si, "Pa"))

    return WaterDewPointLineResult(
        temperature=tuple(temperatures),
        pressure=tuple(pressures),
        failed_points=failed,
        warnings=tuple(warnings),
    )
