"""``eos.hydrate_equilibrium_line``: the Python reference.

The mirror of ``crates/azoth-eos/src/hydrate_equilibrium_line.rs``. It is a grid over
:func:`azoth.eos.reference.hydrate_formation_temperature` and nothing else: ten equally
spaced pressures from a minimum to a maximum, and that model solved at each.

**The seeding NeqSim does between points is not here, and that is measured.** NeqSim sets the
system's temperature to the previous point's answer before each solve; the ten pressures
solved independently on a fresh system agree with its seeded walk to every printed digit
(``validation/neqsim/captures/hydrate_equilibrium_probe.tsv``, ``worst_relative_gap = 0``),
so the seeding is a starting guess rather than part of the answer.

**The count is ten and is not an input.** NeqSim's ``HydrateEquilibriumLine.numberOfPoints``
is a field initialised to ten that no constructor or setter writes, so its own surface cannot
vary it either. A caller wanting another grid wants
:func:`azoth.eos.reference.hydrate_formation_temperature` over pressures of their own.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HydrateEquilibriumLineResult
from azoth.core.units import Q, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos.reference import _hydrate
from azoth.eos.reference.hydrate_formation_temperature import hydrate_formation_temperature

MODEL_ID = "eos.hydrate_equilibrium_line"

#: The points NeqSim's line always solves, whatever the bounds.
POINTS = 10


def hydrate_equilibrium_line(
    components: list[str],
    P_min: Q,
    P_max: Q,
    z: list[float],
    eos: str = "srk",
    hydrate_model: str = _hydrate.PVTSIM,
) -> HydrateEquilibriumLineResult:
    """The formation temperature at ten pressures between a minimum and a maximum.

    Raises:
        InvalidInputError: if either bound is not positive or the maximum is not above the
            minimum.
        OutOfRangeError: from any point's own solve.
        SolverNotConvergedError: from any point's own solve, where NeqSim's loop would catch
            it and report the previous point's temperature at that index instead.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    p_min = input_to_si(spec, "P_min", P_min)
    p_max = input_to_si(spec, "P_max", P_max)
    apply_checks(checks.on_input, {"P_min": p_min, "P_max": p_max}.get, warnings)
    if p_min <= 0.0 or p_max <= 0.0 or p_max <= p_min:
        raise InvalidInputError(
            "P_max",
            f"the grid runs from a positive minimum up to a larger maximum, and this was "
            f"given {p_min} and {p_max} Pa",
        )

    # ``min + dp i`` over ``0..points``, which lands on the maximum at the last point -
    # NeqSim's own stepping, and why it never sets the maximum separately.
    dp = (p_max - p_min) / (POINTS - 1)
    temperature: list[float] = []
    pressure: list[float] = []
    for point in range(POINTS):
        p = p_min + dp * point
        solved = hydrate_formation_temperature(
            components, P=quantity(p, "Pa"), z=z, eos=eos, hydrate_model=hydrate_model
        )
        temperature.append(float(solved.temperature.to("K").magnitude))
        pressure.append(p)
        warnings.extend(solved.warnings)

    return HydrateEquilibriumLineResult(
        temperature=tuple(temperature),
        pressure=tuple(pressure),
        warnings=tuple(warnings),
    )
