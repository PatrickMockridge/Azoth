"""``eos.water_dew_point_multiphase``: the Python reference.

The mirror of ``crates/azoth-eos/src/water_dew_point_multiphase.rs``. NeqSim's
``WaterDewPointTemperatureMultiphaseFlash``: start at ``500`` K, step **down** by ``20`` K
until an aqueous phase appears - up instead if ``500`` K already has one - then bisect to
``1e-4`` K and answer the bracket's midpoint.

**It brackets on a multiphase flash finding an aqueous phase**, which is why the id exists
only now: that needed ``eos.multiflash_stability_test``, ``TPmultiflash``'s own
``stabilityAnalysis``, before the flash answered anything but a single gas phase here.

**The two caps are upstream's own counter**, not one each: the coarse search runs while
``i < 50`` and the bisection while ``i < 100``, over a counter neither resets.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WaterDewPointMultiphaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import reduced_parameters
from azoth.eos.reference.hydrate_inhibitor_wt import AQUEOUS, phase_label
from azoth.eos.reference.tp_multiflash import tp_multiflash

MODEL_ID = "eos.water_dew_point_multiphase"

#: The temperature the coarse search starts from, upstream's ``Thigh``.
START_TEMPERATURE = 500.0

#: The floor the coarse search steps down to, upstream's ``Tlow``.
FLOOR_TEMPERATURE = 200.0

#: The temperature it expands to when ``500`` K already has an aqueous phase.
CEILING_TEMPERATURE = 600.0

#: The coarse step, upstream's ``step``.
STEP = 20.0

#: The bracket width the bisection stops at, upstream's ``tolerance``.
TOLERANCE = 1.0e-4

#: The coarse search's step cap, upstream's ``i < 50``.
COARSE_CAP = 50

#: The bisection's step cap, upstream's ``i < 100`` - **the same counter**.
TOTAL_CAP = 100


def _aqueous(
    mixture: Mixture, temperature: float, pressure: float, z: list[float]
) -> tuple[bool, int]:
    """Whether the fluid at a state has an aqueous phase, and how many phases it has.

    The label is NeqSim's own ``PhaseEos.init`` rule, which
    :mod:`azoth.eos.reference.hydrate_inhibitor_wt` already carries.
    """
    from azoth.core.units import ureg

    result = tp_multiflash(
        mixture,
        ureg.Quantity(temperature, "K"),
        ureg.Quantity(pressure, "Pa"),
        z,
    )
    reduced = reduced_parameters(mixture, temperature, pressure)
    aqueous = any(
        phase_label(mixture, reduced, list(row), z_factor) == AQUEOUS
        for row, z_factor in zip(result.x, result.z_factor, strict=True)
    )
    return aqueous, result.phase_count


def water_dew_point_multiphase(
    components: list[str], z: list[float], P: Q, eos: str = "srk"
) -> WaterDewPointMultiphaseResult:
    """The temperature a gas condenses water at, by NeqSim's bracketing route.

    Raises:
        InvalidInputError: if ``z`` does not match ``components``, if the fluid names no
            ``water``, if that water carries a zero mole fraction, or if no aqueous phase
            appears anywhere in ``200`` to ``500`` K - which upstream answers with the floor.
        OutOfRangeError: from a trial flash.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    p_si = float(input_to_si(spec, "P", P))
    apply_checks(checks.on_input, {"P": p_si}.get, warnings)

    if len(components) != len(z):
        raise InvalidInputError("z", f"{len(components)} components and {len(z)} mole fractions")

    from azoth.eos import components as databank

    mixture, _ = databank.mixture_of(components, eos=eos)
    lowered = [name.strip().lower() for name in components]
    if "water" not in lowered:
        raise InvalidInputError(
            "components",
            "the fluid names no 'water', and a water dew point is the temperature an aqueous "
            "phase appears at - without water there is none to find",
        )
    water = lowered.index("water")
    if not z[water] > 0.0:
        raise InvalidInputError(
            "z",
            f"'water' is present and its mole fraction is {z[water]}, so there is no aqueous "
            f"phase for the bracket to find",
        )

    iterations = 0
    high = START_TEMPERATURE
    low = FLOOR_TEMPERATURE

    aqueous, _ = _aqueous(mixture, high, p_si, z)
    iterations += 1
    if aqueous:
        low = high
        high = CEILING_TEMPERATURE
        _aqueous(mixture, high, p_si, z)
        iterations += 1
    else:
        temperature = high - STEP
        found = False
        while temperature >= low and iterations < COARSE_CAP:
            aqueous, _ = _aqueous(mixture, temperature, p_si, z)
            iterations += 1
            if aqueous:
                low = temperature
                high = temperature + STEP
                found = True
                break
            temperature -= STEP
        if not found:
            raise InvalidInputError(
                "T",
                f"no aqueous phase appears anywhere in {FLOOR_TEMPERATURE} to "
                f"{START_TEMPERATURE} K at this pressure, so there is no water dew point to "
                f"bracket. Upstream returns the floor, which reads as an answer",
            )

    while (high - low) > TOLERANCE and iterations < TOTAL_CAP:
        middle = 0.5 * (high + low)
        aqueous, _ = _aqueous(mixture, middle, p_si, z)
        if aqueous:
            low = middle
        else:
            high = middle
        iterations += 1

    answer = 0.5 * (high + low)
    _, phase_count = _aqueous(mixture, answer, p_si, z)
    return WaterDewPointMultiphaseResult(
        temperature=from_si(answer, "K"),
        bracket_width=from_si(high - low, "K"),
        phase_count=phase_count,
        iterations=iterations,
        warnings=tuple(warnings),
    )
