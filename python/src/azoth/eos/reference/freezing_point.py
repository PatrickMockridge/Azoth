"""``eos.freezing_point`` - the temperature at which a pure substance's solid appears.

The Python reference, the mirror of ``crates/azoth-eos/src/freezing_point.rs``:

```text
solve  (g_fluid(T, P) - g_solid(T, P)) / (R T) = 0   for T
```

Two things the Rust module's docstring states and a reader of *this* file needs in the same
breath, because they are where a port goes wrong silently:

* **the solid is calibrated, not raw.** ``eos.parahydrogen_solid_phase`` is the thesis
  equation as it stands and its reference is fixed elsewhere - by NeqSim's
  ``SystemLeachmanEos``, which shifts it to the para-Leachman liquid at the triple point.
  Unshifted it misses the triple point's own definition, the two Gibbs energies equal there,
  by ``3.12`` J/mol.
* **the fluid's root is NeqSim's rule**: a gas below the triple-point pressure and a liquid at
  and above it, which are different states at the same temperature and pressure.
"""

from __future__ import annotations

from azoth.core.errors import (
    InvalidInputError,
    OutOfRangeError,
    SolverNotConvergedError,
)
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FreezingPointResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference import hydrogen_phase as fluid
from azoth.eos.reference.parahydrogen_solid_phase import parahydrogen_solid_phase

MODEL_ID = "eos.freezing_point"

#: NeqSim's ``ThermodynamicConstantsInterface.R``, which its residual divides by. **Not** the
#: equation's own: the two reference equations carry ``8.31451`` and ``8.3144621``, each a
#: rounded value of its own era, and this is the third.
R = 8.3144621

#: The tolerance on the dimensionless residual and on the bracket.
RESIDUAL_TOLERANCE = 1.0e-10
BRACKET_TOLERANCE = 1.0e-12

#: How far the search may walk, in kelvin, and how many steps it may take.
MINIMUM_TEMPERATURE = 0.5
MAXIMUM_TEMPERATURE = 5000.0
SPAN_GROWTH = 1.8
MAXIMUM_BRACKET_STEPS = 50
MAXIMUM_SOLVER_STEPS = 100

#: The triple point, and the enthalpy of fusion there: the coordinates the calibration is
#: measured at. Para-hydrogen's, from NeqSim's ``ParaHydrogenSolidHelmholtzEquation``.
TRIPLE_POINT_TEMPERATURE = 13.8033
TRIPLE_POINT_PRESSURE = 7042.0
TRIPLE_POINT_ENTHALPY_OF_FUSION = 118.0

#: The substances this has a solid equation for. One name, because
#: ``thermo/util/solid/`` holds two equations and the other is argon's.
SUBSTANCES = ("para-hydrogen",)


def _solid_state(t: float, p: float) -> tuple[float, float]:
    """The raw solid's molar Gibbs energy and entropy at a state, in J/mol and J/(mol*K)."""
    state = parahydrogen_solid_phase(T=from_si(t, "K"), P=from_si(p, "Pa"))
    return state.g.to("J/mol").magnitude, state.s.to("J/(mol*K)").magnitude


def _fluid_state(t: float, p: float, root: str) -> tuple[float, float]:
    """The fluid's molar Gibbs energy and entropy, on the root named, at a state."""
    density = fluid._solve_density(t, p, "para") if root == "gas" else _dense_density(t, p)
    properties = fluid._properties(t, density, "para")
    return properties["g"], properties["s"]


def _dense_density(t: float, p: float) -> float:
    density = fluid._solve_density_dense(t, p, "para")
    if density is None:
        raise OutOfRangeError(
            "P",
            p,
            "the freezing search asks for the fluid's dense root, and this pressure has none "
            "at this temperature within the range the equation is fitted to",
        )
    return density


def calibration() -> tuple[float, float]:
    """The two constants NeqSim's ``SystemLeachmanEos`` computes at the triple point.

    ``(gibbs_shift, entropy_shift)``, in J/mol and J/(mol*K). Measurements of the two models
    rather than tunables: the first comes out ``-3.122683``, which is the offset the capture
    measures between NeqSim's solid and the raw one.
    """
    g_liquid, s_liquid = _fluid_state(TRIPLE_POINT_TEMPERATURE, TRIPLE_POINT_PRESSURE, "liquid")
    g_raw, s_raw = _solid_state(TRIPLE_POINT_TEMPERATURE, TRIPLE_POINT_PRESSURE)
    return (
        g_liquid - g_raw,
        s_liquid - s_raw - TRIPLE_POINT_ENTHALPY_OF_FUSION / TRIPLE_POINT_TEMPERATURE,
    )


def calibrated_solid_gibbs(t: float, p: float, shifts: tuple[float, float]) -> float:
    """The solid's molar Gibbs energy at a state, on the reference the liquid is on."""
    gibbs_shift, entropy_shift = shifts
    g_raw, _ = _solid_state(t, p)
    # A shift of `-S (T - T_tp) + G` in the Helmholtz energy is the same shift in the Gibbs
    # energy: it carries no volume, so `g = a + P v` passes it through.
    return g_raw - entropy_shift * (t - TRIPLE_POINT_TEMPERATURE) + gibbs_shift


def residual(t: float, p: float, root: str, shifts: tuple[float, float]) -> float:
    """The dimensionless residual the solve drives to zero."""
    g_fluid, _ = _fluid_state(t, p, root)
    return (g_fluid - calibrated_solid_gibbs(t, p, shifts)) / (R * t)


def freezing_point(components: list[str], P: Q) -> FreezingPointResult:
    """The freezing-point temperature of para-hydrogen at a pressure.

    Raises:
        InvalidInputError: if ``components`` is not one entry, or names a substance this has
            no solid equation for.
        OutOfRangeError: if ``P`` is not positive.
        SolverNotConvergedError: if the residual will not bracket or the bracket collapses
            without reaching the tolerance.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"P": p_si}.get, warnings)

    if len(components) != 1:
        raise InvalidInputError(
            "components",
            f"{len(components)} entries; a freezing point is a pure substance's, and this "
            f"crate has a solid equation for para-hydrogen alone",
        )
    name = components[0].strip().lower()
    if name not in SUBSTANCES:
        raise InvalidInputError(
            "components",
            f"{components[0]!r}; `thermo/util/solid/` gives an equation for para-hydrogen "
            f"(and argon, which is not a hydrogen phase), so that is the name this takes",
        )

    algorithm = spec["algorithm"]
    # NeqSim's own rule, from `initializeCoexistingFluidPhase`: the fluid is a gas below the
    # triple-point pressure and a liquid at and above it.
    root = "gas" if p_si < TRIPLE_POINT_PRESSURE else "liquid"
    shifts = calibration()

    start = float(algorithm.get("initial_temperature", 14.0))
    start_residual = residual(start, p_si, root, shifts)
    if abs(start_residual) < RESIDUAL_TOLERANCE:
        return _result(start, 1, start_residual, warnings)

    span = max(start * 0.01, 0.1)
    low, low_residual = start, start_residual
    high, high_residual = start, start_residual
    iterations = 1
    bracketed = False
    for _ in range(MAXIMUM_BRACKET_STEPS):
        candidate_low = max(start - span, MINIMUM_TEMPERATURE)
        candidate_high = min(start + span, MAXIMUM_TEMPERATURE)
        # A trial the equation cannot evaluate is skipped rather than fatal, which is NeqSim's
        # `evaluateResidualIfValid`: near the triple point one side of a span can ask for a
        # dense root that does not exist.
        try:
            low, low_residual = candidate_low, residual(candidate_low, p_si, root, shifts)
            iterations += 1
        except OutOfRangeError:
            pass
        try:
            high, high_residual = candidate_high, residual(candidate_high, p_si, root, shifts)
            iterations += 1
        except OutOfRangeError:
            pass
        if low_residual * high_residual <= 0.0:
            bracketed = True
            break
        span *= SPAN_GROWTH

    if not bracketed:
        raise _not_converged(iterations, start_residual, algorithm)

    for _ in range(MAXIMUM_SOLVER_STEPS):
        mid = 0.5 * (low + high)
        value = residual(mid, p_si, root, shifts)
        iterations += 1
        if abs(value) < RESIDUAL_TOLERANCE or high - low < BRACKET_TOLERANCE:
            return _result(mid, iterations, value, warnings)
        if low_residual * value <= 0.0:
            high = mid
        else:
            low, low_residual = mid, value

    raise _not_converged(iterations, 0.5 * (low + high), algorithm)


def _result(
    temperature: float, iterations: int, value: float, warnings: list[Warning]
) -> FreezingPointResult:
    return FreezingPointResult(
        temperature=from_si(temperature, "K"),
        iterations=iterations,
        residual=value,
        warnings=tuple(warnings),
    )


def _not_converged(iterations: int, value: float, algorithm: dict[str, object]) -> Exception:
    return SolverNotConvergedError(
        iterations=iterations,
        residual=value,
        tolerance=float(str(algorithm["tolerance"])),
    )


__all__ = ["freezing_point"]
